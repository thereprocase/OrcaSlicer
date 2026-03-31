// gpu_collision.cpp -- CollisionEvaluator implementations (CPU + GPU)
//
// CPU path: direct port of the inline evaluate() logic from snuggle_nester.hpp.
// GPU path: OpenGL 4.3 compute shader, one work group per individual.

#include "gpu_collision.hpp"
#include "snuggle_nester.hpp"

#include <boost/log/trivial.hpp>
#include <cmath>
#include <cstring>

#ifdef SLIC3R_GUI
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifdef _WIN32
#include <windows.h>
#endif
#include <GL/glew.h>
#endif // SLIC3R_GUI

namespace snuggle {

// =====================================================================
// CpuCollisionEvaluator
// =====================================================================

void CpuCollisionEvaluator::upload_grids(
    const std::vector<PartInfo>& /*parts*/,
    const std::vector<std::vector<VoxelGrid>>& rot_cache)
{
    rot_cache_ = &rot_cache;
}

const VoxelGrid& CpuCollisionEvaluator::get_rotated(size_t part_idx, float angle) const
{
    int bin = (int)std::floor(angle * ROT_CACHE_BINS / (2.0f * 3.14159265f));
    bin = ((bin % ROT_CACHE_BINS) + ROT_CACHE_BINS) % ROT_CACHE_BINS;
    return (*rot_cache_)[part_idx][bin];
}

void CpuCollisionEvaluator::evaluate_batch(
    std::vector<Individual>& individuals,
    const std::vector<PartInfo>& parts,
    float bed_w, float bed_h, float bed_margin)
{
    if (!rot_cache_) return;
    size_t n = parts.size();

    for (auto& ind : individuals) {
        ind.collision_count = 0;
        ind.oob_count = 0;

        std::vector<const VoxelGrid*> rotated(n);
        for (size_t i = 0; i < n; i++)
            rotated[i] = &get_rotated(i, ind.placements[i].zrot);

        for (size_t i = 0; i < n; i++) {
            const auto& pi = ind.placements[i];
            Vec3f off_i = {pi.x, pi.y, 0.0f};

            for (size_t j = i + 1; j < n; j++) {
                const auto& pj = ind.placements[j];
                Vec3f off_j = {pj.x, pj.y, 0.0f};
                size_t c = VoxelGrid::collision_count(*rotated[i], off_i, *rotated[j], off_j);
                ind.collision_count += c;
            }

            const auto& g = *rotated[i];
            float part_min_x = pi.x + g.origin.x;
            float part_min_y = pi.y + g.origin.y;
            float part_max_x = pi.x + g.origin.x + g.nx * g.voxel_size;
            float part_max_y = pi.y + g.origin.y + g.ny * g.voxel_size;

            if (part_min_x < bed_margin)
                ind.oob_count += (size_t)((bed_margin - part_min_x) / g.voxel_size);
            if (part_min_y < bed_margin)
                ind.oob_count += (size_t)((bed_margin - part_min_y) / g.voxel_size);
            if (part_max_x > bed_w - bed_margin)
                ind.oob_count += (size_t)((part_max_x - bed_w + bed_margin) / g.voxel_size);
            if (part_max_y > bed_h - bed_margin)
                ind.oob_count += (size_t)((part_max_y - bed_h + bed_margin) / g.voxel_size);
        }
    }
}

// =====================================================================
// GpuCollisionEvaluator (only when SLIC3R_GUI is defined)
// =====================================================================

#ifdef SLIC3R_GUI

// Metadata struct matching the shader layout (48 bytes, std430)
struct alignas(4) GridMeta {
    uint32_t data_offset;  // byte offset into voxel SSBO
    uint32_t nx, ny, nz;
    float    voxel_size;
    float    origin_x, origin_y, origin_z;
    uint32_t _pad[4];      // pad to 48 bytes
};
static_assert(sizeof(GridMeta) == 48, "GridMeta must be 48 bytes for shader alignment");

// Placement struct matching the shader layout (16 bytes)
struct alignas(4) PlacementGPU {
    float x, y, zrot;
    float _pad;
};
static_assert(sizeof(PlacementGPU) == 16, "PlacementGPU must be 16 bytes");

// ── Compute shader source ────────────────────────────────────
// POC version: one invocation per individual, serial pair iteration.
// A future version would use work-group parallelism across pairs.
static const char* const COLLISION_SHADER_SRC = R"GLSL(
#version 430
layout(local_size_x = 1, local_size_y = 1, local_size_z = 1) in;

struct GridMeta {
    uint data_offset;
    uint nx, ny, nz;
    float voxel_size;
    float origin_x, origin_y, origin_z;
    uint _pad0, _pad1, _pad2, _pad3;
};

struct PlacementData {
    float x, y, zrot;
    float _pad;
};

layout(std430, binding = 0) readonly buffer VoxelData { uint voxel_bits[]; };
layout(std430, binding = 1) readonly buffer MetaData  { GridMeta metas[];  };
layout(std430, binding = 2) readonly buffer Placements { PlacementData placements[]; };
layout(std430, binding = 3) writeonly buffer Results   { uvec2 results[];  };

uniform uint u_n_parts;
uniform float u_bed_w;
uniform float u_bed_h;
uniform float u_bed_margin;
uniform uint u_rot_bins;

uint angle_to_bin(float zrot) {
    int bin = int(floor(zrot * float(u_rot_bins) / 6.2831853));
    return uint(((bin % int(u_rot_bins)) + int(u_rot_bins)) % int(u_rot_bins));
}

bool voxel_get(GridMeta m, int gx, int gy, int gz) {
    if (gx < 0 || gy < 0 || gz < 0) return false;
    if (uint(gx) >= m.nx || uint(gy) >= m.ny || uint(gz) >= m.nz) return false;
    uint idx = uint(gx) + uint(gy) * m.nx + uint(gz) * m.nx * m.ny;
    uint word_idx = m.data_offset / 4u + idx / 32u;
    uint bit_idx = idx % 32u;
    return (voxel_bits[word_idx] & (1u << bit_idx)) != 0u;
}

void main() {
    uint ind_idx = gl_GlobalInvocationID.x;
    uint base = ind_idx * u_n_parts;

    uint total_collisions = 0u;
    uint total_oob = 0u;

    for (uint i = 0u; i < u_n_parts; i++) {
        PlacementData pi = placements[base + i];
        uint rot_bin_i = angle_to_bin(pi.zrot);
        uint meta_idx_i = i * u_rot_bins + rot_bin_i;
        GridMeta mi = metas[meta_idx_i];

        // Bounds check
        float pmin_x = pi.x + mi.origin_x;
        float pmin_y = pi.y + mi.origin_y;
        float pmax_x = pmin_x + float(mi.nx) * mi.voxel_size;
        float pmax_y = pmin_y + float(mi.ny) * mi.voxel_size;

        if (pmin_x < u_bed_margin) total_oob += uint((u_bed_margin - pmin_x) / mi.voxel_size);
        if (pmin_y < u_bed_margin) total_oob += uint((u_bed_margin - pmin_y) / mi.voxel_size);
        if (pmax_x > u_bed_w - u_bed_margin) total_oob += uint((pmax_x - u_bed_w + u_bed_margin) / mi.voxel_size);
        if (pmax_y > u_bed_h - u_bed_margin) total_oob += uint((pmax_y - u_bed_h + u_bed_margin) / mi.voxel_size);

        // Pairwise collision with parts j > i
        for (uint j = i + 1u; j < u_n_parts; j++) {
            PlacementData pj = placements[base + j];
            uint rot_bin_j = angle_to_bin(pj.zrot);
            uint meta_idx_j = j * u_rot_bins + rot_bin_j;
            GridMeta mj = metas[meta_idx_j];

            // World-space bounds for overlap test
            float ai_min_x = pi.x + mi.origin_x;
            float ai_min_y = pi.y + mi.origin_y;
            float aj_min_x = pj.x + mj.origin_x;
            float aj_min_y = pj.y + mj.origin_y;

            float ox_min = max(ai_min_x, aj_min_x);
            float oy_min = max(ai_min_y, aj_min_y);
            float oz_min = max(mi.origin_z, mj.origin_z);
            float ox_max = min(ai_min_x + float(mi.nx) * mi.voxel_size,
                               aj_min_x + float(mj.nx) * mj.voxel_size);
            float oy_max = min(ai_min_y + float(mi.ny) * mi.voxel_size,
                               aj_min_y + float(mj.ny) * mj.voxel_size);
            float oz_max = min(mi.origin_z + float(mi.nz) * mi.voxel_size,
                               mj.origin_z + float(mj.nz) * mj.voxel_size);

            if (ox_min >= ox_max || oy_min >= oy_max || oz_min >= oz_max) continue;

            float vs = max(mi.voxel_size, mj.voxel_size);

            for (float wz = oz_min + vs * 0.5; wz < oz_max; wz += vs) {
                for (float wy = oy_min + vs * 0.5; wy < oy_max; wy += vs) {
                    for (float wx = ox_min + vs * 0.5; wx < ox_max; wx += vs) {
                        int ax = int(floor((wx - ai_min_x) / mi.voxel_size));
                        int ay = int(floor((wy - ai_min_y) / mi.voxel_size));
                        int az = int(floor((wz - mi.origin_z) / mi.voxel_size));

                        if (voxel_get(mi, ax, ay, az)) {
                            int bx = int(floor((wx - aj_min_x) / mj.voxel_size));
                            int by = int(floor((wy - aj_min_y) / mj.voxel_size));
                            int bz = int(floor((wz - mj.origin_z) / mj.voxel_size));

                            if (voxel_get(mj, bx, by, bz))
                                total_collisions++;
                        }
                    }
                }
            }
        }
    }

    results[ind_idx] = uvec2(total_collisions, total_oob);
}
)GLSL";

// ── GL context initialization (Windows) ──────────────────────

GpuCollisionEvaluator::GpuCollisionEvaluator()
{
    available_ = init_context() && compile_shader();
    if (available_) {
        BOOST_LOG_TRIVIAL(info) << "Snuggle GPU: compute shader backend initialized";
    } else {
        BOOST_LOG_TRIVIAL(warning) << "Snuggle GPU: initialization failed, will use CPU fallback";
        cleanup();
    }
}

GpuCollisionEvaluator::~GpuCollisionEvaluator()
{
    cleanup();
}

bool GpuCollisionEvaluator::init_context()
{
#ifdef _WIN32
    WNDCLASSA wc = {};
    wc.lpfnWndProc = DefWindowProcA;
    wc.hInstance = GetModuleHandle(NULL);
    wc.lpszClassName = "SnuggleGPUCollision";
    if (!RegisterClassA(&wc)) {
        // Class may already be registered from a previous instance
        DWORD err = GetLastError();
        if (err != ERROR_CLASS_ALREADY_EXISTS) {
            BOOST_LOG_TRIVIAL(warning) << "Snuggle GPU: RegisterClass failed, error " << err;
            return false;
        }
    }

    HWND hwnd = CreateWindowA("SnuggleGPUCollision", "", 0,
                              0, 0, 1, 1, NULL, NULL, wc.hInstance, NULL);
    if (!hwnd) {
        BOOST_LOG_TRIVIAL(warning) << "Snuggle GPU: CreateWindow failed";
        return false;
    }
    gl_hwnd_ = (void*)hwnd;

    HDC hdc = GetDC(hwnd);
    if (!hdc) {
        BOOST_LOG_TRIVIAL(warning) << "Snuggle GPU: GetDC failed";
        return false;
    }
    gl_dc_ = (void*)hdc;

    PIXELFORMATDESCRIPTOR pfd = {};
    pfd.nSize = sizeof(pfd);
    pfd.nVersion = 1;
    pfd.dwFlags = PFD_SUPPORT_OPENGL;
    pfd.iPixelType = PFD_TYPE_RGBA;
    pfd.cColorBits = 32;

    int fmt = ChoosePixelFormat(hdc, &pfd);
    if (!fmt || !SetPixelFormat(hdc, fmt, &pfd)) {
        BOOST_LOG_TRIVIAL(warning) << "Snuggle GPU: SetPixelFormat failed";
        return false;
    }

    HGLRC ctx = wglCreateContext(hdc);
    if (!ctx) {
        BOOST_LOG_TRIVIAL(warning) << "Snuggle GPU: wglCreateContext failed";
        return false;
    }
    gl_context_ = (void*)ctx;

    if (!wglMakeCurrent(hdc, ctx)) {
        BOOST_LOG_TRIVIAL(warning) << "Snuggle GPU: wglMakeCurrent failed";
        return false;
    }

    // Initialize GLEW
    glewExperimental = GL_TRUE;
    GLenum glew_err = glewInit();
    if (glew_err != GLEW_OK) {
        BOOST_LOG_TRIVIAL(warning) << "Snuggle GPU: glewInit failed: "
                                   << glewGetErrorString(glew_err);
        return false;
    }

    // Check GL 4.3+ (compute shader support)
    GLint major = 0, minor = 0;
    glGetIntegerv(GL_MAJOR_VERSION, &major);
    glGetIntegerv(GL_MINOR_VERSION, &minor);
    BOOST_LOG_TRIVIAL(info) << "Snuggle GPU: OpenGL " << major << "." << minor
                            << " (" << glGetString(GL_RENDERER) << ")";

    if (major < 4 || (major == 4 && minor < 3)) {
        BOOST_LOG_TRIVIAL(warning) << "Snuggle GPU: OpenGL 4.3+ required for compute shaders, got "
                                   << major << "." << minor;
        return false;
    }

    return true;
#else
    // Non-Windows platforms: not implemented yet
    BOOST_LOG_TRIVIAL(info) << "Snuggle GPU: compute context not implemented for this platform";
    return false;
#endif
}

bool GpuCollisionEvaluator::compile_shader()
{
    GLuint shader = glCreateShader(GL_COMPUTE_SHADER);
    if (!shader) {
        BOOST_LOG_TRIVIAL(warning) << "Snuggle GPU: glCreateShader failed";
        return false;
    }

    glShaderSource(shader, 1, &COLLISION_SHADER_SRC, nullptr);
    glCompileShader(shader);

    GLint status = 0;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &status);
    if (status != GL_TRUE) {
        char log[2048] = {};
        glGetShaderInfoLog(shader, sizeof(log), nullptr, log);
        BOOST_LOG_TRIVIAL(warning) << "Snuggle GPU: shader compile failed:\n" << log;
        glDeleteShader(shader);
        return false;
    }

    program_ = glCreateProgram();
    glAttachShader(program_, shader);
    glLinkProgram(program_);

    glGetProgramiv(program_, GL_LINK_STATUS, &status);
    if (status != GL_TRUE) {
        char log[2048] = {};
        glGetProgramInfoLog(program_, sizeof(log), nullptr, log);
        BOOST_LOG_TRIVIAL(warning) << "Snuggle GPU: program link failed:\n" << log;
        glDeleteShader(shader);
        glDeleteProgram(program_);
        program_ = 0;
        return false;
    }

    glDeleteShader(shader);
    BOOST_LOG_TRIVIAL(info) << "Snuggle GPU: compute shader compiled and linked";
    return true;
}

void GpuCollisionEvaluator::upload_grids(
    const std::vector<PartInfo>& parts,
    const std::vector<std::vector<VoxelGrid>>& rot_cache)
{
    if (!available_) return;

    n_parts_ = parts.size();
    size_t total_metas = n_parts_ * ROT_BINS;

    // Pass 1: compute total voxel data size (repacked to uint32 alignment)
    size_t total_voxel_bytes = 0;
    for (size_t p = 0; p < n_parts_; p++) {
        for (int r = 0; r < ROT_BINS; r++) {
            const auto& grid = rot_cache[p][r];
            size_t total_bits = grid.nx * grid.ny * grid.nz;
            // Round up to 4-byte boundary for uint32 word alignment
            size_t grid_bytes = ((total_bits + 31) / 32) * 4;
            total_voxel_bytes += grid_bytes;
        }
    }

    BOOST_LOG_TRIVIAL(info) << "Snuggle GPU: uploading " << n_parts_ << " parts x "
                            << ROT_BINS << " rotations = "
                            << (total_voxel_bytes / (1024 * 1024)) << " MB voxel data";

    // Allocate CPU-side buffers
    std::vector<uint8_t> voxel_data(total_voxel_bytes, 0);
    std::vector<GridMeta> meta_data(total_metas);

    // Pass 2: repack uint8 bit-packed data into uint32 bit-packed data
    // Source: bit N is at byte[N/8], bit (N%8)
    // Dest:   bit N is at uint32[N/32], bit (N%32)
    size_t offset = 0;
    for (size_t p = 0; p < n_parts_; p++) {
        for (int r = 0; r < ROT_BINS; r++) {
            const auto& grid = rot_cache[p][r];
            size_t total_bits = grid.nx * grid.ny * grid.nz;
            size_t grid_bytes = ((total_bits + 31) / 32) * 4;

            // Fill metadata
            size_t meta_idx = p * ROT_BINS + r;
            meta_data[meta_idx].data_offset = (uint32_t)offset;
            meta_data[meta_idx].nx = (uint32_t)grid.nx;
            meta_data[meta_idx].ny = (uint32_t)grid.ny;
            meta_data[meta_idx].nz = (uint32_t)grid.nz;
            meta_data[meta_idx].voxel_size = grid.voxel_size;
            meta_data[meta_idx].origin_x = grid.origin.x;
            meta_data[meta_idx].origin_y = grid.origin.y;
            meta_data[meta_idx].origin_z = grid.origin.z;
            std::memset(meta_data[meta_idx]._pad, 0, sizeof(meta_data[meta_idx]._pad));

            // Repack bits: read from uint8 source, write into uint32 dest
            // Both are linear bit arrays indexed by (x + y*nx + z*nx*ny),
            // but the byte/word packing differs.
            uint32_t* dest_words = reinterpret_cast<uint32_t*>(&voxel_data[offset]);
            const uint8_t* src = grid.bits_data();
            for (size_t bit = 0; bit < total_bits; bit++) {
                bool set = (src[bit / 8] >> (bit % 8)) & 1;
                if (set) {
                    dest_words[bit / 32] |= (1u << (bit % 32));
                }
            }

            offset += grid_bytes;
        }
    }

    // Upload voxel data SSBO (binding 0)
    if (voxel_ssbo_) glDeleteBuffers(1, &voxel_ssbo_);
    glGenBuffers(1, &voxel_ssbo_);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, voxel_ssbo_);
    glBufferData(GL_SHADER_STORAGE_BUFFER, total_voxel_bytes, voxel_data.data(), GL_STATIC_DRAW);

    // Upload metadata SSBO (binding 1)
    if (meta_ssbo_) glDeleteBuffers(1, &meta_ssbo_);
    glGenBuffers(1, &meta_ssbo_);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, meta_ssbo_);
    glBufferData(GL_SHADER_STORAGE_BUFFER, total_metas * sizeof(GridMeta), meta_data.data(), GL_STATIC_DRAW);

    // Pre-allocate placement and results SSBOs (will be resized in evaluate_batch)
    if (!placement_ssbo_) glGenBuffers(1, &placement_ssbo_);
    if (!results_ssbo_) glGenBuffers(1, &results_ssbo_);

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

    GLenum err = glGetError();
    if (err != GL_NO_ERROR) {
        BOOST_LOG_TRIVIAL(warning) << "Snuggle GPU: GL error after upload_grids: 0x"
                                   << std::hex << err << std::dec;
        available_ = false;
    } else {
        BOOST_LOG_TRIVIAL(info) << "Snuggle GPU: grid upload complete";
    }
}

void GpuCollisionEvaluator::evaluate_batch(
    std::vector<Individual>& individuals,
    const std::vector<PartInfo>& parts,
    float bed_w, float bed_h, float bed_margin)
{
    if (!available_ || individuals.empty()) return;

    size_t pop_size = individuals.size();
    size_t n = parts.size();

    // Build placement buffer: pop_size * n_parts PlacementGPU structs
    std::vector<PlacementGPU> placements(pop_size * n);
    for (size_t i = 0; i < pop_size; i++) {
        for (size_t p = 0; p < n; p++) {
            auto& dst = placements[i * n + p];
            dst.x = individuals[i].placements[p].x;
            dst.y = individuals[i].placements[p].y;
            dst.zrot = individuals[i].placements[p].zrot;
            dst._pad = 0.0f;
        }
    }

    // Upload placements (binding 2)
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, placement_ssbo_);
    glBufferData(GL_SHADER_STORAGE_BUFFER,
                 placements.size() * sizeof(PlacementGPU),
                 placements.data(), GL_DYNAMIC_DRAW);

    // Allocate results buffer (binding 3): one uvec2 per individual
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, results_ssbo_);
    glBufferData(GL_SHADER_STORAGE_BUFFER,
                 pop_size * 2 * sizeof(uint32_t),
                 nullptr, GL_DYNAMIC_READ);

    // Bind SSBOs
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, voxel_ssbo_);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, meta_ssbo_);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, placement_ssbo_);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, results_ssbo_);

    // Set uniforms
    glUseProgram(program_);
    glUniform1ui(glGetUniformLocation(program_, "u_n_parts"), (GLuint)n);
    glUniform1f(glGetUniformLocation(program_, "u_bed_w"), bed_w);
    glUniform1f(glGetUniformLocation(program_, "u_bed_h"), bed_h);
    glUniform1f(glGetUniformLocation(program_, "u_bed_margin"), bed_margin);
    glUniform1ui(glGetUniformLocation(program_, "u_rot_bins"), (GLuint)ROT_BINS);

    // Dispatch: one invocation per individual
    glDispatchCompute((GLuint)pop_size, 1, 1);
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

    // Read back results (implicit sync via glGetBufferSubData)
    std::vector<uint32_t> results(pop_size * 2);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, results_ssbo_);
    glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0,
                       results.size() * sizeof(uint32_t),
                       results.data());

    // Write results into individuals
    for (size_t i = 0; i < pop_size; i++) {
        individuals[i].collision_count = results[i * 2];
        individuals[i].oob_count = results[i * 2 + 1];
    }

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
    glUseProgram(0);

    GLenum err = glGetError();
    if (err != GL_NO_ERROR) {
        BOOST_LOG_TRIVIAL(warning) << "Snuggle GPU: GL error in evaluate_batch: 0x"
                                   << std::hex << err << std::dec
                                   << " -- falling back to CPU for remaining evaluations";
        available_ = false;
    }
}

void GpuCollisionEvaluator::cleanup()
{
    // Delete GL objects if context is still current
    if (gl_context_) {
#ifdef _WIN32
        HDC hdc = (HDC)gl_dc_;
        HGLRC ctx = (HGLRC)gl_context_;
        wglMakeCurrent(hdc, ctx);
#endif
        if (program_) { glDeleteProgram(program_); program_ = 0; }
        if (voxel_ssbo_) { glDeleteBuffers(1, &voxel_ssbo_); voxel_ssbo_ = 0; }
        if (meta_ssbo_) { glDeleteBuffers(1, &meta_ssbo_); meta_ssbo_ = 0; }
        if (placement_ssbo_) { glDeleteBuffers(1, &placement_ssbo_); placement_ssbo_ = 0; }
        if (results_ssbo_) { glDeleteBuffers(1, &results_ssbo_); results_ssbo_ = 0; }

#ifdef _WIN32
        wglMakeCurrent(NULL, NULL);
        wglDeleteContext(ctx);
        gl_context_ = nullptr;

        if (gl_hwnd_) {
            HWND hwnd = (HWND)gl_hwnd_;
            if (gl_dc_) {
                ReleaseDC(hwnd, (HDC)gl_dc_);
                gl_dc_ = nullptr;
            }
            DestroyWindow(hwnd);
            gl_hwnd_ = nullptr;
        }
#endif
    }
}

#endif // SLIC3R_GUI

// =====================================================================
// Factory
// =====================================================================

std::unique_ptr<CollisionEvaluator> create_collision_evaluator()
{
#ifdef SLIC3R_GUI
    auto gpu = std::make_unique<GpuCollisionEvaluator>();
    if (gpu->is_available()) {
        BOOST_LOG_TRIVIAL(warning) << "Snuggle: using GPU collision evaluator";
        return gpu;
    }
    BOOST_LOG_TRIVIAL(warning) << "Snuggle: GPU not available, using CPU collision evaluator";
#else
    BOOST_LOG_TRIVIAL(warning) << "Snuggle: headless build, using CPU collision evaluator";
#endif
    return std::make_unique<CpuCollisionEvaluator>();
}

} // namespace snuggle
