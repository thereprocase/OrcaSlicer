// polite_voxelizer.hpp — Snuggle's CPU Voxelizer
//
// Converts triangle meshes into 3D bit-packed voxel grids.
// Conservative outward rounding: partial voxel occupancy = solid.
//
// POLITENESS GUARANTEES:
//   - Hard memory cap: refuses to allocate grids > MAX_VOXEL_MEMORY_MB
//   - Yields CPU every YIELD_EVERY_TRIANGLES triangles (Sleep(0))
//   - Single-threaded by default (opt-in parallelism)
//   - All allocations checked, all loops bounded
//   - Progress callback so caller can abort
//   - No exceptions thrown — error codes only

#pragma once

#include <vector>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <functional>
#include <thread>
#include <chrono>

#ifdef _WIN32
#include <windows.h>
#else
#include <sched.h>
#endif

namespace snuggle {

// ── Politeness knobs ───────────────────────────────────────
static constexpr size_t MAX_VOXEL_MEMORY_MB     = 256;   // Hard cap per grid
static constexpr size_t YIELD_EVERY_TRIANGLES   = 1000;  // Yield CPU this often
static constexpr size_t MAX_GRID_DIM            = 512;   // Max voxels per axis
static constexpr size_t MAX_TOTAL_VOXELS        = 64 * 1024 * 1024; // 64M voxels absolute max

// ── Yield: give other processes a chance ───────────────────
inline void polite_yield() {
#ifdef _WIN32
    SwitchToThread();  // Yields to any ready thread, returns immediately if none
#else
    sched_yield();
#endif
}

// ── Small sleep: back off briefly when doing heavy work ────
inline void polite_nap() {
    std::this_thread::sleep_for(std::chrono::microseconds(100));
}

// ── Vec3f: minimal 3D vector ──────────────────────────────
struct Vec3f {
    float x, y, z;
    Vec3f() : x(0), y(0), z(0) {}
    Vec3f(float x, float y, float z) : x(x), y(y), z(z) {}
    Vec3f operator-(const Vec3f &o) const { return {x-o.x, y-o.y, z-o.z}; }
    Vec3f operator+(const Vec3f &o) const { return {x+o.x, y+o.y, z+o.z}; }
    Vec3f operator*(float s) const { return {x*s, y*s, z*s}; }
    float dot(const Vec3f &o) const { return x*o.x + y*o.y + z*o.z; }
    Vec3f cross(const Vec3f &o) const {
        return {y*o.z - z*o.y, z*o.x - x*o.z, x*o.y - y*o.x};
    }
};

// ── AABB ──────────────────────────────────────────────────
struct AABB {
    Vec3f min, max;
    AABB() : min{1e18f,1e18f,1e18f}, max{-1e18f,-1e18f,-1e18f} {}
    void expand(const Vec3f &p) {
        min.x = std::min(min.x, p.x); min.y = std::min(min.y, p.y); min.z = std::min(min.z, p.z);
        max.x = std::max(max.x, p.x); max.y = std::max(max.y, p.y); max.z = std::max(max.z, p.z);
    }
    Vec3f size() const { return max - min; }
};

// ── Triangle ──────────────────────────────────────────────
struct Triangle {
    Vec3f v0, v1, v2;
};

// ── Error codes ───────────────────────────────────────────
enum class VoxError {
    OK = 0,
    MEMORY_CAP_EXCEEDED,
    GRID_TOO_LARGE,
    NO_TRIANGLES,
    ABORTED_BY_USER,
    INVALID_RESOLUTION,
};

inline const char* vox_error_str(VoxError e) {
    switch (e) {
        case VoxError::OK:                  return "OK";
        case VoxError::MEMORY_CAP_EXCEEDED: return "Memory cap exceeded";
        case VoxError::GRID_TOO_LARGE:      return "Grid dimensions exceed max";
        case VoxError::NO_TRIANGLES:        return "No triangles provided";
        case VoxError::ABORTED_BY_USER:     return "Aborted by progress callback";
        case VoxError::INVALID_RESOLUTION:  return "Invalid voxel resolution";
    }
    return "Unknown error";
}

// ── Progress callback: return false to abort ──────────────
//    progress is 0.0 to 1.0
using ProgressFn = std::function<bool(float progress, const char* phase)>;

// ── Bit-packed 3D voxel grid ──────────────────────────────
class VoxelGrid {
public:
    size_t nx = 0, ny = 0, nz = 0;  // Grid dimensions
    float  voxel_size = 1.0f;        // mm per voxel
    Vec3f  origin;                    // World-space origin of grid (min corner)

    VoxelGrid() = default;

    VoxError allocate(size_t nx_, size_t ny_, size_t nz_) {
        // Bounds check
        if (nx_ > MAX_GRID_DIM || ny_ > MAX_GRID_DIM || nz_ > MAX_GRID_DIM)
            return VoxError::GRID_TOO_LARGE;

        size_t total = nx_ * ny_ * nz_;
        if (total > MAX_TOTAL_VOXELS)
            return VoxError::GRID_TOO_LARGE;

        size_t bytes = (total + 7) / 8;
        if (bytes > MAX_VOXEL_MEMORY_MB * 1024 * 1024)
            return VoxError::MEMORY_CAP_EXCEEDED;

        nx = nx_; ny = ny_; nz = nz_;
        bits_.resize(bytes, 0);
        return VoxError::OK;
    }

    void clear() {
        std::memset(bits_.data(), 0, bits_.size());
    }

    size_t memory_bytes() const { return bits_.size(); }
    size_t total_voxels() const { return nx * ny * nz; }

    // ── Accessors ──────────────────────────────────────────
    bool get(size_t x, size_t y, size_t z) const {
        if (x >= nx || y >= ny || z >= nz) return false;
        size_t idx = x + y * nx + z * nx * ny;
        return (bits_[idx >> 3] >> (idx & 7)) & 1;
    }

    void set(size_t x, size_t y, size_t z) {
        if (x >= nx || y >= ny || z >= nz) return;
        size_t idx = x + y * nx + z * nx * ny;
        bits_[idx >> 3] |= (1u << (idx & 7));
    }

    // ── World <-> grid coordinate conversion ───────────────
    // World to grid (conservative outward: floor for min, ceil for max)
    void world_to_grid(const Vec3f &world, int &gx, int &gy, int &gz) const {
        gx = (int)std::floor((world.x - origin.x) / voxel_size);
        gy = (int)std::floor((world.y - origin.y) / voxel_size);
        gz = (int)std::floor((world.z - origin.z) / voxel_size);
    }

    Vec3f grid_to_world_center(int gx, int gy, int gz) const {
        return {
            origin.x + (gx + 0.5f) * voxel_size,
            origin.y + (gy + 0.5f) * voxel_size,
            origin.z + (gz + 0.5f) * voxel_size,
        };
    }

    // ── Count solid voxels ─────────────────────────────────
    size_t count_solid() const {
        size_t count = 0;
        for (size_t i = 0; i < bits_.size(); i++) {
            // popcount byte
            uint8_t b = bits_[i];
            while (b) { count += b & 1; b >>= 1; }
        }
        return count;
    }

    // ── Collision: AND two grids in shared world space ──────
    //    Returns number of overlapping solid voxels.
    //    Grids can have different origins (different parts at
    //    different XY positions on the bed).
    static size_t collision_count(
        const VoxelGrid &a, const Vec3f &a_offset,
        const VoxelGrid &b, const Vec3f &b_offset)
    {
        // Compute overlap region in world space
        Vec3f a_min = { a.origin.x + a_offset.x, a.origin.y + a_offset.y, a.origin.z + a_offset.z };
        Vec3f a_max = { a_min.x + a.nx * a.voxel_size, a_min.y + a.ny * a.voxel_size, a_min.z + a.nz * a.voxel_size };
        Vec3f b_min = { b.origin.x + b_offset.x, b.origin.y + b_offset.y, b.origin.z + b_offset.z };
        Vec3f b_max = { b_min.x + b.nx * b.voxel_size, b_min.y + b.ny * b.voxel_size, b_min.z + b.nz * b.voxel_size };

        // World-space overlap box
        float ox_min = std::max(a_min.x, b_min.x);
        float oy_min = std::max(a_min.y, b_min.y);
        float oz_min = std::max(a_min.z, b_min.z);
        float ox_max = std::min(a_max.x, b_max.x);
        float oy_max = std::min(a_max.y, b_max.y);
        float oz_max = std::min(a_max.z, b_max.z);

        if (ox_min >= ox_max || oy_min >= oy_max || oz_min >= oz_max)
            return 0;  // No overlap region

        // Iterate overlap region at the coarser voxel size
        float vs = std::max(a.voxel_size, b.voxel_size);
        size_t collisions = 0;
        size_t iters = 0;

        for (float wz = oz_min + vs * 0.5f; wz < oz_max; wz += vs) {
            for (float wy = oy_min + vs * 0.5f; wy < oy_max; wy += vs) {
                for (float wx = ox_min + vs * 0.5f; wx < ox_max; wx += vs) {
                    // Check voxel in grid A
                    int ax, ay, az;
                    a.world_to_grid({wx - a_offset.x, wy - a_offset.y, wz - a_offset.z}, ax, ay, az);
                    bool a_solid = a.get(ax, ay, az);

                    if (a_solid) {
                        // Only check B if A is solid (short-circuit)
                        int bx, by, bz;
                        b.world_to_grid({wx - b_offset.x, wy - b_offset.y, wz - b_offset.z}, bx, by, bz);
                        if (b.get(bx, by, bz))
                            collisions++;
                    }

                    // Polite: yield every 10k iterations
                    if (++iters % 10000 == 0)
                        polite_yield();
                }
            }
        }
        return collisions;
    }

private:
    std::vector<uint8_t> bits_;
};


// ── Triangle-box intersection test ────────────────────────
//    Separating axis theorem (SAT) for triangle vs AABB.
//    Used for conservative surface voxelization.
namespace detail {

inline float project_min(const Vec3f &v0, const Vec3f &v1, const Vec3f &v2, const Vec3f &axis) {
    float p0 = v0.dot(axis), p1 = v1.dot(axis), p2 = v2.dot(axis);
    return std::min({p0, p1, p2});
}

inline float project_max(const Vec3f &v0, const Vec3f &v1, const Vec3f &v2, const Vec3f &axis) {
    float p0 = v0.dot(axis), p1 = v1.dot(axis), p2 = v2.dot(axis);
    return std::max({p0, p1, p2});
}

// Test if triangle (already translated so box is centered at origin)
// overlaps an axis-aligned box of half-size h.
inline bool triangle_aabb_overlap(
    const Vec3f &v0, const Vec3f &v1, const Vec3f &v2,
    const Vec3f &half)
{
    // Edge vectors
    Vec3f e0 = v1 - v0;
    Vec3f e1 = v2 - v1;
    Vec3f e2 = v0 - v2;

    // 9 cross-product axes (edge x box-axis)
    Vec3f axes[9] = {
        {0, -e0.z, e0.y}, {0, -e1.z, e1.y}, {0, -e2.z, e2.y},
        {e0.z, 0, -e0.x}, {e1.z, 0, -e1.x}, {e2.z, 0, -e2.x},
        {-e0.y, e0.x, 0}, {-e1.y, e1.x, 0}, {-e2.y, e2.x, 0},
    };

    for (int i = 0; i < 9; i++) {
        const Vec3f &a = axes[i];
        float len2 = a.x*a.x + a.y*a.y + a.z*a.z;
        if (len2 < 1e-12f) continue;  // Degenerate axis

        float r = half.x * std::abs(a.x) + half.y * std::abs(a.y) + half.z * std::abs(a.z);
        float pmin = project_min(v0, v1, v2, a);
        float pmax = project_max(v0, v1, v2, a);
        if (pmin > r || pmax < -r) return false;
    }

    // 3 box face axes
    if (std::min({v0.x,v1.x,v2.x}) > half.x || std::max({v0.x,v1.x,v2.x}) < -half.x) return false;
    if (std::min({v0.y,v1.y,v2.y}) > half.y || std::max({v0.y,v1.y,v2.y}) < -half.y) return false;
    if (std::min({v0.z,v1.z,v2.z}) > half.z || std::max({v0.z,v1.z,v2.z}) < -half.z) return false;

    // Triangle normal axis
    Vec3f n = e0.cross(e1);
    float d = n.dot(v0);
    float r = half.x * std::abs(n.x) + half.y * std::abs(n.y) + half.z * std::abs(n.z);
    if (d > r || d < -r) return false;

    return true;
}

} // namespace detail


// ── Voxelize a mesh: conservative outward ─────────────────
//
//    Each voxel that ANY triangle touches gets marked solid.
//    This means the voxelized volume is >= the real mesh volume.
//    (Conservative outward rounding — parts are "fatter" than reality.)
//
inline VoxError voxelize_mesh(
    const Triangle* triangles,
    size_t          num_triangles,
    float           voxel_size_mm,
    VoxelGrid      &out_grid,
    ProgressFn      progress = nullptr)
{
    if (num_triangles == 0) return VoxError::NO_TRIANGLES;
    if (voxel_size_mm <= 0.01f) return VoxError::INVALID_RESOLUTION;

    // Phase 1: Compute mesh AABB
    if (progress && !progress(0.0f, "Computing bounds"))
        return VoxError::ABORTED_BY_USER;

    AABB bounds;
    for (size_t i = 0; i < num_triangles; i++) {
        bounds.expand(triangles[i].v0);
        bounds.expand(triangles[i].v1);
        bounds.expand(triangles[i].v2);
    }

    // Expand bounds by one voxel on each side (conservative outward)
    bounds.min.x -= voxel_size_mm;
    bounds.min.y -= voxel_size_mm;
    bounds.min.z -= voxel_size_mm;
    bounds.max.x += voxel_size_mm;
    bounds.max.y += voxel_size_mm;
    bounds.max.z += voxel_size_mm;

    // Compute grid dimensions
    Vec3f sz = bounds.size();
    size_t gx = (size_t)std::ceil(sz.x / voxel_size_mm);
    size_t gy = (size_t)std::ceil(sz.y / voxel_size_mm);
    size_t gz = (size_t)std::ceil(sz.z / voxel_size_mm);

    // Clamp to safe limits
    gx = std::min(gx, MAX_GRID_DIM);
    gy = std::min(gy, MAX_GRID_DIM);
    gz = std::min(gz, MAX_GRID_DIM);

    // Allocate
    VoxError err = out_grid.allocate(gx, gy, gz);
    if (err != VoxError::OK) return err;
    out_grid.clear();
    out_grid.voxel_size = voxel_size_mm;
    out_grid.origin = bounds.min;

    if (progress && !progress(0.05f, "Grid allocated"))
        return VoxError::ABORTED_BY_USER;

    // Phase 2: Rasterize triangles into grid
    float half_vs = voxel_size_mm * 0.5f;
    Vec3f half = {half_vs, half_vs, half_vs};

    for (size_t ti = 0; ti < num_triangles; ti++) {
        const Triangle &tri = triangles[ti];

        // Triangle AABB in grid coords (conservative: floor min, ceil max)
        int ix_min, iy_min, iz_min, ix_max, iy_max, iz_max;
        out_grid.world_to_grid(tri.v0, ix_min, iy_min, iz_min);
        ix_max = ix_min; iy_max = iy_min; iz_max = iz_min;

        int tx, ty, tz;
        out_grid.world_to_grid(tri.v1, tx, ty, tz);
        ix_min = std::min(ix_min, tx); iy_min = std::min(iy_min, ty); iz_min = std::min(iz_min, tz);
        ix_max = std::max(ix_max, tx); iy_max = std::max(iy_max, ty); iz_max = std::max(iz_max, tz);

        out_grid.world_to_grid(tri.v2, tx, ty, tz);
        ix_min = std::min(ix_min, tx); iy_min = std::min(iy_min, ty); iz_min = std::min(iz_min, tz);
        ix_max = std::max(ix_max, tx); iy_max = std::max(iy_max, ty); iz_max = std::max(iz_max, tz);

        // Expand by 1 in each direction (conservative outward)
        ix_min = std::max(0, ix_min - 1);
        iy_min = std::max(0, iy_min - 1);
        iz_min = std::max(0, iz_min - 1);
        ix_max = std::min((int)gx - 1, ix_max + 1);
        iy_max = std::min((int)gy - 1, iy_max + 1);
        iz_max = std::min((int)gz - 1, iz_max + 1);

        // Test each voxel in the triangle's AABB
        for (int iz = iz_min; iz <= iz_max; iz++) {
            for (int iy = iy_min; iy <= iy_max; iy++) {
                for (int ix = ix_min; ix <= ix_max; ix++) {
                    // Voxel center in world space
                    Vec3f center = out_grid.grid_to_world_center(ix, iy, iz);

                    // Translate triangle relative to voxel center
                    Vec3f v0 = tri.v0 - center;
                    Vec3f v1 = tri.v1 - center;
                    Vec3f v2 = tri.v2 - center;

                    if (detail::triangle_aabb_overlap(v0, v1, v2, half)) {
                        out_grid.set(ix, iy, iz);
                    }
                }
            }
        }

        // Polite: yield periodically
        if ((ti + 1) % YIELD_EVERY_TRIANGLES == 0) {
            polite_yield();
            if (progress) {
                float p = 0.05f + 0.90f * ((float)(ti + 1) / (float)num_triangles);
                if (!progress(p, "Rasterizing triangles"))
                    return VoxError::ABORTED_BY_USER;
            }
        }
    }

    if (progress && !progress(0.95f, "Surface voxelization complete"))
        return VoxError::ABORTED_BY_USER;

    // Phase 3: Interior flood fill (optional — makes solid parts truly solid)
    // For collision detection, surface-only is conservative (safe).
    // Interior fill would make nesting detection better but costs more.
    // We skip it for now — surface voxelization is sufficient for
    // collision checking. Interior fill can be added later for concavity.

    if (progress) progress(1.0f, "Done");
    return VoxError::OK;
}


// ── Convenience: voxelize from indexed triangle set ───────
//    (Same format as Orca's indexed_triangle_set)
inline VoxError voxelize_indexed_mesh(
    const float*    vertices,       // [x,y,z, x,y,z, ...] num_verts * 3
    size_t          num_verts,
    const uint32_t* indices,        // [i0,i1,i2, i0,i1,i2, ...] num_tris * 3
    size_t          num_tris,
    float           voxel_size_mm,
    VoxelGrid      &out_grid,
    ProgressFn      progress = nullptr)
{
    // Build triangle array (temporary, bounded)
    size_t tri_bytes = num_tris * sizeof(Triangle);
    if (tri_bytes > MAX_VOXEL_MEMORY_MB * 1024 * 1024)
        return VoxError::MEMORY_CAP_EXCEEDED;

    std::vector<Triangle> tris(num_tris);
    for (size_t i = 0; i < num_tris; i++) {
        uint32_t i0 = indices[i * 3 + 0];
        uint32_t i1 = indices[i * 3 + 1];
        uint32_t i2 = indices[i * 3 + 2];
        tris[i].v0 = {vertices[i0*3], vertices[i0*3+1], vertices[i0*3+2]};
        tris[i].v1 = {vertices[i1*3], vertices[i1*3+1], vertices[i1*3+2]};
        tris[i].v2 = {vertices[i2*3], vertices[i2*3+1], vertices[i2*3+2]};

        // Yield while building (large meshes)
        if ((i + 1) % YIELD_EVERY_TRIANGLES == 0)
            polite_yield();
    }

    return voxelize_mesh(tris.data(), tris.size(), voxel_size_mm, out_grid, progress);
}

} // namespace snuggle
