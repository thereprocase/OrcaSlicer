// snuggle_test — Standalone test harness for the Snuggle radial nester
//
// Exercises polite_voxelizer, auto_snuggle, snuggle_radial, and optionally
// the GPU collision evaluator — without OrcaSlicer's GUI.
//
// Build: cmake -DSNUGGLE_TEST_GPU=ON for GPU tests
// Usage: snuggle_test [flags]
//
// Flags:
//   --filter <prefix>     Run only tests matching prefix (e.g., T03)
//   --stl <path>          Use real STL mesh for integration tests
//   --copies <N>          Number of copies for STL test (default 3)
//   --verbose             Print all assertions, not just failures
//   --no-stress           Skip stress tests (T09-T11)
//   --gpu                 Run GPU evaluator parity tests (needs SLIC3R_GUI)
//   --parts <N>           Custom part count for synthetic test
//   --shape <type>        Shape: cube, l-shape, bar (default cube)
//   --size <mm>           Part size in mm (default 20)
//   --voxel <mm>          Voxel size (default 2.0)
//   --gap <mm>            Gap between parts (default 1.5)
//   --bed <WxH>           Bed size (default 256x256)
//   --rotation <mode>     locked, 15, 45, 90 (default 15)
//   --benchmark           Run timing benchmarks
//   --repeat <N>          Repeat benchmark N times (default 5)
//   --csv <path>          Write results to CSV
//   --adversarial         Generate random adversarial geometries
//   --seed <N>            Random seed for adversarial mode (default 42)
//   --waves <N>           Number of adversarial waves (default 5)
//   --sweep-parts <list>  Comma-separated part counts to sweep
//   --sweep-voxel <list>  Comma-separated voxel sizes to sweep

#include "libslic3r/Arrange/polite_voxelizer.hpp"
#include "libslic3r/Arrange/auto_snuggle.hpp"
#include "libslic3r/Arrange/gpu_collision.hpp"
#include "libslic3r/Arrange/snuggle_nester.hpp"
#include "libslic3r/Arrange/snuggle_radial.hpp"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <string>
#include <vector>
#include <chrono>
#include <functional>
#include <algorithm>
#include <numeric>
#include <sstream>
#include <iostream>
#include <fstream>
#include <random>

using namespace snuggle;

// ═══════════════════════════════════════════════════════════
// Test framework — minimal, no external deps
// ═══════════════════════════════════════════════════════════

struct TestResult {
    std::string name;
    bool        passed;
    std::string message;
    double      elapsed_ms;
};

static std::vector<TestResult> g_results;
static bool g_verbose = false;

#define ASSERT_TRUE(cond, msg) \
    do { if (!(cond)) { \
        throw std::runtime_error(std::string("ASSERT_TRUE failed: ") + (msg) \
                                 + " (" #cond ")"); \
    }} while(0)

#define ASSERT_EQ(a, b, msg) \
    do { if ((a) != (b)) { \
        std::ostringstream _oss; \
        _oss << (msg) << " [got " << (a) << ", expected " << (b) << "]"; \
        throw std::runtime_error(_oss.str()); \
    }} while(0)

#define ASSERT_NEAR(a, b, tol, msg) \
    do { if (std::abs((float)(a) - (float)(b)) > (float)(tol)) { \
        std::ostringstream _oss; \
        _oss << (msg) << " [got " << (a) << ", expected ~" << (b) \
             << " +/-" << (tol) << "]"; \
        throw std::runtime_error(_oss.str()); \
    }} while(0)

static void run_test(const std::string& name,
                     const std::function<void()>& fn)
{
    auto t0 = std::chrono::steady_clock::now();
    TestResult r;
    r.name = name;
    try {
        fn();
        r.passed = true;
        r.message = "OK";
    } catch (const std::exception& e) {
        r.passed = false;
        r.message = e.what();
    } catch (...) {
        r.passed = false;
        r.message = "Unknown exception";
    }
    auto t1 = std::chrono::steady_clock::now();
    r.elapsed_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    g_results.push_back(r);

    printf("[%s] %-50s  %.1f ms\n",
           r.passed ? "PASS" : "FAIL", name.c_str(), r.elapsed_ms);
    if (!r.passed || g_verbose)
        printf("       %s\n", r.message.c_str());
    fflush(stdout);
}

// ═══════════════════════════════════════════════════════════
// Config — parsed from flags
// ═══════════════════════════════════════════════════════════

struct HarnessConfig {
    std::string filter;
    std::string stl_path;
    int         copies      = 3;
    bool        no_stress   = false;
    bool        gpu         = false;
    int         parts       = 0;      // 0 = use test default
    std::string shape       = "cube";
    float       size_mm     = 20.0f;
    float       voxel_mm    = 2.0f;
    float       gap_mm      = 1.5f;
    float       bed_w       = 256.0f;
    float       bed_h       = 256.0f;
    int         rotation    = 15;     // degrees, 0 = locked
    bool        benchmark   = false;
    int         repeat      = 5;
    std::string csv_path;
    bool        adversarial = false;
    int         seed        = 42;
    int         waves       = 5;
    std::vector<int>   sweep_parts;
    std::vector<float> sweep_voxel;
};

static std::vector<int> parse_int_list(const std::string& s) {
    std::vector<int> v;
    std::istringstream ss(s);
    std::string tok;
    while (std::getline(ss, tok, ','))
        v.push_back(std::stoi(tok));
    return v;
}

static std::vector<float> parse_float_list(const std::string& s) {
    std::vector<float> v;
    std::istringstream ss(s);
    std::string tok;
    while (std::getline(ss, tok, ','))
        v.push_back(std::stof(tok));
    return v;
}

static HarnessConfig parse_args(int argc, char** argv) {
    HarnessConfig c;
    for (int i = 1; i < argc; i++) {
        std::string a = argv[i];
        auto next = [&]() -> std::string {
            return (i + 1 < argc) ? argv[++i] : "";
        };
        if (a == "--verbose")       g_verbose = true;
        else if (a == "--filter")   c.filter = next();
        else if (a == "--stl")      c.stl_path = next();
        else if (a == "--copies")   c.copies = std::stoi(next());
        else if (a == "--no-stress") c.no_stress = true;
        else if (a == "--gpu")      c.gpu = true;
        else if (a == "--parts")    c.parts = std::stoi(next());
        else if (a == "--shape")    c.shape = next();
        else if (a == "--size")     c.size_mm = std::stof(next());
        else if (a == "--voxel")    c.voxel_mm = std::stof(next());
        else if (a == "--gap")      c.gap_mm = std::stof(next());
        else if (a == "--bed") {
            std::string b = next();
            auto x = b.find('x');
            if (x != std::string::npos) {
                c.bed_w = std::stof(b.substr(0, x));
                c.bed_h = std::stof(b.substr(x + 1));
            }
        }
        else if (a == "--rotation") c.rotation = std::stoi(next());
        else if (a == "--benchmark") c.benchmark = true;
        else if (a == "--repeat")   c.repeat = std::stoi(next());
        else if (a == "--csv")      c.csv_path = next();
        else if (a == "--adversarial") c.adversarial = true;
        else if (a == "--seed")     c.seed = std::stoi(next());
        else if (a == "--waves")    c.waves = std::stoi(next());
        else if (a == "--sweep-parts") c.sweep_parts = parse_int_list(next());
        else if (a == "--sweep-voxel") c.sweep_voxel = parse_float_list(next());
        else fprintf(stderr, "Unknown flag: %s\n", a.c_str());
    }
    return c;
}

// ═══════════════════════════════════════════════════════════
// Mesh factories
// ═══════════════════════════════════════════════════════════

struct Tri { float v0[3], v1[3], v2[3]; };

static std::vector<Tri> make_box(float w, float h, float d) {
    float V[][3] = {
        {0,0,0},{w,0,0},{w,h,0},{0,h,0},
        {0,0,d},{w,0,d},{w,h,d},{0,h,d}
    };
    auto t = [&](int a, int b, int c) -> Tri {
        return {{V[a][0],V[a][1],V[a][2]},
                {V[b][0],V[b][1],V[b][2]},
                {V[c][0],V[c][1],V[c][2]}};
    };
    return {
        t(0,2,1), t(0,3,2), // bottom
        t(4,5,6), t(4,6,7), // top
        t(0,1,5), t(0,5,4), // front
        t(1,2,6), t(1,6,5), // right
        t(2,3,7), t(2,7,6), // back
        t(3,0,4), t(3,4,7)  // left
    };
}

static std::vector<Tri> make_L(float size, float height) {
    auto a = make_box(size/2, size, height);
    auto b = make_box(size/2, size/2, height);
    for (auto& t : b) {
        t.v0[0] += size/2; t.v1[0] += size/2; t.v2[0] += size/2;
    }
    a.insert(a.end(), b.begin(), b.end());
    return a;
}

static std::vector<Tri> make_bar(float length, float width, float height) {
    return make_box(length, width, height);
}

static std::vector<Tri> make_shape(const std::string& shape, float size) {
    if (shape == "cube")    return make_box(size, size, size);
    if (shape == "l-shape") return make_L(size, size / 2);
    if (shape == "bar")     return make_bar(size, size / 4, size / 4);
    return make_box(size, size, size);
}

static std::vector<Tri> load_stl(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("Cannot open STL: " + path);
    char header[80]; f.read(header, 80);
    uint32_t n = 0; f.read(reinterpret_cast<char*>(&n), 4);
    std::vector<Tri> tris(n);
    for (uint32_t i = 0; i < n; i++) {
        float buf[12]; uint16_t attr;
        f.read(reinterpret_cast<char*>(buf), 48);
        f.read(reinterpret_cast<char*>(&attr), 2);
        tris[i] = {{buf[3],buf[4],buf[5]},{buf[6],buf[7],buf[8]},{buf[9],buf[10],buf[11]}};
    }
    return tris;
}

// ═══════════════════════════════════════════════════════════
// Helpers
// ═══════════════════════════════════════════════════════════

static VoxelGrid voxelize(const std::vector<Tri>& mesh, float voxel_mm) {
    VoxelGrid g;
    std::vector<float> verts;
    std::vector<uint32_t> indices;
    for (size_t i = 0; i < mesh.size(); i++) {
        size_t base = verts.size() / 3;
        for (int v = 0; v < 3; v++) {
            const float* p = (v == 0) ? mesh[i].v0 : (v == 1) ? mesh[i].v1 : mesh[i].v2;
            verts.push_back(p[0]); verts.push_back(p[1]); verts.push_back(p[2]);
        }
        indices.push_back((uint32_t)(base));
        indices.push_back((uint32_t)(base + 1));
        indices.push_back((uint32_t)(base + 2));
    }
    auto err = voxelize_indexed_mesh(verts.data(), verts.size() / 3,
                                      indices.data(), mesh.size(),
                                      voxel_mm, g);
    if (err != VoxError::OK)
        throw std::runtime_error("Voxelization failed: " + std::to_string((int)err));
    return g;
}

static PartInfo make_part(const VoxelGrid& grid, const std::string& name,
                           float height, float zrot = 0.0f) {
    PartInfo p;
    p.grid = grid;
    p.name = name;
    p.max_height_mm = height;
    p.hull_area_mm2 = (float)(grid.nx * grid.ny) * grid.voxel_size * grid.voxel_size;
    p.initial_zrot = zrot;
    return p;
}

static std::vector<VoxelGrid> build_cache(const VoxelGrid& grid, int bins) {
    std::vector<VoxelGrid> c(bins);
    for (int i = 0; i < bins; i++)
        c[i] = grid.rotated_copy((float)i * TWO_PI_F / bins);
    return c;
}

static void validate_no_collisions(const RadialResult& result,
    const std::vector<std::vector<VoxelGrid>>& caches)
{
    for (size_t i = 0; i < result.placements.size(); i++) {
        if (!result.placements[i].placed) continue;
        for (size_t j = i + 1; j < result.placements.size(); j++) {
            if (!result.placements[j].placed) continue;
            const auto& pi = result.placements[i];
            const auto& pj = result.placements[j];
            int ni = (int)caches[i].size(), nj = (int)caches[j].size();
            int bi = ((int)std::floor(pi.zrot * ni / TWO_PI_F) % ni + ni) % ni;
            int bj = ((int)std::floor(pj.zrot * nj / TWO_PI_F) % nj + nj) % nj;
            Vec3f oi = {pi.x, pi.y, 0}, oj = {pj.x, pj.y, 0};
            size_t hits = VoxelGrid::collision_count(caches[i][bi], oi, caches[j][bj], oj);
            if (hits > 0) {
                std::ostringstream s;
                s << "Collision: part " << i << " and " << j << " (" << hits << " voxels)";
                throw std::runtime_error(s.str());
            }
        }
    }
}

static void validate_bed_bounds(const RadialResult& result,
    const std::vector<std::vector<VoxelGrid>>& caches,
    float bed_w, float bed_h, float margin)
{
    for (size_t i = 0; i < result.placements.size(); i++) {
        if (!result.placements[i].placed) continue;
        const auto& p = result.placements[i];
        int n = (int)caches[i].size();
        int b = ((int)std::floor(p.zrot * n / TWO_PI_F) % n + n) % n;
        const auto& g = caches[i][b];
        float x0 = p.x + g.origin.x, y0 = p.y + g.origin.y;
        float x1 = x0 + g.nx * g.voxel_size, y1 = y0 + g.ny * g.voxel_size;
        if (x0 < margin || y0 < margin || x1 > bed_w - margin || y1 > bed_h - margin) {
            std::ostringstream s;
            s << "Part " << i << " OOB: [" << x0 << "," << y0 << "]-["
              << x1 << "," << y1 << "] on " << bed_w << "x" << bed_h;
            throw std::runtime_error(s.str());
        }
    }
}

// ═══════════════════════════════════════════════════════════
// Mock evaluator for batch index test
// ═══════════════════════════════════════════════════════════

class MockEvaluator : public CollisionEvaluator {
public:
    std::vector<RadialCollisionResult> canned;
    bool oob = false;

    void upload_grids(const std::vector<PartInfo>&,
                      const std::vector<std::vector<VoxelGrid>>&) override {}
    void evaluate_batch(std::vector<Individual>&,
                        const std::vector<PartInfo>&,
                        float, float, float) override {}
    void evaluate_radial(
        const std::vector<size_t>&, const std::vector<RadialCandidate>&,
        size_t, const std::vector<RadialCandidate>& candidates,
        float, float, float, float,
        std::vector<RadialCollisionResult>& results) override
    {
        results.resize(candidates.size());
        for (size_t i = 0; i < candidates.size(); i++) {
            if (i < canned.size()) results[i] = canned[i];
            else { oob = true; results[i] = {true}; }
        }
    }
};

// ═══════════════════════════════════════════════════════════
// Tests
// ═══════════════════════════════════════════════════════════

static void t01_auto_config() {
    // Zero parts
    auto c0 = compute_auto_config(0, {}, 256, 256, 1.5f, false, 15, false);
    ASSERT_TRUE(!std::isnan(c0.voxel_mm), "0 parts: no NaN");

    // 1 part, 1mm
    auto c1 = compute_auto_config(1, {1.0f}, 256, 256, 1.5f, false, 15, false);
    ASSERT_NEAR(c1.voxel_mm, 0.5f, 0.01f, "1mm: voxel clamp 0.5");
    ASSERT_EQ(c1.n_directions, 36, "1 part: 36 dirs");

    // 50 parts
    std::vector<float> d50(50, 30.0f);
    auto c50 = compute_auto_config(50, d50, 256, 256, 1.5f, false, 15, false);
    ASSERT_EQ(c50.n_directions, 12, "50 parts: 12 dirs");
    ASSERT_TRUE(c50.n_rotations <= 8, "50 parts: rots<=8");

    // Lock rotation
    auto cl = compute_auto_config(5, {30.f}, 256, 256, 1.5f, true, 15, false);
    ASSERT_EQ(cl.n_rotations, 1, "locked: 1 rot");

    // Vase trap
    auto cv = compute_auto_config(10, {200.f,15.f,15.f,15.f,15.f,15.f,15.f,15.f,15.f,15.f},
                                   256, 256, 1.5f, false, 15, false);
    ASSERT_TRUE(cv.voxel_mm >= 200.f / 120.f, "vase trap: min voxel");

    // GPU bonus (only >5 parts, needs voxel above clamp floor to see effect)
    std::vector<float> d10(10, 80.f);  // 80mm parts → base_voxel = 1.33, above 0.5 clamp
    auto cc = compute_auto_config(10, d10, 256, 256, 1.5f, false, 15, false);
    auto cg = compute_auto_config(10, d10, 256, 256, 1.5f, false, 15, true);
    ASSERT_TRUE(cg.voxel_mm < cc.voxel_mm, "GPU: finer voxel");
}

static void t02_voxelgrid() {
    VoxelGrid g;
    ASSERT_TRUE(g.allocate(4,4,4) == VoxError::OK, "alloc 4x4x4");
    g.set(0,0,0); g.set(3,3,3);
    ASSERT_TRUE(g.get(0,0,0), "get set voxel");
    ASSERT_TRUE(!g.get(1,0,0), "get unset voxel");

    ASSERT_TRUE(VoxelGrid().allocate(513,1,1) == VoxError::GRID_TOO_LARGE, "oversize reject");

    // Collision: same grid overlaps itself
    VoxelGrid s; s.allocate(4,4,4); s.voxel_size = 2.0f; s.origin = {0,0,0};
    for (int z=0;z<4;z++) for (int y=0;y<4;y++) for (int x=0;x<4;x++) s.set(x,y,z);
    ASSERT_EQ(s.count_solid(), (size_t)64, "count_solid");
    Vec3f o0={0,0,0}, o1={100,0,0};
    ASSERT_TRUE(VoxelGrid::collision_count(s,o0,s,o0) > 0, "self-overlap");
    ASSERT_EQ(VoxelGrid::collision_count(s,o0,s,o1), (size_t)0, "separated");

    // Rotation preserves solid count
    VoxelGrid bar; bar.allocate(3,1,1); bar.voxel_size = 2.0f; bar.origin = {0,0,0};
    bar.set(0,0,0); bar.set(1,0,0); bar.set(2,0,0);
    auto rot = bar.rotated_copy(PI_F / 2.0f);
    ASSERT_EQ(rot.count_solid(), bar.count_solid(), "rotation preserves count");
}

static void t03_radial_synthetic(const HarnessConfig& cfg) {
    float sz = cfg.parts > 0 ? cfg.size_mm : 10.0f;
    int n = cfg.parts > 0 ? cfg.parts : 2;
    auto mesh = make_shape(cfg.shape, sz);
    auto grid = voxelize(mesh, cfg.voxel_mm);

    std::vector<PartInfo> parts;
    std::vector<std::vector<VoxelGrid>> caches(n);
    for (int i = 0; i < n; i++) {
        parts.push_back(make_part(grid, "p" + std::to_string(i), sz));
        caches[i] = build_cache(grid, 24);
    }

    RadialConfig rc;
    rc.bed_width_mm = cfg.bed_w; rc.bed_height_mm = cfg.bed_h;
    rc.min_gap_mm = cfg.gap_mm; rc.bed_margin_mm = cfg.gap_mm;
    rc.n_directions = 24; rc.n_rotations = 24;
    rc.step_mm = cfg.voxel_mm; rc.timeout_s = 30.0;
    rc.lock_rotation = (cfg.rotation == 0);

    auto result = radial_arrange(parts, caches, rc);
    printf("  Placed %d/%d in %.1fms\n", result.placed_count, n, result.time_ms);
    ASSERT_EQ(result.placed_count, n, "all parts placed");
    validate_no_collisions(result, caches);
    validate_bed_bounds(result, caches, cfg.bed_w, cfg.bed_h, cfg.gap_mm);
}

static void t04_batch_index() {
    VoxelGrid g; g.allocate(1,1,1); g.voxel_size = 2.0f; g.origin = {0,0,0}; g.set(0,0,0);
    PartInfo p0 = make_part(g, "placed", 2.f);
    PartInfo p1 = make_part(g, "cand", 2.f);

    int N_DIRS = 3, N_ROTS = 4;
    std::vector<std::vector<VoxelGrid>> caches(2);
    caches[0] = std::vector<VoxelGrid>(N_ROTS, g);
    caches[1] = std::vector<VoxelGrid>(N_ROTS, g);

    RadialConfig rc;
    rc.bed_width_mm = 256; rc.bed_height_mm = 256;
    rc.min_gap_mm = 0; rc.bed_margin_mm = 0;
    rc.n_directions = N_DIRS; rc.n_rotations = N_ROTS;
    rc.step_mm = 100.0f; rc.timeout_s = 30.0;

    MockEvaluator mock;
    mock.canned.resize(N_DIRS * N_ROTS, {true});
    mock.canned[6] = {false}; // dir 1, rot 2

    auto result = radial_arrange({p0, p1}, caches, rc, &mock);
    ASSERT_TRUE(!mock.oob, "no OOB access in batch scan");
}

static void t05_degradation() {
    auto mesh = make_box(20, 20, 5);
    auto grid = voxelize(mesh, 2.0f);
    int N = 10;
    std::vector<PartInfo> parts;
    std::vector<std::vector<VoxelGrid>> caches(N);
    for (int i = 0; i < N; i++) {
        parts.push_back(make_part(grid, "p" + std::to_string(i), 5.f, (float)i * PI_F / 6));
        caches[i] = build_cache(grid, 24);
    }

    // Near-zero timeout forces degradation/timeout
    RadialConfig rc;
    rc.bed_width_mm = 256; rc.bed_height_mm = 256;
    rc.min_gap_mm = 1; rc.bed_margin_mm = 3;
    rc.n_directions = 24; rc.n_rotations = 24;
    rc.step_mm = 2.0f; rc.timeout_s = 0.0001;

    auto result = radial_arrange(parts, caches, rc);
    ASSERT_TRUE(result.timed_out || result.placed_count < N,
                "instant timeout: timed out or partial placement");

    // Lock rotation
    rc.timeout_s = 30.0;
    rc.lock_rotation = true; rc.n_rotations = 1;
    auto locked = radial_arrange(parts, caches, rc);
    for (int i = 0; i < N; i++) {
        if (!locked.placements[i].placed) continue;
        ASSERT_NEAR(locked.placements[i].zrot, parts[i].initial_zrot, 1e-5f,
                    "locked: zrot matches initial");
    }
}

static void t06_integration(const HarnessConfig& cfg) {
    std::vector<Tri> mesh;
    if (!cfg.stl_path.empty()) mesh = load_stl(cfg.stl_path);
    else mesh = make_L(60, 60); // synthetic benchy substitute

    auto t0 = std::chrono::steady_clock::now();
    auto grid = voxelize(mesh, cfg.voxel_mm);
    auto t1 = std::chrono::steady_clock::now();
    double vox_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    int N = cfg.copies;
    std::vector<PartInfo> parts;
    std::vector<std::vector<VoxelGrid>> caches(N);
    auto t2 = std::chrono::steady_clock::now();
    for (int i = 0; i < N; i++) {
        parts.push_back(make_part(grid, "part_" + std::to_string(i), 60.f));
        caches[i] = build_cache(grid, 24);
    }
    auto t3 = std::chrono::steady_clock::now();
    double cache_ms = std::chrono::duration<double, std::milli>(t3 - t2).count();

    RadialConfig rc;
    rc.bed_width_mm = cfg.bed_w; rc.bed_height_mm = cfg.bed_h;
    rc.min_gap_mm = cfg.gap_mm; rc.bed_margin_mm = cfg.gap_mm;
    rc.n_directions = 24; rc.n_rotations = 24;
    rc.step_mm = cfg.voxel_mm; rc.timeout_s = 60.0;

    auto result = radial_arrange(parts, caches, rc);
    printf("  vox=%.1fms cache=%.1fms place=%.1fms\n", vox_ms, cache_ms, result.time_ms);

    ASSERT_EQ(result.placed_count, N, "all placed");
    validate_no_collisions(result, caches);
    validate_bed_bounds(result, caches, cfg.bed_w, cfg.bed_h, cfg.gap_mm);

    // Determinism check
    auto result2 = radial_arrange(parts, caches, rc);
    for (int i = 0; i < N; i++) {
        if (!result.placements[i].placed) continue;
        ASSERT_EQ(result.placements[i].x, result2.placements[i].x, "deterministic x");
        ASSERT_EQ(result.placements[i].y, result2.placements[i].y, "deterministic y");
        ASSERT_EQ(result.placements[i].zrot, result2.placements[i].zrot, "deterministic zrot");
    }
}

static void t09_50_cubes(const HarnessConfig& cfg) {
    auto mesh = make_box(20, 20, 20);
    auto grid = voxelize(mesh, cfg.voxel_mm);
    int N = 50;
    std::vector<PartInfo> parts;
    std::vector<std::vector<VoxelGrid>> caches(N);
    for (int i = 0; i < N; i++) {
        parts.push_back(make_part(grid, "cube_" + std::to_string(i), 20.f));
        caches[i] = build_cache(grid, 8);
    }

    RadialConfig rc;
    rc.bed_width_mm = cfg.bed_w; rc.bed_height_mm = cfg.bed_h;
    rc.min_gap_mm = 1; rc.bed_margin_mm = 3;
    rc.n_directions = 12; rc.n_rotations = 8;
    rc.step_mm = cfg.voxel_mm; rc.timeout_s = 60.0;

    auto result = radial_arrange(parts, caches, rc);
    printf("  50 cubes: %d/%d in %.1fms\n", result.placed_count, N, result.time_ms);
    ASSERT_EQ(result.placed_count, N, "all 50 placed");
    validate_no_collisions(result, caches);
}

static void t10_vase_trap() {
    auto big = make_box(200, 200, 200);
    auto small_mesh = make_box(15, 15, 5);
    std::vector<float> dims = {200.f};
    for (int i = 0; i < 15; i++) dims.push_back(15.f);

    auto ac = compute_auto_config(16, dims, 256, 256, 1.5f, false, 15, false);
    ASSERT_TRUE(ac.voxel_mm >= 200.f / 120.f, "vase trap fires");

    auto big_grid = voxelize(big, ac.voxel_mm);
    ASSERT_TRUE(big_grid.nx <= 512, "grid fits in MAX_GRID_DIM");
}

// ═══════════════════════════════════════════════════════════
// Adversarial mode
// ═══════════════════════════════════════════════════════════

static void adversarial_wave(int wave, std::mt19937& rng, const HarnessConfig& cfg) {
    std::uniform_int_distribution<int> d_parts(1, 30);
    std::uniform_real_distribution<float> d_size(1.0f, 150.0f);
    std::uniform_real_distribution<float> d_voxel(0.5f, 4.0f);
    std::uniform_real_distribution<float> d_rot(0.0f, TWO_PI_F);

    int n = d_parts(rng);
    float vx = d_voxel(rng);

    std::vector<PartInfo> parts;
    std::vector<std::vector<VoxelGrid>> caches(n);

    for (int i = 0; i < n; i++) {
        float sz = d_size(rng);
        auto mesh = make_box(sz, sz * 0.5f + 0.1f, sz * 0.3f + 0.1f);
        try {
            auto grid = voxelize(mesh, vx);
            parts.push_back(make_part(grid, "adv_" + std::to_string(i),
                                       sz * 0.3f, d_rot(rng)));
            caches[i] = build_cache(grid, 12);
        } catch (...) {
            // Voxelization can fail on extreme sizes — skip this part
            VoxelGrid empty;
            empty.allocate(1,1,1); empty.voxel_size = vx; empty.origin = {0,0,0};
            parts.push_back(make_part(empty, "adv_fail_" + std::to_string(i), 1.0f));
            caches[i] = {empty};
        }
    }

    RadialConfig rc;
    rc.bed_width_mm = cfg.bed_w; rc.bed_height_mm = cfg.bed_h;
    rc.min_gap_mm = 1.5f; rc.bed_margin_mm = 3.0f;
    rc.n_directions = 16; rc.n_rotations = 12;
    rc.step_mm = vx; rc.timeout_s = 15.0;

    auto result = radial_arrange(parts, caches, rc);

    // Validate: placed parts must not collide and must be in bounds
    validate_no_collisions(result, caches);
    validate_bed_bounds(result, caches, cfg.bed_w, cfg.bed_h, rc.bed_margin_mm);
}

// ═══════════════════════════════════════════════════════════
// Sweep mode
// ═══════════════════════════════════════════════════════════

static void run_sweep(const HarnessConfig& cfg) {
    auto parts_list = cfg.sweep_parts.empty() ? std::vector<int>{2,5,10,20} : cfg.sweep_parts;
    auto voxel_list = cfg.sweep_voxel.empty() ? std::vector<float>{1.0f,2.0f,4.0f} : cfg.sweep_voxel;

    std::ofstream csv;
    if (!cfg.csv_path.empty()) {
        csv.open(cfg.csv_path);
        csv << "parts,voxel_mm,placed,time_ms,collisions\n";
    }

    for (int n : parts_list) {
        for (float vx : voxel_list) {
            auto mesh = make_box(20, 20, 20);
            VoxelGrid grid;
            try { grid = voxelize(mesh, vx); } catch (...) { continue; }

            std::vector<PartInfo> parts;
            std::vector<std::vector<VoxelGrid>> caches(n);
            for (int i = 0; i < n; i++) {
                parts.push_back(make_part(grid, "s" + std::to_string(i), 20.f));
                caches[i] = build_cache(grid, 12);
            }

            RadialConfig rc;
            rc.bed_width_mm = cfg.bed_w; rc.bed_height_mm = cfg.bed_h;
            rc.min_gap_mm = 1.5f; rc.bed_margin_mm = 3.0f;
            rc.n_directions = 16; rc.n_rotations = 12;
            rc.step_mm = vx; rc.timeout_s = 30.0;

            auto result = radial_arrange(parts, caches, rc);

            bool clean = true;
            try { validate_no_collisions(result, caches); }
            catch (...) { clean = false; }

            printf("  N=%2d vx=%.1f: placed %d/%d in %.1fms %s\n",
                   n, vx, result.placed_count, n, result.time_ms,
                   clean ? "CLEAN" : "COLLISION");

            if (csv.is_open())
                csv << n << "," << vx << "," << result.placed_count
                    << "," << result.time_ms << "," << (clean ? 0 : 1) << "\n";
        }
    }
}

// ═══════════════════════════════════════════════════════════
// Main
// ═══════════════════════════════════════════════════════════

int main(int argc, char** argv) {
    auto cfg = parse_args(argc, argv);

    printf("Snuggle Radial Nester — Test Harness\n");
    printf("=====================================\n");
    if (!cfg.stl_path.empty()) printf("STL: %s\n", cfg.stl_path.c_str());
    if (cfg.adversarial) printf("Adversarial: seed=%d waves=%d\n", cfg.seed, cfg.waves);
    if (!cfg.sweep_parts.empty()) printf("Sweep mode\n");
    printf("\n");

    auto maybe = [&](const std::string& name, const std::function<void()>& fn) {
        if (cfg.filter.empty() || name.find(cfg.filter) != std::string::npos)
            run_test(name, fn);
    };

    // Unit tests
    maybe("T01_auto_config",  t01_auto_config);
    maybe("T02_voxelgrid",    t02_voxelgrid);
    maybe("T03_radial",       [&]{ t03_radial_synthetic(cfg); });
    maybe("T04_batch_index",  t04_batch_index);
    maybe("T05_degradation",  t05_degradation);

    // Integration
    maybe("T06_integration",  [&]{ t06_integration(cfg); });

    // Stress
    if (!cfg.no_stress) {
        maybe("T09_50_cubes", [&]{ t09_50_cubes(cfg); });
        maybe("T10_vase_trap", t10_vase_trap);
    }

    // Adversarial
    if (cfg.adversarial) {
        std::mt19937 rng(cfg.seed);
        for (int w = 0; w < cfg.waves; w++) {
            maybe("ADV_wave_" + std::to_string(w),
                  [&]{ adversarial_wave(w, rng, cfg); });
        }
    }

    // Sweep
    if (!cfg.sweep_parts.empty() || !cfg.sweep_voxel.empty()) {
        maybe("SWEEP", [&]{ run_sweep(cfg); });
    }

    // Benchmark
    if (cfg.benchmark) {
        maybe("BENCH", [&]{
            auto mesh = make_shape(cfg.shape, cfg.size_mm);
            auto grid = voxelize(mesh, cfg.voxel_mm);
            int n = cfg.parts > 0 ? cfg.parts : 10;
            std::vector<PartInfo> parts;
            std::vector<std::vector<VoxelGrid>> caches(n);
            for (int i = 0; i < n; i++) {
                parts.push_back(make_part(grid, "b" + std::to_string(i), cfg.size_mm));
                caches[i] = build_cache(grid, 24);
            }
            RadialConfig rc;
            rc.bed_width_mm = cfg.bed_w; rc.bed_height_mm = cfg.bed_h;
            rc.min_gap_mm = cfg.gap_mm; rc.bed_margin_mm = cfg.gap_mm;
            rc.n_directions = 24; rc.n_rotations = 24;
            rc.step_mm = cfg.voxel_mm; rc.timeout_s = 60.0;

            double total = 0;
            for (int r = 0; r < cfg.repeat; r++) {
                auto res = radial_arrange(parts, caches, rc);
                total += res.time_ms;
            }
            printf("  %d reps, avg %.1fms\n", cfg.repeat, total / cfg.repeat);
        });
    }

    // Summary
    printf("\n=====================================\n");
    int passed = 0, failed = 0;
    for (auto& r : g_results) r.passed ? passed++ : failed++;
    printf("Results: %d passed, %d failed\n", passed, failed);
    if (failed > 0) {
        printf("\nFailed:\n");
        for (auto& r : g_results)
            if (!r.passed) printf("  %s: %s\n", r.name.c_str(), r.message.c_str());
    }

    // CSV summary
    if (!cfg.csv_path.empty() && !cfg.sweep_parts.empty()) {
        printf("CSV written to %s\n", cfg.csv_path.c_str());
    }

    return failed > 0 ? 1 : 0;
}
