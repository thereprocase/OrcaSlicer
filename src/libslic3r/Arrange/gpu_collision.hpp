// gpu_collision.hpp -- Collision evaluator interface + GPU/CPU implementations
//
// Abstracts collision + bounds evaluation so the nester can delegate to
// GPU, SIMD, or any other accelerated backend. The nester owns fitness
// scoring (cheap, stays on CPU); this interface owns the expensive
// pairwise overlap and bounds checks.
//
// Lifetime: caller creates the evaluator, passes a non-owning pointer
// to SnuggleNester, and must keep it alive for the duration of run().

#pragma once

#include "polite_voxelizer.hpp"
#include <vector>
#include <memory>
#include <atomic>
#include <string>
#include <boost/log/trivial.hpp>

namespace snuggle {

// Forward declarations (defined in snuggle_nester.hpp)
struct PartInfo;
struct Individual;

// Radial nester types — used by evaluate_radial()
struct RadialCandidate {
    float x, y, zrot;
};

struct RadialCollisionResult {
    bool collides;
};

// Abstract interface -- nester calls this, doesn't know if CPU or GPU
class CollisionEvaluator {
public:
    virtual ~CollisionEvaluator() = default;

    // Upload pre-rotated grids so the evaluator can cache them.
    // rot_cache[part_idx][angle_bin] = rotated VoxelGrid.
    // Called once after the rotation cache is built.
    virtual void upload_grids(
        const std::vector<PartInfo>& parts,
        const std::vector<std::vector<VoxelGrid>>& rot_cache) = 0;

    // GA path: evaluate collision_count and oob_count for every individual.
    virtual void evaluate_batch(
        std::vector<Individual>& pop,
        const std::vector<PartInfo>& parts,
        float bed_width_mm,
        float bed_height_mm,
        float bed_margin_mm = 0.0f) = 0;

    // Radial path: test one part at multiple candidate positions against
    // already-placed parts. Gap enforcement uses 9-probe pattern.
    virtual void evaluate_radial(
        const std::vector<size_t>&            placed_parts,
        const std::vector<RadialCandidate>&   placed_positions,
        size_t                                candidate_part,
        const std::vector<RadialCandidate>&   candidates,
        float                                 gap_mm,
        float                                 bed_w,
        float                                 bed_h,
        float                                 bed_margin,
        std::vector<RadialCollisionResult>&   results) = 0;
};

// CPU fallback -- extracts the collision/bounds logic from SnuggleNester::evaluate()
class CpuCollisionEvaluator : public CollisionEvaluator {
public:
    void upload_grids(const std::vector<PartInfo>& parts,
                      const std::vector<std::vector<VoxelGrid>>& rot_cache) override;
    void evaluate_batch(std::vector<Individual>& individuals,
                        const std::vector<PartInfo>& parts,
                        float bed_w, float bed_h,
                        float bed_margin = 0.0f) override;
    void evaluate_radial(
        const std::vector<size_t>&            placed_parts,
        const std::vector<RadialCandidate>&   placed_positions,
        size_t                                candidate_part,
        const std::vector<RadialCandidate>&   candidates,
        float gap_mm, float bed_w, float bed_h, float bed_margin,
        std::vector<RadialCollisionResult>&   results) override;
private:
    const std::vector<std::vector<VoxelGrid>>* rot_cache_ = nullptr;
    static constexpr int ROT_CACHE_BINS = 360;

    const VoxelGrid& get_rotated(size_t part_idx, float angle) const;
};

#ifdef SLIC3R_GUI

// GPU compute shader backend (OpenGL 4.3 compute shaders)
class GpuCollisionEvaluator : public CollisionEvaluator {
public:
    GpuCollisionEvaluator();
    ~GpuCollisionEvaluator() override;

    bool is_available() const { return available_; }
    bool is_worthwhile() const { return available_ && gpu_cpu_ratio_ > 0.0f && gpu_cpu_ratio_ < 1.0f; }
    float gpu_cpu_ratio() const { return gpu_cpu_ratio_; }
    const std::string& renderer() const { return renderer_name_; }

    // Set cached probe result (skip re-benchmark if renderer matches)
    void set_cached_probe(const std::string& cached_renderer, float cached_ratio) {
        if (cached_renderer == renderer_name_ && cached_ratio > 0.0f) {
            gpu_cpu_ratio_ = cached_ratio;
            BOOST_LOG_TRIVIAL(info) << "Snuggle GPU: using cached probe ratio="
                << cached_ratio << " for " << renderer_name_;
        }
    }

    void upload_grids(const std::vector<PartInfo>& parts,
                      const std::vector<std::vector<VoxelGrid>>& rot_cache) override;
    void evaluate_batch(std::vector<Individual>& individuals,
                        const std::vector<PartInfo>& parts,
                        float bed_w, float bed_h,
                        float bed_margin = 0.0f) override;
    void evaluate_radial(
        const std::vector<size_t>&            placed_parts,
        const std::vector<RadialCandidate>&   placed_positions,
        size_t                                candidate_part,
        const std::vector<RadialCandidate>&   candidates,
        float gap_mm, float bed_w, float bed_h, float bed_margin,
        std::vector<RadialCollisionResult>&   results) override;

private:
    std::atomic<bool> available_{false};
    // GL handles — GA shader
    unsigned int program_ = 0;
    unsigned int voxel_ssbo_ = 0;
    unsigned int meta_ssbo_ = 0;
    unsigned int placement_ssbo_ = 0;
    unsigned int results_ssbo_ = 0;
    // GL handles — radial shader (shares voxel_ssbo_ and meta_ssbo_)
    unsigned int radial_program_ = 0;
    unsigned int radial_candidates_ssbo_ = 0;
    unsigned int radial_placed_ssbo_ = 0;
    unsigned int radial_results_ssbo_ = 0;

    // CPU fallback for evaluate_radial when shader isn't ready
    const std::vector<std::vector<VoxelGrid>>* rot_cache_ptr_ = nullptr;
    int actual_rot_bins_ = 0;  // actual cache size per part (may differ from ROT_BINS)
    float gpu_cpu_ratio_ = 0.0f;  // <1 means GPU is faster. 0 = not yet probed.
    std::string renderer_name_;

    bool compile_radial_shader();
public:
    void probe_gpu_performance();  // micro-benchmark to decide GPU vs CPU
private:
    // Context (platform-specific, stored as opaque pointers)
    void* gl_context_ = nullptr;
    void* gl_dc_ = nullptr;
    void* gl_hwnd_ = nullptr;

    size_t n_parts_ = 0;
    static constexpr int ROT_BINS = 360;

    bool init_context();
    bool compile_shader();
    void cleanup();
};

#endif // SLIC3R_GUI

// Factory: tries GPU first, falls back to CPU
std::unique_ptr<CollisionEvaluator> create_collision_evaluator();

} // namespace snuggle
