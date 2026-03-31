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

namespace snuggle {

// Forward declarations (defined in snuggle_nester.hpp)
struct PartInfo;
struct Individual;

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

    // Evaluate collision_count and oob_count for every individual in the batch.
    // Must write ind.collision_count and ind.oob_count for each individual.
    // Must NOT write ind.fitness (caller handles that).
    virtual void evaluate_batch(
        std::vector<Individual>& pop,
        const std::vector<PartInfo>& parts,
        float bed_width_mm,
        float bed_height_mm) = 0;
};

// CPU fallback -- extracts the collision/bounds logic from SnuggleNester::evaluate()
class CpuCollisionEvaluator : public CollisionEvaluator {
public:
    void upload_grids(const std::vector<PartInfo>& parts,
                      const std::vector<std::vector<VoxelGrid>>& rot_cache) override;
    void evaluate_batch(std::vector<Individual>& individuals,
                        const std::vector<PartInfo>& parts,
                        float bed_w, float bed_h) override;
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

    void upload_grids(const std::vector<PartInfo>& parts,
                      const std::vector<std::vector<VoxelGrid>>& rot_cache) override;
    void evaluate_batch(std::vector<Individual>& individuals,
                        const std::vector<PartInfo>& parts,
                        float bed_w, float bed_h) override;

private:
    bool available_ = false;
    // GL handles
    unsigned int program_ = 0;
    unsigned int voxel_ssbo_ = 0;
    unsigned int meta_ssbo_ = 0;
    unsigned int placement_ssbo_ = 0;
    unsigned int results_ssbo_ = 0;
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
