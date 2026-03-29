#ifndef ARRANGE_HPP
#define ARRANGE_HPP

#include "ExPolygon.hpp"
#include "PrintConfig.hpp"
#include "Print.hpp"

#define BED_SHRINK_SEQ_PRINT 5

namespace Slic3r {

class BoundingBox;

namespace arrangement {

/// A geometry abstraction for a circular print bed. Similarly to BoundingBox.
class CircleBed {
    Point center_;
    double radius_;
public:

    inline CircleBed(): center_(0, 0), radius_(std::nan("")) {}
    explicit inline CircleBed(const Point& c, double r): center_(c), radius_(r) {}

    inline double radius() const { return radius_; }
    inline const Point& center() const { return center_; }
};

/// Representing an unbounded bed.
struct InfiniteBed {
    Point center;
    explicit InfiniteBed(const Point &p = {0, 0}): center{p} {}
};

/// A logical bed representing an object not being arranged. Either the arrange
/// has not yet successfully run on this ArrangePolygon or it could not fit the
/// object due to overly large size or invalid geometry.
static const constexpr int UNARRANGED = -1;

/// Input/Output structure for the arrange() function. The poly field will not
/// be modified during arrangement. Instead, the translation and rotation fields
/// will mark the needed transformation for the polygon to be in the arranged
/// position. These can also be set to an initial offset and rotation.
///
/// The bed_idx field will indicate the logical bed into which the
/// polygon belongs: UNARRANGED means no place for the polygon
/// (also the initial state before arrange), 0..N means the index of the bed.
/// Zero is the physical bed, larger than zero means a virtual bed.
struct ArrangePolygon {
    ExPolygon poly;                 /// The 2D silhouette to be arranged
    Vec2crd   translation{0, 0};    /// The translation of the poly
    double    rotation{0.0};        /// The rotation of the poly in radians
    coord_t   inflation = 0;        /// Arrange with inflated polygon
    int       bed_idx{UNARRANGED};  /// To which logical bed does poly belong...
    int       priority{0};
    //BBS: add locked_plate to indicate whether it is in the locked plate
    int       locked_plate{ -1 };
    bool      is_virt_object{ false };
    bool      is_extrusion_cali_object{ false };
    bool      is_wipe_tower{ false };
    bool      has_tree_support{false};
    //BBS: add row/col for sudoku-style layout
    int       row{0};
    int       col{0};
    std::vector<int> extrude_ids{};      /// extruder_id for least extruder switch
    int filament_temp_type{ -1 };
    int       bed_temp{0};         ///bed temperature for different material judge
    int       print_temp{0};      ///print temperature for different material judge
    int       first_bed_temp{ 0 };      ///first layer bed temperature for different material judge
    int       first_print_temp{ 0 };      ///first layer print temperature for different material judge
    int       vitrify_temp{ 0 };   // max bed temperature for material compatibility, which is usually the filament vitrification temp
    int       itemid{ 0 };         // item id in the vector, used for accessing all possible params like extrude_id
    int       is_applied{ 0 };     // transform has been applied
    double    height{ 0 };         // item height
    double    brim_width{ 0 };     // brim width
    std::string name;

    // If empty, any rotation is allowed (currently unsupported)
    // If only a zero is there, no rotation is allowed
    std::vector<double> allowed_rotations = {0.};

    /// Projected 2D triangle vertices for direct bitmap rasterization (concave mode).
    /// Every 3 consecutive Points form one triangle. Set during prepare(),
    /// consumed by BitmapArranger which scanline-fills each triangle into the
    /// collision bitmap — same as GPU rasterization but on CPU. No Clipper needed.
    Points concave_triangles;
    /// Z coordinates of triangle vertices (parallel to concave_triangles).
    /// Used for 3D-aware nesting to assign triangles to Z slices.
    std::vector<float> concave_z;

    /// Optional setter function which can store arbitrary data in its closure
    std::function<void(const ArrangePolygon&)> setter = nullptr;

    /// Helper function to call the setter with the arrange data arguments
    void apply() {
        if (setter && !is_applied) {
            setter(*this);
            is_applied = 1;
        }
    }

    /// Test if arrange() was called previously and gave a successful result.
    bool is_arranged() const { return bed_idx != UNARRANGED; }

    inline ExPolygon transformed_poly() const
    {
        ExPolygon ret = poly;
        ret.rotate(rotation);
        ret.translate(translation.x(), translation.y());

        return ret;
    }
};

using ArrangePolygons = std::vector<ArrangePolygon>;

struct ArrangeParams {

    /// The minimum distance which is allowed for any
    /// pair of items on the print bed in any direction.
    coord_t min_obj_distance = 0;

    /// The accuracy of optimization.
    /// Goes from 0.0 to 1.0 and scales performance as well
    float accuracy = 1.f;

    /// Allow parallel execution.
    bool parallel = true;

    bool allow_rotations = false;

    bool do_final_align = true;

    //BBS: add specific arrange params
    bool  allow_multi_materials_on_same_plate = true;
    bool  avoid_extrusion_cali_region         = true;
    bool  is_seq_print                        = false;
    bool  align_to_y_axis                     = false;
    bool  use_concave_hulls                   = false;
    bool  allow_multi_plate                   = true;
    bool  consolidate_plates                  = false;
    double rotation_step_rad                  = PI / 4.;
    bool  avoid_purge_pad                     = false;
    int   purge_pad_edge                      = 0;      // 0=front, 1=back, 2=left, 3=right
    float purge_pad_mm                        = 5.f;
    float bitmap_resolution_mm                = 0.5f;
    bool  nesting_3d                          = false;
    float slice_height_mm                     = 12.7f;
    float z_clearance_mm                      = 2.0f;
    int   compaction_mode                     = 3;  // 0=none, 1=coarse bitmap, 2=reverse skyline, 3=both
    bool  gravity_compact                     = false; // slide parts toward Y=0 to tighten layout
    bool  best_fit_compact                    = true;  // try all plates, pick lowest Y (vs first-fit)
    float bed_shrink_x = 1;
    float bed_shrink_y = 1;
    float brim_skirt_distance = 0;
    float clearance_height_to_rod = 0;
    float clearance_height_to_lid = 0;
    float clearance_radius = 0;
    float object_skirt_offset = 0;
    float nozzle_height = 0;
    float printable_height = 256.0;
    Vec2d align_center{ 0.5,0.5 };

    ArrangePolygons excluded_regions;   // regions cant't be used
    ArrangePolygons nonprefered_regions; // regions can be used but not prefered

    /// Progress indicator callback called when an object gets packed.
    /// The unsigned argument is the number of items remaining to pack.
    std::function<void(unsigned, std::string)> progressind = [](unsigned st, std::string str = "") {
        std::cout << "st=" << st << ", " << str << std::endl;
    };

    std::function<void(const ArrangePolygon &)> on_packed;

    /// A predicate returning true if abort is needed.
    std::function<bool(void)>     stopcondition;

    ArrangeParams() = default;
    explicit ArrangeParams(coord_t md) : min_obj_distance(md) {}
    // to json format
    std::string to_json() const{
        std::string ret = "{";
        ret += "\"min_obj_distance\":" + std::to_string(min_obj_distance) + ",";
        ret += "\"accuracy\":" + std::to_string(accuracy) + ",";
        ret += "\"parallel\":" + std::to_string(parallel) + ",";
        ret += "\"allow_rotations\":" + std::to_string(allow_rotations) + ",";
        ret += "\"do_final_align\":" + std::to_string(do_final_align) + ",";
        ret += "\"allow_multi_materials_on_same_plate\":" + std::to_string(allow_multi_materials_on_same_plate) + ",";
        ret += "\"avoid_extrusion_cali_region\":" + std::to_string(avoid_extrusion_cali_region) + ",";
        ret += "\"is_seq_print\":" + std::to_string(is_seq_print) + ",";
        ret += "\"bed_shrink_x\":" + std::to_string(bed_shrink_x) + ",";
        ret += "\"bed_shrink_y\":" + std::to_string(bed_shrink_y) + ",";
        ret += "\"brim_skirt_distance\":" + std::to_string(brim_skirt_distance) + ",";
        ret += "\"clearance_height_to_rod\":" + std::to_string(clearance_height_to_rod) + ",";
        ret += "\"clearance_height_to_lid\":" + std::to_string(clearance_height_to_lid) + ",";
        ret += "\"clearance_radius\":" + std::to_string(clearance_radius) + ",";
        ret += "\"printable_height\":" + std::to_string(printable_height) + ",";
        ret += "\"use_concave_hulls\":" + std::to_string(use_concave_hulls) + ",";
        ret += "\"allow_multi_plate\":" + std::to_string(allow_multi_plate) + ",";
        ret += "\"consolidate_plates\":" + std::to_string(consolidate_plates) + ",";
        ret += "\"rotation_step_rad\":" + std::to_string(rotation_step_rad) + ",";
        ret += "\"align_to_y_axis\":" + std::to_string(align_to_y_axis) + ",";
        ret += "\"avoid_purge_pad\":" + std::to_string(avoid_purge_pad) + ",";
        ret += "\"purge_pad_edge\":" + std::to_string(purge_pad_edge) + ",";
        ret += "\"purge_pad_mm\":" + std::to_string(purge_pad_mm) + ",";
        ret += "\"bitmap_resolution_mm\":" + std::to_string(bitmap_resolution_mm) + ",";
        ret.pop_back();
        ret += "}";
        return ret;
    }

};

void update_arrange_params(ArrangeParams& params, const DynamicPrintConfig* print_cfg, const ArrangePolygons& selected);

void update_selected_items_inflation(ArrangePolygons& selected, const DynamicPrintConfig* print_cfg, ArrangeParams& params);

void update_unselected_items_inflation(ArrangePolygons& unselected, const DynamicPrintConfig* print_cfg, const ArrangeParams& params);

void update_selected_items_axis_align(ArrangePolygons& selected, const DynamicPrintConfig* print_cfg, const ArrangeParams& params);

Points get_shrink_bedpts(const DynamicPrintConfig* print_cfg, const ArrangeParams& params);

/**
 * \brief Arranges the input polygons.
 *
 * When use_concave_hulls is set in ArrangeParams, the BitmapArranger is
 * used instead of libnest2d, enabling tighter packing for irregular
 * shapes (crescents, L-brackets, U-channels).
 *
 * \param items Input vector of ArrangePolygons. The transformation, rotation
 * and bin_idx fields will be changed after the call finished and can be used
 * to apply the result on the input polygon.
 */
template<class TBed> void arrange(ArrangePolygons &items, const ArrangePolygons &excludes, const TBed &bed, const ArrangeParams &params = {});

// A dispatch function that determines the bed shape from a set of points.
template<> void arrange(ArrangePolygons &items, const ArrangePolygons &excludes, const Points &bed, const ArrangeParams &params);

template<> void arrange(ArrangePolygons &items, const ArrangePolygons &excludes, const BoundingBox &bed, const ArrangeParams &params);
extern template void arrange(ArrangePolygons &items, const ArrangePolygons &excludes, const CircleBed &bed, const ArrangeParams &params);
extern template void arrange(ArrangePolygons &items, const ArrangePolygons &excludes, const Polygon &bed, const ArrangeParams &params);
extern template void arrange(ArrangePolygons &items, const ArrangePolygons &excludes, const InfiniteBed &bed, const ArrangeParams &params);

inline void arrange(ArrangePolygons &items, const Points &bed, const ArrangeParams &params = {}) { arrange(items, {}, bed, params); }
inline void arrange(ArrangePolygons &items, const BoundingBox &bed, const ArrangeParams &params = {}) { arrange(items, {}, bed, params); }
inline void arrange(ArrangePolygons &items, const CircleBed &bed, const ArrangeParams &params = {}) { arrange(items, {}, bed, params); }
inline void arrange(ArrangePolygons &items, const Polygon &bed, const ArrangeParams &params = {}) { arrange(items, {}, bed, params); }
inline void arrange(ArrangePolygons &items, const InfiniteBed &bed, const ArrangeParams &params = {}) { arrange(items, {}, bed, params); }

}} // namespace Slic3r::arrangement

#endif // MODELARRANGE_HPP
