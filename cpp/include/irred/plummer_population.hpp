#pragma once
#include "irred/plummer_sphere.hpp"
#include <optional>
#include <span>
#include <vector>
namespace irred::gravity {
inline constexpr std::string_view plummer_population_id =
    "self-gravitating-isotropic-collisionless-Plummer-mass-traces-tracer-SI/v1";
struct PlummerPopulationPreparationPolicy {
  double relative_tolerance=1e-10;
  std::size_t maximum_radii=4096, maximum_native_evaluations=8192;
  std::size_t maximum_payload_bytes=4*1024*1024;
};
struct PlummerPopulationEvaluationPolicy {
  double relative_tolerance=1e-10;
  std::size_t maximum_rows=65536, maximum_velocity_evaluations=65536;
  std::size_t maximum_payload_bytes=32*1024*1024;
};
struct PlummerPopulationRadius {
  double radius_metres=0;
  numerics::Status status=numerics::Status::invalid_input;
  numerics::ScalarResult density_kg_m3, relative_potential_m2_s2, enclosed_mass_kg,
                         projected_surface_density_kg_m2;
  numerics::ScalarResult escape_speed_m_s, one_axis_variance_m2_s2,
                         projected_los_variance_m2_s2;
};
struct PlummerVelocityRequest { std::size_t radius_index; double speed_m_s; };
enum class PlummerEnergySupport { not_evaluated, bound, outside, ambiguous };
struct PlummerVelocityRow {
  numerics::Status status=numerics::Status::invalid_input;
  PlummerVelocityRequest request{};
  double radius_metres=0;
  PlummerEnergySupport support=PlummerEnergySupport::not_evaluated;
  // Wide energy and inherited arithmetic diagnostic preserve the sign test.
  // The diagnostic is empirical, not a certified support enclosure.
  long double binding_energy_m2_s2=0, binding_energy_absolute_error_m2_s2=0;
  numerics::ScalarResult distribution_function_kg_s3_m6,
                         vector_velocity_density_s3_m3, speed_density_s_m;
};
struct PlummerPopulationBatch {
  numerics::Status status=numerics::Status::invalid_input;
  std::optional<PlummerSphere> source;
  std::size_t requested_rows=0, velocity_evaluations=0;
  std::vector<PlummerVelocityRow> rows;
};
class PlummerPopulation {
public:
  PlummerPopulation()=default;
  PlummerPopulation(const PlummerPopulation&)=default;
  PlummerPopulation& operator=(const PlummerPopulation&)=default;
  PlummerPopulation(PlummerPopulation&&) noexcept;
  PlummerPopulation& operator=(PlummerPopulation&&) noexcept;
  numerics::Status status() const noexcept { return status_; }
  const PlummerSphere* source() const noexcept { return source_?&*source_:nullptr; }
  std::span<const PlummerPopulationRadius> radii() const noexcept { return radii_; }
  std::size_t requested_radii() const noexcept { return requested_radii_; }
  std::size_t native_evaluations() const noexcept { return native_evaluations_; }
  PlummerPopulationBatch evaluate(std::span<const PlummerVelocityRequest>,
                                  PlummerPopulationEvaluationPolicy={}) const;
private:
  numerics::Status status_=numerics::Status::invalid_input;
  std::optional<PlummerSphere> source_;
  std::vector<PlummerPopulationRadius> radii_;
  long double log_coefficient_=0, log_coefficient_error_=0;
  std::size_t requested_radii_=0, native_evaluations_=0;
  friend PlummerPopulation prepare_plummer_population(PlummerSphere,
      std::span<const double>,PlummerPopulationPreparationPolicy);
};
std::optional<std::size_t> plummer_population_payload_bound(std::size_t radii) noexcept;
std::optional<std::size_t> plummer_population_batch_payload_bound(std::size_t rows) noexcept;
// Owns the ordered physical/projected radius grid. All required radius states
// must be admitted; failed preparations retain available radius diagnostics.
PlummerPopulation prepare_plummer_population(PlummerSphere,std::span<const double>,
                                             PlummerPopulationPreparationPolicy={});
}
