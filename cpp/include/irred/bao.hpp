#pragma once
#include "irred/background.hpp"
#include "irred/piecewise_background.hpp"
#include "irred/statistics.hpp"
#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <string_view>
#include <vector>
namespace irred::bao {
// Flat late-time distance ratios. Ruler is an empirically free scale, not a
// computed drag-epoch or last-scattering sound horizon.
enum class Observable : std::uint32_t {
  transverse_over_ruler = 0,
  hubble_over_ruler = 1,
  volume_over_ruler = 2
};
struct Query {
  double z;
  Observable observable;
};
struct Ruler {
  double h0_rd_km_s;
  explicit Ruler(double product_km_s) : h0_rd_km_s(product_km_s) {}
};
struct Policy {
  cosmology::Policy background;
  std::size_t maximum_queries = 4096;
  // Conservative dynamic payload peak: returned slots + admitted queries +
  // index map + Background slots; excludes caller input and allocator overhead.
  std::size_t maximum_native_bytes = 1048576;
};
struct Slot {
  Query source{};
  cosmology::Status status = cosmology::Status::invalid_input;
  numerics::Status numerical_status = numerics::Status::invalid_input;
  // Meaningful only when status==ok; exact transverse/volume zero at z0 valid.
  double dimensionless_value = 0;
  std::size_t evaluations = 0;
};
struct Batch {
  cosmology::Status status = cosmology::Status::invalid_input;
  std::size_t evaluations = 0;
  std::vector<Slot> slots;
};
inline constexpr std::string_view equation_id = "P01/flat-free-ruler-BAO/v1";
inline constexpr std::string_view ruler_convention_id =
    "P01/free-H0rd-km-s-no-early-physics/v1";
// z[0,5], H0rd[5000,15000] km/s; owned Background's existing model domain and
// all-observable representability limits apply. Its computational H0 cancels.
// One background callback cap applies across the whole batch. No likelihood,
// joint-probe independence, automatic qualification or runtime CPL ABI claim.
Batch evaluate(const cosmology::Background &, Ruler, std::span<const Query>,
               Policy);

enum class RowRole : std::uint32_t {
  released_fitted_distance_summary = 0,
  synthetic_control = 1,
  unknown = UINT32_MAX
};
enum class CovarianceUnit : std::uint32_t {
  dimensionless_ratio_squared = 0,
  unknown = UINT32_MAX
};
struct DensityInput {
  std::vector<Query> queries;
  std::vector<double> observed, covariance;
  std::vector<std::string> ordered_ids;
  RowRole role = RowRole::unknown;
  CovarianceUnit covariance_unit = CovarianceUnit::unknown;
  std::string table_identity, covariance_identity, ordering_provenance,
      calibration_provenance, dependence_provenance;
};
struct DensityPolicy {
  Policy observables;
  // Required density policy is unset until the caller supplies every field.
  std::size_t maximum_models = 0, maximum_matrix_elements = 0,
              maximum_string_bytes = 0, maximum_native_bytes = 0;
  double maximum_forward_sensitivity = 0;
  numerics::Arithmetic arithmetic =
      static_cast<numerics::Arithmetic>(UINT32_MAX);
};
struct ModelQuery {
  cosmology::Background background;
  Ruler ruler;
};
struct DensitySlot {
  cosmology::Parameters attempted_model;
  double attempted_h0_rd_km_s;
  std::string model_id, equation_id, arithmetic_id;
  cosmology::Status background_status = cosmology::Status::invalid_input;
  numerics::Status numerical_status = numerics::Status::invalid_input;
  statistics::GaussianResult result;
  std::vector<double> predictions, residuals;
  std::size_t evaluations = 0;
};
struct DensityBatch {
  statistics::DensityStatus status = statistics::DensityStatus::invalid_input;
  numerics::Status numerical_status = numerics::Status::invalid_input;
  std::vector<DensitySlot> slots;
  std::size_t evaluations = 0;
};
// Distinct analytic entry: fixed five-bin q and explicitly free H0*r_d.
// The computational H0 used to construct geometry cancels from these ratios.
inline constexpr double piecewise_computational_h0_km_s_mpc = 70;
struct PiecewiseModelPoint {
  std::array<double, 5> q;
  Ruler ruler;
  PiecewiseModelPoint(std::array<double, 5> values, Ruler scale)
      : q(values), ruler(scale) {}
};
struct PiecewisePolicy {
  cosmology::PiecewisePolicy background;
  std::size_t maximum_queries = 4096, maximum_native_bytes = 1048576;
};
struct PiecewiseSlot {
  Query source{};
  cosmology::Status status = cosmology::Status::invalid_input;
  numerics::Status numerical_status = numerics::Status::invalid_input;
  double dimensionless_value = 0;
  std::size_t segment_visits = 0;
};
struct PiecewiseBatch {
  cosmology::Status status = cosmology::Status::invalid_input;
  std::size_t segment_visits = 0;
  std::vector<PiecewiseSlot> slots;
};
PiecewiseBatch evaluate_piecewise(const cosmology::PiecewiseBackground &, Ruler,
                                  std::span<const Query>, PiecewisePolicy);
struct PiecewiseDensityPolicy {
  PiecewisePolicy observables;
  std::size_t maximum_models = 0, maximum_native_bytes = 0;
  double maximum_forward_sensitivity = 0;
  numerics::Arithmetic arithmetic =
      static_cast<numerics::Arithmetic>(UINT32_MAX);
};
struct PiecewiseDensitySlot {
  PiecewiseModelPoint source;
  explicit PiecewiseDensitySlot(PiecewiseModelPoint point) : source(point) {}
  std::string model_id, equation_id, arithmetic_id;
  cosmology::Status background_status = cosmology::Status::invalid_input;
  numerics::Status numerical_status = numerics::Status::invalid_input;
  statistics::GaussianResult result;
  std::vector<double> predictions, residuals;
  std::size_t segment_visits = 0;
};
struct PiecewiseDensityBatch {
  statistics::DensityStatus status = statistics::DensityStatus::invalid_input;
  numerics::Status numerical_status = numerics::Status::invalid_input;
  std::vector<PiecewiseDensitySlot> slots;
  std::size_t segment_visits = 0;
};
class PreparedDensity {
public:
  statistics::DensityStatus status() const noexcept { return status_; }
  numerics::Status numerical_status() const noexcept {
    return numerical_status_;
  }
  const DensityInput &source() const noexcept { return source_; }
  const statistics::Metadata &metadata() const noexcept {
    return gaussian_.metadata();
  }
  DensityBatch evaluate(std::span<const ModelQuery>, DensityPolicy) const;
  // z[0,2.5], finite q_i[-3,2], free H0rd[5000,15000]; no extrapolation,
  // quadrature settings, priors or joint-probe independence claim. Derivative
  // absence at q jumps does not invalidate these geometry-only observables.
  PiecewiseDensityBatch evaluate_piecewise(std::span<const PiecewiseModelPoint>,
                                           PiecewiseDensityPolicy) const;

private:
  template <class Model, class Result, class EvaluationPolicy>
  Result evaluate_common(std::span<const Model>, EvaluationPolicy) const;
  DensityInput source_;
  statistics::Gaussian gaussian_;
  statistics::DensityStatus status_ = statistics::DensityStatus::invalid_input;
  numerics::Status numerical_status_ = numerics::Status::invalid_input;
  friend PreparedDensity prepare_density(const DensityInput &, DensityPolicy);
};
// Original full covariance/order retained. Units and identity declarations are
// mandatory; unknown provenance remains unknown. No source covariance repair.
PreparedDensity prepare_density(const DensityInput &, DensityPolicy);
} // namespace irred::bao
