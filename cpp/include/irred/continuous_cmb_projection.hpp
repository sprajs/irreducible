#pragma once
#include "irred/numerics.hpp"
#include <cstddef>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace irred::projection {
// Flat scalar split-source geometry. All four channels are per conformal Mpc
// per the SAME signed initial mode. CLASS p means sqrt(6)*g*P_internal, with
// its positive E convention. These are finite-support transfers, not spectra.
struct ContinuousCmbSource {
  std::vector<double> k_mpc_inverse, eta_mpc;
  // eta-major, k-minor: channel[eta_index*k.size()+k_index].
  std::vector<double> t0, t1, t2, polarization;
  double observer_eta_mpc = 0;
  std::string producer_id, signed_mode_id, normalization_id;
};
struct ContinuousCmbPreparationPolicy {
  std::size_t maximum_k = 256, maximum_eta = 16384;
  std::size_t maximum_cells = 1048576;
  std::size_t maximum_payload_bytes = 64 * 1024 * 1024;
};
class ContinuousCmbProjection {
public:
  ContinuousCmbProjection() = default;
  ContinuousCmbProjection(const ContinuousCmbProjection &) = delete;
  ContinuousCmbProjection &operator=(const ContinuousCmbProjection &) = delete;
  ContinuousCmbProjection(ContinuousCmbProjection &&) noexcept;
  ContinuousCmbProjection &operator=(ContinuousCmbProjection &&) noexcept;
  numerics::Status status() const noexcept { return status_; }
  const ContinuousCmbSource *source() const noexcept { return source_.get(); }
  std::optional<std::size_t> retained_payload_bound() const noexcept;
private:
  numerics::Status status_ = numerics::Status::invalid_input;
  std::shared_ptr<const ContinuousCmbSource> source_;
  friend ContinuousCmbProjection prepare_continuous_cmb_projection(
      ContinuousCmbSource &&, ContinuousCmbPreparationPolicy);
  friend struct ContinuousCmbProjectionAccess;
};
// Validate/cap before consuming. Successful acquisition moves the original
// buffers exactly once; callers relinquish mutable aliases to those buffers.
ContinuousCmbProjection prepare_continuous_cmb_projection(
    ContinuousCmbSource &&, ContinuousCmbPreparationPolicy = {});
inline constexpr unsigned continuous_temperature = 1, continuous_e_mode = 2;
struct ContinuousCmbPolicy {
  // Estimates concern integration of the declared linear source interpolation.
  // Producer-grid, omitted-support and radial/libm errors remain separate.
  double absolute_tolerance = 1e-9, relative_tolerance = 1e-5;
  std::size_t maximum_multipoles = 128;
  std::size_t maximum_kernel_evaluations = 2000000;
  std::size_t maximum_payload_bytes = 128 * 1024 * 1024;
  unsigned maximum_depth = 20;
};
struct ContinuousCmbRow {
  std::size_t multipole_index = 0, k_index = 0;
  unsigned ell = 0;
  double k_mpc_inverse = 0;
  numerics::Status status = numerics::Status::invalid_input;
  std::optional<double> temperature, e_mode;
  double temperature_quadrature_estimate = 0, e_quadrature_estimate = 0;
  double temperature_arithmetic_estimate = 0, e_arithmetic_estimate = 0;
  // Absence is deliberate: time refinement earns neither of these estimates.
  std::optional<double> radial_error_estimate, source_grid_error_estimate;
  std::size_t kernel_evaluations = 0;
  std::size_t completed_source_cells = 0;
};
struct ContinuousCmbResult {
  numerics::Status status = numerics::Status::invalid_input;
  unsigned requested_outputs = 0;
  // SAME original acquired source, also retained on evaluation refusal.
  std::shared_ptr<const ContinuousCmbSource> source_owner;
  std::vector<ContinuousCmbRow> rows; // multipole-major, original k order
  std::size_t kernel_evaluations = 0, payload_bound = 0;
};
// Hard first slice: T ell<=64, E 2<=ell<=64, 0<=k*(eta0-eta)<=512.
// Within every supplied cell use linear interpolation, preserving all finite
// endpoints. No k interpolation, zero-outside rule, or source differentiation.
// Failure withholds dependent values and preserves estimates/attempted work.
ContinuousCmbResult project_continuous_cmb(
    const ContinuousCmbProjection &, std::span<const unsigned> multipoles,
    unsigned requested_outputs, ContinuousCmbPolicy = {});
inline constexpr std::string_view continuous_cmb_projection_model_id =
    "flat-supplied-continuous-scalar-split-source-finite-support/v1";
inline constexpr std::string_view continuous_cmb_projection_method_id =
    "linear-source-vector-simpson-phase-split-spherical-bessel/v1";
inline constexpr std::string_view continuous_cmb_projection_arithmetic_id =
    "binary64-source-wide-nearest-standard-spherical-bessel/v1";
} // namespace irred::projection
