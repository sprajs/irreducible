#pragma once
#include "irred/numerics.hpp"
#include "irred/quantities.hpp"
#include <cstddef>
#include <cstdint>
#include <span>
#include <string_view>
#include <vector>
namespace irred::cosmology {
enum class Model : std::uint32_t {
  flat_lcdm_late_v1 = 0,
  constant_q_flat_v1 = 1,
  flat_cpl_late_v1 = 2
};
enum class Convention : std::uint32_t {
  geometric_same_redshift,
  released_zhd_zhel
};
enum class Status : std::uint32_t {
  ok,
  invalid_input,
  unsupported_domain,
  incompatible_convention,
  numerical_failure,
  work_limit
};
struct Parameters {
  Model model = Model::flat_lcdm_late_v1;
  double h0_km_s_mpc = 70;
  double omega_m = 0.3;  // active for LCDM/CPL; must be zero for constant-q
  double constant_q = 0; // active only constant-q; zero for LCDM/CPL
  double w0 = -1,
         wa = 0; // CPL retained state; canonical inactive defaults otherwise
};
// No default constructor: CPL admission requires all four explicit parameters.
// Existing prepare(Parameters) accepts only the two legacy models.
struct CplParameters {
  double h0_km_s_mpc, omega_m, w0, wa;
  CplParameters(double h0, double matter, double present_w, double evolution_w)
      : h0_km_s_mpc(h0), omega_m(matter), w0(present_w), wa(evolution_w) {}
};
struct Query {
  double z_expansion;
  double z_observer;
  Convention convention;
};
struct Policy {
  numerics::IntegrationPolicy integration{1e-13, 1e-12, 100000, 30};
  std::size_t maximum_queries = 4096;
  std::size_t maximum_total_evaluations = 2000000;
};
// Values have scientific meaning only when status==ok. Integral error estimates
// are empirical F02 diagnostics, not proven bounds or automatic qualification.
struct Slot {
  Query source{};
  std::string_view luminosity_equation_id, shape_equation_id;
  Status status = Status::invalid_input;
  double expansion_E = 0, h_km_s_mpc = 0;
  double radial_integral = 0, radial_mpc = 0, transverse_mpc = 0;
  double angular_diameter_mpc = 0, luminosity_mpc = 0;
  double dimensionless_luminosity_shape =
      0; // (1+z_observer)*I; no absolute H0 information
  double lookback_seconds = 0, volume_mpc3_per_sr_per_redshift = 0;
  double deceleration_q = 0, jerk = 0;
  double radial_integral_error = 0, lookback_integral_error = 0;
  std::size_t evaluations = 0;
  numerics::Status numerical_status = numerics::Status::invalid_input;
};
struct BatchResult {
  Status status = Status::invalid_input;
  std::vector<Slot> slots;
};
class Background {
public:
  Status status() const noexcept { return status_; }
  const Parameters &parameters() const noexcept { return parameters_; }
  std::string_view model_id() const noexcept;
  static constexpr std::string_view constants_id = irred::constant_set_id;
  static constexpr std::string_view radial_equation_id =
      "P01/flat-radial-comoving-distance/v1";
  BatchResult evaluate_batch(std::span<const Query>, Policy) const;

private:
  Status status_ = Status::invalid_input;
  Parameters parameters_{};
  double hubble_distance_mpc_ = 0, hubble_time_seconds_ = 0;
  friend Background prepare(Parameters);
  friend Background prepare_cpl(CplParameters);
};
// Flat FLRW distance geometry. LCDM is radiation-free, nonnegative
// matter/Lambda components normalized at z0; constant-q is kinematic, no
// early-time assertion. Nonzero physical bundle outputs must be normal
// binary64; underflow/overflow fails the bundle, separately from physical/model
// domain. Analytic z0 zeros valid. First supported late-time interval z=[0,5];
// H0 consumer qualification40..100 is separate from finite-positive arithmetic
// representability. No age/rd/CMB.
Background prepare(Parameters);
// Radiation-free flat CPL, Om[0,1], w0[-2,0], wa[-2,2], z[0,5].
// Background only: no dark-energy perturbation/crossing closure or early epoch.
Background prepare_cpl(CplParameters);
} // namespace irred::cosmology
