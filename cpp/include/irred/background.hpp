#pragma once
#include "irred/numerics.hpp"
#include "irred/quantities.hpp"
#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <string_view>
#include <variant>
#include <vector>
namespace irred::cosmology {
struct FlatFLRW {};
struct LCDM {
  double omega_m;
  explicit LCDM(double x) : omega_m(x) {}
};
struct ConstantQ {
  double q;
  explicit ConstantQ(double x) : q(x) {}
};
struct CPL {
  double omega_m, w0, wa;
  CPL(double m, double w, double a) : omega_m(m), w0(w), wa(a) {}
};
struct FixedFiveBinQ {
  std::array<double, 5> q;
  explicit FixedFiveBinQ(std::array<double, 5> x) : q(x) {}
};
using ExpansionSpec = std::variant<LCDM, ConstantQ, CPL, FixedFiveBinQ>;
inline constexpr std::array<double, 6> piecewise_q_edges{0, .1, .3, .6, 1, 2.5};
enum class Convention : std::uint32_t {
  geometric_same_redshift,
  released_zhd_zhel
};
struct Observer {
  double redshift;
  Convention convention;
};
struct PhysicalScale {
  double h0_km_s_mpc;
  explicit PhysicalScale(double h) : h0_km_s_mpc(h) {}
};
enum class Observable : std::uint32_t {
  radial = 1,
  luminosity_shape = 2,
  clock = 4,
  flat_distances_volume = 8,
  kinematics = 16,
  expansion = 32
};
constexpr std::uint32_t operator|(Observable a, Observable b) {
  return (std::uint32_t)a | (std::uint32_t)b;
}
struct Request {
  double z_expansion;
  std::uint32_t requested;
  std::optional<Observer> observer;
  std::optional<PhysicalScale> physical_scale;
  Request(double z, std::uint32_t outputs, std::optional<Observer> o = {},
          std::optional<PhysicalScale> h = {})
      : z_expansion(z), requested(outputs), observer(o), physical_scale(h) {}
};
enum class Status : std::uint32_t {
  ok,
  invalid_input,
  unsupported_domain,
  incompatible_convention,
  numerical_failure,
  work_limit
};
enum class Availability : std::uint32_t {
  not_requested,
  available,
  unavailable,
  failed
};
template <class T> struct Outcome {
  Availability availability = Availability::not_requested;
  Status status = Status::invalid_input;
  numerics::Status numerical_status = numerics::Status::invalid_input;
  std::optional<T> value;
};
namespace detail {
struct RadialAccess;
}
struct ExpansionValue {
  double expansion_E = 0;
  Outcome<double> h_km_s_mpc;

private:
  long double precise_E_ = 0;
  friend class Expansion;
  friend struct detail::RadialAccess;
};
struct Radial {
  double expansion_E = 0, integral = 0, error_estimate = 0;

private:
  long double precise_E_ = 0, precise_integral_ = 0;
  friend class Expansion;
  friend struct detail::RadialAccess;
};
struct Clock {
  double integral = 0, error_estimate = 0;
  Outcome<double> lookback_seconds;
};
struct Distances {
  double radial_mpc = 0, transverse_mpc = 0, angular_diameter_mpc = 0,
         luminosity_mpc = 0, volume_mpc3_per_sr_per_redshift = 0;
};
enum class QConvention : std::uint32_t {
  not_assessed,
  interior_constant_bin,
  right_limit_at_internal_jump,
  right_limit_at_zero,
  left_limit_at_final_endpoint,
  ordinary_smooth_model
};
enum class JerkAvailability : std::uint32_t {
  not_assessed,
  ordinary_within_bin,
  one_sided_endpoint,
  unavailable_at_jump
};
struct Kinematics {
  double q = 0;
  std::optional<double> jerk, q0_within_piecewise_model;
  QConvention q_convention = QConvention::not_assessed;
  JerkAvailability jerk_availability = JerkAvailability::not_assessed;
  std::size_t bin = 0;
};
struct Work {
  std::size_t callbacks = 0, segment_visits = 0;
};
struct Slot {
  Request source;
  Status admission_status = Status::invalid_input;
  std::optional<std::size_t> node_index;
  Outcome<ExpansionValue> expansion;
  Outcome<Radial> radial;
  Outcome<double> luminosity_shape;
  Outcome<Clock> clock;
  Outcome<Distances> physical;
  Outcome<Kinematics> kinematics;
  explicit Slot(Request r) : source(r) {}
};
struct NodeWork {
  double z_expansion;
  Work work;
};
struct BatchResult {
  Status status = Status::invalid_input;
  std::vector<Slot> slots;
  std::vector<NodeWork> nodes;
  Work work;
};
struct EvaluationPolicy {
  std::optional<numerics::IntegrationPolicy> integration;
  std::size_t maximum_queries = 0, maximum_callbacks = 0,
              maximum_segment_visits = 0, maximum_native_bytes = 0;
};
class Expansion {
public:
  Status status() const noexcept { return status_; }
  const ExpansionSpec &specification() const noexcept { return spec_; }
  std::string_view model_id() const noexcept;
  static constexpr std::string_view constants_id = irred::constant_set_id;
  static constexpr std::string_view radial_equation_id =
      "P01/flat-radial-comoving-distance/v1";
  // Conservative simultaneous dynamic payload for this supported container
  // implementation, excluding borrowed inputs, stack headers, allocator
  // metadata and RSS. Zero requests allocate nothing; overflow has no
  // representable bound.
  static std::optional<std::size_t>
  workspace_payload_bound(std::size_t) noexcept;
  BatchResult evaluate(std::span<const Request>, EvaluationPolicy) const;

private:
  explicit Expansion(ExpansionSpec s) : spec_(std::move(s)) {}
  ExpansionSpec spec_;
  Status status_ = Status::invalid_input;
  std::array<long double, 5> start_E_{};
  friend Expansion prepare(ExpansionSpec, FlatFLRW);
};
// Only these four compiled late-time hypotheses. No implicit inactive fields,
// physical H0, curvature, expressions, priors or early-ruler physics.
Expansion prepare(ExpansionSpec, FlatFLRW);
} // namespace irred::cosmology
