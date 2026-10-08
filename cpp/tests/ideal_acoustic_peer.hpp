#pragma once

// Test-only independent conformal-time oscillator / spatial Einstein TRACE.
// The retained coefficients are shared ancestry.  The polynomial below is an
// explicitly identified high-precision SOURCE CONTROL, never a production H
// replacement or an independent physical mapping.  No momentum Einstein RHS,
// constraint projection, finite-zeta reset, or damping advances this reference.
// Authored from the pinned ideal-acoustic source equations, not copied code.
#include "../src/thermal_conformal_epoch.hpp"
#include "irred/ideal_acoustic.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <limits>
#include <mpfr.h>
#include <optional>
#include <span>
#include <string_view>
#include <vector>

namespace ideal_acoustic_peer {
using Status = irred::numerics::Status;
inline constexpr std::string_view method_id =
    "test-retained-polynomial-eta-TRACE-weighted-oscillator-DP54-kappa0/v1";
inline constexpr std::string_view source_control_id =
    "test-exact-import-emitted-map-private-retained-Lambda-polynomial/v1";

struct Stop {
  Status status;
};
struct Budget {
  Status status = Status::ok;
  std::size_t rhs = 0, rhs_this_k = 0, background = 0, writes = 0;
  std::size_t endpoints = 0, attempts = 0, artifact_bytes = 0;
  std::size_t maximum_rhs = 4000000, maximum_rhs_per_k = 2000000;
  std::size_t maximum_background = 6000000, maximum_writes = 100000000;
  std::size_t maximum_endpoints = 500000;
  std::size_t maximum_payload = 64 * 1024 * 1024;
  std::size_t maximum_artifact_bytes = 128 * 1024 * 1024;
  std::size_t retained_payload_prefix = 0;
  // An arithmetic library's temporary/cache allocation is not proved by the
  // inline scalar layout.  A later whole-call allocation check must earn this.
  bool external_arithmetic_payload_admitted = false;
  void charge(std::size_t &counter, std::size_t n, std::size_t cap) {
    if (status != Status::ok)
      throw Stop{status};
    if (n > cap || counter > cap - n) {
      status = Status::work_limit;
      throw Stop{status};
    }
    counter += n;
  }
  void scalar(std::size_t n = 1) { charge(writes, n, maximum_writes); }
  void stage() {
    // Six physical derivatives and the a clock are distinct RHS owners.
    charge(rhs, 2, maximum_rhs);
    charge(rhs_this_k, 2, maximum_rhs_per_k);
  }
  void query() { charge(background, 1, maximum_background); }
  void endpoint() { charge(endpoints, 1, maximum_endpoints); }
};
inline thread_local Budget *active_budget = nullptr;
struct BudgetScope {
  Budget *previous;
  explicit BudgetScope(Budget &b) : previous(active_budget) {
    active_budget = &b;
  }
  ~BudgetScope() { active_budget = previous; }
};
inline void write(std::size_t n = 1) {
  if (active_budget)
    active_budget->scalar(n);
}

// Fixed inline limbs.  At 400/300 bits, 2^10>10^3 establishes at least the
// requested 120/90 decimal digits.  MPFR internal integer scratch remains an
// explicit external-runtime resource gate.  Every actual scalar construction,
// copy, arithmetic destination and emitted floating destination is counted.
template <unsigned DecimalDigits> class Real {
  static_assert(DecimalDigits == 90 || DecimalDigits == 120);
  static constexpr mpfr_prec_t bits = DecimalDigits * 10 / 3;
  std::array<mp_limb_t, (bits + GMP_NUMB_BITS - 1) / GMP_NUMB_BITS> limbs_{};
  mpfr_t v_;
  void init() {
    mpfr_custom_init(limbs_.data(), bits);
    mpfr_custom_init_set(v_, MPFR_ZERO_KIND, 0, bits, limbs_.data());
  }

public:
  Real() {
    write();
    init();
  }
  explicit Real(long value) {
    write();
    init();
    mpfr_set_si(v_, value, MPFR_RNDN);
  }
  explicit Real(double value) {
    write();
    init();
    mpfr_set_d(v_, value, MPFR_RNDN);
  }
  explicit Real(long double value) {
    write();
    init();
    mpfr_set_ld(v_, value, MPFR_RNDN);
  }
  Real(const Real &o) {
    write();
    init();
    mpfr_set(v_, o.v_, MPFR_RNDN);
  }
  Real(Real &&o) : Real(static_cast<const Real &>(o)) {}
  Real &operator=(const Real &o) {
    if (this != &o) {
      write();
      mpfr_set(v_, o.v_, MPFR_RNDN);
    }
    return *this;
  }
  Real &operator=(Real &&o) { return *this = static_cast<const Real &>(o); }
  bool finite() const { return mpfr_number_p(v_); }
  bool zero() const { return mpfr_zero_p(v_); }
  long double wide() const {
    write();
    const long double out = mpfr_get_ld(v_, MPFR_RNDN);
    if (!finite() || !std::isfinite(out) ||
        (!zero() && std::fpclassify(out) != FP_NORMAL))
      throw Stop{Status::overflow};
    return out;
  }
  long double upward_wide() const {
    write();
    const long double out = mpfr_get_ld(v_, MPFR_RNDU);
    if (!finite() || !std::isfinite(out) ||
        (!zero() && std::fpclassify(out) != FP_NORMAL))
      throw Stop{Status::overflow};
    return out;
  }
  friend Real operator+(const Real &a, const Real &b) {
    Real r;
    write();
    mpfr_add(r.v_, a.v_, b.v_, MPFR_RNDN);
    return r;
  }
  friend Real operator-(const Real &a, const Real &b) {
    Real r;
    write();
    mpfr_sub(r.v_, a.v_, b.v_, MPFR_RNDN);
    return r;
  }
  friend Real operator*(const Real &a, const Real &b) {
    Real r;
    write();
    mpfr_mul(r.v_, a.v_, b.v_, MPFR_RNDN);
    return r;
  }
  friend Real operator/(const Real &a, const Real &b) {
    Real r;
    write();
    mpfr_div(r.v_, a.v_, b.v_, MPFR_RNDN);
    return r;
  }
  friend Real operator-(const Real &a) {
    Real r;
    write();
    mpfr_neg(r.v_, a.v_, MPFR_RNDN);
    return r;
  }
  friend Real abs(const Real &a) {
    Real r;
    write();
    mpfr_abs(r.v_, a.v_, MPFR_RNDN);
    return r;
  }
  friend Real sqrt(const Real &a) {
    Real r;
    write();
    mpfr_sqrt(r.v_, a.v_, MPFR_RNDN);
    return r;
  }
  friend Real sin(const Real &a) {
    Real r;
    write();
    mpfr_sin(r.v_, a.v_, MPFR_RNDN);
    return r;
  }
  friend Real cos(const Real &a) {
    Real r;
    write();
    mpfr_cos(r.v_, a.v_, MPFR_RNDN);
    return r;
  }
  // The response's positive radii use MPFR directed elementary arithmetic.
  // This does not cover the finite DP path, state/slope arithmetic or external
  // library payload; those retain their separate admission/measurement gates.
  friend Real upper_add(const Real &a, const Real &b) {
    Real r;
    write();
    mpfr_add(r.v_, a.v_, b.v_, MPFR_RNDU);
    return r;
  }
  friend Real lower_add(const Real &a, const Real &b) {
    Real r;
    write();
    mpfr_add(r.v_, a.v_, b.v_, MPFR_RNDD);
    return r;
  }
  friend Real lower_subtract(const Real &a, const Real &b) {
    Real r;
    write();
    mpfr_sub(r.v_, a.v_, b.v_, MPFR_RNDD);
    return r;
  }
  friend Real upper_multiply(const Real &a, const Real &b) {
    Real r;
    write();
    mpfr_mul(r.v_, a.v_, b.v_, MPFR_RNDU);
    return r;
  }
  friend Real upper_divide(const Real &a, const Real &b) {
    Real r;
    write();
    mpfr_div(r.v_, a.v_, b.v_, MPFR_RNDU);
    return r;
  }
  friend Real lower_multiply(const Real &a, const Real &b) {
    Real r;
    write();
    mpfr_mul(r.v_, a.v_, b.v_, MPFR_RNDD);
    return r;
  }
  friend Real lower_divide(const Real &a, const Real &b) {
    Real r;
    write();
    mpfr_div(r.v_, a.v_, b.v_, MPFR_RNDD);
    return r;
  }
  friend Real upper_sqrt(const Real &a) {
    Real r;
    write();
    mpfr_sqrt(r.v_, a.v_, MPFR_RNDU);
    return r;
  }
  static Real unit_roundoff() {
    Real r;
    write();
    mpfr_set_ui_2exp(r.v_, 1, -bits, MPFR_RNDN);
    return r;
  }
  friend bool operator<(const Real &a, const Real &b) {
    return mpfr_less_p(a.v_, b.v_);
  }
  friend bool operator>(const Real &a, const Real &b) { return b < a; }
  friend bool operator<=(const Real &a, const Real &b) {
    return mpfr_lessequal_p(a.v_, b.v_);
  }
  friend bool operator>=(const Real &a, const Real &b) {
    return mpfr_greaterequal_p(a.v_, b.v_);
  }
  friend bool operator==(const Real &a, const Real &b) {
    return mpfr_equal_p(a.v_, b.v_);
  }
};
template <class T> T rat(long n, long d) { return T(n) / T(d); }
inline long double up_add(long double a, long double b);

struct Source {
  Status status = Status::invalid_input;
  long double h0 = 0, gamma = 0, baryon = 0, cdm = 0, lambda = 0;
  long double lambda_getter_loss = 0, normalization_estimate = 0;
  std::array<long double, 3> map_estimates{}; // gamma,baryon,CDM
  std::array<long double, 3> map_wide_values{}, map_operation_estimates{};
  bool map_witnesses_available = false;
  bool formal_dust = false; // equation control only; capture never sets this
  // Captured once from the real retained owner and original map.  No remap or
  // fresh closure is made here.  The member estimates are NOT zeroed merely
  // because the independent control uses more precision.
  static Source
  capture(const irred::cosmology::ThermalBackground &background,
          const irred::cosmology::ThermalPhysicalMapping &mapping) {
    Source s;
    const auto witness =
        irred::cosmology::detail::ThermalRetainedCoefficientAccess::capture(
            background);
    if (!witness || mapping.status != Status::ok || !mapping.model ||
        !mapping.scalar_witnesses || background.status() != Status::ok)
      return s;
    const auto &m = background.source();
    const auto &mapped = *mapping.model;
    if (!m.species.empty() || m.omega_massless_nonphoton != 0 ||
        !(m.h0_km_s_mpc > 0) || !(m.omega_gamma > 0) || !(m.omega_b > 0) ||
        !(m.omega_cdm > 0) || !mapped.species.empty() ||
        mapped.omega_massless_nonphoton != 0 ||
        mapped.h0_km_s_mpc != m.h0_km_s_mpc ||
        mapped.omega_gamma != m.omega_gamma || mapped.omega_b != m.omega_b ||
        mapped.omega_cdm != m.omega_cdm)
      return s;
    write(9);
    s.h0 = m.h0_km_s_mpc;
    s.gamma = m.omega_gamma;
    s.baryon = m.omega_b;
    s.cdm = m.omega_cdm;
    s.lambda = witness->lambda_retained;
    s.lambda_getter_loss = witness->lambda_getter_signed_loss;
    s.normalization_estimate = witness->species_normalization_error;
    for (std::size_t j = 0; j < 3; ++j) {
      const auto &w = (*mapping.scalar_witnesses)[j];
      const std::array<double, 3> actual{m.omega_gamma, m.omega_b, m.omega_cdm};
      if (w.emitted_value != actual[j] || !std::isfinite(w.wide_value) ||
          !(w.wide_value > 0) || !std::isfinite(w.wide_operation_estimate) ||
          !(w.wide_operation_estimate >= 0) ||
          !std::isfinite(w.measured_absolute_cast_loss) ||
          !(w.measured_absolute_cast_loss >= 0))
        return Source{};
      write();
      s.map_estimates[j] =
          up_add(w.wide_operation_estimate, w.measured_absolute_cast_loss);
      write(2);
      s.map_wide_values[j] = w.wide_value;
      s.map_operation_estimates[j] = w.wide_operation_estimate;
    }
    if (!std::isfinite(s.lambda) || !(s.lambda >= 0))
      return Source{};
    s.map_witnesses_available = true;
    s.status = Status::ok;
    return s;
  }
};

template <class T> struct Epoch {
  T a, p, h, ell, fg, fb, fc, fl, loading;
};
template <class T> T c_km_s() { return T(299792458L) / T(1000L); }
template <class T> Epoch<T> epoch(const Source &s, const T &a) {
  if (active_budget)
    active_budget->query();
  if (s.status != Status::ok || !a.finite() || !(a > T(0L)))
    throw Stop{Status::outside_domain};
  const T a2 = a * a, a4 = a2 * a2;
  const T radiation(s.gamma), matter = T(s.baryon) + T(s.cdm);
  const T lambda(s.lambda), p = radiation + matter * a + lambda * a4;
  if (!p.finite() || !(p > T(0L)))
    throw Stop{Status::overflow};
  const T h = T(s.h0) / c_km_s<T>() * sqrt(p) / a;
  T loading;
  if (radiation > T(0L))
    loading = T(3L) * T(s.baryon) * a / (T(4L) * radiation);
  else if (!(s.formal_dust && s.gamma == 0 && s.baryon == 0))
    throw Stop{Status::outside_domain};
  Epoch<T> e{a,
             p,
             h,
             (matter * a + T(4L) * lambda * a4) / (T(2L) * p) - T(1L),
             radiation / p,
             T(s.baryon) * a / p,
             T(s.cdm) * a / p,
             lambda * a4 / p,
             loading};
  if (!e.h.finite() || !(e.h > T(0L)))
    throw Stop{Status::overflow};
  return e;
}
// State order a,phi,phi_eta',y,py,delta_c,theta_c.  py=(1+R)y'.
template <class T> using State = std::array<T, 7>;
template <class T>
State<T> rhs(const Source &s, const T &k, const State<T> &v) {
  if (active_budget)
    active_budget->stage();
  const auto e = epoch(s, v[0]);
  const T k2 = k * k, dg = T(4L) * (v[3] + v[1]), h2 = e.h * e.h;
  State<T> f{v[0] * e.h,
             v[2],
             -T(3L) * e.h * v[2] - h2 * (T(2L) * e.ell + T(1L)) * v[1] +
                 h2 * e.fg * dg / T(2L),
             v[4] / (T(1L) + e.loading),
             -k2 * (v[3] + (T(2L) + e.loading) * v[1]) / T(3L),
             -v[6] + T(3L) * v[2],
             -e.h * v[6] + k2 * v[1]};
  for (const auto &x : f)
    if (!x.finite())
      throw Stop{Status::overflow};
  return f;
}

// Four physical directions in the same once-imported source frame.  Gamma,
// baryon and CDM each have dLambda=-dDensity; the fourth is the independently
// owned Lambda-preparation remainder.  These are conditional derivative
// operators, not a source-radius or arithmetic admission by themselves.
enum class SourceParameter : unsigned {
  photon,
  baryon,
  cdm,
  lambda_preparation
};
template <class T> struct ParameterFamily {
  std::array<T, 4> signed_shift, operation_radius, box_radius;
};
template <class T> ParameterFamily<T> parameter_family(const Source &s) {
  if (s.status != Status::ok || !s.map_witnesses_available ||
      s.normalization_estimate != 0)
    throw Stop{Status::conditioning_budget_exceeded};
  ParameterFamily<T> out;
  std::array<T, 3> emitted{T(s.gamma), T(s.baryon), T(s.cdm)};
  const T u = T::unit_roundoff(), one(1L);
  for (unsigned j = 0; j < 3; ++j) {
    const T wide(s.map_wide_values[j]);
    out.signed_shift[j] = wide - emitted[j];
    // Correctly-rounded subtraction amplitude loss is owned once here, with
    // the actual map operation estimate; no measured cast loss is added again.
    const T assembly = upper_multiply(u, upper_add(abs(wide), abs(emitted[j])));
    out.operation_radius[j] =
        upper_add(T(s.map_operation_estimates[j]), assembly);
  }
  std::sort(emitted.begin(), emitted.end(),
            [](const T &a, const T &b) { return a > b; });
  T residue(one), absolute_terms(one);
  for (const auto &density : emitted) {
    residue = residue - density;
    absolute_terms = upper_add(absolute_terms, abs(density));
  }
  out.signed_shift[3] = residue - T(s.lambda);
  absolute_terms = upper_add(absolute_terms, abs(T(s.lambda)));
  const T four_u = upper_multiply(T(4L), u);
  const T denominator = lower_subtract(one, four_u);
  if (!(denominator > T(0L)))
    throw Stop{Status::conditioning_budget_exceeded};
  out.operation_radius[3] =
      upper_multiply(upper_divide(four_u, denominator), absolute_terms);
  for (unsigned j = 0; j < 4; ++j) {
    out.box_radius[j] =
        upper_add(abs(out.signed_shift[j]), out.operation_radius[j]);
    if (!out.box_radius[j].finite() || out.operation_radius[j] < T(0L))
      throw Stop{Status::overflow};
  }
  return out;
}
template <class T> struct PhysicalSourceDirection {
  T h, loading, reciprocal_loading;
};
template <class T>
PhysicalSourceDirection<T> source_direction(const Source &s, const Epoch<T> &e,
                                            SourceParameter parameter) {
  if (!(s.gamma > 0) || !(e.a > T(0L)) || !(e.p > T(0L)))
    throw Stop{Status::outside_domain};
  const T a2 = e.a * e.a, a4 = a2 * a2;
  T pj, rj;
  switch (parameter) {
  case SourceParameter::photon:
    pj = T(1L) - a4;
    rj = -e.loading / T(s.gamma);
    break;
  case SourceParameter::baryon:
    pj = e.a - a4;
    rj = T(3L) * e.a / (T(4L) * T(s.gamma));
    break;
  case SourceParameter::cdm:
    pj = e.a - a4;
    break;
  case SourceParameter::lambda_preparation:
    pj = a4;
    break;
  default:
    throw Stop{Status::invalid_input};
  }
  const T one_plus_r = T(1L) + e.loading;
  return {e.h * pj / (T(2L) * e.p), rj, -rj / (one_plus_r * one_plus_r)};
}
template <class T>
State<T> jacobian_action(const Source &s, const T &k, const State<T> &v,
                         const State<T> &u, const Epoch<T> &e) {
  // COMPLETE seven-state action, including the actual a clock column.
  // The source-control C,D identity is independent of the momentum equation.
  const T a2 = e.a * e.a, a3 = a2 * e.a;
  const T alpha = T(s.h0) / c_km_s<T>(), alpha2 = alpha * alpha;
  const T ha = e.ell * e.h / e.a;
  const T c = T(3L) * alpha2 * (T(s.gamma) / a2 - T(s.lambda) * a2);
  const T d = T(2L) * alpha2 * T(s.gamma) / a2;
  const T ca = -T(6L) * alpha2 * (T(s.gamma) / a3 + T(s.lambda) * e.a);
  const T da = -T(4L) * alpha2 * T(s.gamma) / a3;
  const T ra = e.loading / e.a, one_plus_r = T(1L) + e.loading;
  const T ja = -ra / (one_plus_r * one_plus_r), kk = k * k;
  State<T> out{e.h * (T(1L) + e.ell) * u[0],
               u[2],
               (ca * v[1] - T(3L) * ha * v[2] + da * v[3]) * u[0] + c * u[1] -
                   T(3L) * e.h * u[2] + d * u[3],
               ja * v[4] * u[0] + u[4] / one_plus_r,
               -kk * (ra * v[1] * u[0] + (T(2L) + e.loading) * u[1] + u[3]) /
                   T(3L),
               T(3L) * u[2] - u[6],
               -ha * v[6] * u[0] + kk * u[1] - e.h * u[6]};
  for (const auto &value : out)
    if (!value.finite())
      throw Stop{Status::overflow};
  return out;
}
template <class T>
State<T> parameter_force(const Source &s, const T &k, const State<T> &v,
                         const Epoch<T> &e, SourceParameter parameter) {
  const auto direction = source_direction(s, e, parameter);
  const T alpha = T(s.h0) / c_km_s<T>(), alpha2 = alpha * alpha;
  const T a2 = e.a * e.a;
  T trace;
  switch (parameter) {
  case SourceParameter::photon:
    // The SAME Gamma variation drives C and D.  Preserve 3phi+2y before
    // multiplying by the early a^-2 coefficient.  The explicit a^2 term is
    // its correlated vacuum direction; it must not be rounded away as 1-a^4.
    trace = alpha2 * (T(3L) * v[1] + T(2L) * v[3]) / a2 +
            T(3L) * alpha2 * a2 * v[1];
    break;
  case SourceParameter::baryon:
  case SourceParameter::cdm:
    trace = T(3L) * alpha2 * a2 * v[1];
    break;
  case SourceParameter::lambda_preparation:
    trace = -T(3L) * alpha2 * a2 * v[1];
    break;
  default:
    throw Stop{Status::invalid_input};
  }
  State<T> out{e.a * direction.h,
               T(0L),
               trace - T(3L) * direction.h * v[2],
               direction.reciprocal_loading * v[4],
               -k * k * direction.loading * v[1] / T(3L),
               T(0L),
               -direction.h * v[6]};
  for (const auto &value : out)
    if (!value.finite())
      throw Stop{Status::overflow};
  return out;
}
template <class T>
State<T> sensitivity_rhs(const Source &s, const T &k, const State<T> &v,
                         const State<T> &u, const Epoch<T> &e,
                         SourceParameter parameter) {
  // A future response integration is part of the original attempt, with its
  // own actual RHS owner and scalar charges, not a free post-hoc readout.
  if (active_budget)
    active_budget->stage();
  State<T> out = jacobian_action(s, k, v, u, e);
  const auto source = parameter_force(s, k, v, e, parameter);
  for (std::size_t i = 0; i < out.size(); ++i)
    out[i] = out[i] + source[i];
  return out;
}

// Specific finite Taylor terms for this ideal acoustic polynomial.  All are
// positive MPFR-directed expressions.  The central coefficient/RHS arithmetic
// and finite response integrator are separate errors; this mathematical source
// enclosure alone does not supply their missing admission.
template <class T> struct CoefficientRemainders {
  T h, c, d, loading, reciprocal_loading;
  T h_change, c_change, d_change, loading_change, reciprocal_loading_change;
  T h_lower, h_upper;
};
template <class T>
CoefficientRemainders<T>
coefficient_remainders(const Source &s, const Epoch<T> &e,
                       const ParameterFamily<T> &family,
                       const T &clock_radius) {
  const T zero(0L), one(1L), two(2L), three(3L), four(4L), six(6L);
  if (!clock_radius.finite() || clock_radius < zero)
    throw Stop{Status::overflow};
  const auto &q = family.box_radius;
  const T al = lower_subtract(e.a, clock_radius),
          ah = upper_add(e.a, clock_radius);
  const T gl = lower_subtract(T(s.gamma), q[0]),
          gh = upper_add(T(s.gamma), q[0]);
  const T bl = lower_subtract(T(s.baryon), q[1]),
          bh = upper_add(T(s.baryon), q[1]);
  const T cl = lower_subtract(T(s.cdm), q[2]), ch = upper_add(T(s.cdm), q[2]);
  const T lr = upper_add(upper_add(q[0], q[1]), upper_add(q[2], q[3]));
  const T ll = lower_subtract(T(s.lambda), lr), lh = upper_add(T(s.lambda), lr);
  // This is the mathematical positive source neighborhood needed for clock
  // inversion; it does not admit output requests above the original .01.
  if (!(al > zero) || !(gl > zero) || !(bl > zero) || !(cl > zero) ||
      ll < zero || ah > one)
    throw Stop{Status::conditioning_budget_exceeded};
  const T al2 = lower_multiply(al, al), al3 = lower_multiply(al2, al),
          al4 = lower_multiply(al2, al2);
  const T ah2 = upper_multiply(ah, ah), ah3 = upper_multiply(ah2, ah),
          ah4 = upper_multiply(ah2, ah2);
  const T ml = lower_multiply(lower_add(bl, cl), al);
  const T pl = lower_add(gl, lower_add(ml, lower_multiply(ll, al4)));
  const T ph = upper_add(gh, upper_add(upper_multiply(upper_add(bh, ch), ah),
                                       upper_multiply(lh, ah4)));
  if (!(pl > zero))
    throw Stop{Status::conditioning_budget_exceeded};
  const T alpha_lower =
      lower_divide(lower_multiply(T(s.h0), T(1000L)), T(299792458L));
  const T alpha_upper =
      upper_divide(upper_multiply(T(s.h0), T(1000L)), T(299792458L));
  const T alpha2 = upper_multiply(alpha_upper, alpha_upper);
  CoefficientRemainders<T> out;
  out.h_upper = upper_divide(upper_multiply(alpha_upper, upper_sqrt(ph)), al);
  // sqrt is monotone; lower 1/sqrt(ph) and P_lo ensure a positive H lower.
  out.h_lower = lower_divide(lower_multiply(alpha_lower, pl),
                             upper_multiply(ah, upper_sqrt(ph)));
  if (!(out.h_lower > zero))
    throw Stop{Status::conditioning_budget_exceeded};
  const T pa = upper_add(upper_add(bh, ch),
                         upper_multiply(four, upper_multiply(lh, ah3)));
  const T pdir =
      upper_add(q[0], upper_add(upper_multiply(ah, upper_add(q[1], q[2])),
                                upper_multiply(ah4, q[3])));
  const T b1 = upper_add(upper_multiply(pa, clock_radius), pdir);
  const T aj = upper_add(
      upper_multiply(four, upper_multiply(ah3, upper_add(q[0], q[3]))),
      upper_multiply(upper_add(one, upper_multiply(four, ah3)),
                     upper_add(q[1], q[2])));
  const T aa = upper_multiply(clock_radius, clock_radius);
  const T b2 = upper_add(
      upper_multiply(T(12L), upper_multiply(upper_multiply(lh, ah2), aa)),
      upper_multiply(two, upper_multiply(aj, clock_radius)));
  const T b1p = upper_divide(b1, pl), apl = lower_multiply(al, pl);
  const T hterm = upper_add(
      upper_divide(b2, upper_multiply(two, pl)),
      upper_add(upper_divide(upper_multiply(b1p, b1p), four),
                upper_add(upper_divide(upper_multiply(clock_radius, b1), apl),
                          upper_divide(upper_multiply(two, aa), al2))));
  out.h = upper_multiply(upper_divide(out.h_upper, two), hterm);
  const T cterm = upper_add(
      upper_multiply(upper_add(upper_divide(upper_multiply(six, gh), al4),
                               upper_multiply(two, lh)),
                     aa),
      upper_add(
          upper_divide(upper_multiply(four, upper_multiply(q[0], clock_radius)),
                       al3),
          upper_multiply(
              four, upper_multiply(upper_multiply(ah, lr), clock_radius))));
  out.c =
      upper_multiply(upper_divide(upper_multiply(three, alpha2), two), cterm);
  out.d = upper_multiply(
      alpha2,
      upper_add(
          upper_divide(upper_multiply(six, upper_multiply(gh, aa)), al4),
          upper_divide(upper_multiply(four, upper_multiply(q[0], clock_radius)),
                       al3)));
  const T rfirst =
      upper_add(upper_divide(upper_add(upper_multiply(q[1], ah),
                                       upper_multiply(bh, clock_radius)),
                             gl),
                upper_divide(upper_multiply(upper_multiply(bh, ah), q[0]),
                             lower_multiply(gl, gl)));
  const T rterm = upper_add(
      upper_divide(upper_multiply(q[1], clock_radius), gl),
      upper_add(upper_divide(
                    upper_multiply(upper_add(upper_multiply(q[1], ah),
                                             upper_multiply(bh, clock_radius)),
                                   q[0]),
                    lower_multiply(gl, gl)),
                upper_divide(upper_multiply(upper_multiply(bh, ah),
                                            upper_multiply(q[0], q[0])),
                             lower_multiply(lower_multiply(gl, gl), gl))));
  out.loading = upper_multiply(rat<T>(3, 4), rterm);
  const T rfirst_scaled = upper_multiply(rat<T>(3, 4), rfirst);
  const T rl = lower_divide(lower_multiply(three, lower_multiply(bl, al)),
                            upper_multiply(four, gh));
  const T jrden = lower_add(one, rl);
  out.reciprocal_loading = upper_add(
      upper_divide(upper_multiply(rfirst_scaled, rfirst_scaled),
                   lower_multiply(lower_multiply(jrden, jrden), jrden)),
      upper_divide(out.loading, lower_multiply(jrden, jrden)));
  // First variations are bounded on the same box, including clock cross.
  out.h_change = upper_add(
      upper_add(upper_divide(upper_multiply(out.h_upper, clock_radius), al),
                upper_divide(upper_multiply(out.h_upper, pdir),
                             upper_multiply(two, pl))),
      out.h);
  const T cafirst = upper_multiply(
      six, upper_multiply(alpha2, upper_add(upper_divide(gh, al3),
                                            upper_multiply(lh, ah))));
  const T csource = upper_multiply(
      three, upper_multiply(alpha2, upper_add(upper_divide(q[0], al2),
                                              upper_multiply(lr, ah2))));
  out.c_change = upper_add(
      upper_add(upper_multiply(cafirst, clock_radius), csource), out.c);
  out.d_change = upper_add(
      upper_add(
          upper_divide(upper_multiply(
                           four, upper_multiply(
                                     alpha2, upper_multiply(gh, clock_radius))),
                       al3),
          upper_divide(upper_multiply(two, upper_multiply(alpha2, q[0])), al2)),
      out.d);
  out.loading_change = upper_add(rfirst_scaled, out.loading);
  out.reciprocal_loading_change =
      upper_add(upper_divide(rfirst_scaled, lower_multiply(jrden, jrden)),
                out.reciprocal_loading);
  for (const auto *value :
       {&out.h, &out.c, &out.d, &out.loading, &out.reciprocal_loading,
        &out.h_change, &out.c_change, &out.d_change, &out.loading_change,
        &out.reciprocal_loading_change, &out.h_lower, &out.h_upper})
    if (!value->finite() || *value < zero)
      throw Stop{Status::overflow};
  return out;
}
template <class T>
State<T> source_state_box(const std::array<State<T>, 4> &sensitivity,
                          const State<T> &remainder,
                          const ParameterFamily<T> &family) {
  State<T> box;
  for (std::size_t i = 0; i < box.size(); ++i) {
    if (!remainder[i].finite() || remainder[i] < T(0L))
      throw Stop{Status::conditioning_budget_exceeded};
    box[i] = remainder[i];
    for (unsigned j = 0; j < 4; ++j)
      box[i] = upper_add(
          box[i], upper_multiply(abs(sensitivity[j][i]), family.box_radius[j]));
  }
  return box;
}
template <class T>
State<T> finite_source_remainder_force(const T &k, const State<T> &v,
                                       const State<T> &box,
                                       const CoefficientRemainders<T> &r) {
  const T three(3L), kk = upper_multiply(abs(k), abs(k));
  const T K = upper_divide(kk, three);
  State<T> out;
  out[0] = upper_add(upper_multiply(abs(v[0]), r.h),
                     upper_multiply(box[0], r.h_change));
  out[2] = upper_add(
      upper_multiply(r.c, abs(v[1])),
      upper_add(
          upper_multiply(three, upper_multiply(r.h, abs(v[2]))),
          upper_add(upper_multiply(r.d, abs(v[3])),
                    upper_add(upper_multiply(r.c_change, box[1]),
                              upper_add(upper_multiply(
                                            three,
                                            upper_multiply(r.h_change, box[2])),
                                        upper_multiply(r.d_change, box[3]))))));
  out[3] = upper_add(upper_multiply(r.reciprocal_loading, abs(v[4])),
                     upper_multiply(r.reciprocal_loading_change, box[4]));
  out[4] =
      upper_multiply(K, upper_add(upper_multiply(r.loading, abs(v[1])),
                                  upper_multiply(r.loading_change, box[1])));
  out[6] = upper_add(upper_multiply(r.h, abs(v[6])),
                     upper_multiply(r.h_change, box[6]));
  for (const auto &value : out)
    if (!value.finite())
      throw Stop{Status::overflow};
  return out;
}
template <class T>
State<T> remainder_comparison_rhs(const Source &s, const T &k,
                                  const State<T> &v, const State<T> &radius,
                                  const Epoch<T> &e, const State<T> &force) {
  // Full seven-dimensional comparison, not a constraint norm multiplied by a
  // final endpoint factor.  This central comparison operator still requires
  // its own finite integration/rounding response before admitting a bound.
  if (active_budget)
    active_budget->stage();
  const T alpha = T(s.h0) / c_km_s<T>(), a2 = e.a * e.a, a3 = a2 * e.a,
          alpha2 = alpha * alpha;
  const T ha = e.ell * e.h / e.a;
  const T C = T(3L) * alpha2 * (T(s.gamma) / a2 - T(s.lambda) * a2),
          D = T(2L) * alpha2 * T(s.gamma) / a2;
  const T ca = -T(6L) * alpha2 * (T(s.gamma) / a3 + T(s.lambda) * e.a),
          da = -T(4L) * alpha2 * T(s.gamma) / a3;
  const T ra = e.loading / e.a, J = T(1L) / (T(1L) + e.loading),
          ja = -ra * J * J, K = k * k / T(3L);
  for (std::size_t i = 0; i < radius.size(); ++i)
    if (!radius[i].finite() || radius[i] < T(0L) || !force[i].finite() ||
        force[i] < T(0L))
      throw Stop{Status::conditioning_budget_exceeded};
  State<T> out{e.h * (T(1L) + e.ell) * radius[0] + force[0],
               radius[2] + force[1],
               abs(ca * v[1] - T(3L) * ha * v[2] + da * v[3]) * radius[0] +
                   abs(C) * radius[1] - T(3L) * e.h * radius[2] +
                   abs(D) * radius[3] + force[2],
               abs(ja * v[4]) * radius[0] + J * radius[4] + force[3],
               abs(K * ra * v[1]) * radius[0] +
                   K * (T(2L) + e.loading) * radius[1] + K * radius[3] +
                   force[4],
               T(3L) * radius[2] + radius[6] + force[5],
               abs(ha * v[6]) * radius[0] + k * k * radius[1] -
                   e.h * radius[6] + force[6]};
  for (const auto &value : out)
    if (!value.finite())
      throw Stop{Status::overflow};
  return out;
}
template <class T> struct SourceResponse {
  std::array<State<T>, 4> sensitivity;
  State<T> nonlinear_radius;
};
template <class T>
SourceResponse<T>
source_response_rhs(const Source &s, const T &k, const State<T> &v,
                    const SourceResponse<T> &response, const Epoch<T> &e,
                    const ParameterFamily<T> &family) {
  SourceResponse<T> out;
  for (unsigned j = 0; j < 4; ++j)
    out.sensitivity[j] = sensitivity_rhs(s, k, v, response.sensitivity[j], e,
                                         static_cast<SourceParameter>(j));
  const auto box =
      source_state_box(response.sensitivity, response.nonlinear_radius, family);
  const auto coefficients = coefficient_remainders(s, e, family, box[0]);
  const auto force = finite_source_remainder_force(k, v, box, coefficients);
  out.nonlinear_radius =
      remainder_comparison_rhs(s, k, v, response.nonlinear_radius, e, force);
  return out;
}
template <class T> struct Constraints {
  T momentum, hamiltonian, zeta, momentum_scale, hamiltonian_scale;
};
template <class T>
Constraints<T> constraints(const Source &s, const T &k, const State<T> &v) {
  const auto e = epoch(s, v[0]);
  const T z = v[2] / e.h, dg = T(4L) * (v[3] + v[1]);
  const T db = T(3L) * (v[3] + v[1]), ta = -T(3L) * v[4] / (T(1L) + e.loading);
  const T vc = e.h * v[6] / (k * k), va = e.h * ta / (k * k);
  const T F = e.fc + e.fb + rat<T>(4, 3) * e.fg;
  const T W = e.fc * vc + (e.fb + rat<T>(4, 3) * e.fg) * va;
  const T A = e.fc * v[5] + e.fb * db + e.fg * dg, x2 = k * k / (e.h * e.h);
  return {z + v[1] - rat<T>(3, 2) * W,
          x2 * v[1] + T(3L) * (z + v[1]) + rat<T>(3, 2) * A,
          -v[1] + A / (T(3L) * F),
          abs(z) + abs(v[1]) +
              rat<T>(3, 2) *
                  (e.fc * abs(vc) + (e.fb + rat<T>(4, 3) * e.fg) * abs(va)),
          abs(x2 * v[1]) + T(3L) * abs(z + v[1]) +
              rat<T>(3, 2) *
                  (e.fc * abs(v[5]) + e.fb * abs(db) + e.fg * abs(dg))};
}
template <class T> T normalized(const T &residual, const T &scale) {
  if (!residual.finite() || !scale.finite() || scale < T(0L))
    throw Stop{Status::overflow};
  if (scale.zero()) {
    if (!residual.zero())
      throw Stop{Status::conditioning_budget_exceeded};
    return T(0L);
  }
  return abs(residual) / scale;
}
template <class T> State<T> initial(const Source &s, const T &k, const T &a) {
  const auto e = epoch(s, a);
  const T p0 = -rat<T>(2, 3), m = (T(s.baryon) + T(s.cdm)) * a / T(s.gamma);
  const T m2 = m * m, x2 = k * k / (e.h * e.h);
  const T phi = p0 * (T(1L) - m / T(16L) + rat<T>(7, 160) * m2 - x2 / T(30L));
  const T z = p0 * (-m / T(16L) + rat<T>(7, 80) * m2 - x2 / T(15L));
  const T vc =
      p0 * (rat<T>(1, 2) + m / T(16L) - rat<T>(7, 160) * m2 - x2 / T(120L));
  const T va =
      p0 * (rat<T>(1, 2) + m / T(16L) - rat<T>(7, 160) * m2 - x2 / T(20L));
  const T dc = p0 * (-rat<T>(3, 2) - rat<T>(3, 16) * m + rat<T>(21, 160) * m2 -
                     rat<T>(7, 20) * x2);
  const T dg =
      p0 * (-T(2L) - m / T(4L) + rat<T>(7, 40) * m2 - rat<T>(7, 15) * x2);
  // Independent raw primitive series: no Hamiltonian/momentum projection.
  return {a,
          phi,
          e.h * z,
          dg / T(4L) - phi,
          -(T(1L) + e.loading) * k * k * va / (T(3L) * e.h),
          dc,
          k * k * vc / e.h};
}
template <class T>
State<T> initial_sensitivity(const Source &s, const T &k, const State<T> &raw,
                             const Epoch<T> &e, SourceParameter parameter) {
  // This differentiates the unchanged RAW initial family at fixed ai.  It
  // neither projects Delta nor changes its asymptotic primordial normalization.
  const auto direction = source_direction(s, e, parameter);
  const T m = (T(s.baryon) + T(s.cdm)) * e.a / T(s.gamma);
  const T x2 = k * k / (e.h * e.h), x2j = -T(2L) * x2 * direction.h / e.h;
  T mj;
  switch (parameter) {
  case SourceParameter::photon:
    mj = -m / T(s.gamma);
    break;
  case SourceParameter::baryon:
  case SourceParameter::cdm:
    mj = e.a / T(s.gamma);
    break;
  case SourceParameter::lambda_preparation:
    break;
  default:
    throw Stop{Status::invalid_input};
  }
  const T p0 = -rat<T>(2, 3);
  const T phij =
      p0 * ((-rat<T>(1, 16) + rat<T>(7, 80) * m) * mj - x2j / T(30L));
  const T zj = p0 * ((-rat<T>(1, 16) + rat<T>(7, 40) * m) * mj - x2j / T(15L));
  const T vcj = p0 * ((rat<T>(1, 16) - rat<T>(7, 80) * m) * mj - x2j / T(120L));
  const T vaj = p0 * ((rat<T>(1, 16) - rat<T>(7, 80) * m) * mj - x2j / T(20L));
  const T dcj =
      p0 * ((-rat<T>(3, 16) + rat<T>(21, 80) * m) * mj - rat<T>(7, 20) * x2j);
  const T dgj =
      p0 * ((-rat<T>(1, 4) + rat<T>(7, 20) * m) * mj - rat<T>(7, 15) * x2j);
  const T z = raw[2] / e.h, vc = raw[6] * e.h / (k * k);
  const T va = -T(3L) * e.h * raw[4] / ((T(1L) + e.loading) * k * k);
  State<T> out{T(0L),
               phij,
               direction.h * z + e.h * zj,
               dgj / T(4L) - phij,
               -k * k *
                   (direction.loading * va + (T(1L) + e.loading) * vaj -
                    (T(1L) + e.loading) * va * direction.h / e.h) /
                   (T(3L) * e.h),
               dcj,
               k * k * (vcj - vc * direction.h / e.h) / e.h};
  for (const auto &value : out)
    if (!value.finite())
      throw Stop{Status::overflow};
  return out;
}
template <class T>
SourceResponse<T>
initial_source_response(const Source &s, const T &k, const State<T> &raw,
                        const Epoch<T> &e, const ParameterFamily<T> &family) {
  SourceResponse<T> out;
  for (unsigned j = 0; j < 4; ++j)
    out.sensitivity[j] =
        initial_sensitivity(s, k, raw, e, static_cast<SourceParameter>(j));
  const auto r = coefficient_remainders(s, e, family, T(0L));
  const auto &q = family.box_radius;
  const T gl = lower_subtract(T(s.gamma), q[0]);
  const T Mh =
      upper_add(upper_add(T(s.baryon), T(s.cdm)), upper_add(q[1], q[2]));
  const T mhi = upper_divide(upper_multiply(Mh, e.a), gl);
  const T a2 = upper_multiply(e.a, e.a), a4 = upper_multiply(a2, a2);
  const T mr = upper_add(q[1], q[2]);
  const T Rm = upper_multiply(
      e.a,
      upper_add(upper_divide(upper_multiply(mr, q[0]), lower_multiply(gl, gl)),
                upper_divide(upper_multiply(Mh, upper_multiply(q[0], q[0])),
                             lower_multiply(lower_multiply(gl, gl), gl))));
  const T Dm = upper_add(upper_divide(upper_multiply(e.a, mr), T(s.gamma)),
                         upper_divide(upper_multiply(mhi, q[0]), gl));
  const T Bm = upper_add(Dm, Rm);
  const T Pdir = upper_add(
      q[0], upper_add(upper_multiply(e.a, mr), upper_multiply(a4, q[3])));
  const T Px = upper_divide(Pdir, gl);
  const T x2hi = upper_divide(upper_multiply(abs(k), abs(k)),
                              lower_multiply(r.h_lower, r.h_lower));
  const T Rx2 = upper_multiply(x2hi, upper_multiply(Px, Px));
  const T Dx2 = upper_divide(upper_multiply(x2hi, Pdir), gl);
  const T p0abs = upper_divide(T(2L), T(3L));
  const auto uq = [](long n, long d) { return upper_divide(T(n), T(d)); };
  const auto mode_remainder = [&](const T &c1abs, const T &c2abs,
                                  const T &cxabs) {
    return upper_multiply(
        p0abs,
        upper_add(upper_multiply(
                      upper_add(c1abs, upper_multiply(
                                           T(2L), upper_multiply(c2abs, mhi))),
                      Rm),
                  upper_add(upper_multiply(c2abs, upper_multiply(Bm, Bm)),
                            upper_multiply(cxabs, Rx2))));
  };
  const auto mode_change = [&](const T &c1abs, const T &c2abs, const T &cxabs,
                               const T &remainder) {
    return upper_add(
        upper_multiply(
            p0abs,
            upper_add(
                upper_multiply(
                    upper_add(c1abs, upper_multiply(
                                         T(2L), upper_multiply(c2abs, mhi))),
                    Dm),
                upper_multiply(cxabs, Dx2))),
        remainder);
  };
  const T Rphi = mode_remainder(uq(1, 16), uq(7, 160), uq(1, 30));
  const T Rz = mode_remainder(uq(1, 16), uq(7, 80), uq(1, 15));
  const T Rvc = mode_remainder(uq(1, 16), uq(7, 160), uq(1, 120));
  const T Rva = mode_remainder(uq(1, 16), uq(7, 160), uq(1, 20));
  const T Rdc = mode_remainder(uq(3, 16), uq(21, 160), uq(7, 20));
  const T Rdg = mode_remainder(uq(1, 4), uq(7, 40), uq(7, 15));
  const T Bz = mode_change(uq(1, 16), uq(7, 80), uq(1, 15), Rz);
  const T Bvc = mode_change(uq(1, 16), uq(7, 160), uq(1, 120), Rvc);
  const T Bva = mode_change(uq(1, 16), uq(7, 160), uq(1, 20), Rva);
  const auto mode_magnitude = [&](bool constant, long wave_denominator,
                                  long matter_denominator) {
    T sum = upper_add(upper_multiply(uq(1, 16), mhi),
                      upper_add(upper_multiply(uq(7, matter_denominator),
                                               upper_multiply(mhi, mhi)),
                                upper_multiply(uq(1, wave_denominator), x2hi)));
    if (constant)
      sum = upper_add(uq(1, 2), sum);
    return upper_multiply(p0abs, sum);
  };
  const T zabs = mode_magnitude(false, 15, 80),
          vcabs = mode_magnitude(true, 120, 160);
  const T vaabs = mode_magnitude(true, 20, 160);
  const T inverse = upper_divide(T(1L), r.h_lower);
  const T inverse_change =
      upper_divide(r.h_change, lower_multiply(r.h_lower, r.h_lower));
  const T inverse_remainder = upper_divide(
      upper_add(
          r.h, upper_divide(upper_multiply(r.h_change, r.h_change), r.h_lower)),
      lower_multiply(r.h_lower, r.h_lower));
  out.nonlinear_radius[1] = Rphi;
  out.nonlinear_radius[2] = upper_add(
      upper_multiply(r.h_upper, Rz),
      upper_add(upper_multiply(zabs, r.h), upper_multiply(r.h_change, Bz)));
  out.nonlinear_radius[3] = upper_add(upper_divide(Rdg, T(4L)), Rphi);
  out.nonlinear_radius[5] = Rdc;
  out.nonlinear_radius[6] = upper_multiply(
      upper_multiply(abs(k), abs(k)),
      upper_add(upper_multiply(inverse, Rvc),
                upper_add(upper_multiply(vcabs, inverse_remainder),
                          upper_multiply(Bvc, inverse_change))));
  const T Rhi = upper_divide(
      upper_multiply(T(3L), upper_multiply(upper_add(T(s.baryon), q[1]), e.a)),
      lower_multiply(T(4L), gl));
  const T one_plus_R = upper_add(T(1L), Rhi);
  const T productabs = upper_multiply(one_plus_R, vaabs);
  const T product_remainder =
      upper_add(upper_multiply(one_plus_R, Rva),
                upper_add(upper_multiply(vaabs, r.loading),
                          upper_multiply(r.loading_change, Bva)));
  const T product_change =
      upper_add(upper_multiply(one_plus_R, Bva),
                upper_add(upper_multiply(vaabs, r.loading_change),
                          upper_multiply(r.loading_change, Bva)));
  out.nonlinear_radius[4] = upper_multiply(
      upper_divide(upper_multiply(abs(k), abs(k)), T(3L)),
      upper_add(upper_multiply(inverse, product_remainder),
                upper_add(upper_multiply(productabs, inverse_remainder),
                          upper_multiply(product_change, inverse_change))));
  for (const auto &value : out.nonlinear_radius)
    if (!value.finite() || value < T(0L))
      throw Stop{Status::overflow};
  // Exact raw-initial family Taylor remainder only.  The separate m^3,m*x2,
  // R*x2,wave/Lambda truncation and local initialization arithmetic remain in
  // their original start/arithmetic shares, never relabelled as source zero.
  return out;
}
template <class T>
State<T> fixed_a_sensitivity(const State<T> &flow,
                             const State<T> &sensitivity) {
  // Linear event inversion is not the complete finite endpoint error.  Its
  // positivity interval and nonlinear event remainder are still required.
  if (!flow[0].finite() || !(flow[0] > T(0L)))
    throw Stop{Status::conditioning_budget_exceeded};
  const T delay = sensitivity[0] / flow[0];
  State<T> out;
  for (std::size_t i = 1; i < out.size(); ++i)
    out[i] = sensitivity[i] - flow[i] * delay;
  // The requested a is fixed by definition of this differentiated output.
  // No central physical state or its clock is reset by this response map.
  for (const auto &value : out)
    if (!value.finite())
      throw Stop{Status::overflow};
  return out;
}
template <class T>
std::array<T, 9> outputs(const Source &s, const T &k, const State<T> &v) {
  const auto e = epoch(s, v[0]);
  const T theta = -T(3L) * v[4] / (T(1L) + e.loading),
          temperature = v[3] + v[1];
  return {v[5] + T(3L) * e.h * v[6] / (k * k),
          T(3L) * temperature + T(3L) * e.h * theta / (k * k),
          temperature,
          v[1],
          v[5],
          T(3L) * temperature,
          theta / k,
          v[6] / k,
          v[2] / e.h};
}
template <class T>
std::array<T, 9>
output_sensitivity(const Source &s, const T &k, const State<T> &v,
                   const State<T> &fixed_a_u, const Epoch<T> &e,
                   SourceParameter parameter) {
  if (!fixed_a_u[0].zero())
    throw Stop{Status::invalid_input};
  const auto direction = source_direction(s, e, parameter);
  const T temperature = fixed_a_u[3] + fixed_a_u[1];
  const T j = T(1L) / (T(1L) + e.loading), hj = direction.h;
  const T theta = -T(3L) * j * v[4];
  const T theta_j =
      -T(3L) * (direction.reciprocal_loading * v[4] + j * fixed_a_u[4]);
  std::array<T, 9> out{
      fixed_a_u[5] + T(3L) * (hj * v[6] + e.h * fixed_a_u[6]) / (k * k),
      T(3L) * temperature + T(3L) * (hj * theta + e.h * theta_j) / (k * k),
      temperature,
      fixed_a_u[1],
      fixed_a_u[5],
      T(3L) * temperature,
      theta_j / k,
      fixed_a_u[6] / k,
      fixed_a_u[2] / e.h - v[2] * hj / (e.h * e.h)};
  for (const auto &value : out)
    if (!value.finite())
      throw Stop{Status::overflow};
  return out;
}

template <class T> struct Trial {
  State<T> high, low;
};
template <class T>
Trial<T> trial(const Source &s, const T &k, const State<T> &v, const T &h) {
  std::array<State<T>, 7> f;
  f[0] = rhs(s, k, v);
  const auto stage = [&](std::size_t j, std::initializer_list<T> coefficients) {
    State<T> w(v);
    std::size_t q = 0;
    for (const auto &c : coefficients) {
      for (std::size_t i = 0; i < 7; ++i)
        w[i] = w[i] + h * c * f[q][i];
      ++q;
    }
    f[j] = rhs(s, k, w);
  };
  stage(1, {rat<T>(1, 5)});
  stage(2, {rat<T>(3, 40), rat<T>(9, 40)});
  stage(3, {rat<T>(44, 45), -rat<T>(56, 15), rat<T>(32, 9)});
  stage(4, {rat<T>(19372, 6561), -rat<T>(25360, 2187), rat<T>(64448, 6561),
            -rat<T>(212, 729)});
  stage(5, {rat<T>(9017, 3168), -rat<T>(355, 33), rat<T>(46732, 5247),
            rat<T>(49, 176), -rat<T>(5103, 18656)});
  stage(6, {rat<T>(35, 384), T(0L), rat<T>(500, 1113), rat<T>(125, 192),
            -rat<T>(2187, 6784), rat<T>(11, 84)});
  Trial<T> result{v, v};
  const std::array<T, 7> bh{rat<T>(35, 384),
                            T(0L),
                            rat<T>(500, 1113),
                            rat<T>(125, 192),
                            -rat<T>(2187, 6784),
                            rat<T>(11, 84),
                            T(0L)};
  const std::array<T, 7> bl{rat<T>(5179, 57600),    T(0L),
                            rat<T>(7571, 16695),    rat<T>(393, 640),
                            -rat<T>(92097, 339200), rat<T>(187, 2100),
                            rat<T>(1, 40)};
  for (std::size_t i = 0; i < 7; ++i)
    for (std::size_t j = 0; j < 7; ++j) {
      result.high[i] = result.high[i] + h * bh[j] * f[j][i];
      result.low[i] = result.low[i] + h * bl[j] * f[j][i];
    }
  return result;
}

struct Point {
  double requested_a = 0;
  long double actual_a = 0, eta_elapsed = 0, hcal = 0, endpoint_a_error = 0;
  long double conformal_age_mpc = 0, initial_age_lambda_error_mpc = 0;
  std::array<long double, 9> values{}, emission_losses{};
  long double momentum = 0, hamiltonian = 0;
  long double normalized_momentum = 0, normalized_hamiltonian = 0;
  // These are full-output source/round/endpoint contributions, not a residual
  // norm times an arbitrary endpoint multiplier.  Missing is NOT zero.
  std::optional<std::array<long double, 9>> source_age_endpoint_error;
  std::optional<std::array<long double, 9>> arithmetic_constraint_error;
};
struct MeshEndpoint {
  std::array<long double, 7> state{};
  long double eta_elapsed = 0, momentum = 0, hamiltonian = 0;
  long double normalized_momentum = 0, normalized_hamiltonian = 0;
};
struct Run {
  Status status = Status::invalid_input;
  bool initial_state_available = false;
  unsigned decimal_digits = 0, resolution = 0;
  long double initial_a = 0;
  long double initial_eta_mpc = 0, initial_age_lambda_error_mpc = 0;
  std::array<long double, 7> initial_state{};
  long double initial_momentum = 0, initial_hamiltonian = 0, initial_zeta = 0;
  long double maximum_normalized_momentum = 0,
              maximum_normalized_hamiltonian = 0;
  std::vector<Point> points;
  std::vector<MeshEndpoint> mesh;
};
inline void payload(const Run &r, const Budget &b) {
  // Inline scalar scratch is separately a conservative 256 scalar objects at
  // the larger precision.  Actual MPFR integer temporaries/caches are NOT
  // inferred from this bound and retain their explicit whole-call gate.
  if (b.retained_payload_prefix > b.maximum_payload)
    throw Stop{Status::work_limit};
  const std::size_t fixed = sizeof(Run) + 256 * sizeof(Real<120>);
  const std::size_t available = b.maximum_payload - b.retained_payload_prefix;
  if (fixed > available)
    throw Stop{Status::work_limit};
  if (r.points.capacity() > (available - fixed) / sizeof(Point) ||
      r.mesh.capacity() >
          (available - fixed - r.points.capacity() * sizeof(Point)) /
              sizeof(MeshEndpoint))
    throw Stop{Status::work_limit};
}
template <class T> long double emitted_loss(const T &v, long double emitted) {
  return abs(v - T(emitted)).upward_wide();
}
template <class T>
void record_mesh(const Source &s, const T &k, const State<T> &v, const T &eta,
                 Run &r, Budget &b) {
  b.endpoint();
  if (r.mesh.size() == r.mesh.capacity()) {
    const std::size_t next = std::min(
        b.maximum_endpoints, std::max(std::size_t(16), r.mesh.capacity() * 2));
    if (next <= r.mesh.size())
      throw Stop{Status::work_limit};
    // Refuse before reserve; include old+new blocks during vector growth.
    if (b.retained_payload_prefix > b.maximum_payload)
      throw Stop{Status::work_limit};
    const std::size_t available = b.maximum_payload - b.retained_payload_prefix;
    const std::size_t other = sizeof(Run) + 256 * sizeof(Real<120>) +
                              r.points.capacity() * sizeof(Point);
    if (other > available ||
        r.mesh.capacity() > (available - other) / sizeof(MeshEndpoint))
      throw Stop{Status::work_limit};
    if (next > (available - other) / sizeof(MeshEndpoint) - r.mesh.capacity())
      throw Stop{Status::work_limit};
    if (r.mesh.size() > std::numeric_limits<std::size_t>::max() / 12)
      throw Stop{Status::work_limit};
    // Existing twelve-field records are copied into the new allocation.  This
    // counted prefix is checked before reserve, together with old+new payload.
    write(12 * r.mesh.size());
    r.mesh.reserve(next);
  }
  const auto c = constraints(s, k, v);
  write(12);
  MeshEndpoint w;
  for (std::size_t i = 0; i < 7; ++i) {
    write();
    w.state[i] = v[i].wide();
  }
  write(3);
  w.eta_elapsed = eta.wide();
  w.momentum = c.momentum.wide();
  w.hamiltonian = c.hamiltonian.wide();
  write(2);
  w.normalized_momentum =
      normalized(c.momentum, c.momentum_scale).upward_wide();
  w.normalized_hamiltonian =
      normalized(c.hamiltonian, c.hamiltonian_scale).upward_wide();
  write(2);
  r.maximum_normalized_momentum =
      std::max(r.maximum_normalized_momentum, w.normalized_momentum);
  r.maximum_normalized_hamiltonian =
      std::max(r.maximum_normalized_hamiltonian, w.normalized_hamiltonian);
  write(12);
  r.mesh.push_back(w);
  payload(r, b);
  if (w.normalized_hamiltonian > 1e-6L || w.normalized_momentum > 1e-6L)
    throw Stop{Status::conditioning_budget_exceeded};
}

// Literal fixed-mesh DP54 source draft.  Base eta step=.02/(H+k); its half and
// quarter are the three time resolutions.  No FSAL reuse crosses an endpoint.
// An a-event is bracketed by charged trial bisections, with no a/state reset.
// The final actual-a discrepancy and its full physical output translation must
// be admitted separately before this raw reference can qualify a prediction.
template <unsigned Digits>
void run_into(Run &r, const Source &s, double k_input, long double initial_a,
              std::span<const double> requested_a, unsigned resolution,
              Budget &b) {
  BudgetScope scope(b);
  using T = Real<Digits>;
  // Populate the final owned campaign slot directly.  No intermediate Run
  // floating-header copy/move or assumed NRVO is used by the ten-attempt path.
  r.decimal_digits = Digits;
  r.resolution = resolution;
  ++b.attempts;
  try {
    write(16);
    r.initial_a = initial_a; // named initialized record fields + ai
    if (s.status != Status::ok || !std::isfinite(k_input) || !(k_input > 0) ||
        !(initial_a > 0) || resolution > 2 || requested_a.empty())
      throw Stop{Status::invalid_input};
    for (std::size_t j = 0; j < requested_a.size(); ++j)
      if (!std::isfinite(requested_a[j]) || !(requested_a[j] >= initial_a) ||
          requested_a[j] > .01 || (j && !(requested_a[j] > requested_a[j - 1])))
        throw Stop{Status::outside_domain};
    const std::size_t fixed = sizeof(Run) + 256 * sizeof(T);
    if (b.retained_payload_prefix > b.maximum_payload ||
        fixed > b.maximum_payload - b.retained_payload_prefix ||
        requested_a.size() >
            (b.maximum_payload - b.retained_payload_prefix - fixed) /
                sizeof(Point))
      throw Stop{Status::work_limit};
    r.points.reserve(requested_a.size());
    const T k(k_input), start(initial_a);
    if (k_input < 1e-7 || k_input > .01)
      throw Stop{Status::outside_domain};
    const auto started = epoch(s, start);
    const T m_start = (T(s.baryon) + T(s.cdm)) * start / T(s.gamma);
    if (m_start > rat<T>(1, 10000) || started.loading > rat<T>(1, 10000) ||
        k * k / (started.h * started.h) > rat<T>(1, 1000000) ||
        started.fl > rat<T>(1, 100000) * rat<T>(1, 100000))
      throw Stop{Status::outside_domain};
    State<T> v = initial(s, k, start);
    T eta;
    // Independent SOURCE-control age: radiation+matter analytic age is an
    // upper value.  Nonnegative Lambda has the explicit integral difference
    // below.  Source/map and clock propagation are still required separately.
    const T radiation(s.gamma), matter = T(s.baryon) + T(s.cdm);
    const T eta_i = T(2L) * c_km_s<T>() / T(s.h0) * start /
                    (sqrt(radiation + matter * start) + sqrt(radiation));
    const T a2 = start * start, a4 = a2 * a2;
    const T age_lambda_error =
        c_km_s<T>() * T(s.lambda) * a4 * start /
        (T(10L) * T(s.h0) * radiation * radiation * sqrt(radiation));
    write(2);
    r.initial_eta_mpc = eta_i.wide();
    r.initial_age_lambda_error_mpc = age_lambda_error.upward_wide();
    const auto ci = constraints(s, k, v);
    for (std::size_t j = 0; j < 7; ++j) {
      write();
      r.initial_state[j] = v[j].wide();
    }
    write(3);
    r.initial_momentum = ci.momentum.wide();
    r.initial_hamiltonian = ci.hamiltonian.wide();
    r.initial_zeta = ci.zeta.wide();
    r.initial_state_available = true;
    record_mesh(s, k, v, eta, r, b);
    for (double requested : requested_a) {
      const T target(requested);
      while (v[0] < target) {
        const auto e = epoch(s, v[0]);
        const T base =
            rat<T>(1, 50) / (e.h + k) / T(static_cast<long>(1U << resolution));
        T h = base;
        Trial<T> t = trial(s, k, v, h);
        if (t.high[0] > target) {
          T lo, hi(h);
          bool landed = false;
          // Fixed finite ceiling; all attempted event stages count.
          for (unsigned q = 0; q < 160; ++q) {
            h = (lo + hi) / T(2L);
            t = trial(s, k, v, h);
            const T distance = abs(t.high[0] - target);
            const T small = rat<T>(1, 1000000),
                    event_tolerance = target * small * small * small * small;
            if (distance <= event_tolerance && t.high[0] <= target) {
              landed = true;
              break;
            }
            if (t.high[0] > target)
              hi = h;
            else
              lo = h;
          }
          if (!landed)
            throw Stop{Status::conditioning_budget_exceeded};
        }
        if (!(t.high[0] > v[0]) || !t.high[0].finite())
          throw Stop{Status::overflow};
        v = t.high;
        eta = eta + h;
        record_mesh(s, k, v, eta, r, b);
        // If the actual event is just below the requested scale factor, keep
        // that state and its error; do not take a vanishing unobserved step.
        const T small = rat<T>(1, 1000000);
        if (abs(v[0] - target) <= target * small * small * small * small)
          break;
      }
      // Raw source-control age/phase guard.  The complete source/time/endpoint
      // phase radius remains unavailable together with prediction admission.
      if (k * (eta_i + eta + age_lambda_error) > T(20L))
        throw Stop{Status::outside_domain};
      const auto out = outputs(s, k, v);
      const auto e = epoch(s, v[0]);
      const auto c = constraints(s, k, v);
      write(29);
      Point p;
      write(5);
      p.requested_a = requested;
      p.actual_a = v[0].wide();
      p.eta_elapsed = eta.wide();
      p.hcal = e.h.wide();
      p.endpoint_a_error = abs(v[0] - target).upward_wide();
      write(2);
      p.conformal_age_mpc = (eta_i + eta).wide();
      p.initial_age_lambda_error_mpc = age_lambda_error.upward_wide();
      for (std::size_t j = 0; j < 9; ++j) {
        write();
        p.values[j] = out[j].wide();
        write();
        p.emission_losses[j] = emitted_loss(out[j], p.values[j]);
      }
      write(2);
      p.momentum = c.momentum.wide();
      p.hamiltonian = c.hamiltonian.wide();
      write(2);
      p.normalized_momentum =
          normalized(c.momentum, c.momentum_scale).upward_wide();
      p.normalized_hamiltonian =
          normalized(c.hamiltonian, c.hamiltonian_scale).upward_wide();
      write(29);
      r.points.push_back(p);
      payload(r, b);
    }
    // Central computation completed.  This is a raw control status, not a
    // scientific/source/resource admission; missing optional errors stay
    // absent.
    r.status = Status::ok;
  } catch (const Stop &failure) {
    r.status = failure.status; // first failed prefix/partial endpoints retained
  }
}
template <unsigned Digits>
Run run(const Source &s, double k_input, long double initial_a,
        std::span<const double> requested_a, unsigned resolution, Budget &b) {
  Run result;
  run_into<Digits>(result, s, k_input, initial_a, requested_a, resolution, b);
  return result; // standalone equation control; campaign never uses this move
}

struct Campaign {
  std::array<Run, 10> attempts;
  Status status = Status::invalid_input;
  // Return attempt8: ai/16,120digits,quarter mesh.  Attempt9 is 90digits with
  // exactly the same start/mesh/endpoints.  No smaller start replaces native
  // ai.
  std::vector<std::array<long double, 9>> time_precision, initial_refinement;
  bool source_admitted = false, resource_admitted = false;
};
// Literal outward measurements of stored wide witnesses.  Their arithmetic
// postconditions are checked; they do not supply missing physical source bands.
inline long double up_add(long double a, long double b) {
  if (!std::isfinite(a) || !std::isfinite(b) || !(a >= 0) || !(b >= 0))
    throw Stop{Status::overflow};
  if (a == 0 || b == 0) {
    write();
    return a == 0 ? b : a;
  }
  write();
  const long double raw = a + b;
  if (!std::isfinite(raw) || std::fpclassify(raw) != FP_NORMAL)
    throw Stop{Status::overflow};
  write();
  const long double out =
      std::nextafter(raw, std::numeric_limits<long double>::infinity());
  if (!std::isfinite(out))
    throw Stop{Status::overflow};
  return out;
}
inline long double up_times_eight(long double a) {
  if (!std::isfinite(a) || !(a >= 0))
    throw Stop{Status::overflow};
  if (a == 0) {
    write();
    return 0;
  }
  write();
  const long double raw = 8 * a;
  if (!std::isfinite(raw) || std::fpclassify(raw) != FP_NORMAL)
    throw Stop{Status::overflow};
  write();
  const long double out =
      std::nextafter(raw, std::numeric_limits<long double>::infinity());
  if (!std::isfinite(out))
    throw Stop{Status::overflow};
  return out;
}
inline long double up_difference(long double a, long double b) {
  if (!std::isfinite(a) || !std::isfinite(b))
    throw Stop{Status::overflow};
  if (a == b) {
    write();
    return 0;
  }
  write();
  const long double raw = a - b;
  if (!std::isfinite(raw) || std::fpclassify(raw) != FP_NORMAL)
    throw Stop{Status::overflow};
  write();
  const long double magnitude = std::abs(raw);
  write();
  const long double out =
      std::nextafter(magnitude, std::numeric_limits<long double>::infinity());
  if (!std::isfinite(out))
    throw Stop{Status::overflow};
  return out;
}
inline long double measured_difference(const Point &a, const Point &b,
                                       std::size_t j) {
  return up_add(
      up_add(up_difference(a.values[j], b.values[j]), a.emission_losses[j]),
      b.emission_losses[j]);
}
inline Campaign campaign(const Source &s, double k, long double native_ai,
                         std::span<const double> a, Budget &b) {
  BudgetScope scope(b);
  Campaign c;
  const std::size_t outer_prefix = b.retained_payload_prefix;
  struct PrefixRestore {
    Budget &budget;
    std::size_t prefix;
    ~PrefixRestore() { budget.retained_payload_prefix = prefix; }
  } restore{b, outer_prefix};
  b.rhs_this_k = 0;
  if (!(native_ai > 0) || a.empty())
    return c;
  try {
    if (outer_prefix > b.maximum_payload)
      throw Stop{Status::work_limit};
    for (unsigned start = 0; start < 3; ++start)
      for (unsigned time = 0; time < 3; ++time) {
        const unsigned index = 3 * start + time;
        std::size_t prefix = outer_prefix;
        if (sizeof(Campaign) > b.maximum_payload - prefix)
          throw Stop{Status::work_limit};
        prefix += sizeof(Campaign);
        for (unsigned prior = 0; prior < index; ++prior) {
          const auto &saved = c.attempts[prior];
          if (saved.points.capacity() >
              (b.maximum_payload - prefix) / sizeof(Point))
            throw Stop{Status::work_limit};
          prefix += saved.points.capacity() * sizeof(Point);
          if (saved.mesh.capacity() >
              (b.maximum_payload - prefix) / sizeof(MeshEndpoint))
            throw Stop{Status::work_limit};
          prefix += saved.mesh.capacity() * sizeof(MeshEndpoint);
        }
        b.retained_payload_prefix = prefix;
        const long double ai =
            std::ldexp(native_ai, -static_cast<int>(start + 2));
        write();
        run_into<120>(c.attempts[index], s, k, ai, a, time, b);
        if (c.attempts[index].status != Status::ok) {
          c.status = c.attempts[index].status;
          return c;
        }
      }
    write();
    const long double ai = std::ldexp(native_ai, -4);
    std::size_t prefix = outer_prefix;
    if (prefix > b.maximum_payload ||
        sizeof(Campaign) > b.maximum_payload - prefix)
      throw Stop{Status::work_limit};
    prefix += sizeof(Campaign);
    for (unsigned prior = 0; prior < 9; ++prior) {
      const auto &saved = c.attempts[prior];
      if (saved.points.capacity() >
          (b.maximum_payload - prefix) / sizeof(Point))
        throw Stop{Status::work_limit};
      prefix += saved.points.capacity() * sizeof(Point);
      if (saved.mesh.capacity() >
          (b.maximum_payload - prefix) / sizeof(MeshEndpoint))
        throw Stop{Status::work_limit};
      prefix += saved.mesh.capacity() * sizeof(MeshEndpoint);
    }
    b.retained_payload_prefix = prefix;
    run_into<90>(c.attempts[9], s, k, ai, a, 2, b);
    if (c.attempts[9].status != Status::ok) {
      c.status = c.attempts[9].status;
      return c;
    }
    // All ten retained raw histories, both refinement arrays and the full
    // current scratch are simultaneously owned.  Preflight before resize.
    prefix = outer_prefix;
    const std::size_t fixed = sizeof(Campaign) + 256 * sizeof(Real<120>);
    if (prefix > b.maximum_payload || fixed > b.maximum_payload - prefix)
      throw Stop{Status::work_limit};
    prefix += fixed;
    for (const auto &saved : c.attempts) {
      if (saved.points.capacity() >
          (b.maximum_payload - prefix) / sizeof(Point))
        throw Stop{Status::work_limit};
      prefix += saved.points.capacity() * sizeof(Point);
      if (saved.mesh.capacity() >
          (b.maximum_payload - prefix) / sizeof(MeshEndpoint))
        throw Stop{Status::work_limit};
      prefix += saved.mesh.capacity() * sizeof(MeshEndpoint);
    }
    if (a.size() >
        (b.maximum_payload - prefix) / (2 * sizeof(std::array<long double, 9>)))
      throw Stop{Status::work_limit};
    if (a.size() > b.maximum_writes / 18)
      throw Stop{Status::work_limit};
    write(18 * a.size()); // actual zero initialization of the two error vectors
    c.time_precision.resize(a.size());
    c.initial_refinement.resize(a.size());
    if (c.time_precision.capacity() >
        (b.maximum_payload - prefix) / sizeof(std::array<long double, 9>))
      throw Stop{Status::work_limit};
    prefix += c.time_precision.capacity() * sizeof(std::array<long double, 9>);
    if (c.initial_refinement.capacity() >
        (b.maximum_payload - prefix) / sizeof(std::array<long double, 9>))
      throw Stop{Status::work_limit};
    for (std::size_t i = 0; i < a.size(); ++i)
      for (std::size_t j = 0; j < 9; ++j) {
        const auto &r0 = c.attempts[6].points[i], &r1 = c.attempts[7].points[i];
        const auto &r2 = c.attempts[8].points[i],
                   &r90 = c.attempts[9].points[i];
        write();
        c.time_precision[i][j] =
            up_add(up_times_eight(std::max(measured_difference(r2, r1, j),
                                           measured_difference(r1, r0, j))),
                   measured_difference(r2, r90, j));
        write();
        c.initial_refinement[i][j] = up_times_eight(
            std::max(measured_difference(c.attempts[8].points[i],
                                         c.attempts[5].points[i], j),
                     measured_difference(c.attempts[5].points[i],
                                         c.attempts[2].points[i], j)));
      }
    c.resource_admitted = b.external_arithmetic_payload_admitted;
    c.status = Status::conditioning_budget_exceeded;
    // Full source/constraint/endpoint propagation and the external-runtime
    // payload gate are required before any reference numerical acceptance.
  } catch (const Stop &failure) {
    c.status = failure.status;
  }
  return c;
}
} // namespace ideal_acoustic_peer
