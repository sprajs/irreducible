#pragma once
#include "irred/ideal_acoustic.hpp"
#include "ideal_acoustic_response.hpp"
#include <algorithm>
#include <array>
#include <cfenv>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>

namespace irred::cosmology::detail {
// Private model-specific response, not another nominal H/background owner.
// Y=(Delta,S,dV,Vc,phi,eta); R=(four signed six-vector sensitivities,
// six finite SOURCE remainders, six ARITHMETIC radii). Parent owns RK stages,
// every added derivative-owner charge, prefix records and final admission.
inline constexpr std::size_t ideal_acoustic_response_width = 36;
using IdealAcousticResponseState = std::array<long double, 36>;

struct IdealAcousticTransportAccounting {
  void *context = nullptr;
  bool (*state_writes)(void *, std::size_t) noexcept = nullptr;
  bool (*diagnostic)(void *) noexcept = nullptr;
};
struct IdealAcousticSourceUncertainty {
  numerics::Status status = numerics::Status::invalid_input;
  // Gamma,b,c,Lambda-preparation. Physical density shifts are wide-emitted;
  // the fourth shift is (1-sum(emitted densities))-actual retained Lambda.
  std::array<long double, 4> signed_shift{}, operation_radius{};
  // Original mapper fourth coordinate is retained, NOT our fourth direction.
  bool original_other_zero_witness = false;
  long double original_other_shift = 0, original_other_radius = 0;
};
struct IdealAcousticTransportArithmetic {
  numerics::Status status = numerics::Status::invalid_input;
  // Absolute actual-center minus exact retained smooth-shadow bounds, with
  // all literal center/gradient assembly loss already included by its owner.
  // Order x2,F,B,L,T,cs2. Getter loss belongs here, never another source shift.
  std::array<long double, 6> coefficient_radius{};
  // Same-actual-a arithmetic only. SOURCE quotient/readout derivatives use
  // this band; the full family clock change above is owned by ARITHMETIC.
  std::array<long double, 6> coefficient_at_actual_a_radius{};
  long double hcal_radius = 0, shadow_p_radius = 0, shadow_p_n_radius = 0;
  long double sound_denominator_radius = 0;
  std::array<std::array<long double, 6>, 4> gradient_radius{};
  std::array<long double, 4> eta_gradient_radius{};
  // Additional finite nominal/stage/response defects supplied by the caller.
  // Normal RHS/readout assembly under the explicit profile is added here once;
  // caller must not independently inject that same literal RHS loss again.
  std::array<long double, 6> additional_rhs_radius{};
  bool conditional_wide_clock_profile = false;
};
struct IdealAcousticTransportFrame {
  numerics::Status status = numerics::Status::invalid_input;
  IdealAcousticCoefficients actual{};
  IdealAcousticSourceCenter center{};
  IdealAcousticSourceDirections directions{};
  // Diagnostic centers BORROW the actual H/x. Their nonzero differences to the
  // smooth shadow are in arithmetic; no second H query/sqrt/reclosure occurs.
  long double hcal = 0, x = 0;
  IdealAcousticTransportArithmetic arithmetic;
};
struct IdealAcousticTransportInitialBounds {
  numerics::Status status = numerics::Status::invalid_input;
  // Full smooth-source initial FAMILY difference from its smooth nominal
  // initial map, including Delta-only projection, with correlated parameters.
  // eta can put its entire initial source/early-age error in component5.
  std::array<long double, 6> source_family_radius{};
  // Actual rounded initial map / age estimate / conditional clock seed.
  std::array<long double, 6> arithmetic_radius{};
};

namespace ideal_acoustic_transport_internal {
using W = long double;
using S = numerics::Status;
inline bool normal_or_zero(W x) noexcept {
  return std::isfinite(x) && (x == 0 || std::isnormal(x));
}
inline bool radius(W x) noexcept { return normal_or_zero(x) && x >= 0; }
inline IdealAcousticRadiusFailure record_radius_failure(
    unsigned coordinate, IdealAcousticRadiusChannel channel, W provisional,
    std::optional<W> allowance, W completed, S arithmetic_status) noexcept {
  IdealAcousticRadiusFailure out;
  out.coordinate = static_cast<IdealAcousticStateCoordinate>(coordinate);
  out.channel = channel;
  if (std::isfinite(provisional)) out.provisional_radius = provisional;
  if (allowance && std::isfinite(*allowance))
    out.local_assembly_allowance = *allowance;
  if (std::isfinite(completed)) out.completed_radius = completed;
  out.arithmetic_status = arithmetic_status;
  out.completed_finite = std::isfinite(completed);
  out.completed_normal_or_zero = normal_or_zero(completed);
  out.completed_nonnegative = completed >= 0;
  return out;
}
struct Arithmetic {
  S status = S::ok;
  W checked(W x) noexcept {
    if (!normal_or_zero(x))
      status = S::outside_domain;
    return x;
  }
  W add(W a, W b) noexcept { return checked(a + b); }
  W sub(W a, W b) noexcept { return checked(a - b); }
  W mul(W a, W b) noexcept {
    const W r = a * b;
    if (a != 0 && b != 0 && r == 0)
      status = S::outside_domain;
    return checked(r);
  }
  W div(W a, W b) noexcept {
    if (b == 0) {
      status = S::conditioning_budget_exceeded;
      return 0;
    }
    const W r = a / b;
    if (a != 0 && r == 0)
      status = S::outside_domain;
    return checked(r);
  }
  W magnitude(W a) noexcept { return checked(std::abs(a)); }
  W up(W a) noexcept {
    if (!radius(a)) {
      status = S::outside_domain;
      return 0;
    }
    if (a == 0)
      return 0;
    return checked(std::nextafter(a, std::numeric_limits<W>::infinity()));
  }
  W signed_upper(W x) noexcept {
    if (!normal_or_zero(x)) {
      status = S::outside_domain;
      return 0;
    }
    if (x == 0)
      return 0;
    return checked(std::nextafter(x, std::numeric_limits<W>::infinity()));
  }
  W plus(W a, W b) noexcept {
    if (!radius(a) || !radius(b)) {
      status = S::outside_domain;
      return 0;
    }
    if (a == 0 || b == 0)
      return a == 0 ? b : a;
    return up(add(a, b));
  }
  // A signed RK provisional is not a radius until its once-owned assembly
  // allowance has been added. Complete that signed sum outward; the caller
  // still owns the final nonnegative-radius guard and all prior status checks.
  W assemble_signed_radius(W provisional, W allowance) noexcept {
    if (!radius(allowance)) {
      status = S::outside_domain;
      return 0;
    }
    return signed_upper(add(provisional, allowance));
  }
  W times(W a, W b) noexcept {
    if (!radius(a) || !radius(b)) {
      status = S::outside_domain;
      return 0;
    }
    if (a == 0 || b == 0)
      return 0;
    return up(mul(a, b));
  }
  W over(W a, W b) noexcept {
    if (!radius(a) || !(b > 0) || !normal_or_zero(b)) {
      status = S::conditioning_budget_exceeded;
      return 0;
    }
    if (a == 0)
      return 0;
    return up(div(a, b));
  }
  W lower(W a, W b) noexcept {
    if (!(a > 0) || !normal_or_zero(a) || !radius(b)) {
      status = S::conditioning_budget_exceeded;
      return 0;
    }
    const W r = std::nextafter(sub(a, b), -std::numeric_limits<W>::infinity());
    if (!(r > 0) || !normal_or_zero(r))
      status = S::conditioning_budget_exceeded;
    return r;
  }
  W product_lower(W a, W b) noexcept {
    const W r = std::nextafter(mul(a, b), 0.L);
    if (!(a > 0 && b > 0 && r > 0) || !normal_or_zero(r))
      status = S::conditioning_budget_exceeded;
    return r;
  }
  W ratio_lower(W a, W b) noexcept {
    const W r = std::nextafter(div(a, b), 0.L);
    if (!(a > 0 && b > 0 && r > 0) || !normal_or_zero(r))
      status = S::conditioning_budget_exceeded;
    return r;
  }
  // Conditional selected normal Wide profile, including literal arithmetic
  // and sqrt/log/exp when owned by that profile. Not a universal libm bound.
  W loss(W positive_term_scale) noexcept {
    return times(64 * std::numeric_limits<W>::epsilon(), positive_term_scale);
  }
};
inline S account(const IdealAcousticTransportAccounting &a,
                 std::size_t writes = 0) noexcept {
  if (!a.state_writes || !a.diagnostic)
    return S::invalid_input;
  if (!a.diagnostic(a.context) ||
      (writes && !a.state_writes(a.context, writes)))
    return S::work_limit;
  return S::ok;
}
inline bool valid(const IdealAcousticTransportFrame &f,
                  const IdealAcousticSourceUncertainty &u) noexcept {
  if (f.status != S::ok || u.status != S::ok || f.arithmetic.status != S::ok ||
      f.directions.status != S::ok ||
      !f.arithmetic.conditional_wide_clock_profile ||
      !u.original_other_zero_witness || u.original_other_shift != 0 ||
      u.original_other_radius != 0 || std::fegetround() != FE_TONEAREST ||
      std::numeric_limits<W>::digits < 64 ||
      std::numeric_limits<W>::max_exponent < 16384)
    return false;
  if (!(f.center.a > 0 && f.center.a <= W(.01)) || !(f.center.photon > 0) ||
      !(f.center.baryon > 0) || !(f.center.cdm > 0) ||
      !(f.center.shadow_p > 0) || !(f.center.shadow_p_n >= 0) || !(f.hcal > 0))
    return false;
  for (W v :
       {f.center.a, f.center.photon, f.center.baryon, f.center.cdm,
        f.center.shadow_p, f.center.shadow_p_n, f.hcal, f.actual.x2, f.actual.L,
        f.actual.F, f.actual.B, f.actual.loading_over_one_plus_loading,
        f.actual.sound_speed_squared, f.actual.acceleration_defect})
    if (!normal_or_zero(v))
      return false;
  if (!(f.actual.x2 > 0 && f.actual.F > 0 && f.actual.B > 0 &&
        f.actual.loading_over_one_plus_loading > 0 &&
        f.actual.loading_over_one_plus_loading < 1 &&
        f.actual.sound_speed_squared > 0))
    return false;
  // SourceCenter's coefficient slots are borrowed ACTUAL diagnostic centers,
  // with their nonzero exact-shadow discrepancy in the arithmetic bridge.
  if (f.center.x2 != f.actual.x2 || f.center.F != f.actual.F ||
      f.center.B != f.actual.B || f.center.L != f.actual.L ||
      f.center.loading_fraction != f.actual.loading_over_one_plus_loading ||
      f.center.sound_speed_squared != f.actual.sound_speed_squared)
    return false;
  for (W v : f.arithmetic.coefficient_radius)
    if (!(v > 0) || !radius(v))
      return false;
  for (std::size_t i=0;i<6;++i)
    if (!(f.arithmetic.coefficient_at_actual_a_radius[i]>0) ||
        !radius(f.arithmetic.coefficient_at_actual_a_radius[i]) ||
        f.arithmetic.coefficient_at_actual_a_radius[i]>
            f.arithmetic.coefficient_radius[i])
      return false;
  for (W v :
       {f.arithmetic.hcal_radius, f.arithmetic.shadow_p_radius,
        f.arithmetic.shadow_p_n_radius, f.arithmetic.sound_denominator_radius})
    if (!(v > 0) || !radius(v))
      return false;
  for (std::size_t j = 0; j < 4; ++j) {
    if (!normal_or_zero(u.signed_shift[j]) || !(u.operation_radius[j] > 0) ||
        !radius(u.operation_radius[j]) ||
        !(f.arithmetic.eta_gradient_radius[j] > 0) ||
        !radius(f.arithmetic.eta_gradient_radius[j]))
      return false;
    const auto &g = f.directions.gradient[j];
    const std::array<W, 6> values{
        g.x2, g.F, g.B, g.L, g.loading_fraction, g.sound_speed_squared};
    for (std::size_t k = 0; k < 6; ++k) {
      if (!normal_or_zero(values[k]) ||
          !radius(f.arithmetic.gradient_radius[j][k]))
        return false;
      if ((k < 4 || j < 2) && !(f.arithmetic.gradient_radius[j][k] > 0))
        return false;
    }
  }
  for (W v : f.arithmetic.additional_rhs_radius)
    if (!radius(v))
      return false;
  return true;
}
// Fixed-a finite source quotient bands. These are diagnostic/radius assembly,
// charged as such, not a generic all-floating-operation state-write counter.
struct Bands {
  std::array<W, 4> amplitude, pj;
  std::array<W, 6> change, remainder;
  W inverse_h, inverse_h_arithmetic, eta_remainder;
};
inline S bands(const IdealAcousticTransportFrame &f,
               const IdealAcousticSourceUncertainty &u, Bands &o,
               Arithmetic &a) noexcept {
  const W a2 = a.mul(f.center.a, f.center.a), a4 = a.mul(a2, a2);
  o.pj = {a.sub(1, a4), a.sub(f.center.a, a4), a.sub(f.center.a, a4), a4};
  W dp = 0, dn_f = 0, dn_b = 0, dn_q = 0, dd = 0, dn_t = 0, dn_cs = 0;
  const std::array<W, 4> pnj{-4 * a4, f.center.a - 4 * a4, f.center.a - 4 * a4,
                             4 * a4};
  for (std::size_t j = 0; j < 4; ++j) {
    o.amplitude[j] =
        a.plus(a.magnitude(u.signed_shift[j]), u.operation_radius[j]);
    // Own the fresh a4/pj assembly, even if a4 is lost in represented 1-a4.
    const W pj_bound =
        a.plus(a.magnitude(o.pj[j]),
               a.loss(a.plus(j == 0 ? 1 : (j == 3 ? 0 : f.center.a), a4)));
    dp = a.plus(dp, a.times(pj_bound, o.amplitude[j]));
    dn_q = a.plus(dn_q, a.times(a.plus(a.magnitude(pnj[j]),
                                       a.loss(a.plus(f.center.a, 4 * a4))),
                                o.amplitude[j]));
  }
  const W rg = o.amplitude[0], rb = o.amplitude[1], rc = o.amplitude[2];
  if (!(f.center.photon > rg) || !(f.center.baryon > rb) ||
      !(f.center.cdm > rc))
    return S::conditioning_budget_exceeded;
  dn_f =
      a.plus(a.times(a.up(4.L / 3), rg), a.times(f.center.a, a.plus(rb, rc)));
  dn_b = a.plus(a.times(a.up(4.L / 3), rg), a.times(f.center.a, rb));
  dn_q = a.over(dn_q, 2);
  dd = a.plus(a.times(4, rg), a.times(a.times(3, f.center.a), rb));
  dn_t = a.times(a.times(3, f.center.a), rb);
  dn_cs = a.times(a.up(4.L / 3), rg);
  const W d = a.add(a.mul(4, f.center.photon),
                    a.mul(a.mul(3, f.center.baryon), f.center.a));
  const W plo = a.lower(f.center.shadow_p, f.arithmetic.shadow_p_radius);
  const W dlo = a.lower(d, f.arithmetic.sound_denominator_radius);
  const W p_family_lo = a.lower(plo, dp), d_family_lo = a.lower(dlo, dd);
  if (a.status != S::ok)
    return a.status;
  const std::array<W, 6> denominators{plo, plo, plo, plo, dlo, dlo};
  const std::array<W, 6> changes{dp, dp, dp, dp, dd, dd};
  const std::array<W, 6> numerators{0, dn_f, dn_b, dn_q, dn_t, dn_cs};
  const std::array<W, 6> centers{f.actual.x2,
                                 f.actual.F,
                                 f.actual.B,
                                 f.actual.L + 1,
                                 f.actual.loading_over_one_plus_loading,
                                 f.actual.sound_speed_squared};
  for (std::size_t j = 0; j < 6; ++j) {
    const W fhi =
        a.plus(a.magnitude(centers[j]), f.arithmetic.coefficient_at_actual_a_radius[j]);
    const W numerator = a.plus(numerators[j], a.times(fhi, changes[j]));
    const W lower = j < 4 ? p_family_lo : d_family_lo;
    o.change[j] = a.over(numerator, lower);
    o.remainder[j] = a.times(a.over(changes[j], denominators[j]), o.change[j]);
  }
  o.inverse_h = a.div(1, f.hcal);
  const W hlo = a.lower(f.hcal, f.arithmetic.hcal_radius);
  o.inverse_h_arithmetic =
      a.plus(a.over(f.arithmetic.hcal_radius, a.product_lower(f.hcal, hlo)),
             a.loss(o.inverse_h));
  // P^-1/2 Taylor remainder, without creating another nominal H/sqrt.
  const W ih_hi = a.over(
      a.times(a.plus(o.inverse_h, o.inverse_h_arithmetic), plo), p_family_lo);
  const W t = a.over(dp, p_family_lo);
  o.eta_remainder = a.times(a.up(3.L / 8), a.times(ih_hi, a.times(t, t)));
  return a.status;
}
inline std::array<W, 6> rhs_scale(std::span<const W, 6> y,
                                  const IdealAcousticCoefficients &c,
                                  W inverse_h, Arithmetic &a) noexcept {
  const W dv = a.magnitude(y[2]), vc = a.magnitude(y[3]),
          phi = a.magnitude(y[4]);
  return {
      a.plus(
          a.times(c.x2, vc),
          a.plus(a.times(a.times(4.5L, c.B), dv),
                 a.times(a.times(3, a.magnitude(c.acceleration_defect)), vc))),
      a.times(c.x2, dv),
      a.plus(a.times(a.plus(a.magnitude(c.L), c.loading_over_one_plus_loading),
                     dv),
             a.times(c.sound_speed_squared,
                     a.plus(a.magnitude(y[0]), a.magnitude(y[1])))),
      a.plus(a.times(a.plus(a.magnitude(c.L), 1), vc), phi),
      a.plus(phi, a.plus(a.times(a.times(1.5L, c.F), vc),
                         a.times(a.times(1.5L, c.B), dv))),
      inverse_h};
}
inline std::array<W, 6> matrix_force(const std::array<W, 6> &coefficient,
                                     std::span<const W, 6> y,
                                     Arithmetic &a) noexcept {
  return {a.plus(a.times(coefficient[0], y[3]),
                 a.times(a.times(4.5L, coefficient[2]), y[2])),
          a.times(coefficient[0], y[2]),
          a.plus(a.times(a.plus(coefficient[3], coefficient[4]), y[2]),
                 a.times(coefficient[5], a.plus(y[0], y[1]))),
          a.times(coefficient[3], y[3]),
          a.plus(a.times(a.times(1.5L, coefficient[1]), y[3]),
                 a.times(a.times(1.5L, coefficient[2]), y[2])),
          0};
}
inline S response_sizes(std::span<const W, 36> r,
                        const IdealAcousticSourceUncertainty &u,
                        std::array<W, 6> &linear, std::array<W, 6> &total,
                        Arithmetic &a) noexcept {
  for (std::size_t i = 0; i < 6; ++i) {
    if (!radius(r[24 + i]) || !radius(r[30 + i]))
      return S::conditioning_budget_exceeded;
    W signed_sum = 0, uncertain = 0, term_scale = 0;
    for (std::size_t j = 0; j < 4; ++j) {
      if (!normal_or_zero(r[6 * j + i]))
        return S::outside_domain;
      signed_sum = a.add(signed_sum, a.mul(u.signed_shift[j], r[6 * j + i]));
      uncertain = a.plus(
          uncertain, a.times(u.operation_radius[j], a.magnitude(r[6 * j + i])));
      term_scale = a.plus(term_scale, a.times(a.magnitude(u.signed_shift[j]),
                                              a.magnitude(r[6 * j + i])));
    }
    linear[i] =
        a.plus(a.magnitude(signed_sum), a.plus(uncertain, a.loss(term_scale)));
    total[i] = a.plus(linear[i], a.plus(r[24 + i], r[30 + i]));
  }
  return a.status;
}
} // namespace ideal_acoustic_transport_internal

// Definitions below use caller-owned buffers. Successful status is availability
// of these response estimates only; full prediction/campaign admission is the
// caller's distinct fixed shares, refinement, source and clock/output gates.
// Central/signed RHS only. Fresh literal assembly leaves are returned as
// tau, then injected at their actual discrete RK4 sites by the local shadow.
// Inherited coefficient/gradient/clock uncertainty stays in the radius forcing.
inline numerics::Status ideal_acoustic_transport_rhs(
    std::span<const long double, 6> y, std::span<const long double, 24> z,
    const IdealAcousticTransportFrame &f,
    const IdealAcousticSourceUncertainty &u, std::span<long double, 24> dz,
    std::span<long double, 6> tau,
    const IdealAcousticTransportAccounting &work) noexcept {
  using namespace ideal_acoustic_transport_internal;
  if (!valid(f, u)) return S::conditioning_budget_exceeded;
  for (W v : y) if (!normal_or_zero(v)) return S::outside_domain;
  for (W v : z) if (!normal_or_zero(v)) return S::outside_domain;
  S status = account(work, 6);
  if (status != S::ok) return status;
  Arithmetic a;
  Bands b;
  status = bands(f, u, b, a);
  if (status != S::ok) return status;
  std::array<W, 6> absolute_y;
  for (std::size_t i = 0; i < 6; ++i) absolute_y[i] = a.magnitude(y[i]);
  const auto nominal_scale = rhs_scale(y, f.actual, b.inverse_h, a);
  for (std::size_t i = 0; i < 6; ++i) tau[i] = a.loss(nominal_scale[i]);
  for (std::size_t j = 0; j < 4; ++j) {
    const std::span<const W, 6> zj(z.data() + 6*j, 6);
    if (!work.state_writes(work.context, 6)) return S::work_limit;
    // Keep the separately rounded original derivative/source expressions and
    // their final addition. Neither intermediate five-state vector is stored.
    // Scalar expression scratch retains the existing diagnostic ownership.
    const auto assign = [&]<unsigned I>() noexcept {
      dz[6*j+I] = a.add(ideal_acoustic_derivative_coordinate<I>(
                            std::span<const W,5>(zj.data(),5), f.actual),
                        ideal_acoustic_source_force_coordinate<I>(
                            std::span<const W,5>(y.data(),5),
                            f.directions.gradient[j]));
    };
    assign.template operator()<0>();
    assign.template operator()<1>();
    assign.template operator()<2>();
    assign.template operator()<3>();
    assign.template operator()<4>();
    dz[6*j+5] = -a.div(a.mul(b.inverse_h, b.pj[j]),
                        a.mul(2, f.center.shadow_p));
    const auto sensitivity_scale = rhs_scale(zj, f.actual, 0, a);
    const auto &g = f.directions.gradient[j];
    const std::array<W, 6> absolute_gradient{
        a.magnitude(g.x2), a.magnitude(g.F), a.magnitude(g.B),
        a.magnitude(g.L), a.magnitude(g.loading_fraction),
        a.magnitude(g.sound_speed_squared)};
    auto force_scale = matrix_force(absolute_gradient, absolute_y, a);
    force_scale[5] = a.magnitude(dz[6*j+5]);
    if (!work.state_writes(work.context, 6)) return S::work_limit;
    for (std::size_t i = 0; i < 6; ++i)
      tau[i] = a.plus(tau[i], a.times(b.amplitude[j],
          a.loss(a.plus(sensitivity_scale[i], force_scale[i]))));
  }
  for (W v : dz) if (!normal_or_zero(v)) return S::outside_domain;
  for (W v : tau) if (!radius(v)) return S::outside_domain;
  return a.status;
}

namespace ideal_acoustic_transport_internal {
// Kplus is evaluated directly without subtracting a negative drift from a
// completed derivative. positive[] is the separately earned diagonal bound.
inline std::array<W,6> positive_comparison(std::span<const W,6> r,
    const IdealAcousticCoefficients &c, std::span<const W,6> positive,
    Arithmetic &a) noexcept {
  const W signed_coefficient = a.add(-c.x2, a.mul(3,c.acceleration_defect));
  const W coefficient = a.plus(a.magnitude(signed_coefficient),
      a.loss(a.plus(c.x2,a.times(3,a.magnitude(c.acceleration_defect)))));
  return {
      a.plus(a.times(positive[0],r[0]),a.plus(a.times(coefficient,r[3]),
                                            a.times(a.times(4.5L,c.B),r[2]))),
      a.plus(a.times(positive[1],r[1]),a.times(c.x2,r[2])),
      a.plus(a.times(positive[2],r[2]),
             a.times(c.sound_speed_squared,a.plus(r[0],r[1]))),
      a.plus(a.times(positive[3],r[3]),r[4]),
      a.plus(a.times(positive[4],r[4]),
             a.plus(a.times(a.times(1.5L,c.F),r[3]),
                    a.times(a.times(1.5L,c.B),r[2]))),
      a.times(positive[5],r[5])};
}
inline std::array<W,6> arithmetic_matrix_force(
    const IdealAcousticTransportFrame &f,std::span<const W,6> r,
    Arithmetic &a) noexcept {
  auto out=matrix_force(f.arithmetic.coefficient_radius,r,a);
  out[0]=a.plus(out[0],a.times(a.times(3,a.magnitude(f.actual.acceleration_defect)),r[3]));
  return out;
}
} // namespace ideal_acoustic_transport_internal

// Complete nonnegative P in R'=D R+P, with the original inherited forcing and
// measurement law. Fresh central RHS/combine pulses are owned separately.
inline numerics::Status ideal_acoustic_transport_forcing(
    std::span<const long double,6> y,std::span<const long double,36> r,
    const IdealAcousticTransportFrame &f,const IdealAcousticSourceUncertainty &u,
    std::span<const long double,6> positive,std::span<long double,12> out,
    const IdealAcousticTransportAccounting &work) noexcept {
  using namespace ideal_acoustic_transport_internal;
  if (!valid(f,u)) return S::conditioning_budget_exceeded;
  for (W v:y) if (!normal_or_zero(v)) return S::outside_domain;
  for (W v:positive) if (!radius(v)) return S::outside_domain;
  S status=account(work,12);
  if (status!=S::ok) return status;
  Arithmetic a; Bands b;
  status=bands(f,u,b,a);
  if (status!=S::ok) return status;
  std::array<W,6> linear,total,absolute_y,full;
  status=response_sizes(r,u,linear,total,a);
  if (status!=S::ok) return status;
  for (std::size_t i=0;i<6;++i) {
    absolute_y[i]=a.magnitude(y[i]);
    full[i]=a.plus(absolute_y[i],total[i]);
  }
  const auto change=matrix_force(b.change,total,a);
  const auto remainder=matrix_force(b.remainder,absolute_y,a);
  auto arithmetic=arithmetic_matrix_force(f,full,a);
  arithmetic[5]=b.inverse_h_arithmetic;
  std::array<W,6> gradient{};
  for (std::size_t j=0;j<4;++j) {
    const auto inherited=matrix_force(f.arithmetic.gradient_radius[j],absolute_y,a);
    for (std::size_t i=0;i<6;++i)
      gradient[i]=a.plus(gradient[i],a.times(b.amplitude[j],
          i==5 ? f.arithmetic.eta_gradient_radius[j] : inherited[i]));
  }
  const auto source=positive_comparison(
      std::span<const W,6>(r.data()+24,6),f.actual,positive,a);
  const auto rounding=positive_comparison(
      std::span<const W,6>(r.data()+30,6),f.actual,positive,a);
  const auto source_scale=rhs_scale(std::span<const W,6>(r.data()+24,6),f.actual,0,a);
  const auto rounding_scale=rhs_scale(std::span<const W,6>(r.data()+30,6),f.actual,0,a);
  for (std::size_t i=0;i<6;++i) {
    const W finite_source=i==5 ? b.eta_remainder : a.plus(change[i],remainder[i]);
    out[i]=a.plus(source[i],a.plus(finite_source,a.loss(source_scale[i])));
    out[6+i]=a.plus(rounding[i],a.plus(arithmetic[i],a.plus(gradient[i],
        a.plus(f.arithmetic.additional_rhs_radius[i],a.loss(rounding_scale[i])))));
  }
  for (W v:out) if (!radius(v)) return S::outside_domain;
  return a.status;
}

// Coupled12 local shadow: e bounds the SUM of primitive central/signed pulse
// error and additional arithmetic-response variation. Its SOURCE cross is
// present at every actual intermediate RK impulse stage. This positive action
// encloses an exact-real conditional scale law; it is not rounded-map Lipschitz.
inline numerics::Status ideal_acoustic_transport_local_action(
    std::span<const long double,12> v,const IdealAcousticTransportFrame &f,
    const IdealAcousticSourceUncertainty &u,std::span<long double,12> out,
    const IdealAcousticTransportAccounting &work) noexcept {
  using namespace ideal_acoustic_transport_internal;
  if (!valid(f,u)) return S::conditioning_budget_exceeded;
  for (W x:v) if (!radius(x)) return S::outside_domain;
  S status=account(work,12);
  if (status!=S::ok) return status;
  Arithmetic a; Bands b;
  status=bands(f,u,b,a);
  if (status!=S::ok) return status;
  const std::span<const W,6> c(v.data(),6),e(v.data()+6,6);
  const auto sc=rhs_scale(c,f.actual,0,a),se=rhs_scale(e,f.actual,0,a);
  const auto xc=matrix_force(b.change,c,a),xe=matrix_force(b.change,e,a);
  const auto re=matrix_force(b.remainder,e,a);
  const auto ac=arithmetic_matrix_force(f,c,a),ae=arithmetic_matrix_force(f,e,a);
  auto je=se;
  std::array<W,6> ge{};
  for (std::size_t j=0;j<4;++j) {
    const auto &g=f.directions.gradient[j];
    const std::array<W,6> absolute_gradient{
        a.magnitude(g.x2),a.magnitude(g.F),a.magnitude(g.B),a.magnitude(g.L),
        a.magnitude(g.loading_fraction),a.magnitude(g.sound_speed_squared)};
    const auto central=matrix_force(absolute_gradient,e,a);
    const auto inherited=matrix_force(f.arithmetic.gradient_radius[j],e,a);
    for (std::size_t i=0;i<6;++i) {
      je[i]=a.plus(je[i],a.times(b.amplitude[j],central[i]));
      ge[i]=a.plus(ge[i],a.times(b.amplitude[j],inherited[i]));
    }
  }
  const W xi=64*std::numeric_limits<W>::epsilon();
  const W one_xi=a.plus(1,xi),two_xi=a.plus(2,xi);
  for (std::size_t i=0;i<6;++i) {
    out[i]=a.plus(a.plus(sc[i],xc[i]),
        a.plus(a.loss(sc[i]),a.plus(a.times(one_xi,xe[i]),re[i])));
    out[6+i]=a.plus(je[i],a.plus(ac[i],a.plus(a.times(two_xi,ae[i]),
        a.plus(ge[i],a.loss(se[i])))));
  }
  for (W x:out) if (!radius(x)) return S::outside_domain;
  return a.status;
}

// Narrow finite PC endpoint-map bridge, distinct from the RK4 pulse gain.
// q is PURE local Q_E, not flow+Q. Conservative overlap is retained rather
// than assuming equality. The fixed borrowed frame has no Y/Z clock feedback.
inline numerics::Status ideal_acoustic_transport_endpoint_bridge(
    std::span<const long double,6> q,const IdealAcousticTransportFrame &f,
    const IdealAcousticSourceUncertainty &u,long double h,
    std::span<const long double,6> mu_lower,std::span<long double,12> out,
    const IdealAcousticTransportAccounting &work) noexcept {
  using namespace ideal_acoustic_transport_internal;
  if (!valid(f,u)) return S::conditioning_budget_exceeded;
  if (!(h>0) || !radius(h)) return S::invalid_input;
  for (W v:q) if (!radius(v)) return S::outside_domain;
  for (W v:mu_lower) if (!radius(v)) return S::outside_domain;
  S status=account(work,12);
  if (status!=S::ok) return status;
  Arithmetic a; Bands b;
  status=bands(f,u,b,a);
  if (status!=S::ok) return status;
  const auto x=matrix_force(b.change,q,a),r=matrix_force(b.remainder,q,a);
  const auto arithmetic=arithmetic_matrix_force(f,q,a);
  std::array<W,6> gradient{};
  for (std::size_t j=0;j<4;++j) {
    const auto g=matrix_force(f.arithmetic.gradient_radius[j],q,a);
    for (std::size_t i=0;i<6;++i)
      gradient[i]=a.plus(gradient[i],a.times(b.amplitude[j],g[i]));
  }
  const W xi=64*std::numeric_limits<W>::epsilon();
  const W one_xi=a.plus(1,xi),two_xi=a.plus(2,xi),half_h=a.div(h,2);
  for (std::size_t i=0;i<6;++i) {
    W denominator=1;
    if (mu_lower[i]!=0)
      denominator=a.checked(std::nextafter(
          a.add(1,a.product_lower(half_h,mu_lower[i])),0.L));
    if (!(denominator>0)) return S::conditioning_budget_exceeded;
    out[i]=a.over(a.times(half_h,a.plus(a.times(one_xi,x[i]),r[i])),denominator);
    out[6+i]=a.over(a.times(half_h,
        a.plus(a.times(two_xi,arithmetic[i]),gradient[i])),denominator);
  }
  for (W v:out) if (!radius(v)) return S::outside_domain;
  return a.status;
}

inline numerics::Status ideal_acoustic_transport_initial(
    std::span<const long double, 6> y, const IdealAcousticTransportFrame &f,
    const IdealAcousticSourceUncertainty &u,
    const IdealAcousticTransportInitialBounds &initial,
    std::span<long double, 36> r,
    const IdealAcousticTransportAccounting &work) noexcept {
  using namespace ideal_acoustic_transport_internal;
  if (!valid(f, u) || initial.status != S::ok)
    return S::conditioning_budget_exceeded;
  for (W v : y)
    if (!normal_or_zero(v))
      return S::outside_domain;
  for (W v : initial.source_family_radius)
    if (!radius(v))
      return S::conditioning_budget_exceeded;
  for (W v : initial.arithmetic_radius)
    if (!radius(v))
      return S::conditioning_budget_exceeded;
  for (std::size_t i = 0; i < 6; ++i)
    if (y[i] != 0 && !(initial.arithmetic_radius[i] > 0))
      return S::conditioning_budget_exceeded;
  S status = account(work);
  if (status != S::ok)
    return status;
  Arithmetic a;
  Bands b;
  status = bands(f, u, b, a);
  if (status != S::ok)
    return status;
  const W p0 = a.div(-2, 3),
          m = a.div(a.mul(a.add(f.center.baryon, f.center.cdm), f.center.a),
                    f.center.photon);
  const W x2 = f.actual.x2;
  const W S0 = -a.div(a.mul(p0, a.mul(x2, x2)), 96),
          dV0 = -a.div(a.mul(p0, x2), 24);
  const std::array<W, 4> mj{-a.div(m, f.center.photon),
                            a.div(f.center.a, f.center.photon),
                            a.div(f.center.a, f.center.photon), 0};
  if (a.status != S::ok)
    return a.status;
  if (!work.state_writes(work.context, 36))
    return S::work_limit;
  for (std::size_t j = 0; j < 4; ++j) {
    const auto &g = f.directions.gradient[j];
    const W phi_j =
        a.mul(p0, a.sub(a.mul(a.add(-1.L / 16, a.mul(7.L / 80, m)), mj[j]),
                        a.div(g.x2, 30)));
    const W vc_j =
        a.mul(p0, a.sub(a.mul(a.sub(1.L / 16, a.mul(7.L / 80, m)), mj[j]),
                        a.div(g.x2, 120)));
    const W S_j = -a.div(a.mul(p0, a.mul(x2, g.x2)), 48);
    const W dV_j = -a.div(a.mul(p0, g.x2), 24);
    // Derivative of the ACTUAL declared Delta-only projection; no reset of
    // Vc/phi/S/dV, finite-zeta normalization or native momentum initializer.
    const W numerator =
        a.sub(a.add(-a.mul(2.L / 3, a.add(a.mul(g.x2, y[4]), a.mul(x2, phi_j))),
                    a.add(a.mul(g.B, S0), a.mul(f.actual.B, S_j))),
              a.add(a.mul(3, a.add(a.mul(g.B, dV0), a.mul(f.actual.B, dV_j))),
                    a.mul(y[0], g.F)));
    r[6 * j] = a.div(numerator, f.actual.F);
    r[6 * j + 1] = S_j;
    r[6 * j + 2] = dV_j;
    r[6 * j + 3] = vc_j;
    r[6 * j + 4] = phi_j;
    r[6 * j + 5] = 0; // Entire initial eta source uncertainty is in remainder5.
  }
  for (std::size_t i = 0; i < 6; ++i) {
    W linear_bound = 0;
    for (std::size_t j = 0; j < 4; ++j)
      linear_bound = a.plus(linear_bound,
                            a.times(b.amplitude[j], a.magnitude(r[6 * j + i])));
    // A conservative full-family initial enclosure plus |linear| is an
    // explicit finite remainder, not a claim that sensitivities alone suffice.
    r[24 + i] = a.plus(initial.source_family_radius[i], linear_bound);
    r[30 + i] = initial.arithmetic_radius[i];
  }
  return a.status;
}

// Final fixed-a translation of the exact accumulated RK path. The caller's
// log-distance includes its initial/stage N-update lineage and target-log
// allowance; the borrowed frame covers that entire positive-family clock
// neighborhood. Stored exp/log equality alone does not establish zero distance.
// A finite integral inequality gives |deltaY|_inf <= d*K*|Y|_inf/(1-d*K).
// Eta has its separate inhomogeneous inverse-H bound and units. This owns no
// new H/query/log/exp and changes only the caller-owned ARITHMETIC radii.
inline numerics::Status ideal_acoustic_transport_fixed_a(
    std::span<const long double, 6> y,
    std::span<long double, 36> r,
    const IdealAcousticTransportFrame &f,
    const IdealAcousticSourceUncertainty &u, long double log_distance,
    const IdealAcousticTransportAccounting &work) noexcept {
  using namespace ideal_acoustic_transport_internal;
  if (!valid(f,u) || !(log_distance>0) || !radius(log_distance))
    return S::conditioning_budget_exceeded;
  const S charged=account(work);
  if (charged!=S::ok) return charged;
  Arithmetic a;
  Bands b;
  S status=bands(f,u,b,a);
  if (status!=S::ok) return status;
  std::array<W,6> linear,total,upper;
  status=response_sizes(r,u,linear,total,a);
  if (status!=S::ok) return status;
  for (unsigned i=0;i<6;++i) {
    if (!normal_or_zero(y[i])) return S::outside_domain;
    upper[i]=a.plus(a.magnitude(y[i]),total[i]);
  }
  std::array<W,6> coefficient;
  const std::array<W,6> center{f.actual.x2,f.actual.F,f.actual.B,
      a.magnitude(f.actual.L),f.actual.loading_over_one_plus_loading,
      f.actual.sound_speed_squared};
  for (unsigned i=0;i<6;++i)
    coefficient[i]=a.plus(center[i],a.plus(b.change[i],f.arithmetic.coefficient_radius[i]));
  const auto [q,F,B,L,T,cs]=coefficient;
  const std::array<W,5> row_sum{
      a.plus(q,a.plus(a.times(4.5L,B),a.times(3,a.magnitude(f.actual.acceleration_defect)))),
      q,a.plus(a.plus(L,T),a.times(2,cs)),a.plus(L,2),
      a.plus(1,a.times(1.5L,a.plus(F,B)))};
  const W K=*std::max_element(row_sum.begin(),row_sum.end());
  const W magnitude=*std::max_element(upper.begin(),upper.begin()+5);
  const W growth=a.times(log_distance,K);
  const W physical=a.over(a.times(growth,magnitude),a.lower(1,growth));
  const W a2=a.times(f.center.a,f.center.a),a4=a.times(a2,a2);
  const W vacuum=a.plus(a.plus(b.amplitude[0],b.amplitude[1]),
                          a.plus(b.amplitude[2],b.amplitude[3]));
  const W dp=a.plus(b.amplitude[0],a.plus(
      a.times(f.center.a,a.plus(b.amplitude[1],b.amplitude[2])),a.times(a4,vacuum)));
  const W phi=a.plus(f.center.shadow_p,f.arithmetic.shadow_p_radius);
  const W plo=a.lower(f.center.shadow_p,f.arithmetic.shadow_p_radius);
  const W inverse_h_family=a.times(a.plus(b.inverse_h,b.inverse_h_arithmetic),
                                  a.over(phi,a.lower(plo,dp)));
  const W eta=a.over(a.times(log_distance,inverse_h_family),a.lower(1,log_distance));
  if (a.status!=S::ok) return a.status;
  if (!work.state_writes(work.context,6)) return S::work_limit;
  for (unsigned i=0;i<5;++i) r[30+i]=a.plus(r[30+i],physical);
  r[35]=a.plus(r[35],eta);
  return a.status;
}

inline numerics::Status ideal_acoustic_transport_readout(
    std::span<const long double, 6> y, std::span<const long double, 36> r,
    const IdealAcousticTransportFrame &f,
    const IdealAcousticSourceUncertainty &u, std::span<long double, 10> source,
    std::span<long double, 10> arithmetic,
    const IdealAcousticTransportAccounting &work) noexcept {
  using namespace ideal_acoustic_transport_internal;
  if (!valid(f, u))
    return S::conditioning_budget_exceeded;
  // Only readout requires the actual sqrt(x2) already evaluated by its owner.
  // RHS/initialization do not request or interpret a new x/sqrt witness.
  if (!(f.x > 0) || !normal_or_zero(f.x))
    return S::conditioning_budget_exceeded;
  for (W v : y)
    if (!normal_or_zero(v))
      return S::outside_domain;
  S status = account(work);
  if (status != S::ok)
    return status;
  Arithmetic a;
  Bands b;
  status = bands(f, u, b, a);
  if (status != S::ok)
    return status;
  std::array<W, 6> linear, total, absolute_y, full;
  status = response_sizes(r, u, linear, total, a);
  if (status != S::ok)
    return status;
  for (std::size_t i = 0; i < 6; ++i) {
    absolute_y[i] = a.magnitude(y[i]);
    full[i] = a.plus(absolute_y[i], total[i]);
  }
  // Literal squared-identity postcheck of the ONCE borrowed native x witness.
  // Its own normal multiplication/subtraction allowance is arithmetic, not
  // another source coefficient, and needs no independent sqrt certificate.
  const W x_squared = a.mul(f.x, f.x);
  const W squared_loss = a.plus(a.magnitude(a.sub(x_squared, f.actual.x2)),
                                a.loss(a.plus(x_squared, f.actual.x2)));
  const W x2_arithmetic =
      a.plus(f.arithmetic.coefficient_radius[0], squared_loss);
  const W x2_actual_a_arithmetic =
      a.plus(f.arithmetic.coefficient_at_actual_a_radius[0], squared_loss);
  const W x2lo = a.lower(f.actual.x2, x2_actual_a_arithmetic);
  const W x2_full_lo = a.lower(f.actual.x2,x2_arithmetic);
  const W x2familylo = a.lower(x2lo, b.change[0]);
  if (a.status != S::ok)
    return a.status;
  // Exact-real harmonic lower: 2t/(x+t/x)<=sqrt(t), for positive t,x.
  // Outward rational evaluation supplies lower bounds without another sqrt.
  const W xlo = a.ratio_lower(a.product_lower(2, x2lo),
                              a.plus(f.x, a.over(x2lo, f.x))),
          x_full_lo = a.ratio_lower(a.product_lower(2,x2_full_lo),
                                    a.plus(f.x,a.over(x2_full_lo,f.x))),
          xfamilylo = a.ratio_lower(a.product_lower(2, x2familylo),
                                    a.plus(f.x, a.over(x2familylo, f.x)));
  if (!(xlo > 0 && x_full_lo>0 && xfamilylo > 0) || !normal_or_zero(xlo) ||
      !normal_or_zero(x_full_lo) ||
      !normal_or_zero(xfamilylo))
    return S::conditioning_budget_exceeded;
  const W x_arithmetic = a.plus(
      a.over(x2_arithmetic,
             a.lower(a.add(f.x, x_full_lo), a.loss(a.add(f.x, x_full_lo)))),
      a.loss(f.x));
  const W x_actual_a_arithmetic = a.plus(
      a.over(x2_actual_a_arithmetic,
             a.lower(a.add(f.x,xlo),a.loss(a.add(f.x,xlo)))),a.loss(f.x));
  const W x_change = a.over(b.change[0], a.product_lower(2, xfamilylo));
  const W x_remainder = a.plus(
      a.over(b.remainder[0], a.product_lower(2, xfamilylo)),
      a.over(a.times(b.change[0], b.change[0]),
             a.product_lower(
                 8, a.product_lower(xfamilylo,
                                    a.product_lower(xfamilylo, xfamilylo)))));
  const auto linear_map = [&](std::span<const W, 6> v) {
    return std::array<W, 10>{
        v[0],
        a.plus(a.plus(v[0], v[1]), a.times(3, v[2])),
        a.plus(a.over(a.plus(v[0], v[1]), 3), v[3]),
        v[4],
        a.plus(v[0], a.times(3, v[3])),
        a.plus(a.plus(v[0], a.times(3, v[3])), v[1]),
        a.times(f.x, a.plus(v[3], v[2])),
        a.times(f.x, v[3]),
        a.plus(v[4], a.plus(a.times(a.times(1.5L, f.actual.F), v[3]),
                            a.times(a.times(1.5L, f.actual.B), v[2]))),
        v[5]};
  };
  const auto finite_source =
      linear_map(std::span<const W, 6>(r.data() + 24, 6));
  const auto finite_arithmetic =
      linear_map(std::span<const W, 6>(r.data() + 30, 6));
  const auto literal_scale = linear_map(absolute_y);
  std::array<W, 10> signed_sum{}, uncertain{}, sum_scale{}, gradient_loss{};
  // Literal signed first-order readout, including coefficient derivatives.
  // No |sum| replacement by sum| | for the actual mapper shifts.
  for (std::size_t j = 0; j < 4; ++j) {
    const auto &g = f.directions.gradient[j];
    const W *z = r.data() + 6 * j;
    const W xj = a.div(g.x2, a.mul(2, f.x));
    if (!work.state_writes(work.context, 10))
      return S::work_limit;
    const std::array<W, 10> zj{
        z[0],
        a.add(a.sub(z[0], z[1]), a.mul(3, z[2])),
        a.sub(a.div(a.sub(z[0], z[1]), 3), z[3]),
        z[4],
        a.sub(z[0], a.mul(3, z[3])),
        a.sub(a.sub(z[0], a.mul(3, z[3])), z[1]),
        a.add(a.mul(f.x, a.add(z[3], z[2])), a.mul(xj, a.add(y[3], y[2]))),
        a.add(a.mul(f.x, z[3]), a.mul(xj, y[3])),
        a.add(-z[4],
              a.mul(1.5L, a.add(a.add(a.mul(f.actual.F, z[3]),
                                      a.mul(f.actual.B, z[2])),
                                a.add(a.mul(g.F, y[3]), a.mul(g.B, y[2]))))),
        z[5]};
    std::array<W, 6> az;
    for (std::size_t i = 0; i < 6; ++i)
      az[i] = a.magnitude(z[i]);
    auto scale = linear_map(az);
    scale[6] = a.plus(scale[6], a.times(a.magnitude(xj),
                                        a.plus(absolute_y[3], absolute_y[2])));
    scale[7] = a.plus(scale[7], a.times(a.magnitude(xj), absolute_y[3]));
    scale[8] =
        a.plus(scale[8],
               a.times(1.5L, a.plus(a.times(a.magnitude(g.F), absolute_y[3]),
                                    a.times(a.magnitude(g.B), absolute_y[2]))));
    const W xj_error = a.plus(
        a.over(f.arithmetic.gradient_radius[j][0], a.product_lower(2, xlo)),
        a.plus(a.over(a.times(a.magnitude(g.x2), x_actual_a_arithmetic),
                      a.product_lower(2, a.product_lower(f.x, xlo))),
               a.loss(a.magnitude(xj))));
    for (std::size_t i = 0; i < 10; ++i) {
      signed_sum[i] = a.add(signed_sum[i], a.mul(u.signed_shift[j], zj[i]));
      uncertain[i] = a.plus(uncertain[i],
                            a.times(u.operation_radius[j], a.magnitude(zj[i])));
      sum_scale[i] =
          a.plus(sum_scale[i],
                 a.times(a.magnitude(u.signed_shift[j]), a.magnitude(zj[i])));
      gradient_loss[i] =
          a.plus(gradient_loss[i], a.times(b.amplitude[j], a.loss(scale[i])));
    }
    gradient_loss[6] = a.plus(
        gradient_loss[6],
        a.times(b.amplitude[j],
                a.times(xj_error, a.plus(absolute_y[3], absolute_y[2]))));
    gradient_loss[7] =
        a.plus(gradient_loss[7],
               a.times(b.amplitude[j], a.times(xj_error, absolute_y[3])));
    gradient_loss[8] =
        a.plus(gradient_loss[8],
               a.times(a.times(1.5L, b.amplitude[j]),
                       a.plus(a.times(f.arithmetic.gradient_radius[j][1],
                                      absolute_y[3]),
                              a.times(f.arithmetic.gradient_radius[j][2],
                                      absolute_y[2]))));
  }
  if (a.status != S::ok)
    return a.status;
  if (!work.state_writes(work.context, 20))
    return S::work_limit;
  for (std::size_t i = 0; i < 10; ++i) {
    source[i] = a.plus(a.plus(a.magnitude(signed_sum[i]), uncertain[i]),
                       finite_source[i]);
    arithmetic[i] =
        a.plus(finite_arithmetic[i],
               a.plus(gradient_loss[i],
                      a.plus(a.loss(sum_scale[i]), a.loss(literal_scale[i]))));
  }
  // Full finite coefficient/state cross terms. Physical source coefficient
  // change times ARITHMETIC state error remains SOURCE, matching the RHS split.
  const W velocity = a.plus(absolute_y[3], absolute_y[2]),
          velocity_error = a.plus(total[3], total[2]);
  if (!work.state_writes(work.context, 6))
    return S::work_limit;
  source[6] = a.plus(source[6], a.plus(a.times(x_remainder, velocity),
                                       a.times(x_change, velocity_error)));
  source[7] = a.plus(source[7], a.plus(a.times(x_remainder, absolute_y[3]),
                                       a.times(x_change, total[3])));
  source[8] = a.plus(
      source[8],
      a.times(1.5L, a.plus(a.plus(a.times(b.remainder[1], absolute_y[3]),
                                  a.times(b.remainder[2], absolute_y[2])),
                           a.plus(a.times(b.change[1], total[3]),
                                  a.times(b.change[2], total[2])))));
  arithmetic[6] =
      a.plus(arithmetic[6], a.times(x_arithmetic, a.plus(full[3], full[2])));
  arithmetic[7] = a.plus(arithmetic[7], a.times(x_arithmetic, full[3]));
  arithmetic[8] = a.plus(
      arithmetic[8],
      a.times(1.5L,
              a.plus(a.times(f.arithmetic.coefficient_radius[1], full[3]),
                     a.times(f.arithmetic.coefficient_radius[2], full[2]))));
  // Six second writes above are separately precharged before mutation.
  for (W v : source)
    if (!radius(v))
      return S::outside_domain;
  for (W v : arithmetic)
    if (!radius(v))
      return S::outside_domain;
  return a.status;
}
} // namespace irred::cosmology::detail
