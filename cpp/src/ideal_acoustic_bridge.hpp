#pragma once
#include "ideal_acoustic_transport.hpp"
#include "thermal_conformal_epoch.hpp"
#include "thermal_ruler.hpp"

namespace irred::cosmology::detail {
struct IdealAcousticBridgeClock {
  numerics::Status status = numerics::Status::invalid_input;
  long double actual_N = 0;
  // Includes the selected initial log profile and ALL actual N updates/stage
  // assembly. The separate fixed-requested-a endpoint displacement is not here.
  long double accumulated_N_absolute_radius = 0;
  bool conditional_profile = false;
};

namespace ideal_acoustic_bridge_internal {
using W = long double;
using S = numerics::Status;
using Arithmetic = ideal_acoustic_transport_internal::Arithmetic;
// Only this model's diagnostic primitive graph. No nominal coefficients,
// background query, source map, derivative owner or runtime allocator lives
// here.
struct Ball {
  W center = 0, radius = 0;
};
struct Graph {
  Arithmetic a;
  Ball add(Ball x, Ball y) noexcept {
    return {
        a.add(x.center, y.center),
        a.plus(a.plus(x.radius, y.radius),
               a.loss(a.plus(a.magnitude(x.center), a.magnitude(y.center))))};
  }
  Ball subtract(Ball x, Ball y) noexcept {
    return {
        a.sub(x.center, y.center),
        a.plus(a.plus(x.radius, y.radius),
               a.loss(a.plus(a.magnitude(x.center), a.magnitude(y.center))))};
  }
  Ball multiply(Ball x, Ball y) noexcept {
    return {a.mul(x.center, y.center),
            a.plus(a.plus(a.times(a.magnitude(x.center), y.radius),
                          a.times(a.magnitude(y.center), x.radius)),
                   a.plus(a.times(x.radius, y.radius),
                          a.loss(a.times(a.magnitude(x.center),
                                         a.magnitude(y.center)))))};
  }
  Ball divide(Ball n, Ball d) noexcept {
    const W lower = a.lower(d.center, d.radius);
    if (a.status != S::ok)
      return {};
    const W q = a.div(n.center, d.center), product = a.mul(q, d.center);
    const W residual =
        a.plus(a.magnitude(a.sub(n.center, product)),
               a.loss(a.plus(a.magnitude(n.center), a.magnitude(product))));
    return {q, a.plus(a.over(a.plus(n.radius,
                                    a.plus(a.times(a.magnitude(q), d.radius),
                                           residual)),
                             lower),
                      a.loss(a.magnitude(q)))};
  }
  W difference(W actual, Ball exact) noexcept {
    // Measured difference plus BOTH finite graph and subtraction ancestry.
    return a.plus(
        exact.radius,
        a.plus(a.magnitude(a.sub(actual, exact.center)),
               a.loss(a.plus(a.magnitude(actual), a.magnitude(exact.center)))));
  }
  W magnitude(Ball x) noexcept {
    return a.plus(a.magnitude(x.center), x.radius);
  }
};
inline bool positive(W x) noexcept {
  return x > 0 && ideal_acoustic_transport_internal::normal_or_zero(x);
}
} // namespace ideal_acoustic_bridge_internal

// Borrow the already charged actual stage, loading/sound and original primitive
// frame. This only assembles diagnostics. SOURCE mapper shifts/op radii are NOT
// imported as another physical force. They are used ONLY to enclose the
// arithmetic a-clock change across the entire positive physical source family.
// SOURCE/gradient target is the actual borrowed stage a; ARITHMETIC clock owns
// the full family change to the exact N-scale, including physical×clock cross.
// Exact stored primitive inputs define the retained smooth shadow. No second
// H, x, sqrt, exp, log, getter, capture, remap or background query is
// performed. The explicit selected normal Wide log/exp law is conditional, not
// a universal libm certificate. Parent owns finite RK response, full
// seeds/early age, fixed-a endpoint translation, original counters and final
// admission.
inline numerics::Status ideal_acoustic_transport_frame(
    const ThermalConformalEpoch &epoch, const ThermalBaryonLoading &loading,
    const ThermalBaryonSound &sound,
    const IdealAcousticCoefficients &borrowed_actual,
    const IdealAcousticSourceCenter &borrowed_center,
    const IdealAcousticSourceDirections &borrowed_directions,
    long double retained_lambda,
    const IdealAcousticSourceUncertainty &uncertainty,
    const IdealAcousticBridgeClock &clock, IdealAcousticTransportFrame &out,
    const IdealAcousticTransportAccounting &work) noexcept {
  using namespace ideal_acoustic_bridge_internal;
  const auto charged = ideal_acoustic_transport_internal::account(work);
  if (charged != S::ok)
    return charged;
  // Small diagnostic metadata copies allow these borrowed slots to come from
  // out itself. They belong to this charged diagnostic, not physical vectors.
  const auto actual = borrowed_actual;
  const auto center = borrowed_center;
  const auto directions = borrowed_directions;
  out = {};
  if (epoch.status != S::ok || !epoch.shadow || !epoch.forward ||
      loading.status != S::ok || sound.status != S::ok ||
      directions.status != S::ok || clock.status != S::ok ||
      uncertainty.status != S::ok || !uncertainty.original_other_zero_witness ||
      uncertainty.original_other_shift != 0 ||
      uncertainty.original_other_radius != 0 || !clock.conditional_profile ||
      std::fegetround() != FE_TONEAREST ||
      std::numeric_limits<W>::digits < 64 ||
      std::numeric_limits<W>::max_exponent < 16384)
    return S::conditioning_budget_exceeded;
  const auto &shadow = *epoch.shadow;
  const auto normal = ideal_acoustic_transport_internal::normal_or_zero;
  const auto radius = ideal_acoustic_transport_internal::radius;
  if (!epoch.raw_scaled_query || epoch.raw_scaled_query->status != S::ok ||
      epoch.raw_scaled_query->a4_e2 != epoch.p || epoch.fr != epoch.fg ||
      shadow.source_law_id != "mapped-photon-cold-matter-retained-polynomial/"
                              "conditional-diagnostic/v1" ||
      !(center.a > 0 && center.a <= W(.01)) || !positive(center.photon) ||
      !positive(center.baryon) || !positive(center.cdm) ||
      !(retained_lambda >= 0) || !normal(retained_lambda) ||
      !positive(epoch.hcal) || !positive(epoch.x2) ||
      !positive(shadow.actual_hcal_arithmetic_estimate) ||
      !positive(shadow.actual_x2_arithmetic_estimate) ||
      !positive(clock.accumulated_N_absolute_radius) ||
      !normal(clock.actual_N) || !positive(loading.ratio_today) ||
      !positive(loading.arithmetic_estimate) || !positive(loading.numerator) ||
      !positive(loading.denominator) || !positive(sound.loading) ||
      !positive(sound.loading_estimate) || !positive(sound.denominator) ||
      !positive(sound.sound_speed_squared_over_c_squared) ||
      !positive(sound.sound_speed_squared_estimate_over_c_squared) ||
      center.shadow_p != shadow.p_shadow ||
      center.shadow_p_n != shadow.p_shadow_n || center.x2 != actual.x2 ||
      center.F != actual.F || center.B != actual.B || center.L != actual.L ||
      center.loading_fraction != actual.loading_over_one_plus_loading ||
      center.sound_speed_squared != actual.sound_speed_squared ||
      actual.x2 != epoch.x2 || actual.L != epoch.g ||
      actual.F != epoch.enthalpy_fraction ||
      actual.B != epoch.fb + 4 * epoch.fg / 3 ||
      actual.acceleration_defect != epoch.acceleration_defect ||
      actual.loading_over_one_plus_loading !=
          sound.loading / sound.denominator ||
      actual.sound_speed_squared != sound.sound_speed_squared_over_c_squared ||
      sound.loading != loading.ratio_today * center.a ||
      sound.denominator != 1 + sound.loading ||
      loading.numerator != 3 * center.baryon ||
      loading.denominator != 4 * center.photon)
    return S::conditioning_budget_exceeded;
  for (W v : {center.a, center.shadow_p, center.shadow_p_n, actual.x2, actual.F,
              actual.B, actual.L, actual.loading_over_one_plus_loading,
              actual.sound_speed_squared, actual.acceleration_defect})
    if (!normal(v))
      return S::outside_domain;
  if (!positive(center.shadow_p) || !(center.shadow_p_n >= 0))
    return S::conditioning_budget_exceeded;

  Graph g;
  auto &a = g.a;
  std::array<W, 4> amplitude{};
  for (unsigned j = 0; j < 4; ++j) {
    if (!normal(uncertainty.signed_shift[j]) ||
        !positive(uncertainty.operation_radius[j]))
      return S::conditioning_budget_exceeded;
    amplitude[j] = a.plus(a.magnitude(uncertainty.signed_shift[j]),
                          uncertainty.operation_radius[j]);
  }
  const W vacuum_amplitude = a.plus(a.plus(amplitude[0], amplitude[1]),
                                    a.plus(amplitude[2], amplitude[3]));
  if (!(a.lower(center.photon, amplitude[0]) > 0) ||
      !(a.lower(center.baryon, amplitude[1]) > 0) ||
      !(a.lower(center.cdm, amplitude[2]) > 0) ||
      !(a.lower(retained_lambda, vacuum_amplitude) > 0) || a.status != S::ok)
    return S::conditioning_budget_exceeded;
  const W eps_profile = 64 * std::numeric_limits<W>::epsilon();
  const W selected_log_floor =
      a.times(eps_profile, a.plus(1, a.magnitude(clock.actual_N)));
  if (clock.accumulated_N_absolute_radius < selected_log_floor ||
      clock.accumulated_N_absolute_radius >= .25L)
    return S::conditioning_budget_exceeded;
  // exp(t)<=1/(1-|t|), |t|<1; actual exp has the selected relative 64eps
  // error. This finite ratio band owns their cross term without another exp.
  const W nlo = a.lower(1, clock.accumulated_N_absolute_radius);
  const W elo = a.lower(1, eps_profile);
  const W relative_a =
      a.over(a.plus(clock.accumulated_N_absolute_radius, eps_profile),
             a.product_lower(nlo, elo));
  if (!(relative_a > 0 && relative_a < .25L) || !radius(relative_a))
    return S::conditioning_budget_exceeded;
  const W scale_lower = a.lower(1, relative_a);
  const W absolute_a = a.times(center.a, relative_a);
  if (!(a.lower(center.a, absolute_a) > 0))
    return S::conditioning_budget_exceeded;
  // The positive diagnostic neighborhood alone may extend above .01. Its
  // complete finite bound remains within a<=1, without admitting new outputs.
  if (a.plus(center.a, absolute_a) > 1)
    return S::conditioning_budget_exceeded;
  const Ball scale{center.a, 0}, photon{center.photon, 0},
      baryon{center.baryon, 0}, cdm{center.cdm, 0}, vacuum{retained_lambda, 0};
  const Ball one{1, 0}, two{2, 0}, three{3, 0}, four{4, 0};
  const Ball a2 = g.multiply(scale, scale), a4 = g.multiply(a2, a2);
  const Ball matter = g.multiply(g.add(baryon, cdm), scale);
  const Ball ba = g.multiply(baryon, scale), va = g.multiply(vacuum, a4);
  const Ball P = g.add(photon, g.add(matter, va));
  const Ball PN = g.add(matter, g.multiply(four, va));
  const Ball D = g.add(g.multiply(four, photon), g.multiply(three, ba));
  const Ball photon_enthalpy = g.divide(g.multiply(four, photon), three);
  const Ball F = g.divide(g.add(matter, photon_enthalpy), P);
  const Ball B = g.divide(g.add(ba, photon_enthalpy), P);
  const Ball half_q = g.divide(PN, g.multiply(two, P));
  const Ball L = g.subtract(half_q, one);
  const Ball T = g.divide(g.multiply(three, ba), D);
  const Ball cs = g.divide(photon_enthalpy, D);

  // Separate full-family clock graph at the SAME actual stage a. Physical
  // primitive radii enter the factors below only multiplied by clock support;
  // they are not added again as a map/background/source forcing amplitude.
  const Ball family_photon{center.photon, amplitude[0]},
      family_baryon{center.baryon, amplitude[1]},
      family_cdm{center.cdm, amplitude[2]},
      family_vacuum{retained_lambda, vacuum_amplitude};
  const Ball fm = g.multiply(g.add(family_baryon, family_cdm), scale),
             fb = g.multiply(family_baryon, scale),
             fv = g.multiply(family_vacuum, a4);
  const Ball fP = g.add(family_photon, g.add(fm, fv));
  const Ball fD = g.add(g.multiply(four, family_photon), g.multiply(three, fb));
  const Ball fe = g.divide(g.multiply(four, family_photon), three);
  const Ball fF = g.divide(g.add(fm, fe), fP), fB = g.divide(g.add(fb, fe), fP),
             fq = g.divide(g.add(fm, g.multiply(four, fv)), fP),
             fT = g.divide(g.multiply(three, fb), fD), fcs = g.divide(fe, fD);
  const W r2 = a.times(relative_a, relative_a), r3 = a.times(r2, relative_a),
          r4 = a.times(r2, r2);
  const W tail4 = a.plus(a.times(6, r2), a.plus(a.times(4, r3), r4));
  const W P_clock_change = a.plus(
      a.times(a.plus(g.magnitude(fm), a.times(4, g.magnitude(fv))), relative_a),
      a.times(g.magnitude(fv), tail4));
  const W D_clock_change = a.times(a.times(3, g.magnitude(fb)), relative_a);
  const W fP_lo = a.lower(fP.center, fP.radius),
          fP_clock_lo = a.lower(fP.center, a.plus(fP.radius, P_clock_change)),
          fD_clock_lo = a.lower(fD.center, a.plus(fD.radius, D_clock_change));
  const W P_lo = a.lower(P.center, P.radius), P_hi = a.plus(P.center, P.radius),
          fP_hi = a.plus(fP.center, fP.radius);
  if (a.status != S::ok)
    return a.status;
  const Ball F_clock_linear = g.subtract(g.multiply(fm, g.subtract(one, fF)),
                                         g.multiply(four, g.multiply(fF, fv)));
  const Ball B_clock_linear = g.subtract(g.subtract(fb, g.multiply(fB, fm)),
                                         g.multiply(four, g.multiply(fB, fv)));
  const Ball L_clock_linear =
      g.add(g.multiply(fm, g.subtract(one, fq)),
            g.multiply(fv, g.subtract(Ball{16, 0}, g.multiply(four, fq))));
  const Ball L_clock_tail = g.multiply(fv, g.subtract(four, fq));
  const Ball T_clock_linear =
      g.multiply(g.multiply(three, fb), g.subtract(one, fT));
  const Ball cs_clock_linear = g.multiply(fcs, g.multiply(three, fb));
  const W clock_F =
      a.over(a.plus(a.times(g.magnitude(F_clock_linear), relative_a),
                    a.times(g.magnitude(g.multiply(fF, fv)), tail4)),
             fP_clock_lo);
  const W clock_B =
      a.over(a.plus(a.times(g.magnitude(B_clock_linear), relative_a),
                    a.times(g.magnitude(g.multiply(fB, fv)), tail4)),
             fP_clock_lo);
  const W clock_L =
      a.over(a.plus(a.times(g.magnitude(L_clock_linear), relative_a),
                    a.times(g.magnitude(L_clock_tail), tail4)),
             a.product_lower(2, fP_clock_lo));
  const W clock_T =
      a.over(a.times(g.magnitude(T_clock_linear), relative_a), fD_clock_lo);
  const W clock_cs =
      a.over(a.times(g.magnitude(cs_clock_linear), relative_a), fD_clock_lo);
  // For nonnegative retained Gamma/matter/Lambda, exact |g_H|<=1. Thus H
  // changes by at most rho/(1-rho); x2 changes by at most
  // (2rho-rho²)/(1-rho)², enclosed here by the positive 2rho+rho² numerator.
  // Retain arithmetic×clock cross products. H and x2
  // centers remain the actual borrowed shared witnesses, without a new sqrt.
  const W h_relative = a.over(relative_a, scale_lower);
  const W x_relative = a.over(a.times(relative_a, a.plus(2, relative_a)),
                              a.product_lower(scale_lower, scale_lower));
  // Ratios >=1 enclose source-family magnitudes without another H/sqrt.
  // sqrt(r)<=r for r>=1, while x2 is exactly inverse-linear in P at fixed a.
  const W h_family_ratio = a.over(fP_hi, P_lo),
          x_family_ratio = a.over(P_hi, fP_lo);
  if (!(h_family_ratio >= 1 && x_family_ratio >= 1))
    return S::conditioning_budget_exceeded;
  const W h_family_upper =
      a.times(a.plus(epoch.hcal, shadow.actual_hcal_arithmetic_estimate),
              h_family_ratio);
  const W x_family_upper = a.times(
      a.plus(epoch.x2, shadow.actual_x2_arithmetic_estimate), x_family_ratio);
  const W nominal_h_lower =
      a.lower(epoch.hcal, shadow.actual_hcal_arithmetic_estimate);
  const W inverse_h_family_upper = a.over(x_family_ratio, nominal_h_lower);
  const W eta_clock = a.times(inverse_h_family_upper, h_relative);
  // The held inverse-H kernel uses R/[Hactual(Hactual-R)]. In addition to the
  // absolute H clock jump, choose R>=Hactual²*eta_clock; its kernel is then
  // >=eta_clock across the physical family, without a nominal-denominator-only
  // omission of physical×clock cross. R0+R is superadditive in that kernel.
  const W h_clock =
      std::max(a.times(h_family_upper, h_relative),
               a.times(a.times(epoch.hcal, epoch.hcal), eta_clock));
  const W h_radius = a.plus(shadow.actual_hcal_arithmetic_estimate, h_clock);
  const W x_radius = a.plus(shadow.actual_x2_arithmetic_estimate,
                            a.times(x_family_upper, x_relative));
  // Gradient target stays at actual a. Whole-family clock forcing above
  // already owns that physical×clock contribution, so gradients do not repeat
  // it.
  const Ball H{epoch.hcal, shadow.actual_hcal_arithmetic_estimate},
      X{epoch.x2, shadow.actual_x2_arithmetic_estimate};
  if (!(a.lower(epoch.hcal, h_radius) > 0) ||
      !(a.lower(epoch.x2, x_radius) > 0) || a.status != S::ok)
    return a.status == S::ok ? S::conditioning_budget_exceeded : a.status;

  out.actual = actual;
  out.center = center;
  out.directions = directions;
  out.hcal = epoch.hcal;
  // x stays unavailable here. The caller borrows the native readout's one
  // already computed sqrt(x2) witness only when requesting response readouts.
  out.arithmetic.hcal_radius = h_radius;
  out.arithmetic.shadow_p_radius = g.difference(center.shadow_p, P);
  out.arithmetic.shadow_p_n_radius = g.difference(center.shadow_p_n, PN);
  const W d_center =
      a.add(a.mul(4, center.photon), a.mul(a.mul(3, center.baryon), center.a));
  out.arithmetic.sound_denominator_radius = g.difference(d_center, D);
  const std::array<Ball, 6> coefficients{X, F, B, L, T, cs};
  const std::array<W, 6> actual_values{actual.x2,
                                       actual.F,
                                       actual.B,
                                       actual.L,
                                       actual.loading_over_one_plus_loading,
                                       actual.sound_speed_squared};
  const std::array<W, 6> clock_coefficients{a.times(x_family_upper, x_relative),
                                            clock_F,
                                            clock_B,
                                            clock_L,
                                            clock_T,
                                            clock_cs};
  for (unsigned i = 0; i < 6; ++i) {
    out.arithmetic.coefficient_at_actual_a_radius[i] =
        g.difference(actual_values[i],coefficients[i]);
    out.arithmetic.coefficient_radius[i] = a.plus(
        out.arithmetic.coefficient_at_actual_a_radius[i], clock_coefficients[i]);
  }

  // SAME finite nominal polynomial at actual a in every gradient.
  // Structural zero b/c-independent sound directions stay exact zero. The
  // rational 4/3 construction and all freshly reused expressions have their
  // own operation ancestry, including cancellation in fn/bn/ln.
  const Ball zero{};
  const Ball four_thirds = g.divide(four, three);
  const std::array<Ball, 4> pj{g.subtract(one, a4), g.subtract(scale, a4),
                               g.subtract(scale, a4), a4};
  const Ball three_a = g.multiply(three, scale);
  const std::array<Ball, 4> dj{four, three_a, zero, zero};
  const std::array<Ball, 4> tnj{zero, three_a, zero, zero};
  const std::array<Ball, 4> csnj{four_thirds, zero, zero, zero};
  const Ball Fa4 = g.multiply(F, a4), Ba4 = g.multiply(B, a4);
  const std::array<Ball, 4> fn{
      g.add(g.subtract(four_thirds, F), Fa4),
      g.add(g.multiply(scale, g.subtract(one, F)), Fa4),
      g.add(g.multiply(scale, g.subtract(one, F)), Fa4), g.subtract(zero, Fa4)};
  const std::array<Ball, 4> bn{
      g.add(g.subtract(four_thirds, B), Ba4),
      g.add(g.multiply(scale, g.subtract(one, B)), Ba4),
      g.add(g.subtract(zero, g.multiply(B, scale)), Ba4),
      g.subtract(zero, Ba4)};
  const Ball q_vacuum = g.multiply(a4, g.subtract(half_q, two));
  const Ball q_matter =
      g.multiply(scale, g.subtract(g.divide(one, two), half_q));
  const std::array<Ball, 4> ln{
      g.add(g.subtract(zero, half_q), q_vacuum), g.add(q_matter, q_vacuum),
      g.add(q_matter, q_vacuum), g.multiply(a4, g.subtract(two, half_q))};
  const Ball inverse_h = g.divide(one, H);
  for (unsigned j = 0; j < 4; ++j) {
    const std::array<Ball, 6> gradient{
        g.divide(g.subtract(zero, g.multiply(X, pj[j])), P),
        g.divide(fn[j], P),
        g.divide(bn[j], P),
        g.divide(ln[j], P),
        j < 2 ? g.divide(g.subtract(tnj[j], g.multiply(T, dj[j])), D) : zero,
        j < 2 ? g.divide(g.subtract(csnj[j], g.multiply(cs, dj[j])), D) : zero};
    const auto &used = directions.gradient[j];
    const std::array<W, 6> used_values{used.x2,
                                       used.F,
                                       used.B,
                                       used.L,
                                       used.loading_fraction,
                                       used.sound_speed_squared};
    for (unsigned i = 0; i < 6; ++i) {
      if (!normal(used_values[i]) || (j >= 2 && i >= 4 && used_values[i] != 0))
        return S::conditioning_budget_exceeded;
      out.arithmetic.gradient_radius[j][i] =
          g.difference(used_values[i], gradient[i]);
    }
    const W used_eta = -a.div(a.mul(a.div(1, epoch.hcal), pj[j].center),
                              a.mul(2, center.shadow_p));
    const Ball eta_gradient = g.divide(
        g.subtract(zero, g.multiply(inverse_h, pj[j])), g.multiply(two, P));
    out.arithmetic.eta_gradient_radius[j] =
        g.difference(used_eta, eta_gradient);
  }
  if (a.status != S::ok)
    return a.status;
  for (W v : out.arithmetic.coefficient_radius)
    if (!positive(v))
      return S::conditioning_budget_exceeded;
  for (W v : {out.arithmetic.hcal_radius, out.arithmetic.shadow_p_radius,
              out.arithmetic.shadow_p_n_radius,
              out.arithmetic.sound_denominator_radius})
    if (!positive(v))
      return S::conditioning_budget_exceeded;
  for (unsigned j = 0; j < 4; ++j) {
    if (!positive(out.arithmetic.eta_gradient_radius[j]))
      return S::conditioning_budget_exceeded;
    for (unsigned i = 0; i < 6; ++i)
      if (!radius(out.arithmetic.gradient_radius[j][i]) ||
          ((i < 4 || j < 2) && !positive(out.arithmetic.gradient_radius[j][i])))
        return S::conditioning_budget_exceeded;
  }
  // Actual Rbg remains in actual and is injected ONCE by transport_rhs. It is
  // not copied into an independent Dp/max-constraint/additional force here.
  // Additional RHS defects stay caller-owned and are not certified by this
  // coefficient bridge. No complete seed/RK/output admission follows from ok.
  out.arithmetic.conditional_wide_clock_profile = true;
  out.arithmetic.status = S::ok;
  out.status = S::ok;
  return S::ok;
}
} // namespace irred::cosmology::detail
