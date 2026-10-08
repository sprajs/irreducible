#pragma once
#include "ideal_acoustic_bridge.hpp"

namespace irred::cosmology::detail {
namespace ideal_acoustic_initial_internal {
using W = long double;
using S = numerics::Status;
using Ball = ideal_acoustic_bridge_internal::Ball;
using Graph = ideal_acoustic_bridge_internal::Graph;

// Exact declared finite-start map. Its Taylor truncation is not an error
// estimate here: the independent start refinement owns that approximation.
struct Map {
  Ball matter;
  std::array<Ball, 5> y;
};
inline Map map(Graph &g, Ball matter, Ball q, Ball F, Ball B) noexcept {
  const Ball zero{}, one{1, 0}, three{3, 0};
  const Ball p = g.divide(Ball{-2, 0}, three);
  const Ball m16 = g.divide(matter, Ball{16, 0});
  const Ball m2 = g.multiply(matter, matter);
  const Ball mterm = g.divide(g.multiply(Ball{7, 0}, m2), Ball{160, 0});
  const Ball phi = g.multiply(p, g.subtract(g.add(g.subtract(one, m16), mterm),
                                            g.divide(q, Ball{30, 0})));
  const Ball vc =
      g.multiply(p, g.subtract(g.subtract(g.add(Ball{.5L, 0}, m16), mterm),
                               g.divide(q, Ball{120, 0})));
  const Ball entropy =
      g.subtract(zero, g.divide(g.multiply(p, g.multiply(q, q)), Ball{96, 0}));
  const Ball velocity =
      g.subtract(zero, g.divide(g.multiply(p, q), Ball{24, 0}));
  const Ball numerator = g.subtract(
      g.add(g.multiply(g.multiply(p, q), phi), g.multiply(B, entropy)),
      g.multiply(three, g.multiply(B, velocity)));
  return {matter, {g.divide(numerator, F), entropy, velocity, vc, phi}};
}

// Finite polynomial differences and the actual Delta-only quotient projection.
// Every positive product includes the corresponding finite cross term.
inline std::array<W, 5> variation(Graph &g, W p, W dm, W dq, W dF, W dB,
                                  W phi_m_slope, W vc_m_slope, W q_hi, W B_hi,
                                  W F_lo,
                                  const std::array<W, 5> &initial_hi) noexcept {
  auto &a = g.a;
  const W quadratic = a.over(a.times(7, a.times(dm, dm)), 160);
  const W dphi = a.times(
      p, a.plus(a.times(phi_m_slope, dm), a.plus(quadratic, a.over(dq, 30))));
  const W dvc = a.times(
      p, a.plus(a.times(vc_m_slope, dm), a.plus(quadratic, a.over(dq, 120))));
  const W ds = a.over(
      a.times(p, a.plus(a.times(a.times(2, q_hi), dq), a.times(dq, dq))), 96);
  const W ddv = a.over(a.times(p, dq), 24);
  const W wave =
      a.times(p, a.plus(a.times(q_hi, dphi),
                        a.plus(a.times(initial_hi[4], dq), a.times(dq, dphi))));
  const W entropy = a.plus(a.times(B_hi, ds),
                           a.plus(a.times(initial_hi[1], dB), a.times(dB, ds)));
  const W velocity =
      a.times(3, a.plus(a.times(B_hi, ddv),
                        a.plus(a.times(initial_hi[2], dB), a.times(dB, ddv))));
  const W numerator = a.plus(wave, a.plus(entropy, velocity));
  const W delta =
      a.over(a.plus(numerator, a.times(initial_hi[0], dF)), a.lower(F_lo, dF));
  return {delta, ds, ddv, dvc, dphi};
}

// Recover the separately owned positive clock increment from the bridge's
// stored base+clock sum, including the finite subtraction ancestry.
inline W clock_increment(Graph &g, W full, W base) noexcept {
  auto &a = g.a;
  if (!(full >= base && base > 0)) {
    a.status = S::conditioning_budget_exceeded;
    return 0;
  }
  return a.plus(a.magnitude(a.sub(full, base)), a.loss(a.plus(full, base)));
}
} // namespace ideal_acoustic_initial_internal

// Input frame MUST be the actual borrowed-stage output of
// ideal_acoustic_transport_frame. SOURCE is the finite primitive family at
// that actual a. ARITHMETIC is the measured nominal-map discrepancy plus the
// entire positive-family clock change, including source-times-clock products.
// No second H, background query, source capture/map, sqrt, log or exp occurs.
// The caller supplies complete positive early-age estimates; eta SOURCE is
// entirely a remainder seed, while transport_initial sets eta sensitivities0.
// Success provides initial bounds only, not campaign/output admission.
inline numerics::Status ideal_acoustic_initial_bounds(
    std::span<const long double, 6> y, const IdealAcousticTransportFrame &f,
    const IdealAcousticSourceUncertainty &u, long double retained_lambda,
    long double initial_age_source_radius,
    long double initial_age_arithmetic_radius,
    IdealAcousticTransportInitialBounds &out,
    const IdealAcousticTransportAccounting &work) noexcept {
  using namespace ideal_acoustic_initial_internal;
  using ideal_acoustic_bridge_internal::positive;
  using ideal_acoustic_transport_internal::normal_or_zero;
  out.status = S::invalid_input;
  if (!ideal_acoustic_transport_internal::valid(f, u) ||
      !(retained_lambda >= 0) || !normal_or_zero(retained_lambda) ||
      !positive(initial_age_source_radius) ||
      !positive(initial_age_arithmetic_radius) || !positive(y[5]))
    return out.status = S::conditioning_budget_exceeded;
  for (W v : y)
    if (!normal_or_zero(v))
      return out.status = S::outside_domain;
  const S charged = ideal_acoustic_transport_internal::account(work);
  if (charged != S::ok)
    return out.status = charged;
  Graph g;
  auto &a = g.a;
  ideal_acoustic_transport_internal::Bands bands;
  const S band_status =
      ideal_acoustic_transport_internal::bands(f, u, bands, a);
  if (band_status != S::ok)
    return out.status = band_status;
  const auto &base = f.arithmetic.coefficient_at_actual_a_radius;
  const auto &full = f.arithmetic.coefficient_radius;
  const Ball scale{f.center.a, 0}, photon{f.center.photon, 0};
  const Ball matter = g.divide(
      g.multiply(g.add(Ball{f.center.baryon, 0}, Ball{f.center.cdm, 0}), scale),
      photon);
  const Ball q{f.actual.x2, base[0]}, F{f.actual.F, base[1]},
      B{f.actual.B, base[2]};
  const Map nominal = map(g, matter, q, F, B);
  const Ball p = g.divide(Ball{-2, 0}, Ball{3, 0});
  const W p_hi = g.magnitude(p), m_hi = g.magnitude(matter);
  const W q_hi = g.magnitude(q), B_hi = g.magnitude(B);
  const W F_lo = a.lower(F.center, F.radius);
  const W vacuum_amplitude =
      a.plus(a.plus(bands.amplitude[0], bands.amplitude[1]),
             a.plus(bands.amplitude[2], bands.amplitude[3]));
  if (!(a.lower(retained_lambda, vacuum_amplitude) > 0))
    return out.status = S::conditioning_budget_exceeded;
  std::array<W, 5> nominal_hi;
  for (std::size_t i = 0; i < 5; ++i)
    nominal_hi[i] = g.magnitude(nominal.y[i]);
  const W dm = a.over(a.plus(a.times(f.center.a, a.plus(bands.amplitude[1],
                                                        bands.amplitude[2])),
                             a.times(m_hi, bands.amplitude[0])),
                      a.lower(f.center.photon, bands.amplitude[0]));
  const Ball m_slope = g.divide(g.multiply(Ball{7, 0}, matter), Ball{80, 0});
  const W phi_slope = g.magnitude(g.add(Ball{-1.L / 16, 0}, m_slope));
  const W vc_slope = g.magnitude(g.subtract(Ball{1.L / 16, 0}, m_slope));
  const auto source =
      variation(g, p_hi, dm, bands.change[0], bands.change[1], bands.change[2],
                phi_slope, vc_slope, q_hi, B_hi, F_lo, nominal_hi);

  const W clock_q = clock_increment(g, full[0], base[0]);
  const W clock_F = clock_increment(g, full[1], base[1]);
  const W clock_B = clock_increment(g, full[2], base[2]);
  if (!(full[0] > base[0]))
    return out.status = S::conditioning_budget_exceeded;
  // The actual bridge stores Cq >= Xfamily_hi*(2rho+rho^2)/(1-rho)^2
  // with Xfamily_hi >= actual x2. Thus this outward recovered rho encloses
  // that same selected clock support. Arbitrary hand-written frame radii do
  // not establish this provenance. No missing clock is inferred as zero.
  const W rho = a.over(clock_q, a.product_lower(2, f.actual.x2));
  if (!(rho > 0 && rho < .25L))
    return out.status = S::conditioning_budget_exceeded;
  const W family_m_hi = a.plus(m_hi, dm);
  const W family_q_hi = a.plus(q_hi, bands.change[0]);
  const W family_B_hi = a.plus(B_hi, bands.change[2]);
  const W family_F_lo = a.lower(F_lo, bands.change[1]);
  std::array<W, 5> family_hi;
  for (std::size_t i = 0; i < 5; ++i)
    family_hi[i] = a.plus(nominal_hi[i], source[i]);
  const W clock_dm = a.times(family_m_hi, rho);
  const W family_slope =
      a.plus(a.up(1.L / 16), a.over(a.times(7, family_m_hi), 80));
  const auto clock =
      variation(g, p_hi, clock_dm, clock_q, clock_F, clock_B, family_slope,
                family_slope, family_q_hi, family_B_hi, family_F_lo, family_hi);

  // Bound the integral-age displacement by the complete family inverse-H
  // magnitude and finite log-distance. |g_H|<=1 yields one extra 1/(1-rho)
  // over this recovered support; no new inverse-H source or H query is made.
  const W a2 = a.times(f.center.a, f.center.a), a4 = a.times(a2, a2);
  const W dp = a.plus(bands.amplitude[0],
                      a.plus(a.times(f.center.a, a.plus(bands.amplitude[1],
                                                        bands.amplitude[2])),
                             a.times(a4, vacuum_amplitude)));
  const W P_hi = a.plus(f.center.shadow_p, f.arithmetic.shadow_p_radius);
  const W P_lo = a.lower(f.center.shadow_p, f.arithmetic.shadow_p_radius);
  const W family_P_lo = a.lower(P_lo, dp);
  const W inverse_h_family =
      a.times(a.plus(bands.inverse_h, bands.inverse_h_arithmetic),
              a.over(P_hi, family_P_lo));
  const W scale_lo = a.lower(1, rho);
  const W eta_clock = a.over(a.times(inverse_h_family, rho),
                             a.product_lower(scale_lo, scale_lo));
  std::array<W, 5> arithmetic;
  for (std::size_t i = 0; i < 5; ++i)
    arithmetic[i] = a.plus(g.difference(y[i], nominal.y[i]), clock[i]);
  const W eta_arithmetic = a.plus(initial_age_arithmetic_radius,
                                  a.plus(eta_clock, a.loss(a.magnitude(y[5]))));
  if (a.status != S::ok)
    return out.status = a.status;
  for (std::size_t i = 0; i < 5; ++i)
    if (!positive(source[i]) || !positive(arithmetic[i]))
      return out.status = S::conditioning_budget_exceeded;
  if (!positive(eta_arithmetic))
    return out.status = S::conditioning_budget_exceeded;
  // Only final owned error-result destinations consume state_element_writes;
  // all graph/interval assembly above belongs to the charged diagnostic.
  if (!work.state_writes(work.context, 12))
    return out.status = S::work_limit;
  for (std::size_t i = 0; i < 5; ++i) {
    out.source_family_radius[i] = source[i];
    out.arithmetic_radius[i] = arithmetic[i];
  }
  out.source_family_radius[5] = initial_age_source_radius;
  out.arithmetic_radius[5] = eta_arithmetic;
  return out.status = S::ok;
}
} // namespace irred::cosmology::detail
