#pragma once
#include "ideal_acoustic_equations.hpp"
#include "irred/numerics.hpp"
#include <array>
#include <cmath>
#include <initializer_list>

namespace irred::cosmology::detail {
// Model-specific diagnostic directions for the selected, fixed no-relic
// physical source. They do not provide a second nominal background or H.
// The caller owns the original mapper radii and once-captured Lambda facts.
struct IdealAcousticSourceCenter {
  long double a, photon, baryon, cdm;
  long double shadow_p, shadow_p_n;
  long double x2, F, B, L, loading_fraction, sound_speed_squared;
};
struct IdealAcousticSourceGradient {
  long double x2=0, F=0, B=0, L=0;
  long double loading_fraction=0, sound_speed_squared=0;
};
struct IdealAcousticSourceDirections {
  numerics::Status status=numerics::Status::invalid_input;
  // photon, baryon, CDM with correlated dLambda=-sum(ddensity), then
  // the independent retained Lambda-preparation arithmetic remainder.
  std::array<IdealAcousticSourceGradient,4> gradient{};
};
inline IdealAcousticSourceDirections ideal_acoustic_source_directions(
    const IdealAcousticSourceCenter &c) noexcept {
  using S=numerics::Status; using W=long double;
  IdealAcousticSourceDirections out;
  if (!(c.a>0 && c.a<=W(.01)) || !(c.photon>0) || !(c.baryon>=0) ||
      !(c.cdm>=0) || !(c.shadow_p>0) || !(c.shadow_p_n>=0) || !(c.x2>=0))
    return out;
  for (W v:{c.a,c.photon,c.baryon,c.cdm,c.shadow_p,c.shadow_p_n,
            c.x2,c.F,c.B,c.L,c.loading_fraction,c.sound_speed_squared})
    if (!std::isfinite(v)) return out;
  const W a2=c.a*c.a,a4=a2*a2,D=4*c.photon+3*c.baryon*c.a;
  if (!(D>0) || !std::isnormal(D) || !std::isnormal(c.shadow_p)) {
    out.status=S::outside_domain; return out;
  }
  const std::array<W,4> pj{1-a4,c.a-a4,c.a-a4,a4};
  const std::array<W,4> dj{4,3*c.a,0,0},tnj{0,3*c.a,0,0},csnj{4.L/3,0,0,0};
  const W half_q=c.shadow_p_n/(2*c.shadow_p);
  // Literal correlated affine grouping retains the vacuum support even when
  // 1-a^4 rounds to one. Do not replace these signed directions by unrelated
  // fraction maxima or drop their a^4 term after a radiation cancellation.
  const std::array<W,4> fn{
      (4.L/3-c.F)+c.F*a4,c.a*(1-c.F)+c.F*a4,
      c.a*(1-c.F)+c.F*a4,-c.F*a4};
  const std::array<W,4> bn{
      (4.L/3-c.B)+c.B*a4,c.a*(1-c.B)+c.B*a4,
      -c.B*c.a+c.B*a4,-c.B*a4};
  const std::array<W,4> ln{
      -half_q+a4*(half_q-2),c.a*(.5L-half_q)+a4*(half_q-2),
      c.a*(.5L-half_q)+a4*(half_q-2),a4*(2-half_q)};
  for (unsigned j=0;j<4;++j) {
    auto &g=out.gradient[j];
    g.x2=-c.x2*pj[j]/c.shadow_p;
    g.F=fn[j]/c.shadow_p;
    g.B=bn[j]/c.shadow_p;
    g.L=ln[j]/c.shadow_p;
    g.loading_fraction=(tnj[j]-c.loading_fraction*dj[j])/D;
    g.sound_speed_squared=(csnj[j]-c.sound_speed_squared*dj[j])/D;
    for (W v:{g.x2,g.F,g.B,g.L,g.loading_fraction,g.sound_speed_squared})
      if (!std::isfinite(v)) { out.status=S::overflow; return out; }
  }
  out.status=S::ok; return out;
}
// The signed source force A_j Y preserves radiation/adiabatic cancellations.
// Acceleration-defect differences belong to actual-versus-shadow arithmetic,
// not an extra physical source direction: every exact ideal source has Rbg=0.
template<unsigned I>
inline long double ideal_acoustic_source_force_coordinate(
    std::span<const long double,5> y,const IdealAcousticSourceGradient &g) noexcept {
  static_assert(I<5);
  if constexpr (I==0)
    return -g.x2*y[3]+4.5L*g.B*y[2];
  else if constexpr (I==1)
    return g.x2*y[2];
  else if constexpr (I==2)
    return (g.L-g.loading_fraction)*y[2]+g.sound_speed_squared*(y[0]-y[1]);
  else if constexpr (I==3)
    return g.L*y[3];
  else
    return 1.5L*g.F*y[3]+1.5L*g.B*y[2];
}
inline std::array<long double,5> ideal_acoustic_source_force(
    std::span<const long double,5> y,const IdealAcousticSourceGradient &g) noexcept {
  return {ideal_acoustic_source_force_coordinate<0>(y,g),
          ideal_acoustic_source_force_coordinate<1>(y,g),
          ideal_acoustic_source_force_coordinate<2>(y,g),
          ideal_acoustic_source_force_coordinate<3>(y,g),
          ideal_acoustic_source_force_coordinate<4>(y,g)};
}
struct IdealAcousticQuotientRemainder {
  numerics::Status status=numerics::Status::invalid_input;
  long double absolute_estimate=0;
};
// Exact real identity: f(n+dn,d+dd)-f-f_linear =
// -(dd/d)*(dn-f*dd)/(d+dd). Its finite-radius bound is useful only
// together with the actual center/gradient operation and cast estimates.
// This function supplies the mathematical quotient term, not a complete
// source/arithmetic certificate or a zero-on-missing admission route.
inline IdealAcousticQuotientRemainder ideal_acoustic_quotient_remainder(
    long double f,long double denominator,long double numerator_radius,
    long double denominator_radius) noexcept {
  using S=numerics::Status;
  IdealAcousticQuotientRemainder out;
  if (!(denominator>0) || !(numerator_radius>=0) || !(denominator_radius>=0) ||
      !std::isfinite(f) || !std::isfinite(denominator) ||
      !std::isfinite(numerator_radius) || !std::isfinite(denominator_radius)) return out;
  const auto lower=denominator-denominator_radius;
  if (!(lower>0) || !std::isnormal(lower)) {
    out.status=S::conditioning_budget_exceeded; return out;
  }
  out.absolute_estimate=(denominator_radius/denominator)*
      (numerator_radius+std::abs(f)*denominator_radius)/lower;
  if (!std::isfinite(out.absolute_estimate) ||
      (out.absolute_estimate!=0 && !std::isnormal(out.absolute_estimate))) {
    out.status=S::outside_domain; return out;
  }
  out.status=S::ok; return out;
}
} // namespace irred::cosmology::detail
