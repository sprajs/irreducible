// Original direct-z/sqrt(a) GL8 polynomial controls and explicit 3x3 cofactors.
// Massive facts retain thermal_observables_peer mpmath60/90 digits direct
// momentum GL32/48 tails128/160 ancestry and CLASS3.3.0
// 0ceb7a9a4c1e444ef5d5d56a8328a0640be91b18 explicit constants/species.
// That reference's refinement <2.86e-9 of ratio allocation and CLASS <.00177;
// common equations/constants are shared ancestry, not independent measurements.
// Ratio comparison 1e-10+5e-10|ref|; density components absolute1e-8.
// New polynomial ratio and propagated density refinements consume <=5%.
#include "test_bao_thermal_fixture.hpp"
#include <iostream>
#include <numbers>
using namespace thermal_bao_test;
namespace {
constexpr W ratio_absolute = 1e-10L, ratio_relative = 5e-10L,
            density_absolute = 1e-8L;
W allowance(W ref) { return ratio_absolute + ratio_relative * std::abs(ref); }
template <class F> W gl(F f, W end, unsigned n) {
  constexpr W x[]{
      .183434642495649804939476142360184L, .525532409916328985817739049189246L,
      .796666477413626739591553936475831L, .960289856497536231683560868569473L};
  constexpr W w[]{
      .362683783378361982965150449277195L, .313706645877887287337962201986601L,
      .222381034453374470544355994426240L, .101228536290376259152531354309962L};
  W sum = 0, h = end / (2 * n);
  for (unsigned i = 0; i < n; ++i)
    for (unsigned j = 0; j < 4; ++j) {
      W mid = (2 * i + 1) * h;
      sum += h * w[j] * (f(mid - h * x[j]) + f(mid + h * x[j]));
    }
  return sum;
}
W critical(W H) {
  const W pi = std::numbers::pi_v<W>, c = 299792458.L, ev = 1.602176634e-19L,
          hbar = 6.62607015e-34L / (2 * pi),
          mpc = 648000000000.L * 149597870700.L / pi;
  const W hs = H * 1000 / mpc;
  return 3 * hs * hs * c * c / (8 * pi * 6.67430e-11L) / ev *
         std::pow(hbar * c / ev, 3);
}
std::array<W, 3> polynomial(cosmology::ThermalObservableRequest q, W z,
                            unsigned n) {
  auto m = q.model;
  const W pi = std::numbers::pi_v<W>, kb = 1.380649e-23L / 1.602176634e-19L,
          h2 = std::pow(W(m.h0_km_s_mpc) / 100, 2),
          b = m.physical_baryon_density / h2,
          g = pi * pi / 15 * std::pow(kb * m.tcmb_kelvin, 4) /
              critical(m.h0_km_s_mpc),
          matter = b + m.physical_cdm_density / h2;
  W radiation = g + m.physical_massless_nonphoton_density / h2;
  for (auto sp : m.species) {
    need(sp.mass_ev == 0, "polynomial only massless FD limit");
    radiation += sp.statistical_weight * 7 * pi * pi / 240 *
                 std::pow(kb * sp.temperature_today_kelvin, 4) /
                 critical(m.h0_km_s_mpc);
  }
  W lambda = 1 - radiation - matter, scale = 299792.458L / m.h0_km_s_mpc;
  auto P = [&](W a) { return radiation + matter * a + lambda * a * a * a * a; };
  auto E = [&](W z0) {
    W a = 1 / (1 + z0);
    return std::sqrt(P(a)) / (a * a);
  };
  W rs = scale / std::sqrt(3.L) *
         gl(
             [&](W u) {
               W a = u * u;
               return 2 * u / std::sqrt(P(a)) /
                      std::sqrt(1 + 3 * b * a / (4 * g));
             },
             1 / std::sqrt(1 + W(q.z_drag)), n),
    dm = scale * gl([&](W t) { return 1 / E(t); }, z, n), dh = scale / E(z);
  return {dm / rs, dh / rs, std::cbrt(dm * dm * z * dh) / rs};
}
struct Density {
  W q, ld, norm, logp;
};
Density cofactors(const bao::DensityInput &s, std::span<const W> mu) {
  const auto &c = s.covariance;
  const W a = c[0], b = c[1], cc = c[2], d = c[4], e = c[5], f = c[8];
  const W inverse[]{d * f - e * e,  cc * e - b * f,  b * e - cc * d,
                    cc * e - b * f, a * f - cc * cc, b * cc - a * e,
                    b * e - cc * d, b * cc - a * e,  a * d - b * b};
  const W det = a * inverse[0] + b * inverse[1] + cc * inverse[2];
  need(det > 0, "independent cofactor SPD determinant");
  W r[3], q = 0;
  for (unsigned j = 0; j < 3; ++j)
    r[j] = W(s.observed[j]) - mu[j];
  for (unsigned j = 0; j < 3; ++j)
    for (unsigned k = 0; k < 3; ++k)
      q += r[j] * inverse[j * 3 + k] * r[k] / det;
  W ld = std::log(det), norm = 3 * std::log(2 * std::numbers::pi_v<W>);
  return {q, ld, norm, -(q + ld + norm) / 2};
}
void density_check(const bao::ThermalDensitySlot &s, Density ref) {
  if (!s.result)
    std::cerr << "missing density: status=" << unsigned(s.numerical_status)
              << " prepare=" << unsigned(s.preparation_status) << " projection="
              << s.projection_log_density_error_estimate.value_or(-1)
              << " callbacks=" << s.callbacks << '\n';
  need(s.result &&
           s.result->density.status == statistics::DensityStatus::finite,
       "strict normalized density available");
  near(s.result->quadratic, ref.q, density_absolute,
       "independent cofactor quadratic");
  near(s.result->log_determinant, ref.ld, density_absolute,
       "independent determinant");
  near(s.result->normalization, ref.norm, density_absolute,
       "ratio measure normalization");
  near(s.result->density.log_value, ref.logp, density_absolute,
       "independent normalized density");
  need(s.projection_log_density_error_estimate &&
           *s.projection_log_density_error_estimate <= 1e-8,
       "unchanged density projection allocation");
}
void massive_facts() {
  constexpr W ref[]{
      23.8845464757980689796700765015910639718009504300807168774140453259228452189617087589041495L,
      17.4466655220062015394759710599562292487011950165058252947423953870961312129398339569739476L,
      21.5104160838493437484767027864037581211333042395728507780758064035644645503060388743665683L};
  auto req = massive();
  auto d = prepared(ref);
  auto got = evaluate(d, req);
  need(got.slots.size() == 1, "massive one retained slot");
  auto &s = got.slots[0];
  need(s.predictions.size() == 3, "massive predictions");
  for (unsigned j = 0; j < 3; ++j)
    near(s.predictions[j], ref[j], allowance(ref[j]),
         "inherited independent high precision massive fact");
  density_check(s, cofactors(d.source(), ref));
  // CLASS independently supplied rs(z_drag), not its predicted drag epoch.
  W rs = 142.40563343054799L, dm = 3401.2939701017422L,
    dh = 299792.458L / 67.4L / 1.7902809950671763L;
  const W external[]{dm / rs, dh / rs, std::cbrt(dm * dm * dh) / rs};
  for (unsigned j = 0; j < 3; ++j)
    near(s.predictions[j], external[j], allowance(external[j]),
         "matched external CLASS ratio");
  density_check(s, cofactors(d.source(), external));
  auto perm = evaluate(prepared(ref, 1, true), req);
  density_check(perm.slots[0],
                cofactors(prepared(ref, 1, true).source(),
                          std::array<W, 3>{ref[2], ref[1], ref[0]}));
  near(perm.slots[0].result->density.log_value, s.result->density.log_value,
       1e-8, "row covariance permutation");
  auto p = policy();
  p.predictions.absolute_tolerance_mpc = 1e-10;
  p.predictions.relative_tolerance = 2e-12;
  p.predictions.absolute_tolerance_ratio = 1e-12;
  p.predictions.relative_tolerance_ratio = 5e-12;
  p.predictions.thermal.absolute_tolerance = 1e-14;
  p.predictions.thermal.relative_tolerance = 2e-14;
  // Frozen tighter-quality refinement originally exhausted20-million per-row
  // work (retained local witness). Increase explicit resource ceilings, leaving
  // the original acceptance and the tighter numerical settings unchanged.
  p.predictions.maximum_callbacks_per_point = 100000000;
  p.predictions.maximum_total_callbacks = 500000000;
  p.predictions.thermal.maximum_total_callbacks = 500000000;
  p.maximum_total_callbacks = 500000000;
  auto refined = evaluate(d, req, p);
  density_check(refined.slots[0], cofactors(d.source(), ref));
  near(refined.slots[0].result->density.log_value, s.result->density.log_value,
       1e-8, "production refinement retains density budget");
  // Tiny declared observation covariance is not inflated to absorb provider
  // diagnostics. Default1e-8 projection refuses while ratios remain available.
  auto in = d.source();
  for (auto &v : in.covariance)
    v *= 1e-8;
  auto tiny = bao::prepare_density(std::move(in),
                                   {16, 256, 4096, 8 * 1024 * 1024, 1e-8,
                                    numerics::Arithmetic::longdouble_cpu_v1});
  auto refused = evaluate(tiny, req);
  need(refused.slots[0].predictions.size() == 3 && !refused.slots[0].result &&
           refused.slots[0].numerical_status == S::conditioning_budget_exceeded,
       "density level default refusal leaves original predictions");
}
void analytic_and_scaling() {
  for (bool species : {false, true}) {
    auto req = massless();
    if (species)
      req.model.species.push_back({0, 1.9, 2});
    for (double z : {0., 1e-8, .7, 4.}) {
      auto ref = polynomial(req, z, 512), coarse = polynomial(req, z, 256);
      for (unsigned j = 0; j < 3; ++j)
        near(ref[j], coarse[j], .05L * allowance(ref[j]),
             "independent polynomial ratio refinement");
      auto d = prepared(ref, z);
      auto reference = cofactors(d.source(), ref),
           old = cofactors(d.source(), coarse);
      near(reference.q, old.q, .05L * density_absolute,
           "cofactor quadratic refinement");
      near(reference.logp, old.logp, .05L * density_absolute,
           "cofactor log density refinement");
      auto got = evaluate(d, req);
      need(got.slots[0].predictions.size() == 3,
           "analytic predictions available");
      for (unsigned j = 0; j < 3; ++j)
        near(got.slots[0].predictions[j], ref[j], allowance(ref[j]),
             "massless polynomial independent ratios");
      density_check(got.slots[0], reference);
    }
  }
  auto req = massless();
  auto ref = polynomial(req, .7, 512);
  auto d = prepared(ref, .7);
  auto before = evaluate(d, req);
  req.model.h0_km_s_mpc *= 2;
  req.model.physical_baryon_density *= 4;
  req.model.physical_cdm_density *= 4;
  req.model.physical_massless_nonphoton_density *= 4;
  req.model.tcmb_kelvin *= std::sqrt(2.);
  auto fractions = evaluate(d, req);
  density_check(fractions.slots[0], cofactors(d.source(), ref));
  for (unsigned j = 0; j < 3; ++j)
    near(fractions.slots[0].predictions[j], ref[j], allowance(ref[j]),
         "fixed fractions explicit scaling");
  req = massless();
  req.model.h0_km_s_mpc *= 2;
  auto changed = polynomial(req, .7, 512);
  auto physical = evaluate(d, req);
  need(physical.slots[0].predictions.size() == 3,
       "fixed physical density ratios");
  for (unsigned j = 0; j < 3; ++j)
    near(physical.slots[0].predictions[j], changed[j], allowance(changed[j]),
         "fixed physical density independent scaling");
  need(std::abs(physical.slots[0].predictions[0] -
                before.slots[0].predictions[0]) > .1,
       "no false fixed physical H0 cancellation");
  density_check(physical.slots[0], cofactors(d.source(), changed));
  req = massless();
  req.z_drag = 1200;
  auto drag = evaluate(d, req);
  need(drag.slots[0].predictions[0] > before.slots[0].predictions[0],
       "conditional drag dependence");
}
} // namespace
int main() {
  try {
    massive_facts();
    analytic_and_scaling();
    std::cout << "PASS " << checks << " thermal BAO peer controls\n";
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
