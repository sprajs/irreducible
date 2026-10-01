// Independent synthetic controls; no copied reference code/assets. SI
// c=299792458 m/s exact, shared defining constant is not independent physics
// evidence. Distances: composite Simpson directly in z, unlike production
// log-redshift. Ruler: u=sqrt(a), independently derived dz/H -> 2u
// du/sqrt(r+m*u²+l*u8). Frozen budgets: distances 1e-9+2e-11|ref| Mpc;
// ratios1e-11+5e-11|ref|; reference refinement <=5% of each budget. No
// inference/released-data claim.
#include "irred/early_late.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace irred::cosmology;
namespace {
using W = long double;
unsigned checks = 0;
void need(bool b, const char *s) {
  ++checks;
  if (!b)
    throw std::runtime_error(s);
}
template <class F> W simpson(F f, W end, unsigned n) {
  W sum = f(0) + f(end);
  for (unsigned i = 1; i < n; ++i)
    sum += (i % 2 ? 4 : 2) * f(end * i / n);
  return sum * end / (3 * n);
}
W ez(EarlyFlatModel m, W z) {
  W y = 1 + z;
  return std::sqrt(W(m.omega_r) * y * y * y * y + W(m.omega_m) * y * y * y + 1 -
                   W(m.omega_m) - m.omega_r);
}
W dc(EarlyFlatModel m, W z, unsigned n) {
  return 299792.458L / m.h0_km_s_mpc *
         simpson([&](W x) { return 1 / ez(m, x); }, z, n);
}
W rs(SoundHorizonRequest r, unsigned n) {
  W b = 3 * W(r.model.omega_b) / (4 * r.model.omega_gamma),
    end = 1 / std::sqrt(1 + W(r.z_drag)),
    l = 1 - W(r.model.omega_m) - r.model.omega_r;
  return 299792.458L / r.model.h0_km_s_mpc / std::sqrt(3.L) *
         simpson(
             [&](W u) {
               W a = u * u;
               return 2 * u /
                      std::sqrt(W(r.model.omega_r) + W(r.model.omega_m) * a +
                                l * a * a * a * a) /
                      std::sqrt(1 + b * a);
             },
             end, n);
}
W val(const EarlyLateRow &r, EarlyLateOutput o) {
  const auto &v = r.outputs[unsigned(o)];
  need(v.status == irred::numerics::Status::ok && v.value.has_value(),
       "requested finite output");
  return *v.value;
}
void near(W a, W b, W abs, W rel, const char *s) {
  need(std::abs(a - b) <= abs + rel * std::abs(b), s);
}
EarlyLatePolicy policy() {
  return {1e-9,
          2e-11,
          1e-11,
          5e-11,
          1000000,
          4000000,
          40,
          64,
          4 * 1024 * 1024,
          {1e-9, 2e-11, 1000000, 40, 64, 4000000, 4 * 1024 * 1024}};
}
} // namespace
int main() {
  try {
    need(std::numeric_limits<W>::digits >= 64, "wide reference arithmetic");
    unsigned mask = (1u << early_late_output_count) - 1;
    const double zs[]{0, .001, .1, .7, 2, 4};
    for (auto model : {EarlyFlatModel{70, .3, 9e-5, .05, 5e-5},
                       EarlyFlatModel{61, .2, .1, .03, .04},
                       EarlyFlatModel{70, 0, 1, 0, .1}}) {
      SoundHorizonRequest req{model, 1059, "synthetic supplied drag"};
      auto got = evaluate_early_late(req, zs, mask, policy());
      need(got.status == irred::numerics::Status::ok && got.rows.size() == 6 &&
               got.ruler.has_value(),
           "batch rows and single ruler");
      W ruler = rs(req, 8192), coarse = rs(req, 4096);
      near(ruler, coarse, .05e-9L, .05L * 2e-11L,
           "ruler refinement allocation");
      near(*got.ruler->sound_horizon_mpc, ruler, 1e-9L, 2e-11L,
           "independent u ruler");
      size_t callbacks = got.ruler->callbacks;
      for (size_t i = 0; i < 6; ++i) {
        W z = zs[i], e = ez(model, z), h = model.h0_km_s_mpc * e,
          dh = 299792.458L / h, dm = dc(model, z, 8192), dl = (1 + z) * dm,
          dv = std::cbrt(z * dh * dm * dm);
        callbacks += got.rows[i].callbacks;
        near(dm, dc(model, z, 4096), .05e-9L, .05L * 2e-11L,
             "z refinement allocation");
        near(val(got.rows[i], EarlyLateOutput::e), e, 0, 2e-11L,
             "independent E");
        near(val(got.rows[i], EarlyLateOutput::h_km_s_mpc), h, 0, 2e-11L,
             "independent H");
        near(val(got.rows[i], EarlyLateOutput::dh_mpc), dh, 0, 2e-11L,
             "independent DH");
        for (auto pair : {std::pair{EarlyLateOutput::dm_mpc, dm},
                          std::pair{EarlyLateOutput::dl_mpc, dl},
                          std::pair{EarlyLateOutput::dv_mpc, dv}})
          near(val(got.rows[i], pair.first), pair.second, 1e-9L, 2e-11L,
               "independent distance");
        for (auto pair : {std::pair{EarlyLateOutput::dm_over_rs, dm / ruler},
                          std::pair{EarlyLateOutput::dh_over_rs, dh / ruler},
                          std::pair{EarlyLateOutput::dv_over_rs, dv / ruler}})
          near(val(got.rows[i], pair.first), pair.second, 1e-11L, 5e-11L,
               "conditional ratio");
        if (model.omega_r == 1)
          near(dm, 299792.458L / model.h0_km_s_mpc * z / (1 + z), 1e-9L, 2e-11L,
               "exact radiation distance");
      }
      need(callbacks == got.callbacks, "single ruler callback accounting");
      req.model.h0_km_s_mpc *= 2;
      auto scaled = evaluate_early_late(req, zs, mask, policy());
      need(scaled.rows.size() == 6, "scaled rows");
      for (size_t i = 0; i < 6; ++i) {
        near(val(scaled.rows[i], EarlyLateOutput::dm_mpc) * 2,
             val(got.rows[i], EarlyLateOutput::dm_mpc), 1e-9L, 2e-11L,
             "H0 inverse distance");
        near(val(scaled.rows[i], EarlyLateOutput::dm_over_rs),
             val(got.rows[i], EarlyLateOutput::dm_over_rs), 1e-11L, 5e-11L,
             "H0 ratio cancellation");
      }
    }
    // Astropy8.0.1 FlatLambdaCDM named comparison, independently acquired BSD-3
    // library numerical facts, not observational truth. Exact massless mapping:
    // H0=70 Om0=.3 Ob0=.05 Tcmb0=2.7255K Neff=3.046 m_nu=0eV;
    // Or=Ogamma0+Onu0, so no implicit conversion enters the product API.
    // Official contract:
    // https://docs.astropy.org/en/stable/api/astropy.cosmology.FlatLambdaCDM.html
    const double az[]{.001, .1, .7, 2., 4.};
    const W ae[]{1.0004505196856988L, 1.0484939797831612L, 1.4746279987544053L,
                 2.9676304579685473L, 6.184923457366982L};
    const W adm[]{4.28178520714217L, 418.45080980593247L, 2505.543413272741L,
                  5179.011131409537L, 7168.467673285191L};
    const W adl[]{4.286066992349312L, 460.2958907865258L, 4259.42380256366L,
                  15537.033394228609L, 35842.33836642596L};
    SoundHorizonRequest ar{
        {70, .3, 8.538168828261095e-5, .05, 5.046888425947121e-5},
        1059,
        "synthetic supplied drag; Astropy mapping compares background only"};
    auto ab = evaluate_early_late(ar, az,
                                  early_late_mask(EarlyLateOutput::e) |
                                      early_late_mask(EarlyLateOutput::dm_mpc) |
                                      early_late_mask(EarlyLateOutput::dl_mpc),
                                  policy());
    need(ab.rows.size() == 5 && !ab.ruler, "Astropy background only");
    for (size_t i = 0; i < 5; ++i) {
      const W directed_z_reference = dc(ar.model, az[i], 8192);
      near(directed_z_reference, dc(ar.model, az[i], 4096), .05e-9L,
           .05L * 2e-11L,
           "Astropy-coordinate independent refinement allocation");
      near(adm[i], directed_z_reference, .05e-9L, .05L * 2e-11L,
           "Astropy facts reference allocation checked independently");
      near(val(ab.rows[i], EarlyLateOutput::e), ae[i], 0, 2e-11L,
           "Astropy E named agreement");
      near(val(ab.rows[i], EarlyLateOutput::dm_mpc), adm[i], 1e-9L, 2e-11L,
           "Astropy DM named agreement");
      near(val(ab.rows[i], EarlyLateOutput::dl_mpc), adl[i], 1e-9L, 2e-11L,
           "Astropy DL named agreement");
    }
    SoundHorizonRequest source{
        {70, .3, 9e-5, .05, 5e-5}, 1059, "synthetic origin"};
    double z = .5;
    auto limited = evaluate_early_late(
        source, {&z, 1}, early_late_mask(EarlyLateOutput::e), policy());
    need(!limited.ruler && limited.callbacks == 0,
         "state only no ruler or quadrature");
    need(!limited.rows[0].outputs[unsigned(EarlyLateOutput::dm_mpc)].value,
         "unrequested absent");
    for (double bad : {-1., std::numeric_limits<double>::quiet_NaN(),
                       std::numeric_limits<double>::infinity()}) {
      auto failure = evaluate_early_late(source, {&bad, 1}, mask, policy());
      need(
          failure.rows.empty() ||
              !failure.rows[0].outputs[unsigned(EarlyLateOutput::dm_mpc)].value,
          "invalid redshift absent distance");
    }
    std::cout << "PASS " << checks << " early late independent peer controls\n";
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
