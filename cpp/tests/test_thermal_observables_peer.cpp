// Original independent controls, not copied implementation/assets. FD
// reference: mpmath1.3.0 60/90 digits, direct momentum Gauss-Legendre32/48 per
// panel, q tails128/160. Distances use direct z; ruler uses u=sqrt(a), whereas
// production uses log(1+z) and direct a. Empirical reference refinements
// <=2.86e-9 of the fixed allocation; not rigorous universal bounds.
// CLASS v3.3.0 commit0ceb7a9a4c1e444ef5d5d56a8328a0640be91b18:
// public background_init/background_at_z, explicit FD
// temperature/state/constant matching; qmax80 bins4096/8192; fixed1e-12
// background integration. Its maximum refinement fraction <.00177. Sound
// horizon evaluated at supplied1059.95, never CLASS's predicted drag. Common
// equations/constants are shared ancestry. Mpc allocation1e-8+2e-10|ref|;
// ratios1e-10+5e-10|ref|; E/H2e-10rel. These named synthetic controls establish
// neither full Planck nor inference.
#include "irred/thermal_observables.hpp"
#include <algorithm>
#include <array>
#include <cfenv>
#include <cmath>
#include <iostream>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <utility>
namespace {
using namespace irred::cosmology;
using S = irred::numerics::Status;
using W = long double;
unsigned checks = 0;
ThermalObservablePolicy scientific_policy() {
  ThermalObservablePolicy p;
#ifdef IRRED_TEST_NESTED_CC
  p.thermal.momentum_method=ThermalMomentumMethod::nested_clenshaw_curtis;
#endif
  return p;
}
void need(bool b, const char *s) {
  ++checks;
  if (!b)
    throw std::runtime_error(s);
}
void near(W got, W ref, unsigned output, const char *s) {
  const W abs = output < 2   ? 0
                : output < 6 ? 1e-8L
                             : 1e-10L,
          rel = output < 6 ? 2e-10L : 5e-10L;
  need(std::abs(got - ref) <= abs + rel * std::abs(ref), s);
}
W value(const EarlyLateRow &r, unsigned i) {
  need(r.outputs[i].status == S::ok && r.outputs[i].value.has_value(),
       "requested independently checked output");
  return *r.outputs[i].value;
}
constexpr unsigned all = (1u << early_late_output_count) - 1;
constexpr unsigned mask(EarlyLateOutput x) { return early_late_mask(x); }
constexpr W light = 299792.458L;
// Independently derived natural/SI conversion; same defining constants are
// explicit shared ancestry, not independent measurements.
W critical(W H) {
  const W pi = std::numbers::pi_v<W>, c = 299792458.L, ev = 1.602176634e-19L,
          hbar = 6.62607015e-34L / (2 * pi),
          mpc = 648000000000.L * 149597870700.L / pi;
  const W hs = H * 1000 / mpc;
  return 3 * hs * hs * c * c / (8 * pi * 6.67430e-11L) / ev *
         std::pow(hbar * c / ev, 3);
}
W photon(W T, W H) {
  const W kb = 1.380649e-23L / 1.602176634e-19L;
  return std::numbers::pi_v<W> * std::numbers::pi_v<W> / 15 *
         std::pow(kb * T, 4) / critical(H);
}
struct Polynomial {
  W H, g, r, m, b, l;
};
Polynomial polynomial(const ThermalPhysicalModel &m) {
  W h2 = std::pow(W(m.h0_km_s_mpc) / 100, 2),
    g = photon(m.tcmb_kelvin, m.h0_km_s_mpc),
    r = g + W(m.physical_massless_nonphoton_density) / h2,
    b = W(m.physical_baryon_density) / h2,
    mat = b + W(m.physical_cdm_density) / h2;
  for (const auto &s : m.species) {
    need(s.mass_ev == 0, "analytic polynomial only massless species");
    const W kb = 1.380649e-23L / 1.602176634e-19L, pi = std::numbers::pi_v<W>;
    r += W(s.statistical_weight) * 7 * pi * pi / 240 *
         std::pow(kb * s.temperature_today_kelvin, 4) / critical(m.h0_km_s_mpc);
  }
  return {m.h0_km_s_mpc, g, r, mat, b, 1 - r - mat};
}
W P(Polynomial m, W a) { return m.r + m.m * a + m.l * a * a * a * a; }
W E(Polynomial m, W z) {
  W a = 1 / (1 + z);
  return std::sqrt(P(m, a)) / (a * a);
}
template <class F> W simpson(F f, W end, unsigned n) {
  W sum = f(0) + f(end);
  for (unsigned i = 1; i < n; ++i)
    sum += (i % 2 ? 4 : 2) * f(end * i / n);
  return sum * end / (3 * n);
}
W distance(Polynomial m, W z, unsigned n) {
  return light / m.H * simpson([&](W t) { return 1 / E(m, t); }, z, n);
}
W ruler(Polynomial m, W z, unsigned n) {
  const W u = 1 / std::sqrt(1 + z), B = 3 * m.b / (4 * m.g);
  return light / m.H / std::sqrt(3.L) *
         simpson(
             [&](W t) {
               return 2 * t / std::sqrt(P(m, t * t) * (1 + B * t * t));
             },
             u, n);
}
ThermalObservableRequest massive() {
  return {{67.4,
           .02237,
           .12,
           2.725501088928264,
           1.703989276e-5,
           {{.06, 1.9517599631877276, 2}}},
          1059.95,
          "synthetic supplied drag",
          "matched explicit CLASS constants/FD source"};
}
ThermalObservableRequest massless(bool species = false) {
  ThermalObservableRequest r{{70, .0245, .1225, 2.7255, 1e-5, {}},
                             1059.95,
                             "synthetic supplied drag",
                             "analytic massless synthetic source"};
  if (species)
    r.model.species.push_back({0, 1.9, 2});
  return r;
}

void frozen_references() {
  const double z[]{0, .1, 1, 10, 100, 1000};
  constexpr W ref[][9]{
      {1.0L, 67.400000000000005684341886080801486968994140625L,
       4447.9593175074180224967143013235740679294494208784079674923474950367101930486015681684408L,
       0.0L, 0.0L, 0.0L, 0.0L,
       31.2344337113419664373801703190901001266483840871950426660383345031989127926431849906550066L,
       0.0L},
      {1.0508317697470469758828886502139495483928467453966669421488227221161684860007351822155196L,
       70.8260612809509721477937390219757189322347231555228542580395188034845338348891265204734712L,
       4232.79867012216165350901249076826506775032893317655927732001808007959328940722113590403649L,
       434.112781390693674650668271241952780242256757408271290686199565316962299482363971667820113L,
       477.524059529763044525545124286224565747599629567476999847946291440292828482017132721527409L,
       430.471350231587317427908110075324341750494120764758346756914823645762581767069886134242033L,
       3.04842421562262741964605527603530369891996327367201662257287528565984957049521498797640134L,
       29.7235338810327596205383343914216768891389634536138633476476190552843252176670828122262332L,
       3.02285337919302175006282882309230937218474540525886084140386229574989586008503938068533778L},
      {1.79028099506720536641443960930127133573577613679043077518746346973421247661863719073640372L,
       120.66493906752965187290247778183790379592867575424462145962625335969897021288457671261577L,
       2484.50345491180620919869912720153472725463308422231031268133898056874613997721691892223955L,
       3401.29397008685395878728137850573159264202624104125296481101178768313136695288439513759216L,
       6802.58794017370791757456275701146318528405248208250592962202357536626273390576879027518432L,
       3063.20442777473352212771664734893147277679841618276260508387511631419712254974286964028735L,
       23.8845464757980689796700765015910639718009504300807168774140453259228452189617087589041495L,
       17.4466655220062015394759710599562292487011950165058252947423953870961312129398339569739476L,
       21.5104160838493437484767027864037581211333042395728507780758064035644645503060388743665683L},
      {20.5198761070776970170829739282123226485556946030658795982285519873100366498450883184597486L,
       1383.03964961703689559338369541192125900607280437719819732345189937835579290532083204089131L,
       216.7634587215286322947516603761513694557617060039110163690744690194638914023494946537876L,
       9631.25246285576943792799746314559750824527280379327651080646919311755259531156255947042233L,
       105943.777091413463817207972094601572590698000841726041618871161124293078548427188154174646L,
       5858.46539250452175494677413369562151354339295360146436203600396941894861487261680087139236L,
       67.6325242958488200624954033174923309775109595422009075575010835622995412157656022268492582L,
       1.52215508263076763653103629212654051896691006527065872941962063889697669325398352034999019L,
       41.139281160268459997731683214778889352147316759217902720356884862475101524747587470215408L},
      {578.348445108067350207242526872617745312703082346357369910004629002911357051910040748136405L,
       38980.6852002837426914984375887048856770731908467704342822029453666242274633079871913410651L,
       7.6907949785812844852972127250743535057225271806310769472653984517168084086285904218023813L,
       12817.5728195271664088650256141115758644407867880652253608927170380700534859615243509259768L,
       1294574.85477224380729536758702526916230851946559458776145016442084507540208211395944352366L,
       5017.96475939935217090062930633577779586411576622337884105833096323237316192204113013865605L,
       90.0074843301783745981128236838782290096304248368469618495607810389020537239205873230910098L,
       0.0540062551832496991413372230187018723791635308788480953871346970720996018497870327392990043L,
       35.2371225668360028389783986747969850646596369732923101199618031922081277136366304898698238L},
      {20477.937796413383064646103857784073441003811095579693658278215258963868718004934212794966L,
       1380213.00747826213496074695672442829926331935206869874276103103328358359335646405271977761L,
       0.217207384929475550402770729514275056539868095769219919369668533811648407217285693610731587L,
       13844.4942628161921692820469795741306247911024305129994896634695269245026653928340280936173L,
       13858338.7570790083614513290265537047554158935329435124891531329964514271680582268621217109L,
       3465.84864321470604560304142289956765539157432010925930311898122434097144279039281855050365L,
       97.2187260384638956789216671753772242507763293400046573244976115380506834502587692678878327L,
       0.00152527241863253116849860751968515857596671896421172848750298445706667545934890048015394571L,
       24.3378619210704350214654910791451724072926051070333610164692656348987772059652276829066098L},
  };
  constexpr W rs =
      142.405633430525686388027211133418785193859200471881493404261868025356164067826055195413854L;
  auto input = massive();
  auto mapping = map_thermal_physical_model(input.model);
  need(mapping.status == S::ok && mapping.model, "physical mapping admitted");
  near(mapping.model->omega_gamma,
       photon(input.model.tcmb_kelvin, input.model.h0_km_s_mpc), 0,
       "independent SI photon mapping");
  near(mapping.model->omega_b,
       W(input.model.physical_baryon_density) /
           std::pow(W(input.model.h0_km_s_mpc) / 100, 2),
       0, "physical baryon Omega h squared mapping");
  near(mapping.model->species[0].temperature_today_ev,
       W(input.model.species[0].temperature_today_kelvin) *
           (1.380649e-23L / 1.602176634e-19L),
       0, "explicit Kelvin eV conversion");
  auto p=scientific_policy();
  auto op = prepare_thermal_observables(input,p);
  need(op.status() == S::ok, "massive retained source admitted");
#ifdef IRRED_TEST_NESTED_CC
  need(op.background().momentum_method()==ThermalMomentumMethod::nested_clenshaw_curtis,
       "independent observable case actually retains opt-in CC");
#endif
  auto got = op.evaluate(z, all,p);
  need(got.status == S::ok && got.rows.size() == 6 && got.ruler &&
           got.ruler->status == S::ok && got.ruler->value,
       "frozen reference full batch");
  near(*got.ruler->value, rs, 3, "independent supplied ruler");
  for (unsigned i = 0; i < 6; ++i)
    for (unsigned k = 0; k < 9; ++k) {
      near(value(got.rows[i], k), ref[i][k], k,
           "highprecision direct-z/sqrt-a FD reference");
      const auto &v = got.rows[i].outputs[k];
      need(std::isfinite(v.error_estimate) && v.error_estimate >= 0,
           "finite diagnostic");
      const W allowance = (k < 2   ? 0
                           : k < 6 ? 1e-8L
                                   : 1e-10L) +
                          (k < 6 ? 2e-10L : 5e-10L) * std::abs(W(*v.value));
      need(v.error_estimate <= allowance,
           "public output diagnostic within frozen budget");
    }
  need(got.callbacks == got.outer_callbacks + got.momentum_callbacks &&
           got.outer_callbacks > 0 && got.momentum_callbacks > 0,
       "coarse plus nested callback accounting");
  constexpr W class_e[]{1L,
                        1.0508317697470424L,
                        1.7902809950671763L,
                        20.519876107077128L,
                        578.34844510805067L,
                        20477.937796412822L};
  constexpr W class_dm[]{0L,
                         434.11278136081603L,
                         3401.2939701017422L,
                         9631.2524628307419L,
                         12817.572819501434L,
                         13844.494262790595L};
  constexpr W class_dl[]{0L,
                         477.52405949688705L,
                         6802.5879402034834L,
                         105943.77709113804L,
                         1294574.8547696434L,
                         13858338.757053372L};
  for (unsigned i = 0; i < 6; ++i) {
    near(value(got.rows[i], 0), class_e[i], 0, "external CLASS E");
    near(value(got.rows[i], 1), 67.4L * class_e[i], 1, "external CLASS H");
    near(value(got.rows[i], 3), class_dm[i], 3, "external CLASS DM");
    near(value(got.rows[i], 4), class_dl[i], 4, "external CLASS DL");
  }
  near(*got.ruler->value, 142.40563343054799L, 3,
       "external CLASS supplied-redshift rs");
}
void analytic_and_scaling() {
  for (bool species : {false, true}) {
    auto req = massless(species);
    auto model = polynomial(req.model);
    auto op = prepare_thermal_observables(req);
    need(op.status() == S::ok, "massless analytic source admitted");
    const double z[]{0, 1e-8, .1, .7, 2, 4};
    auto got = op.evaluate(z, all);
    need(got.status == S::ok && got.rows.size() == 6 && got.ruler &&
             got.ruler->value,
         "analytic retained outputs");
    W rs = ruler(model, req.z_drag, 8192);
    near(rs, ruler(model, req.z_drag, 4096), 3, "ruler mesh refinement");
    need(std::abs(rs - ruler(model, req.z_drag, 4096)) <=
             .05L * (1e-8L + 2e-10L * std::abs(rs)),
         "ruler reference 5 percent allocation");
    near(*got.ruler->value, rs, 3, "analytic polynomial independent u ruler");
    for (unsigned i = 0; i < 6; ++i) {
      W ez = E(model, z[i]), dh = light / model.H / ez,
        dm = distance(model, z[i], 8192), dl = (1 + W(z[i])) * dm,
        dv = std::cbrt(dm * dm * z[i] * dh);
      need(std::abs(dm - distance(model, z[i], 4096)) <=
               .05L * (1e-8L + 2e-10L * std::abs(dm)),
           "distance reference 5 percent allocation");
      const W expected[]{ez, model.H * ez, dh,      dm,     dl,
                         dv, dm / rs,      dh / rs, dv / rs};
      for (unsigned k = 0; k < 9; ++k)
        near(value(got.rows[i], k), expected[k], k,
             "analytic massless independent reference");
    }
  }
  auto req = massless();
  const double z[]{.1, 1, 4};
  auto before = prepare_thermal_observables(req).evaluate(z, all);
  // Preserve all fractions: physical densities multiply by H0 squared; photon
  // temperature multiplies by sqrt(H0 ratio), rather than remaining fixed.
  req.model.h0_km_s_mpc *= 2;
  req.model.physical_baryon_density *= 4;
  req.model.physical_cdm_density *= 4;
  req.model.physical_massless_nonphoton_density *= 4;
  req.model.tcmb_kelvin *= std::sqrt(2.);
  auto fixedfractions = prepare_thermal_observables(req).evaluate(z, all);
  need(fixedfractions.rows.size() == 3, "fixed-fraction scaling admitted");
  for (unsigned i = 0; i < 3; ++i) {
    near(value(fixedfractions.rows[i], 0), value(before.rows[i], 0), 0,
         "fixed fractions E unchanged");
    near(2 * value(fixedfractions.rows[i], 3), value(before.rows[i], 3), 3,
         "fixed fractions inverse H0 DM");
    near(value(fixedfractions.rows[i], 6), value(before.rows[i], 6), 6,
         "fixed fractions H0 ratio cancellation");
  }
  req = massless();
  req.model.h0_km_s_mpc *= 2;
  auto fixedphysical = prepare_thermal_observables(req).evaluate(
      z, mask(EarlyLateOutput::e) | mask(EarlyLateOutput::dm_mpc));
  need(fixedphysical.rows.size() == 3,
       "fixed physical density scaling admitted");
  auto ref = polynomial(req.model);
  for (unsigned i = 0; i < 3; ++i) {
    near(value(fixedphysical.rows[i], 0), E(ref, z[i]), 0,
         "physical density scaling E reference");
    near(value(fixedphysical.rows[i], 3), distance(ref, z[i], 8192), 3,
         "physical density scaling distance reference");
  }
  need(std::abs(value(fixedphysical.rows[1], 0) - value(before.rows[1], 0)) >
           .1L,
       "fixed physical densities alter fractions and E");
  // Exact low-level states bypass physical-source fraction rounding. P=a4E2
  // has an analytic radiation endpoint even with an explicit massless species.
  auto radiation = prepare_thermal_background({70, 1, 0, 0, 0, {}});
  for (W a : {0.L, .001L, .5L, 1.L}) {
    auto p = radiation.scaled_expansion(a);
    need(p.status == S::ok, "exact radiation scaled domain");
    near(p.a4_e2, 1, 0, "exact radiation scaled polynomial");
  }
  auto dust = prepare_thermal_background({70, .25, 0, .25, .5, {}});
  for (W a : {0.L, .001L, .5L, 1.L}) {
    auto p = dust.scaled_expansion(a);
    need(p.status == S::ok, "Lambda zero scaled domain");
    near(p.a4_e2, .25L + .75L * a, 0, "analytic matter radiation polynomial");
  }
}
void zero_endpoint() {
  for (auto model : {ThermalFlatModel{70, 0, 0, 0, 0, {}},
                     ThermalFlatModel{70, 0, 0, .25, .75, {}}}) {
    auto background = prepare_thermal_background(model);
    need(background.status() == S::ok,
         "zero-radiation physical source admitted");
    auto scaled = background.scaled_expansion(0);
    need(scaled.status == S::ok && scaled.a4_e2 == 0 &&
             scaled.error_estimate == 0 && scaled.callbacks == 0,
         "analytic zero scaled endpoint is a valid coordinate");
    const double a[]{0};
    auto e = background.evaluate(a, thermal_e);
    need(e.rows.size() == 1 && !e.rows[0].e.value &&
             e.rows[0].e.status == S::outside_domain,
         "zero scaled coordinate does not admit E0 or its inverse");
  }
  auto vacuum = prepare_thermal_background({70, 0, 0, 0, 0, {}});
  // P(a)=a^4 is positive for every positive a; a numerical zero is not its
  // mathematical radiation endpoint. The wide arithmetic contract rejects
  // either subnormal coordinate or diagnostic rather than inventing a floor.
  for (W a : {1e-1230L, 1e-1235L, 1e-2000L}) {
    auto p = vacuum.scaled_expansion(a);
    need(p.status == S::outside_domain,
         "positive wide P or diagnostic underflow cannot become exact zero");
  }
  const W admitted_a = 1e-1200L;
  auto wide = vacuum.scaled_expansion(admitted_a);
  need(wide.status == S::ok && std::isnormal(wide.a4_e2) &&
           std::isnormal(wide.error_estimate),
       "representable wide coordinate and diagnostic remain available");
  near(std::sqrt(wide.a4_e2) / (admitted_a * admitted_a), 1, 0,
       "pure Lambda analytic E reconstructed from admitted wide state");
  auto req = massless();
  req.model.tcmb_kelvin = 0;
  need(prepare_thermal_observables(req).status() != S::ok,
       "ruler/distance consumer requires a strictly positive photon source");
  ThermalFlatModel relic{70, 5e-5, 2e-5, .05, .25, {{.06, .000168, 2}}};
  auto background = prepare_thermal_background(relic);
  auto at_zero = background.scaled_expansion(0);
  const W pi = std::numbers::pi_v<W>;
  const W expected =
      W(relic.omega_gamma) + relic.omega_massless_nonphoton +
      W(relic.species[0].statistical_weight) * 7 * pi * pi / 240 *
          std::pow(W(relic.species[0].temperature_today_ev), 4) / critical(70);
  need(at_zero.status == S::ok && at_zero.a4_e2 > 0 && at_zero.callbacks == 0,
       "analytic FD radiation endpoint has no momentum callbacks");
  near(at_zero.a4_e2, expected, 0, "independent FD scaled origin convention");
}
void boundaries_and_lifetimes() {
  auto req = massive();
  auto op = prepare_thermal_observables(req);
  const double z[]{0, .1, 1};
  req.model.species.clear();
  req.drag_origin.clear();
  need(op.source().model.species.size() == 1 &&
           !op.source().drag_origin.empty(),
       "input acquisition owns physical species/origins");
  auto copy = op;
  need(copy.source().model.species.data() != op.source().model.species.data(),
       "copies own distinct physical species");
  auto moved = std::move(op);
  need(op.status() != S::ok && op.evaluate(z, all).rows.empty(),
       "move invalidates source");
  need(moved.status() == S::ok, "move target retains admission");
  auto *alias = &moved;
  moved = std::move(*alias);
  need(moved.status() == S::ok, "self move preserves owner");
  ThermalObservables assigned;
  assigned = std::move(moved);
  need(moved.status() != S::ok && assigned.status() == S::ok,
       "move assignment transfers complete owner");
  auto got = assigned.evaluate(z, mask(EarlyLateOutput::e));
  need(got.rows.size() == 3 && !got.ruler && got.outer_callbacks == 0,
       "state only omits outer/ruler work");
  near(value(got.rows[2], 0), 1.7902809950671473L, 0,
       "moved retained independent E");
  for (const auto &row : got.rows)
    for (unsigned i = 1; i < 9; ++i)
      need(!row.outputs[i].value, "unrequested values absent");
  auto generic = prepare_thermal_observables(massless());
  ThermalObservablePolicy limited;
  limited.maximum_callbacks_per_point = 4;
  auto partial = generic.evaluate(
      z, mask(EarlyLateOutput::e) | mask(EarlyLateOutput::dm_mpc), limited);
  need(partial.rows.size() == 3 && partial.rows[2].outputs[0].value &&
           !partial.rows[2].outputs[3].value,
       "distance callback failure preserves admitted E");
  limited = {};
  limited.maximum_total_callbacks = 0;
  auto nointegral = generic.evaluate(z, mask(EarlyLateOutput::e), limited);
  need(nointegral.rows.size() == 3 && nointegral.rows[2].outputs[0].value,
       "zero total quota permits no quadrature state");
  auto refused = generic.evaluate(z, mask(EarlyLateOutput::dm_mpc), limited);
  need(refused.rows.empty() || !refused.rows.back().outputs[3].value,
       "zero callback budget withholds integral");
  limited = {};
  limited.maximum_native_bytes = 0;
  need(generic.evaluate(z, all, limited).rows.empty(),
       "payload admission before rows");
  limited = {};
  limited.maximum_points = 2;
  need(generic.evaluate(z, all, limited).rows.empty(),
       "point admission before rows");
  limited = {};
  limited.absolute_tolerance_mpc = 0;
  limited.relative_tolerance = 1e-30;
  auto impossible = generic.evaluate(z, mask(EarlyLateOutput::dm_mpc), limited);
  need(impossible.rows.empty() || !impossible.rows.back().outputs[3].value,
       "unattainable numerical budget refused");
  for (double bad : {-1., std::numeric_limits<double>::infinity(),
                     std::numeric_limits<double>::quiet_NaN()}) {
    auto r = generic.evaluate({&bad, 1}, all);
    need(r.rows.empty() || !r.rows[0].outputs[3].value,
         "invalid redshift withholds usable distance");
  }
  for (unsigned mode = 0; mode < 7; ++mode) {
    auto r = massless();
    if (mode == 0)
      r.model.h0_km_s_mpc = 0;
    if (mode == 1)
      r.model.tcmb_kelvin = 0;
    if (mode == 2)
      r.model.physical_baryon_density = -1;
    if (mode == 3)
      r.model.physical_cdm_density = 2;
    if (mode == 4)
      r.z_drag = -1;
    if (mode == 5)
      r.drag_origin.clear();
    if (mode == 6)
      r.source_origin.clear();
    need(prepare_thermal_observables(r).status() != S::ok,
         "invalid physical/source/drag admission");
  }
  auto extreme_source = massless();
  extreme_source.model = {1e100, 0, 0, 1e50, 0, {}};
  auto extreme_owner = prepare_thermal_observables(extreme_source);
  need(extreme_owner.status() == S::ok,
       "independent extreme physical source admitted");
  const double extreme_z[]{1e105};
  auto extreme = extreme_owner.evaluate(
      extreme_z, mask(EarlyLateOutput::e) | mask(EarlyLateOutput::h_km_s_mpc) |
                     mask(EarlyLateOutput::dh_mpc));
  need(extreme.rows.size() == 1, "extreme output row retained");
  near(value(extreme.rows[0], 0),
       E(polynomial(extreme_source.model), extreme_z[0]), 0,
       "independent extreme dimensionless expansion");
  need(!extreme.rows[0].outputs[1].value && !extreme.rows[0].outputs[2].value,
       "dimensionless E survives dimensional overflow/diagnostic underflow");
  need(extreme.rows[0].outputs[1].status == S::outside_domain &&
           extreme.rows[0].outputs[2].status == S::conditioning_budget_exceeded,
       "representable DH with unrepresentable positive diagnostic refused");
  const double value_underflow_z[]{1e110};
  auto value_underflow = extreme_owner.evaluate(
      value_underflow_z, mask(EarlyLateOutput::e) |
                             mask(EarlyLateOutput::h_km_s_mpc) |
                             mask(EarlyLateOutput::dh_mpc));
  need(value_underflow.rows.size() == 1, "underflow row retained");
  near(value(value_underflow.rows[0], 0),
       E(polynomial(extreme_source.model), value_underflow_z[0]), 0,
       "independent E at dimensional value underflow");
  need(value_underflow.rows[0].outputs[1].status == S::outside_domain &&
           !value_underflow.rows[0].outputs[1].value &&
           value_underflow.rows[0].outputs[2].status == S::outside_domain &&
           !value_underflow.rows[0].outputs[2].value,
       "unrepresentable dimensional values separately refused");
  limited = {};
  limited.maximum_total_callbacks = 0;
  auto failed_ruler = generic.evaluate(
      z, mask(EarlyLateOutput::e) | mask(EarlyLateOutput::dm_over_rs), limited);
  need(failed_ruler.ruler && !failed_ruler.ruler->value &&
           failed_ruler.rows[0].outputs[0].value &&
           !failed_ruler.rows[0].outputs[6].value,
       "failed supplied ruler preserves independent E");
  auto tiny_source = massless();
  tiny_source.model = {
      std::numeric_limits<double>::denorm_min(), 0, 0, 1e-162, 0, {}};
  auto tiny_owner = prepare_thermal_observables(tiny_source);
  need(tiny_owner.status() == S::ok,
       "explicit tiny H0 source maps to finite positive photon fraction");
  const double present[]{0};
  auto tiny_present = tiny_owner.evaluate(
      present, mask(EarlyLateOutput::e) | mask(EarlyLateOutput::h_km_s_mpc) |
                   mask(EarlyLateOutput::dh_mpc));
  need(tiny_present.rows.size() == 1 && value(tiny_present.rows[0], 0) == 1,
       "tiny source retains exact dimensionless present normalization");
  need(tiny_present.rows[0].outputs[1].status == S::outside_domain &&
           !tiny_present.rows[0].outputs[1].value &&
           tiny_present.rows[0].outputs[2].status == S::outside_domain &&
           !tiny_present.rows[0].outputs[2].value,
       "exact present identity cannot expose positive subnormal H");
  tiny_source.model.h0_km_s_mpc = std::numeric_limits<double>::min();
  tiny_source.model.tcmb_kelvin = 1e-154;
  auto normal_present =
      prepare_thermal_observables(tiny_source)
          .evaluate(present, mask(EarlyLateOutput::e) |
                                 mask(EarlyLateOutput::h_km_s_mpc));
  need(normal_present.rows.size() == 1 &&
           value(normal_present.rows[0], 0) == 1 &&
           value(normal_present.rows[0], 1) ==
               std::numeric_limits<double>::min(),
       "adjacent representable normal H identity remains available");
  const auto prior = std::fegetround();
  need(std::fesetround(FE_DOWNWARD) == 0, "set hostile rounding");
  auto rounding = generic.evaluate(z, all);
  need(std::fesetround(prior) == 0, "restore rounding");
  need(rounding.status != S::ok && rounding.rows.empty(),
       "unsupported arithmetic refuses processing");
}
} // namespace
int main() {
  try {
    need(std::numeric_limits<W>::digits >= 64, "wide peer arithmetic");
    frozen_references();
    analytic_and_scaling();
    zero_endpoint();
    boundaries_and_lifetimes();
    std::cout << "PASS " << checks << " thermal observables peer controls\n";
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
