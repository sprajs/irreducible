// Original ownership, admission and finite-measure controls. Independent
// direct-momentum/Radau mathematics lives in the separately authored peer.
#include "irred/recombination_drag.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace irred::cosmology;
using S = irred::numerics::Status;
namespace {
unsigned checks = 0;
void need(bool condition, const char *message) {
  ++checks;
  if (!condition) throw std::runtime_error(message);
}
double value(const HydrogenHistoryValue &v) {
  need(v.status == S::ok && v.value.has_value(), "accepted scalar group");
  return *v.value;
}
auto groups(const HydrogenHistoryRow &r) {
  return std::array{r.electron_fraction, r.drag_depth,
                   r.matter_temperature_kelvin, r.thomson_depth,
                   r.thomson_opacity_per_redshift, r.visibility_per_redshift,
                   r.finite_endpoint_survival};
}
PureHydrogenRequest request() {
  return {{67.4, .0224, .12, 2.7255, 0,
           {{.06, 1.95, 2}, {0, 1.95, 2}, {0, 1.95, 2}}},
          1600, 300, HydrogenTemperatureModel::evolved_compton_adiabatic};
}
PureHydrogenPolicy extended() {
  PureHydrogenPolicy p;
  p.maximum_total_work = 700000000; // Explicit caller choice, not a default.
  return p;
}
PureHydrogenHistory prepare(const PureHydrogenRequest &r) {
  auto h = prepare_pure_hydrogen_history(r, extended());
  auto w = h.work();
  std::cerr << "mass-profile";
  for (auto s : r.model.species) std::cerr << ' ' << s.mass_ev;
  std::cerr << " status=" << int(h.status()) << " total=" << w.total()
            << " momentum=" << w.momentum_callbacks
            << " background=" << w.background_evaluations
            << " rhs=" << w.rhs_evaluations
            << " temperature_rhs=" << w.temperature_rhs_evaluations << '\n';
  need(h.status() == S::ok && h.source() && h.background(), "retained relic history");
  need(w.total() <= extended().maximum_total_work &&
           w.momentum_callbacks > 4000000 && w.background_evaluations == 32769 &&
           w.equilibrium_solves == 1,
       "nested work charged without a hidden cap increase");
  return h;
}
}
int main() {
  try {
    auto r = request();
    auto refused = prepare_pure_hydrogen_history(r);
    need(PureHydrogenPolicy{}.maximum_total_work == 4000000,
         "unchanged combined default cap");
    need(refused.status() == S::work_limit &&
             refused.work().total() <= 4000000 &&
             refused.work().momentum_callbacks > 0 &&
             refused.work().rhs_evaluations == 0,
         "actual nested-provider refusal preserved");
    // Domain controls use a deliberately tiny caller cap after admission, so
    // admitted endpoints are distinguished from physical-domain refusals.
    auto cap = extended(); cap.maximum_total_work = 1;
    for (auto s : {ThermalPhysicalSpecies{.001, 1.8, 1},
                   ThermalPhysicalSpecies{.3, 2, 2}}) {
      auto edge = r; edge.model.species = {s};
      need(prepare_pure_hydrogen_history(edge, cap).status() == S::work_limit,
           "positive-profile inclusive domain endpoints");
    }
    for (auto s : {ThermalPhysicalSpecies{-.001, 1.95, 2},
                   ThermalPhysicalSpecies{.300001, 1.95, 2},
                   ThermalPhysicalSpecies{.06, 1.799, 2},
                   ThermalPhysicalSpecies{.06, 2.001, 2},
                   ThermalPhysicalSpecies{.06, 1.95, .999},
                   ThermalPhysicalSpecies{.06, 1.95, 2.001}}) {
      auto bad = r; bad.model.species = {s};
      need(prepare_pure_hydrogen_history(bad, cap).status() == S::outside_domain,
           "bounded explicit massive profile");
    }
    for (unsigned field = 0; field != 3; ++field) {
      auto bad = r; auto &s = bad.model.species[0];
      if (field == 0) s.mass_ev = std::numeric_limits<double>::quiet_NaN();
      if (field == 1) s.temperature_today_kelvin = std::numeric_limits<double>::infinity();
      if (field == 2) s.statistical_weight = std::numeric_limits<double>::quiet_NaN();
      need(prepare_pure_hydrogen_history(bad, cap).status() == S::nonfinite_input,
           "typed nonfinite species");
    }
    auto bad = r; bad.model.species.assign(4, {.06, 1.95, 1});
    need(prepare_pure_hydrogen_history(bad, cap).status() == S::outside_domain,
         "at most three positive species");
    bad.model.species.assign(17, {0, 1.95, 2});
    need(prepare_pure_hydrogen_history(bad, cap).status() == S::work_limit,
         "unchanged total species bound");
    // Massless temperatures and populated-state weights retain provider
    // admission; the new massive-profile bounds are not applied to them.
    bad = r; bad.model.species = {{0, 1, .5}};
    need(prepare_pure_hydrogen_history(bad, cap).status() == S::work_limit,
         "massless admission remains distinct");
    auto h = prepare(r);
    need(h.model_identity() == evolved_hydrogen_relic_history_id,
         "distinct evolved positive-mass identity");
    double z[]{1600, 1599.999, 1500, 1100, 1000, 600, 300};
    auto b = h.evaluate(z, 127);
    need(b.status == S::ok && b.rows.size() == 7, "ordered seven-group batch");
    for (const auto &row : b.rows) {
      for (auto v : groups(row)) need(value(v) >= 0, "nonnegative scalar");
      double x = value(row.electron_fraction), T = value(row.matter_temperature_kelvin),
             tau = value(row.thomson_depth), q = value(row.thomson_opacity_per_redshift),
             g = value(row.visibility_per_redshift), Sboundary = value(row.finite_endpoint_survival);
      need(x > 0 && x < 1 && T > 0 && T <= r.model.tcmb_kelvin * (1 + row.redshift),
           "physical evolved state");
      need(q > 0 && g > 0 && Sboundary > 0 && Sboundary <= 1,
           "positive visibility and survival");
      need(std::abs(Sboundary - std::exp(-tau)) <= 4e-15 * Sboundary &&
               std::abs(g - q * Sboundary) <= 4e-15 * g,
           "finite-endpoint exponential projection");
    }
    need(value(b.rows.back().thomson_depth) == 0 &&
             value(b.rows.back().finite_endpoint_survival) == 1,
         "exact finite late endpoint");
    need(h.conditional_unit_depth_redshift().status == S::ok,
         "conditional unit drag depth accepted separately");
    auto massless = r;
    for (auto &s : massless.model.species) s.mass_ev = 0;
    auto old = prepare_pure_hydrogen_history(massless);
    need(old.status() == S::ok && old.model_identity() == evolved_hydrogen_history_id,
         "unchanged legacy massless identity and default usability");
    auto initial = old.evaluate(std::span(z, 1), 5);
    need(value(initial.rows[0].electron_fraction) == value(b.rows[0].electron_fraction) &&
             value(initial.rows[0].matter_temperature_kelvin) == value(b.rows[0].matter_temperature_kelvin),
         "initial Saha and radiation boundary do not depend on relic H");
    need(old.method_identity() == h.method_identity(), "unchanged ODE arithmetic method identity");
    auto before = h.work().total();
    r.model.species.clear(); r.model.physical_baryon_density = .03;
    need(h.source()->model.species.size() == 3 && h.source()->model.physical_baryon_density == .0224,
         "immutable source lifetime");
    double mixed[]{1000, 299, std::numeric_limits<double>::quiet_NaN(), 300};
    auto m = h.evaluate(mixed, 127);
    need(m.rows[0].visibility_per_redshift.status == S::ok &&
             m.rows[1].visibility_per_redshift.status == S::outside_domain &&
             m.rows[2].visibility_per_redshift.status == S::nonfinite_input &&
             m.rows[3].finite_endpoint_survival.status == S::ok,
         "typed partial rows");
    auto mask = h.evaluate(std::span(z + 4, 1), 4);
    need(mask.rows[0].matter_temperature_kelvin.status == S::ok &&
             !mask.rows[0].electron_fraction.value && !mask.rows[0].visibility_per_redshift.value,
         "unrequested projections skipped");
    need(h.evaluate(z, 128).status == S::invalid_input &&
             h.evaluate(z, 127, 1).status == S::work_limit &&
             h.evaluate(z, 127, 4096, 1).status == S::work_limit,
         "mask and query resource refusal");
    need(h.work().total() == before, "queries perform no new provider work");
    auto copy = h;
    auto moved = std::move(copy);
    need(!copy.source() && moved.source() &&
             value(moved.evaluate(std::span(z + 4, 1), 1).rows[0].electron_fraction) ==
                 value(b.rows[4].electron_fraction),
         "copy/move own complete relic history");
    auto split = request(); split.model.species = {{.06,1.95,1},{.06,1.95,1},{0,1.95,2},{0,1.95,2}};
    auto sh = prepare(split); auto sb = sh.evaluate(z,127);
    for (std::size_t i = 0; i != b.rows.size(); ++i) {
      auto a = groups(b.rows[i]), c = groups(sb.rows[i]);
      for (unsigned j = 0; j != 7; ++j)
        need(std::abs(value(a[j])-value(c[j])) <= 2e-12 * std::max(std::abs(value(a[j])), 1e-20),
             "split populated-state weights preserve all seven scalars");
    }
    need(sh.work().momentum_callbacks > h.work().momentum_callbacks,
         "split species count actual additional momentum work");
    auto three = request(); three.model.species = {{.001,1.8,1},{.1,1.95,2},{.3,2,2}};
    auto th = prepare(three);
    for (auto row : th.evaluate(z,127).rows)
      for (auto v : groups(row)) value(v);
    auto prescribed = request(); prescribed.model.species = {{.3,2,2}};
    prescribed.temperature_model = HydrogenTemperatureModel::prescribed_radiation;
    auto ph = prepare(prescribed);
    need(ph.model_identity() == pure_hydrogen_relic_history_id,
         "distinct prescribed positive-mass identity");
    auto pb = ph.evaluate(std::span(z + 4, 1),127);
    need(pb.rows[0].electron_fraction.status == S::ok && pb.rows[0].drag_depth.status == S::ok &&
             value(pb.rows[0].matter_temperature_kelvin) == prescribed.model.tcmb_kelvin * 1001 &&
             pb.rows[0].visibility_per_redshift.status == S::outside_domain,
         "prescribed scientific boundary retains accepted legacy groups");
    // Original composite GL4 control of the public finite-interval measure.
    constexpr long double qx[]{.3399810435848562648L,.8611363115940525752L};
    constexpr long double qw[]{.6521451548625461426L,.3478548451374538574L};
    std::vector<double> queries; queries.reserve(8192);
    const long double dz = 1300.L/2048;
    for(unsigned i=0;i!=2048;++i)
      for(unsigned j=0;j!=2;++j)
        for(int sign : {-1,1}) queries.push_back((double)(300+(i+.5L)*dz+sign*qx[j]*dz/2));
    long double integral = 0;
    for(unsigned part=0;part!=2;++part) {
      auto flux = h.evaluate(std::span(queries).subspan(4096*part,4096),32);
      need(flux.status == S::ok && flux.rows.size() == 4096,
           "visibility quadrature streams admitted query batches");
      for(unsigned i=0;i!=1024;++i)
        for(unsigned j=0;j!=2;++j)
          for(unsigned k=0;k!=2;++k)
            integral += dz/2*qw[j]*value(flux.rows[4*i+2*j+k].visibility_per_redshift);
    }
    long double residual = std::abs(integral + value(b.rows.front().finite_endpoint_survival)-1);
    std::cerr << "visibility mass residual=" << (double)residual << '\n';
    need(residual <= 2e-7L, "visibility mass plus finite boundary without normalization");
    std::cout << "hydrogen relic owner " << checks << " checks passed\n";
  } catch(const std::exception &e) { std::cerr << e.what() << '\n'; return 1; }
}
