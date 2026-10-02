// Original bounded thermal/opacity ownership and physical-state controls.
// Separate peer tests own the independently authored stiff ODE reference.
#include "../src/hydrogen_rates.hpp"
#include "irred/recombination_drag.hpp"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace irred::cosmology;
using S = irred::numerics::Status;
namespace {
unsigned checks = 0;
void need(bool b, const char *s) {
  ++checks;
  if (!b)
    throw std::runtime_error(s);
}
double value(const HydrogenHistoryValue &v) {
  if (v.status != S::ok || !v.value)
    std::cerr << "group status=" << int(v.status)
              << " estimate=" << v.absolute_error_estimate << '\n';
  need(v.status == S::ok && v.value.has_value(),
       "accepted requested coordinate");
  return *v.value;
}
PureHydrogenRequest request() {
  return {{67.4, .0224, .12, 2.7255, 0, {}},
          1600,
          300,
          HydrogenTemperatureModel::evolved_compton_adiabatic};
}
} // namespace
int main() {
  try {
    // Original rational antiderivative: e(z)=2-2z integrates to3/4 on[0,1/2],
    // while a linear blend of full-cell cumulative error gives only1/2.
    need(
        detail::hydrogen_opacity_cell_integral(2, 0, 1, .5L) == .75L,
        "partial decreasing opacity-error integral preserves positive weights");
    need(detail::hydrogen_opacity_cell_integral(2, 0, 1, .5L) +
                 detail::hydrogen_opacity_cell_integral(1, 0, .5L, 1) ==
             1,
         "finite-cell primitive is additive without normalization");
    need(detail::hydrogen_opacity_cell_integral(0, 2, 1, .5L) == .25L,
         "partial increasing opacity-error rational control");
    auto source = request();
    auto owner = prepare_pure_hydrogen_history(source);
    std::cerr << "prepare=" << int(owner.status())
              << " work=" << owner.work().total()
              << " temperature_rhs=" << owner.work().temperature_rhs_evaluations
              << '\n';
    need(owner.status() == S::ok && owner.source() && owner.background(),
         "owned evolved state");
    need(owner.model_identity() == evolved_hydrogen_history_id &&
             !owner.method_identity().empty(),
         "explicit model and numerical identities");
    need(owner.work().equilibrium_solves == 1 &&
             owner.work().temperature_rhs_evaluations > 0 &&
             owner.work().temperature_rhs_evaluations ==
                 owner.work().rhs_evaluations &&
             owner.work().total() <= 4000000,
         "bounded actual coupled work");
    double zs[]{1600, 1599.999, 1599.99, 1599.9, 1500, 1300,
                1100, 1000,     800,     600,    300};
    auto batch = owner.evaluate(zs, 127);
    need(batch.status == S::ok && batch.rows.size() == 11 &&
             batch.requested == 127,
         "ordered coarse query");
    for (const auto &r : batch.rows) {
      double x = value(r.electron_fraction),
             T = value(r.matter_temperature_kelvin),
             tau = value(r.thomson_depth),
             q = value(r.thomson_opacity_per_redshift),
             g = value(r.visibility_per_redshift),
             survive = value(r.finite_endpoint_survival);
      need(x > 0 && x < 1 && T > 0 &&
               T <= source.model.tcmb_kelvin * (1 + r.redshift),
           "positive coupled physical state");
      need(tau >= 0 && q > 0 && g > 0 && survive > 0 && survive <= 1,
           "positive finite interval opacity");
      need(std::abs(survive - std::exp(-tau)) <= 4e-15 * survive,
           "survival exponential");
      need(std::abs(g - q * survive) <= 4e-15 * g,
           "visibility has per-redshift measure");
      std::cerr << r.redshift << " T=" << T << " T_error="
                << r.matter_temperature_kelvin.absolute_error_estimate
                << " tau=" << tau
                << " tau_error=" << r.thomson_depth.absolute_error_estimate
                << '\n';
    }
    need(value(batch.rows.front().matter_temperature_kelvin) ==
             source.model.tcmb_kelvin * 1601,
         "declared initial thermal boundary");
    need(value(batch.rows.back().thomson_depth) == 0 &&
             value(batch.rows.back().finite_endpoint_survival) == 1,
         "late zero depth and surviving boundary");
    auto before = owner.work().total();
    double bad[]{1000, std::numeric_limits<double>::quiet_NaN(), 299, 1601,
                 300};
    auto partial = owner.evaluate(bad, 127);
    need(partial.status == S::ok && partial.rows.size() == 5,
         "mixed row batch preserved");
    need(partial.rows[1].visibility_per_redshift.status == S::nonfinite_input &&
             partial.rows[2].matter_temperature_kelvin.status ==
                 S::outside_domain &&
             partial.rows[3].electron_fraction.status == S::outside_domain &&
             partial.rows[4].thomson_depth.status == S::ok,
         "typed partial failures");
    need(owner.work().total() == before,
         "queries reuse retained evolution/background");
    auto only = owner.evaluate(std::span(zs + 7, 1), 4);
    need(only.rows[0].electron_fraction.status == S::invalid_input &&
             !only.rows[0].thomson_depth.value &&
             only.rows[0].matter_temperature_kelvin.status == S::ok,
         "scalar masks skip projection");
    need(owner.evaluate(zs, 128).status == S::invalid_input &&
             owner.evaluate(zs, 0).status == S::invalid_input,
         "unknown and empty masks");
    need(owner.evaluate(zs, 127, 1).status == S::work_limit &&
             owner.evaluate(zs, 127, 4096, 1).status == S::work_limit,
         "query cap truthful refusal");
    auto cap = PureHydrogenPolicy{};
    cap.maximum_total_work = 1;
    auto failed = prepare_pure_hydrogen_history(source, cap);
    need(failed.status() == S::work_limit && failed.work().total() <= 1,
         "preparation actual work cap");
    cap = {};
    cap.maximum_native_bytes = 1;
    need(prepare_pure_hydrogen_history(source, cap).status() == S::work_limit,
         "payload before preparation");
    auto outside = source;
    outside.initial_redshift = 1600.001;
    need(prepare_pure_hydrogen_history(outside).status() == S::outside_domain,
         "bounded evolved qualification scope");
    cap = {};
    cap.absolute_temperature_tolerance_kelvin = 1e-30;
    cap.relative_temperature_tolerance = 0;
    auto strict = prepare_pure_hydrogen_history(source, cap)
                      .evaluate(std::span(zs + 7, 1), 5);
    need(strict.rows[0].matter_temperature_kelvin.status ==
                 S::conditioning_budget_exceeded &&
             !strict.rows[0].matter_temperature_kelvin.value &&
             strict.rows[0].electron_fraction.status == S::ok,
         "temperature refusal preserves accepted ion fraction");
    cap = {};
    cap.absolute_survival_tolerance = 0;
    cap.relative_survival_tolerance = 1e-30;
    strict = prepare_pure_hydrogen_history(source, cap)
                 .evaluate(std::span(zs + 7, 1), 65);
    need(strict.rows[0].finite_endpoint_survival.status ==
                 S::conditioning_budget_exceeded &&
             strict.rows[0].electron_fraction.status == S::ok,
         "exponential tau sensitivity refusal");
    auto copy = owner;
    source.model.tcmb_kelvin = 2.7;
    source.initial_redshift = 1550;
    need(copy.source()->initial_redshift == 1600 &&
             copy.source()->model.tcmb_kelvin == 2.7255,
         "owned source survives caller mutation");
    auto moved = std::move(copy);
    need(!copy.source() && copy.status() == S::invalid_input &&
             moved.status() == S::ok,
         "move leaves explicit invalid source");
    auto prescribed = request();
    prescribed.temperature_model =
        HydrogenTemperatureModel::prescribed_radiation;
    auto old = prepare_pure_hydrogen_history(prescribed);
    auto oldrows = old.evaluate(std::span(zs + 7, 1), 127);
    need(old.model_identity() == pure_hydrogen_history_id &&
             old.work().temperature_rhs_evaluations == 0,
         "distinct prescribed identity");
    need(value(oldrows.rows[0].matter_temperature_kelvin) ==
                 prescribed.model.tcmb_kelvin * 1001 &&
             oldrows.rows[0].visibility_per_redshift.status ==
                 S::outside_domain &&
             oldrows.rows[0].electron_fraction.status == S::ok,
         "prescribed temperature with unqualified photon groups refused");
    // Original composite Gauss-Legendre quadrature over the public per-z
    // measure. This checks conservation without renormalizing visibility or
    // sharing tau cells.
    constexpr double nodes[]{-.8611363115940525752, -.3399810435848562648,
                             .3399810435848562648, .8611363115940525752};
    constexpr double weights[]{.3478548451374538574, .6521451548625461426,
                               .6521451548625461426, .3478548451374538574};
    std::vector<double> queries;
    queries.reserve(2048);
    constexpr unsigned cells = 512;
    double width = 1300. / cells;
    for (unsigned i = 0; i < cells; ++i)
      for (double n : nodes)
        queries.push_back(300 + width * (i + .5) + width * n / 2);
    auto flux = owner.evaluate(queries, 32);
    need(flux.status == S::ok, "bounded visibility quadrature batch");
    long double mass = 0;
    for (size_t i = 0; i < flux.rows.size(); ++i)
      mass += weights[i % 4] * value(flux.rows[i].visibility_per_redshift) *
              width / 2;
    long double target = 1 - value(batch.rows.front().finite_endpoint_survival);
    std::cerr << "visibility mass residual=" << double(std::abs(mass - target))
              << '\n';
    need(std::abs(mass - target) <= 2e-7L,
         "finite visibility mass plus survival");
    std::cout << checks << " hydrogen thermal checks\n";
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
