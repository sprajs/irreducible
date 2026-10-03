#include "irred/effective_fluid.hpp"
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace irred::cosmology;
using S = irred::numerics::Status;
namespace {
bool same_source(const ThermalFlatModel &a, const ThermalFlatModel &b) {
  if (a.h0_km_s_mpc != b.h0_km_s_mpc || a.omega_gamma != b.omega_gamma ||
      a.omega_massless_nonphoton != b.omega_massless_nonphoton || a.omega_b != b.omega_b ||
      a.omega_cdm != b.omega_cdm || a.species.size() != b.species.size()) return false;
  for (std::size_t i = 0; i < a.species.size(); ++i)
    if (a.species[i].mass_ev != b.species[i].mass_ev ||
        a.species[i].temperature_today_ev != b.species[i].temperature_today_ev ||
        a.species[i].statistical_weight != b.species[i].statistical_weight) return false;
  return true;
}
}
int main() {
  // Explicit synthetic physical source. This is not an external CLASS fit.
  const ThermalPhysicalModel physical{67.4, .0224, .12, 2.7255, 1.68e-5, {{.06, 1.95, 2}}};
  ThermalObservablePolicy p;
  const auto map = map_thermal_physical_model(physical, p.thermal);
  if (map.status != S::ok || !map.model || !map.scalar_witnesses) throw std::runtime_error("complete physical map");
  const EffectiveFluidRequest request{*map.model, 1e9, 1. / 3500, 1. / 1060,
      "synthetic-complete-positive-mass-FD-physical-source", "fixed-synthetic-scale-not-predicted-drag"};
  auto zero_request = request; zero_request.density_at_transition = 0;
  const auto zero = prepare_effective_fluid(zero_request, p), owner = prepare_effective_fluid(request, p);
  if (owner.status() != S::ok || zero.status() != S::ok ||
      !same_source(owner.source(), *map.model) || !same_source(zero.source(), *map.model) ||
      owner.source().species.size() != 1 || owner.source().species[0].mass_ev != .06 ||
      owner.source().species[0].temperature_today_ev != map.model->species[0].temperature_today_ev ||
      owner.source().species[0].statistical_weight != 2) throw std::runtime_error("complete relic source retained");
  auto failed_request = request;
  failed_request.density_at_transition = std::numeric_limits<double>::max();
  const auto failed = prepare_effective_fluid(failed_request, p);
  const auto refused = failed.evaluate({}, effective_fluid_ruler, p);
  if (failed.status() != S::outside_domain || failed.reference_background().status() != S::ok ||
      failed.request().density_at_transition != failed_request.density_at_transition ||
      !same_source(failed.source(), *map.model) || !same_source(failed.reference_background().source(), *map.model) ||
      failed.source().species.size() != 1 || failed.source().species[0].mass_ev != .06 ||
      failed.source().species[0].temperature_today_ev != map.model->species[0].temperature_today_ev ||
      failed.source().species[0].statistical_weight != 2 || failed.preparation_callbacks() == 0 ||
      failed.preparation_callbacks() != owner.preparation_callbacks() ||
      failed.preparation_callbacks() > p.thermal.maximum_total_callbacks ||
      refused.status != S::outside_domain || !refused.rows.empty() || refused.ruler || refused.callbacks != 0)
    throw std::runtime_error("late closure refusal retains complete earned positive-FD normalization prefix");
  const std::array<double, 4> a{1, request.ruler_endpoint, request.transition_scale, request.ruler_endpoint};
  const auto baseline = zero.evaluate(a, thermal_e | thermal_h | effective_fluid_ruler, p),
             got = owner.evaluate(a, thermal_e | thermal_h | effective_fluid_ruler, p);
  if (!baseline.ruler || !baseline.ruler->value || !got.ruler || !got.ruler->value ||
      !(got.ruler->status == S::ok && baseline.ruler->status == S::ok) ||
      !(*got.ruler->value < *baseline.ruler->value)) throw std::runtime_error("matched massive-FD finite ruler response");
  for (std::size_t i = 0; i < a.size(); ++i) {
    if (!got.rows[i].e.value || !baseline.rows[i].e.value ||
        *got.rows[i].e.value < *baseline.rows[i].e.value ||
        got.rows[i].scale_factor != a[i]) throw std::runtime_error("matched massive-FD ordered expansion");
  }
  if (got.callbacks != got.outer_callbacks + got.momentum_callbacks || got.momentum_callbacks == 0 ||
      got.callbacks > p.maximum_total_callbacks || owner.preparation_callbacks() == 0 ||
      got.rows[1].e.value != got.rows[3].e.value) throw std::runtime_error("actual nested work/order");
  std::cout << effective_fluid_model_id << '\n' << effective_fluid_ruler_id << '\n'
            << std::hexfloat << "a_end=" << request.ruler_endpoint
            << " Xc=" << request.density_at_transition << " ac=" << request.transition_scale
            << " rs_ref=" << *baseline.ruler->value << " rs_new=" << *got.ruler->value
            << " estimate=" << got.ruler->error_estimate << std::defaultfloat
            << " prep=" << owner.preparation_callbacks() << " refused_prep=" << failed.preparation_callbacks()
            << " outer=" << got.outer_callbacks
            << " momentum=" << got.momentum_callbacks << '\n';
}
