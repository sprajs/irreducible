#include "irred/effective_fluid.hpp"
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
using namespace irred::cosmology;
using S = irred::numerics::Status;
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
      owner.source().species.size() != 1 || owner.source().species[0].mass_ev != .06 ||
      owner.source().species[0].temperature_today_ev != map.model->species[0].temperature_today_ev ||
      owner.source().species[0].statistical_weight != 2) throw std::runtime_error("complete relic source retained");
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
            << " prep=" << owner.preparation_callbacks() << " outer=" << got.outer_callbacks
            << " momentum=" << got.momentum_callbacks << '\n';
}
