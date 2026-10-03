#include "irred/effective_fluid.hpp"
#include "irred/quantities.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace irred::cosmology;
using S = irred::numerics::Status;
namespace {
using W = long double;
unsigned checks = 0;
void require(bool x, const char *why) { ++checks; if (!x) throw std::runtime_error(why); }
// Independently authored log-density and fixed direct-a composite Simpson.
// Empty-species dyadic reference is independently polynomial, not an FD oracle.
W density(W xc, W ac, W a) {
  return a == 0 ? 2 * xc : 2 * xc / (1 + std::exp(4.5L * (std::log(a) - std::log(ac))));
}
W reference(const EffectiveFluidRequest &r, unsigned n) {
  const auto &m = r.reference;
  const W lambda = 1.L - m.omega_gamma - m.omega_massless_nonphoton - m.omega_b - m.omega_cdm;
  const W today = density(r.density_at_transition, r.transition_scale, 1),
      step = W(r.ruler_endpoint) / n, loading = 3.L * m.omega_b / (4.L * m.omega_gamma);
  W sum = 0;
  for (unsigned i = 0; i <= n; ++i) {
    const W a = step * i, a2 = a * a;
    const W p = m.omega_gamma + W(m.omega_massless_nonphoton) +
        (W(m.omega_b) + m.omega_cdm) * a +
        (lambda + density(r.density_at_transition, r.transition_scale, a) - today) * a2 * a2;
    const W f = 1 / std::sqrt(3 * p * (1 + loading * a));
    sum += (i == 0 || i == n ? 1 : (i & 1 ? 4 : 2)) * f;
  }
  return W(irred::speed_of_light_m_per_s) / 1000 / m.h0_km_s_mpc * step * sum / 3;
}
}
int main() {
  for (double xc : {0., .25, .5}) for (double endpoint : {.25, .5, 1.}) {
    const EffectiveFluidRequest r{{70, .5, .125, .0625, .0625, {}}, xc, .5, endpoint,
      "independent-dyadic-polynomial-control", "matched-finite-endpoint"};
    const auto owner = prepare_effective_fluid(r);
    require(owner.status() == S::ok, "peer owner domain");
    const auto got = owner.evaluate({}, effective_fluid_ruler);
    require(got.ruler && got.ruler->status == S::ok && got.ruler->value, "peer ruler available");
    const W coarse = reference(r, 4096), fine = reference(r, 8192), allocation = 1e-8L + 2e-10L * std::abs(fine);
    require(std::abs(coarse - fine) <= .05L * allocation, "EACH reference refinement <=5% frozen ruler allocation");
    require(std::abs(W(*got.ruler->value) - fine) <= allocation, "independent direct-a ruler comparison");
    for (W a : {.125L, .25L, .5L, .75L, 1.L}) {
      const auto state = owner.fluid_state(a);
      const W expected = density(xc, .5L, a);
      require(state.status == S::ok && std::abs(state.density - expected) <= 2e-15L * (1 + expected),
              "independent log-density arithmetic");
      const W da = 0x1p-18L,
          derivative = (density(xc, .5L, a * std::exp(da)) - density(xc, .5L, a * std::exp(-da))) / (2 * da);
      require(std::abs(state.log_scale_derivative - derivative) <= 2e-10L * (1 + std::abs(derivative)),
              "independent symmetric log-step continuity derivative");
    }
  }
  std::cout << "effective_fluid peer checks=" << checks << " direct_reference_nodes=110610\n";
}
