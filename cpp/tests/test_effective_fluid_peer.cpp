#include "irred/effective_fluid.hpp"
#include "irred/quantities.hpp"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace irred::cosmology;
using S = irred::numerics::Status;
namespace {
using W = long double;
unsigned checks = 0;
void require(bool x, const char *why) { ++checks; if (!x) throw std::runtime_error(why); }
// Independently authored log-density and fixed direct-a composite Simpson.
// Empty-species dyadic reference is independently polynomial, not an FD oracle.
constexpr W operation_floor = 128 * std::numeric_limits<W>::epsilon();
struct Density { W value, estimate; };
Density density_record(W xc, W ac, W a) {
  if (a == 0) return {2 * xc, operation_floor * std::abs(2 * xc)};
  const W la = std::log(a), lc = std::log(ac);
  const W value = 2 * xc / (1 + std::exp(4.5L * (la - lc)));
  // Conditional empirical log/exp/operator allowance; no rigorous libm claim.
  return {value, operation_floor * (1 + 9 * (std::abs(la) + std::abs(lc))) * std::abs(value)};
}
W density(W xc, W ac, W a) { return density_record(xc, ac, a).value; }
struct Reference { W value, arithmetic_estimate, reporting_cast_loss; };
Reference reference(const EffectiveFluidRequest &r, unsigned n) {
  const auto &m = r.reference;
  require((n == 4096 || n == 8192) && r.transition_scale == .5 &&
          (r.density_at_transition == 0 || r.density_at_transition == .25 || r.density_at_transition == .5) &&
          (r.ruler_endpoint == .25 || r.ruler_endpoint == .5 || r.ruler_endpoint == 1) &&
          m.h0_km_s_mpc == 70 && m.omega_gamma == .5 && m.omega_massless_nonphoton == .125 &&
          m.omega_b == .0625 && m.omega_cdm == .0625 && m.species.empty(),
          "empirical reference has only its frozen dyadic domain");
  const W lambda = 1.L - m.omega_gamma - m.omega_massless_nonphoton - m.omega_b - m.omega_cdm;
  const auto today = density_record(r.density_at_transition, r.transition_scale, 1);
  const W step = W(r.ruler_endpoint) / n, loading = 3.L * m.omega_b / (4.L * m.omega_gamma),
      lambda_estimate = operation_floor * (1 + m.omega_gamma + m.omega_massless_nonphoton + m.omega_b + m.omega_cdm);
  require(step * n == W(r.ruler_endpoint) && loading == .09375L,
          "exact finite dyadic endpoint/mesh/loading source witnesses");
  W sum = 0, estimated_sum = 0;
  for (unsigned i = 0; i <= n; ++i) {
    const W a = step * i, a2 = a * a;
    const auto fluid = density_record(r.density_at_transition, r.transition_scale, a);
    const W p = m.omega_gamma + W(m.omega_massless_nonphoton) +
        (W(m.omega_b) + m.omega_cdm) * a +
        (lambda + fluid.value - today.value) * a2 * a2;
    const W f = 1 / std::sqrt(3 * p * (1 + loading * a)); // Original nominal expression/order.
    const W weight = i == 0 || i == n ? 1 : (i & 1 ? 4 : 2);
    sum += weight * f; // Original positive naive Simpson accumulation.
    const W a4 = a2 * a2, load = 1 + loading * a,
        pe = a4 * (lambda_estimate + fluid.estimate + today.estimate) +
          operation_floor * (m.omega_gamma + W(m.omega_massless_nonphoton) +
          std::abs((W(m.omega_b) + m.omega_cdm) * a) +
          (std::abs(lambda) + std::abs(fluid.value) + std::abs(today.value)) * a4),
        le = operation_floor * load;
    require(p > pe && load > le && std::isfinite(f), "positive reference dependency intervals");
    const W rp = pe / (std::sqrt(p - pe) * (std::sqrt(p) + std::sqrt(p - pe))),
        rl = le / (std::sqrt(load - le) * (std::sqrt(load) + std::sqrt(load - le)));
    estimated_sum += weight * f * (rp + rl + rp * rl + operation_floor);
  }
  const W value = W(irred::speed_of_light_m_per_s) / 1000 / m.h0_km_s_mpc * step * sum / 3,
      scale = W(irred::speed_of_light_m_per_s) / 1000 / m.h0_km_s_mpc * step / 3,
      product = W(n + 1) * std::numeric_limits<W>::epsilon(), gamma = product / (1 - product);
  require(product < 1 && gamma < 1 && estimated_sum > 0 && std::isfinite(value),
          "bounded positive reference sum/scale arithmetic");
  const W arithmetic = std::abs(scale) * (estimated_sum / (1 - gamma) + gamma * std::abs(sum)) +
      operation_floor * std::abs(value);
  const W cast = std::abs(value - W(static_cast<double>(value)));
  require(std::isnormal(arithmetic) && (cast == 0 || std::isnormal(cast)), "reference arithmetic/cast reporting resolves");
  return {value, arithmetic, cast};
}
bool adequate(W center, W reported, W reference_center, W reference_error, W allocation) {
  return std::isfinite(center) && std::isfinite(reported) && reported >= 0 &&
      std::isfinite(reference_error) && reference_error >= 0 && reference_error <= .05L * allocation &&
      std::abs(center - reference_center) <= allocation &&
      std::abs(center - reference_center) <= reported + reference_error;
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
    const auto coarse = reference(r, 4096), fine = reference(r, 8192);
    const W allocation = 1e-8L + 2e-10L * std::abs(fine.value), refinement = std::abs(coarse.value - fine.value),
        partial = refinement + coarse.arithmetic_estimate + fine.arithmetic_estimate +
          coarse.reporting_cast_loss + fine.reporting_cast_loss,
        reference_error = partial + operation_floor * (std::abs(coarse.value) + std::abs(fine.value) + partial);
    require(std::isnormal(reference_error) && reference_error <= .05L * allocation,
            "EACH complete empirical reference refinement/arithmetic/cast <=5% frozen allocation");
    require(adequate(*got.ruler->value, got.ruler->error_estimate, fine.value, reference_error, allocation),
            "native discrepancy covered by actual native diagnostic plus complete empirical reference error");
    require(!adequate(fine.value + allocation / 2, 0, fine.value, reference_error, allocation),
            "underreported zero diagnostic refuses original half-allocation witness");
    require(adequate(fine.value + reference_error / 2, 0, fine.value, reference_error, allocation),
            "within-reference-radius center is admitted");
    std::cout << std::hexfloat << "peer Xc=" << xc << " endpoint=" << endpoint
        << " refinement=" << refinement << " reference_arithmetic=" << coarse.arithmetic_estimate + fine.arithmetic_estimate
        << " reference_cast=" << coarse.reporting_cast_loss + fine.reporting_cast_loss
        << " Eref=" << reference_error << " native_estimate=" << got.ruler->error_estimate << std::defaultfloat << '\n';
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
