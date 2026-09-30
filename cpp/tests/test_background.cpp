// Owned-background comparison suite. Analytic EdS/deSitter/constant-q
// definitions and fixed-panel Simpson quadrature differ from production
// adaptive Simpson. These test named late-time geometries, not empirical FLRW
// truth. Dimensionless integral budget 2e-15+2e-10*|I|; propagated mu <=1e-8
// mag for z>=1e-4. Reference refinement uncertainty is checked separately.
#include "irred/background.hpp"
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
using namespace irred::cosmology;
namespace {
int checks = 0;
long double maximum_integral_error = 0, maximum_mu_error = 0;
void check(bool value, const char *why) {
  ++checks;
  if (!value)
    throw std::runtime_error(why);
}
void near(long double value, long double expected, long double absolute,
          long double relative, const char *why) {
  check(std::abs(value - expected) <= absolute + relative * std::abs(expected),
        why);
}
std::vector<Request> queries(double h0 = 70) {
  std::vector<Request> out;
  const double z[] = {0, 1e-10, 1e-4, .01, .1, .5, 1, 2, 5};
  for (double value : z)
    out.emplace_back(value, 63,
                     Observer{value, Convention::geometric_same_redshift},
                     PhysicalScale{h0});
  return out;
}
long double q_integral(long double z, long double q) {
  return q == 0 ? std::log1p(z) : -std::expm1(-q * std::log1p(z)) / q;
}
template <class F> long double fixed(F f, long double z, unsigned panels) {
  const auto h = z / panels;
  long double sum = f(0) + f(z);
  for (unsigned i = 1; i < panels; ++i)
    sum += (i % 2 ? 4 : 2) * f(h * i);
  return h * sum / 3;
}
void check_slot(const Slot &slot, long double integral, long double clock,
                double h0, long double expected_q, long double expected_jerk) {
  check(slot.radial.value && slot.clock.value && slot.physical.value &&
            slot.kinematics.value,
        "requested groups available");
  const auto dh = 299792.458L / h0; // exact SI c, explicit km/s/Mpc scaling
  // Pinned independent Machin/Decimal SI-IAU reference, generated before core.
  const auto ht =
      (70.L / h0) / 2.2685455026110555162663814014470700829931485893862e-18L;
  maximum_integral_error =
      std::max(maximum_integral_error,
               std::abs(static_cast<long double>(slot.radial.value->integral) -
                        integral));
  near(slot.radial.value->integral, integral, 2e-15L, 2e-10L,
       "radial integral budget");
  near(slot.physical.value->radial_mpc, dh * integral, dh * 2e-15L, 2e-10L,
       "distance unit and budget");
  near(*slot.clock.value->lookback_seconds.value, ht * clock, ht * 2e-15L, 2e-10L,
       "clock units and budget");
  near(slot.kinematics.value->q, expected_q, 2e-14L, 2e-14L, "q definition");
  near(*slot.kinematics.value->jerk, expected_jerk, 2e-14L, 2e-14L,
       "jerk definition");
  near(slot.physical.value->luminosity_mpc,
       slot.physical.value->angular_diameter_mpc *
           std::pow(1 + slot.source.z_expansion, 2),
       dh * 2e-15L, 2e-10L, "same-redshift reciprocity");
  if (slot.source.z_expansion >= 1e-4) {
    const auto target = dh * integral * (1 + slot.source.observer->redshift);
    const auto error = std::abs(
        5 / std::log(10.L) *
        std::log(static_cast<long double>(slot.physical.value->luminosity_mpc) /
                 target));
    maximum_mu_error = std::max(maximum_mu_error, error);
    check(error <= 1e-8, "magnitude-consumer error budget");
  }
  if (slot.source.z_expansion == 0)
    check(slot.physical.value->radial_mpc == 0 &&
              slot.physical.value->luminosity_mpc == 0 &&
              slot.physical.value->angular_diameter_mpc == 0 &&
              *slot.clock.value->lookback_seconds.value == 0 &&
              slot.physical.value->volume_mpc3_per_sr_per_redshift == 0,
          "z0 valid exact zeros, no log or quadrature");
}
} // namespace
int main() {
  EvaluationPolicy policy;
  policy.integration =
      irred::numerics::IntegrationPolicy{1e-12, 1e-12, 100000, 30};
  policy.maximum_queries = 100;
  policy.maximum_callbacks = 2000000;
  policy.maximum_segment_visits = 1000;
  policy.maximum_native_bytes = 1 << 20;
  auto grid = queries();
  for (double h0 : {40., 70., 100.}) {
    grid = queries(h0);
    for (double q :
         {-2., -1.0000000001, -1., -.9999999999, -1e-10, 0., 1e-10, .5, 2.}) {
      const auto background = prepare(ConstantQ{q}, FlatFLRW{});
      check(background.status() == Status::ok, "constant-q preparation");
      const auto batch = background.evaluate(grid, policy);
      check(batch.status == Status::ok && batch.slots.size() == grid.size(),
            "coarse stable order");
      for (std::size_t i = 0; i < grid.size(); ++i) {
        const auto z = grid[i].z_expansion;
        check_slot(batch.slots[i], q_integral(z, q), q_integral(z, 1 + q), h0,
                   q, q * (2 * q + 1));
        near(batch.slots[i].expansion.value->expansion_E,
             std::pow(1 + static_cast<long double>(z), 1 + q), 2e-14, 2e-14,
             "constant-q E");
      }
    }
  }
  grid = queries(70);
  for (double omega : {0., .3, 1.}) {
    const auto background = prepare(LCDM{omega}, FlatFLRW{});
    const auto batch = background.evaluate(grid, policy);
    check(batch.status == Status::ok, "LCDM batch");
    for (std::size_t i = 0; i < grid.size(); ++i) {
      const auto z = grid[i].z_expansion;
      const auto inverse_E = [omega](long double x) {
        return 1 / std::sqrt(omega * std::pow(1 + x, 3) + 1 - omega);
      };
      const auto reference = fixed(inverse_E, z, 8192);
      const auto coarse = fixed(inverse_E, z, 4096);
      check(std::abs(coarse - reference) <=
                (2e-15L + 2e-10L * std::abs(reference)) / 10,
            "independent fixed-panel reference refinement");
      const auto clock =
          fixed([&](long double x) { return inverse_E(x) / (1 + x); }, z, 8192);
      const auto e2 =
          omega * std::pow(1 + static_cast<long double>(z), 3) + 1 - omega;
      check_slot(
          batch.slots[i], reference, clock, 70,
          1.5L * omega * std::pow(1 + static_cast<long double>(z), 3) / e2 - 1,
          1);
      if (omega == 0)
        near(reference, z, 2e-15, 2e-10, "deSitter distance");
      if (omega == 1)
        near(reference,
             2 * (1 - 1 / std::sqrt(1 + static_cast<long double>(z))), 2e-15,
             2e-10, "EdS distance");
    }
  }
  const auto background = prepare(LCDM{.3}, FlatFLRW{});
  std::array<Request, 3> mixed{
      {Request{.2, 63, Observer{.3, Convention::released_zhd_zhel},
               PhysicalScale{70}},
       Request{.2, 63, Observer{.3, Convention::geometric_same_redshift},
               PhysicalScale{70}},
       Request{0, 63, Observer{0, Convention::geometric_same_redshift},
               PhysicalScale{70}}}};
  const auto batch = background.evaluate(mixed, policy);
  check(batch.slots[0].physical.value &&
            batch.slots[1].physical.status == Status::incompatible_convention &&
            batch.slots[1].radial.value && batch.slots[2].physical.value,
        "observer failure isolated from radial");
  near(batch.slots[0].physical.value->luminosity_mpc /
           batch.slots[0].physical.value->transverse_mpc,
       1.3L, 2e-14, 2e-14, "released observer prefactor");
  near(batch.slots[0].physical.value->angular_diameter_mpc /
           batch.slots[0].physical.value->transverse_mpc,
       1 / 1.2L, 2e-14, 2e-14, "DA expansion redshift");
  // Old copied per-slot equation string checks replaced by explicit typed
  // observer convention provenance, with the same distance formula controls.
  check(batch.slots[0].source.observer->convention !=
            batch.slots[2].source.observer->convention,
        "convention retained in source");
  auto cap = policy;
  cap.maximum_queries = 2;
  check(background.evaluate(grid, cap).status == Status::work_limit &&
            background.evaluate(grid, cap).slots.empty(),
        "batch limit before slot allocation");
  cap = policy;
  cap.maximum_callbacks = 3;
  const auto limited = background.evaluate(grid, cap);
  check(limited.work.callbacks <= 3 && limited.slots[0].radial.value &&
            limited.slots.back().radial.status == Status::work_limit,
        "global radial-clock budget");
  // Inactive fields no longer exist in the active model variant. Domain safety
  // remains tested on active parameters, while H0 is a projection-only input.
  check(prepare(LCDM{-.1}, FlatFLRW{}).status() == Status::unsupported_domain,
        "active LCDM domain");
  std::array<Request, 1> invalid_h{
      {Request{.5, 63, Observer{.5, Convention::geometric_same_redshift},
               PhysicalScale{-70}}}};
  auto h = background.evaluate(invalid_h, policy);
  check(h.slots[0].radial.value &&
            h.slots[0].physical.status == Status::invalid_input,
        "invalid physical scale does not poison dimensionless integral");
  std::array<Request, 3> extremes{
      {Request{1e-300, 63,
               Observer{1e-300, Convention::geometric_same_redshift},
               PhysicalScale{70}},
       Request{1e-100, 63,
               Observer{1e-100, Convention::geometric_same_redshift},
               PhysicalScale{70}},
       Request{std::numeric_limits<double>::denorm_min(), 63,
               Observer{std::numeric_limits<double>::denorm_min(),
                        Convention::geometric_same_redshift},
               PhysicalScale{70}}}};
  const auto extreme = background.evaluate(extremes, policy);
  check(extreme.slots[0].physical.status == Status::numerical_failure &&
            extreme.slots[0].radial.value,
        "physical volume underflow isolated");
  check(extreme.slots[1].physical.value.has_value(),
        "tiny normal physical output accepted");
  check(extreme.slots[2].radial.status == Status::work_limit &&
            extreme.slots[2].radial.numerical_status ==
                irred::numerics::Status::work_limit &&
            extreme.nodes[*extreme.slots[2].node_index].work.callbacks == 6,
        "smallest interval truthful failed radial+clock work");
  std::printf("{\"suite\":\"background_owner\",\"checks\":%d,\"max_integral_"
              "error\":%.17Lg,\"max_mu_error\":%.17Lg,\"passed\":true}\n",
              checks, maximum_integral_error, maximum_mu_error);
}
