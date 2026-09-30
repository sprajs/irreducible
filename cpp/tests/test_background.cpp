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
std::array<Query, 9> queries() {
  std::array<Query, 9> out{};
  const double z[] = {0, 1e-10, 1e-4, .01, .1, .5, 1, 2, 5};
  for (std::size_t i = 0; i < out.size(); ++i)
    out[i] = {z[i], z[i], Convention::geometric_same_redshift};
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
  check(slot.status == Status::ok, "valid query bundle");
  const auto dh = 299792.458L / h0; // exact SI c, explicit km/s/Mpc scaling
  // Pinned independent Machin/Decimal SI-IAU reference, generated before core.
  const auto ht =
      (70.L / h0) / 2.2685455026110555162663814014470700829931485893862e-18L;
  maximum_integral_error = std::max(
      maximum_integral_error,
      std::abs(static_cast<long double>(slot.radial_integral) - integral));
  near(slot.radial_integral, integral, 2e-15L, 2e-10L,
       "radial integral budget");
  near(slot.radial_mpc, dh * integral, dh * 2e-15L, 2e-10L,
       "distance unit and budget");
  near(slot.lookback_seconds, ht * clock, ht * 2e-15L, 2e-10L,
       "clock units and budget");
  near(slot.deceleration_q, expected_q, 2e-14L, 2e-14L, "q definition");
  near(slot.jerk, expected_jerk, 2e-14L, 2e-14L, "jerk definition");
  near(slot.luminosity_mpc,
       slot.angular_diameter_mpc * std::pow(1 + slot.source.z_expansion, 2),
       dh * 2e-15L, 2e-10L, "same-redshift reciprocity");
  if (slot.source.z_expansion >= 1e-4) {
    const auto target = dh * integral * (1 + slot.source.z_observer);
    const auto error = std::abs(
        5 / std::log(10.L) *
        std::log(static_cast<long double>(slot.luminosity_mpc) / target));
    maximum_mu_error = std::max(maximum_mu_error, error);
    check(error <= 1e-8, "magnitude-consumer error budget");
  }
  if (slot.source.z_expansion == 0)
    check(slot.radial_mpc == 0 && slot.luminosity_mpc == 0 &&
              slot.angular_diameter_mpc == 0 && slot.lookback_seconds == 0 &&
              slot.volume_mpc3_per_sr_per_redshift == 0 &&
              slot.evaluations == 0,
          "z0 valid exact zeros, no log or quadrature");
}
} // namespace
int main() {
  const auto grid = queries();
  for (double h0 : {40., 70., 100.}) {
    for (double q :
         {-2., -1.0000000001, -1., -.9999999999, -1e-10, 0., 1e-10, .5, 2.}) {
      const auto background = prepare({Model::constant_q_flat_v1, h0, 0, q});
      check(background.status() == Status::ok, "constant-q preparation");
      const auto batch = background.evaluate_batch(grid, {});
      check(batch.status == Status::ok && batch.slots.size() == grid.size(),
            "coarse stable order");
      for (std::size_t i = 0; i < grid.size(); ++i) {
        const auto z = grid[i].z_expansion;
        check_slot(batch.slots[i], q_integral(z, q), q_integral(z, 1 + q), h0,
                   q, q * (2 * q + 1));
        near(batch.slots[i].expansion_E,
             std::pow(1 + static_cast<long double>(z), 1 + q), 2e-14, 2e-14,
             "constant-q E");
      }
    }
  }
  for (double omega : {0., .3, 1.}) {
    const auto background = prepare({Model::flat_lcdm_late_v1, 70, omega, 0});
    const auto batch = background.evaluate_batch(grid, {});
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
  const auto background = prepare({Model::flat_lcdm_late_v1, 70, .3, 0});
  std::array<Query, 3> mixed{{{.2, .3, Convention::released_zhd_zhel},
                              {.2, .3, Convention::geometric_same_redshift},
                              {0, 0, Convention::geometric_same_redshift}}};
  const auto batch = background.evaluate_batch(mixed, {});
  check(batch.slots[0].status == Status::ok &&
            batch.slots[1].status == Status::incompatible_convention &&
            batch.slots[2].status == Status::ok,
        "mixed conventions stable failure slots");
  near(batch.slots[0].luminosity_mpc / batch.slots[0].transverse_mpc, 1.3L,
       2e-14, 2e-14, "released observer prefactor");
  near(batch.slots[0].angular_diameter_mpc / batch.slots[0].transverse_mpc,
       1 / 1.2L, 2e-14, 2e-14, "DA uses expansion redshift");
  check(batch.slots[0].luminosity_equation_id !=
                batch.slots[2].luminosity_equation_id &&
            batch.slots[0].shape_equation_id !=
                batch.slots[2].shape_equation_id,
        "convention-specific equations");
  Policy cap;
  cap.maximum_queries = 2;
  check(background.evaluate_batch(grid, cap).status == Status::work_limit &&
            background.evaluate_batch(grid, cap).slots.empty(),
        "batch limit before slot allocation");
  cap = {};
  cap.maximum_total_evaluations = 3;
  const auto limited = background.evaluate_batch(grid, cap);
  std::size_t used = 0;
  for (const auto &slot : limited.slots)
    used += slot.evaluations;
  check(used <= 3 && limited.slots[0].status == Status::ok &&
            limited.slots.back().status == Status::work_limit,
        "shared radial-clock work budget");
  check(prepare({Model::flat_lcdm_late_v1, 70, .3, 1}).status() ==
            Status::invalid_input,
        "inactive model field rejects");
  check(prepare({Model::flat_lcdm_late_v1, -70, .3, 0}).status() ==
            Status::invalid_input,
        "positive expansion H0");
  std::array<Query, 3> extremes{
      {{1e-300, 1e-300, Convention::geometric_same_redshift},
       {1e-100, 1e-100, Convention::geometric_same_redshift},
       {std::numeric_limits<double>::denorm_min(),
        std::numeric_limits<double>::denorm_min(),
        Convention::geometric_same_redshift}}};
  const auto extreme = background.evaluate_batch(extremes, {});
  check(extreme.slots[0].status == Status::numerical_failure,
        "nonzero physical volume underflow fails bundle");
  check(extreme.slots[1].status == Status::ok,
        "tiny but normal physical bundle accepted");
  check(extreme.slots[2].status == Status::work_limit &&
            extreme.slots[2].numerical_status ==
                irred::numerics::Status::work_limit &&
            extreme.slots[2].evaluations == 3,
        "smallest positive interval reports actual quadrature work failure");
  std::printf("{\"suite\":\"background_owner\",\"checks\":%d,\"max_integral_"
              "error\":%.17Lg,\"max_mu_error\":%.17Lg,\"passed\":true}\n",
              checks, maximum_integral_error, maximum_mu_error);
}
