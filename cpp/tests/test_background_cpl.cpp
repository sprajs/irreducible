// CPL background only: flat matter + conserved CPL homogeneous density.
// Independent reference variable a=1/(1+z), composite Simpson panel refinement;
// analytic pure constant-w and EdS limits. Shared libm ancestry disclosed.
// Budgets fixed before candidate: I/J abs2e-15+rel2e-10; refinement<=1/10;
// analytic E/q/j abs2e-13+rel2e-12; magnitude<=1e-8 for z>=1e-4.
// No external assets, likelihood/perturbation/early-universe qualification.
#include "irred/background.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <type_traits>
using namespace irred::cosmology;
namespace {
int checks = 0;
long double maximum_fraction = 0;
void check(bool v, const char *s) {
  ++checks;
  if (!v)
    throw std::runtime_error(s);
}
void near(long double x, long double y, long double budget, const char *s) {
  check(std::abs(x - y) <= budget, s);
  maximum_fraction = std::max(maximum_fraction, std::abs(x - y) / budget);
}
long double density_a(CplParameters p, long double a) {
  return p.omega_m / (a * a * a) +
         (1 - static_cast<long double>(p.omega_m)) *
             std::pow(a, -3 * (1 + static_cast<long double>(p.w0) + p.wa)) *
             std::exp(-3 * p.wa * (1 - a));
}
long double panels(CplParameters p, double z, bool clock, unsigned n) {
  auto lower = 1 / (1 + static_cast<long double>(z));
  auto h = (1 - lower) / n;
  auto f = [&](long double a) {
    return 1 / (std::sqrt(density_a(p, a)) * (clock ? a : a * a));
  };
  long double sum = f(lower) + f(1);
  for (unsigned i = 1; i < n; ++i)
    sum += (i % 2 ? 4 : 2) * f(lower + i * h);
  return h * sum / 3;
}
Slot slot(Background b, double z, Policy p = {}) {
  std::array<Query, 1> queries{{{z, z, Convention::geometric_same_redshift}}};
  auto batch = b.evaluate_batch(queries, p);
  check(batch.status == Status::ok && batch.slots.size() == 1, "batch");
  return batch.slots[0];
}
} // namespace
int main() {
  try {
    static_assert(!std::is_default_constructible_v<CplParameters>);
    for (double w : {-2., -1. - 1e-7, -1., -1. + 1e-7, -1. / 3 - 1e-7, -1. / 3,
                     -1. / 3 + 1e-7, 0.})
      for (double z : {0., 1e-4, .1, 1., 5.}) {
        auto s = slot(prepare_cpl({70, 0, w, 0}), z);
        check(s.status == Status::ok, "constantw finite");
        auto L = std::log1p(static_cast<long double>(z)),
             power = 1.5L * (1 + static_cast<long double>(w));
        auto I = power == 1 ? L : std::expm1((1 - power) * L) / (1 - power);
        auto J = power == 0 ? L : -std::expm1(-power * L) / power;
        near(s.radial_integral, I, 2e-15L + 2e-10L * std::abs(I),
             "constantw radial");
        // Clock's H0 factor recovered from independently related Hubble
        // distance/c.
        auto time = s.radial_integral == 0
                        ? 0.L
                        : s.radial_mpc / s.radial_integral *
                              irred::parsec_in_metres() * 1e6L /
                              irred::speed_of_light_m_per_s;
        if (z > 0)
          near(s.lookback_seconds / time, J, 2e-15L + 2e-10L * std::abs(J),
               "constantw clock");
        near(s.expansion_E, std::exp(power * L),
             2e-13L + 2e-12L * std::exp(power * L), "constantw E");
        near(s.deceleration_q, power - 1, 2e-13L + 2e-12L * std::abs(power - 1),
             "constantw q");
        near(s.jerk, (power - 1) * (2 * power - 1),
             2e-13L + 2e-12L * std::abs((power - 1) * (2 * power - 1)),
             "constantw jerk");
      }
    for (double om : {0., .3, 1.})
      for (double z : {0., 1e-4, .1, 1., 5.}) {
        auto c = slot(prepare_cpl({70, om, -1, 0}), z),
             l = slot(prepare({Model::flat_lcdm_late_v1, 70, om, 0}), z);
        check(c.status == Status::ok && l.status == Status::ok,
              "Lambda limit finite");
        check(c.radial_integral == l.radial_integral &&
                  c.luminosity_mpc == l.luminosity_mpc &&
                  c.lookback_seconds == l.lookback_seconds,
              "Lambda distances bit parity");
        near(c.deceleration_q, l.deceleration_q, 2e-13L, "Lambda q");
        near(c.jerk, 1, 2e-13L, "Lambda jerk");
      }
    for (double om : {0., 1e-12, .3, 1 - 1e-12, 1.})
      for (double w : {-2., -.9, 0.})
        for (double wa : {-2., .7, 2.})
          for (double z : {.1, 1., 5.}) {
            CplParameters p{70, om, w, wa};
            auto s = slot(prepare_cpl(p), z);
            check(s.status == Status::ok, "CPL grid finite");
            auto i1 = panels(p, z, false, 8192),
                 i2 = panels(p, z, false, 16384), j1 = panels(p, z, true, 8192),
                 j2 = panels(p, z, true, 16384);
            auto bi = 2e-15L + 2e-10L * std::abs(i2),
                 bj = 2e-15L + 2e-10L * std::abs(j2);
            near(i1, i2, bi / 10, "reference radial refinement");
            near(j1, j2, bj / 10, "reference clock refinement");
            near(s.radial_integral, i2, bi, "CPL a radial");
            auto time = s.radial_mpc / s.radial_integral *
                        irred::parsec_in_metres() * 1e6L /
                        irred::speed_of_light_m_per_s;
            near(s.lookback_seconds / time, j2, bj, "CPL a clock");
            near(5 / std::log(10.L) * std::log(s.radial_integral / i2), 0,
                 1e-8L, "CPL magnitude");
            if (om == 1) {
              near(s.deceleration_q, .5L, 2e-13L, "allmatter q");
              near(s.jerk, 1, 2e-13L, "allmatter jerk");
            }
          }
    for (auto p : {CplParameters{70, .3, -.9, .7}, CplParameters{70, 0, -2, 2},
                   CplParameters{70, .8, 0, -2}})
      for (double z : {.25, 1., 3., 4.9}) {
        auto candidate = slot(prepare_cpl(p), z);
        auto derivatives = [&](long double h) {
          auto S = [&](long double zz) { return density_a(p, 1 / (1 + zz)); };
          const long double zz = z, sm2 = S(zz - 2 * h), sm = S(zz - h),
                            s0 = S(zz), sp = S(zz + h), sp2 = S(zz + 2 * h),
                            u = 1 + zz;
          auto d1 = (sm2 - 8 * sm + 8 * sp - sp2) / (12 * h);
          auto d2 = (-sp2 + 16 * sp - 30 * s0 + 16 * sm - sm2) / (12 * h * h);
          return std::array<long double, 2>{
              -1 + .5L * u * d1 / s0, 1 - u * d1 / s0 + .5L * u * u * d2 / s0};
        };
        auto coarse = derivatives(1e-3L), medium = derivatives(5e-4L),
             fine = derivatives(2.5e-4L);
        near(coarse[0], medium[0], 2e-6L, "derivative coarse q convergence");
        near(medium[0], fine[0], 2e-6L, "derivative fine q convergence");
        near(medium[1], fine[1], 2e-6L, "derivative fine jerk convergence");
        near(candidate.deceleration_q, fine[0], 2e-6L,
             "CPL q independent derivatives");
        near(candidate.jerk, fine[1], 2e-6L,
             "CPL jerk independent derivatives");
        Policy refined{};
        refined.integration = {1e-14, 1e-13, 100000, 30};
        auto second = slot(prepare_cpl(p), z, refined);
        refined.integration = {1e-15, 1e-14, 100000, 30};
        auto third = slot(prepare_cpl(p), z, refined);
        near(candidate.radial_integral, third.radial_integral,
             2e-15L + 2e-10L * third.radial_integral,
             "production integral refinement");
        near(5 / std::log(10.L) *
                 std::log(second.radial_integral / third.radial_integral),
             0, 1e-8L, "production magnitude refinement");
      }
    auto zero = slot(prepare_cpl({70, .3, -.9, .7}), 0);
    check(zero.status == Status::ok && zero.expansion_E == 1 &&
              zero.h_km_s_mpc == 70,
          "z0 expansion valid");
    check(zero.radial_mpc == 0 && zero.luminosity_mpc == 0 &&
              zero.lookback_seconds == 0 &&
              zero.volume_mpc3_per_sr_per_redshift == 0,
          "z0 analytic zeros");
    auto background = prepare_cpl({70, .3, -.9, .7});
    std::array<Query, 2> observer_queries{
        {{1, 1, Convention::geometric_same_redshift},
         {1, .3, Convention::released_zhd_zhel}}};
    auto observers = background.evaluate_batch(observer_queries, {});
    check(observers.slots[0].status == Status::ok &&
              observers.slots[1].status == Status::ok,
          "observer modes finite");
    near(observers.slots[1].luminosity_mpc / observers.slots[0].luminosity_mpc,
         .65L, 2e-13L, "released observer prefactor");
    check(observers.slots[0].angular_diameter_mpc ==
              observers.slots[1].angular_diameter_mpc,
          "angular expansion redshift retained");
    check(observers.slots[0].luminosity_equation_id !=
              observers.slots[1].luminosity_equation_id,
          "convention identity distinct");
    Policy capped{};
    capped.maximum_queries = 1;
    auto cap_failure = background.evaluate_batch(observer_queries, capped);
    check(cap_failure.status == Status::work_limit && cap_failure.slots.empty(),
          "batch cap before output allocation");
    observer_queries[0] = {-1, -1, Convention::geometric_same_redshift};
    observer_queries[1] = {5.0001, 5.0001, Convention::geometric_same_redshift};
    auto domains = background.evaluate_batch(observer_queries, {});
    check(domains.slots[0].status == Status::unsupported_domain &&
              domains.slots[1].status == Status::unsupported_domain,
          "z domain boundaries");
    check(prepare_cpl({70, -1e-15, -1, 0}).status() ==
              Status::unsupported_domain,
          "matter lower bound");
    check(prepare_cpl({70, 1 + 1e-15, -1, 0}).status() ==
              Status::unsupported_domain,
          "matter upper bound");
    check(prepare_cpl({std::numeric_limits<double>::infinity(), .3, -1, 0})
                  .status() == Status::invalid_input,
          "H0 infinity");
    check(
        prepare_cpl({std::numeric_limits<double>::min(), .3, -1, 0}).status() ==
            Status::numerical_failure,
        "H0 extreme representation");
    auto a = slot(prepare_cpl({40, .3, -.9, .7}), 1),
         b = slot(prepare_cpl({100, .3, -.9, .7}), 1);
    near(a.radial_mpc / b.radial_mpc, 2.5L, 2e-12L, "H0 distance");
    near(a.lookback_seconds / b.lookback_seconds, 2.5L, 2e-12L, "H0 clock");
    near(a.volume_mpc3_per_sr_per_redshift / b.volume_mpc3_per_sr_per_redshift,
         15.625L, 2e-11L, "H0 volume");
    auto attempted = prepare_cpl({70, .3, 1, 0});
    check(attempted.status() == Status::unsupported_domain &&
              attempted.model_id() == "P01/flat-cpl-radiation-free/v1" &&
              attempted.parameters().w0 == 1,
          "failed attempted identity");
    check(prepare({Model::flat_cpl_late_v1, 70, .3, 0, -1, 0}).status() ==
              Status::invalid_input,
          "legacy cannot admit CPL");
    check(prepare({Model::flat_lcdm_late_v1, 70, .3, 0, -.9, 0}).status() ==
              Status::invalid_input,
          "legacy inactive w rejected");
    check(prepare_cpl({70, .3, std::numeric_limits<double>::quiet_NaN(), 0})
                  .status() == Status::invalid_input,
          "NaN w");
    check(prepare_cpl({70, .3, -1, 2.0001}).status() ==
              Status::unsupported_domain,
          "wa range");
    auto good = prepare_cpl({70, .3, -.9, .7});
    Policy tiny{};
    tiny.maximum_total_evaluations = 3;
    auto failed = slot(good, 1, tiny);
    check(failed.status == Status::work_limit && failed.evaluations <= 3,
          "global work limit");
    check(slot(good, 1e-200).status == Status::numerical_failure,
          "physical underflow");
    std::printf("CPL owner: %d checks PASS; max budget fraction %.9Lg\n",
                checks, maximum_fraction);
  } catch (const std::exception &e) {
    std::fprintf(stderr, "%s\n", e.what());
    return 1;
  }
}
