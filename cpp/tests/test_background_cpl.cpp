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
long double density_a(CPL p, long double a) {
  return p.omega_m / (a * a * a) +
         (1 - static_cast<long double>(p.omega_m)) *
             std::pow(a, -3 * (1 + static_cast<long double>(p.w0) + p.wa)) *
             std::exp(-3 * p.wa * (1 - a));
}
long double panels(CPL p, double z, bool clock, unsigned n) {
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
EvaluationPolicy policy() {
  EvaluationPolicy p;
  p.integration = irred::numerics::IntegrationPolicy{1e-12, 1e-12, 100000, 30};
  p.maximum_queries = 100;
  p.maximum_callbacks = 1000000;
  p.maximum_segment_visits = 1000;
  p.maximum_native_bytes = 1 << 20;
  return p;
}
Slot slot(Expansion b, double z, EvaluationPolicy p = policy(),
          double h0 = 70) {
  std::array<Request, 1> requests{
      {Request{z, 63, Observer{z, Convention::geometric_same_redshift},
               PhysicalScale{h0}}}};
  auto batch = b.evaluate(requests, p);
  check(batch.status == Status::ok && batch.slots.size() == 1, "batch");
  return batch.slots[0];
}
} // namespace
int main() {
  try {
    static_assert(!std::is_default_constructible_v<CPL>);
    for (double w : {-2., -1. - 1e-7, -1., -1. + 1e-7, -1. / 3 - 1e-7, -1. / 3,
                     -1. / 3 + 1e-7, 0.})
      for (double z : {0., 1e-4, .1, 1., 5.}) {
        auto s = slot(prepare(CPL{0, w, 0}, FlatFLRW{}), z);
        check(s.radial.status == Status::ok, "constantw finite");
        auto L = std::log1p(static_cast<long double>(z)),
             power = 1.5L * (1 + static_cast<long double>(w));
        auto I = power == 1 ? L : std::expm1((1 - power) * L) / (1 - power);
        auto J = power == 0 ? L : -std::expm1(-power * L) / power;
        near(s.radial.value->integral, I, 2e-15L + 2e-10L * std::abs(I),
             "constantw radial");
        // Clock's H0 factor recovered from independently related Hubble
        // distance/c.
        auto time = s.radial.value->integral == 0
                        ? 0.L
                        : s.physical.value->radial_mpc /
                              s.radial.value->integral *
                              irred::parsec_in_metres() * 1e6L /
                              irred::speed_of_light_m_per_s;
        if (z > 0)
          near(*s.clock.value->lookback_seconds.value / time, J,
               2e-15L + 2e-10L * std::abs(J), "constantw clock");
        near(s.expansion.value->expansion_E, std::exp(power * L),
             2e-13L + 2e-12L * std::exp(power * L), "constantw E");
        near(s.kinematics.value->q, power - 1,
             2e-13L + 2e-12L * std::abs(power - 1), "constantw q");
        near(s.kinematics.value->jerk.value(), (power - 1) * (2 * power - 1),
             2e-13L + 2e-12L * std::abs((power - 1) * (2 * power - 1)),
             "constantw jerk");
      }
    for (double om : {0., .3, 1.})
      for (double z : {0., 1e-4, .1, 1., 5.}) {
        auto c = slot(prepare(CPL{om, -1, 0}, FlatFLRW{}), z),
             l = slot(prepare(LCDM{om}, FlatFLRW{}), z);
        check(c.radial.status == Status::ok &&
                  l.radial.status == Status::ok,
              "Lambda limit finite");
        check(c.radial.value->integral == l.radial.value->integral &&
                  c.physical.value->luminosity_mpc ==
                      l.physical.value->luminosity_mpc &&
                  *c.clock.value->lookback_seconds.value ==
                      *l.clock.value->lookback_seconds.value,
              "Lambda distances bit parity");
        near(c.kinematics.value->q, l.kinematics.value->q, 2e-13L, "Lambda q");
        near(c.kinematics.value->jerk.value(), 1, 2e-13L, "Lambda jerk");
      }
    for (double om : {0., 1e-12, .3, 1 - 1e-12, 1.})
      for (double w : {-2., -.9, 0.})
        for (double wa : {-2., .7, 2.})
          for (double z : {.1, 1., 5.}) {
            CPL p{om, w, wa};
            auto s = slot(prepare(p, FlatFLRW{}), z);
            check(s.radial.status == Status::ok, "CPL grid finite");
            auto i1 = panels(p, z, false, 8192),
                 i2 = panels(p, z, false, 16384), j1 = panels(p, z, true, 8192),
                 j2 = panels(p, z, true, 16384);
            auto bi = 2e-15L + 2e-10L * std::abs(i2),
                 bj = 2e-15L + 2e-10L * std::abs(j2);
            near(i1, i2, bi / 10, "reference radial refinement");
            near(j1, j2, bj / 10, "reference clock refinement");
            near(s.radial.value->integral, i2, bi, "CPL a radial");
            auto time = s.physical.value->radial_mpc /
                        s.radial.value->integral * irred::parsec_in_metres() *
                        1e6L / irred::speed_of_light_m_per_s;
            near(*s.clock.value->lookback_seconds.value / time, j2, bj,
                 "CPL a clock");
            near(5 / std::log(10.L) * std::log(s.radial.value->integral / i2),
                 0, 1e-8L, "CPL magnitude");
            if (om == 1) {
              near(s.kinematics.value->q, .5L, 2e-13L, "allmatter q");
              near(s.kinematics.value->jerk.value(), 1, 2e-13L,
                   "allmatter jerk");
            }
          }
    for (auto p : {CPL{.3, -.9, .7}, CPL{0, -2, 2}, CPL{.8, 0, -2}})
      for (double z : {.25, 1., 3., 4.9}) {
        auto candidate = slot(prepare(p, FlatFLRW{}), z);
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
        near(candidate.kinematics.value->q, fine[0], 2e-6L,
             "CPL q independent derivatives");
        near(candidate.kinematics.value->jerk.value(), fine[1], 2e-6L,
             "CPL jerk independent derivatives");
        auto refined = policy();
        refined.integration =
            irred::numerics::IntegrationPolicy{1e-14, 1e-13, 100000, 30};
        auto second = slot(prepare(p, FlatFLRW{}), z, refined);
        refined.integration =
            irred::numerics::IntegrationPolicy{1e-15, 1e-14, 100000, 30};
        auto third = slot(prepare(p, FlatFLRW{}), z, refined);
        near(candidate.radial.value->integral, third.radial.value->integral,
             2e-15L + 2e-10L * third.radial.value->integral,
             "production integral refinement");
        near(5 / std::log(10.L) *
                 std::log(second.radial.value->integral /
                          third.radial.value->integral),
             0, 1e-8L, "production magnitude refinement");
      }
    auto zero = slot(prepare(CPL{.3, -.9, .7}, FlatFLRW{}), 0);
    check(zero.radial.status == Status::ok &&
              zero.expansion.value->expansion_E == 1 &&
              zero.expansion.value->h_km_s_mpc.value.value() == 70,
          "z0 expansion valid");
    check(zero.physical.value->radial_mpc == 0 &&
              zero.physical.value->luminosity_mpc == 0 &&
              *zero.clock.value->lookback_seconds.value == 0 &&
              zero.physical.value->volume_mpc3_per_sr_per_redshift == 0,
          "z0 analytic zeros");
    auto background = prepare(CPL{.3, -.9, .7}, FlatFLRW{});
    std::array<Request, 2> observer_queries{
        {Request{.5, 63, Observer{.5, Convention::geometric_same_redshift},
                 PhysicalScale{70}},
         Request{.5, 63, Observer{.7, Convention::released_zhd_zhel},
                 PhysicalScale{70}}}};
    auto observers = background.evaluate(observer_queries, policy());
    check(observers.slots[0].physical.value &&
              observers.slots[1].physical.value,
          "observer geometry finite");
    near(observers.slots[1].physical.value->luminosity_mpc /
             observers.slots[0].physical.value->luminosity_mpc,
         1.7L / 1.5L, 2e-13L, "observer prefactor");
    check(observers.slots[0].physical.value->angular_diameter_mpc ==
              observers.slots[1].physical.value->angular_diameter_mpc,
          "DA observer invariant");
    check(observers.slots[0].source.observer->convention !=
              observers.slots[1].source.observer->convention,
          "explicit equation convention identity");
    auto capped = policy();
    capped.maximum_queries = 1;
    auto cap_failure = background.evaluate(observer_queries, capped);
    check(cap_failure.status == Status::work_limit && cap_failure.slots.empty(),
          "query cap before output");
    observer_queries[0].z_expansion = -.1;
    observer_queries[1].z_expansion = 5.0001;
    auto domains = background.evaluate(observer_queries, policy());
    check(domains.slots[0].radial.status == Status::unsupported_domain &&
              domains.slots[1].radial.status == Status::unsupported_domain,
          "redshift domain");
    check(prepare(CPL{-1e-15, -1, 0}, FlatFLRW{}).status() ==
              Status::unsupported_domain,
          "matter lower bound");
    check(prepare(CPL{1 + 1e-15, -1, 0}, FlatFLRW{}).status() ==
              Status::unsupported_domain,
          "matter upper bound");
    auto valid = prepare(CPL{.3, -1, 0}, FlatFLRW{});
    auto bad_scale =
        slot(valid, 1, policy(), std::numeric_limits<double>::infinity());
    check(bad_scale.radial.value &&
              bad_scale.physical.status == Status::invalid_input,
          "invalid H0 isolated from expansion");
    auto tiny_scale =
        slot(valid, 1, policy(), std::numeric_limits<double>::min());
    check(tiny_scale.radial.value &&
              tiny_scale.physical.status == Status::numerical_failure,
          "tiny H0 representability failure belongs to physical projection");
    auto a = slot(prepare(CPL{.3, -.9, .7}, FlatFLRW{}), 1, policy(), 40),
         b = slot(prepare(CPL{.3, -.9, .7}, FlatFLRW{}), 1, policy(), 100);
    near(a.physical.value->radial_mpc / b.physical.value->radial_mpc, 2.5L,
         2e-12L, "H0 distance");
    near(*a.clock.value->lookback_seconds.value /
             *b.clock.value->lookback_seconds.value,
         2.5L, 2e-12L, "H0 clock");
    near(a.physical.value->volume_mpc3_per_sr_per_redshift /
             b.physical.value->volume_mpc3_per_sr_per_redshift,
         15.625L, 2e-11L, "H0 volume");
    auto attempted = prepare(CPL{.3, 1, 0}, FlatFLRW{});
    check(attempted.status() == Status::unsupported_domain &&
              attempted.model_id() == "P01/flat-cpl-radiation-free/v1" &&
              std::get<CPL>(attempted.specification()).w0 == 1,
          "failed attempted identity");
    // Removed old model-number admission/inactive fields: active variants
    // retain distinct hypotheses without compatibility projection packets.
    static_assert(!std::is_default_constructible_v<CPL>);
    check(prepare(CPL{.3, std::numeric_limits<double>::quiet_NaN(), 0},
                  FlatFLRW{})
                  .status() == Status::invalid_input,
          "NaN w");
    check(prepare(CPL{.3, -1, 2.0001}, FlatFLRW{}).status() ==
              Status::unsupported_domain,
          "wa range");
    auto good = prepare(CPL{.3, -.9, .7}, FlatFLRW{});
    auto tiny = policy();
    tiny.maximum_callbacks = 3;
    auto failed = slot(good, 1, tiny);
    check(failed.radial.status == Status::work_limit &&
              failed.radial.numerical_status ==
                  irred::numerics::Status::work_limit,
          "global work limit");
    check(slot(good, 1e-200).physical.status == Status::numerical_failure,
          "physical underflow");
    std::printf("CPL owner: %d checks PASS; max budget fraction %.9Lg\n",
                checks, maximum_fraction);
  } catch (const std::exception &e) {
    std::fprintf(stderr, "%s\n", e.what());
    return 1;
  }
}
