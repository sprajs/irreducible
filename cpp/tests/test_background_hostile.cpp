// Independent P01 tests: closed EdS/deSitter integrals; long-double fixed-panel
// Simpson reference written here, NOT production F02 quadrature. Three panel
// counts4096/8192/16384 challenge reference refinement; logs/exp share libm.
// Dimensions/constants independently SI: Mpc=1e6*(648000/pi)*149597870700 m.
// Budgets predeclared: integral abs2e-15+rel2e-10; reference refinement<=1/10;
// same scaled distance budget, clocks same relative dimensionless budget;
// delta magnitude<=1e-8 at z>=1e-4. Synthetic cases, not W01
// target/reproduction.
#include "irred/background.hpp"
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <string_view>
using namespace irred::cosmology;
namespace {
int checks = 0;
double maximum_integral = 0, maximum_mag = 0, maximum_reference_delta = 0;
long double radial_fraction = 0, clock_fraction = 0, distance_fraction = 0,
            reference_fraction = 0;
void check(bool b, const char *n) {
  ++checks;
  if (!b)
    throw std::runtime_error(n);
}
void near(long double v, long double ref, long double budget, const char *n) {
  const auto fraction = std::abs(v - ref) / budget;
  const std::string_view name(n);
  if (name.find("radial integral") != name.npos ||
      name.find("defining integral") != name.npos ||
      name.find("independent integration") != name.npos)
    radial_fraction = std::max(radial_fraction, fraction);
  if (name.find("clock") != name.npos)
    clock_fraction = std::max(clock_fraction, fraction);
  if (name.find("physical distance") != name.npos)
    distance_fraction = std::max(distance_fraction, fraction);
  check(std::abs(v - ref) <= budget, n);
}
long double scale(long double value) {
  return 2e-15L + 2e-10L * std::abs(value);
}
long double reference(Parameters p, long double z, bool clock,
                      unsigned panels) {
  auto f = [&](long double x) {
    const auto u = 1 + x;
    long double e;
    if (p.model == Model::flat_lcdm_late_v1)
      e = std::sqrt(static_cast<long double>(p.omega_m) * u * u * u +
                    (1 - static_cast<long double>(p.omega_m)));
    else
      e = std::pow(u, 1 + static_cast<long double>(p.constant_q));
    return 1 / (clock ? u * e : e);
  };
  long double h = z / panels, total = f(0) + f(z);
  for (unsigned i = 1; i < panels; ++i)
    total += (i % 2 ? 4 : 2) * f(i * h);
  return total * h / 3;
}
Slot one(const Background &b, Query q, Policy p = {}) {
  auto result = b.evaluate_batch(std::span<const Query>(&q, 1), p);
  check(result.status == Status::ok && result.slots.size() == 1,
        "valid result batch shape");
  return result.slots[0];
}
long double mpc() {
  return 1e6L * 648000.L / std::numbers::pi_v<long double> * 149597870700.L;
}
} // namespace
int main() {
  try {
    const auto dh = 299792.458L / 70;
    const auto th = mpc() / (70.L * 1000);
    Parameters p{};
    auto base = prepare(p);
    check(base.status() == Status::ok &&
              base.constants_id == irred::constant_set_id,
          "model constants retained");
    for (double om : {0., 1.}) {
      p.omega_m = om;
      auto b = prepare(p);
      for (double z : {0., 1e-4, .1, 1., 5.}) {
        auto s = one(b, {z, z, Convention::geometric_same_redshift});
        check(s.status == Status::ok, "analytic cosmology slot");
        const long double u = 1 + static_cast<long double>(z);
        const auto i = om == 0 ? static_cast<long double>(z)
                               : -2 * std::expm1(-std::log(u) / 2);
        const auto t = om == 0 ? std::log(u)
                               : -(2.L / 3) * std::expm1(-1.5L * std::log(u));
        near(s.radial_integral, i, scale(i), "analytic radial integral");
        near(s.lookback_seconds / th, t, scale(t), "analytic clock");
        near(s.radial_mpc, dh * i, dh * scale(i),
             "independent physical distance");
        near(s.luminosity_mpc, u * u * s.angular_diameter_mpc,
             1e-8L + 2e-10L * std::abs(s.luminosity_mpc),
             "geometric reciprocity");
        near(s.deceleration_q, om == 0 ? -1 : .5, 1e-15, "analytic q");
        near(s.jerk, 1, 1e-15, "LCDM jerk");
        if (z == 0)
          check(s.evaluations == 0 && s.radial_mpc == 0 &&
                    s.volume_mpc3_per_sr_per_redshift == 0,
                "zero exact no callbacks");
      }
    }
    for (double q :
         {-2., -1. - 1e-8, -1., -1. + 1e-8, -1e-8, 0., 1e-8, .5, 2.}) {
      p.model = Model::constant_q_flat_v1;
      p.omega_m = 0;
      p.constant_q = q;
      auto b = prepare(p);
      check(b.status() == Status::ok, "constant q model");
      for (double z : {1e-4, .7, 5.}) {
        auto s = one(b, {z, z, Convention::geometric_same_redshift});
        check(s.status == Status::ok, "constant q finite");
        auto a = reference(p, z, false, 4096), c = reference(p, z, false, 8192),
             d = reference(p, z, false, 16384);
        maximum_reference_delta = std::max(
            maximum_reference_delta, static_cast<double>(std::abs(c - d)));
        reference_fraction =
            std::max(reference_fraction, std::abs(c - d) / scale(d));
        check(std::abs(c - d) <= scale(d) / 10 &&
                  std::abs(c - d) <= std::abs(a - c) + 1e-17L,
              "independent reference refinement");
        maximum_integral =
            std::max(maximum_integral,
                     static_cast<double>(std::abs(s.radial_integral - d)));
        near(s.radial_integral, d, scale(d),
             "q defining integral independent panels");
        auto t0 = reference(p, z, true, 4096), t1 = reference(p, z, true, 8192),
             t = reference(p, z, true, 16384);
        reference_fraction =
            std::max(reference_fraction, std::abs(t1 - t) / scale(t));
        check(std::abs(t1 - t) <= scale(t) / 10 &&
                  std::abs(t1 - t) <= std::abs(t0 - t1) + 1e-17L,
              "clock reference refinement");
        near(s.lookback_seconds / th, t, scale(t), "q clock limit both sides");
        auto dm = 5 / std::log(10.L) * std::log(s.radial_integral / d);
        maximum_mag = std::max(maximum_mag, static_cast<double>(std::abs(dm)));
        near(dm, 0, 1e-8, "magnitude error propagation");
        near(s.deceleration_q, q, 1e-15, "q parameter returned");
        near(s.jerk, static_cast<long double>(q) * (2 * q + 1), 1e-14,
             "constant q jerk");
      }
    }
    p = Parameters{};
    for (double om : {.01, .3, .99}) {
      p.omega_m = om;
      auto b = prepare(p);
      for (double z : {1e-4, .5, 2., 5.}) {
        auto s = one(b, {z, z, Convention::geometric_same_redshift});
        check(s.status == Status::ok, "LCDM fixed panelcase");
        auto a = reference(p, z, false, 4096), c = reference(p, z, false, 8192),
             d = reference(p, z, false, 16384);
        reference_fraction =
            std::max(reference_fraction, std::abs(c - d) / scale(d));
        check(std::abs(c - d) <= scale(d) / 10 &&
                  std::abs(c - d) <= std::abs(a - c) + 1e-17L,
              "LCDM reference refined");
        maximum_reference_delta = std::max(
            maximum_reference_delta, static_cast<double>(std::abs(c - d)));
        near(s.radial_integral, d, scale(d), "LCDM independent integration");
        auto dm = 5 / std::log(10.L) * std::log(s.radial_integral / d);
        maximum_mag = std::max(maximum_mag, static_cast<double>(std::abs(dm)));
        near(dm, 0, 1e-8, "LCDM magnitude budget");
      }
    }
    p = Parameters{};
    auto b = prepare(p);
    auto geo = one(b, {.5, .5, Convention::geometric_same_redshift}),
         obs = one(b, {.5, .7, Convention::released_zhd_zhel});
    check(obs.status == Status::ok && obs.source.z_observer == .7 &&
              obs.luminosity_equation_id != geo.luminosity_equation_id,
          "observer convention source identity");
    near(obs.luminosity_mpc / geo.luminosity_mpc, 1.7L / 1.5, 1e-14,
         "observer prefactor distinct");
    check(obs.angular_diameter_mpc == geo.angular_diameter_mpc &&
              obs.radial_mpc == geo.radial_mpc,
          "observer redshift doesn't alter radial or DA");
    p.h0_km_s_mpc = 100;
    auto higher =
        one(prepare(p), {.5, .5, Convention::geometric_same_redshift});
    near(higher.radial_mpc / geo.radial_mpc, .7L, 1e-14, "H0 distance scale");
    near(higher.lookback_seconds / geo.lookback_seconds, .7L, 1e-14,
         "H0 clock scale");
    near(higher.volume_mpc3_per_sr_per_redshift /
             geo.volume_mpc3_per_sr_per_redshift,
         .343L, 1e-14, "H0 volume scale");
    check(
        one(b, {1e-200, 1e-200, Convention::geometric_same_redshift}).status ==
            Status::numerical_failure,
        "physical volume underflow rejected");
    check(one(b, {std::numeric_limits<double>::denorm_min(),
                  std::numeric_limits<double>::denorm_min(),
                  Convention::geometric_same_redshift})
                  .status == Status::work_limit,
          "unbisectable minimumpositive interval truthful worklimit");
    check(one(b, {.5, std::numeric_limits<double>::max(),
                  Convention::released_zhd_zhel})
                  .status == Status::numerical_failure,
          "finite observer prefactor overflow is failure");
    check(one(b, {.5, .6, Convention::geometric_same_redshift}).status ==
              Status::incompatible_convention,
          "geometric mismatched redshift");
    check(one(b, {.5, -1, Convention::released_zhd_zhel}).status ==
              Status::incompatible_convention,
          "observer boundary");
    check(one(b, {-1, -1, Convention::released_zhd_zhel}).status ==
              Status::unsupported_domain,
          "negative expansion redshift");
    check(one(b, {5.1, 5.1, Convention::geometric_same_redshift}).status ==
              Status::unsupported_domain,
          "expansion domain");
    check(one(b, {NAN, .1, Convention::released_zhd_zhel}).status ==
              Status::invalid_input,
          "nonfinite query");
    check(one(b, {.1, .1, static_cast<Convention>(99)}).status ==
              Status::incompatible_convention,
          "unknown convention");
    p = Parameters{};
    p.model = static_cast<Model>(99);
    check(prepare(p).status() == Status::invalid_input, "unknown model");
    p = Parameters{};
    p.constant_q = 1;
    check(prepare(p).status() == Status::invalid_input,
          "inactive parameter cannot hide");
    p = Parameters{};
    p.omega_m = 1.1;
    check(prepare(p).status() == Status::unsupported_domain, "fraction domain");
    p = Parameters{};
    p.h0_km_s_mpc = 0;
    check(prepare(p).status() == Status::invalid_input, "H0 positive");
    p.h0_km_s_mpc = std::numeric_limits<double>::denorm_min();
    check(prepare(p).status() == Status::numerical_failure,
          "H0 conversion representability");
    Policy policy{};
    policy.maximum_queries = 0;
    Query query{.5, .5, Convention::geometric_same_redshift};
    check(b.evaluate_batch(std::span<const Query>(&query, 1), policy).status ==
              Status::work_limit,
          "batch resource limit before allocation");
    policy = Policy{};
    policy.integration.max_evaluations = 2;
    check(b.evaluate_batch(std::span<const Query>(&query, 1), policy).status ==
              Status::invalid_input,
          "invalid integration policy");
    policy = Policy{};
    policy.maximum_total_evaluations = 3;
    std::array<Query, 3> mixed{{{.5, .5, Convention::geometric_same_redshift},
                                {0, 0, Convention::geometric_same_redshift},
                                {.1, .1, Convention::geometric_same_redshift}}};
    auto work = b.evaluate_batch(mixed, policy);
    check(work.status == Status::ok &&
              work.slots[0].status == Status::work_limit &&
              work.slots[1].status == Status::ok &&
              work.slots[2].status == Status::work_limit,
          "global work budget mixedslots and exactzero");
    std::size_t total = 0;
    for (auto &s : work.slots)
      total += s.evaluations;
    check(total <= 3, "actual counted callbacks budget");
    check(b.evaluate_batch({}, Policy{}).status == Status::ok,
          "empty batch defined");
    std::printf("{\"budget_fractions\":{\"radial\":%.17Lg,\"clock\":%.17Lg,"
                "\"distance\":%.17Lg,\"magnitude\":%.17g,\"reference_"
                "refinement\":%.17Lg}}\n",
                radial_fraction, clock_fraction, distance_fraction,
                maximum_mag / 1e-8, reference_fraction);
    std::printf(
        "{\"suite\":\"independent_background\",\"checks\":%d,\"maximum_"
        "integral_error\":%.17g,\"maximum_magnitude_error\":%.17g,\"maximum_"
        "reference_refinement_delta\":%.17g,\"passed\":true}\n",
        checks, maximum_integral, maximum_mag, maximum_reference_delta);
    return 0;
  } catch (const std::exception &e) {
    std::fprintf(stderr, "%s\n", e.what());
    return 1;
  }
}
