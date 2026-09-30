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
// Local independent integrand parameters, not a compatibility model API.
struct ReferenceParameters {
  bool lcdm = true;
  double omega_m = .3, q = 0;
};
long double reference(ReferenceParameters p, long double z, bool clock,
                      unsigned panels) {
  auto f = [&](long double x) {
    const auto u = 1 + x;
    long double e;
    if (p.lcdm)
      e = std::sqrt(static_cast<long double>(p.omega_m) * u * u * u +
                    (1 - static_cast<long double>(p.omega_m)));
    else
      e = std::pow(u, 1 + static_cast<long double>(p.q));
    return 1 / (clock ? u * e : e);
  };
  long double h = z / panels, total = f(0) + f(z);
  for (unsigned i = 1; i < panels; ++i)
    total += (i % 2 ? 4 : 2) * f(i * h);
  return total * h / 3;
}
EvaluationPolicy policy() {
  return {{irred::numerics::IntegrationPolicy{1e-13, 1e-12, 100000, 30}},
          64,
          1000000,
          2000,
          1000000};
}
Expansion model(ReferenceParameters p) {
  return p.lcdm ? prepare(LCDM(p.omega_m), FlatFLRW{})
                : prepare(ConstantQ(p.q), FlatFLRW{});
}
Request full(double z, double observer, double h = 70,
             Convention c = Convention::geometric_same_redshift) {
  return {z, 63, Observer{observer, c}, PhysicalScale(h)};
}
Slot one(const Expansion &e, Request q, EvaluationPolicy p = policy()) {
  auto result = e.evaluate(std::span<const Request>(&q, 1), p);
  check(result.status == Status::ok && result.slots.size() == 1,
        "valid result batch shape");
  return result.slots[0];
}
void available(const Slot &s) {
  check(s.radial.value && s.expansion.value && s.clock.value &&
            s.physical.value && s.luminosity_shape.value && s.kinematics.value,
        "all requested geometric groups available");
}
long double mpc() {
  return 1e6L * 648000.L / std::numbers::pi_v<long double> * 149597870700.L;
}
} // namespace
int main() {
  try {
    const auto dh = 299792.458L / 70, th = mpc() / (70.L * 1000);
    ReferenceParameters p;
    check(model(p).status() == Status::ok &&
              Expansion::constants_id == irred::constant_set_id,
          "constants identity");
    for (double om : {0., 1.}) {
      p.omega_m = om;
      auto b = model(p);
      for (double z : {0., 1e-4, .1, 1., 5.}) {
        auto s = one(b, full(z, z));
        available(s);
        const long double u = 1 + (long double)z;
        const auto i =
            om == 0 ? (long double)z : -2 * std::expm1(-std::log(u) / 2);
        const auto t = om == 0 ? std::log(u)
                               : -(2.L / 3) * std::expm1(-1.5L * std::log(u));
        near(s.radial.value->integral, i, scale(i), "analytic radial integral");
        near(*s.clock.value->lookback_seconds.value / th, t, scale(t),
             "analytic clock");
        near(s.physical.value->radial_mpc, dh * i, dh * scale(i),
             "independent physical distance");
        near(s.physical.value->luminosity_mpc,
             u * u * s.physical.value->angular_diameter_mpc,
             1e-8L + 2e-10L * std::abs(s.physical.value->luminosity_mpc),
             "geometric reciprocity");
        near(s.kinematics.value->q, om == 0 ? -1 : .5, 1e-15, "analytic q");
        near(*s.kinematics.value->jerk, 1, 1e-15, "LCDM jerk");
        if (z == 0)
          check(s.physical.value->radial_mpc == 0 &&
                    s.physical.value->volume_mpc3_per_sr_per_redshift == 0,
                "zero physical outputs");
      }
    }
    for (double q :
         {-2., -1. - 1e-8, -1., -1. + 1e-8, -1e-8, 0., 1e-8, .5, 2.}) {
      p = {false, 0, q};
      auto b = model(p);
      check(b.status() == Status::ok, "constant q model");
      for (double z : {1e-4, .7, 5.}) {
        auto s = one(b, full(z, z));
        available(s);
        auto a = reference(p, z, false, 4096), c = reference(p, z, false, 8192),
             d = reference(p, z, false, 16384);
        maximum_reference_delta =
            std::max(maximum_reference_delta, (double)std::abs(c - d));
        reference_fraction =
            std::max(reference_fraction, std::abs(c - d) / scale(d));
        check(std::abs(c - d) <= scale(d) / 10 &&
                  std::abs(c - d) <= std::abs(a - c) + 1e-17L,
              "independent reference refinement");
        maximum_integral = std::max(
            maximum_integral, (double)std::abs(s.radial.value->integral - d));
        near(s.radial.value->integral, d, scale(d),
             "q defining integral independent panels");
        auto t0 = reference(p, z, true, 4096), t1 = reference(p, z, true, 8192),
             t = reference(p, z, true, 16384);
        reference_fraction =
            std::max(reference_fraction, std::abs(t1 - t) / scale(t));
        check(std::abs(t1 - t) <= scale(t) / 10 &&
                  std::abs(t1 - t) <= std::abs(t0 - t1) + 1e-17L,
              "clock reference refinement");
        near(*s.clock.value->lookback_seconds.value / th, t, scale(t),
             "q clock limit both sides");
        auto dm = 5 / std::log(10.L) * std::log(s.radial.value->integral / d);
        maximum_mag = std::max(maximum_mag, (double)std::abs(dm));
        near(dm, 0, 1e-8, "magnitude error propagation");
        near(s.kinematics.value->q, q, 1e-15, "q parameter returned");
        near(*s.kinematics.value->jerk, (long double)q * (2 * q + 1), 1e-14,
             "constant q jerk");
      }
    }
    p = {};
    for (double om : {.01, .3, .99}) {
      p.omega_m = om;
      auto b = model(p);
      for (double z : {1e-4, .5, 2., 5.}) {
        auto s = one(b, full(z, z));
        available(s);
        auto a = reference(p, z, false, 4096), c = reference(p, z, false, 8192),
             d = reference(p, z, false, 16384);
        reference_fraction =
            std::max(reference_fraction, std::abs(c - d) / scale(d));
        check(std::abs(c - d) <= scale(d) / 10 &&
                  std::abs(c - d) <= std::abs(a - c) + 1e-17L,
              "LCDM reference refined");
        maximum_reference_delta =
            std::max(maximum_reference_delta, (double)std::abs(c - d));
        near(s.radial.value->integral, d, scale(d),
             "LCDM independent integration");
        auto dm = 5 / std::log(10.L) * std::log(s.radial.value->integral / d);
        maximum_mag = std::max(maximum_mag, (double)std::abs(dm));
        near(dm, 0, 1e-8, "LCDM magnitude budget");
      }
    }
    p = {};
    auto b = model(p);
    auto geo = one(b, full(.5, .5)),
         obs = one(b, full(.5, .7, 70, Convention::released_zhd_zhel));
    available(obs);
    check(obs.source.observer->redshift == .7, "observer attempt retained");
    near(obs.physical.value->luminosity_mpc /
             geo.physical.value->luminosity_mpc,
         1.7L / 1.5, 1e-14, "observer prefactor distinct");
    check(obs.physical.value->angular_diameter_mpc ==
                  geo.physical.value->angular_diameter_mpc &&
              obs.physical.value->radial_mpc == geo.physical.value->radial_mpc,
          "observer doesn't change radial or DA");
    auto higher = one(b, full(.5, .5, 100));
    available(higher);
    near(higher.physical.value->radial_mpc / geo.physical.value->radial_mpc,
         .7L, 1e-14, "H0 distance scale");
    near(*higher.clock.value->lookback_seconds.value /
             *geo.clock.value->lookback_seconds.value,
         .7L, 1e-14, "H0 clock scale");
    near(higher.physical.value->volume_mpc3_per_sr_per_redshift /
             geo.physical.value->volume_mpc3_per_sr_per_redshift,
         .343L, 1e-14, "H0 volume scale");
    auto tiny = one(b, full(1e-200, 1e-200));
    check(tiny.physical.status == Status::numerical_failure &&
              !tiny.physical.value,
          "physical volume underflow rejected");
    auto denorm = one(b, full(std::numeric_limits<double>::denorm_min(),
                              std::numeric_limits<double>::denorm_min()));
    check(denorm.radial.status == Status::work_limit,
          "unbisectable interval work limit");
    auto extreme_observer = one(b, full(.5, std::numeric_limits<double>::max(),
                                        70, Convention::released_zhd_zhel));
    check(extreme_observer.luminosity_shape.value &&
              extreme_observer.physical.status == Status::numerical_failure,
          "physical observer overflow isolated from finite shape");
    for (auto r :
         {full(.5, .6), full(.5, -1, 70, Convention::released_zhd_zhel),
          full(.1, .1, 70, (Convention)99)})
      check(one(b, r).luminosity_shape.status ==
                Status::incompatible_convention,
            "observer convention admission");
    for (auto z : {-1., 5.1})
      check(one(b, full(z, z)).radial.status == Status::unsupported_domain,
            "redshift domain");
    check(one(b, full(NAN, .1)).radial.status == Status::invalid_input,
          "nonfinite z");
    check(prepare(LCDM(1.1), FlatFLRW{}).status() == Status::unsupported_domain,
          "fraction domain");
    check(prepare(ConstantQ(NAN), FlatFLRW{}).status() == Status::invalid_input,
          "nonfinite active field");
    auto badscale = one(b, full(.5, .5, 0));
    check(badscale.physical.status == Status::invalid_input &&
              badscale.radial.value && badscale.luminosity_shape.value,
          "scale failure isolated");
    auto pp = policy();
    pp.maximum_queries = 0;
    auto rq = full(.5, .5);
    check(b.evaluate(std::span<const Request>(&rq, 1), pp).status ==
              Status::work_limit,
          "preallocation query cap");
    pp = policy();
    pp.integration->max_evaluations = 2;
    check(b.evaluate(std::span<const Request>(&rq, 1), pp).status ==
              Status::invalid_input,
          "invalid numerical policy");
    pp = policy();
    pp.maximum_callbacks = 3;
    std::array<Request, 3> mixed{full(.5, .5), full(0, 0), full(.1, .1)};
    auto work = b.evaluate(mixed, pp);
    check(work.status == Status::ok &&
              work.slots[0].radial.status == Status::work_limit &&
              work.slots[1].radial.value &&
              work.slots[2].radial.status == Status::work_limit &&
              work.work.callbacks <= 3,
          "one global callback cap including failed work");
    // E-only is observer/H0/integration independent, including q-bin edges.
    EvaluationPolicy zero{{}, 16, 0, 0, 100000};
    Request eonly{.5, 32};
    auto ev = b.evaluate(std::span<const Request>(&eonly, 1), zero);
    check(ev.status == Status::ok && ev.slots[0].expansion.value &&
              ev.work.callbacks == 0 && ev.work.segment_visits == 0 &&
              !ev.slots[0].radial.value,
          "E only zero work");
    eonly.physical_scale = PhysicalScale(0);
    ev = b.evaluate(std::span<const Request>(&eonly, 1), zero);
    check(ev.slots[0].expansion.value &&
              ev.slots[0].expansion.value->h_km_s_mpc.status ==
                  Status::invalid_input,
          "invalid H scale doesn't poison E");
    for (uint32_t mask : {0u, 64u}) {
      Request bad{.5, mask};
      auto v = b.evaluate(std::span<const Request>(&bad, 1), policy());
      check(v.status == Status::invalid_input && v.slots.empty() &&
                v.work.callbacks == 0,
            "mask structural rejection");
    }
    std::array<Request, 3> duplicates{
        {Request{.5, 1},
         Request{.5, 2, Observer{.7, Convention::released_zhd_zhel}},
         Request{.5, 1}}};
    auto reused = b.evaluate(duplicates, policy());
    check(reused.nodes.size() == 1 && reused.slots.size() == 3 &&
              reused.slots[0].node_index == reused.slots[2].node_index &&
              reused.slots[0].radial.value->integral ==
                  reused.slots[2].radial.value->integral,
          "exact z reuse preserves rows");
    auto single =
        b.evaluate(std::span<const Request>(duplicates.data(), 1), policy());
    check(single.work.callbacks == reused.work.callbacks,
          "duplicate radial work charged once");
    check(b.evaluate({}, policy()).status == Status::ok, "empty batch defined");
    auto empty_bound = Expansion::workspace_payload_bound(0);
    check(empty_bound && *empty_bound == 0,
          "empty workspace has no heap payload");
    check(
        !Expansion::workspace_payload_bound(std::numeric_limits<size_t>::max()),
        "workspace multiplication overflow rejected without allocation");
    for (size_t count : {3u, 7u}) {
      std::vector<Request> requests;
      for (size_t i = 0; i < count; ++i)
        requests.emplace_back(.1 + i * .01, 32);
      auto bound = Expansion::workspace_payload_bound(count);
      check(bound && *bound > 0, "nonpower workspace bound available");
      EvaluationPolicy exact{{}, count, 0, 0, *bound};
      auto admitted = b.evaluate(requests, exact);
      check(admitted.status == Status::ok && admitted.slots.size() == count &&
                admitted.nodes.size() == count && admitted.work.callbacks == 0,
            "exact authoritative workspace bound admitted");
      --exact.maximum_native_bytes;
      auto rejected = b.evaluate(requests, exact);
      check(rejected.status == Status::work_limit && rejected.slots.empty(),
            "workspace cap minus one rejected before payload allocation");
    }
    for (double h : {0., -1., std::numeric_limits<double>::infinity(),
                     std::numeric_limits<double>::denorm_min()}) {
      Request clock_only{.5, 4, {}, PhysicalScale(h)};
      auto retained = one(b, clock_only);
      check(retained.clock.value &&
                retained.clock.availability == Availability::available &&
                retained.clock.value->integral > 0,
            "dimensionless clock survives dimensional scale failure");
      const auto &seconds = retained.clock.value->lookback_seconds;
      check(seconds.availability == Availability::failed && !seconds.value &&
                seconds.numerical_status != irred::numerics::Status::ok,
            "nested physical time retains explicit failure");
      check(retained.radial.availability == Availability::not_requested &&
                retained.kinematics.availability == Availability::not_requested,
            "clock request does not force unrelated groups");
    }
    Request huge_h{
        .5, 4, {}, PhysicalScale(std::numeric_limits<double>::max())};
    auto tiny_time = one(b, huge_h);
    check(tiny_time.clock.value &&
              tiny_time.clock.value->lookback_seconds.value &&
              *tiny_time.clock.value->lookback_seconds.value > 0,
          "positive representable physical time at huge valid H0");
    std::printf(
        "{\"suite\":\"independent_background_current\",\"checks\":%d,\"radial_"
        "budget_fraction\":%.17Lg,\"clock_budget_fraction\":%.17Lg,\"distance_"
        "budget_fraction\":%.17Lg,\"magnitude_budget_fraction\":%.17g,"
        "\"reference_refinement_fraction\":%.17Lg,\"passed\":true}\n",
        checks, radial_fraction, clock_fraction, distance_fraction,
        maximum_mag / 1e-8, reference_fraction);
    return 0;
  } catch (const std::exception &e) {
    std::fprintf(stderr, "%s\n", e.what());
    return 1;
  }
}
