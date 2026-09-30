// Native owner: analytic antiderivatives + separately split scale-factor
// Simpson. Budgets fixed before implementation: I/H0-clock abs2e-15+rel2e-10;
// E abs2e-14+rel2e-12; reference refinement<=10%; magnitude abs1e-8.
// Shared libm/equations are disclosed; this is not independent-author evidence.
#include "fixtures/piecewise_historical.hpp"
#include "irred/piecewise_background.hpp"
#include <array>
#include <bit>
#include <cmath>
#include <cstdio>
#include <limits>
#include <numbers>
#include <stdexcept>
using namespace irred::cosmology;
namespace {
int checks = 0;
long double worst = 0;
void check(bool b, const char *why) {
  ++checks;
  if (!b)
    throw std::runtime_error(why);
}
void near(double x, long double y) {
  auto budget = 2e-15L + 2e-10L * std::abs(y);
  worst = std::max(worst, std::abs(x - y) / budget);
  check(std::abs(x - y) <= budget, "integral budget");
}
EvaluationPolicy policy() {
  EvaluationPolicy p;
  p.maximum_queries = 100;
  p.maximum_segment_visits = 10000;
  p.maximum_callbacks = 100000;
  p.maximum_native_bytes = 1 << 20;
  p.integration = irred::numerics::IntegrationPolicy{1e-12, 1e-12, 100000, 30};
  return p;
}
Expansion fixed(std::array<double, 5> q) {
  return prepare(FixedFiveBinQ(q), FlatFLRW{});
}
Slot one(const Expansion &b, double z) {
  std::array<Request, 1> q{
      Request(z, 63, Observer{z, Convention::geometric_same_redshift},
              PhysicalScale(70))};
  return b.evaluate(q, policy()).slots[0];
}
long double E_at(const std::array<double, 5> &q, long double z) {
  long double logarithm = 0;
  for (int k = 0; k < 5; ++k) {
    long double lo = piecewise_q_edges[k],
                hi = std::min(z, (long double)piecewise_q_edges[k + 1]);
    if (hi > lo)
      logarithm += (1 + (long double)q[k]) * std::log1p((hi - lo) / (1 + lo));
  }
  return std::exp(logarithm);
}
long double reference(const std::array<double, 5> &q, double z, bool clock,
                      int panels) {
  long double sum = 0;
  for (int k = 0; k < 5; ++k) {
    long double lo = piecewise_q_edges[k],
                hi = std::min((long double)z,
                              (long double)piecewise_q_edges[k + 1]);
    if (hi <= lo)
      continue;
    auto a0 = 1 / (1 + hi), a1 = 1 / (1 + lo), h = (a1 - a0) / panels;
    for (int i = 0; i <= panels; ++i) {
      auto a = a0 + i * h;
      auto e = E_at(q, 1 / a - 1);
      sum += (i == 0 || i == panels ? 1
              : i % 2               ? 4
                                    : 2) *
             h / 3 / (a * a * e) * (clock ? a : 1);
    }
  }
  return sum;
}

// Current same-platform audit mode; historical bit receipts remain immutable.
void transcript() {
  for (auto q : {std::array<double, 5>{0, 0, 0, 0, 0},
                 std::array<double, 5>{-3, 2, -3, 2, -3},
                 std::array<double, 5>{-1, -1, -1, -1, -1}}) {
    auto b = fixed(q);
    for (double edge : piecewise_q_edges)
      for (double z : {std::nextafter(edge, -INFINITY), edge,
                       std::nextafter(edge, INFINITY)})
        if (z >= 0 && z <= 2.5) {
          auto s = one(b, z);
          std::printf("%a %u ", z, (unsigned)s.admission_status);
          std::printf("%u %a %u %a %u %a ", (unsigned)s.expansion.availability,
                      s.expansion.value ? s.expansion.value->expansion_E : 0.,
                      (unsigned)s.radial.availability,
                      s.radial.value ? s.radial.value->integral : 0.,
                      (unsigned)s.clock.availability,
                      s.clock.value ? s.clock.value->integral : 0.);
          std::printf("%u %u\n", (unsigned)s.kinematics.value->q_convention,
                      (unsigned)s.kinematics.value->jerk_availability);
        }
  }
}
} // namespace
int main(int argc, char **argv) {
  try {
    if (argc == 2 && std::string_view(argv[1]) == "--piecewise-transcript") {
      transcript();
      return 0;
    }
    for (const auto &reference : irred::test_reference::piecewise_historical) {
      auto b = fixed(reference.q);
      auto s = one(b, reference.z);
      check(s.radial.value && s.clock.value && s.kinematics.value,
            "historical fixture geometry");
      near(s.radial.value->integral, reference.I);
      check(std::abs(s.expansion.value->expansion_E - reference.E) <=
                2e-14 + 2e-12 * std::abs(reference.E),
            "historical E");
      check(s.kinematics.value->q == reference.assigned_q,
            "historical right-bin assignment");
    }
    auto mpc = 1e6L * 648000 / std::numbers::pi_v<long double> * 149597870700.L;
    auto time = mpc / 70000;
    for (double q :
         {-3., -1., -1 - 1e-10, -1 + 1e-10, -1e-10, 0., 1e-10, .5, 2.}) {
      std::array<double, 5> values;
      values.fill(q);
      auto b = fixed(values);
      check(b.status() == Status::ok, "constant prepare");
      for (double z : {0., .0001, .1, .3, .6, 1., 2.5}) {
        auto s = one(b, z);
        check(s.radial.value && s.clock.value && s.kinematics.value,
              "constant geometry");
        auto l = std::log1p((long double)z);
        auto I = q == 0 ? l : -std::expm1(-q * l) / q;
        auto T = q == -1 ? l
                         : -std::expm1(-(1 + (long double)q) * l) /
                               (1 + (long double)q);
        if (q == -3)
          I = (long double)z + (long double)z * z + (long double)z * z * z / 3;
        if (q == -1)
          I = z;
        if (q == 0)
          T = (long double)z / (1 + (long double)z);
        near(s.radial.value->integral, I);
        near(*s.clock.value->lookback_seconds.value / time, T);
        check(std::abs(s.expansion.value->expansion_E -
                       std::exp((1 + (long double)q) * l)) <=
                  2e-14L + 2e-12L * std::exp((1 + (long double)q) * l),
              "E budget");
        check(s.kinematics.value->jerk.has_value(),
              "no derivative jump for equal bins");
        check(std::abs(*s.kinematics.value->jerk -
                       (long double)q * (2 * (long double)q + 1)) <= 1e-12L,
              "constant jerk");
      }
    }
    for (auto q : {std::array<double, 5>{-.4, -.4, -.2, .1, .3},
                   std::array<double, 5>{-1, 0, -1, 0, -1},
                   std::array<double, 5>{-3, 2, -3, 2, -3}}) {
      auto b = fixed(q);
      for (double z : {.0001, .1, .3, .6, 1., 2.5}) {
        auto s = one(b, z);
        check(s.radial.value && s.clock.value && s.kinematics.value,
              "piecewise geometry");
        for (bool clock : {false, true}) {
          auto a = reference(q, z, clock, 1024),
               c = reference(q, z, clock, 2048);
          check(std::abs(a - c) < .1L * (2e-15L + 2e-10L * std::abs(c)),
                "split reference refinement");
          near(clock ? *s.clock.value->lookback_seconds.value / time
                     : s.radial.value->integral,
               c);
        }
      }
      for (std::size_t k = 1; k < 5; ++k) {
        auto edge = piecewise_q_edges[k];
        auto s = one(b, edge), l = one(b, std::nextafter(edge, 0.)),
             r = one(b, std::nextafter(edge, 3.));
        check(s.kinematics.value->q == q[k] &&
                  l.kinematics.value->q == q[k - 1] &&
                  r.kinematics.value->q == q[k],
              "right-bin assignment");
        near(s.radial.value->integral, l.radial.value->integral);
        near(s.radial.value->integral, r.radial.value->integral);
        near(*s.clock.value->lookback_seconds.value / time,
             *l.clock.value->lookback_seconds.value / time);
        near(*s.clock.value->lookback_seconds.value / time,
             *r.clock.value->lookback_seconds.value / time);
        check(std::abs(s.expansion.value->expansion_E -
                       r.expansion.value->expansion_E) < 1e-12,
              "right E continuity");
        check(std::abs(s.expansion.value->expansion_E -
                       l.expansion.value->expansion_E) < 1e-12,
              "E continuity");
        if (q[k] != q[k - 1])
          check(!s.kinematics.value->jerk &&
                    s.kinematics.value->q_convention ==
                        QConvention::right_limit_at_internal_jump,
                "undefined ordinary jerk");
      }
    }

    Kinematics unassessed;
    check(unassessed.q_convention == QConvention::not_assessed &&
              unassessed.jerk_availability == JerkAvailability::not_assessed,
          "default semantics");
    auto b = fixed({0, 0, 0, 0, 0});
    auto z0 = one(b, 0);
    check(z0.kinematics.value->q0_within_piecewise_model == 0 &&
              z0.kinematics.value->q_convention ==
                  QConvention::right_limit_at_zero,
          "q0 model convention");
    check(one(b, 2.5).kinematics.value->q_convention ==
              QConvention::left_limit_at_final_endpoint,
          "last endpoint");
    check(one(b, 1e-50).radial.value->integral > 0, "small normal z");
    // Current group isolation deliberately replaces old bundle-wide failure.
    auto range = one(b, std::numeric_limits<double>::min());
    check(range.radial.value && range.kinematics.value &&
              range.physical.availability == Availability::failed,
          "tiny volume failure preserves radial/kinematics");
    auto sub = fixed({std::numeric_limits<double>::denorm_min(), 0, 0, 0, 0});
    check(sub.status() == Status::ok &&
              one(sub, 0).kinematics.value->q ==
                  std::numeric_limits<double>::denorm_min(),
          "retained subnormal q");
    check(one(b, 2.500001).admission_status == Status::unsupported_domain,
          "no extrapolation");
    check(fixed({-3.1, 0, 0, 0, 0}).status() == Status::unsupported_domain,
          "q domain");
    check(fixed({NAN, 0, 0, 0, 0}).status() == Status::invalid_input,
          "nonfinite q");
    for (double h : {0., std::numeric_limits<double>::min()}) {
      std::array<Request, 1> q{
          Request(.1, 63, Observer{.1, Convention::geometric_same_redshift},
                  PhysicalScale(h))};
      auto r = b.evaluate(q, policy()).slots[0];
      check(r.radial.value && r.kinematics.value &&
                r.physical.availability == Availability::failed,
            "scale fails only requested dimensional projection");
      check(r.clock.value && r.clock.value->lookback_seconds.availability ==
                                 Availability::failed,
            "nested time failure preserves clock");
    }
    std::array<Request, 4> invalid{
        Request(.1, 63, Observer{.2, Convention::geometric_same_redshift},
                PhysicalScale(70)),
        Request(.1, 63, Observer{-1, Convention::released_zhd_zhel},
                PhysicalScale(70)),
        Request(.1, 63, Observer{.1, (Convention)99}, PhysicalScale(70)),
        Request(.1, 63, Observer{NAN, Convention::released_zhd_zhel},
                PhysicalScale(70))};
    for (const auto &s : b.evaluate(invalid, policy()).slots) {
      check(s.luminosity_shape.availability == Availability::failed &&
                s.kinematics.value,
            "observer failure preserves unrelated kinematics");
    }
    std::array<Request, 1> overflow{
        Request(1, 63,
                Observer{std::numeric_limits<double>::max(),
                         Convention::released_zhd_zhel},
                PhysicalScale(70))};
    auto r = b.evaluate(overflow, policy()).slots[0];
    check(r.physical.availability == Availability::failed && r.kinematics.value,
          "physical overflow isolation");
    check(!Expansion::workspace_payload_bound(SIZE_MAX).has_value(),
          "count-only overflow guard without an invalid span");
    std::array<Request, 2> qs{
        Request(.3, 63, Observer{.4, Convention::released_zhd_zhel},
                PhysicalScale(70)),
        Request(1, 63, Observer{1, Convention::geometric_same_redshift},
                PhysicalScale(70))};
    auto limitedpolicy = policy();
    limitedpolicy.maximum_segment_visits = 2;
    auto limited = b.evaluate(qs, limitedpolicy);
    check(limited.slots[0].radial.value &&
              limited.slots[1].radial.status == Status::work_limit &&
              limited.slots[1].kinematics.value,
          "global segment cap isolated");
    check(*limited.slots[0].luminosity_shape.value !=
              *one(b, .3).luminosity_shape.value,
          "observer prefactor retained");
    std::printf("piecewise current owner %d PASS max_fraction %.17Lg\n", checks,
                worst);
  } catch (const std::exception &e) {
    std::fprintf(stderr, "%s\n", e.what());
    return 1;
  }
}
