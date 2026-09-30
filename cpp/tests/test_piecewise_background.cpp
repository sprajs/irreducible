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
PiecewiseSlot one(const PiecewiseBackground &b, double z) {
  std::array<Query, 1> q{{{z, z, Convention::geometric_same_redshift}}};
  return b.evaluate_batch(q, {}).slots[0];
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
void transcript() {
  for (int model = 0; model < 3; ++model)
    for (double value : {0., .3, 1.}) {
      auto b = model == 2
                   ? prepare_cpl({70, value, -.9, .4})
                   : prepare({static_cast<Model>(model), 70,
                              model == 1 ? 0 : value, model == 1 ? value : 0});
      for (double z : {0., .0001, .1, .3, .6, 1., 2.5, 5.}) {
        std::array<Query, 1> qs{{{z, z, Convention::geometric_same_redshift}}};
        auto batch = b.evaluate_batch(qs, {});
        auto &s = batch.slots[0];
        std::printf("%u %u %zu ", (unsigned)s.status,
                    (unsigned)s.numerical_status, s.evaluations);
        for (double x :
             {s.expansion_E, s.h_km_s_mpc, s.radial_integral, s.radial_mpc,
              s.transverse_mpc, s.angular_diameter_mpc, s.luminosity_mpc,
              s.dimensionless_luminosity_shape, s.lookback_seconds,
              s.volume_mpc3_per_sr_per_redshift, s.deceleration_q, s.jerk,
              s.radial_integral_error, s.lookback_integral_error})
          std::printf("%a ", x);
        std::printf("\n");
      }
    }
}
// Same-platform before/after audit mode; no cross-platform libm bit promise.
void transcript_double_bits(double v) {
  std::printf("%016llx ",
              (unsigned long long)std::bit_cast<unsigned long long>(v));
}
void piecewise_transcript() {
  for (auto q : {std::array<double, 5>{0, 0, 0, 0, 0},
                 std::array<double, 5>{-3, 2, -3, 2, -3},
                 std::array<double, 5>{-1, -1, -1, -1, -1}}) {
    auto b = prepare_piecewise_q({70, q});
    std::vector<Query> qs;
    for (double e : piecewise_q_edges)
      for (double z :
           {std::nextafter(e, -INFINITY), e, std::nextafter(e, INFINITY)})
        qs.push_back({z, z, Convention::geometric_same_redshift});
    qs.push_back({1e-200, 1e-200, Convention::geometric_same_redshift});
    for (auto cap : {size_t(1000), size_t(5)}) {
      auto r = b.evaluate_batch(qs, {100, cap});
      std::printf("%u %zu\n", (unsigned)r.status, r.segments_processed);
      for (auto &s : r.slots) {
        std::printf("%u %u %zu %zu %u %u ", (unsigned)s.status,
                    (unsigned)s.numerical_status, s.bin, s.segments_processed,
                    (unsigned)s.q_convention, (unsigned)s.jerk_availability);
        auto &g = s.geometry;
        for (double x : {g.expansion_E, g.h_km_s_mpc, g.radial_integral,
                         g.radial_mpc, g.transverse_mpc, g.angular_diameter_mpc,
                         g.luminosity_mpc, g.dimensionless_luminosity_shape,
                         g.lookback_seconds, g.volume_mpc3_per_sr_per_redshift})
          transcript_double_bits(x);
        for (auto v : {s.assigned_q, s.jerk, s.q0_within_piecewise_model}) {
          std::printf("%u ", (unsigned)v.has_value());
          if (v)
            transcript_double_bits(*v);
        }
        std::puts("");
      }
    }
  }
}

} // namespace
int main(int argc, char **argv) {
  try {
    if (argc == 2 && std::string_view(argv[1]) == "--piecewise-transcript") {
      piecewise_transcript();
      return 0;
    }
    if (argc == 2 && std::string_view(argv[1]) == "--legacy-transcript") {
      transcript();
      return 0;
    }
    for (const auto &reference : irred::test_reference::piecewise_historical) {
      auto b = prepare_piecewise_q({70, reference.q});
      auto s = one(b, reference.z);
      check(s.status == Status::ok, "historical fixture geometry");
      near(s.geometry.radial_integral, reference.I);
      check(std::abs(s.geometry.expansion_E - reference.E) <=
                2e-14 + 2e-12 * std::abs(reference.E),
            "historical E");
      check(s.assigned_q == reference.assigned_q,
            "historical right-bin assignment");
    }
    auto mpc = 1e6L * 648000 / std::numbers::pi_v<long double> * 149597870700.L;
    auto time = mpc / 70000;
    for (double q :
         {-3., -1., -1 - 1e-10, -1 + 1e-10, -1e-10, 0., 1e-10, .5, 2.}) {
      std::array<double, 5> values;
      values.fill(q);
      auto b = prepare_piecewise_q({70, values});
      check(b.status() == Status::ok, "constant prepare");
      for (double z : {0., .0001, .1, .3, .6, 1., 2.5}) {
        auto s = one(b, z);
        check(s.status == Status::ok, "constant geometry");
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
        near(s.geometry.radial_integral, I);
        near(s.geometry.lookback_seconds / time, T);
        check(std::abs(s.geometry.expansion_E -
                       std::exp((1 + (long double)q) * l)) <=
                  2e-14L + 2e-12L * std::exp((1 + (long double)q) * l),
              "E budget");
        check(s.jerk.has_value(), "no derivative jump for equal bins");
        check(std::abs(*s.jerk - (long double)q * (2 * (long double)q + 1)) <=
                  1e-12L,
              "constant jerk");
      }
    }
    for (auto q : {std::array<double, 5>{-.4, -.4, -.2, .1, .3},
                   std::array<double, 5>{-1, 0, -1, 0, -1},
                   std::array<double, 5>{-3, 2, -3, 2, -3}}) {
      auto b = prepare_piecewise_q({70, q});
      for (double z : {.0001, .1, .3, .6, 1., 2.5}) {
        auto s = one(b, z);
        check(s.status == Status::ok, "piecewise geometry");
        for (bool clock : {false, true}) {
          auto a = reference(q, z, clock, 1024),
               c = reference(q, z, clock, 2048);
          check(std::abs(a - c) < .1L * (2e-15L + 2e-10L * std::abs(c)),
                "split reference refinement");
          near(clock ? s.geometry.lookback_seconds / time
                     : s.geometry.radial_integral,
               c);
        }
      }
      for (std::size_t k = 1; k < 5; ++k) {
        auto edge = piecewise_q_edges[k];
        auto s = one(b, edge), l = one(b, std::nextafter(edge, 0.)),
             r = one(b, std::nextafter(edge, 3.));
        check(s.assigned_q == q[k] && l.assigned_q == q[k - 1] &&
                  r.assigned_q == q[k],
              "right-bin assignment");
        near(s.geometry.radial_integral, l.geometry.radial_integral);
        near(s.geometry.radial_integral, r.geometry.radial_integral);
        near(s.geometry.lookback_seconds / time,
             l.geometry.lookback_seconds / time);
        near(s.geometry.lookback_seconds / time,
             r.geometry.lookback_seconds / time);
        check(std::abs(s.geometry.expansion_E - r.geometry.expansion_E) < 1e-12,
              "right E continuity");
        check(std::abs(s.geometry.expansion_E - l.geometry.expansion_E) < 1e-12,
              "E continuity");
        if (q[k] != q[k - 1])
          check(!s.jerk &&
                    s.q_convention == QConvention::right_limit_at_internal_jump,
                "undefined ordinary jerk");
      }
    }
    PiecewiseSlot unassessed;
    check(unassessed.q_convention == QConvention::not_assessed &&
              unassessed.jerk_availability == JerkAvailability::not_assessed,
          "default derivative semantics unassessed");
    auto b = prepare_piecewise_q({70, {0, 0, 0, 0, 0}});
    auto z0 = one(b, 0);
    check(z0.q0_within_piecewise_model == 0 &&
              z0.q_convention == QConvention::right_limit_at_zero,
          "q0 imposed equality");
    check(one(b, 2.5).q_convention == QConvention::left_limit_at_final_endpoint,
          "last endpoint");
    auto tiny = one(b, 1e-50);
    check(tiny.status == Status::ok && tiny.geometry.radial_integral > 0,
          "positive small normal z");
    auto range = one(b, std::numeric_limits<double>::min());
    check(range.status == Status::numerical_failure && !range.assigned_q &&
              !range.jerk,
          "tiny range failure no derived payload");
    check(range.q_convention == QConvention::not_assessed &&
              range.jerk_availability == JerkAvailability::not_assessed,
          "failed range derivative semantics unassessed");
    auto subq = prepare_piecewise_q(
        {70, {std::numeric_limits<double>::denorm_min(), 0, 0, 0, 0}});
    check(subq.status() == Status::ok &&
              one(subq, 0).assigned_q ==
                  std::numeric_limits<double>::denorm_min(),
          "retained subnormal q");
    check(one(b, 2.500001).status == Status::unsupported_domain,
          "no extrapolation");
    check(prepare_piecewise_q({70, {-3.1, 0, 0, 0, 0}}).status() ==
              Status::unsupported_domain,
          "q domain");
    check(prepare_piecewise_q({70, {NAN, 0, 0, 0, 0}}).status() ==
              Status::invalid_input,
          "nonfinite q");
    check(prepare_piecewise_q({0, {0, 0, 0, 0, 0}}).status() ==
              Status::invalid_input,
          "nonpositive H0");
    check(prepare_piecewise_q(
              {std::numeric_limits<double>::min(), {0, 0, 0, 0, 0}})
                  .status() == Status::numerical_failure,
          "H0 representational failure");
    std::array<Query, 4> invalid_queries{
        {{.1, .2, Convention::geometric_same_redshift},
         {.1, -1, Convention::released_zhd_zhel},
         {.1, .1, static_cast<Convention>(99)},
         {.1, NAN, Convention::released_zhd_zhel}}};
    auto invalid_batch = b.evaluate_batch(invalid_queries, {});
    for (std::size_t i = 0; i < invalid_queries.size(); ++i) {
      check(invalid_batch.slots[i].status ==
                (i == 3 ? Status::invalid_input
                        : Status::incompatible_convention),
            "query convention/domain");
      check(!invalid_batch.slots[i].assigned_q && !invalid_batch.slots[i].jerk,
            "failed query no derivative payload");
    }
    std::array<Query, 1> enormous_observer{
        {{1, std::numeric_limits<double>::max(),
          Convention::released_zhd_zhel}}};
    auto overflow = b.evaluate_batch(enormous_observer, {}).slots[0];
    check(overflow.status == Status::numerical_failure &&
              overflow.numerical_status == irred::numerics::Status::overflow &&
              !overflow.assigned_q,
          "geometry overflow cause");
    auto huge = std::span<const Query>((const Query *)nullptr,
                                       std::numeric_limits<std::size_t>::max());
    check(b.evaluate_batch(huge, {}).status == Status::work_limit,
          "pre-copy resource guard");
    std::array<Query, 2> qs{{{.3, .4, Convention::released_zhd_zhel},
                             {1, 1, Convention::geometric_same_redshift}}};
    auto limited = b.evaluate_batch(qs, {10, 2});
    check(limited.slots[0].status == Status::ok &&
              limited.slots[1].status == Status::work_limit,
          "global segment cap");
    check(limited.slots[0].geometry.dimensionless_luminosity_shape !=
              one(b, .3).geometry.dimensionless_luminosity_shape,
          "observer prefactor retained");
    std::printf("{\"suite\":\"piecewise_owner\",\"checks\":%d,\"maximum_"
                "integral_budget_fraction\":%.17Lg,\"passed\":true}\n",
                checks, worst);
  } catch (const std::exception &e) {
    std::fprintf(stderr, "%s\n", e.what());
    return 1;
  }
}
