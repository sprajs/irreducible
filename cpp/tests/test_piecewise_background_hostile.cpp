// Independent split composite Simpson integrals, rather than production
// segment antiderivatives. Samples do not certify arbitrary q(z) or priors.
#include "irred/piecewise_background.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
using namespace irred::cosmology;
namespace {
int checks = 0;
long double radial_ratio = 0, clock_ratio = 0, expansion_ratio = 0,
            refinement_ratio = 0;
void check(bool b, const char *why) {
  ++checks;
  if (!b)
    throw std::runtime_error(why);
}
long double expansion(long double z, const std::array<double, 5> &q) {
  long double e = 1;
  for (int i = 0; i < 5; ++i) {
    auto lo = (long double)piecewise_q_edges[i];
    if (z <= lo)
      break;
    auto hi = std::min(z, (long double)piecewise_q_edges[i + 1]);
    e *= std::pow((1 + hi) / (1 + lo), 1 + (long double)q[i]);
  }
  return e;
}
long double integral(double z, const std::array<double, 5> &q, int n,
                     bool clock) {
  long double sum = 0;
  for (int i = 0; i < 5; ++i) {
    long double lo = piecewise_q_edges[i],
                hi = std::min((long double)z,
                              (long double)piecewise_q_edges[i + 1]);
    if (hi <= lo)
      break;
    auto h = (hi - lo) / n;
    long double local = 0;
    for (int j = 0; j <= n; ++j) {
      auto x = lo + j * h;
      auto f = 1 / expansion(x, q);
      if (clock)
        f /= 1 + x;
      local += (j == 0 || j == n ? 1 : j % 2 ? 4 : 2) * f;
    }
    sum += h * local / 3;
  }
  return sum;
}
} // namespace
int main() {
  try {
    for (auto q : {std::array<double, 5>{0, 0, 0, 0, 0},
                   std::array<double, 5>{-1, -1, -1, -1, -1},
                   std::array<double, 5>{-.4, -.4, -.2, .1, .3},
                   std::array<double, 5>{-3, 2, -1e-10, 1e-10, -1 + 1e-10}}) {
      auto p = prepare_piecewise_q({70, q});
      check(p.status() == Status::ok, "valid provider");
      for (double z : {0., .03, .1, .3, .6, 1., 1.7, 2.5}) {
        Query query{z, z, Convention::geometric_same_redshift};
        auto b = p.evaluate_batch(std::span(&query, 1), {});
        auto &s = b.slots[0];
        check(s.status == Status::ok, "valid geometry");
        auto r = integral(z, q, 4096, false),
             fine = integral(z, q, 8192, false);
        check(std::abs(r - fine) <= .1L * (2e-15L + 2e-10L * std::abs(fine)),
              "reference refinement");
        check(std::abs(s.geometry.radial_integral - fine) <
                  2e-15L + 2e-10L * std::abs(fine),
              "independent radial");
        check(std::abs(s.geometry.expansion_E - expansion(z, q)) <
                  2e-14L + 2e-12L * expansion(z, q),
              "independent expansion");
        auto clock_coarse = integral(z, q, 4096, true);
        auto clock_fine = integral(z, q, 8192, true);
        check(std::abs(clock_coarse - clock_fine) <=
                  .1L * (2e-15L + 2e-10L * std::abs(clock_fine)),
              "independent clock refinement");
        // Hubble time independently reconstructed from the same declared IAU
        // parsec definition, not from the production geometry helper.
        const auto mpc_metres =
            1e6L * (648000.L / std::acos(-1.L)) * 149597870700.L;
        const auto expected_seconds = clock_fine * mpc_metres / 70000.L;
        check(std::abs(s.geometry.lookback_seconds - expected_seconds) <
                  (2e-15L + 2e-10L * std::abs(clock_fine)) * mpc_metres /
                      70000.L,
              "independent lookback units");
        radial_ratio =
            std::max(radial_ratio, std::abs(s.geometry.radial_integral - fine) /
                                       (2e-15L + 2e-10L * std::abs(fine)));
        expansion_ratio =
            std::max(expansion_ratio,
                     std::abs(s.geometry.expansion_E - expansion(z, q)) /
                         (2e-14L + 2e-12L * expansion(z, q)));
        clock_ratio =
            std::max(clock_ratio,
                     std::abs(s.geometry.lookback_seconds - expected_seconds) /
                         ((2e-15L + 2e-10L * std::abs(clock_fine)) *
                          mpc_metres / 70000.L));
        refinement_ratio = std::max(
            {refinement_ratio,
             std::abs(r - fine) / (.1L * (2e-15L + 2e-10L * std::abs(fine))),
             std::abs(clock_coarse - clock_fine) /
                 (.1L * (2e-15L + 2e-10L * std::abs(clock_fine)))});
        auto k = s.bin;
        bool jump =
            z > 0 && z < 2.5 && z == piecewise_q_edges[k] && q[k] != q[k - 1];
        check(s.jerk.has_value() != jump,
              "jerk availability follows actual jump");
        if (s.jerk)
          check(std::abs(*s.jerk - q[k] * (2 * q[k] + 1)) < 2e-14,
                "jerk ordinary value");
        if (z == 0)
          check(s.geometry.radial_mpc == 0 && s.geometry.luminosity_mpc == 0 &&
                    s.q0_within_piecewise_model == q[0],
                "zero geometry and one-sided q0");
      }
    }
    auto p = prepare_piecewise_q({70, {-1, 0, 1, -.5, .2}});
    for (int i = 1; i < 5; ++i)
      for (double z :
           {std::nextafter(piecewise_q_edges[i], 0.), piecewise_q_edges[i],
            std::nextafter(piecewise_q_edges[i], 3.)}) {
        Query q{z, z, Convention::geometric_same_redshift};
        auto s = p.evaluate_batch(std::span(&q, 1), {}).slots[0];
        check(s.status == Status::ok, "nextafter boundary geometry");
        check(s.assigned_q == p.parameters().q[s.bin], "retained side q");
      }
    Query queries[]{{2., 2., Convention::geometric_same_redshift},
                    {2., 2., Convention::geometric_same_redshift}};
    PiecewisePolicy cap;
    cap.maximum_segment_visits = 5;
    auto b = p.evaluate_batch(queries, cap);
    check(b.slots[0].status == Status::ok &&
              b.slots[1].status == Status::work_limit,
          "global segment cap");
    cap.maximum_queries = 0;
    check(p.evaluate_batch(queries, cap).slots.empty(),
          "count cap before allocation");
    Query normal{.4, .4, Convention::geometric_same_redshift};
    auto a = p.evaluate_batch(std::span(&normal, 1), {}).slots[0];
    auto scaled = prepare_piecewise_q({100, p.parameters().q})
                      .evaluate_batch(std::span(&normal, 1), {})
                      .slots[0];
    check(std::abs(scaled.geometry.radial_mpc / a.geometry.radial_mpc - .7) <
              2e-15,
          "H0 inverse distance scaling");
    check(std::abs(scaled.geometry.h_km_s_mpc / a.geometry.h_km_s_mpc -
                   100. / 70) < 2e-15,
          "H0 expansion scaling");
    for (double z : {1e-200, std::numeric_limits<double>::denorm_min()}) {
      Query tiny{z, z, Convention::geometric_same_redshift};
      auto failed = p.evaluate_batch(std::span(&tiny, 1), {}).slots[0];
      check(failed.status != Status::ok && !failed.assigned_q && !failed.jerk,
            "unrepresentable geometry has no derivative payload");
    }
    check(prepare_piecewise_q({70, {0, 0, 3, 0, 0}}).status() ==
              Status::unsupported_domain,
          "q domain");
    check(prepare_piecewise_q({0, {0, 0, 0, 0, 0}}).status() != Status::ok,
          "invalid H0");
    for (double derivative :
         {std::numeric_limits<double>::denorm_min(),
          -std::numeric_limits<double>::denorm_min(), 0.0, -0.0}) {
      auto tiny_q = prepare_piecewise_q({70, {derivative, 0, 0, 0, 0}});
      Query query{.05, .05, Convention::geometric_same_redshift};
      auto s = tiny_q.evaluate_batch(std::span(&query, 1), {}).slots[0];
      check(s.status == Status::ok && s.assigned_q == derivative && s.jerk,
            "dimensionless subnormal q does not erase valid geometry");
      check(*s.jerk == derivative * (2 * derivative + 1),
            "dimensionless subnormal jerk retained");
    }
    for (Query query : {Query{.2, .3, Convention::geometric_same_redshift},
                        Query{.2, -1, Convention::released_zhd_zhel},
                        Query{.2, .2, static_cast<Convention>(99)},
                        Query{3, 3, Convention::geometric_same_redshift},
                        Query{std::numeric_limits<double>::quiet_NaN(), .2,
                              Convention::released_zhd_zhel}}) {
      auto s = p.evaluate_batch(std::span(&query, 1), {}).slots[0];
      check(s.status != Status::ok && !s.assigned_q && !s.jerk &&
                s.geometry.radial_mpc == 0 &&
                s.q_convention == QConvention::not_assessed &&
                s.jerk_availability == JerkAvailability::not_assessed,
            "invalid source has no finite geometric or derivative payload");
    }
    Query released{.4, .45, Convention::released_zhd_zhel};
    auto released_slot = p.evaluate_batch(std::span(&released, 1), {}).slots[0];
    check(released_slot.status == Status::ok &&
              released_slot.source.convention == Convention::released_zhd_zhel,
          "released observer convention retained");
    check(std::abs(released_slot.geometry.luminosity_mpc /
                       a.geometry.luminosity_mpc -
                   1.45 / 1.4) < 2e-15 &&
              released_slot.geometry.angular_diameter_mpc ==
                  a.geometry.angular_diameter_mpc,
          "observer luminosity prefactor distinct from geometric DA");
    Query tiny_normal{1e-90, 1e-90, Convention::geometric_same_redshift};
    check(p.evaluate_batch(std::span(&tiny_normal, 1), {}).slots[0].status ==
              Status::ok,
          "tiny but normal physical bundle remains valid");
    PiecewiseBackground default_provider;
    check(default_provider.status() == Status::invalid_input &&
              default_provider.parameters().h0_km_s_mpc == 0,
          "default provider inspectable and unprepared");
    auto no_owner = default_provider.evaluate_batch(std::span(&normal, 1), {});
    check(no_owner.status == Status::invalid_input && no_owner.slots.empty(),
          "default provider yields no scientific slots");
    PiecewiseSlot default_slot;
    check(!default_slot.assigned_q && !default_slot.jerk &&
              default_slot.q_convention == QConvention::not_assessed &&
              default_slot.jerk_availability == JerkAvailability::not_assessed,
          "default derivative assessment unset");
    auto empty = p.evaluate_batch({}, {});
    check(empty.status == Status::ok && empty.slots.empty() &&
              empty.segments_processed == 0,
          "empty valid batch is completed without work");
    auto copied = p;
    auto copied_result =
        copied.evaluate_batch(std::span(&normal, 1), {}).slots[0];
    check(copied_result.geometry.radial_integral ==
                  a.geometry.radial_integral &&
              copied.parameters().q == p.parameters().q,
          "copied immutable provider preserves numerical/source identity");
    std::printf("piecewise independent checks %d PASS radial %.18Lg clock "
                "%.18Lg expansion %.18Lg refinement %.18Lg\n",
                checks, radial_ratio, clock_ratio, expansion_ratio,
                refinement_ratio);
    return 0;
  } catch (const std::exception &e) {
    std::fprintf(stderr, "FAIL %s after %d\n", e.what(), checks);
    return 1;
  }
}
