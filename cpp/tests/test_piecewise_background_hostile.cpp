// Independent split composite Simpson integrals, rather than production
// segment antiderivatives. Samples do not certify arbitrary q(z) or priors.
#include "irred/background.hpp"
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
EvaluationPolicy policy() { return {{}, 64, 0, 2000, 1000000}; }
Request full(double z, double obs, double h = 70,
             Convention c = Convention::geometric_same_redshift) {
  return {z, 63, Observer{obs, c}, PhysicalScale(h)};
}
Slot one(const Expansion &p, Request r, EvaluationPolicy cap = policy()) {
  auto b = p.evaluate(std::span(&r, 1), cap);
  check(b.status == Status::ok && b.slots.size() == 1,
        "admitted current batch");
  return b.slots[0];
}
void geometry(const Slot &s) {
  check(s.radial.value && s.expansion.value && s.clock.value &&
            s.physical.value && s.luminosity_shape.value && s.kinematics.value,
        "independent requested groups available");
}
} // namespace
int main() {
  try {
    for (auto q : {std::array<double, 5>{0, 0, 0, 0, 0},
                   std::array<double, 5>{-1, -1, -1, -1, -1},
                   std::array<double, 5>{-.4, -.4, -.2, .1, .3},
                   std::array<double, 5>{-3, 2, -1e-10, 1e-10, -1 + 1e-10}}) {
      auto p = prepare(FixedFiveBinQ(q), FlatFLRW{});
      check(p.status() == Status::ok, "valid provider");
      for (double z : {0., .03, .1, .3, .6, 1., 1.7, 2.5}) {
        auto query = full(z, z);
        auto b = p.evaluate(std::span(&query, 1), policy());
        auto &s = b.slots[0];
        geometry(s);
        auto r = integral(z, q, 4096, false),
             fine = integral(z, q, 8192, false);
        check(std::abs(r - fine) <= .1L * (2e-15L + 2e-10L * std::abs(fine)),
              "reference refinement");
        check(std::abs(s.radial.value->integral - fine) <
                  2e-15L + 2e-10L * std::abs(fine),
              "independent radial");
        check(std::abs(s.expansion.value->expansion_E - expansion(z, q)) <
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
        check(std::abs(*s.clock.value->lookback_seconds.value -
                       expected_seconds) <
                  (2e-15L + 2e-10L * std::abs(clock_fine)) * mpc_metres /
                      70000.L,
              "independent lookback units");
        radial_ratio =
            std::max(radial_ratio, std::abs(s.radial.value->integral - fine) /
                                       (2e-15L + 2e-10L * std::abs(fine)));
        expansion_ratio =
            std::max(expansion_ratio, std::abs(s.expansion.value->expansion_E -
                                               expansion(z, q)) /
                                          (2e-14L + 2e-12L * expansion(z, q)));
        clock_ratio = std::max(clock_ratio,
                               std::abs(*s.clock.value->lookback_seconds.value -
                                        expected_seconds) /
                                   ((2e-15L + 2e-10L * std::abs(clock_fine)) *
                                    mpc_metres / 70000.L));
        refinement_ratio = std::max(
            {refinement_ratio,
             std::abs(r - fine) / (.1L * (2e-15L + 2e-10L * std::abs(fine))),
             std::abs(clock_coarse - clock_fine) /
                 (.1L * (2e-15L + 2e-10L * std::abs(clock_fine)))});
        auto k = s.kinematics.value->bin;
        bool jump =
            z > 0 && z < 2.5 && z == piecewise_q_edges[k] && q[k] != q[k - 1];
        check(s.kinematics.value->jerk.has_value() != jump,
              "jerk availability follows actual jump");
        if (s.kinematics.value->jerk)
          check(std::abs(*s.kinematics.value->jerk - q[k] * (2 * q[k] + 1)) <
                    2e-14,
                "jerk ordinary value");
        if (z == 0)
          check(s.physical.value->radial_mpc == 0 &&
                    s.physical.value->luminosity_mpc == 0 &&
                    s.kinematics.value->q0_within_piecewise_model == q[0],
                "zero geometry and one-sided q0");
      }
    }
    const std::array<double, 5> source_q{-1, 0, 1, -.5, .2};
    auto p = prepare(FixedFiveBinQ(source_q), FlatFLRW{});
    for (int i = 1; i < 5; ++i)
      for (double z :
           {std::nextafter(piecewise_q_edges[i], 0.), piecewise_q_edges[i],
            std::nextafter(piecewise_q_edges[i], 3.)}) {
        auto s = one(p, full(z, z));
        geometry(s);
        check(s.kinematics.value->q == source_q[s.kinematics.value->bin],
              "retained side q");
      }
    std::array<Request, 2> duplicate{full(2, 2), full(2, 2)};
    auto cap = policy();
    cap.maximum_segment_visits = 5;
    auto reused = p.evaluate(duplicate, cap);
    check(reused.slots[0].radial.value && reused.slots[1].radial.value &&
              reused.nodes.size() == 1 && reused.work.segment_visits == 5,
          "exact duplicate reuses one segment budget");
    std::array<Request, 2> distinct{full(2, 2), full(1.9, 1.9)};
    auto limited = p.evaluate(distinct, cap);
    check(limited.slots[0].radial.value &&
              limited.slots[1].radial.status == Status::work_limit &&
              limited.work.segment_visits == 5,
          "global distinct-node segment cap");
    cap.maximum_queries = 0;
    check(p.evaluate(duplicate, cap).slots.empty(),
          "count cap before allocation");
    auto a = one(p, full(.4, .4)), scaled = one(p, full(.4, .4, 100));
    geometry(a);
    geometry(scaled);
    check(std::abs(scaled.physical.value->radial_mpc /
                       a.physical.value->radial_mpc -
                   .7) < 2e-15,
          "H0 inverse distance scaling");
    check(std::abs(*scaled.expansion.value->h_km_s_mpc.value /
                       *a.expansion.value->h_km_s_mpc.value -
                   100. / 70) < 2e-15,
          "H0 expansion scaling");
    for (double z : {1e-200, std::numeric_limits<double>::denorm_min()}) {
      auto s = one(p, full(z, z));
      check(!s.physical.value && s.kinematics.value && s.kinematics.value->jerk,
            "failed tiny physical group leaves kinematics available");
    }
    check(prepare(FixedFiveBinQ({0, 0, 3, 0, 0}), FlatFLRW{}).status() ==
              Status::unsupported_domain,
          "q domain");
    auto badh = one(p, full(.2, .2, 0));
    check(!badh.physical.value && badh.expansion.value &&
              !badh.expansion.value->h_km_s_mpc.value,
          "H0 dimensional failure isolated from E");
    for (double d : {std::numeric_limits<double>::denorm_min(),
                     -std::numeric_limits<double>::denorm_min(), 0., -0.}) {
      auto tiny = prepare(FixedFiveBinQ({d, 0, 0, 0, 0}), FlatFLRW{});
      auto s = one(tiny, full(.05, .05));
      geometry(s);
      check(s.kinematics.value->q == d && s.kinematics.value->jerk,
            "subnormal dimensionless q retained");
      check(*s.kinematics.value->jerk == d * (2 * d + 1),
            "subnormal jerk retained");
    }
    for (auto r :
         {full(.2, .3), full(.2, -1, 70, Convention::released_zhd_zhel),
          full(.2, .2, 70, (Convention)99)}) {
      auto s = one(p, r);
      check(!s.physical.value && !s.luminosity_shape.value &&
                s.kinematics.value && s.radial.value,
            "invalid observer isolated from expansion/redshift groups");
    }
    for (double z : {3., std::numeric_limits<double>::quiet_NaN()}) {
      auto s = one(p, full(z, .2));
      check(!s.radial.value && !s.kinematics.value && !s.physical.value,
            "invalid expansion source has no payload");
    }
    auto released = one(p, full(.4, .45, 70, Convention::released_zhd_zhel));
    geometry(released);
    check(released.source.observer->convention == Convention::released_zhd_zhel,
          "observer convention retained");
    check(std::abs(released.physical.value->luminosity_mpc /
                       a.physical.value->luminosity_mpc -
                   1.45 / 1.4) < 2e-15 &&
              released.physical.value->angular_diameter_mpc ==
                  a.physical.value->angular_diameter_mpc,
          "observer prefactor distinct from DA");
    geometry(one(p, full(1e-90, 1e-90)));
    Slot unset(Request{0, 1});
    check(unset.kinematics.availability == Availability::not_requested &&
              !unset.kinematics.value,
          "unset derivatives unassessed");
    auto empty = p.evaluate({}, policy());
    check(empty.status == Status::ok && empty.slots.empty() &&
              empty.work.segment_visits == 0,
          "empty completed zero work");
    auto copied = p;
    auto c = one(copied, full(.4, .4));
    check(c.radial.value->integral == a.radial.value->integral &&
              std::get<FixedFiveBinQ>(copied.specification()).q == source_q,
          "copied immutable source and value");
    // E/kinematics alone uses no radial or clock segments, at the actual
    // jump.
    Request edge{.3, 48};
    EvaluationPolicy zero{{}, 4, 0, 0, 100000};
    auto only = p.evaluate(std::span(&edge, 1), zero);
    check(only.slots[0].expansion.value && only.slots[0].kinematics.value &&
              !only.slots[0].kinematics.value->jerk &&
              only.work.segment_visits == 0,
          "jump E and right-limit q zero segment work");
    Request scale_failure{.3, 32, {}, PhysicalScale(0)};
    auto isolated = p.evaluate(std::span(&scale_failure, 1), zero);
    check(isolated.slots[0].expansion.value &&
              isolated.slots[0].expansion.value->h_km_s_mpc.availability ==
                  Availability::failed &&
              isolated.work.segment_visits == 0,
          "nested H failure doesn't erase E");
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
