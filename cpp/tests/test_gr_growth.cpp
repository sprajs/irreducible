// Original GR ODE RK4 references, separate from production Heath quadrature.
// D_x=V,V_x=-(2-1.5Om(a))V+1.5Om(a)D. Integrate g=D/a,w=V/a
// from a0<=min(1e-7,a/1000), q0<=1e-15 with the growing initial series
// g=1-2q0/11,w=1-8q0/11. 8192/16384 step and a0 refinement occupy
// <=5% of the frozen1e-8 relative named-case comparison allocation.
#include "irred/gr_growth.hpp"
#include <array>
#include <cfenv>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace irred::cosmology;
namespace {
using W = long double;
unsigned checks = 0;
W maximum_refinement_fraction = 0, maximum_difference = 0;
void need(bool b, const char *s) {
  ++checks;
  if (!b)
    throw std::runtime_error(s);
}
void near(W got, W ref, W relative, const char *s) {
  need(std::abs(got - ref) <= relative * std::abs(ref), s);
}
std::array<W, 2> ode(double matter, double a, unsigned n, W shift = 1) {
  W a0 = std::min(1e-7L, W(a) / 1000) * shift, x0 = std::log(a0),
    h = (std::log(W(a)) - x0) / n, q0 = (1 - W(matter)) / matter * a0 * a0 * a0;
  std::array<W, 2> y{1 - 2 * q0 / 11, 1 - 8 * q0 / 11};
  auto rhs = [&](W x, std::array<W, 2> g) {
    W scale = std::exp(x), q = (1 - W(matter)) / matter * scale * scale * scale,
      om = 1 / (1 + q);
    return std::array<W, 2>{g[1] - g[0],
                            -(3 - 1.5L * om) * g[1] + 1.5L * om * g[0]};
  };
  auto add = [](std::array<W, 2> y, std::array<W, 2> k, W step) {
    return std::array<W, 2>{y[0] + step * k[0], y[1] + step * k[1]};
  };
  for (unsigned i = 0; i < n; ++i) {
    W x = x0 + i * h;
    auto k1 = rhs(x, y), k2 = rhs(x + h / 2, add(y, k1, h / 2)),
         k3 = rhs(x + h / 2, add(y, k2, h / 2)), k4 = rhs(x + h, add(y, k3, h));
    for (unsigned j = 0; j < 2; ++j)
      y[j] += h * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]) / 6;
  }
  return {W(a) * y[0], y[1] / y[0]};
}
} // namespace
int main() {
  try {
    for (double m : {1., .9, .3, .01, 1e-6}) {
      auto background = prepare(LCDM(m), FlatFLRW{});
      auto owner = prepare_gr_growth(background);
      need(owner.status() == irred::numerics::Status::ok &&
               owner.background() &&
               owner.background()->model_id() == background.model_id(),
           "owned compatible LCDM identity");
      const double a[]{1e-8, .01, .2, .7, 1};
      auto result = owner.evaluate(a, 3);
      need(result.status == irred::numerics::Status::ok &&
               result.rows.size() == 5,
           "admitted ordered growth batch");
      for (unsigned i = 0; i < 5; ++i) {
        auto ref = ode(m, a[i], 16384), old = ode(m, a[i], 8192),
             early = ode(m, a[i], 16384, .1L);
        need(result.rows[i].d.value && result.rows[i].f.value,
             "D and f accepted");
        for (unsigned j = 0; j < 2; ++j) {
          W fraction = std::abs(ref[j] - old[j]) / (1e-8L * std::abs(ref[j]));
          maximum_refinement_fraction =
              std::max(maximum_refinement_fraction, fraction);
          near(ref[j], old[j], .05L * 1e-8, "ODE refinement gate");
          near(ref[j], early[j], .05L * 1e-8,
               "growing initial endpoint refinement");
          W got = j ? *result.rows[i].f.value : *result.rows[i].d.value;
          maximum_difference = std::max(
              maximum_difference, std::abs(got - ref[j]) / std::abs(ref[j]));
          if (std::abs(got - ref[j]) > 1e-8L * std::abs(ref[j]))
            std::cerr << "m=" << m << " a=" << a[i] << " group=" << j
                      << " got=" << (double)got << " ref=" << (double)ref[j]
                      << "\n";
          near(got, ref[j], 1e-8, "independent GR ODE");
        }
        if (m == 1) {
          need(*result.rows[i].d.value == a[i] &&
                   *result.rows[i].f.value == 1 &&
                   result.rows[i].callbacks == 0,
               "EdS exact growing mode");
        }
        need(result.rows[i].d.absolute_error_estimate <=
                     1e-10 * *result.rows[i].d.value &&
                 result.rows[i].f.absolute_error_estimate <=
                     1e-10 * *result.rows[i].f.value,
             "diagnostic within requested budget");
      }
    }
    auto owner = prepare_gr_growth(prepare(LCDM(.3), FlatFLRW{}));
    const double a[]{.5,   -1,   std::numeric_limits<double>::quiet_NaN(),
                     1e-9, 1.01, .5};
    auto mixed = owner.evaluate(a, 3);
    need(mixed.rows.size() == 6 && mixed.rows[0].d.value &&
             !mixed.rows[1].d.value && !mixed.rows[2].d.value &&
             !mixed.rows[3].d.value && !mixed.rows[4].d.value &&
             mixed.rows[5].d.value,
         "mixed domain failures remain independent");
    for (unsigned mask : {1u, 2u}) {
      auto single = owner.evaluate(std::span(a, 1), mask);
      need(bool(single.rows[0].d.value) == bool(mask & 1) &&
               bool(single.rows[0].f.value) == bool(mask & 2),
           "unrequested output omitted");
    }
    GrowthPolicy p;
    p.maximum_total_callbacks = 0;
    need(!owner.evaluate(std::span(a, 1), 3, p).rows[0].d.value,
         "global callback refusal");
    p = {};
    p.maximum_callbacks_per_point = 4;
    auto limited = owner.evaluate(std::span(a, 1), 3, p);
    need(!limited.rows[0].d.value && limited.callbacks <= 4,
         "point callback refusal");
    p = {};
    p.maximum_points = 0;
    need(owner.evaluate(std::span(a, 1), 3, p).rows.empty(), "point admission");
    p = {};
    p.maximum_native_bytes = 1;
    need(owner.evaluate(std::span(a, 1), 3, p).rows.empty(),
         "payload admission");
    p = {};
    p.relative_tolerance = 1e-30;
    need(!owner.evaluate(std::span(a, 1), 3, p).rows[0].d.value,
         "unattainable budget refusal");
    need(!growth_payload_bound(SIZE_MAX), "payload overflow");
    for (auto spec :
         {ExpansionSpec(CPL(.3, -1, 0)), ExpansionSpec(ConstantQ(-.5)),
          ExpansionSpec(LCDM(0)), ExpansionSpec(LCDM(1e-7))})
      need(prepare_gr_growth(prepare(spec, FlatFLRW{})).status() !=
               irred::numerics::Status::ok,
           "unsupported closure or growing-normalization domain rejected");
    auto copy = owner;
    auto moved = std::move(owner);
    auto *self = &moved;
    moved = std::move(*self);
    need(owner.evaluate(std::span(a, 1), 3).rows.empty() &&
             moved.evaluate(std::span(a, 1), 3).rows[0].d.value &&
             copy.evaluate(std::span(a, 1), 3).rows[0].d.value,
         "owner copy move lifetime");
    auto rounding = std::fegetround();
    need(std::fesetround(FE_DOWNWARD) == 0, "set rounding");
    auto bad = moved.evaluate(std::span(a, 1), 3);
    need(std::fesetround(rounding) == 0, "restore rounding");
    need(bad.rows.empty(), "unsupported arithmetic refuses processing");
    std::cout << "PASS " << checks
              << " GR growth controls; maximum relative difference="
              << double(maximum_difference) << " reference allocation fraction="
              << double(maximum_refinement_fraction) << '\n';
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
