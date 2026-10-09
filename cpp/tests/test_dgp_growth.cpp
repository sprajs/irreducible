// Independent fixed-mesh explicit midpoint/Richardson in original D,D_N,
// direct source E/beta. Production evolves g,g_N using adaptive RK4 and stable
// mu.
#include "irred/dgp_growth.hpp"
#include <algorithm>
#include <array>
#include <cfenv>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace irred::cosmology;
namespace {
using W = long double;
using State = std::array<W, 2>;
using S = irred::numerics::Status;
size_t checks = 0;
W maxref = 0, maxdifference = 0;
void need(bool yes, const char *why) {
  ++checks;
  if (!yes)
    throw std::runtime_error(why);
}
void near(W a, W b, W tolerance, const char *why) {
  need(std::abs(a - b) <= tolerance * std::abs(b), why);
}
State midpoint(W m, W a, size_t count, W initial = 1e-8L) {
  const W x = (1 - m) * std::pow(initial, 1.5L) / (2 * std::sqrt(m));
  State y{initial * (1 - 11 * x / 12), initial * (1 - 55 * x / 24)};
  const W start = std::log(initial), h = (std::log(a) - start) / count;
  auto rhs = [&](W n, State state) {
    const W aa = std::exp(n), rho = m / (aa * aa * aa), rc = (1 - m) / 2,
            s = std::sqrt(rho + rc * rc), e = s + rc, om = rho / (e * e),
            hn = -1.5L * rho / (s * e),
            beta = m == 1 ? -std::numeric_limits<W>::infinity()
                          : 1 - 2 * e / (1 - m) * (1 + hn / 3),
            mu = m == 1 ? 1 : 1 + 1 / (3 * beta);
    return State{state[1], -(2 + hn) * state[1] + 1.5L * mu * om * state[0]};
  };
  for (size_t i = 0; i < count; ++i) {
    const W n = start + i * h;
    const auto k = rhs(n, y);
    const auto mid = rhs(n + h / 2, {y[0] + h * k[0] / 2, y[1] + h * k[1] / 2});
    y = {y[0] + h * mid[0], y[1] + h * mid[1]};
  }
  return y;
}
State reference(W m, W a, size_t n, W initial = 1e-8L) {
  auto coarse = midpoint(m, a, n, initial),
       fine = midpoint(m, a, n * 2, initial);
  State y{(4 * fine[0] - coarse[0]) / 3, (4 * fine[1] - coarse[1]) / 3};
  return {y[0], y[1] / y[0]};
}
void refused(const DGPValue &v, S cause, const char *why) {
  need(v.status == cause && !v.value, why);
}
} // namespace
int main() {
  try {
    for (double m : {.05, .3, .8, 1.}) {
      auto owner = prepare_dgp_growth({m, 70});
      need(owner.status() == S::ok, "DGP admitted own model");
      const double scales[]{1e-4, .01, .5, 1};
      auto got = owner.evaluate(scales, 63);
      need(got.status == S::ok && got.rows.size() == 4, "complete DGP outputs");
      DGPPolicy refined;
      refined.relative_tolerance = 1e-10;
      auto fine = owner.evaluate(scales, dgp_d | dgp_f, refined);
      need(fine.status == S::ok, "production refined policy");
      for (size_t i = 0; i < 4; ++i) {
        const auto &r = got.rows[i];
        auto ref = reference(m, scales[i], 65536),
             old = reference(m, scales[i], 32768),
             early = reference(m, scales[i], 65536, 1e-9L);
        const W om = *r.omega_m.value, e = *r.e.value, a = scales[i];
        near(e * e - (1 - m) * e, m / (a * a * a), 3e-14,
             "modified Friedmann constraint");
        near(*r.h.value, 70 * e, 2e-15, "H units");
        near(*r.mu.value, 2 * (1 + 2 * om * om) / (3 * (1 + om * om)), 2e-15,
             "source mu stable identity");
        need(*r.mu.value >= 2. / 3 && *r.mu.value <= 1,
             "mu physical dust bounds");
        if (a == 1)
          near(e, 1, 2e-15, "flat present closure");
        for (size_t j = 0; j < 2; ++j) {
          W actual = j ? *r.f.value : *r.d.value,
            prod = j ? *fine.rows[i].f.value : *fine.rows[i].d.value;
          maxref = std::max(maxref, std::abs(ref[j] - old[j]) /
                                        (2e-7L * std::abs(ref[j])));
          near(ref[j], old[j], .05L * 2e-7,
               "independent reference mesh refinement");
          near(ref[j], early[j], .05L * 2e-7,
               "independent initial endpoint refinement");
          maxdifference = std::max(maxdifference, std::abs(actual - ref[j]) /
                                                      std::abs(ref[j]));
          near(actual, ref[j], 2e-7,
               "independent midpoint original-D reference");
          near(actual, prod, 1e-8, "production step policy refinement");
        }
        need(r.d.absolute_error_estimate <= 1e-8 * *r.d.value &&
                 r.f.absolute_error_estimate <= 1e-8 * *r.f.value,
             "returned D/f budget");
        if (m == 1)
          need(*r.d.value == a && *r.f.value == 1 && r.trials == 0 &&
                   r.callbacks == 0,
               "exact EdS infinite-rc endpoint");
      }
      const double tiny[]{1e-4};
      auto early = owner.evaluate(tiny, 63);
      const W x = (1 - W(m)) * 1e-6L / (2 * std::sqrt(W(m)));
      if (m != 1) {
        near(*early.rows[0].d.value / 1e-4, 1 - 11 * x / 12, 2e-10,
             "early growing coefficient");
        near(*early.rows[0].f.value, 1 - 11 * x / 8, 2e-10,
             "early f coefficient");
      }
    }
    auto owner = prepare_dgp_growth({.3, 70});
    const double mixed[]{.5,   -1,   std::numeric_limits<double>::quiet_NaN(),
                         1e-5, 1.01, .5};
    auto batch = owner.evaluate(mixed, 63);
    need(batch.rows.size() == 6, "mixed order retained");
    need(batch.rows[0].d.value == batch.rows[5].d.value,
         "duplicate order deterministic");
    refused(batch.rows[1].d, S::outside_domain, "negative a");
    refused(batch.rows[2].f, S::nonfinite_input, "NaN a");
    refused(batch.rows[3].e, S::outside_domain, "unsupported early domain");
    for (DGPParameters p :
         {DGPParameters{.049, 70}, {1.01, 70}, {.3, 0}, {.3, 1001}})
      need(prepare_dgp_growth(p).status() == S::outside_domain,
           "unsupported model parameters");
    need(prepare_dgp_growth({std::numeric_limits<double>::quiet_NaN(), 70})
                 .status() == S::nonfinite_input,
         "nonfinite parameter");
    const double a[]{1, .5};
    need(owner.evaluate(a, 0).status == S::invalid_input &&
             owner.evaluate(a, 64).status == S::invalid_input,
         "invalid masks");
    auto only = owner.evaluate(a, dgp_e | dgp_mu);
    need(only.status == S::ok && only.callbacks == 0 && !only.rows[0].d.value &&
             !only.rows[0].f.value,
         "background avoids ODE");
    DGPPolicy p;
    p.maximum_points = 1;
    need(owner.evaluate(a, 63, p).status == S::work_limit,
         "point admission before allocation");
    p = {};
    p.maximum_native_bytes = 1;
    need(owner.evaluate(a, 63, p).rows.empty(), "byte quota before allocation");
    p = {};
    p.maximum_steps_per_point = 1;
    auto exhausted = owner.evaluate(a, 63, p);
    need(exhausted.status == S::work_limit && exhausted.rows[0].e.value &&
             !exhausted.rows[0].d.value && exhausted.trials == 2 &&
             exhausted.callbacks == 24,
         "growth work failure retains background and exact counters");
    p = {};
    p.maximum_total_steps = 1;
    auto total = owner.evaluate(a, dgp_d, p);
    need(total.trials == 1 && total.status == S::work_limit,
         "shared total work budget");
    p = {};
    p.relative_tolerance = 1e-25;
    auto tight = owner.evaluate(a, dgp_d | dgp_f, p);
    refused(tight.rows[0].d, S::conditioning_budget_exceeded,
            "unattainable tolerance preserved");
    p.relative_tolerance = std::numeric_limits<double>::infinity();
    need(owner.evaluate(a, 63, p).status == S::invalid_input,
         "nonfinite policy");
    need(!dgp_payload_bound(std::numeric_limits<size_t>::max()),
         "payload arithmetic overflow");
    auto copied = owner;
    auto moved = std::move(copied);
    need(copied.status() == S::invalid_input && !copied.parameters() &&
             moved.status() == S::ok,
         "move invalidates owner");
    auto *self = &moved;
    moved = std::move(*self);
    need(moved.status() == S::ok, "self move");
    const int rounding = std::fegetround();
    std::fesetround(FE_DOWNWARD);
    need(owner.evaluate(a, 63).status == S::invalid_input, "rounding refusal");
    std::fesetround(rounding);
    std::cout << "DGP checks=" << checks
              << " max_reference_budget_fraction=" << (double)maxref
              << " max_relative_difference=" << (double)maxdifference << "\n";
    return 0;
  } catch (const std::exception &e) {
    std::cerr << "DGP failure after " << checks << " checks: " << e.what()
              << "\n";
    return 1;
  }
}
