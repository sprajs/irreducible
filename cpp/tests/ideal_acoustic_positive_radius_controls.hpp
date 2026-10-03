#pragma once
#include "../src/ideal_acoustic_positive_radius.hpp"
#include <cstdlib>
#include <iostream>

namespace ideal_acoustic_positive_test {
namespace detail = irred::cosmology::detail;
using W = long double;
using S = irred::numerics::Status;
using Vector = detail::IdealAcousticRadiusVector;
using Diagonal = detail::IdealAcousticRadiusDiagonal;

inline void require(bool good, const char *message) {
  if (!good) {
    std::cerr << "FAIL ideal acoustic positive radius: " << message << '\n';
    std::exit(1);
  }
}
inline void close(W actual, W expected, const char *message) {
  const W allowance = 4096 * std::numeric_limits<W>::epsilon() *
                      (1 + std::abs(actual) + std::abs(expected));
  require(std::isfinite(actual) && std::isfinite(expected) &&
              std::abs(actual - expected) <= allowance,
          message);
}
struct Counter {
  std::size_t writes = 0, diagnostics = 0;
  std::size_t maximum_writes = 100000, maximum_diagnostics = 100000;
  static bool write(void *opaque, std::size_t count) noexcept {
    auto &self = *static_cast<Counter *>(opaque);
    if (self.writes > self.maximum_writes ||
        count > self.maximum_writes - self.writes)
      return false;
    self.writes += count;
    return true;
  }
  static bool diagnostic(void *opaque) noexcept {
    auto &self = *static_cast<Counter *>(opaque);
    if (self.diagnostics >= self.maximum_diagnostics)
      return false;
    ++self.diagnostics;
    return true;
  }
  detail::IdealAcousticTransportAccounting accounting() noexcept {
    return {this, write, diagnostic};
  }
};
inline Diagonal synthetic_diagonal(W first, W second = 0) {
  Diagonal diagonal{};
  diagonal.status = S::ok;
  diagonal.mu_lower[0] = diagonal.mu_upper[0] = first;
  diagonal.mu_lower[1] = diagonal.mu_upper[1] = second;
  return diagonal;
}
inline Vector scalar_step(W base, W h, W mu_start, W mu_endpoint,
                          W p_start, W p_endpoint, Counter &counter) {
  Vector initial{}, start{}, endpoint{}, result{};
  initial[0] = base;
  start[0] = p_start;
  endpoint[0] = p_endpoint;
  require(detail::ideal_acoustic_radius_corrector(
              initial, start, endpoint, h, synthetic_diagonal(mu_start),
              synthetic_diagonal(mu_endpoint), result, counter.accounting()) ==
              S::ok,
          "scalar PC corrector admitted");
  return result;
}
inline void diagonal_controls() {
  Counter counter;
  Diagonal diagonal;
  const detail::IdealAcousticCoefficients actual{.01L, -1, 1, 1, .25L, .25L, 0};
  require(detail::ideal_acoustic_radius_diagonal(actual, diagonal,
                                                counter.accounting()) == S::ok,
          "actual stored coefficient diagonal admitted");
  require(counter.writes == 18 && counter.diagnostics == 1,
          "diagonal has eighteen owned scalar slots and one graph diagnostic");
  require(diagonal.mu_lower[2] <= 1.25L && diagonal.mu_upper[2] >= 1.25L &&
              diagonal.mu_lower[3] <= 2 && diagonal.mu_upper[3] >= 2 &&
              diagonal.positive[2] == 0 && diagonal.positive[3] == 0 &&
              diagonal.mu_lower[4] == 1 && diagonal.mu_upper[4] == 1,
          "negative stored diagonal enclosed and exact phi damping retained");
  for (std::size_t i : {std::size_t(0), std::size_t(1), std::size_t(5)})
    require(diagonal.mu_lower[i] == 0 && diagonal.mu_upper[i] == 0 &&
                diagonal.positive[i] == 0,
            "structural diagonal zeros remain exact");
  auto growing = actual;
  growing.L = 2;
  require(detail::ideal_acoustic_radius_diagonal(growing, diagonal,
                                                counter.accounting()) == S::ok &&
              diagonal.mu_lower[2] == 0 && diagonal.mu_upper[2] == 0 &&
              diagonal.positive[2] >= 1.75L && diagonal.positive[3] >= 1,
          "positive actual diagonal remains in the positive force");
  auto uncertain_zero = actual;
  uncertain_zero.L = uncertain_zero.loading_over_one_plus_loading;
  require(detail::ideal_acoustic_radius_diagonal(uncertain_zero, diagonal,
                                                counter.accounting()) == S::ok &&
              diagonal.mu_lower[2] == 0 && diagonal.mu_upper[2] > 0 &&
              diagonal.positive[2] > 0,
          "uncertain zero damping lower is a valid coefficient bound");
}
inline void decay_and_nonautonomous_controls() {
  Counter counter;
  const W h = .25L;
  const W decay = scalar_step(1, h, 1, 1, 0, 0, counter)[0];
  close(decay, (1 - h / 2) / (1 + h / 2), "pure decay rational PC factor");
  require(decay > 0 && decay < std::exp(-h),
          "positive PC decay is not a continuous upper envelope");
  const auto forcing = [](W t) {
    return 2 + 2 * t + 4 * t * t + t * t * t + t * t * t * t;
  };
  const W step = .125L;
  const W exact = 1 + step + step * step * step;
  const W defect = step * step * step / (2 * (1 + step / 2 + step * step / 2));
  close(scalar_step(1, step, 1, 1 + step, forcing(0), forcing(step), counter)[0],
        exact + defect, "nonautonomous scalar exact local order-three defect");
  W errors[2]{};
  for (unsigned level = 0; level < 2; ++level) {
    const unsigned count = level == 0 ? 8 : 16;
    const W dt = 1.L / count;
    W value = 1;
    for (unsigned j = 0; j < count; ++j) {
      const W t = W(j) / count;
      value = scalar_step(value, dt, 1 + t, 1 + t + dt, forcing(t),
                          forcing(t + dt), counter)[0];
    }
    errors[level] = std::abs(value - 3);
  }
  require(errors[0] > 0 && errors[1] > 0 && errors[1] < errors[0] / 3,
          "nonautonomous radius mesh refinement retains conditional order two");
}
inline Vector coupled_step(const Vector &base, W h, bool forced,
                           Counter &counter) {
  Vector start{}, predictor{}, endpoint{}, result{};
  start[0] = base[1] + (forced ? 1 : 0);
  start[1] = 2 * base[0] + (forced ? 2 : 0);
  const auto diagonal = synthetic_diagonal(2, 1);
  require(detail::ideal_acoustic_radius_predictor(
              base, start, h, diagonal, predictor, counter.accounting()) == S::ok,
          "coupled complete predictor admitted");
  endpoint[0] = predictor[1] + (forced ? 1 : 0);
  endpoint[1] = 2 * predictor[0] + (forced ? 2 : 0);
  require(detail::ideal_acoustic_radius_corrector(
              base, start, endpoint, h, diagonal, diagonal, result,
              counter.accounting()) == S::ok,
          "coupled complete endpoint predictor force admitted");
  return result;
}
inline void coupled_controls() {
  Counter counter;
  Vector equilibrium{};
  equilibrium[0] = 1;
  equilibrium[1] = 2;
  const W h = .125L;
  const auto unchanged = coupled_step(equilibrium, h, false, counter);
  close(unchanged[0], 1, "coupled neutral equilibrium first component");
  close(unchanged[1], 2, "coupled neutral equilibrium second component");
  const auto forced = coupled_step(equilibrium, h, true, counter);
  close(forced[0], 1 + h - h * h * h / ((1 + h) * (1 + h)),
        "coupled forced first exact PC defect");
  close(forced[1], 2 + 2 * h - 2 * h * h * h / ((1 + 2 * h) * (1 + h / 2)),
        "coupled forced second exact PC defect");
  require(forced[0] < 1 + h && forced[1] < 2 + 2 * h,
          "positive coupling and forcing do not certify continuous domination");
  W errors[2]{};
  for (unsigned level = 0; level < 2; ++level) {
    const unsigned count = level == 0 ? 8 : 16;
    Vector value = equilibrium;
    for (unsigned j = 0; j < count; ++j)
      value = coupled_step(value, 1.L / count, true, counter);
    errors[level] = std::abs(value[0] - 2) + std::abs(value[1] - 4);
  }
  require(errors[0] > 0 && errors[1] > 0 && errors[1] < errors[0] / 3,
          "coupled radius mesh refinement retains conditional order two");
}
inline void discrete_impulse_controls() {
  Counter counter;
  const W h = .5L, delta = .125L;
  Vector zero{}, force1{}, force2{}, force3{}, force4{}, v2{}, v3{}, v4{}, result{};
  std::array<W, 6> b2{}, no_assembly{};
  b2[0] = delta;
  require(detail::ideal_acoustic_radius_stage(zero, b2, h / 2, v2,
                                             counter.accounting()) == S::ok,
          "stage-two assembly enters aggregate arithmetic pulse once");
  force2[6] = v2[6]; // Scalar positive |A|=1, distinct from signed A=-1.
  require(detail::ideal_acoustic_radius_stage(force2, no_assembly, h / 2, v3,
                                             counter.accounting()) == S::ok,
          "stage-three positive causal impulse");
  force3[6] = v3[6];
  require(detail::ideal_acoustic_radius_stage(force3, no_assembly, h, v4,
                                             counter.accounting()) == S::ok,
          "stage-four positive causal impulse");
  force4[6] = v4[6];
  const std::array<const Vector *, 4> forces{&force1, &force2, &force3, &force4};
  const std::array<W, 4> weights{h / 6, h / 3, h / 3, h / 6};
  require(detail::ideal_acoustic_radius_endpoint(zero, forces, weights,
                                                no_assembly, result,
                                                counter.accounting()) == S::ok,
          "four once-owned RK impulse forces reach final endpoint");
  const W signed_v3 = -h * delta / 2;
  const W signed_v4 = -h * signed_v3;
  const W signed_endpoint = h * (-2 * delta - 2 * signed_v3 - signed_v4) / 6;
  close(signed_endpoint, -13 * delta / 96, "independent signed RK4 impulse gain");
  close(result[6], 21 * delta / 96, "positive local RK4 impulse gain");
  require(result[6] >= std::abs(signed_endpoint),
          "positive pulse graph covers signed discrete endpoint error");
  std::array<W, 6> bend{};
  bend[0] = delta;
  const std::array<const Vector *, 4> zero_forces{&zero, &zero, &zero, &zero};
  require(detail::ideal_acoustic_radius_endpoint(zero, zero_forces, weights,
                                                bend, result,
                                                counter.accounting()) == S::ok,
          "final assembly endpoint admitted");
  require(result[6] == delta && result[0] == 0,
          "final central assembly outside damping denominator and no source donation");
  // This primitive-only pulse is not the separate PC-map endpoint bridge.
  // Its absent SOURCE component is an adversarial witness: physical endpoint
  // Y/Z perturbations can still change P_endpoint. Parent owns that bridge's
  // actual-frame controls; no continuous/PC-map certificate follows here.
  Vector source_force{}, source_stage{};
  source_force[0] = v2[6];
  require(detail::ideal_acoustic_radius_stage(source_force, no_assembly, h / 2,
                                             source_stage,
                                             counter.accounting()) == S::ok &&
              source_stage[0] >= h * delta / 2 && source_stage[6] == 0,
          "joint stage preserves SOURCE feedback from an intermediate arithmetic pulse");
}
inline void refusal_controls() {
  Counter counter;
  Vector base{}, force{}, result{};
  std::array<W, 6> assembly{};
  const auto diagonal = synthetic_diagonal(1);
  require(detail::ideal_acoustic_radius_predictor(base, force, .125L, diagonal,
                                                 result, counter.accounting()) ==
              S::ok,
          "exact-zero predictor remains zero");
  for (W value : result)
    require(value == 0, "no manufactured radius on an exact-zero graph");
  require(detail::ideal_acoustic_radius_corrector(
              base, force, force, 2, diagonal, diagonal, result,
              counter.accounting()) == S::conditioning_budget_exceeded,
          "nonpositive actual numerator-factor guard refuses rather than clips");
  auto uncertain = synthetic_diagonal(0);
  uncertain.mu_upper[0] = 1;
  base[0] = 1;
  require(detail::ideal_acoustic_radius_corrector(
              base, force, force, .125L, uncertain, uncertain, result,
              counter.accounting()) == S::ok && result[0] >= 1,
          "uncertain zero lower damping never becomes excessive attenuation");
  for (W invalid : {std::numeric_limits<W>::quiet_NaN(),
                    std::numeric_limits<W>::infinity(), -1.L}) {
    force[0] = invalid;
    require(detail::ideal_acoustic_radius_stage(force, assembly, .125L, result,
                                               counter.accounting()) != S::ok,
            "invalid radius inputs remain refused");
  }
  force = {};
  force[0] = std::numeric_limits<W>::min();
  require(detail::ideal_acoustic_radius_stage(force, assembly, .5L, result,
                                             counter.accounting()) ==
              S::outside_domain,
          "positive produced subnormal is not silently zeroed");
  Counter limited;
  limited.maximum_writes = 3;
  force = {};
  result.fill(-1);
  require(detail::ideal_acoustic_radius_stage(force, assembly, .125L, result,
                                             limited.accounting()) == S::work_limit &&
              limited.writes == 3 && limited.diagnostics == 1 && result[2] == 0 &&
              result[3] == -1,
          "first exhausted output prefix preserved without prospective padding");
  Counter no_diagnostic;
  no_diagnostic.maximum_diagnostics = 0;
  require(detail::ideal_acoustic_radius_stage(force, assembly, .125L, result,
                                             no_diagnostic.accounting()) ==
              S::work_limit && no_diagnostic.writes == 0,
          "diagnostic refusal precedes any output assignment");
}
inline void run() {
  diagonal_controls();
  decay_and_nonautonomous_controls();
  coupled_controls();
  discrete_impulse_controls();
  refusal_controls();
}
} // namespace ideal_acoustic_positive_test
