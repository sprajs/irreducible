// Original independent characteristic-function inversion. This integrates
// exp(lambda*(exp(it)-1)-sigma^2*t^2/2), not the native Poisson mixture sum.
// GL8 refinement128/256 consumes<=5% of the unchanged log-density allocation.
#include "irred/detector_selection.hpp"
#include <cmath>
#include <iostream>
#include <numbers>
#include <stdexcept>
using W = long double;
using namespace irred::detector;
namespace {
unsigned checks = 0;
void need(bool v, const char *why) {
  ++checks;
  if (!v)
    throw std::runtime_error(why);
}
void near(W a, W b, W allocation) {
  need(std::abs(a - b) <= allocation, "independent frozen allocation");
}
template <class F> W gl(F f, W end, unsigned panels) {
  const W nodes[]{
      .183434642495649804939476142360184L, .525532409916328985817739049189246L,
      .796666477413626739591553936475831L, .960289856497536231683560868569473L};
  const W weights[]{
      .362683783378361982965150449277195L, .313706645877887287337962201986601L,
      .222381034453374470544355994426240L, .101228536290376259152531354309962L};
  W total = 0, h = end / (2 * panels);
  for (unsigned n = 0; n < panels; ++n)
    for (unsigned j = 0; j < 4; ++j) {
      const W mid = (2 * n + 1) * h;
      total += h * weights[j] * (f(mid - h * nodes[j]) + f(mid + h * nodes[j]));
    }
  return total;
}
struct Reference {
  W below, density;
};
Reference invert(W lambda, W sd, W gain, W threshold, W observed,
                 unsigned panels) {
  const W end = 24 / sd, pi = std::numbers::pi_v<W>;
  auto amplitude = [&](W t) {
    return std::exp(lambda * (std::cos(t) - 1) - sd * sd * t * t / 2);
  };
  const W below =
      .5L - gl(
                [&](W t) {
                  return amplitude(t) *
                         std::sin(lambda * std::sin(t) - threshold * t) / t;
                },
                end, panels) /
                pi;
  const W density =
      gain *
      gl(
          [&](W t) {
            return amplitude(t) * std::cos(lambda * std::sin(t) - observed * t);
          },
          end, panels) /
      pi;
  // exp(-288) tail bound is negligible for these sigma>=.75 controls;
  // floating/quadrature ancestry/refinement remains empirical.
  return {below, density};
}
void cases() {
  for (W lambda : {0.L, 2.L, 15.L, 64.L})
    for (W sd : {.75L, 2.L}) {
      Input model{PhotonLaw::poisson_arrivals,
                  double(lambda),
                  1,
                  0,
                  0,
                  1,
                  double(sd),
                  2,
                  3};
      const W threshold = lambda + .25L, observed = lambda + 1.L;
      const W threshold_adu = threshold / 2 + 3,
              observed_adu = observed / 2 + 3;
      const Observation rows[] = {{double(threshold_adu),
                                   false,
                                   {},
                                   {},
                                   SelectionMeasure::joint_detection_record},
                                  {double(threshold_adu),
                                   true,
                                   {},
                                   double(observed_adu),
                                   SelectionMeasure::joint_detection_record},
                                  {double(threshold_adu),
                                   true,
                                   {},
                                   double(observed_adu),
                                   SelectionMeasure::selected_only}};
      auto coarse = invert(lambda, sd, 2, threshold, observed, 128),
           fine = invert(lambda, sd, 2, threshold, observed, 256);
      need(fine.below > 0 && fine.below < 1 && fine.density > 0,
           "positive independent probability/density");
      const W refs[]{std::log(fine.below), std::log(fine.density),
                     std::log(fine.density / (1 - fine.below))};
      const W olds[]{std::log(coarse.below), std::log(coarse.density),
                     std::log(coarse.density / (1 - coarse.below))};
      const auto result =
          likelihood(model, rows, {3, 1024 * 1024, 256, 5e-13, 2e-11});
      need(result.status == irred::numerics::Status::ok &&
               result.rows.size() == 3,
           "native likelihood batch");
      for (unsigned j = 0; j < 3; ++j) {
        const W budget = 5e-13L + 2e-11L * (1 + std::abs(refs[j]));
        near(refs[j], olds[j], .05L * budget);
        need(result.rows[j].status == irred::numerics::Status::ok &&
                 result.rows[j].log_value,
             "native finite density");
        near(*result.rows[j].log_value, refs[j], budget);
        near(*result.rows[j].detection_probability, 1 - fine.below,
             5e-13L + 2e-11L * (1 - fine.below));
      }
    }
}
} // namespace
int main() {
  try {
    cases();
    std::cout << "PASS " << checks
              << " independent Fourier detector controls\n";
  } catch (const std::exception &e) {
    std::cerr << "FAIL " << checks << ": " << e.what() << '\n';
    return 1;
  }
}
