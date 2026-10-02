#include "irred/photometry_calibration.hpp"
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <numbers>
#include <stdexcept>
using namespace irred::photometry;
namespace {
unsigned checks = 0;
void check(bool x, const char *m) {
  ++checks;
  if (!x)
    throw std::runtime_error(m);
}
void near(double x, long double y, long double scale, bool covariance = false) {
  check(y == 0 || x != 0, "nonzero reference");
  check(std::isfinite(x), "finite");
  check(std::abs(x - y) <= 1e-280L + (covariance ? 2e-8L : 2e-10L) * scale,
        "frozen allocation");
}
} // namespace
int main() {
  try {
    std::array<double, 2> wave{1, 2}, lum{1, 1}, low{0.25, 0.25},
        high{0.75, 0.75};
    std::array<CalibrationBand, 2> bands{
        {{"a", wave, 1, 1}, {"b", wave, 2, 3}}};
    std::array<std::span<const double>, 2> t0{low, high}, t1{high, low};
    std::array<CalibrationState, 2> states{{{"s0", 1, t0}, {"s1", 3, t1}}};
    CalibrationInput in{{wave, lum},
                        1,
                        0,
                        bands,
                        states,
                        "synthetic-source",
                        "synthetic-optics",
                        "relative-discrete-law",
                        "shared-state"};
    auto run = [&]() { return evaluate_calibration(in); };
    const long double k = 1 / (4 * std::numbers::pi_v<long double>),
                      ph = 1.5L / (6.62607015e-34L * 299792458);
    auto r = run();
    check(r.status == irred::numerics::Status::ok && r.moments.has_value(),
          "law available");
    check(r.axes.size() == 4 && r.attempts.size() == 4, "ordering dimensions");
    std::array<long double, 4> a{k, ph * k, 6 * k, 6 * ph * k},
        mu{0.625L, 0.625L, 0.375L, 0.375L};
    std::array<int, 4> sign{1, 1, -1, -1};
    for (unsigned i = 0; i < 4; ++i) {
      near(r.moments->mean[i], a[i] * mu[i], a[i] * mu[i]);
      for (unsigned j = 0; j < 4; ++j)
        near(r.moments->covariance[i * 4 + j],
             sign[i] * sign[j] * 0.046875L * a[i] * a[j],
             0.046875L * a[i] * a[j], true);
    }
    auto original = *r.moments;
    states[0].relative_mass = 2;
    states[1].relative_mass = 6;
    r = run();
    check(r.status == irred::numerics::Status::ok, "mass scaling");
    for (unsigned i = 0; i < 4; ++i)
      near(r.moments->mean[i], original.mean[i], original.mean[i]);
    std::swap(states[0], states[1]);
    r = run();
    check(r.status == irred::numerics::Status::ok, "state permutation");
    for (unsigned i = 0; i < 16; ++i)
      near(r.moments->covariance[i], original.covariance[i],
           std::abs(original.covariance[i]), true);
    std::swap(states[0], states[1]);
    in.luminosity_distance_metre = 2;
    r = run();
    for (unsigned i = 0; i < 4; ++i)
      near(r.moments->mean[i], original.mean[i] / 4, original.mean[i] / 4);
    for (unsigned i = 0; i < 16; ++i)
      near(r.moments->covariance[i], original.covariance[i] / 16,
           std::abs(original.covariance[i]) / 16, true);
    in.luminosity_distance_metre = 1;
    t1 = t0;
    r = run();
    check(r.status == irred::numerics::Status::ok, "exact constant inputs");
    for (double x : r.moments->covariance)
      check(x == 0, "constant covariance zero");
    in.states = std::span(states).first(1);
    r = run();
    check(r.status == irred::numerics::Status::ok, "one state");
    in.states = states;
    t1 = {high, low};
    bands[0].collecting_area_square_metre = 0;
    r = run();
    check(r.status == irred::numerics::Status::ok, "constant zero area");
    check(r.moments->mean[0] == 0 && r.moments->covariance[0] == 0,
          "zero axis");
    bands[0].collecting_area_square_metre = 1;
    lum = {0, 0};
    r = run();
    check(r.status == irred::numerics::Status::ok, "zero source");
    for (double x : r.moments->covariance)
      check(x == 0, "zero source covariance");
    lum = {1, 1};
    low = {0.5, 0.5};
    high = {std::nextafter(0.5, 1.0), std::nextafter(0.5, 1.0)};
    r = run();
    check(r.status == irred::numerics::Status::conditioning_budget_exceeded &&
              !r.moments,
          "unresolved spread refused");
    auto state_policy = CalibrationPolicy{};
    state_policy.aggregation = CalibrationAggregation::state_resolved_only;
    auto state_law = evaluate_calibration(in, state_policy);
    check(state_law.status == irred::numerics::Status::ok && !state_law.moments &&
              state_law.attempts.size() == 4 && state_law.normalized_state_mass.size() == 2,
          "explicit state law does not request unresolved population moments");
    check(state_law.normalized_state_mass[0] == .25L &&
              state_law.normalized_state_mass[1] == .75L,
          "retained normalized masses share one owner");
    low = {0.25, 0.25};
    high = {0.75, 0.75};
    high[1] = 2;
    r = run();
    check(r.status != irred::numerics::Status::ok && !r.moments &&
              r.attempts.size() == 4,
          "partial failures retained");
    state_law = evaluate_calibration(in, state_policy);
    check(state_law.status != irred::numerics::Status::ok && !state_law.moments &&
              state_law.attempts.size() == 4 && state_law.normalized_state_mass.size() == 2,
          "state-only preserves every required refusal and normalized mass");
    high[1] = 0.75;
    for (double mass : {0.0, -1.0, std::numeric_limits<double>::denorm_min(),
                        std::numeric_limits<double>::infinity()}) {
      states[0].relative_mass = mass;
      check(!run().moments, "hostile mass");
    }
    states[0].relative_mass = 2;
    states[1].id = states[0].id;
    check(!run().moments, "duplicate state");
    states[1].id = "s1";
    bands[1].id = "a";
    check(!run().moments, "duplicate band");
    bands[1].id = "b";
    in.dependence_origin = "";
    check(!run().moments, "missing dependence");
    in.dependence_origin = "shared-state";
    auto p = CalibrationPolicy{};
    p.maximum_output_bytes = 0;
    check(evaluate_calibration(in, p).status ==
              irred::numerics::Status::work_limit,
          "byte quota");
    p = {};
    p.maximum_states = 1;
    check(evaluate_calibration(in, p).status ==
              irred::numerics::Status::work_limit,
          "state quota");
    p = {};
    p.maximum_total_knots = 1;
    check(evaluate_calibration(in, p).status ==
              irred::numerics::Status::work_limit,
          "knot quota");
    p = {};
    p.requested_outputs = incident_flux;
    check(!evaluate_calibration(in, p).moments,
          "fixed flux not stochastic output");

    // Same joint state moves both bands together, giving positive cross-band
    // covariance.
    t0 = {low, low};
    t1 = {high, high};
    r = run();
    check(r.status == irred::numerics::Status::ok, "correlated bands");
    near(r.moments->covariance[2], 0.046875L * a[0] * a[2],
         0.046875L * a[0] * a[2], true);
    t0 = {low, high};
    t1 = {high, low};
    // Per-band area/exposure scaling and band permutation preserve the full
    // law.
    bands[0].collecting_area_square_metre = 2;
    bands[0].observer_exposure_second = 3;
    r = run();
    check(r.status == irred::numerics::Status::ok, "area exposure scaling");
    near(r.moments->mean[0], 6 * original.mean[0], 6 * original.mean[0]);
    near(r.moments->covariance[2], 6 * original.covariance[2],
         6 * std::abs(original.covariance[2]), true);
    bands[0].collecting_area_square_metre = 1;
    bands[0].observer_exposure_second = 1;
    std::swap(bands[0], bands[1]);
    std::swap(t0[0], t0[1]);
    std::swap(t1[0], t1[1]);
    r = run();
    check(r.status == irred::numerics::Status::ok, "band permutation");
    near(r.moments->mean[0], original.mean[2], original.mean[2]);
    std::swap(bands[0], bands[1]);
    std::swap(t0[0], t0[1]);
    std::swap(t1[0], t1[1]);
    p = {};
    p.requested_outputs = collected_energy;
    r = evaluate_calibration(in, p);
    check(r.status == irred::numerics::Status::ok && r.axes.size() == 2,
          "energy only");
    lum = {1e-160, 1e-160};
    r = evaluate_calibration(in, p);
    check(!r.moments, "positive covariance underflow refused");
    lum = {1, 1};
    in.luminosity_distance_metre = 1e200;
    r = run();
    check(!r.moments && r.attempts.size() == 4, "parent positive underflow");
    in.luminosity_distance_metre = 1;
    bands[0].collecting_area_square_metre = 1e308;
    bands[0].observer_exposure_second = 1e308;
    r = run();
    check(!r.moments && r.attempts.size() == 4, "parent overflow");
    bands[0].collecting_area_square_metre = 1;
    bands[0].observer_exposure_second = 1;
    wave = {2, 1};
    check(!run().moments, "hostile axis");
    wave = {1, 2};
    p = {};
    p.maximum_bands = 1;
    check(evaluate_calibration(in, p).status ==
              irred::numerics::Status::work_limit,
          "band quota");
    // Independent polynomial antiderivatives: L=lambda, T=lambda-1; E
    // integral=5/6, photon integral=17/12.
    lum = {1, 2};
    low = {0, 1};
    t0 = {low, low};
    in.states = std::span(states).first(1);
    r = run();
    check(r.status == irred::numerics::Status::ok, "linear controls");
    near(r.moments->mean[0], 5 * k / 6, 5 * k / 6);
    near(r.moments->mean[1], 17 * k / (12 * 6.62607015e-34L * 299792458),
         17 * k / (12 * 6.62607015e-34L * 299792458));
    std::cout << "photometry calibration checks=" << checks << " failures=0\n";
  } catch (const std::exception &e) {
    std::cerr << "checks=" << checks << " failure=" << e.what() << '\n';
    return 1;
  }
}
