// Independent fixed GL8 in x=sqrt(a/a_drag), unlike production adaptive
// Simpson in a/a_drag. Analytic Lambda=0 and radiation-only controls are
// derived directly from the declared equation, not package-generated goldens.
// Photon-baryon convention: Eisenstein & Hu astro-ph/9709112 equations5,6.
#include "irred/sound_horizon.hpp"
#include <array>
#include <bit>
#include <cfenv>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace irred;
using namespace irred::cosmology;
namespace {
unsigned checks = 0;
long double max_error = 0, max_refinement = 0;
void check(bool b, const char *w) {
  ++checks;
  if (!b)
    throw std::runtime_error(w);
}
SoundHorizonPolicy policy() {
  return {1e-9, 2e-11, 200000, 40, 32, 1000000, 1u << 20};
}
SoundHorizonRequest typical() {
  return {{70, 0.3, 9e-5, 0.05, 5e-5},
          1059,
          "explicit synthetic drag epoch, not predicted"};
}
SoundHorizonRow eval(SoundHorizonRequest x, SoundHorizonPolicy p = policy()) {
  auto b = evaluate_sound_horizon(std::span(&x, 1), p);
  check(b.status == numerics::Status::ok && b.rows.size() == 1,
        "materialized batch");
  check(b.callbacks == b.rows[0].callbacks, "batch exact callback accounting");
  return b.rows[0];
}
void finite(const SoundHorizonRow &r) {
  check(r.status == numerics::Status::ok, "finite status");
  check(r.sound_horizon_mpc && std::isfinite(*r.sound_horizon_mpc) &&
            *r.sound_horizon_mpc > 0,
        "positive finite horizon");
  check(std::isfinite(r.error_estimate_mpc) && r.error_estimate_mpc >= 0,
        "finite error diagnostic");
}
long double prefactor(const SoundHorizonRequest &r) {
  return 299792458.0L / (1000.L * r.model.h0_km_s_mpc * std::sqrt(3.L));
}
long double reference(const SoundHorizonRequest &r, unsigned panels) {
  constexpr std::array<long double, 4> nodes{
      .183434642495649804939476142360184L, .525532409916328985817739049189246L,
      .796666477413626739591553936475831L, .960289856497536231683560868569473L};
  constexpr std::array<long double, 4> weights{
      .362683783378361982965150449277195L, .313706645877887287337962201986601L,
      .222381034453374470544355994426240L, .101228536290376259152531354309962L};
  const long double ad = 1 / (1 + static_cast<long double>(r.z_drag)),
                    om = r.model.omega_m, orr = r.model.omega_r,
                    lambda = 1 - om - orr,
                    B = 3 * static_cast<long double>(r.model.omega_b) /
                        (4 * r.model.omega_gamma);
  long double sum = 0;
  for (unsigned j = 0; j < panels; ++j) {
    const long double mid = (j + 0.5L) / panels, half = 0.5L / panels;
    for (unsigned k = 0; k < 4; ++k)
      for (int sign : {-1, 1}) {
        const long double x = mid + sign * half * nodes[k], a = ad * x * x,
                          a2 = a * a;
        const long double f =
            2 * ad * x /
            (std::sqrt(orr + om * a + lambda * a2 * a2) * std::sqrt(1 + B * a));
        sum += half * weights[k] * f;
      }
  }
  return prefactor(r) * sum;
}
long double analytic_lambda_zero(const SoundHorizonRequest &r) {
  const long double ad = 1 / (1 + static_cast<long double>(r.z_drag)),
                    m = r.model.omega_m, rad = r.model.omega_r,
                    B = 3 * static_cast<long double>(r.model.omega_b) /
                        (4 * r.model.omega_gamma);
  if (B == 0)
    return prefactor(r) * 2 * ad / (std::sqrt(rad + m * ad) + std::sqrt(rad));
  const long double den = std::sqrt(m) + std::sqrt(B * rad),
                    num = std::sqrt(m) * std::sqrt(1 + B * ad) +
                          std::sqrt(B) * std::sqrt(rad + m * ad);
  return prefactor(r) * 2 / std::sqrt(m * B) * std::log1p((num - den) / den);
}
void control(SoundHorizonRequest r, bool analytic = false) {
  auto row = eval(r);
  finite(row);
  const auto a = reference(r, 64), b = reference(r, 128), c = reference(r, 256);
  const long double budget = 1e-9L + 2e-11L * std::fabs(c),
                    refdiff = std::max(std::fabs(a - b), std::fabs(b - c));
  max_refinement = std::max(max_refinement, refdiff);
  if (refdiff > 0.05L * budget)
    std::cerr << "reference z=" << r.z_drag
              << " refinements=" << (double)std::fabs(a - b) << ","
              << (double)std::fabs(b - c)
              << " allowance=" << (double)(.05L * budget) << "\n";
  check(refdiff <= 0.05L * budget, "frozen independent reference allocation");
  const auto error =
      std::fabs(static_cast<long double>(*row.sound_horizon_mpc) - c);
  max_error = std::max(max_error, error);
  check(error <= budget, "frozen horizon observable allocation");
  check(row.error_estimate_mpc <=
            policy().absolute_tolerance_mpc +
                policy().relative_tolerance * std::fabs(*row.sound_horizon_mpc),
        "reported diagnostic within requested allowance");
  if (analytic)
    check(std::fabs(c - analytic_lambda_zero(r)) <= 0.05L * budget,
          "independent analytic reference agreement");
}
void rejected(SoundHorizonRequest r, numerics::Status cause) {
  auto x = eval(std::move(r));
  check(x.status == cause && !x.sound_horizon_mpc,
        "rejected precise cause no payload");
  check(x.callbacks == 0, "model rejection zero callback");
}
} // namespace
int main() {
  try {
    control(typical());
    for (double z : {0.0, 10.0, 100.0, 2000.0}) {
      auto r = typical();
      r.z_drag = z;
      control(r);
    }
    for (double baryon : {0.0, 0.1, 0.5}) {
      SoundHorizonRequest r{
          {70, 0.75, 0.25, baryon, 0.1}, 2, "Lambda-zero analytic control"};
      control(r, true);
    }
    SoundHorizonRequest radiation{
        {70, 0, 1, 0, 1}, 1000, "radiation-only analytic control"};
    control(radiation, true);
    auto raw = eval(radiation);
    finite(raw);
    const long double exact =
        prefactor(radiation) / (1 + static_cast<long double>(radiation.z_drag));
    check(std::fabs(*raw.sound_horizon_mpc - exact) <= 1e-9L + 2e-11L * exact,
          "radiation exact limit");
    auto base = eval(typical());
    auto x = typical();
    x.model.h0_km_s_mpc *= 2;
    auto doubled = eval(x);
    finite(doubled);
    check(std::fabs(*doubled.sound_horizon_mpc - *base.sound_horizon_mpc / 2) <=
              1e-9,
          "inverse H0 scaling");
    x = typical();
    x.z_drag *= 2;
    auto later = eval(x);
    finite(later);
    check(*later.sound_horizon_mpc < *base.sound_horizon_mpc,
          "drag-redshift monotonicity");
    x = typical();
    x.model.omega_b *= 2;
    auto loaded = eval(x);
    finite(loaded);
    check(*loaded.sound_horizon_mpc < *base.sound_horizon_mpc,
          "baryon loading monotonicity");
    x = typical();
    x.model.omega_m = 1;
    x.model.omega_r = std::numeric_limits<double>::denorm_min();
    x.model.omega_gamma = x.model.omega_r;
    rejected(x, numerics::Status::outside_domain);
    x = typical();
    x.model.omega_r = 1;
    x.model.omega_m = std::numeric_limits<double>::denorm_min();
    x.model.omega_b = 0;
    x.model.omega_gamma = 1;
    rejected(x, numerics::Status::outside_domain);
    x = typical();
    x.model.omega_b = x.model.omega_m + 0.1;
    rejected(x, numerics::Status::outside_domain);
    x = typical();
    x.model.omega_gamma = 2 * x.model.omega_r;
    rejected(x, numerics::Status::outside_domain);
    x = typical();
    x.model.omega_r = 0;
    rejected(x, numerics::Status::outside_domain);
    x = typical();
    x.drag_origin.clear();
    rejected(x, numerics::Status::outside_domain);
    x = typical();
    x.z_drag = -1;
    rejected(x, numerics::Status::outside_domain);
    x = typical();
    x.model.h0_km_s_mpc = std::numeric_limits<double>::quiet_NaN();
    rejected(x, numerics::Status::nonfinite_input);
    auto tight = policy();
    tight.absolute_tolerance_mpc = 1e-30;
    tight.relative_tolerance = 0;
    auto impossible = eval(radiation, tight);
    check(impossible.status == numerics::Status::conditioning_budget_exceeded &&
              !impossible.sound_horizon_mpc,
          "unattainable tolerance precise failure");
    check(impossible.callbacks > 0, "failed computed work retained");
    auto limited = policy();
    limited.maximum_callbacks_per_point = 3;
    limited.maximum_total_callbacks = 3;
    std::array two{typical(), typical()};
    auto failed = evaluate_sound_horizon(two, limited);
    check(failed.status == numerics::Status::ok && failed.rows.size() == 2,
          "owned exhausted model attempts");
    check(failed.callbacks == 3 && failed.rows[0].callbacks == 3 &&
              failed.rows[1].callbacks == 0,
          "global actual failed callback charging");
    for (auto &r : failed.rows)
      check(r.status == numerics::Status::work_limit && !r.sound_horizon_mpc,
            "work failure no payload");
    auto five = policy();
    five.maximum_total_callbacks = 5;
    std::array exact_two{radiation, radiation};
    auto partial = evaluate_sound_horizon(exact_two, five);
    check(partial.callbacks == 5 && partial.rows.size() == 2,
          "five callback global cap");
    finite(partial.rows[0]);
    check(partial.rows[1].status == numerics::Status::work_limit &&
              !partial.rows[1].sound_horizon_mpc,
          "later model budget not reset");
    auto bound = sound_horizon_payload_bound(std::span(&radiation, 1));
    check(bound.has_value(), "payload bound");
    auto capped = policy();
    capped.maximum_native_bytes = *bound;
    check(evaluate_sound_horizon(std::span(&radiation, 1), capped).status ==
              numerics::Status::ok,
          "exact payload admitted");
    capped.maximum_native_bytes = *bound - 1;
    auto denied = evaluate_sound_horizon(std::span(&radiation, 1), capped);
    check(denied.status == numerics::Status::work_limit &&
              denied.rows.empty() && denied.callbacks == 0,
          "cap minus one prework failure");
    auto empty_policy = policy();
    empty_policy.maximum_points = empty_policy.maximum_total_callbacks =
        empty_policy.maximum_native_bytes = 0;
    auto empty = evaluate_sound_horizon({}, empty_policy);
    check(empty.status == numerics::Status::ok && empty.rows.empty() &&
              empty.callbacks == 0,
          "empty admitted zero dynamic payload");
    check(sound_horizon_payload_bound({}) == 0, "empty helper zero scope");
    SoundHorizonBatch owned;
    {
      auto r = typical();
      r.drag_origin = std::string(300, 'x');
      owned = evaluate_sound_horizon(std::span(&r, 1), policy());
      r.drag_origin = "changed";
      r.model.omega_m = 0.8;
    }
    check(owned.rows[0].source.drag_origin == std::string(300, 'x'),
          "source string owns lifetime");
    check(std::bit_cast<std::uint64_t>(owned.rows[0].source.model.omega_m) ==
              std::bit_cast<std::uint64_t>(0.3),
          "source model bits own lifetime");
    x = radiation;
    x.model.h0_km_s_mpc = std::numeric_limits<double>::denorm_min();
    auto overflow = eval(x);
    check(overflow.status != numerics::Status::ok &&
              !overflow.sound_horizon_mpc,
          "true final overflow rejected");
    x = radiation;
    x.model.h0_km_s_mpc = std::numeric_limits<double>::max();
    x.z_drag = std::numeric_limits<double>::max();
    auto underflow = eval(x);
    check(underflow.status != numerics::Status::ok &&
              !underflow.sound_horizon_mpc,
          "positive final underflow not zero");
    const int saved = std::fegetround();
    std::fesetround(FE_DOWNWARD);
    auto unsupported =
        evaluate_sound_horizon(std::span(&radiation, 1), policy());
    std::fesetround(saved);
    check(unsupported.status == numerics::Status::invalid_input &&
              unsupported.rows.empty(),
          "unsupported rounding rejected");
    std::cout << "PASS " << checks
              << " independent sound-horizon controls max_error="
              << (double)max_error
              << " max_refinement=" << (double)max_refinement << '\n';
  } catch (const std::exception &e) {
    std::fesetround(FE_TONEAREST);
    std::cerr << "FAIL check " << checks << ": " << e.what() << '\n';
    return 1;
  }
}
