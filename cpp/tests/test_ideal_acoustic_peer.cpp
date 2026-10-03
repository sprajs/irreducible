#include "../src/ideal_acoustic_equations.hpp"
#include "irred/ideal_acoustic.hpp"
#include "irred/quantities.hpp"
#if defined(IRRED_IDEAL_ACOUSTIC_MPFR_REFERENCE)
#include "ideal_acoustic_peer.hpp"
#endif
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <span>
#include <streambuf>
#include <string_view>
#include <utility>

namespace {
using W = long double;
using S = irred::numerics::Status;
namespace native = irred::cosmology;
constexpr W c_km_s = W(irred::speed_of_light_m_per_s) / 1000;
constexpr std::array<std::string_view, 9> fields{
    "Delta_c", "Delta_b",        "Theta0",         "phi",  "delta_c",
    "delta_b", "theta_A_over_k", "theta_c_over_k", "phi_N"};

void need(bool condition, const char *message) {
  if (!condition) {
    std::cerr << "FAIL ideal acoustic peer: " << message << '\n';
    std::exit(1);
  }
}
void close(W actual, W expected, W tolerance, const char *message) {
  if (!std::isfinite(actual) || !std::isfinite(expected) ||
      !std::isfinite(tolerance) || !(tolerance >= 0) ||
      std::abs(actual - expected) >
          tolerance * std::max(1.L, std::abs(expected))) {
    std::cerr << std::setprecision(std::numeric_limits<W>::max_digits10)
              << actual << " vs " << expected << '\n';
    need(false, message);
  }
}

// Portable mathematical controls only. This polynomial is explicitly a formal
// source equation, never a replacement for the production retained background.
// There is no momentum Einstein source in this independent eta RHS.
struct FormalSource {
  W h0, gamma, baryon, cdm, lambda;
  bool dust = false;
};
struct FormalEpoch {
  W a, h, ell, fg, fb, fc, loading;
};
using Y = std::array<W, 7>; // a,phi,phi_eta,y,py,delta_c,theta_c
FormalEpoch formal_epoch(FormalSource s, W a) {
  const W a2 = a * a, a4 = a2 * a2;
  const W matter = s.baryon + s.cdm;
  const W p = s.gamma + matter * a + s.lambda * a4;
  need(a > 0 && p > 0, "formal epoch positive domain");
  const W loading = s.gamma > 0 ? 3 * s.baryon * a / (4 * s.gamma) : 0;
  need(s.gamma > 0 || (s.dust && s.gamma == 0 && s.baryon == 0),
       "radiation-free control is explicitly dust");
  return {a,
          s.h0 / c_km_s * std::sqrt(p) / a,
          (matter * a + 4 * s.lambda * a4) / (2 * p) - 1,
          s.gamma / p,
          s.baryon * a / p,
          s.cdm * a / p,
          loading};
}
Y formal_rhs(FormalSource s, W k, const Y &v) {
  const auto e = formal_epoch(s, v[0]);
  return {v[0] * e.h,
          v[2],
          -3 * e.h * v[2] - (2 * e.ell + 1) * e.h * e.h * v[1] +
              2 * e.h * e.h * e.fg * (v[3] + v[1]),
          v[4] / (1 + e.loading),
          -k * k * (v[3] + (2 + e.loading) * v[1]) / 3,
          -v[6] + 3 * v[2],
          -e.h * v[6] + k * k * v[1]};
}
std::array<W, 9> formal_outputs(FormalSource s, W k, const Y &v) {
  const auto e = formal_epoch(s, v[0]);
  const W temperature = v[3] + v[1], theta = -3 * v[4] / (1 + e.loading);
  return {v[5] + 3 * e.h * v[6] / (k * k),
          3 * temperature + 3 * e.h * theta / (k * k),
          temperature,
          v[1],
          v[5],
          3 * temperature,
          theta / k,
          v[6] / k,
          v[2] / e.h};
}
struct Residuals {
  W momentum, hamiltonian, momentum_N, hamiltonian_N;
};
Residuals formal_residuals(FormalSource s, W k, const Y &v) {
  const auto e = formal_epoch(s, v[0]);
  const auto f = formal_rhs(s, k, v);
  const W z = v[2] / e.h, z_N = f[2] / (e.h * e.h) - e.ell * z;
  const W theta = -3 * v[4] / (1 + e.loading);
  const W theta_prime =
      -3 * f[4] / (1 + e.loading) +
      3 * v[4] * e.h * e.loading / ((1 + e.loading) * (1 + e.loading));
  const W vc = e.h * v[6] / (k * k), va = e.h * theta / (k * k);
  const W vc_N = e.ell * vc + f[6] / (k * k);
  const W va_N = e.ell * va + theta_prime / (k * k);
  const W q = 2 * (e.ell + 1), B = e.fb + 4 * e.fg / 3;
  const W fc_N = (1 - q) * e.fc, fb_N = (1 - q) * e.fb, fg_N = -q * e.fg;
  const W B_N = (1 - q) * B - 4 * e.fg / 3;
  const W dg = 4 * (v[3] + v[1]), db = 3 * (v[3] + v[1]);
  const W dg_N = 4 * (f[3] + f[1]) / e.h, db_N = 3 * (f[3] + f[1]) / e.h;
  const W A = e.fc * v[5] + e.fb * db + e.fg * dg;
  const W A_N = fc_N * v[5] + e.fc * f[5] / e.h + fb_N * db + e.fb * db_N +
                fg_N * dg + e.fg * dg_N;
  const W x2 = k * k / (e.h * e.h);
  return {z + v[1] - 1.5L * (e.fc * vc + B * va),
          x2 * v[1] + 3 * (z + v[1]) + 1.5L * A,
          z_N + z - 1.5L * (fc_N * vc + e.fc * vc_N + B_N * va + B * va_N),
          -2 * e.ell * x2 * v[1] + x2 * z + 3 * (z_N + z) + 1.5L * A_N};
}

// For |u|<=1 this alternating tracer series decreases term by term; the
// first omitted term bounds the exact analytic tail. Floating evaluation is
// checked separately and is not called an outward interval certificate.
std::pair<W, W> tracer_series(W u, unsigned terms) {
  need(std::abs(u) <= 1 && terms > 0, "tracer series declared domain");
  W term = u * u / 12, sum = 0;
  for (unsigned n = 1; n <= terms; ++n) {
    sum += term;
    term *= -u * u * W(n) / (W(n + 1) * W(2 * n + 2) * W(2 * n + 3));
  }
  return {sum, std::abs(term)};
}
Y radiation_state(W a, W k, W beta, W &phi_second) {
  constexpr W p0 = -2.L / 3;
  const W x = k * a / beta, u = x / std::sqrt(3.L), h = beta / a;
  W term = p0, phi = 0, z = 0, phi_xx = 0;
  for (unsigned n = 0; n < 28; ++n) {
    phi += term;
    z += W(2 * n) * term;
    if (n)
      phi_xx += W(2 * n) * W(2 * n - 1) * term / (x * x);
    term *= -u * u * W(n + 2) / (W(n + 1) * W(2 * n + 4) * W(2 * n + 5));
  }
  const W vc = 3 * p0 * (1 - std::sin(u) / u) / (u * u);
  const W va = (z + phi) / 2, dg = -2 * (z + phi) - 2 * x * x * phi / 3;
  const auto integral = tracer_series(u, 28);
  const W dc = 3 * phi - 4.5L * p0 - 9 * p0 * integral.first;
  phi_second = k * k * phi_xx;
  return {a,  phi,           h * z, dg / 4 - phi, -k * k * va / (3 * h),
          dc, k * k * vc / h};
}

void portable_controls() {
  const W tolerance = 4096 * std::numeric_limits<W>::epsilon();
  const FormalSource radiation{c_km_s, 1, 0, 0, 0};
  W phi_second = 0;
  const auto v = radiation_state(.25L, 1, 1, phi_second);
  const auto f = formal_rhs(radiation, 1, v);
  const auto out = formal_outputs(radiation, 1, v);
  const W u = .25L / std::sqrt(3.L), p0 = -2.L / 3;
  const W exact_phi = 3 * p0 * (std::sin(u) - u * std::cos(u)) / (u * u * u);
  close(out[3], exact_phi, tolerance,
        "radiation potential independent closed form");
  close(f[2], phi_second, tolerance,
        "radiation analytic metric satisfies original TRACE");
  const auto residual = formal_residuals(radiation, 1, v);
  close(residual.momentum, 0, tolerance, "radiation momentum constraint");
  close(residual.hamiltonian, 0, tolerance, "radiation Hamiltonian constraint");
  close(out[5], 3 * out[2], tolerance, "radiation baryon entropy");
  close(out[0], out[4] + 12 * out[7], tolerance, "CDM comoving readout");
  close(out[1], out[5] + 12 * out[6], tolerance, "baryon comoving readout");
  close(out[8], .25L * v[2], tolerance, "independent derivative readout");
  const auto i14 = tracer_series(u, 14), i28 = tracer_series(u, 28);
  need(std::abs(i14.first - i28.first) <=
           i14.second + tolerance * std::abs(i28.first),
       "tracer analytic tail and separate floating sensitivity");

  const FormalSource loaded{c_km_s, 1, 4.L / 3, .5L, 0};
  const auto e = formal_epoch(loaded, 1);
  const W phi = -.6L, amplitude = .2L, sine_amplitude = -.1L, phase = .3L,
          k = 2;
  const W omega = k / std::sqrt(3 * (1 + e.loading));
  const W oscillation =
      amplitude * std::cos(phase) + sine_amplitude * std::sin(phase);
  Y oscillator{
      1,
      phi,
      0,
      -(2 + e.loading) * phi + oscillation,
      (1 + e.loading) * omega *
          (-amplitude * std::sin(phase) + sine_amplitude * std::cos(phase)),
      0,
      0};
  auto acoustic = formal_rhs(loaded, k, oscillator);
  close(acoustic[4], -k * k * oscillation / 3, tolerance,
        "loaded frozen oscillator pressure sign");
  close(oscillator[3] + 2 * phi, -e.loading * phi + oscillation, tolerance,
        "loaded gravitational zero point");
  const W theta = -3 * oscillator[4] / (1 + e.loading);
  const W theta_prime =
      -3 * acoustic[4] / (1 + e.loading) +
      3 * oscillator[4] * e.h * e.loading / ((1 + e.loading) * (1 + e.loading));
  close(theta_prime,
        -e.h * e.loading / (1 + e.loading) * theta +
            k * k * ((oscillator[3] + phi) / (1 + e.loading) + phi),
        tolerance, "actual loading derivative and weighted common Euler");
  oscillator[3] = -(2 + e.loading) * phi;
  oscillator[4] = 0;
  close(formal_rhs(loaded, k, oscillator)[4], 0, tolerance,
        "loaded forced oscillator equilibrium");

  // Arbitrary off-constraint state: the signed original residual propagation
  // is checked through differentiated primitive definitions, without resets.
  const Y injected{1, -.6L, .17L, 1.2L, -.08L, 1.3L, -.04L};
  const auto r = formal_residuals(loaded, .7L, injected);
  const W x2 = .7L * .7L / (e.h * e.h);
  need(std::abs(r.momentum) > .01L && std::abs(r.hamiltonian) > .01L,
       "injected constraints are meaningful nonzero controls");
  close(r.momentum_N, -(e.ell + 2) * r.momentum, tolerance,
        "signed momentum residual propagation");
  close(r.hamiltonian_N, -(1 + 2 * e.ell) * r.hamiltonian + x2 * r.momentum,
        tolerance, "signed Hamiltonian residual propagation");

  const FormalSource mixed{c_km_s, 1, .125L, .375L, 0};
  const W a = .01L, m = .5L * a;
  const auto em = formal_epoch(mixed, a);
  const auto raw_matter = [&](bool second_order) {
    const W phi_m = -2.L / 3 + m / 24 - (second_order ? 7 * m * m / 240 : 0);
    const W z = m / 24 - (second_order ? 7 * m * m / 120 : 0), V = -1 - phi_m;
    return Y{a,
             phi_m,
             em.h * z,
             1,
             -(1 + em.loading) * V / (3 * em.h),
             3 * (1 + phi_m),
             V / em.h};
  };
  const auto first = formal_residuals(mixed, 1, raw_matter(false));
  const auto second = formal_residuals(mixed, 1, raw_matter(true));
  close(first.momentum, 7 * m * m / (48 * (1 + m)), tolerance,
        "preserved first-order TRACE start residual");
  close(second.momentum, -21 * m * m * m / (160 * (1 + m)), tolerance,
        "second-order start leaves signed cubic residual");
  need(std::abs(second.momentum) < std::abs(first.momentum),
       "higher matter order reduces this initial residual");

  const FormalSource dust{c_km_s, 0, 0, 1, 0, true};
  const W ad = .25L;
  const auto ed = formal_epoch(dust, ad);
  const W delta = 2 * ad / 5, V = -2.L / 5;
  const Y dust_state{ad, -3.L / 5, 0, 0, 0, delta - 3 * V, V / ed.h};
  const auto fd = formal_rhs(dust, 1, dust_state);
  close(fd[1], 0, tolerance, "EdS constant potential");
  close(fd[2], 0, tolerance, "EdS original TRACE");
  close(dust_state[5] + 3 * ed.h * dust_state[6], 2 * ad / 5, tolerance,
        "EdS growing comoving density");
  close(fd[5], -dust_state[6], tolerance, "EdS density continuity");

  // This binds the actual native coefficient to an independently differentiated
  // Hamiltonian identity. The preserved 3B alternative must fail that identity.
  std::array<W, 5> ny{.08L, -.02L, .03L, -.4L, -.6L}, dy{};
  const W F = e.fc + e.fb + 4 * e.fg / 3, B = e.fb + 4 * e.fg / 3;
  native::detail::IdealAcousticCoefficients coeff{
      x2, e.ell, F, B, e.loading / (1 + e.loading), 1 / (3 * (1 + e.loading)),
      0};
  native::detail::ideal_acoustic_derivative(ny, coeff, dy);
  const W q = 2 * (e.ell + 1), FN = (1 - q) * F - 4 * e.fg / 3;
  const W BN = (1 - q) * B - 4 * e.fg / 3;
  const W C = x2 * ny[4] + 1.5L * (F * ny[0] - B * ny[1] + 3 * B * ny[2]);
  const auto derivative = [&](const auto &d) {
    return -2 * e.ell * x2 * ny[4] + x2 * d[4] +
           1.5L * (FN * ny[0] + F * d[0] - BN * ny[1] - B * d[1] +
                   3 * BN * ny[2] + 3 * B * d[2]);
  };
  close(derivative(dy), -(1 + 2 * e.ell) * C, tolerance,
        "actual native 9B/2 Hamiltonian coefficient");
  auto rejected = dy;
  rejected[0] -= 1.5L * B * ny[2];
  close(derivative(rejected) - derivative(dy), -9 * F * B * ny[2] / 4,
        tolerance, "rejected 3B source remains visible");
  need(std::abs(derivative(rejected) + (1 + 2 * e.ell) * C) > .001L,
       "wrong coefficient control cannot accidentally pass");
}

native::IdealAcousticTransfer native_owner() {
  native::IdealAcousticRequest request{
      {70, .02, .10, 2.7, 0, {}},
      1e-10,
      "synthetic:"
      "e8bf08cadf571557a9e4e2f2b59121dc58a64362288a4b95f4b491d41e35ab36"};
  auto owner = native::prepare_ideal_acoustic(std::move(request));
  need(owner.status() == S::ok && owner.source() && owner.physical_mapping() &&
           owner.background(),
       "actual original physical source preparation");
  need(owner.preparation_work().physical_mappings == 1 &&
           owner.preparation_work().background_preparations == 1,
       "one original map and retained background");
  auto moved = std::move(owner);
  need(owner.status() != S::ok && !owner.source() &&
           !owner.physical_mapping() && !owner.background(),
       "move invalidates prior source accessors");
  return moved;
}
void native_refusal_controls(native::IdealAcousticTransfer &owner) {
  const double ai = owner.source()->initial_scale_factor;
  const std::array<double, 2> duplicate_k{1e-4, 1e-4}, early_a{ai, 2 * ai};
  const auto grid =
      owner.evaluate(duplicate_k, early_a, native::acoustic_all_outputs);
  need(grid.rows.size() == 4 && grid.trajectories.size() == 2,
       "ordered early grid retains repeated k");
  need(grid.shared_dependency_status == S::conditioning_budget_exceeded,
       "missing complete source dependency refuses admission");
  for (std::size_t i = 0; i < 2; ++i) {
    need(grid.trajectories[i].attempts_started == 5,
         "five native attempts retained");
    for (std::size_t j = 0; j < 2; ++j) {
      const auto &row = grid.rows[2 * i + j];
      need(row.original_k_index == i && row.original_a_index == j &&
               row.wavenumber_mpc_inverse == duplicate_k[i] &&
               row.requested_scale_factor == early_a[j],
           "original k-major a-minor lineage");
      need(row.epoch.status == S::ok && row.epoch.hcal_mpc_inverse > 0,
           "finite raw state has its matching epoch");
      for (const auto &field : row.outputs)
        need(field.computed && std::isfinite(*field.computed) && !field.value &&
                 field.status == S::conditioning_budget_exceeded &&
                 !field.error.common_source_background_age &&
                 !field.error.arithmetic_storage_constraint &&
                 !field.error.absolute_error_estimate,
             "raw central witness cannot invent admitted errors or value");
      close(*row.outputs[5].computed, 3 * *row.outputs[2].computed,
            4096 * std::numeric_limits<W>::epsilon(),
            "actual native baryon entropy readout");
    }
  }
  const unsigned mask = native::acoustic_intrinsic_temperature |
                        native::acoustic_metric_log_derivative;
  const auto selected =
      owner.evaluate(std::span(duplicate_k).first(1), early_a, mask);
  need(selected.rows.size() == 2, "selected grid shape");
  for (std::size_t j = 0; j < 2; ++j)
    for (unsigned f = 0; f < 9; ++f) {
      if (mask & (1u << f))
        need(selected.rows[j].outputs[f].computed ==
                 grid.rows[j].outputs[f].computed,
             "requested mask retains the same central state");
      else
        need(!selected.rows[j].outputs[f].computed &&
                 !selected.rows[j].outputs[f].value,
             "unrequested coordinates remain absent");
    }
  const std::array<double, 2> unsorted{2 * ai, ai}, repeated_a{ai, ai};
  need(owner.evaluate(duplicate_k, unsorted, mask).status == S::outside_domain,
       "a axis is not silently sorted");
  need(owner.evaluate(duplicate_k, repeated_a, mask).status ==
           S::outside_domain,
       "a axis duplicates refuse");
  need(owner.evaluate(duplicate_k, early_a, 0).status == S::invalid_input,
       "empty output mask refuses");
  need(owner.evaluate(duplicate_k, early_a, 512).status == S::invalid_input,
       "unknown output mask refuses");
  // The owner needs three initial quadrature callbacks. An exhausted age or
  // background allowance is a resource refusal before a generic invalid
  // integration policy; all original requested rows remain withheld.
  for (unsigned category = 0; category < 2; ++category)
    for (std::size_t limit = 0; limit < 3; ++limit) {
      auto capped = native::IdealAcousticPolicy{};
      if (category == 0)
        capped.maximum_background_evaluations = limit;
      else
        capped.maximum_age_evaluations = limit;
      const auto stopped = owner.evaluate(duplicate_k, early_a, mask, capped);
      need(stopped.status == S::work_limit &&
               stopped.shared_dependency_status == S::work_limit,
           "exhausted native quadrature allowance is work_limit");
      need(stopped.evaluation_work.background_evaluations == 0 &&
               stopped.evaluation_work.age_evaluations == 0,
           "subminimum native age/background caps issue no queries");
      need(stopped.rows.size() == 4 && stopped.trajectories.size() == 2,
           "native refused row and trajectory lineage retained");
      for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j) {
          const auto &row = stopped.rows[2 * i + j];
          need(row.original_k_index == i && row.original_a_index == j &&
                   row.wavenumber_mpc_inverse == duplicate_k[i] &&
                   row.requested_scale_factor == early_a[j] &&
                   row.epoch.status == S::work_limit,
               "native resource refusal retains the ordered unavailable epoch");
          for (unsigned f = 0; f < 9; ++f)
            if (mask & (1u << f))
              need(!row.outputs[f].computed && !row.outputs[f].value &&
                       row.outputs[f].status == S::work_limit,
                   "native resource refusal withholds fields");
        }
    }

  for (unsigned variant = 0; variant < 4; ++variant) {
    native::IdealAcousticRequest bad{
        {70, .02, .10, 2.7, 0, {}}, 1e-10, "invalid-model-control"};
    if (variant == 0)
      bad.model.physical_baryon_density = 0;
    if (variant == 1)
      bad.model.tcmb_kelvin = 0;
    if (variant == 2)
      bad.model.physical_massless_nonphoton_density = 1e-5;
    if (variant == 3)
      bad.model.species.push_back({0, 1, 1});
    auto refused = native::prepare_ideal_acoustic(std::move(bad));
    need(refused.status() == S::outside_domain && refused.source(),
         "excluded physical source retained with refusal");
  }
  native::IdealAcousticRequest preserved{
      {70, .02, .10, 2.7, 0, {}}, 1e-10, "preflight-original"};
  auto tiny = native::IdealAcousticPolicy{};
  tiny.maximum_native_bytes = 0;
  auto refused = native::prepare_ideal_acoustic(std::move(preserved), tiny);
  need(refused.status() != S::ok &&
           preserved.source_origin == "preflight-original" &&
           preserved.model.physical_baryon_density == .02,
       "preflight refusal preserves the caller request");
}

#if defined(IRRED_IDEAL_ACOUSTIC_MPFR_REFERENCE)
namespace peer = ideal_acoustic_peer;
using T = peer::Real<120>;
void mp_close(const T &actual, const T &expected, const char *message) {
  const T scale = abs(expected) > T(1L) ? abs(expected) : T(1L);
  need(actual.finite() && expected.finite() &&
           abs(actual - expected) <= T(1e-90L) * scale,
       message);
}
struct AnalyticRadiation {
  peer::State<T> state;
  std::array<T, 9> output;
  T phi_second, integral_tail;
};
AnalyticRadiation mp_radiation(const peer::Source &s, const T &a, const T &k) {
  const T beta = T(s.h0) / peer::c_km_s<T>(), h = beta / a;
  const T x = k * a / beta, u2 = x * x / T(3L), p0 = -peer::rat<T>(2, 3);
  T term(p0), phi, z, phi_xx;
  for (long n = 0; n < 40; ++n) {
    phi = phi + term;
    z = z + T(2 * n) * term;
    if (n)
      phi_xx = phi_xx + T(2 * n) * T(2 * n - 1) * term / (x * x);
    term = -term * u2 * T(n + 2) / (T(n + 1) * T(2 * n + 4) * T(2 * n + 5));
  }
  T vc_term = p0 / T(2L), vc;
  for (long n = 0; n < 40; ++n) {
    vc = vc + vc_term;
    vc_term = -vc_term * u2 / (T(2 * n + 4) * T(2 * n + 5));
  }
  T i_term = u2 / T(12L), integral;
  for (long n = 1; n <= 40; ++n) {
    integral = integral + i_term;
    i_term = -i_term * u2 * T(n) / (T(n + 1) * T(2 * n + 2) * T(2 * n + 3));
  }
  const T va = (z + phi) / T(2L);
  const T dg = -T(2L) * (z + phi) - peer::rat<T>(2, 3) * x * x * phi;
  const T dc = T(3L) * phi - peer::rat<T>(9, 2) * p0 - T(9L) * p0 * integral;
  const T temperature = dg / T(4L), tc = k * k * vc / h;
  return {{a, phi, h * z, temperature - phi, -k * k * va / (T(3L) * h), dc, tc},
          {dc + T(3L) * vc, T(3L) * temperature + T(3L) * va, temperature, phi,
           dc, T(3L) * temperature, x * va, x * vc, z},
          k * k * phi_xx,
          abs(i_term)};
}
peer::Source formal_peer_radiation() {
  peer::Source source;
  source.status = S::ok;
  source.h0 = 299792458.L;
  source.gamma = 1;
  return source; // Deliberately outside the positive-baryon/CDM factory.
}
void mpfr_controls(native::IdealAcousticTransfer &owner) {
  peer::Budget controls;
  {
    peer::BudgetScope scope(controls);
    const T nan(std::numeric_limits<W>::quiet_NaN());
    need(!nan.finite() && !(nan <= T(0L)) && !(nan >= T(0L)),
         "MPFR NaN cannot pass ordered comparison gates");
    const auto source = formal_peer_radiation();
    const T a = peer::rat<T>(1, 4), k(1000L);
    const auto analytic = mp_radiation(source, a, k);
    need(analytic.integral_tail < T(1e-100L),
         "explicit radiation tracer series tail");
    const auto out = peer::outputs(source, k, analytic.state);
    for (std::size_t j = 0; j < out.size(); ++j)
      mp_close(out[j], analytic.output[j], "all9 MPFR radiation readouts");
    const auto f = peer::rhs(source, k, analytic.state);
    mp_close(f[2], analytic.phi_second,
             "MPFR independent analytic radiation TRACE");
    const auto exact_constraints = peer::constraints(source, k, analytic.state);
    mp_close(exact_constraints.momentum, T(0L), "analytic radiation momentum");
    mp_close(exact_constraints.hamiltonian, T(0L),
             "analytic radiation Hamiltonian");

    const auto raw = peer::initial(source, k, a);
    const auto raw_constraints = peer::constraints(source, k, raw);
    const auto e = peer::epoch(source, a);
    const T x2 = k * k / (e.h * e.h), p0 = -peer::rat<T>(2, 3);
    mp_close(raw_constraints.hamiltonian, -p0 * x2 * x2 / T(30L),
             "unprojected radiation initial constraint remains nonzero");
    mp_close(raw_constraints.zeta, T(1L) - p0 * x2 / T(12L),
             "finite-start curvature is not reset to asymptotic zeta");
    need(abs(raw_constraints.hamiltonian) > T(1e-20L) &&
             abs(raw_constraints.zeta - T(1L)) > T(1e-5L),
         "raw initial controls have a resolved violation");

    peer::Source loaded(source);
    loaded.baryon = 4.L / 3;
    loaded.cdm = .5L;
    const auto le = peer::epoch(loaded, a);
    const T phi = -peer::rat<T>(3, 5);
    peer::State<T> equilibrium{a,     phi,   T(0L), -(T(2L) + le.loading) * phi,
                               T(0L), T(0L), T(0L)};
    mp_close(peer::rhs(loaded, k, equilibrium)[4], T(0L),
             "MPFR loaded acoustic gravitational zero point");
    equilibrium[4] = peer::rat<T>(1, 10);
    const auto lf = peer::rhs(loaded, k, equilibrium);
    const T theta = -T(3L) * equilibrium[4] / (T(1L) + le.loading);
    const T theta_prime = -T(3L) * lf[4] / (T(1L) + le.loading) +
                          T(3L) * equilibrium[4] * le.h * le.loading /
                              ((T(1L) + le.loading) * (T(1L) + le.loading));
    mp_close(theta_prime,
             -le.h * le.loading / (T(1L) + le.loading) * theta +
                 k * k * ((equilibrium[3] + phi) / (T(1L) + le.loading) + phi),
             "MPFR weighted common Euler with actual loading derivative");

    const T beta = T(source.h0) / peer::c_km_s<T>();
    const T step = a / beta / T(16L);
    const auto exact_end = mp_radiation(source, a + beta * step, k);
    std::array<T, 3> errors;
    for (unsigned level = 0; level < 3; ++level) {
      auto state = analytic.state;
      const unsigned count = 1u << level;
      for (unsigned j = 0; j < count; ++j)
        state =
            peer::trial(source, k, state, step / T(static_cast<long>(count)))
                .high;
      errors[level] = abs(state[1] - exact_end.state[1]);
    }
    need(errors[1] < errors[0] && errors[2] < errors[1],
         "three independent eta DP meshes refine toward analytic radiation");

    const auto captured =
        peer::Source::capture(*owner.background(), *owner.physical_mapping());
    need(captured.status == S::ok && !captured.formal_dust,
         "actual retained source capture has physical ancestry");
    auto mismatch = *owner.physical_mapping();
    mismatch.model->omega_b *= 2;
    need(peer::Source::capture(*owner.background(), mismatch).status != S::ok,
         "different map cannot supply retained-background error ancestry");
    mismatch = *owner.physical_mapping();
    (*mismatch.scalar_witnesses)[0].emitted_value = 0;
    need(peer::Source::capture(*owner.background(), mismatch).status != S::ok,
         "map witness emitted value must match the retained source");
    need(peer::up_difference(1, 1) == 0 && peer::up_times_eight(0) == 0 &&
             peer::up_add(1, 1) > 2,
         "stored-witness outward measurements preserve exact zero");
  }
  // These are separate bounded mathematical/refusal requests. The original
  // two-row campaign below has one cumulative budget and no retry reset.
  const auto source = formal_peer_radiation();
  const std::array<double, 1> target{.002};
  peer::Budget cap;
  cap.maximum_writes = 0;
  auto refused = peer::run<120>(source, .001, .001L, target, 0, cap);
  need(refused.status == S::work_limit && refused.points.empty() &&
           cap.writes == 0,
       "zero scalar cap produces a retained run refusal");
  cap = peer::Budget{};
  cap.maximum_payload = 0;
  refused = peer::run<120>(source, .001, .001L, target, 0, cap);
  need(refused.status == S::work_limit && cap.background == 0 &&
           refused.mesh.empty(),
       "tiny payload refuses before arithmetic allocation");
  cap = peer::Budget{};
  cap.maximum_rhs = 0;
  refused = peer::run<120>(source, .001, .001L, target, 0, cap);
  need(refused.status == S::work_limit && refused.mesh.size() == 1 &&
           refused.points.empty(),
       "RHS cap retains initialized raw prefix");
  cap = peer::Budget{};
  cap.maximum_endpoints = 1;
  refused = peer::run<120>(source, .001, .001L, target, 0, cap);
  need(refused.status == S::work_limit && cap.endpoints == 1 &&
           refused.mesh.size() == 1,
       "endpoint refusal retains the original initial state");
  cap = peer::Budget{};
  const std::array<double, 2> duplicate_a{.002, .002};
  refused = peer::run<120>(source, .001, .001L, duplicate_a, 0, cap);
  need(refused.status == S::outside_domain && cap.background == 0,
       "peer a axis duplicates refuse before source queries");
  cap = peer::Budget{};
  refused = peer::run<90>(source, std::numeric_limits<double>::quiet_NaN(),
                          .001L, target, 0, cap);
  need(refused.status == S::invalid_input && cap.background == 0,
       "nonfinite peer k refuses with no invented initial state");
  cap = peer::Budget{};
  cap.maximum_payload = 0;
  const auto failed_campaign = peer::campaign(source, .001, .001L, target, cap);
  need(failed_campaign.status == S::work_limit &&
           !failed_campaign.source_admitted &&
           !failed_campaign.resource_admitted &&
           cap.retained_payload_prefix == 0,
       "ten-run payload preflight preserves refusal and caller prefix");
}

class LineBuffer : public std::streambuf {
  std::array<char, 4096> bytes_{};

public:
  LineBuffer() { setp(bytes_.data(), bytes_.data() + bytes_.size() - 1); }
  const char *data() const { return bytes_.data(); }
  std::size_t size() const {
    return static_cast<std::size_t>(pptr() - pbase());
  }
};
class Record : public std::ostream {
  LineBuffer buffer_;

public:
  Record() : std::ostream(nullptr) {
    rdbuf(&buffer_);
    *this << std::boolalpha
          << std::setprecision(std::numeric_limits<W>::max_digits10);
  }
  const char *data() const { return buffer_.data(); }
  std::size_t size() const { return buffer_.size(); }
};
void emit(const Record &record, peer::Budget &budget, bool flush = false) {
  need(record.good(), "fixed record encoder overflow refuses");
  const std::size_t bytes = record.size() + 1;
  if (bytes > budget.maximum_artifact_bytes ||
      budget.artifact_bytes > budget.maximum_artifact_bytes - bytes)
    throw peer::Stop{S::work_limit};
  budget.artifact_bytes += bytes;
  std::cout.write(record.data(), static_cast<std::streamsize>(record.size()));
  std::cout.put('\n');
  if (flush)
    std::cout.flush();
  need(bool(std::cout), "retained diagnostic output failed");
}
void retain_run(const peer::Run &run, std::size_t k_index, std::size_t attempt,
                peer::Budget &budget) {
  Record record;
  record << "{\"kind\":\"peer-attempt\",\"k_index\":" << k_index
         << ",\"attempt\":" << attempt
         << ",\"status\":" << static_cast<int>(run.status)
         << ",\"digits\":" << run.decimal_digits
         << ",\"resolution\":" << run.resolution
         << ",\"initial_a\":" << run.initial_a << ",\"initial_state\":";
  // Availability is earned after all seven raw fields and M/H/zeta have been
  // stored. A refused mesh recorder must not erase that completed IC witness.
  if (!run.initial_state_available)
    record << "null";
  else {
    record << '[';
    for (std::size_t j = 0; j < 7; ++j) {
      if (j)
        record << ',';
      record << run.initial_state[j];
    }
    record << ']';
  }
  record << ",\"initial_M\":";
  if (!run.initial_state_available)
    record << "null";
  else
    record << run.initial_momentum;
  record << ",\"initial_H\":";
  if (!run.initial_state_available)
    record << "null";
  else
    record << run.initial_hamiltonian;
  record << ",\"initial_zeta\":";
  if (!run.initial_state_available)
    record << "null";
  else
    record << run.initial_zeta;
  record << ",\"initial_eta_mpc\":";
  if (!run.initial_state_available)
    record << "null";
  else
    record << run.initial_eta_mpc;
  record << ",\"initial_age_Lambda_bound\":";
  if (!run.initial_state_available)
    record << "null";
  else
    record << run.initial_age_lambda_error_mpc;
  record << ",\"maximum_normalized_M\":" << run.maximum_normalized_momentum
         << ",\"maximum_normalized_H\":" << run.maximum_normalized_hamiltonian
         << ",\"source_admitted\":false}";
  emit(record, budget, true);
  for (std::size_t i = 0; i < run.mesh.size(); ++i) {
    const auto &point = run.mesh[i];
    Record mesh;
    mesh << "{\"kind\":\"peer-mesh\",\"k_index\":" << k_index
         << ",\"attempt\":" << attempt << ",\"index\":" << i << ",\"state\":[";
    for (std::size_t j = 0; j < 7; ++j) {
      if (j)
        mesh << ',';
      mesh << point.state[j];
    }
    mesh << "],\"eta_elapsed\":" << point.eta_elapsed
         << ",\"M\":" << point.momentum << ",\"H\":" << point.hamiltonian
         << ",\"normalized_M\":" << point.normalized_momentum
         << ",\"normalized_H\":" << point.normalized_hamiltonian << '}';
    emit(mesh, budget);
  }
}

void retain_native(const native::IdealAcousticBatch &batch,
                   peer::Budget &budget) {
  Record returned;
  returned << "{\"kind\":\"native-request-returned\",\"status\":"
           << static_cast<int>(batch.status) << ",\"dependency_status\":"
           << static_cast<int>(batch.shared_dependency_status)
           << ",\"attempts\":" << batch.evaluation_work.attempts
           << ",\"rhs\":" << batch.evaluation_work.rhs_evaluations
           << ",\"P_calls\":" << batch.evaluation_work.background_evaluations
           << ",\"age_calls\":" << batch.evaluation_work.age_evaluations
           << ",\"writes\":" << batch.evaluation_work.state_element_writes
           << ",\"predictions_admitted\":false}";
  emit(returned, budget, true);
  for (const auto &trajectory : batch.trajectories)
    for (std::size_t j = 0; j < trajectory.attempts.size(); ++j) {
      const auto &attempt = trajectory.attempts[j];
      Record record;
      record << "{\"kind\":\"native-attempt\",\"k_index\":"
             << trajectory.original_k_index << ",\"attempt\":" << j
             << ",\"started\":" << (j < trajectory.attempts_started)
             << ",\"status\":" << static_cast<int>(attempt.status)
             << ",\"initial_a\":" << attempt.initial_scale_factor
             << ",\"log_step\":" << attempt.maximum_log_step << ",\"raw_IC\":";
      if (!attempt.unprojected_initial_state)
        record << "null";
      else {
        record << '[';
        for (std::size_t q = 0; q < 5; ++q) {
          if (q)
            record << ',';
          record << (*attempt.unprojected_initial_state)[q];
        }
        record << ']';
      }
      record << ",\"projected_IC\":";
      if (!attempt.projected_initial_state)
        record << "null";
      else {
        record << '[';
        for (std::size_t q = 0; q < 5; ++q) {
          if (q)
            record << ',';
          record << (*attempt.projected_initial_state)[q];
        }
        record << ']';
      }
      record << ",\"Delta_projection\":";
      if (attempt.initial_delta_projection)
        record << *attempt.initial_delta_projection;
      else
        record << "null";
      record << ",\"maximum_normalized_H\":"
             << attempt.maximum_normalized_hamiltonian_residual
             << ",\"maximum_absolute_H\":"
             << attempt.maximum_absolute_hamiltonian_residual
             << ",\"maximum_direct_H\":"
             << attempt.maximum_direct_hamiltonian_residual
             << ",\"maximum_H_assembly_discrepancy\":"
             << attempt.maximum_hamiltonian_assembly_discrepancy
             << ",\"maximum_closure_defect\":"
             << attempt.maximum_absolute_closure_defect
             << ",\"maximum_acceleration_defect\":"
             << attempt.maximum_absolute_acceleration_defect
             << ",\"rhs\":" << attempt.work.rhs_evaluations
             << ",\"P_calls\":" << attempt.work.background_evaluations
             << ",\"writes\":" << attempt.work.state_element_writes << '}';
      emit(record, budget);
    }
}

void retained_source_control(native::IdealAcousticTransfer &owner,
                             const peer::Source &source, double a, double k,
                             peer::Budget &budget) {
  peer::BudgetScope scope(budget);
  // One actual P query, once at this common requested a. This retains its raw
  // receipt and literal helper/source-control difference. A measured difference
  // is not a propagated source radius or an independently mapped cosmology.
  budget.query();
  const auto actual = native::detail::thermal_conformal_epoch(
      *owner.background(), W(a), W(k), source.gamma, native::ThermalPolicy{});
  Record record;
  record << "{\"kind\":\"retained-helper-source-control\",\"a\":" << a
         << ",\"k\":" << k
         << ",\"helper_status\":" << static_cast<int>(actual.status)
         << ",\"raw_query_available\":" << actual.raw_scaled_query.has_value();
  if (actual.status != S::ok) {
    record << ",\"source_control\":null,\"full_output_source_bound\":null}";
    emit(record, budget, true);
    need(false, "actual retained helper control query must complete");
  }
  const auto control = peer::epoch(source, T(a));
  const std::array<W, 8> helper{actual.p,  actual.hcal, actual.g,  actual.fg,
                                actual.fb, actual.fc,   actual.fl, actual.fr};
  const std::array<T, 8> high_state{control.p,  control.h,  control.ell,
                                    control.fg, control.fb, control.fc,
                                    control.fl, control.fg};
  std::array<W, 8> high{};
  for (std::size_t j = 0; j < high.size(); ++j) {
    peer::write();
    high[j] = high_state[j].wide();
  }
  const std::array<std::string_view, 8> names{"P",  "Hcal", "ell", "fg",
                                              "fb", "fc",   "fl",  "fr"};
  record << ",\"same_a_differences\":{ ";
  for (std::size_t j = 0; j < helper.size(); ++j) {
    if (j)
      record << ',';
    record << '"' << names[j]
           << "\":" << peer::up_difference(helper[j], high[j]);
  }
  record << "},\"source_control_emission_losses\":[";
  for (std::size_t j = 0; j < high.size(); ++j) {
    if (j)
      record << ',';
    record << peer::emitted_loss(high_state[j], high[j]);
  }
  record
      << "],\"actual_P_error_estimate\":" << actual.p_error
      << ",\"actual_closure_defect\":" << actual.closure_defect
      << ",\"actual_acceleration_defect\":" << actual.acceleration_defect
      << ",\"source_normalization_estimate\":" << source.normalization_estimate
      << ",\"Lambda_getter_signed_loss\":" << source.lambda_getter_loss
      << ",\"map_estimates\":[" << source.map_estimates[0] << ','
      << source.map_estimates[1] << ',' << source.map_estimates[2] << ']'
      << ",\"source_control_loading\":" << control.loading.wide()
      << ",\"actual_loading_radius\":null,\"forward_available\":"
      << actual.forward.has_value()
      << ",\"shadow_available\":" << actual.shadow.has_value()
      << ",\"full_output_source_bound\":null,\"comparison_qualified\":false}";
  emit(record, budget, true);
  need(actual.raw_scaled_query && actual.momentum_callbacks == 0 &&
           !actual.forward && !actual.shadow,
       "same retained helper has a raw receipt and unavailable complete source "
       "bundles");
}

void original_campaign(native::IdealAcousticTransfer &owner) {
  const std::array<double, 2> ks{1e-4, .01};
  const std::array<double, 1> a{.01};
  const native::IdealAcousticPolicy native_policy;
  peer::Budget budget; // Unchanged whole-request caps across both original k.
  peer::Source source;
  {
    peer::BudgetScope scope(budget);
    source =
        peer::Source::capture(*owner.background(), *owner.physical_mapping());
  }
  need(source.status == S::ok, "original once-captured peer source");
  Record started;
  started
      << "{\"kind\":\"original-request-started\",\"request_sha256\":"
         "\"e8bf08cadf571557a9e4e2f2b59121dc58a64362288a4b95f4b491d41e35ab36\","
         "\"k\":[0.0001,0.01],\"a\":[0.01],\"native_ai\":1e-10,"
         "\"model_id\":\""
      << native::ideal_acoustic_model_id << "\",\"peer_method_id\":\""
      << peer::method_id << "\",\"peer_source_control_id\":\""
      << peer::source_control_id
      << "\",\"H0\":" << owner.source()->model.h0_km_s_mpc
      << ",\"omega_b_physical\":"
      << owner.source()->model.physical_baryon_density
      << ",\"omega_c_physical\":" << owner.source()->model.physical_cdm_density
      << ",\"Tcmb\":" << owner.source()->model.tcmb_kelvin
      << ",\"outputs\":" << native::acoustic_all_outputs
      << ",\"native_max_rhs\":" << native_policy.maximum_rhs_batch
      << ",\"native_max_rhs_per_k\":"
      << native_policy.maximum_rhs_per_wavenumber
      << ",\"native_max_P\":" << native_policy.maximum_background_evaluations
      << ",\"native_max_age\":" << native_policy.maximum_age_evaluations
      << ",\"native_max_writes\":" << native_policy.maximum_state_element_writes
      << ",\"native_max_payload\":" << native_policy.maximum_native_bytes
      << ",\"peer_max_rhs\":" << budget.maximum_rhs
      << ",\"peer_max_rhs_per_k\":" << budget.maximum_rhs_per_k
      << ",\"peer_max_P\":" << budget.maximum_background
      << ",\"peer_max_writes\":" << budget.maximum_writes
      << ",\"peer_max_endpoints\":" << budget.maximum_endpoints
      << ",\"peer_max_payload\":" << budget.maximum_payload
      << ",\"peer_max_artifact\":" << budget.maximum_artifact_bytes
      << ",\"source_admitted\":false,\"external_runtime_admitted\":false}";
  emit(started, budget, true);
  const auto actual =
      owner.evaluate(ks, a, native::acoustic_all_outputs, native_policy);
  retain_native(actual, budget);
  need(actual.rows.size() == ks.size() * a.size(),
       "original native ordered nine-field grid");
  need(actual.shared_dependency_status == S::conditioning_budget_exceeded,
       "native full source gate remains explicitly unavailable");
  need(actual.trajectories.size() == ks.size(),
       "original native trajectories retained");
  for (const auto &trajectory : actual.trajectories) {
    need(trajectory.attempts_started == 5,
         "all five original native attempts complete");
    for (std::size_t j = 0; j < trajectory.attempts.size(); ++j) {
      const auto &attempt = trajectory.attempts[j];
      const W expected_start =
          std::ldexp(W(owner.source()->initial_scale_factor),
                     -static_cast<int>(j < 3 ? 0 : j - 2));
      const W expected_step = std::ldexp(W(native_policy.maximum_log_step),
                                         -static_cast<int>(j < 3 ? j : 2));
      need(attempt.status == S::ok &&
               attempt.initial_scale_factor ==
                   static_cast<double>(expected_start) &&
               attempt.maximum_log_step == static_cast<double>(expected_step) &&
               attempt.unprojected_initial_state &&
               attempt.projected_initial_state &&
               attempt.initial_delta_projection,
           "native original starts, meshes and projection ancestry retained");
    }
  }
  retained_source_control(owner, source, a[0], ks[0], budget);
  for (std::size_t i = 0; i < ks.size(); ++i) {
    Record attempt_started;
    attempt_started << "{\"kind\":\"peer-campaign-started\",\"k_index\":" << i
                    << ",\"k\":" << ks[i] << ",\"attempts_requested\":10}";
    emit(attempt_started, budget, true);
    const auto campaign = peer::campaign(
        source, ks[i], owner.source()->initial_scale_factor, a, budget);
    for (std::size_t r = 0; r < 10; ++r)
      retain_run(campaign.attempts[r], i, r, budget);
    Record summary;
    summary << "{\"kind\":\"peer-campaign-returned\",\"k_index\":" << i
            << ",\"status\":" << static_cast<int>(campaign.status)
            << ",\"rhs\":" << budget.rhs << ",\"P_calls\":" << budget.background
            << ",\"writes\":" << budget.writes
            << ",\"endpoints\":" << budget.endpoints
            << ",\"artifact_bytes\":" << budget.artifact_bytes
            << ",\"source_admitted\":" << campaign.source_admitted
            << ",\"resource_admitted\":" << campaign.resource_admitted << '}';
    emit(summary, budget, true);
    need(!campaign.source_admitted && !campaign.resource_admitted &&
             campaign.status != S::ok,
         "missing full peer error and runtime gates cannot admit a prediction");
    for (std::size_t r = 0; r < 10; ++r) {
      const auto &run = campaign.attempts[r];
      need(run.status == S::ok && run.initial_state_available &&
               run.points.size() == a.size(),
           "mandatory raw reference attempt must complete under original caps");
      need(run.decimal_digits == (r == 9 ? 90u : 120u) &&
               run.resolution == (r == 9 ? 2u : unsigned(r % 3)),
           "original precision and time attempt order");
      const W expected_start =
          std::ldexp(W(owner.source()->initial_scale_factor),
                     -static_cast<int>(r == 9 ? 4 : r / 3 + 2));
      need(run.initial_a == expected_start,
           "original peer starts are preserved");
      for (const auto &point : run.points)
        need(!point.source_age_endpoint_error &&
                 !point.arithmetic_constraint_error,
             "absent full reference errors stay absent");
    }
    need(campaign.time_precision.size() == a.size() &&
             campaign.initial_refinement.size() == a.size(),
         "two resolved time and start comparisons exist for every endpoint");
    const auto &row = actual.rows[i];
    const auto &reference = campaign.attempts[8].points[0];
    need(row.original_k_index == i && row.original_a_index == 0 &&
             row.requested_scale_factor == a[0] &&
             row.wavenumber_mpc_inverse == ks[i],
         "original native comparison lineage");
    need(reference.actual_a <= W(a[0]) && reference.endpoint_a_error >= 0 &&
             std::abs(reference.actual_a - W(a[0])) <=
                 reference.endpoint_a_error +
                     4 * std::numeric_limits<W>::epsilon() * W(a[0]),
         "raw actual endpoint discrepancy is retained");
    for (std::size_t j = 0; j < 9; ++j) {
      const auto &field = row.outputs[j];
      need(field.computed && std::isfinite(*field.computed) && !field.value &&
               field.status == S::conditioning_budget_exceeded &&
               !field.error.common_source_background_age &&
               !field.error.arithmetic_storage_constraint &&
               !field.error.absolute_error_estimate,
           "all9 original native raw witnesses remain unqualified");
      const W epsilon =
          W(native_policy.absolute_tolerance) +
          W(native_policy.relative_tolerance) * std::abs(*field.computed);
      const W time = campaign.time_precision[0][j],
              initial = campaign.initial_refinement[0][j];
      const W gap = std::abs(*field.computed - reference.values[j]);
      Record comparison;
      comparison << "{\"kind\":\"partial-nine-field-comparison\",\"k_index\":"
                 << i << ",\"a_index\":0,\"field\":\"" << fields[j]
                 << "\",\"native_status\":" << static_cast<int>(field.status)
                 << ",\"native_computed\":" << *field.computed
                 << ",\"peer_computed\":" << reference.values[j]
                 << ",\"gap\":" << gap << ",\"time_precision\":" << time
                 << ",\"initial\":" << initial
                 << ",\"native_a\":" << row.epoch.state_scale_factor
                 << ",\"peer_a\":" << reference.actual_a
                 << ",\"native_Hcal\":" << row.epoch.hcal_mpc_inverse
                 << ",\"peer_Hcal\":" << reference.hcal
                 << ",\"native_eta_mpc\":" << row.epoch.conformal_age_mpc
                 << ",\"peer_eta_mpc\":" << reference.conformal_age_mpc
                 << ",\"peer_initial_age_Lambda_bound\":"
                 << reference.initial_age_lambda_error_mpc
                 << ",\"native_time\":" << field.error.time_refinement
                 << ",\"native_initial\":" << field.error.initial_refinement
                 << ",\"peer_normalized_M\":" << reference.normalized_momentum
                 << ",\"peer_normalized_H\":"
                 << reference.normalized_hamiltonian
                 << ",\"peer_endpoint_defect\":" << reference.endpoint_a_error
                 << ",\"epsilon\":" << epsilon
                 << ",\"full_source_error\":null,\"full_arithmetic_endpoint_"
                    "error\":null,"
                    "\"comparison_qualified\":false}";
      emit(comparison, budget);
      need(
          std::isfinite(time) && time >= 0 && time <= epsilon / 20 &&
              std::isfinite(initial) && initial >= 0 && initial <= epsilon / 20,
          "known reference shares satisfy their original separate allocations");
      need(std::isfinite(gap) && gap + time + initial <= epsilon,
           "raw independent nine-field agreement and known error shares");
    }
  }
  need(budget.attempts == 20 && budget.rhs <= budget.maximum_rhs &&
           budget.background <= budget.maximum_background &&
           budget.writes <= budget.maximum_writes &&
           budget.endpoints <= budget.maximum_endpoints,
       "all20 original attempts use one unreset whole-request budget");
}
#endif
} // namespace

int main() {
  portable_controls();
  auto owner = native_owner();
  native_refusal_controls(owner);
#if defined(IRRED_IDEAL_ACOUSTIC_MPFR_REFERENCE)
  try {
    mpfr_controls(owner);
    original_campaign(owner);
  } catch (const ideal_acoustic_peer::Stop &failure) {
    std::cerr << "FAIL optional ideal acoustic reference refused: "
              << static_cast<int>(failure.status)
              << "; no prediction admitted\n";
    return 1;
  }
#endif
  std::cout << "ideal acoustic portable TRACE/algebra controls completed; "
               "scientific source admission remains unavailable\n";
}
