#include "irred/continuous_cmb_projection.hpp"
#include <array>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <numbers>
#include <type_traits>
#include <utility>

namespace {
using namespace irred::projection;
using S = irred::numerics::Status;
using W = long double;
using C = std::complex<W>;
void need(bool ok, const char *why) {
  if (!ok) { std::cerr << "FAIL " << why << '\n'; std::exit(1); }
}
ContinuousCmbSource constant_source(unsigned channel, double amplitude = 1,
                                    double observer = 12, double end = 3) {
  ContinuousCmbSource s;
  s.k_mpc_inverse = {.5, 1}; s.eta_mpc = {1, end};
  s.observer_eta_mpc = observer;
  s.producer_id = "synthetic-constant-split-source";
  s.signed_mode_id = "synthetic-positive-amplitude";
  s.normalization_id = "dimensionless-unit-initial-mode";
  s.t0.assign(4, 0); s.t1.assign(4, 0); s.t2.assign(4, 0); s.polarization.assign(4, 0);
  const std::array<std::vector<double> *, 4> channels{&s.t0, &s.t1, &s.t2, &s.polarization};
  channels[channel]->assign(4, amplitude);
  return s;
}
ContinuousCmbPolicy accurate() {
  ContinuousCmbPolicy p; p.absolute_tolerance = 2e-12; p.relative_tolerance = 2e-10; return p;
}
// Independent algorithm: stable-domain forward spherical recurrence with
// analytic trigonometric seeds, rather than std special-function seeds. System
// trig/constants ancestry remains shared; this is an algorithm control.
W forward_j(unsigned l, W x) {
  if (x == 0) return l ? 0 : 1;
  W a = std::sin(x) / x;
  if (!l) return a;
  W b = (std::sin(x) - x * std::cos(x)) / (x * x);
  for (unsigned n = 1; n < l; ++n) {
    const W next = (2 * n + 1) * b / x - a; a = b; b = next;
  }
  return b;
}
// Independently authored angular projection integrates the finite constant
// eta source analytically, then uses composite Simpson in mu. No Bessel calls.
// Outward exp(+ix mu), source T0+i mu T1-P2(mu)T2, CLASS-positive E.
W angular_constant(unsigned l, unsigned channel, W k, W observer,
                   W eta_start, W eta_end, unsigned panels) {
  const std::array<C, 4> il{C(1, 0), C(0, 1), C(-1, 0), C(0, -1)};
  C sum{};
  for (unsigned i = 0; i <= panels; ++i) {
    const W mu = -1 + 2 * W(i) / panels;
    W p = 1, d = 0, dd = 0, oldp = 0, oldd = 0, olddd = 0;
    for (unsigned n = 1; n <= l; ++n) {
      const W np = ((2 * n - 1) * mu * p - (n - 1) * oldp) / n;
      const W nd = ((2 * n - 1) * (p + mu * d) - (n - 1) * oldd) / n;
      const W ndd = ((2 * n - 1) * (2 * d + mu * dd) - (n - 1) * olddd) / n;
      oldp = p; oldd = d; olddd = dd; p = np; d = nd; dd = ndd;
    }
    const C phase = mu == 0 ? C(eta_end - eta_start, 0) :
        (std::exp(C(0, k * (observer - eta_start) * mu)) -
         std::exp(C(0, k * (observer - eta_end) * mu))) / C(0, k * mu);
    C source = 1;
    if (channel == 1) source = C(0, mu);
    if (channel == 2) source = -(3 * mu * mu - 1) / 2;
    if (channel == 3) source = (1 - mu * mu) * (1 - mu * mu) * dd;
    const W weight = (i == 0 || i == panels) ? 1 : ((i % 2) ? 4 : 2);
    sum += weight * phase * source * (channel == 3 ? 1.L : p);
  }
  sum *= 2.L / (3 * panels);
  C result = sum / (2.L * il[l % 4]);
  if (channel == 3) {
    const W product = W(l - 1) * l * (l + 1) * (l + 2);
    result *= -std::sqrt(3.L * product / 8) / product;
  }
  need(std::abs(result.imag()) < 1e-12L, "independent angular parity/signed convention");
  return result.real();
}
void analytic_and_angular() {
  const std::array<unsigned, 3> ell{4, 0, 1}; // original order is deliberate
  auto owner = prepare_continuous_cmb_projection(constant_source(1));
  need(owner.status() == S::ok, "constant Doppler source acquired");
  const auto result = project_continuous_cmb(owner, ell, continuous_temperature, accurate());
  need(result.status == S::ok && result.rows.size() == 6, "finite constant Doppler transfers");
  for (const auto &r : result.rows) {
    const W k = r.k_mpc_inverse;
    const W exact = (forward_j(r.ell, 11 * k) - forward_j(r.ell, 9 * k)) / k;
    need(r.temperature && !r.e_mode && r.ell == ell[r.multipole_index] &&
         r.k_index < 2 && !r.radial_error_estimate && !r.source_grid_error_estimate,
         "requested masks/order/absent external error ownership");
    need(std::abs(W(*r.temperature) - exact) < 2e-10L,
         "constant Doppler exact finite endpoint integral");
  }
  auto mono = prepare_continuous_cmb_projection(constant_source(0));
  const std::array<unsigned, 1> l1{1};
  const auto m = project_continuous_cmb(mono, l1, continuous_temperature, accurate());
  for (const auto &r : m.rows) {
    const W k = r.k_mpc_inverse;
    need(r.temperature && std::abs(W(*r.temperature) -
        (forward_j(0, 9 * k) - forward_j(0, 11 * k)) / k) < 2e-10L,
        "constant monopole l1 analytic endpoint orientation");
  }
  const std::array<unsigned, 4> angular_ells{2, 4, 8, 32};
  for (unsigned channel = 0; channel < 4; ++channel) {
    auto source = constant_source(channel, 1, 34, 3);
    auto angular_owner = prepare_continuous_cmb_projection(std::move(source));
    const unsigned mask = channel == 3 ? continuous_e_mode : continuous_temperature;
    const auto native = project_continuous_cmb(angular_owner, angular_ells, mask, accurate());
    need(native.status == S::ok, "low/turning-region split-channel projection");
    for (const auto &r : native.rows) {
      const W a = angular_constant(r.ell, channel, r.k_mpc_inverse, 34, 1, 3, 8192);
      const W b = angular_constant(r.ell, channel, r.k_mpc_inverse, 34, 1, 3, 16384);
      const W value = channel == 3 ? *r.e_mode : *r.temperature;
      need(std::abs(a - b) < 2e-7L && std::abs(value - b) < 3e-7L,
           "independent angular/refinement T/E sign and normalization");
    }
  }
}
void regular_small_phase_and_sign() {
  ContinuousCmbSource s = constant_source(2, 1, 1e-12, 1e-12);
  s.k_mpc_inverse = {1}; s.eta_mpc = {0, 1e-12};
  s.t0 = {0, 0}; s.t1 = {0, 0}; s.t2 = {1, 1}; s.polarization = {1, 1};
  auto p = prepare_continuous_cmb_projection(std::move(s));
  const std::array<unsigned, 1> l2{2}, l0{0};
  const auto endpoint = project_continuous_cmb(p, l2, continuous_temperature | continuous_e_mode, accurate());
  need(endpoint.status == S::ok && endpoint.rows[0].temperature && endpoint.rows[0].e_mode,
       "regular finite-support l2 endpoints");
  const W leading = W(1e-12) / 5;
  need(std::abs(W(*endpoint.rows[0].temperature) - leading) < leading * 1e-15L &&
       std::abs(W(*endpoint.rows[0].e_mode) - leading) < leading * 1e-15L,
       "l2 quadrupole/E regular one-fifth limit with finite tiny support");
  const auto tiny = project_continuous_cmb(p, l0, continuous_temperature, accurate());
  const W tiny_expected = W(1e-12) * W(1e-12) * W(1e-12) / 45;
  need(tiny.status == S::ok && tiny.rows[0].temperature && *tiny.rows[0].temperature > 0 &&
       std::abs(W(*tiny.rows[0].temperature) - tiny_expected) < tiny_expected * 1e-15L,
       "positive tiny quadrupole preserved without endpoint cancellation");
  auto positive = prepare_continuous_cmb_projection(constant_source(0, 1));
  auto negative = prepare_continuous_cmb_projection(constant_source(0, -3));
  const auto a = project_continuous_cmb(positive, l2, continuous_temperature, accurate());
  const auto b = project_continuous_cmb(negative, l2, continuous_temperature, accurate());
  need(a.status == S::ok && b.status == S::ok &&
       std::abs(W(*b.rows[0].temperature) + 3 * W(*a.rows[0].temperature)) < 2e-10L,
       "signed source amplitude linearity");
}
void manufactured_source_grid_refinement() {
  // Exact smooth manufactured antiderivative is a DIFFERENT source from each
  // PL table. Its converging mismatch tests producer-grid refinement rather
  // than pretending that the native time-quadrature estimate covers that error.
  const std::array<unsigned, 1> ell{2};
  std::array<W, 2> previous_t{}, previous_e{};
  for (unsigned panels : {32, 64, 128}) {
    ContinuousCmbSource s;
    s.k_mpc_inverse = {.5, 1}; s.observer_eta_mpc = 8;
    s.producer_id = "synthetic-smooth-source-sampled-as-PL";
    s.signed_mode_id = "positive-manufactured-scalar";
    s.normalization_id = "q=.001-in-explicit-Mpc-powers";
    for (unsigned i = 0; i <= panels; ++i) {
      const double eta = 2 + 2. * i / panels, chi = 8 - eta;
      s.eta_mpc.push_back(eta);
      for (unsigned ki = 0; ki < 2; ++ki) {
        s.t0.push_back(.001 * chi * chi * chi * chi);
        s.t1.push_back(0); s.t2.push_back(0);
        s.polarization.push_back(.001 * chi * chi * chi * chi * chi * chi);
      }
    }
    auto owner = prepare_continuous_cmb_projection(std::move(s));
    const auto result = project_continuous_cmb(owner, ell,
        continuous_temperature | continuous_e_mode, accurate());
    need(result.status == S::ok, "manufactured PL-source projection");
    for (const auto &r : result.rows) {
      const W k = r.k_mpc_inverse;
      const W endpoint = 1296 * forward_j(3, 6 * k) - 256 * forward_j(3, 4 * k);
      const W te = std::abs(W(*r.temperature) - .001L * endpoint / k);
      const W ee = std::abs(W(*r.e_mode) - .003L * endpoint / (k * k * k));
      if (panels > 32) need(te < .35L * previous_t[r.k_index] &&
                            ee < .35L * previous_e[r.k_index],
                            "separate PL source-grid refinement against smooth antiderivative");
      previous_t[r.k_index] = te; previous_e[r.k_index] = ee;
      need(!r.source_grid_error_estimate, "source refinement does not invent a per-table grid radius");
    }
  }
}
void lifetime_and_refusals() {
  static_assert(!std::is_copy_constructible_v<ContinuousCmbProjection>);
  auto owner = prepare_continuous_cmb_projection(constant_source(0));
  const auto *source = owner.source();
  auto moved = std::move(owner);
  need(owner.status() == S::invalid_input && !owner.source() && moved.source() == source,
       "move transfers the original immutable source");
  const std::array<unsigned, 2> ell{2, 3};
  auto result = project_continuous_cmb(moved, ell, continuous_temperature, accurate());
  need(result.source_owner.get() == source && result.status == S::ok,
       "result retains the same original source without readback");
  moved = ContinuousCmbProjection{};
  need(result.source_owner.get() == source && result.source_owner->producer_id == "synthetic-constant-split-source",
       "result remains valid after prepared owner release");
  auto bad = constant_source(0); bad.t1.pop_back();
  const auto bad_size = bad.t0.size();
  need(prepare_continuous_cmb_projection(std::move(bad)).status() == S::invalid_input &&
       bad.t0.size() == bad_size, "invalid structure refuses before consumption");
  bad = constant_source(0); bad.eta_mpc[1] = bad.eta_mpc[0];
  need(prepare_continuous_cmb_projection(std::move(bad)).status() == S::outside_domain,
       "duplicate time axis refusal");
  bad = constant_source(0); bad.observer_eta_mpc = 10000;
  need(prepare_continuous_cmb_projection(std::move(bad)).status() == S::outside_domain,
       "unadmitted phase refusal");
  bad = constant_source(0); bad.t2[0] = std::numeric_limits<double>::quiet_NaN();
  need(prepare_continuous_cmb_projection(std::move(bad)).status() == S::nonfinite_input,
       "nonfinite required source refusal");
  auto capped_owner = prepare_continuous_cmb_projection(constant_source(0));
  auto cap = accurate(); cap.maximum_kernel_evaluations = 2;
  const auto capped = project_continuous_cmb(capped_owner, ell, continuous_temperature, cap);
  need(capped.status == S::work_limit && capped.kernel_evaluations == 2 &&
       capped.source_owner.get() == capped_owner.source(), "failed node work and source retained");
  for (const auto &r : capped.rows) need(!r.temperature && !r.e_mode, "work refusal withholds all failed transfers");
  const std::array<unsigned, 1> invalid_e{1}, invalid_l{65};
  need(project_continuous_cmb(capped_owner, invalid_e, continuous_e_mode).status == S::outside_domain &&
       project_continuous_cmb(capped_owner, invalid_l, continuous_temperature).status == S::outside_domain,
       "unavailable polarization/multipole domain");
  cap = accurate(); cap.maximum_payload_bytes = 1;
  need(project_continuous_cmb(capped_owner, ell, continuous_temperature, cap).status == S::work_limit,
       "whole retained evaluation payload refusal before nodes");
  cap = accurate(); cap.maximum_depth = 0;
  need(project_continuous_cmb(capped_owner, ell, continuous_temperature, cap).status == S::conditioning_budget_exceeded,
       "unresolved phase refuses original refinement cap");
}
} // namespace
int main() {
  analytic_and_angular(); regular_small_phase_and_sign();
  manufactured_source_grid_refinement(); lifetime_and_refusals();
  std::cout << "continuous_cmb_projection_contract finite-support/analytic/angular/ownership/refusal controls passed\n";
}
