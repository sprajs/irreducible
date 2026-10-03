#include "irred/massless_fd_transfer.hpp"
#include "payload_accounting.hpp"
#include "thermal_conformal_epoch.hpp"
#include <algorithm>
#include <cfenv>
#include <cmath>
#include <limits>
#include <utility>

namespace irred::cosmology {
namespace {
using S = numerics::Status;
using W = long double;
using Epoch = detail::ThermalConformalEpoch;
// Delta_c, relative entropy, Vr-Vc, Vc, phi, sigma/x^2, eta, F3...FL.
using State = std::array<W, 133>;
constexpr W p0 = -10.L / 19, phi0 = -14.L / 19;
constexpr W eps_w = std::numeric_limits<W>::epsilon();
constexpr unsigned all_fields = massless_fd_comoving_cdm |
                                massless_fd_spatial_potential |
                                massless_fd_lapse_potential;
bool arithmetic() {
  return std::fegetround() == FE_TONEAREST &&
         std::numeric_limits<W>::digits >= 64 &&
         std::numeric_limits<W>::max_exponent >= 16384;
}
bool valid_policy(const MasslessFDTransferPolicy &p) {
  return arithmetic() && std::isfinite(p.absolute_tolerance) &&
         std::isfinite(p.relative_tolerance) && p.absolute_tolerance >= 0 &&
         p.relative_tolerance >= 0 &&
         (p.absolute_tolerance > 0 || p.relative_tolerance > 0) &&
         p.absolute_tolerance <= 1e-7 && p.relative_tolerance <= 1e-4 &&
         std::isfinite(p.maximum_log_step) && p.maximum_log_step > 0 &&
         p.maximum_log_step <= .04 && std::isfinite(p.maximum_phase_step) &&
         p.maximum_phase_step > 0 && p.maximum_phase_step <= .04 &&
         std::isfinite(p.maximum_constraint_residual) &&
         p.maximum_constraint_residual > 0 &&
         p.maximum_constraint_residual <= 1e-7 && p.maximum_points <= 16 &&
         p.maximum_rhs_per_point <= 1000000 && p.maximum_total_rhs <= 4000000 &&
         p.maximum_scalar_updates_per_point <= 128000000 &&
         p.maximum_total_scalar_updates <= 512000000 &&
         p.maximum_background_queries <= 4000000 &&
         p.maximum_age_quadrature_evaluations <= 200000 &&
         p.maximum_native_bytes <= 16 * 1024 * 1024;
}
struct Budget {
  const MasslessFDTransferPolicy &policy;
  MasslessFDTransferWork &total;
  MasslessFDTransferWork *row = nullptr;
  bool add(std::size_t MasslessFDTransferWork::*member, std::size_t n,
           std::size_t global_cap, std::size_t row_cap) {
    if (total.*member > global_cap || n > global_cap - total.*member ||
        (row && (row->*member > row_cap || n > row_cap - row->*member)))
      return false;
    total.*member += n;
    if (row)
      row->*member += n;
    return true;
  }
  bool rhs() {
    return add(&MasslessFDTransferWork::rhs, 1, policy.maximum_total_rhs,
               policy.maximum_rhs_per_point);
  }
  bool scalar(std::size_t n) {
    return add(&MasslessFDTransferWork::scalar_updates, n,
               policy.maximum_total_scalar_updates,
               policy.maximum_scalar_updates_per_point);
  }
  bool background() {
    return add(&MasslessFDTransferWork::background_queries, 1,
               policy.maximum_background_queries,
               policy.maximum_background_queries);
  }
  bool quadrature() {
    return add(&MasslessFDTransferWork::age_quadrature_evaluations, 1,
               policy.maximum_age_quadrature_evaluations,
               policy.maximum_age_quadrature_evaluations);
  }
  void momentum(std::size_t n) {
    total.momentum_callbacks += n;
    if (row)
      row->momentum_callbacks += n;
  }
};
MasslessFDTransferWork difference(const MasslessFDTransferWork &a,
                                  const MasslessFDTransferWork &b) {
  return {a.rhs - b.rhs, a.scalar_updates - b.scalar_updates,
          a.background_queries - b.background_queries,
          a.age_quadrature_evaluations - b.age_quadrature_evaluations,
          a.momentum_callbacks - b.momentum_callbacks};
}
// Private numerical-input witnesses; never another admitted physical model.
struct Control {
  unsigned kind = 0;
  int sign = 0;
  W radiation_shift = 0, lambda_shift = 0;
  bool round_core = false;
};
Epoch epoch(const ThermalBackground &b, W a, W k, W radiation,
            const Control &control, Budget &budget) {
  Epoch e;
  if (!budget.background()) {
    e.status = S::work_limit;
    return e;
  }
  e = detail::thermal_conformal_epoch(b, a, k, radiation,
                                      budget.policy.thermal);
  budget.momentum(e.momentum_callbacks);
  if (e.status != S::ok)
    return e;
  // Retain the central shared P. Witness perturbations propagate its declared
  // numerical envelope; they are labelled controls and cannot alter source().
  const W a2 = a * a, a4 = a2 * a2, original_p = e.p;
  W changed_p = original_p;
  if (control.kind == 1)
    changed_p += control.radiation_shift * (1 - a4);
  if (control.kind == 2)
    changed_p += control.sign * e.p_error;
  if (!(changed_p > e.p_error)) {
    e.status = S::conditioning_budget_exceeded;
    return e;
  }
  const W lambda = W(*b.omega_lambda()) + control.lambda_shift;
  if (lambda < 0 || !(radiation + control.radiation_shift > 0)) {
    e.status = S::outside_domain;
    return e;
  }
  e.hcal *= std::sqrt(changed_p / original_p);
  e.p = changed_p;
  e.fc *= original_p / changed_p;
  e.fr = (radiation + control.radiation_shift) / changed_p;
  e.fl = lambda * a4 / changed_p;
  e.x2 = (k / e.hcal) * (k / e.hcal);
  e.g = -1 + e.fc / 2 + 2 * e.fl;
  e.closure_defect = e.fc + e.fr + e.fl - 1;
  e.background_identity_defect = e.g - 1 + 1.5L * e.fc + 2 * e.fr;
  if (!(e.x2 > 0) || !std::isnormal(e.x2) || !std::isfinite(e.g))
    e.status = S::outside_domain;
  return e;
}
struct Age {
  S status = S::invalid_input;
  W value = 0, error = 0;
};
Age early_age(const ThermalBackground &b, W start, W radiation,
              W radiation_error, const Control &control) {
  Age out;
  const W r = radiation + control.radiation_shift;
  const W m = b.source().omega_cdm;
  const W lambda = W(*b.omega_lambda()) + control.lambda_shift;
  if (!(r > radiation_error) || lambda < 0)
    return out;
  const W c_h0 = detail::thermal_conformal_c_km_s / b.source().h0_km_s_mpc;
  out.value = 2 * c_h0 * start / (std::sqrt(r + m * start) + std::sqrt(r));
  const W a2 = start * start, a4 = a2 * a2;
  const W lambda_tail = c_h0 * lambda * a4 * start / (10 * r * std::sqrt(r));
  out.error = lambda_tail + out.value * radiation_error / (2 * r) +
              64 * eps_w * out.value;
  if (control.kind == 4)
    out.value += control.sign * out.error;
  out.status = std::isnormal(out.value) && std::isnormal(out.error) &&
                       out.value > out.error
                   ? S::ok
                   : S::outside_domain;
  return out;
}
// Independent endpoint algorithm: adaptive Simpson in a, including the actual
// finite Big-Bang P(0) limit, not elapsed conformal time from the initializer.
struct AgeQuadrature {
  const ThermalBackground &background;
  Budget &budget;
  S status = S::ok;
  W maximum_relative_background_error = 0;
  W f(W a) {
    if (!budget.quadrature() || !budget.background()) {
      status = S::work_limit;
      return 0;
    }
    auto q = background.scaled_expansion(a, budget.policy.thermal);
    budget.momentum(q.callbacks);
    if (q.status != S::ok) {
      status = q.status;
      return 0;
    }
    if (!(q.a4_e2 > q.error_estimate)) {
      status = S::conditioning_budget_exceeded;
      return 0;
    }
    maximum_relative_background_error =
        std::max(maximum_relative_background_error,
                 q.error_estimate / (q.a4_e2 - q.error_estimate));
    return detail::thermal_conformal_c_km_s /
           (background.source().h0_km_s_mpc * std::sqrt(q.a4_e2));
  }
  Age integrate(W lo, W hi, W flo, W fm, W fhi, W whole, W tolerance,
                unsigned depth) {
    Age out;
    if (status != S::ok) {
      out.status = status;
      return out;
    }
    const W mid = (lo + hi) / 2;
    const W left_mid = (lo + mid) / 2, right_mid = (mid + hi) / 2;
    if (!(left_mid > lo && mid > left_mid && right_mid > mid &&
          hi > right_mid)) {
      out.status = status = S::conditioning_budget_exceeded;
      return out;
    }
    const W fl = f(left_mid), fr = f(right_mid);
    if (status != S::ok) {
      out.status = status;
      return out;
    }
    const W left = (mid - lo) * (flo + 4 * fl + fm) / 6;
    const W right = (hi - mid) * (fm + 4 * fr + fhi) / 6;
    const W refined = left + right, empirical = 8 * std::abs(refined - whole);
    if (empirical <= tolerance) {
      out = {S::ok, refined, empirical};
      return out;
    }
    if (depth >= 30) {
      out.status = status = S::conditioning_budget_exceeded;
      return out;
    }
    auto l = integrate(lo, mid, flo, fl, fm, left, tolerance / 2, depth + 1);
    if (l.status != S::ok)
      return l;
    auto r = integrate(mid, hi, fm, fr, fhi, right, tolerance / 2, depth + 1);
    if (r.status != S::ok)
      return r;
    return {S::ok, l.value + r.value, l.error + r.error};
  }
  Age at(W a) {
    const W f0 = f(0), fm = f(a / 2), f1 = f(a);
    if (status != S::ok)
      return {status, 0, 0};
    const W whole = a * (f0 + 4 * fm + f1) / 6;
    auto out =
        integrate(0, a, f0, fm, f1, whole, 1e-7L + 1e-11L * std::abs(whole), 0);
    if (out.status == S::ok)
      out.error += std::abs(out.value) * maximum_relative_background_error +
                   64 * eps_w * std::abs(out.value);
    return out;
  }
};
W constraint(const State &y, const Epoch &e, W *raw = nullptr) {
  const W F = e.fc + 4 * e.fr / 3;
  const std::array<W, 4> terms{e.x2 * y[4], 1.5L * F * y[0], -2 * e.fr * y[1],
                               6 * e.fr * y[2]};
  W sum = 0, scale = 0;
  for (W v : terms) {
    sum += v;
    scale += std::abs(v);
  }
  if (raw)
    *raw = sum;
  return scale == 0 ? (sum == 0 ? 0 : std::numeric_limits<W>::infinity())
                    : std::abs(sum) / scale;
}
bool round_core(State &y, Budget &budget) {
  if (!budget.scalar(6))
    return false;
  for (unsigned j = 0; j < 6; ++j)
    y[j] = static_cast<double>(y[j]);
  return true;
}
W clock_phase_upper(W k, W previous_eta, W next_eta) {
  return k * (std::abs(next_eta - previous_eta) +
              64 * eps_w * (std::abs(next_eta) + std::abs(previous_eta)));
}
struct Run {
  S status = S::invalid_input;
  State y{};
  std::array<W, 3> outputs{};
  MasslessFDInitialWitness initial;
  Epoch endpoint_epoch;
  bool state_initialized = false;
  W state_a = 0, metric_a = 0;
  W maximum_stage_phase = 0, maximum_phase_increment = 0;
  W maximum_constraint = 0, maximum_closure = 0, maximum_background_defect = 0;
  W inherited_age_error = 0;
};
Run evolve(const ThermalBackground &b, W k, W start, W target, W radiation,
           W radiation_error, unsigned L, unsigned divisor,
           const Control &control, Budget &budget) {
  Run out;
  const unsigned n = L + 5;
  auto e = epoch(b, start, k, radiation, control, budget);
  if (e.status != S::ok) {
    out.status = e.status;
    return out;
  }
  auto age = early_age(b, start, radiation, radiation_error, control);
  if (age.status != S::ok) {
    out.status = age.status;
    return out;
  }
  const W r =
      W(b.source().omega_cdm) * start / (radiation + control.radiation_shift);
  const W u = k * (age.value + age.error), F = e.fc + 4 * e.fr / 3;
  if (start < W(4e-20) / 4 || !(start < target / 16) || r > 1e-8L ||
      u > 1e-4L || e.fl > 1e-12L || !(F > 0) || k * age.error > 1e-8L) {
    out.status = S::outside_domain;
    return out;
  }
  // n initialized physical components and five explicit overwrites. Reserved
  // State padding and diagnostic witness copies are not hierarchy updates.
  if (!budget.scalar(n + 5)) {
    out.status = S::work_limit;
    return out;
  }
  std::fill_n(out.y.begin(), n, W(0));
  out.y[4] = phi0;
  out.y[5] = p0 * (e.hcal * age.value) * (e.hcal * age.value) / 15;
  out.y[6] = age.value;
  const W psi = phi0 - 6 * e.fr * out.y[5];
  out.y[0] = -e.x2 * phi0 / (1.5L * F);
  out.y[3] = (out.y[0] + 1.5L * p0) / 3;
  out.state_initialized = true;
  out.endpoint_epoch = e;
  out.state_a = out.metric_a = start;
  auto &w = out.initial;
  w.status = S::ok;
  w.scale_factor = static_cast<double>(start);
  w.eta_mpc = age.value;
  w.eta_error_mpc = age.error;
  w.hcal_mpc_inverse = e.hcal;
  w.leading_phi = phi0;
  w.leading_psi = p0;
  w.leading_delta_c = -1.5L * p0;
  w.leading_delta_r = -2 * p0;
  w.leading_theta = k * k * age.value * p0 / 2;
  w.leading_sigma = p0 * (k * age.value) * (k * age.value) / 15;
  w.projected_phi = phi0;
  w.projected_psi = psi;
  w.projected_vc = out.y[3];
  w.projected_delta = out.y[0];
  w.projected_phi_n = -psi + 1.5L * F * out.y[3];
  w.scaled_shear = out.y[5];
  w.lapse_correction = psi - p0;
  w.velocity_correction = out.y[3] - e.hcal * age.value * p0 / 2;
  w.hamiltonian_before =
      e.x2 * phi0 + 1.5L * F * (-1.5L * p0 + 3 * e.hcal * age.value * p0 / 2);
  constraint(out.y, e, &w.hamiltonian_after);
  w.radiation_fraction = e.fr;
  w.closure_defect = e.closure_defect;
  // Omitted-order scaling guard, not a remainder or forward-error certificate.
  const W omitted = r + u * u + e.fl + radiation_error / radiation +
                    age.error / age.value + 64 * eps_w;
  if (std::abs(w.lapse_correction) > 32 * omitted ||
      std::abs(w.velocity_correction) > 32 * omitted) {
    w.status = S::conditioning_budget_exceeded;
    out.status = S::conditioning_budget_exceeded;
    return out;
  }
  const W initial_eta = age.value;
  out.inherited_age_error = age.error;
  W N = std::log(start), end = std::log(target);
  W active_step = 0;
  auto rhs = [&](W at, const State &y, State &dy) -> S {
    if (!budget.rhs())
      return S::work_limit;
    auto q = epoch(b, std::exp(at), k, radiation, control, budget);
    if (q.status != S::ok)
      return q.status;
    out.maximum_closure =
        std::max(out.maximum_closure, std::abs(q.closure_defect));
    out.maximum_background_defect = std::max(
        out.maximum_background_defect, std::abs(q.background_identity_defect));
    const W inherited = age.error + 64 *
                                        std::numeric_limits<double>::epsilon() *
                                        std::abs(y[6] - initial_eta);
    // Check the fixed complete mesh at every attempted RK stage. Refuse an
    // unsafe stage; do not hide it by retrying on an unrecorded finer mesh.
    if (!(y[6] > inherited) || active_step * 129 > q.hcal * (y[6] - inherited))
      return S::conditioning_budget_exceeded;
    const W stage_phase = active_step * std::sqrt(q.x2);
    out.maximum_stage_phase = std::max(out.maximum_stage_phase, stage_phase);
    if (stage_phase > W(budget.policy.maximum_phase_step))
      return S::conditioning_budget_exceeded;
    if (!budget.scalar(n))
      return S::work_limit;
    const W x = std::sqrt(q.x2), total = q.fc + 4 * q.fr / 3;
    const W lapse = y[4] - 6 * q.fr * y[5];
    dy[0] = -q.x2 * y[3] + 6 * q.fr * y[2];
    dy[1] = q.x2 * y[2];
    dy[2] = q.g * y[2] + (y[0] - y[1]) / 3 - q.x2 * y[5];
    dy[3] = (q.g - 1) * y[3] + lapse;
    dy[4] = -lapse + 1.5L * total * y[3] + 2 * q.fr * y[2];
    dy[5] = 4 * (y[3] + y[2]) / 15 - 3 * y[7] / (10 * x) + 2 * q.g * y[5];
    dy[6] = 1 / q.hcal;
    for (unsigned l = 3; l < L; ++l) {
      const W previous = l == 3 ? 2 * q.x2 * y[5] : y[l + 3];
      dy[l + 4] = x * (l * previous - (l + 1) * y[l + 5]) / (2 * l + 1);
    }
    if (!(q.hcal * y[6] > 0))
      return S::outside_domain;
    // Original MB51 approximate outgoing endpoint, with actual Big-Bang age.
    dy[L + 4] = x * y[L + 3] - (L + 1) * y[L + 4] / (q.hcal * y[6]);
    for (unsigned j = 0; j < n; ++j)
      if (!std::isfinite(dy[j]))
        return S::overflow;
    return S::ok;
  };
  State k1, k2, k3, k4, temporary;
  auto assemble = [&](const State &slope, W h) -> bool {
    if (!budget.scalar(n))
      return false;
    for (unsigned j = 0; j < n; ++j)
      temporary[j] = out.y[j] + h * slope[j];
    return !control.round_core || round_core(temporary, budget);
  };
  while (N < end) {
    const W a = std::exp(N);
    e = epoch(b, a, k, radiation, control, budget);
    if (e.status != S::ok) {
      out.status = e.status;
      return out;
    }
    const W inherited = age.error + 64 *
                                        std::numeric_limits<double>::epsilon() *
                                        std::abs(out.y[6] - initial_eta);
    const W lower = out.y[6] - inherited;
    if (!(lower > 0)) {
      out.status = S::conditioning_budget_exceeded;
      return out;
    }
    // Every L uses this SAME L128-based complete mesh, then half/quarter.
    const W base_step = std::min(
        {W(budget.policy.maximum_log_step),
         W(budget.policy.maximum_phase_step) *
             std::exp(-W(budget.policy.maximum_log_step)) / std::sqrt(e.x2),
         e.hcal * lower / 130}); // margin for the separately checked stages
    const W h = std::min(base_step / divisor, end - N);
    if (!(h > 0) || !(N + h > N)) {
      out.status = S::conditioning_budget_exceeded;
      return out;
    }
    active_step = h;
    S status = rhs(N, out.y, k1);
    if (status == S::ok && !assemble(k1, h / 2))
      status = S::work_limit;
    if (status == S::ok)
      status = rhs(N + h / 2, temporary, k2);
    if (status == S::ok && !assemble(k2, h / 2))
      status = S::work_limit;
    if (status == S::ok)
      status = rhs(N + h / 2, temporary, k3);
    if (status == S::ok && !assemble(k3, h))
      status = S::work_limit;
    if (status == S::ok)
      status = rhs(N + h, temporary, k4);
    if (status != S::ok) {
      out.status = status;
      return out;
    }
    if (!budget.scalar(n)) {
      out.status = S::work_limit;
      return out;
    }
    const W previous_eta = out.y[6];
    for (unsigned j = 0; j < n; ++j)
      out.y[j] += h * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]) / 6;
    // The combined vector already belongs to this new epoch. Preserve that
    // fact before any cast, phase or background admission can fail.
    N = std::min(N + h, end);
    out.state_a = N == end ? target : std::exp(N);
    // Initial-age uncertainty is common to both clocks and cancels from this
    // increment. Stage H uncertainty has matched +/- input witnesses; this
    // subtraction/assembly allowance belongs to wide clock arithmetic.
    const W phase_upper = clock_phase_upper(k, previous_eta, out.y[6]);
    out.maximum_phase_increment =
        std::max(out.maximum_phase_increment, phase_upper);
    if (phase_upper > W(budget.policy.maximum_phase_step)) {
      out.status = S::conditioning_budget_exceeded;
      return out;
    }
    if (control.round_core && !round_core(out.y, budget)) {
      out.status = S::work_limit;
      return out;
    }
    e = epoch(b, out.state_a, k, radiation, control, budget);
    if (e.status != S::ok) {
      out.status = e.status;
      return out;
    }
    out.endpoint_epoch = e;
    out.metric_a = out.state_a;
    out.maximum_constraint =
        std::max(out.maximum_constraint, constraint(out.y, e));
    if (out.maximum_constraint > budget.policy.maximum_constraint_residual) {
      out.status = S::conditioning_budget_exceeded;
      return out;
    }
    out.inherited_age_error = inherited;
  }
  e = epoch(b, target, k, radiation, control, budget);
  if (e.status != S::ok) {
    out.status = e.status;
    return out;
  }
  out.outputs = {out.y[0], out.y[4], out.y[4] - 6 * e.fr * out.y[5]};
  out.status = S::ok;
  for (W v : out.outputs)
    if (!std::isfinite(v))
      out.status = S::overflow;
  return out;
}
void refuse(MasslessFDTransferRow &row, unsigned outputs, S cause) {
  if (outputs & massless_fd_comoving_cdm)
    row.comoving_cdm.status = cause;
  if (outputs & massless_fd_spatial_potential)
    row.spatial_potential.status = cause;
  if (outputs & massless_fd_lapse_potential)
    row.lapse_potential.status = cause;
}
void accept(MasslessFDTransferValue &out, W value, std::array<W, 5> errors,
            const MasslessFDTransferPolicy &policy) {
  const double rounded = static_cast<double>(value);
  const W cast = std::abs(value - W(rounded));
  errors[4] += cast + 64 * eps_w * std::max(1.L, std::abs(value));
  const W epsilon = W(policy.absolute_tolerance) +
                    W(policy.relative_tolerance) * std::abs(value);
  W sum = 0;
  for (W error : errors) {
    if (!(error >= 0) || !std::isfinite(error)) {
      out.status = S::overflow;
      return;
    }
    sum += error;
  }
  auto outward = [](W error) {
    double q = static_cast<double>(error);
    return W(q) < error
               ? std::nextafter(q, std::numeric_limits<double>::infinity())
               : q;
  };
  out.time_refinement = outward(errors[0]);
  out.initial_refinement = outward(errors[1]);
  out.hierarchy_refinement = outward(errors[2]);
  out.background_age_sensitivity = outward(errors[3]);
  out.arithmetic_cast_sensitivity = outward(errors[4]);
  const std::array<W, 5> published{
      W(out.time_refinement), W(out.initial_refinement),
      W(out.hierarchy_refinement), W(out.background_age_sensitivity),
      W(out.arithmetic_cast_sensitivity)};
  W published_sum = 0;
  for (W error : published)
    if (error > 0)
      published_sum = std::nextafter(published_sum + error,
                                     std::numeric_limits<W>::infinity());
  out.absolute_error_estimate = outward(published_sum);
  if (!std::isfinite(rounded) || (value != 0 && !std::isnormal(rounded)) ||
      !std::isnormal(out.absolute_error_estimate)) {
    out.status = S::outside_domain;
    return;
  }
  const W published_epsilon =
      W(policy.absolute_tolerance) +
      W(policy.relative_tolerance) * std::abs(W(rounded));
  if (sum > 5 * epsilon / 6 ||
      std::any_of(errors.begin(), errors.end(),
                  [&](W e) { return e > epsilon / 6; }) ||
      W(out.absolute_error_estimate) >
          5 * std::min(epsilon, published_epsilon) / 6 ||
      std::any_of(published.begin(), published.end(), [&](W e) {
        return e > std::min(epsilon, published_epsilon) / 6;
      })) {
    out.status = S::conditioning_budget_exceeded;
    return;
  }
  out.value = rounded;
  out.status = S::ok;
}
} // namespace

MasslessFDTransfer::MasslessFDTransfer(MasslessFDTransfer &&other) noexcept {
  *this = std::move(other);
}
MasslessFDTransfer &
MasslessFDTransfer::operator=(MasslessFDTransfer &&other) noexcept {
  if (this == &other)
    return *this;
  status_ = other.status_;
  background_ = std::move(other.background_);
  initial_ = other.initial_;
  other.status_ = S::invalid_input;
  other.background_.reset();
  other.initial_ = 0;
  return *this;
}
MasslessFDTransfer prepare_massless_fd_transfer(const ThermalBackground &b,
                                                double initial) {
  MasslessFDTransfer out;
  if (!arithmetic())
    return out;
  if (b.status() != S::ok) {
    out.status_ = b.status();
    return out;
  }
  if (!std::isfinite(initial)) {
    out.status_ = S::nonfinite_input;
    return out;
  }
  const auto &m = b.source();
  if (initial < 4e-20 || initial > 1e-6 || m.omega_gamma != 0 ||
      m.omega_massless_nonphoton != 0 || m.omega_b != 0 || !(m.omega_cdm > 0) ||
      m.species.empty() || m.species.size() > 16) {
    out.status_ = S::outside_domain;
    return out;
  }
  for (const auto &s : m.species)
    if (s.mass_ev != 0 || !(s.temperature_today_ev > 0) ||
        !(s.statistical_weight > 0) || !std::isfinite(s.temperature_today_ev) ||
        !std::isfinite(s.statistical_weight)) {
      out.status_ = S::outside_domain;
      return out;
    }
  out.background_ = b;
  out.initial_ = initial;
  out.status_ = S::ok;
  return out;
}
std::optional<std::size_t>
massless_fd_transfer_payload_bound(std::size_t points,
                                   std::size_t species) noexcept {
  if (points > 16 || species == 0 || species > 16)
    return {};
  irred::detail::PayloadAccounting bytes(sizeof(MasslessFDTransferBatch) +
                                         sizeof(MasslessFDTransfer));
  bytes.add(points, sizeof(MasslessFDTransferRow));
  bytes.add(24,
            sizeof(State));   // anchors, current control, RK stages and copies
  bytes.add(4096, sizeof(W)); // bounded age quadrature/scratch/initial ledger
  bytes.embedded(thermal_background_payload_bound(0, species),
                 sizeof(ThermalBackground));
  return bytes.result();
}
MasslessFDTransferBatch
MasslessFDTransfer::evaluate(std::span<const double> ks, double target,
                             unsigned outputs,
                             MasslessFDTransferPolicy policy) const {
  MasslessFDTransferBatch out;
  out.scale_factor = target;
  out.requested_outputs = outputs;
  if (status_ != S::ok) {
    out.status = status_;
    return out;
  }
  if (!valid_policy(policy) || !outputs || (outputs & ~all_fields) ||
      ks.empty())
    return out;
  if (policy.thermal.momentum_method != background_->momentum_method())
    return out;
  if (!std::isfinite(target)) {
    out.status = S::nonfinite_input;
    return out;
  }
  if (target < 1e-4 || target > 1) {
    out.status = S::outside_domain;
    return out;
  }
  auto bytes = massless_fd_transfer_payload_bound(
      ks.size(), background_->source().species.size());
  if (!bytes || ks.size() > policy.maximum_points ||
      *bytes > policy.maximum_native_bytes) {
    out.status = S::work_limit;
    return out;
  }
  out.rows.reserve(ks.size());
  out.status = S::ok;
  Budget budget{policy, out.work};
  if (!budget.background()) {
    out.status = S::work_limit;
    return out;
  }
  const auto radiation_state = background_->scaled_expansion(0, policy.thermal);
  budget.momentum(radiation_state.callbacks);
  if (radiation_state.status != S::ok) {
    out.status = radiation_state.status;
    return out;
  }
  const W radiation = radiation_state.a4_e2,
          radiation_error = radiation_state.error_estimate;
  if (!(radiation > radiation_error) || !std::isnormal(radiation_error)) {
    out.status = S::conditioning_budget_exceeded;
    return out;
  }
  Age reference_age;
  bool age_acquired = false;
  for (double k_input : ks) {
    out.rows.emplace_back();
    auto &row = out.rows.back();
    row.wavenumber_mpc_inverse = k_input;
    budget.row = &row.work;
    if (!std::isfinite(k_input)) {
      refuse(row, outputs, S::nonfinite_input);
      continue;
    }
    if (k_input < 1e-7 || k_input > 1e-2) {
      refuse(row, outputs, S::outside_domain);
      continue;
    }
    const W k = k_input;
    if (!age_acquired) {
      reference_age = AgeQuadrature{*background_, budget}.at(target);
      age_acquired = true;
    }
    if (reference_age.status != S::ok) {
      refuse(row, outputs, reference_age.status);
      continue;
    }
    row.conformal_age_mpc = reference_age.value;
    row.conformal_age_error_mpc = reference_age.error;
    if (k * reference_age.error > 1e-6L ||
        k * (reference_age.value + reference_age.error) > 20) {
      refuse(row, outputs, S::outside_domain);
      continue;
    }
    auto trial = [&](W start, unsigned L, unsigned divisor,
                     const Control &control = {}) {
      const auto before = row.work;
      auto run = evolve(*background_, k, start, target, radiation,
                        radiation_error, L, divisor, control, budget);
      auto &witness = row.attempts[row.attempts_recorded++];
      witness.status = run.status;
      witness.endpoint_available = run.state_initialized;
      witness.lapse_available =
          run.state_initialized && run.state_a == run.metric_a;
      witness.hierarchy_l = L;
      witness.time_divisor = divisor;
      witness.control_kind = control.kind;
      witness.control_sign = control.sign;
      witness.initial = run.initial;
      if (run.state_initialized) {
        witness.endpoint[0] = run.y[0];
        witness.endpoint[1] = run.y[4];
        if (witness.lapse_available)
          witness.endpoint[2] = run.y[4] - 6 * run.endpoint_epoch.fr * run.y[5];
      }
      witness.endpoint_scale_factor = run.state_a;
      witness.metric_epoch_scale_factor = run.metric_a;
      witness.endpoint_scaled_shear = run.y[5];
      witness.endpoint_eta_mpc = run.y[6];
      witness.maximum_stage_phase_bound = run.maximum_stage_phase;
      witness.maximum_phase_increment_upper = run.maximum_phase_increment;
      witness.maximum_constraint_residual =
          static_cast<double>(run.maximum_constraint);
      witness.work = difference(row.work, before);
      row.maximum_constraint_residual = std::max(
          row.maximum_constraint_residual, witness.maximum_constraint_residual);
      row.maximum_closure_defect = std::max(
          row.maximum_closure_defect, static_cast<double>(run.maximum_closure));
      row.maximum_background_identity_defect =
          std::max(row.maximum_background_identity_defect,
                   static_cast<double>(run.maximum_background_defect));
      return run;
    };
    // Finest, original/half starts, coarse/half time, L32/L64, joint coarse.
    const std::array<W, 8> starts{
        W(initial_) / 4, W(initial_),     W(initial_) / 2, W(initial_) / 4,
        W(initial_) / 4, W(initial_) / 4, W(initial_) / 4, W(initial_)};
    const std::array<unsigned, 8> levels{128, 128, 128, 128, 128, 32, 64, 32};
    const std::array<unsigned, 8> divisors{4, 4, 4, 1, 2, 4, 4, 1};
    std::array<Run, 8> runs;
    S cause = S::ok;
    for (unsigned j = 0; j < 8; ++j) {
      runs[j] = trial(starts[j], levels[j], divisors[j]);
      if (j == 0)
        row.initial_states[2] = runs[j].initial;
      if (j == 1 || j == 2)
        row.initial_states[j - 1] = runs[j].initial;
      if (runs[j].status != S::ok) {
        cause = runs[j].status;
        break;
      }
    }
    if (cause != S::ok) {
      refuse(row, outputs, cause);
      continue;
    }
    const W age_time = 8 * std::max(std::abs(runs[0].y[6] - runs[4].y[6]),
                                    std::abs(runs[4].y[6] - runs[3].y[6]));
    const W age_difference = std::abs(runs[0].y[6] - reference_age.value);
    row.conformal_age_mpc = runs[0].y[6];
    row.conformal_age_error_mpc =
        reference_age.error + age_time + age_difference;
    if (k * row.conformal_age_error_mpc > 1e-6L ||
        k * (row.conformal_age_mpc + row.conformal_age_error_mpc) > 20) {
      refuse(row, outputs, S::conditioning_budget_exceeded);
      continue;
    }
    std::array<std::array<W, 5>, 3> errors{};
    for (unsigned field = 0; field < 3; ++field) {
      auto &e = errors[field];
      auto v = [&](unsigned j) { return runs[j].outputs[field]; };
      e[0] = 8 * std::max(std::abs(v(0) - v(4)), std::abs(v(4) - v(3)));
      e[1] = 8 * std::max(std::abs(v(1) - v(2)), std::abs(v(2) - v(0)));
      e[2] = 8 * std::max(std::abs(v(5) - v(6)), std::abs(v(6) - v(0)));
      if (std::abs(v(7) - v(0)) >
          e[0] + e[1] + e[2] + 128 * eps_w * std::max(1.L, std::abs(v(0))))
        cause = S::conditioning_budget_exceeded;
    }
    if (cause != S::ok) {
      refuse(row, outputs, cause);
      continue;
    }
    // Every original input axis gets matched FINEST +/- witnesses. Their
    // sufficiency under the original caps must be earned, never presumed.
    const double lambda = *background_->omega_lambda();
    const W lambda_cast =
        lambda == 0
            ? 0
            : std::max(W(std::nextafter(
                           lambda, std::numeric_limits<double>::infinity())) -
                           lambda,
                       W(lambda) - std::nextafter(lambda, 0.0)) /
                  2;
    for (unsigned kind = 1; kind <= 4 && cause == S::ok; ++kind) {
      for (int sign : {-1, 1}) {
        Control control;
        control.kind = kind;
        control.sign = sign;
        if (kind == 1) {
          control.radiation_shift = sign * radiation_error;
          control.lambda_shift = -control.radiation_shift;
        }
        if (kind == 3)
          control.lambda_shift = sign * lambda_cast;
        auto run = trial(W(initial_) / 4, 128, 4, control);
        if (run.status != S::ok) {
          cause = run.status;
          break;
        }
        for (unsigned field = 0; field < 3; ++field)
          errors[field][3] +=
              8 * std::abs(run.outputs[field] - runs[0].outputs[field]);
      }
    }
    if (cause == S::ok) {
      Control reduced;
      reduced.kind = 5;
      reduced.round_core = true;
      auto run = trial(W(initial_) / 4, 128, 4, reduced);
      cause = run.status;
      if (cause == S::ok)
        for (unsigned field = 0; field < 3; ++field)
          errors[field][4] =
              8 * std::abs(run.outputs[field] - runs[0].outputs[field]);
    }
    if (cause != S::ok) {
      refuse(row, outputs, cause);
      continue;
    }
    if (outputs & massless_fd_comoving_cdm)
      accept(row.comoving_cdm, runs[0].outputs[0], errors[0], policy);
    if (outputs & massless_fd_spatial_potential)
      accept(row.spatial_potential, runs[0].outputs[1], errors[1], policy);
    if (outputs & massless_fd_lapse_potential)
      accept(row.lapse_potential, runs[0].outputs[2], errors[2], policy);
  }
  return out;
}
} // namespace irred::cosmology
