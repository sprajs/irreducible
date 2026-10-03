#pragma once

// Test-only MB48 angular transport / Einstein TRACE reference. The retained
// thermal epoch is shared ancestry; metric, state, time and angular algorithms
// are independent of the native momentum/hierarchy implementation. Equations
// are implemented from the accepted source-only v3 derivation, not copied code.
#include "../src/thermal_conformal_epoch.hpp"
#include "irred/massless_fd_transfer.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace massless_fd_peer {
using Status = irred::numerics::Status;
using Background = irred::cosmology::ThermalBackground;
using Epoch = irred::cosmology::detail::ThermalConformalEpoch;
inline constexpr long double p = -10.L / 19;
inline constexpr long double phi_star = 7 * p / 5, s_star = 2 * p / 5;
inline constexpr long double b_star = -2 * p;
inline constexpr const char *method_id =
    "test-MB48-centered-actual-Gram-angular-eta-DP54-original-TRACE-kappa2/v3";

// All controls of one requested k share this budget; a failed/rejected attempt
// does not release work. Ledger bytes also accumulate across its trials.
struct Budget {
  Status status = Status::ok;
  std::size_t rhs = 0, scalar_updates = 0, background_queries = 0;
  std::size_t age_quadrature_evaluations = 0, momentum_callbacks = 0;
  std::size_t ledger_bytes = 0, attempted_stages = 0, logging_refusals = 0;
  std::size_t maximum_rhs = 2000000, maximum_scalar_updates = 512000000;
  std::size_t maximum_background_queries = 4000000;
  std::size_t maximum_age_quadrature_evaluations = 200000;
  std::size_t maximum_payload_bytes = 16 * 1024 * 1024;
  std::size_t maximum_ledger_bytes = 256 * 1024 * 1024;
  std::size_t *campaign_bytes = nullptr;
  std::size_t maximum_campaign_bytes = 2ULL * 1024 * 1024 * 1024;

  bool charge(std::size_t &used, std::size_t amount, std::size_t limit) {
    if (status != Status::ok)
      return false;
    if (amount > limit || used > limit - amount) {
      status = Status::work_limit;
      return false;
    }
    used += amount;
    return true;
  }
  bool scalars(std::size_t amount) {
    return charge(scalar_updates, amount, maximum_scalar_updates);
  }
  bool stage() {
    ++attempted_stages;
    return charge(rhs, 1, maximum_rhs);
  }
};

template <class T> using C = std::complex<T>;
template <class T> T P2(T mu) { return (3 * mu * mu - 1) / 2; }
template <class T> bool finite(C<T> z) {
  return std::isfinite(z.real()) && std::isfinite(z.imag());
}

template <class T> struct Rule {
  Status status = Status::invalid_input;
  std::vector<T> mu, weight;
  T normalization = 0, m2 = 0, m4 = 0, gram_defect = 0;
  T odd_moment = 0, squared_p2 = 0;
};

// A normalized half-mean is the ratio of two sums of the actual stored
// weights. Thus no rounded node or Gram moment is replaced by 1/3 or 1/5.
template <class T, class F>
auto mean(const Rule<T> &r, F f, Budget *budget = nullptr) {
  using V = decltype(f(T{}));
  V sum{}, correction{};
  if (budget &&
      !budget->scalars(r.mu.size() * (std::is_same_v<V, C<T>> ? 2 : 1)))
    return V{};
  for (std::size_t j = 0; j < r.mu.size(); ++j) {
    const V v = r.weight[j] * f(r.mu[j]) - correction;
    const V next = sum + v;
    correction = (next - sum) - v;
    sum = next;
  }
  return sum / r.normalization;
}

template <class T> void inspect_rule(Rule<T> &r) {
  r.status = Status::invalid_input;
  if (r.mu.empty() || r.mu.size() != r.weight.size() || r.mu.size() > 256)
    return;
  r.normalization = 0;
  for (std::size_t j = 0; j < r.mu.size(); ++j) {
    if (!std::isfinite(r.mu[j]) || !std::isfinite(r.weight[j]) ||
        !(r.weight[j] > 0) || std::abs(r.mu[j]) > 1)
      return;
    const auto other = r.mu.size() - 1 - j;
    if (r.mu[j] != -r.mu[other] || r.weight[j] != r.weight[other])
      return;
    r.normalization += r.weight[j];
  }
  if (!(r.normalization > 0) || !std::isfinite(r.normalization))
    return;
  r.m2 = mean(r, [](T u) { return u * u; });
  r.m4 = mean(r, [](T u) { return u * u * u * u; });
  r.odd_moment = mean(r, [](T u) { return u; });
  r.squared_p2 = mean(r, [](T u) { return P2(u) * P2(u); });
  const T determinant = r.m2 - r.odd_moment * r.odd_moment;
  if (!(determinant > T(64) * std::numeric_limits<T>::epsilon()))
    return;
  r.gram_defect = std::max({std::abs(r.m2 - T(1) / 3),
                            std::abs(r.m4 - T(1) / 5), std::abs(r.odd_moment)});
  r.status = Status::ok;
}

template <class T> Rule<T> make_rule(unsigned nodes, Budget &budget) {
  Rule<T> r;
  if (nodes != 64 && nodes != 128 && nodes != 256)
    return r;
  if (sizeof(T) * nodes * 4 > budget.maximum_payload_bytes) {
    budget.status = Status::work_limit;
    return r;
  }
  r.mu.resize(nodes);
  r.weight.resize(nodes);
  const T pi = std::acos(T(-1));
  for (unsigned i = 0; i < nodes / 2; ++i) {
    T z = std::cos(pi * (T(i) + T(.75)) / (T(nodes) + T(.5)));
    T derivative = 0;
    bool converged = false;
    for (unsigned iteration = 0; iteration < 64; ++iteration) {
      if (!budget.scalars(2 * nodes + 4))
        return r;
      T l0 = 1, l1 = z;
      for (unsigned l = 2; l <= nodes; ++l) {
        const T next = ((2 * T(l) - 1) * z * l1 - (T(l) - 1) * l0) / T(l);
        l0 = l1;
        l1 = next;
      }
      derivative = T(nodes) * (z * l1 - l0) / (z * z - 1);
      const T next = z - l1 / derivative;
      if (std::abs(next - z) <= 4 * std::numeric_limits<T>::epsilon()) {
        z = next;
        converged = true;
        break;
      }
      z = next;
    }
    if (!converged) {
      budget.status = Status::conditioning_budget_exceeded;
      return r;
    }
    // Re-evaluate derivative at the stored node before forming its weight.
    if (!budget.scalars(2 * nodes + 4))
      return r;
    T l0 = 1, l1 = z;
    for (unsigned l = 2; l <= nodes; ++l) {
      const T next = ((2 * T(l) - 1) * z * l1 - (T(l) - 1) * l0) / T(l);
      l0 = l1;
      l1 = next;
    }
    derivative = T(nodes) * (z * l1 - l0) / (z * z - 1);
    const T weight = 1 / ((1 - z * z) * derivative * derivative);
    r.mu[i] = -z;
    r.mu[nodes - 1 - i] = z;
    r.weight[i] = r.weight[nodes - 1 - i] = weight;
  }
  if (!budget.scalars(5 * nodes + 7))
    return r;
  inspect_rule(r);
  if (r.status != Status::ok)
    budget.status = r.status;
  return r;
}

template <class T> struct Angular {
  T A = 0;
  C<T> B{};
  std::vector<C<T>> Q;
};
template <class T> struct Translation {
  C<T> alpha{}, beta{};
  T maximum_node_change = 0, shear_translation = 0, parity_defect = 0;
  std::size_t nodes_translated = 0;
};

template <class T>
Angular<T> angular_rhs(const Rule<T> &r, const Angular<T> &y, T k,
                       T metric_source, Budget *budget = nullptr) {
  Angular<T> d;
  d.Q.resize(y.Q.size());
  C<T> m{};
  for (std::size_t j = 0; j < y.Q.size(); ++j) {
    if (budget && !budget->scalars(2))
      return d;
    m += r.weight[j] * r.mu[j] * r.mu[j] * y.Q[j];
  }
  m /= r.normalization;
  d.A = metric_source + k * k * r.m2 * y.B.real();
  d.B = -y.A - k * k * m / r.m2;
  for (std::size_t j = 0; j < y.Q.size(); ++j) {
    if (budget && !budget->scalars(2))
      return d;
    const T u = r.mu[j];
    d.Q[j] = y.B * (u * u - r.m2) - C<T>(0, k) * (u * y.Q[j] - u * m / r.m2);
  }
  return d;
}

template <class T>
Translation<T> recenter(const Rule<T> &r, Angular<T> &y, T k,
                        Budget *budget = nullptr) {
  Translation<T> out;
  C<T> q0{}, q1{};
  for (std::size_t j = 0; j < y.Q.size(); ++j) {
    if (budget && !budget->scalars(4))
      return out;
    q0 += r.weight[j] * y.Q[j];
    q1 += r.weight[j] * r.mu[j] * y.Q[j];
  }
  if (budget && !budget->scalars(9))
    return out;
  q0 /= r.normalization;
  q1 /= r.normalization;
  const T determinant = r.m2 - r.odd_moment * r.odd_moment;
  out.alpha = (r.m2 * q0 - r.odd_moment * q1) / determinant;
  out.beta = (q1 - r.odd_moment * q0) / determinant;
  out.parity_defect =
      std::max(std::abs(out.alpha.imag()), std::abs(out.beta.real()));
  if (budget && out.parity_defect > T(1e-7)) {
    budget->status = Status::conditioning_budget_exceeded;
    return out;
  }
  const T old_A = y.A;
  const C<T> old_B = y.B;
  if (budget && !budget->scalars(3))
    return out;
  y.A += (k * k * out.alpha).real();
  y.B -= C<T>(0, k) * out.beta;
  for (std::size_t j = 0; j < y.Q.size(); ++j) {
    if (budget && !budget->scalars(6))
      return out;
    const T u = r.mu[j];
    const C<T> old_D = old_A + C<T>(0, k * u) * old_B + k * k * y.Q[j];
    y.Q[j] -= out.alpha + out.beta * u;
    ++out.nodes_translated;
    const C<T> new_D = y.A + C<T>(0, k * u) * y.B + k * k * y.Q[j];
    out.maximum_node_change =
        std::max(out.maximum_node_change, std::abs(new_D - old_D));
  }
  const T odd_p2 = mean(r, [](T u) { return u * P2(u); }, budget);
  out.shear_translation =
      (k * k / 2 * (out.alpha * (3 * r.m2 - 1) / T(2) + out.beta * odd_p2))
          .real();
  return out;
}

template <class T> struct Hybrid {
  T delta = 0, theta = 0, Vr = 0, sigma = 0, psi = 0;
};
template <class T>
Hybrid<T> hybrid(const Rule<T> &r, const Angular<T> &y, T k, T H, T phi, T fr,
                 Budget *budget = nullptr) {
  C<T> q2{};
  for (std::size_t j = 0; j < y.Q.size(); ++j) {
    if (budget && !budget->scalars(2))
      return {};
    q2 += r.weight[j] * P2(r.mu[j]) * y.Q[j];
  }
  q2 /= r.normalization;
  Hybrid<T> out;
  out.sigma = -k * k * q2.real() / 2;
  out.psi = phi + 3 * fr * H * H * q2.real();
  out.delta = y.A - 4 * out.psi;
  out.theta = -k * k * y.B.real() / 4;
  out.Vr = -H * y.B.real() / 4;
  return out;
}

template <class T> struct Ideal {
  T phi = T(phi_star), psi = 0, delta_c = -3 * T(p) / 2;
  T delta_r = -2 * T(p), Vc = 0, Vr = 0, sigma = 0, Delta = 0;
  T Z = 0, eta = 0, eta_error = 0, H = 0, fr = 0, fc = 0;
};
template <class T> struct Initial {
  Ideal<T> ideal;
  Angular<T> angular;
  Hybrid<T> actual;
  Translation<T> translation;
  T Z = 0, E = 0, W = 0;
};
template <class T>
Initial<T> initialize(const Rule<T> &r, const Epoch &epoch, T k, T eta,
                      T eta_error, Budget *budget = nullptr) {
  Initial<T> out;
  auto &i = out.ideal;
  i.H = T(epoch.hcal);
  i.fr = T(epoch.fr);
  i.fc = T(epoch.fc);
  i.eta = eta;
  i.eta_error = eta_error;
  const T x2 = T(epoch.x2), F = i.fc + 4 * i.fr / 3;
  i.sigma = T(p) * k * k * eta * eta / 15;
  i.psi = i.phi - 6 * i.fr * T(p) * (i.H * eta) * (i.H * eta) / 15;
  i.Delta = -x2 * i.phi / (T(1.5) * F);
  i.Vc = i.Vr = T(p) / 2 + i.Delta / 3;
  i.Z = -i.psi + T(1.5) * F * i.Vc;
  out.Z = i.Z;
  out.angular.A = i.delta_r + 4 * i.psi;
  out.angular.B = -4 * i.Vr / i.H;
  out.angular.Q.resize(r.mu.size());
  for (std::size_t j = 0; j < r.mu.size(); ++j) {
    if (budget && !budget->scalars(2))
      return out;
    // Stable equivalent of -10 sigma P2/k^2, including tiny k.
    out.angular.Q[j] = -2 * T(p) * eta * eta * P2(r.mu[j]) / 3;
  }
  out.translation = recenter(r, out.angular, k, budget);
  out.actual = hybrid(r, out.angular, k, i.H, i.phi, i.fr, budget);
  out.E = out.Z + out.actual.psi + x2 * i.phi / 3 +
          (i.fc * i.delta_c + i.fr * out.actual.delta) / 2;
  out.W = out.Z + out.actual.psi -
          T(1.5) * (i.fc * i.Vc + 4 * i.fr * out.actual.Vr / 3);
  return out;
}

template <class T>
std::array<T, 2> constraint_rhs(T g, T x, T kappa, T E, T W, T fE = 0,
                                T fW = 0) {
  return {-(2 * g + 1 + kappa) * E + x * x * W / 3 + fE,
          -kappa * E - (g + 2) * W + fW};
}

// The logger retains every attempted RHS, including rejected DP stages. If
// writing is refused, its bounded last refusal/counters remain in Result;
// the existing evidence file is retained without truncation.
struct Ledger {
  std::ofstream stream;
  Status status = Status::ok;
  explicit Ledger(const std::string &path) {
    if (!path.empty()) {
      stream.open(path, std::ios::out | std::ios::app);
      if (!stream || stream.tellp() != std::streampos(0))
        status = Status::invalid_input;
    } else {
      status = Status::invalid_input;
    }
  }
  bool write(const std::string &record, Budget &budget) {
    if (status != Status::ok || record.size() > 8192 ||
        record.size() > budget.maximum_ledger_bytes ||
        budget.ledger_bytes > budget.maximum_ledger_bytes - record.size() ||
        (budget.campaign_bytes &&
         (record.size() > budget.maximum_campaign_bytes ||
          *budget.campaign_bytes >
              budget.maximum_campaign_bytes - record.size()))) {
      ++budget.logging_refusals;
      status = budget.status = Status::work_limit;
      return false;
    }
    stream.write(record.data(), static_cast<std::streamsize>(record.size()));
    stream.flush();
    if (!stream) {
      ++budget.logging_refusals;
      status = budget.status = Status::invalid_input;
      return false;
    }
    budget.ledger_bytes += record.size();
    if (budget.campaign_bytes)
      *budget.campaign_bytes += record.size();
    return true;
  }
};

struct Controls {
  unsigned nodes = 256;
  long double maximum_log_step = .01L, maximum_phase_step = .01L;
  long double absolute_state_tolerance = 2e-12L / 16;
  long double relative_state_tolerance = 2e-11L / 16;
  long double kappa = 2;
  int background_direction = 0, age_direction = 0;
  // Uncertainty directions are distinct controls, never replacement H models.
  // 0=all jointly, 1=P only, 2=retained R normalization, 3=Lambda accessor.
  unsigned background_component = 0;
  bool finite_recentering = true;
  std::string run_id = "main";
};
struct Diagnostics {
  long double maximum_E = 0, maximum_W = 0, maximum_trace = 0;
  long double maximum_slip = 0, maximum_conservation = 0;
  long double maximum_correction = 0, maximum_recenter = 0;
  long double maximum_closure = 0, maximum_background_identity = 0;
  long double integrated_forcing = 0, maximum_subspace_defect = 0;
  long double initial_E = 0, initial_W = 0, initial_zeta_translation = 0;
};
struct Result {
  Status status = Status::invalid_input;
  std::array<long double, 3> value{}; // Delta_c,phi,psi
  // time, angle, initial, background/age, arithmetic/translation
  std::array<std::array<long double, 5>, 3> contributions{};
  std::array<long double, 3> error{};
  Diagnostics diagnostics;
  long double eta = 0, eta_error = 0, endpoint_log_error = 0;
  std::size_t rhs = 0, scalar_updates = 0, background_queries = 0;
  std::size_t age_quadrature_evaluations = 0, ledger_bytes = 0;
};

inline Epoch shared_epoch(const Background &b, long double a, long double k,
                          long double radiation, long double radiation_error,
                          int direction, Budget &budget,
                          unsigned component = 0) {
  Epoch e;
  if (!budget.charge(budget.background_queries, 1,
                     budget.maximum_background_queries)) {
    e.status = budget.status;
    return e;
  }
  irred::cosmology::ThermalPolicy policy;
  policy.momentum_method = b.momentum_method();
  e = irred::cosmology::detail::thermal_conformal_epoch(b, a, k, radiation,
                                                        policy);
  budget.momentum_callbacks += e.momentum_callbacks;
  if (e.status != Status::ok)
    return e;
  if (direction) {
    // A diagnostic uncertainty control of the SAME retained P, explicitly
    // separate from the central shared-H run, never a replacement background.
    const long double controlled_p =
        e.p + ((component == 0 || component == 1) ? direction * e.p_error : 0);
    if (!(controlled_p > 0)) {
      e.status = Status::outside_domain;
      return e;
    }
    const long double ratio = e.p / controlled_p;
    e.hcal *= std::sqrt(controlled_p / e.p);
    e.fc *= ratio;
    e.fr = (radiation + ((component == 0 || component == 2)
                             ? direction * radiation_error
                             : 0)) /
           controlled_p;
    e.fl *= ratio;
    if (component == 0 || component == 3)
      e.fl += direction * e.lambda_cast_error * a * a * a * a / controlled_p;
    e.p = controlled_p;
    e.x2 = (k / e.hcal) * (k / e.hcal);
    e.g = -1 + e.fc / 2 + 2 * e.fl;
    e.closure_defect = e.fc + e.fr + e.fl - 1;
    e.background_identity_defect = e.g - 1 + 1.5L * e.fc + 2 * e.fr;
    if (!(e.fr > 0) || !(e.fc > 0) || !(e.fl >= 0) || !std::isfinite(e.hcal) ||
        !std::isfinite(e.g))
      e.status = Status::conditioning_budget_exceeded;
  }
  return e;
}

struct Age {
  Status status = Status::invalid_input;
  long double value = 0, error = 0;
};
inline Age initial_age(const Background &b, long double a,
                       long double radiation, long double radiation_error,
                       Budget &budget) {
  Age out;
  if (!budget.charge(budget.age_quadrature_evaluations, 1,
                     budget.maximum_age_quadrature_evaluations)) {
    out.status = budget.status;
    return out;
  }
  const auto e = shared_epoch(b, a, 0, radiation, radiation_error, 0, budget);
  if (e.status != Status::ok) {
    out.status = e.status;
    return out;
  }
  const long double c_over_h =
      static_cast<long double>(irred::speed_of_light_m_per_s) /
      (1000 * b.source().h0_km_s_mpc);
  const long double matter = b.source().omega_cdm;
  const long double lambda_upper =
      static_cast<long double>(*b.omega_lambda()) + e.lambda_cast_error;
  out.value = c_over_h * 2 * a /
              (std::sqrt(radiation + matter * a) + std::sqrt(radiation));
  const long double lambda_error = c_over_h * lambda_upper * std::pow(a, 5) /
                                   (10 * std::pow(radiation, 1.5L));
  out.error = lambda_error + out.value * radiation_error / (2 * radiation) +
              64 * std::numeric_limits<long double>::epsilon() * out.value;
  out.status =
      std::isfinite(out.value) && out.value > out.error && e.fl <= 1e-12L
          ? Status::ok
          : Status::outside_domain;
  return out;
}

// Independent endpoint age in log-a with adaptive Simpson refinement, sharing
// only retained H. The early analytic radiation/matter primitive has an
// explicit positive-Lambda ceiling; it is not an alternative H in the ODE.
inline Age endpoint_age(const Background &b, long double target,
                        long double radiation, long double radiation_error,
                        Budget &budget) {
  const long double early = std::min(
      {target, 1e-8L,
       1e-12L * radiation / static_cast<long double>(b.source().omega_cdm)});
  Age out = initial_age(b, early, radiation, radiation_error, budget);
  if (out.status != Status::ok || early == target)
    return out;
  Status status = Status::ok;
  long double maximum_background_ratio = 0;
  auto integrand = [&](long double n) {
    if (!budget.charge(budget.age_quadrature_evaluations, 1,
                       budget.maximum_age_quadrature_evaluations)) {
      status = budget.status;
      return 0.L;
    }
    const auto e =
        shared_epoch(b, std::exp(n), 0, radiation, radiation_error, 0, budget);
    if (e.status != Status::ok) {
      status = e.status;
      return 0.L;
    }
    maximum_background_ratio =
        std::max(maximum_background_ratio, e.p_error / (2 * e.p));
    return 1 / e.hcal;
  };
  const long double lo = std::log(early), hi = std::log(target);
  const long double left = integrand(lo), middle = integrand((lo + hi) / 2),
                    right = integrand(hi);
  if (status != Status::ok) {
    out.status = status;
    return out;
  }
  const long double initial = (hi - lo) * (left + 4 * middle + right) / 6;
  const long double tolerance = 2e-12L + 2e-12L * std::abs(initial);
  long double error = 0;
  auto refine = [&](auto &&self, long double a, long double z, long double fa,
                    long double fm, long double fz, long double old,
                    long double allocation, unsigned depth) -> long double {
    const long double center = (a + z) / 2;
    const long double f1 = integrand((a + center) / 2);
    const long double f2 = integrand((center + z) / 2);
    if (status != Status::ok)
      return 0;
    const long double first = (center - a) * (fa + 4 * f1 + fm) / 6;
    const long double second = (z - center) * (fm + 4 * f2 + fz) / 6;
    const long double change = first + second - old;
    if (!budget.scalars(12)) {
      status = budget.status;
      return 0;
    }
    if (std::abs(change) / 15 <= allocation) {
      error += 8 * std::abs(change) / 15;
      return first + second + change / 15;
    }
    if (!depth) {
      status = Status::conditioning_budget_exceeded;
      return 0;
    }
    return self(self, a, center, fa, f1, fm, first, allocation / 2, depth - 1) +
           self(self, center, z, fm, f2, fz, second, allocation / 2, depth - 1);
  };
  out.value +=
      refine(refine, lo, hi, left, middle, right, initial, tolerance, 30);
  out.error += error + maximum_background_ratio * std::abs(out.value) +
               128 * std::numeric_limits<long double>::epsilon() *
                   budget.age_quadrature_evaluations * std::abs(out.value);
  out.status = status;
  return out;
}

template <class T> struct State {
  // N,h,Z,w,e,c,b-b*,U-U*. Only N is O(log a); residual coordinates are small.
  std::array<T, 7> v{};
  std::vector<C<T>> R;
};
// Charge actual scalar destinations, including baseline copies. Capacity and
// allocator padding are not scalar updates. A refusal precedes every write.
template <class T> std::size_t active_scalars(const State<T> &y) {
  return y.v.size() + 2 * y.R.size();
}
template <class T>
bool copy_state(State<T> &to, const State<T> &from, Budget &budget) {
  if (!budget.scalars(active_scalars(from)))
    return false;
  to = from;
  return true;
}
template <class T>
bool resize_residual(State<T> &y, std::size_t nodes, Budget &budget) {
  const auto added = nodes > y.R.size() ? nodes - y.R.size() : 0;
  if (!budget.scalars(2 * added))
    return false;
  y.R.resize(nodes);
  return true;
}
template <class T, class F>
bool assign_initial_scalars(State<T> &y, Budget &budget, F write) {
  if (!budget.scalars(y.v.size()))
    return false;
  write(y.v);
  return true;
}
template <class T> T centered_cdm(const State<T> &y) {
  return 3 * y.v[3] + y.v[4] + 3 * y.v[1];
}
template <class T> struct Readout {
  T phi = 0, psi = 0, s = 0, v = 0, b = 0, delta_r = 0;
  T E = 0, W = 0, q2 = 0;
  C<T> q0{}, q1{};
  T original_trace = 0, modified_trace = 0, slip_residual = 0;
  T radiation_conservation = 0, radiation_euler = 0;
  T E_normalized = 0, W_normalized = 0, trace_normalized = 0;
  T slip_normalized = 0, conservation_normalized = 0;
  T parity_defect = 0, quadrupole_imaginary = 0;
};
template <class T>
Readout<T> readout(const State<T> &y, const Rule<T> &r, const Epoch &e,
                   Budget *budget = nullptr) {
  Readout<T> o;
  C<T> q2R{}, q0{}, q1{};
  for (std::size_t j = 0; j < r.mu.size(); ++j) {
    if (budget && !budget->scalars(8))
      return o;
    q2R += r.weight[j] * P2(r.mu[j]) * y.R[j];
    const C<T> U = -2 * T(p) * P2(r.mu[j]) / 3 + y.R[j];
    q0 += r.weight[j] * U;
    q1 += r.weight[j] * r.mu[j] * U;
    o.parity_defect =
        std::max(o.parity_defect,
                 std::abs(y.R[j] - std::conj(y.R[r.mu.size() - 1 - j])));
  }
  q2R /= r.normalization;
  o.quadrupole_imaginary = std::abs(q2R.imag());
  o.q0 = q0 / r.normalization;
  o.q1 = q1 / r.normalization;
  const T q2_correction = -2 * T(p) / 3 * (r.squared_p2 - T(1) / 5);
  o.v = T(s_star) * T(e.closure_defect - e.fc - e.fl) -
        3 * T(e.fr) * (q2R.real() + q2_correction);
  o.s = T(s_star) + o.v;
  o.phi = T(phi_star) + y.v[1];
  o.psi = T(p) + y.v[1] - o.v;
  o.b = T(b_star) + y.v[6];
  o.delta_r = -2 * T(p) + 4 * y.v[1] + y.v[5];
  o.q2 = -2 * T(p) / 15 + q2_correction + q2R.real();
  o.slip_residual = o.psi - o.phi - 3 * T(e.fr) * o.q2;
  o.E = y.v[2] + (1 + T(1.5) * T(e.fc) + 2 * T(e.fr)) * y.v[1] - o.v +
        (T(e.fc) * y.v[4] + T(e.fr) * y.v[5]) / 2 +
        T(p) * T(e.fc / 4 + e.fl - e.closure_defect) + T(e.x2) * o.phi / 3;
  o.W = y.v[2] + o.psi -
        T(1.5) * (T(e.fc) * (T(p) / 2 + y.v[3]) - T(e.fr) * o.b / 3);
  auto normalize = [](T numerator, T denominator) {
    return denominator > 0
               ? std::abs(numerator) / denominator
               : (numerator == 0 ? T(0) : std::numeric_limits<T>::infinity());
  };
  const T delta_c = y.v[4] + 3 * o.phi + 3;
  o.E_normalized = normalize(o.E, std::abs(y.v[2]) + std::abs(o.psi) +
                                      std::abs(T(e.x2) * o.phi / 3) +
                                      std::abs(T(e.fc) * delta_c / 2) +
                                      std::abs(T(e.fr) * o.delta_r / 2));
  o.W_normalized =
      normalize(o.W, std::abs(y.v[2]) + std::abs(o.psi) +
                         std::abs(T(1.5) * T(e.fc) * (T(p) / 2 + y.v[3])) +
                         std::abs(T(e.fr) * o.b / 2));
  o.slip_normalized =
      normalize(o.slip_residual, std::abs(o.psi) + std::abs(o.phi) +
                                     std::abs(3 * T(e.fr) * o.q2));
  return o;
}

template <class T>
bool record_stage(Ledger &ledger, Budget &budget, const Controls &control,
                  long double eta, const State<T> &y, const Epoch &e,
                  const Readout<T> &o, long double sN, long double forcingE,
                  long double forcingW, Status stage_status,
                  const char *kind = "rhs", const Rule<T> *rule = nullptr) {
  std::ostringstream text;
  text << std::setprecision(std::numeric_limits<long double>::max_digits10)
       << "method=" << method_id << " arithmetic="
       << (std::numeric_limits<T>::digits > 53 ? "wide" : "binary64")
       << " run=" << control.run_id << " stage=" << budget.attempted_stages
       << " kind=" << kind << " status=" << static_cast<int>(stage_status)
       << " diagnostics_valid=" << (stage_status == Status::ok ? 1 : 0)
       << " eta=" << eta << " N=" << y.v[0] << " E=" << o.E << " W=" << o.W
       << " original_trace=" << o.original_trace
       << " modified_trace=" << o.modified_trace
       << " correction=" << -control.kappa * o.E << " sN=" << sN
       << " forcing_units=d/dN fE_quadrature=" << forcingE
       << " fW_quadrature=" << forcingW
       << " fE_Rbg=" << e.background_identity_defect * y.v[2]
       << " fW_Rbg=" << -e.background_identity_defect * o.psi
       << " complete_derivative_response_qualified=0"
       << " slip=" << o.slip_residual
       << " continuity=" << o.radiation_conservation
       << " Euler=" << o.radiation_euler << " normalized_E=" << o.E_normalized
       << " normalized_W=" << o.W_normalized
       << " normalized_trace=" << o.trace_normalized
       << " normalized_slip=" << o.slip_normalized
       << " normalized_conservation=" << o.conservation_normalized
       << " parity=" << o.parity_defect
       << " quadrupole_imaginary=" << o.quadrupole_imaginary
       << " meanU_re=" << o.q0.real() << " meanU_im=" << o.q0.imag()
       << " dipoleU_re=" << o.q1.real() << " dipoleU_im=" << o.q1.imag()
       << " chi=" << e.closure_defect << " Rbg=" << e.background_identity_defect
       << " P=" << e.p << " Perror=" << e.p_error
       << " lambda_cast=" << e.lambda_cast_error << " rhs=" << budget.rhs
       << " scalars=" << budget.scalar_updates
       << " background=" << budget.background_queries
       << " age_evaluations=" << budget.age_quadrature_evaluations
       << " momentum_callbacks=" << budget.momentum_callbacks;
  if (rule)
    text << " nodes=" << rule->mu.size() << " m2=" << rule->m2
         << " m4=" << rule->m4 << " normalization=" << rule->normalization
         << " Gramodd=" << rule->odd_moment
         << " Gramdefect=" << rule->gram_defect;
  text << '\n';
  return ledger.write(text.str(), budget);
}

template <class T> struct NodeTranslationWitness {
  C<T> before{}, after{}, translated_difference{};
  T scale = 0;
  C<T> remainder_change{};
};
// Every node is checked and retained in bounded blocks. Maxima alone cannot
// substitute for the original per-node preservation witness. Both raw D
// subtraction and the small-coordinate translation identity are recorded.
template <class T, class F>
bool record_gram_nodes(Ledger &ledger, Budget &budget, const Controls &control,
                       long double eta, std::size_t nodes, const char *kind,
                       F witness) {
  long double maximum_raw = 0, maximum_translated = 0;
  for (std::size_t first = 0; first < nodes; first += 24) {
    std::ostringstream record;
    record << std::setprecision(std::numeric_limits<long double>::max_digits10)
           << "Gram-node-preservation run=" << control.run_id
           << " stage=" << budget.attempted_stages << " kind=" << kind
           << " eta=" << eta << " first=" << first;
    const auto last = std::min(nodes, first + 24);
    for (std::size_t j = first; j < last; ++j) {
      // The callbacks below form 24 active real scalar destinations including
      // their low-basis/remainder reconstructions, returned tuple, subtraction,
      // scale and maxima. Immutable input reads and text encoding add no state
      // updates; encoded bytes have their own independent cap.
      if (!budget.scalars(24))
        return false;
      const auto node = witness(j);
      const auto raw = node.after - node.before;
      const T normalized =
          node.scale > 0 ? std::abs(node.translated_difference) / node.scale
                         : (node.translated_difference == C<T>{}
                                ? T(0)
                                : std::numeric_limits<T>::infinity());
      maximum_raw =
          std::max(maximum_raw, static_cast<long double>(std::abs(raw)));
      maximum_translated =
          std::max(maximum_translated, static_cast<long double>(normalized));
      record << " node=" << j << ',' << raw.real() << ',' << raw.imag() << ','
             << node.translated_difference.real() << ','
             << node.translated_difference.imag() << ',' << node.scale << ','
             << normalized << ',' << node.remainder_change.real() << ','
             << node.remainder_change.imag();
    }
    record << " checked=" << last - first << '\n';
    if (!ledger.write(record.str(), budget))
      return false;
  }
  std::ostringstream summary;
  summary << std::setprecision(std::numeric_limits<long double>::max_digits10)
          << "Gram-node-summary run=" << control.run_id
          << " stage=" << budget.attempted_stages << " kind=" << kind
          << " eta=" << eta << " checked=" << nodes
          << " maximum_raw_D_change=" << maximum_raw
          << " maximum_normalized_translation=" << maximum_translated << '\n';
  return ledger.write(summary.str(), budget);
}

template <class T>
Status rhs(const Background &background, long double radiation,
           long double radiation_error, T k, long double eta, const State<T> &y,
           State<T> &d, const Rule<T> &r, const Controls &control,
           Budget &budget, Ledger &ledger, Diagnostics &diagnostics) {
  Epoch e;
  Readout<T> o;
  if (!budget.stage()) {
    record_stage(ledger, budget, control, eta, y, e, o, 0, 0, 0, budget.status,
                 "rhs-refusal");
    return budget.status;
  }
  e = shared_epoch(background, std::exp(static_cast<long double>(y.v[0])), k,
                   radiation, radiation_error, control.background_direction,
                   budget, control.background_component);
  if (e.status != Status::ok) {
    record_stage(ledger, budget, control, eta, y, e, o, 0, 0, 0, e.status,
                 "background-refusal");
    return e.status;
  }
  o = readout(y, r, e, &budget);
  if (!resize_residual(d, r.mu.size(), budget)) {
    record_stage(ledger, budget, control, eta, y, e, o, 0, 0, 0, budget.status,
                 "rhs-buffer-refusal", &r);
    return budget.status;
  }
  C<T> M{}, odd{};
  for (std::size_t j = 0; j < r.mu.size(); ++j) {
    if (!budget.scalars(4))
      break;
    const T u = r.mu[j];
    const C<T> U = -2 * T(p) * P2(u) / 3 + y.R[j];
    M += r.weight[j] * u * u * U;
    odd += r.weight[j] * u * P2(u) * U;
  }
  M /= r.normalization;
  odd /= r.normalization;
  const T x = std::sqrt(T(e.x2)), g = T(e.g), fr = T(e.fr);
  const T am = T(1.5) * (r.m4 - r.m2 * r.m2);
  const T sN = -2 * o.v - 3 * fr * am * y.v[6] +
               2 * T(s_star) * T(e.closure_defect - e.fc - e.fl) +
               6 * T(p) * fr * (am - T(2) / 15) +
               (C<T>(0, 3 * fr * x) * odd).real();
  if (!budget.scalars(d.v.size())) {
    record_stage(ledger, budget, control, eta, y, e, o, sN, 0, 0, budget.status,
                 "rhs-derivative-assignment-refusal", &r);
    return budget.status;
  }
  d.v[0] = 1;
  d.v[1] = y.v[2];
  // Stable transcription of ORIGINAL trace, including the actual closure.
  d.v[2] = sN - (3 + g) * y.v[2] +
           T(3 * (e.fr - e.fl) - e.closure_defect) * y.v[1] +
           T(-e.fr + 3 * e.fl + e.closure_defect) * o.v + fr * y.v[5] / 2 -
           x * x * o.s / 3 - T(p) * T(e.closure_defect + 3 * e.fl) -
           T(control.kappa) * o.E;
  d.v[3] = (g - 1) * y.v[3] + y.v[1] - o.v + T(p) * T(e.fc / 2 + 2 * e.fl) / 2;
  d.v[4] = -x * x * (T(p) / 2 + y.v[3]);
  d.v[5] = x * x * r.m2 * o.b;
  d.v[6] = g * y.v[6] - y.v[5] - 8 * y.v[1] + 4 * o.v -
           2 * T(p) * T(e.fc / 2 + 2 * e.fl) - x * x * M.real() / r.m2;
  const T original_trace_rhs = fr * o.delta_r / 2 - (y.v[2] - sN) -
                               (2 + g) * y.v[2] - (2 * g + 1) * o.psi -
                               x * x * o.s / 3;
  o.original_trace = d.v[2] - original_trace_rhs;
  o.modified_trace = o.original_trace + T(control.kappa) * o.E;
  o.radiation_conservation = d.v[5] - x * x * o.b / 3;
  const T sigma_readout = -x * x * o.q2 / 2;
  o.radiation_euler =
      -d.v[6] / 4 - (-g * o.b / 4 + o.delta_r / 4 - sigma_readout + o.psi);
  auto normalized = [](T numerator, T denominator) {
    return denominator > 0
               ? std::abs(numerator) / denominator
               : (numerator == 0 ? T(0) : std::numeric_limits<T>::infinity());
  };
  o.trace_normalized =
      normalized(o.original_trace,
                 std::abs(d.v[2]) + std::abs(fr * o.delta_r / 2) +
                     std::abs(y.v[2] - sN) + std::abs((2 + g) * y.v[2]) +
                     std::abs((2 * g + 1) * o.psi) + std::abs(x * x * o.s / 3));
  o.conservation_normalized =
      std::max(normalized(o.radiation_conservation,
                          std::abs(d.v[5] + 4 * y.v[2]) +
                              std::abs(x * x * o.b / 3) + std::abs(4 * y.v[2])),
               normalized(o.radiation_euler,
                          std::abs(d.v[6] / 4) + std::abs(g * o.b / 4) +
                              std::abs(o.delta_r / 4) +
                              std::abs(sigma_readout) + std::abs(o.psi)));
  for (std::size_t j = 0; j < r.mu.size(); ++j) {
    if (!budget.scalars(2))
      break;
    const T u = r.mu[j], u2 = u * u;
    const C<T> U = -2 * T(p) * P2(u) / 3 + y.R[j];
    d.R[j] = 2 * g * y.R[j] + y.v[6] * (u2 - r.m2) -
             2 * T(p) * T(e.fc / 2 + 2 * e.fl) * (u2 - T(1) / 3) -
             2 * T(p) * (T(1) / 3 - r.m2) - C<T>(0, x) * (u * U - u * M / r.m2);
  }
  const T sigma = sigma_readout;
  const T fE = fr * x * x * o.b * (r.m2 - T(1) / 3) / 2;
  const T fW = -2 * fr * sigma * (1 - 1 / (3 * r.m2));
  diagnostics.maximum_E =
      std::max(diagnostics.maximum_E, static_cast<long double>(o.E_normalized));
  diagnostics.maximum_W =
      std::max(diagnostics.maximum_W, static_cast<long double>(o.W_normalized));
  diagnostics.maximum_trace = std::max(
      diagnostics.maximum_trace, static_cast<long double>(o.trace_normalized));
  diagnostics.maximum_correction =
      std::max(diagnostics.maximum_correction,
               std::abs(control.kappa * static_cast<long double>(o.E)));
  diagnostics.maximum_slip = std::max(
      diagnostics.maximum_slip, static_cast<long double>(o.slip_normalized));
  diagnostics.maximum_subspace_defect =
      std::max(diagnostics.maximum_subspace_defect,
               std::max({static_cast<long double>(std::abs(o.q0)),
                         static_cast<long double>(std::abs(o.q1)),
                         static_cast<long double>(o.parity_defect),
                         static_cast<long double>(o.quadrupole_imaginary),
                         std::abs(static_cast<long double>(o.q0.imag())),
                         std::abs(static_cast<long double>(o.q1.real()))}));
  diagnostics.maximum_conservation =
      std::max(diagnostics.maximum_conservation,
               static_cast<long double>(o.conservation_normalized));
  diagnostics.maximum_closure =
      std::max(diagnostics.maximum_closure, std::abs(e.closure_defect));
  diagnostics.maximum_background_identity =
      std::max(diagnostics.maximum_background_identity,
               std::abs(e.background_identity_defect));
  if (!budget.scalars(7)) {
    record_stage(ledger, budget, control, eta, y, e, o, sN, fE, fW,
                 budget.status, "scalar-refusal", &r);
    return budget.status;
  }
  for (auto &v : d.v)
    v *= T(e.hcal);
  for (auto &v : d.R) {
    if (!budget.scalars(2))
      break;
    v *= T(e.hcal);
  }
  Status status = budget.status;
  for (T v : d.v)
    if (!std::isfinite(v))
      status = Status::overflow;
  for (auto v : d.R)
    if (!finite(v))
      status = Status::overflow;
  if (!record_stage(ledger, budget, control, eta, y, e, o, sN, fE, fW, status,
                    "rhs", &r))
    return budget.status;
  return status;
}

// Source stage: implementation below is unexecuted until a root compute lease.
// Coarse/fine controls share the same original budgets and file byte guards.
template <class T>
Result run(const Background &background, long double radiation,
           long double radiation_error, long double k, long double start,
           long double target, const Controls &control, Budget &budget,
           Ledger &ledger) {
  Result out;
  long double reached_eta = 0, reached_N = 0;
  bool reached_state_known = false;
  auto finish = [&](Status status) {
    out.status = status;
    std::ostringstream endpoint_record;
    endpoint_record << std::setprecision(
                           std::numeric_limits<long double>::max_digits10)
                    << "run-end run=" << control.run_id
                    << " status=" << static_cast<int>(status) << " k=" << k
                    << " start=" << start << " requested_target=" << target
                    << " reached_state_known=" << (reached_state_known ? 1 : 0)
                    << " reached_N=" << reached_N
                    << " reached_eta=" << reached_eta
                    << " eta_error=" << out.eta_error
                    << " endpoint_log_error=" << out.endpoint_log_error
                    << " rhs=" << budget.rhs
                    << " attempted_stages=" << budget.attempted_stages
                    << " scalars=" << budget.scalar_updates
                    << " background=" << budget.background_queries
                    << " age_evaluations=" << budget.age_quadrature_evaluations
                    << " momentum_callbacks=" << budget.momentum_callbacks
                    << " logging_refusals=" << budget.logging_refusals;
    if (status == Status::ok)
      for (unsigned field = 0; field < out.value.size(); ++field)
        endpoint_record << " X" << field << '=' << out.value[field];
    endpoint_record << '\n';
    if (!ledger.write(endpoint_record.str(), budget) && status == Status::ok)
      out.status = budget.status;
    out.rhs = budget.rhs;
    out.scalar_updates = budget.scalar_updates;
    out.background_queries = budget.background_queries;
    out.age_quadrature_evaluations = budget.age_quadrature_evaluations;
    out.ledger_bytes = budget.ledger_bytes;
    return out;
  };
  if (budget.status != Status::ok || ledger.status != Status::ok)
    return finish(budget.status != Status::ok ? budget.status : ledger.status);
  if (!(k > 0) || !std::isfinite(k) || !(start > 0) || !(target > 16 * start) ||
      target < 1e-4L || target > 1 || !std::isfinite(target) ||
      !(control.maximum_log_step > 0 && control.maximum_log_step <= .04L) ||
      !(control.maximum_phase_step > 0 && control.maximum_phase_step <= .04L) ||
      !(control.kappa == 1.5L || control.kappa == 2 || control.kappa == 3) ||
      !(control.absolute_state_tolerance > 0 &&
        control.absolute_state_tolerance <= 2e-12L) ||
      !(control.relative_state_tolerance > 0 &&
        control.relative_state_tolerance <= 2e-11L))
    return finish(Status::outside_domain);
  // Seven DP slopes, trial/base/embedded states, rule and bounded log record.
  const std::size_t payload =
      sizeof(State<T>) * 12 +
      control.nodes * (24 * sizeof(C<T>) + 4 * sizeof(T)) + 16384;
  if (payload > budget.maximum_payload_bytes)
    return finish(Status::work_limit);
  auto rule = make_rule<T>(control.nodes, budget);
  if (rule.status != Status::ok || budget.status != Status::ok)
    return finish(budget.status != Status::ok ? budget.status : rule.status);
  const auto age =
      initial_age(background, start, radiation, radiation_error, budget);
  if (age.status != Status::ok)
    return finish(age.status);
  auto epoch = shared_epoch(background, start, k, radiation, radiation_error,
                            control.background_direction, budget,
                            control.background_component);
  if (epoch.status != Status::ok)
    return finish(epoch.status);
  if (background.source().omega_cdm * start / radiation > 1e-8L ||
      epoch.fl > 1e-12L || k * (age.value + age.error) > 1e-4L)
    return finish(Status::outside_domain);
  const T initial_eta = T(age.value + control.age_direction * age.error);
  auto initial =
      initialize(rule, epoch, T(k), initial_eta, T(age.error), &budget);
  if (budget.status != Status::ok) {
    std::ostringstream refused;
    refused << std::setprecision(std::numeric_limits<long double>::max_digits10)
            << "initial-Gram-refusal run=" << control.run_id
            << " status=" << static_cast<int>(budget.status)
            << " diagnostics_valid=0 alpha_re="
            << initial.translation.alpha.real()
            << " alpha_im=" << initial.translation.alpha.imag()
            << " beta_re=" << initial.translation.beta.real()
            << " beta_im=" << initial.translation.beta.imag()
            << " nodes_translated=" << initial.translation.nodes_translated
            << " scalars=" << budget.scalar_updates << '\n';
    ledger.write(refused.str(), budget);
    return finish(budget.status);
  }
  // State construction writes seven active zeros, then the source assignments
  // and vector value initialization are separately charged actual writes.
  if (!budget.scalars(7))
    return finish(budget.status);
  State<T> y;
  if (!resize_residual(y, rule.mu.size(), budget) ||
      !assign_initial_scalars(y, budget, [&](auto &v) {
        v[0] = std::log(T(start));
        v[1] = 0;
        v[2] = initial.Z; // ORIGINAL source Z, never a momentum reset
        v[3] = initial.ideal.Delta / 3;
        v[4] = 0;
        // Exact identities of the same translated primitive witness. Retain
        // the tiny projected velocity without an O(1) subtraction.
        v[5] = T(k * k) * initial.translation.alpha.real() -
               4 * (initial.actual.psi - initial.ideal.psi);
        v[6] = -4 * initial.ideal.Delta / 3 -
               (C<T>(0, T(k * epoch.hcal)) * initial.translation.beta).real();
      }))
    return finish(budget.status);
  const T H2 = T(epoch.hcal) * T(epoch.hcal);
  for (std::size_t j = 0; j < y.R.size(); ++j) {
    if (!budget.scalars(2))
      return finish(budget.status);
    y.R[j] = H2 * initial.angular.Q[j] + 2 * T(p) * P2(rule.mu[j]) / 3;
  }
  out.diagnostics.initial_E = initial.E;
  out.diagnostics.initial_W = initial.W;
  out.diagnostics.initial_zeta_translation =
      epoch.fr * (initial.actual.delta - initial.ideal.delta_r) /
      (3 * (epoch.fc + 4 * epoch.fr / 3));
  std::ostringstream witness;
  witness << std::setprecision(std::numeric_limits<long double>::max_digits10)
          << "initial-projected-witness run=" << control.run_id
          << " source_asymptotic_zeta=1 source_p=" << p
          << " eta=" << initial.ideal.eta
          << " eta_error=" << initial.ideal.eta_error
          << " H=" << initial.ideal.H << " phi=" << initial.ideal.phi
          << " psi=" << initial.ideal.psi
          << " delta_c=" << initial.ideal.delta_c
          << " delta_r=" << initial.ideal.delta_r << " Vc=" << initial.ideal.Vc
          << " Vr=" << initial.ideal.Vr << " sigma=" << initial.ideal.sigma
          << " Delta=" << initial.ideal.Delta << " Z=" << initial.ideal.Z
          << " sampled_shear_factor=" << 5 * rule.squared_p2
          << " actual_sigma=" << initial.actual.sigma
          << " actual_psi=" << initial.actual.psi
          << " actual_delta_r=" << initial.actual.delta
          << " actual_Vr=" << initial.actual.Vr
          << " alpha_re=" << initial.translation.alpha.real()
          << " alpha_im=" << initial.translation.alpha.imag()
          << " beta_re=" << initial.translation.beta.real()
          << " beta_im=" << initial.translation.beta.imag()
          << " D_node_change=" << initial.translation.maximum_node_change
          << " E=" << initial.E << " W=" << initial.W
          << " finite_zeta_translation="
          << out.diagnostics.initial_zeta_translation
          << " Z_reset=0 background_component=" << control.background_component
          << " background_direction=" << control.background_direction
          << " age_direction=" << control.age_direction << '\n';
  if (!ledger.write(witness.str(), budget))
    return finish(budget.status);
  if (!record_gram_nodes<T>(
          ledger, budget, control, initial_eta, rule.mu.size(), "initial",
          [&](std::size_t j) {
            const T u = rule.mu[j];
            const C<T> old_Q =
                -2 * T(p) * initial_eta * initial_eta * P2(u) / 3;
            const T old_A = initial.ideal.delta_r + 4 * initial.ideal.psi;
            const C<T> old_B = -4 * initial.ideal.Vr / initial.ideal.H;
            const C<T> old_D =
                old_A + C<T>(0, T(k) * u) * old_B + T(k * k) * old_Q;
            const C<T> new_D = initial.angular.A +
                               C<T>(0, T(k) * u) * initial.angular.B +
                               T(k * k) * initial.angular.Q[j];
            return NodeTranslationWitness<T>{
                old_D, new_D,
                (initial.angular.A - old_A) +
                    C<T>(0, T(k) * u) * (initial.angular.B - old_B) +
                    T(k * k) * (initial.angular.Q[j] - old_Q),
                std::abs(old_A) + std::abs(T(k) * u * old_B) +
                    std::abs(T(k * k) * old_Q),
                initial.angular.Q[j] - old_Q};
          }))
    return finish(budget.status);
  const auto initial_readout = readout(y, rule, epoch, &budget);
  if (!record_stage(ledger, budget, control, initial_eta, y, epoch,
                    initial_readout, 0, 0, 0, budget.status,
                    "initial-Gram-readout-originalZ", &rule))
    return finish(budget.status);
  const T target_N = std::log(T(target));
  long double eta = initial_eta;
  reached_eta = eta;
  reached_N = y.v[0];
  reached_state_known = true;
  long double proposed_step = control.maximum_log_step / epoch.hcal;
  const long double epsilon = std::numeric_limits<T>::epsilon();
  if (!budget.scalars(7 * 7))
    return finish(budget.status);
  std::array<State<T>, 7> slopes;
  for (auto &slope : slopes)
    if (!resize_residual(slope, rule.mu.size(), budget))
      return finish(budget.status);
  if (!budget.scalars(2 * active_scalars(y)))
    return finish(budget.status);
  State<T> trial = y, embedded = y;
  constexpr long double node[7]{0, 1.L / 5, 3.L / 10, 4.L / 5, 8.L / 9, 1, 1};
  constexpr long double a[7][7]{
      {},
      {1.L / 5},
      {3.L / 40, 9.L / 40},
      {44.L / 45, -56.L / 15, 32.L / 9},
      {19372.L / 6561, -25360.L / 2187, 64448.L / 6561, -212.L / 729},
      {9017.L / 3168, -355.L / 33, 46732.L / 5247, 49.L / 176, -5103.L / 18656},
      {35.L / 384, 0, 500.L / 1113, 125.L / 192, -2187.L / 6784, 11.L / 84}};
  constexpr long double fifth[7]{
      35.L / 384, 0, 500.L / 1113, 125.L / 192, -2187.L / 6784, 11.L / 84, 0};
  constexpr long double fourth[7]{
      5179.L / 57600, 0,       7571.L / 16695, 393.L / 640, -92097.L / 339200,
      187.L / 2100,   1.L / 40};
  auto combine = [&](State<T> &to, long double h,
                     const long double *coefficients, unsigned count) {
    if (!copy_state(to, y, budget))
      return false;
    for (unsigned l = 0; l < count; ++l) {
      if (coefficients[l] == 0)
        continue;
      if (!budget.scalars(7 + 2 * y.R.size()))
        return false;
      const T coefficient = T(h * coefficients[l]);
      for (unsigned j = 0; j < y.v.size(); ++j)
        to.v[j] += coefficient * slopes[l].v[j];
      for (std::size_t j = 0; j < y.R.size(); ++j)
        to.R[j] += coefficient * slopes[l].R[j];
    }
    return true;
  };
  while (y.v[0] < target_N) {
    const long double gap = static_cast<long double>(target_N - y.v[0]);
    if (gap <= 64 * epsilon *
                   std::max(1.L, std::abs(static_cast<long double>(target_N))))
      break;
    epoch =
        shared_epoch(background, std::exp(static_cast<long double>(y.v[0])), k,
                     radiation, radiation_error, control.background_direction,
                     budget, control.background_component);
    if (epoch.status != Status::ok)
      return finish(epoch.status);
    const auto base_readout = readout(y, rule, epoch, &budget);
    if (budget.status != Status::ok) {
      record_stage(ledger, budget, control, eta, y, epoch, base_readout, 0, 0,
                   0, budget.status, "base-readout-refusal", &rule);
      return finish(budget.status);
    }
    const long double h =
        std::min({proposed_step, control.maximum_log_step / epoch.hcal,
                  control.maximum_phase_step / k, .8L * gap / epoch.hcal});
    if (!(h > 0) || eta + h == eta)
      return finish(Status::conditioning_budget_exceeded);
    bool guard_failure = false;
    for (unsigned stage = 0; stage < 7; ++stage) {
      if (!combine(trial, h, a[stage], stage)) {
        record_stage(ledger, budget, control, eta, y, epoch, base_readout, 0, 0,
                     0, budget.status, "stage-assembly-refusal", &rule);
        return finish(budget.status);
      }
      if (trial.v[0] > target_N || trial.v[0] < y.v[0] ||
          static_cast<long double>(trial.v[0] - y.v[0]) >
              control.maximum_log_step) {
        guard_failure = true;
        break;
      }
      const auto cause = rhs(background, radiation, radiation_error, T(k),
                             eta + node[stage] * h, trial, slopes[stage], rule,
                             control, budget, ledger, out.diagnostics);
      if (cause != Status::ok)
        return finish(cause);
    }
    if (guard_failure) {
      if (!record_stage(ledger, budget, control, eta, y, epoch, base_readout, 0,
                        0, 0, Status::outside_domain, "stage-guard-rejection",
                        &rule))
        return finish(budget.status);
      proposed_step = h / 2;
      continue;
    }
    if (!combine(trial, h, fifth, 7) || !combine(embedded, h, fourth, 7)) {
      record_stage(ledger, budget, control, eta, y, epoch, base_readout, 0, 0,
                   0, budget.status, "final-combination-refusal", &rule);
      return finish(budget.status);
    }
    long double error = 0;
    auto scaled_error = [&](long double coarse, long double fine) {
      const long double scale = control.absolute_state_tolerance +
                                control.relative_state_tolerance *
                                    std::max(std::abs(coarse), std::abs(fine));
      return std::abs(coarse - fine) / scale;
    };
    if (!budget.scalars(7 + 2 * y.R.size())) {
      record_stage(ledger, budget, control, eta, y, epoch, base_readout, 0, 0,
                   0, budget.status, "estimator-refusal", &rule);
      return finish(budget.status);
    }
    for (unsigned j = 0; j < y.v.size(); ++j)
      error = std::max(error, scaled_error(embedded.v[j], trial.v[j]));
    for (std::size_t j = 0; j < y.R.size(); ++j) {
      error = std::max(error,
                       scaled_error(embedded.R[j].real(), trial.R[j].real()));
      error = std::max(error,
                       scaled_error(embedded.R[j].imag(), trial.R[j].imag()));
    }
    if (!std::isfinite(error))
      return finish(Status::overflow);
    if (error > 1) {
      proposed_step = h * std::max(.1L, .8L * std::pow(error, -.2L));
      continue; // attempted slopes and diagnostics stay recorded/charged
    }
    if (!copy_state(y, trial, budget)) {
      record_stage(ledger, budget, control, eta, y, epoch, base_readout, 0, 0,
                   0, budget.status, "accepted-copy-refusal", &rule);
      return finish(budget.status);
    }
    eta += h;
    reached_eta = eta;
    reached_N = y.v[0];
    const auto current =
        shared_epoch(background, std::exp(static_cast<long double>(y.v[0])), k,
                     radiation, radiation_error, control.background_direction,
                     budget, control.background_component);
    if (current.status != Status::ok)
      return finish(current.status);
    if (control.finite_recentering) {
      // Scaled Gram translation. Do not reconstruct O(1) F or reset E/W.
      C<T> q0{}, q1{};
      for (std::size_t j = 0; j < y.R.size(); ++j) {
        if (!budget.scalars(4))
          return finish(budget.status);
        const C<T> U = -2 * T(p) * P2(rule.mu[j]) / 3 + y.R[j];
        q0 += rule.weight[j] * U;
        q1 += rule.weight[j] * rule.mu[j] * U;
      }
      if (!budget.scalars(9)) {
        record_stage(ledger, budget, control, eta, y, current, base_readout, 0,
                     0, 0, budget.status, "Gram-coefficient-refusal", &rule);
        return finish(budget.status);
      }
      q0 /= rule.normalization;
      q1 /= rule.normalization;
      const T determinant = rule.m2 - rule.odd_moment * rule.odd_moment;
      const C<T> alpha = (rule.m2 * q0 - rule.odd_moment * q1) / determinant;
      const C<T> beta = (q1 - rule.odd_moment * q0) / determinant;
      out.diagnostics.maximum_subspace_defect =
          std::max(out.diagnostics.maximum_subspace_defect,
                   std::max(std::abs(static_cast<long double>(alpha.imag())),
                            std::abs(static_cast<long double>(beta.real()))));
      const auto before = readout(y, rule, current, &budget);
      std::ostringstream began;
      began << std::setprecision(std::numeric_limits<long double>::max_digits10)
            << "Gram-event-begin run=" << control.run_id
            << " stage=" << budget.attempted_stages << " eta=" << eta
            << " N=" << y.v[0]
            << " coordinates=U,HcalB alpha_re=" << alpha.real()
            << " alpha_im=" << alpha.imag() << " beta_re=" << beta.real()
            << " beta_im=" << beta.imag() << " planned_nodes=" << y.R.size()
            << " E_before=" << before.E << " W_before=" << before.W
            << " scalars=" << budget.scalar_updates
            << " status=" << static_cast<int>(budget.status) << '\n';
      if (!ledger.write(began.str(), budget))
        return finish(budget.status);
      if (budget.status != Status::ok)
        return finish(budget.status);
      if (out.diagnostics.maximum_subspace_defect > 1e-7L)
        return finish(Status::conditioning_budget_exceeded);
      for (std::size_t j = 0; j < y.R.size(); ++j) {
        if (!budget.scalars(2)) {
          std::ostringstream refusal;
          refusal << "Gram-node-assignment-refusal run=" << control.run_id
                  << " stage=" << budget.attempted_stages
                  << " nodes_translated=" << j << " diagnostics_valid=0 status="
                  << static_cast<int>(budget.status)
                  << " scalars=" << budget.scalar_updates << '\n';
          ledger.write(refusal.str(), budget);
          return finish(budget.status);
        }
        y.R[j] -= alpha + beta * rule.mu[j];
      }
      if (!budget.scalars(1)) {
        record_stage(ledger, budget, control, eta, y, current, before, 0, 0, 0,
                     budget.status, "Gram-dipole-assignment-refusal", &rule);
        return finish(budget.status);
      }
      y.v[6] -= (C<T>(0, std::sqrt(T(current.x2))) * beta).real();
      const auto after = readout(y, rule, current, &budget);
      const T delta_correction =
          T(current.x2) * alpha.real() - 4 * (after.psi - before.psi);
      if (!budget.scalars(1)) {
        record_stage(ledger, budget, control, eta, y, current, after, 0, 0, 0,
                     budget.status, "Gram-density-assignment-refusal", &rule);
        return finish(budget.status);
      }
      y.v[5] += delta_correction;
      out.diagnostics.maximum_recenter = std::max(
          out.diagnostics.maximum_recenter,
          std::max({std::abs(static_cast<long double>(delta_correction)),
                    std::abs(static_cast<long double>(after.psi - before.psi)),
                    std::abs(static_cast<long double>(std::sqrt(T(current.x2)) *
                                                      beta.imag()))}));
      auto corrected = readout(y, rule, current, &budget);
      if (!budget.scalars(6))
        return finish(budget.status);
      const T delta_psi = corrected.psi - before.psi;
      const T delta_b = y.v[6] - trial.v[6];
      const T delta_Cr = y.v[5] - trial.v[5];
      const T impulse_E = delta_psi + T(current.fr) * delta_Cr / 2;
      const T impulse_W = delta_psi + T(current.fr) * delta_b / 2;
      const T delta_sigma = -T(current.x2) * (corrected.q2 - before.q2) / 2;
      std::ostringstream event;
      event << std::setprecision(std::numeric_limits<long double>::max_digits10)
            << "Gram-event run=" << control.run_id
            << " stage=" << budget.attempted_stages << " eta=" << eta
            << " N=" << y.v[0] << " coordinates=U,HcalB"
            << " alpha_re=" << alpha.real() << " alpha_im=" << alpha.imag()
            << " beta_re=" << beta.real() << " beta_im=" << beta.imag()
            << " delta_Cr=" << delta_Cr << " delta_b=" << delta_b
            << " delta_psi=" << delta_psi << " delta_sigma=" << delta_sigma
            << " impulse_E=" << impulse_E << " impulse_W=" << impulse_W
            << " E_subtraction_witness=" << corrected.E - before.E
            << " W_subtraction_witness=" << corrected.W - before.W
            << " unchanged_Z=" << y.v[2]
            << " constraint_reset=0 rhs_forcing=0\n";
      if (!ledger.write(event.str(), budget) ||
          !record_gram_nodes<T>(
              ledger, budget, control, eta, rule.mu.size(), "accepted-event",
              [&](std::size_t j) {
                const T u = rule.mu[j], x = std::sqrt(T(current.x2));
                const C<T> old_U = -2 * T(p) * P2(u) / 3 + trial.R[j];
                const C<T> new_U = -2 * T(p) * P2(u) / 3 + y.R[j];
                const C<T> old_D = before.delta_r + 4 * before.psi +
                                   C<T>(0, x * u) * before.b +
                                   T(current.x2) * old_U;
                const C<T> new_D = corrected.delta_r + 4 * corrected.psi +
                                   C<T>(0, x * u) * corrected.b +
                                   T(current.x2) * new_U;
                return NodeTranslationWitness<T>{
                    old_D, new_D,
                    delta_Cr + 4 * delta_psi + C<T>(0, x * u) * delta_b +
                        T(current.x2) * (y.R[j] - trial.R[j]),
                    std::abs(before.delta_r + 4 * before.psi) +
                        std::abs(x * u * before.b) +
                        std::abs(T(current.x2) * old_U),
                    y.R[j] - trial.R[j]};
              }))
        return finish(budget.status);
      if (!record_stage(ledger, budget, control, eta, y, current, corrected, 0,
                        0, 0, budget.status, "Gram-translation", &rule))
        return finish(budget.status);
    }
    if (budget.status != Status::ok)
      return finish(budget.status);
    const long double forcing = out.diagnostics.maximum_conservation +
                                out.diagnostics.maximum_background_identity;
    out.diagnostics.integrated_forcing += h * current.hcal * forcing;
    proposed_step =
        h *
        (error == 0 ? 2 : std::clamp(.9L * std::pow(error, -.2L), .2L, 2.L));
  }
  const auto endpoint =
      shared_epoch(background, std::exp(static_cast<long double>(y.v[0])), k,
                   radiation, radiation_error, control.background_direction,
                   budget, control.background_component);
  if (endpoint.status != Status::ok)
    return finish(endpoint.status);
  const auto final = readout(y, rule, endpoint, &budget);
  out.value = {static_cast<long double>(centered_cdm(y)),
               static_cast<long double>(final.phi),
               static_cast<long double>(final.psi)};
  out.eta = eta;
  out.eta_error = age.error + 128 * epsilon * budget.rhs * std::abs(eta);
  out.endpoint_log_error =
      std::abs(static_cast<long double>(target_N - y.v[0]));
  const auto checked_age =
      endpoint_age(background, target, radiation, radiation_error, budget);
  if (checked_age.status != Status::ok)
    return finish(checked_age.status);
  out.eta_error += checked_age.error + std::abs(out.eta - checked_age.value);
  if (k * (checked_age.value + checked_age.error) > 20 ||
      out.diagnostics.maximum_E > 1e-7L || out.diagnostics.maximum_W > 1e-7L ||
      out.diagnostics.maximum_trace > 1e-7L ||
      out.diagnostics.maximum_conservation > 1e-7L ||
      out.diagnostics.maximum_slip > 1e-7L ||
      out.diagnostics.maximum_subspace_defect > 1e-7L)
    return finish(Status::conditioning_budget_exceeded);
  for (auto value : out.value)
    if (!std::isfinite(value))
      return finish(Status::overflow);
  return finish(budget.status);
}

inline Result campaign(const Background &background, long double k,
                       long double native_start, long double target,
                       const std::string &ledger_path, Budget &budget) {
  Result refusal;
  auto snapshot = [&](Result result) {
    result.rhs = budget.rhs;
    result.scalar_updates = budget.scalar_updates;
    result.background_queries = budget.background_queries;
    result.age_quadrature_evaluations = budget.age_quadrature_evaluations;
    result.ledger_bytes = budget.ledger_bytes;
    return result;
  };
  if (background.status() != Status::ok) {
    refusal.status = background.status();
    return snapshot(refusal);
  }
  const auto &source = background.source();
  if (source.omega_gamma != 0 || source.omega_massless_nonphoton != 0 ||
      source.omega_b != 0 || !(source.omega_cdm > 0) ||
      source.species.empty() || source.species.size() > 16 ||
      !(native_start >= static_cast<long double>(4e-20) &&
        native_start <= 1e-6L))
    return snapshot(refusal);
  for (const auto &species : source.species)
    if (species.mass_ev != 0)
      return snapshot(refusal);
  if (!budget.charge(budget.background_queries, 1,
                     budget.maximum_background_queries)) {
    refusal.status = budget.status;
    return snapshot(refusal);
  }
  irred::cosmology::ThermalPolicy thermal;
  thermal.momentum_method = background.momentum_method();
  const auto radiation = background.scaled_expansion(0, thermal);
  budget.momentum_callbacks += radiation.callbacks;
  if (radiation.status != Status::ok || !(radiation.a4_e2 > 0)) {
    refusal.status = radiation.status != Status::ok ? radiation.status
                                                    : Status::outside_domain;
    return snapshot(refusal);
  }
  Ledger ledger(ledger_path);
  if (ledger.status != Status::ok) {
    refusal.status = ledger.status;
    return snapshot(refusal);
  }
  std::ostringstream header;
  header << std::setprecision(std::numeric_limits<long double>::max_digits10)
         << "campaign method=" << method_id << " k=" << k
         << " native_start=" << native_start << " target=" << target
         << " H0=" << source.h0_km_s_mpc << " cdm=" << source.omega_cdm
         << " radiation=" << radiation.a4_e2
         << " radiation_error=" << radiation.error_estimate
         << " source_method=" << static_cast<int>(background.momentum_method())
         << " preparation_callbacks=" << background.preparation_callbacks()
         << " species_count=" << source.species.size() << '\n';
  for (std::size_t j = 0; j < source.species.size(); ++j) {
    const auto &s = source.species[j];
    header << "species index=" << j << " mass=" << s.mass_ev
           << " temperature=" << s.temperature_today_ev
           << " weight=" << s.statistical_weight << '\n';
  }
  if (!ledger.write(header.str(), budget)) {
    refusal.status = budget.status;
    return snapshot(refusal);
  }
  Controls settings;
  auto solve = [&](const std::string &name, unsigned nodes, unsigned refinement,
                   unsigned start_divisor, long double kappa = 2, int bg = 0,
                   int age = 0, bool center = true,
                   unsigned background_component = 0) {
    settings.run_id = name;
    settings.nodes = nodes;
    const long double factor = std::pow(2.L, refinement);
    settings.maximum_log_step = settings.maximum_phase_step = .04L / factor;
    settings.absolute_state_tolerance = 2e-12L / (factor * factor);
    settings.relative_state_tolerance = 2e-11L / (factor * factor);
    settings.kappa = kappa;
    settings.background_direction = bg;
    settings.background_component = background_component;
    settings.age_direction = age;
    settings.finite_recentering = center;
    return run<long double>(
        background, radiation.a4_e2, radiation.error_estimate, k,
        native_start / start_divisor, target, settings, budget, ledger);
  };
  std::array<Result, 3> angles, times, starts;
  for (unsigned j = 0; j < 3; ++j) {
    angles[j] = solve("angle" + std::to_string(j), 64U << j, 2, 16);
    if (angles[j].status != Status::ok)
      return snapshot(angles[j]);
  }
  Result main = angles[2];
  for (unsigned j = 0; j < 2; ++j) {
    times[j] = solve("time" + std::to_string(j), 256, j, 16);
    if (times[j].status != Status::ok)
      return snapshot(times[j]);
    starts[j] = solve("start" + std::to_string(j), 256, 2, 4U << j);
    if (starts[j].status != Status::ok)
      return snapshot(starts[j]);
  }
  times[2] = starts[2] = main;
  const auto kappa_low = solve("kappa1.5", 256, 2, 16, 1.5L);
  if (kappa_low.status != Status::ok)
    return snapshot(kappa_low);
  const auto kappa_high = solve("kappa3", 256, 2, 16, 3);
  if (kappa_high.status != Status::ok)
    return snapshot(kappa_high);
  std::array<Result, 3> background_low, background_high;
  for (unsigned component = 1; component <= 3; ++component) {
    background_low[component - 1] =
        solve("background" + std::to_string(component) + "-minus", 256, 2, 16,
              2, -1, 0, true, component);
    if (background_low[component - 1].status != Status::ok)
      return snapshot(background_low[component - 1]);
    background_high[component - 1] =
        solve("background" + std::to_string(component) + "-plus", 256, 2, 16, 2,
              1, 0, true, component);
    if (background_high[component - 1].status != Status::ok)
      return snapshot(background_high[component - 1]);
  }
  const auto age_low = solve("age-minus", 256, 2, 16, 2, 0, -1);
  if (age_low.status != Status::ok)
    return snapshot(age_low);
  const auto age_high = solve("age-plus", 256, 2, 16, 2, 0, 1);
  if (age_high.status != Status::ok)
    return snapshot(age_high);
  const auto no_recentering =
      solve("no-finite-recenter", 256, 2, 16, 2, 0, 0, false);
  if (no_recentering.status != Status::ok)
    return snapshot(no_recentering);
  settings.run_id = "binary64";
  settings.finite_recentering = true;
  const auto binary64 =
      run<double>(background, radiation.a4_e2, radiation.error_estimate, k,
                  native_start / 16, target, settings, budget, ledger);
  if (binary64.status != Status::ok)
    return snapshot(binary64);
  auto difference = [](const Result &a, const Result &b, unsigned field) {
    return std::abs(a.value[field] - b.value[field]);
  };
  for (unsigned field = 0; field < 3; ++field) {
    auto &terms = main.contributions[field];
    terms[0] = 8 * std::max(difference(times[0], times[1], field),
                            difference(times[1], times[2], field)) +
               8 * std::max(difference(kappa_low, main, field),
                            difference(kappa_high, main, field));
    terms[1] = 8 * std::max(difference(angles[0], angles[1], field),
                            difference(angles[1], angles[2], field));
    terms[2] = 8 * std::max(difference(starts[0], starts[1], field),
                            difference(starts[1], starts[2], field));
    terms[3] = 8 * std::max(difference(age_low, main, field),
                            difference(age_high, main, field));
    for (unsigned component = 0; component < 3; ++component)
      terms[3] +=
          8 * std::max(difference(background_low[component], main, field),
                       difference(background_high[component], main, field));
    terms[3] += main.endpoint_log_error * (1 + std::abs(main.value[field]));
    terms[3] += k * main.eta_error * (1 + std::abs(main.value[field]));
    terms[4] = 8 * (difference(binary64, main, field) +
                    difference(no_recentering, main, field)) +
               8 * main.diagnostics.integrated_forcing *
                   (1 + std::abs(main.value[field])) +
               128 * std::numeric_limits<long double>::epsilon() *
                   budget.scalar_updates * (1 + std::abs(main.value[field]));
    const long double epsilon = 1e-7L + 1e-4L * std::abs(main.value[field]);
    for (auto term : terms) {
      main.error[field] += term;
      if (!std::isfinite(term) || term > epsilon / 30)
        main.status = Status::conditioning_budget_exceeded;
    }
    if (main.error[field] > epsilon / 6)
      main.status = Status::conditioning_budget_exceeded;
  }
  std::ostringstream summary;
  summary << std::setprecision(std::numeric_limits<long double>::max_digits10)
          << "campaign-summary k=" << k
          << " status=" << static_cast<int>(main.status)
          << " rhs=" << budget.rhs << " scalars=" << budget.scalar_updates;
  for (unsigned field = 0; field < 3; ++field)
    summary << " X" << field << '=' << main.value[field] << " error" << field
            << '=' << main.error[field];
  summary << '\n';
  if (!ledger.write(summary.str(), budget))
    main.status = budget.status;
  return snapshot(main);
}

} // namespace massless_fd_peer
