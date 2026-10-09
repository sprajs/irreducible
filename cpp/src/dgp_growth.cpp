#include "irred/dgp_growth.hpp"
#include "payload_accounting.hpp"
#include <algorithm>
#include <array>
#include <cfenv>
#include <cmath>
#include <limits>
namespace irred::cosmology {
namespace {
using W = long double;
using S = numerics::Status;
using State = std::array<W, 2>;
bool arithmetic() {
  return std::fegetround() == FE_TONEAREST &&
         std::numeric_limits<W>::digits >= 64;
}
struct Background {
  W e, om, mu, hn;
};
Background background(W m, W a) {
  const W x = (1 - m) * std::pow(a, 1.5L) / (2 * std::sqrt(m));
  const W factor = std::hypot(1.L, x) + x;
  const W om = 1 / (factor * factor), om2 = om * om;
  return {std::sqrt(m) * std::pow(a, -1.5L) * factor, om,
          1 - (1 - om2) / (3 * (1 + om2)), -3 * om / (1 + om)};
}
State rhs(W m, W n, State y) {
  const auto b = background(m, std::exp(n));
  return {y[1], -(4 + b.hn) * y[1] - (3 + b.hn - 1.5L * b.mu * b.om) * y[0]};
}
State add(State y, State k, W h) { return {y[0] + h * k[0], y[1] + h * k[1]}; }
State rk4(W m, W n, State y, W h) {
  const auto k1 = rhs(m, n, y), k2 = rhs(m, n + h / 2, add(y, k1, h / 2)),
             k3 = rhs(m, n + h / 2, add(y, k2, h / 2)),
             k4 = rhs(m, n + h, add(y, k3, h));
  return {y[0] + h * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) / 6,
          y[1] + h * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) / 6};
}
void value(DGPValue &out, W x, W error, double tolerance) {
  const double rounded = static_cast<double>(x);
  const W estimate = error + std::abs(x - W(rounded)) +
                     64 * std::numeric_limits<double>::epsilon() * std::abs(x);
  if (!std::isfinite(x) || !std::isfinite(estimate) ||
      !std::isnormal(rounded) || x <= 0) {
    out.status = S::overflow;
    return;
  }
  if (estimate > W(tolerance) * std::abs(x)) {
    out.status = S::conditioning_budget_exceeded;
    return;
  }
  out = {S::ok, rounded, static_cast<double>(estimate)};
}
S growth(W m, W a, DGPRow &row, DGPPolicy p, size_t remaining, W &d, W &f,
         W &de, W &fe) {
  constexpr W initial = 1e-8L;
  const W x = (1 - m) * std::pow(initial, 1.5L) / (2 * std::sqrt(m));
  State y{1 - 11 * x / 12, -11 * x / 8}, error{16 * x * x, 16 * x * x};
  W n = std::log(initial), end = std::log(a), length = end - n, h = .1L;
  const auto cap =
      std::min({remaining, p.maximum_steps_per_point, size_t(200000)});
  while (n < end) {
    if (row.trials >= cap)
      return S::work_limit;
    h = std::min(h, end - n);
    if (!(h > 0) || n + h == n)
      return S::conditioning_budget_exceeded;
    const auto coarse = rk4(m, n, y, h), half = rk4(m, n, y, h / 2),
               fine = rk4(m, n + h / 2, half, h / 2);
    ++row.trials;
    row.callbacks += 12;
    State delta{std::abs(fine[0] - coarse[0]) / 15,
                std::abs(fine[1] - coarse[1]) / 15};
    const W norm = std::max(delta[0], delta[1]);
    const W allocation = W(p.relative_tolerance) * h / length / 128;
    if (!std::isfinite(fine[0]) || !std::isfinite(fine[1]))
      return S::overflow;
    if (norm > allocation) {
      h /= 2;
      ++row.rejected_trials;
      continue;
    }
    y = fine;
    n += h;
    for (unsigned j = 0; j < 2; ++j)
      error[j] += delta[j] + 64 * std::numeric_limits<W>::epsilon() *
                                 std::max(1.L, std::abs(y[j]));
    if (norm < allocation / 32)
      h = std::min(.2L, h * 2);
  }
  if (!(y[0] > error[0]))
    return S::conditioning_budget_exceeded;
  d = a * y[0];
  f = 1 + y[1] / y[0];
  de = a * error[0];
  fe = (error[1] + std::abs(y[1] / y[0]) * error[0]) / (y[0] - error[0]);
  return S::ok;
}
} // namespace
std::optional<size_t> dgp_payload_bound(size_t points) noexcept {
  irred::detail::PayloadAccounting b(sizeof(DGPBatch) + 1024);
  b.add(points, sizeof(DGPRow));
  return b.result();
}
DGPGrowth prepare_dgp_growth(DGPParameters p) {
  DGPGrowth out;
  if (!std::isfinite(p.omega_m0) || !std::isfinite(p.h0_km_s_mpc)) {
    out.status_ = S::nonfinite_input;
    return out;
  }
  if (!arithmetic())
    return out;
  if (p.omega_m0 < .05 || p.omega_m0 > 1 || p.h0_km_s_mpc <= 0 ||
      p.h0_km_s_mpc > 1000) {
    out.status_ = S::outside_domain;
    return out;
  }
  out.parameters_ = p;
  out.status_ = S::ok;
  return out;
}
DGPGrowth::DGPGrowth(DGPGrowth &&o) noexcept { *this = std::move(o); }
DGPGrowth &DGPGrowth::operator=(DGPGrowth &&o) noexcept {
  if (this != &o) {
    parameters_ = o.parameters_;
    status_ = o.status_;
    o.parameters_.reset();
    o.status_ = S::invalid_input;
  }
  return *this;
}
DGPBatch DGPGrowth::evaluate(std::span<const double> points, unsigned outputs,
                             DGPPolicy p) const {
  DGPBatch out;
  if (status_ != S::ok) {
    out.status = status_;
    return out;
  }
  if (!arithmetic() || !outputs || (outputs & ~63u) ||
      !std::isfinite(p.relative_tolerance) || p.relative_tolerance <= 0)
    return out;
  const auto bytes = dgp_payload_bound(points.size());
  if (points.size() > 65536 || points.size() > p.maximum_points || !bytes ||
      *bytes > p.maximum_native_bytes || *bytes > size_t(1024) * 1024 * 1024) {
    out.status = S::work_limit;
    return out;
  }
  out.requested = outputs;
  out.status = S::ok;
  out.rows.reserve(points.size());
  size_t remaining = std::min(p.maximum_total_steps, size_t(4000000));
  for (double a : points) {
    out.rows.push_back({});
    auto &r = out.rows.back();
    r.scale_factor = a;
    S cause = !std::isfinite(a)
                  ? S::nonfinite_input
                  : (a < 1e-4 || a > 1 ? S::outside_domain : S::ok);
    auto refuse = [&](unsigned mask, DGPValue &v, S s) {
      if (outputs & mask)
        v.status = s;
    };
    if (cause != S::ok) {
      refuse(dgp_e, r.e, cause);
      refuse(dgp_h, r.h, cause);
      refuse(dgp_omega_m, r.omega_m, cause);
      refuse(dgp_mu, r.mu, cause);
      refuse(dgp_d, r.d, cause);
      refuse(dgp_f, r.f, cause);
      out.status = cause;
      continue;
    }
    const auto b = background(parameters_->omega_m0, a);
    if (outputs & dgp_e)
      value(r.e, b.e, 0, p.relative_tolerance);
    if (outputs & dgp_h)
      value(r.h, b.e * parameters_->h0_km_s_mpc, 0, p.relative_tolerance);
    if (outputs & dgp_omega_m)
      value(r.omega_m, b.om, 0, p.relative_tolerance);
    if (outputs & dgp_mu)
      value(r.mu, b.mu, 0, p.relative_tolerance);
    W d = a, f = 1, de = 0, fe = 0;
    if ((outputs & (dgp_d | dgp_f)) && parameters_->omega_m0 != 1) {
      if (p.relative_tolerance < 128 * std::numeric_limits<double>::epsilon())
        cause = S::conditioning_budget_exceeded;
      else
        cause = growth(parameters_->omega_m0, a, r, p, remaining, d, f, de, fe);
      remaining -= r.trials;
      out.trials += r.trials;
      out.rejected_trials += r.rejected_trials;
      out.callbacks += r.callbacks;
    }
    if (cause == S::ok) {
      if (outputs & dgp_d)
        value(r.d, d, de, p.relative_tolerance);
      if (outputs & dgp_f)
        value(r.f, f, fe, p.relative_tolerance);
    } else {
      refuse(dgp_d, r.d, cause);
      refuse(dgp_f, r.f, cause);
    }
    for (const auto *v : {&r.e, &r.h, &r.omega_m, &r.mu, &r.d, &r.f})
      if (v->status != S::ok && v->status != S::invalid_input)
        out.status = v->status;
  }
  return out;
}
} // namespace irred::cosmology
