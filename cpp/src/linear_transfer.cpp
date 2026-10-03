#include "irred/linear_transfer.hpp"
#include "irred/quantities.hpp"
#include "payload_accounting.hpp"
#include <algorithm>
#include <cfenv>
#include <cmath>
#include <limits>
namespace irred::cosmology {
namespace {
using S = numerics::Status;
using W = long double;
// Delta_c, relative entropy, V_r-V_c, V_c, phi. Never subtract O(1) states
// to form the O(k^2) comoving density.
using State = std::array<W, 5>;
constexpr W c_km_s = W(irred::speed_of_light_m_per_s) / 1000;
bool arithmetic() {
  return std::fegetround() == FE_TONEAREST &&
         std::numeric_limits<W>::digits >= 64 &&
         std::numeric_limits<W>::max_exponent >= 16384;
}
bool policy(const TransferPolicy &p) {
  return arithmetic() && std::isfinite(p.absolute_tolerance) &&
         std::isfinite(p.relative_tolerance) && p.absolute_tolerance >= 0 &&
         p.relative_tolerance > 0 && std::isfinite(p.maximum_log_step) &&
         p.maximum_log_step > 0 && p.maximum_log_step <= .2 &&
         std::isfinite(p.maximum_constraint_residual) &&
         p.maximum_constraint_residual > 0;
}
struct Epoch {
  S status = S::invalid_input;
  W fm = 0, fr = 0, fl = 0, x2 = 0, log_h_derivative = 0;
};
Epoch epoch(const ThermalBackground &b, W a, W k) {
  Epoch out;
  ThermalPolicy p;
  p.momentum_method = b.momentum_method();
  const auto e = b.scaled_expansion(a, p);
  out.status = e.status;
  if (e.status != S::ok)
    return out;
  const auto &m = b.source();
  const W a2 = a * a, a4 = a2 * a2, scaled = e.a4_e2;
  const W hubble_conformal = W(m.h0_km_s_mpc) / c_km_s * std::sqrt(scaled) / a;
  out.fm = W(m.omega_cdm) * a / scaled;
  out.fr = (W(m.omega_gamma) + m.omega_massless_nonphoton) / scaled;
  out.fl = W(*b.omega_lambda()) * a4 / scaled;
  out.x2 = (k / hubble_conformal) * (k / hubble_conformal);
  out.log_h_derivative = -1 + out.fm / 2 + 2 * out.fl;
  if (!std::isfinite(out.x2) || hubble_conformal <= 0)
    out.status = S::overflow;
  return out;
}
State derivative(const State &y, const Epoch &e) {
  const W F = e.fm + 4 * e.fr / 3;
  return {-e.x2 * y[3] + 6 * e.fr * y[2], e.x2 * y[2],
          e.log_h_derivative * y[2] + (y[0] - y[1]) / 3,
          (e.log_h_derivative - 1) * y[3] + y[4],
          -y[4] + 1.5L * F * y[3] + 2 * e.fr * y[2]};
}
W residual(const State &y, const Epoch &e) {
  const std::array<W, 4> terms{e.x2 * y[4], 1.5L * (e.fm + 4 * e.fr / 3) * y[0],
                               -2 * e.fr * y[1], 6 * e.fr * y[2]};
  W total = 0, scale = 0;
  for (auto t : terms) {
    total += t;
    scale += std::abs(t);
  }
  if (scale == 0)
    return total == 0 ? 0 : std::numeric_limits<W>::infinity();
  return std::abs(total) / scale;
}
struct Run {
  S status = S::invalid_input;
  State y{};
  W constraint = 0;
};
Run evolve(const ThermalBackground &b, W k, W start, W target, W step,
           std::size_t &remaining, std::size_t &used) {
  Run out;
  auto e = epoch(b, start, k);
  if (e.status != S::ok) {
    out.status = e.status;
    return out;
  }
  const auto &m = b.source();
  const W radiation = W(m.omega_gamma) + m.omega_massless_nonphoton;
  if (e.fl > 1e-10L) {
    out.status = S::outside_domain;
    return out;
  }
  if (radiation > 0) {
    const W r = W(m.omega_cdm) * start / radiation;
    if (r > 1e-4L || e.x2 > 1e-6L) {
      out.status = S::outside_domain;
      return out;
    }
    const W phi0 = -2.L / 3;
    out.y = {0, 0, -phi0 * e.x2 / 24, phi0 * (.5L + r / 16 - e.x2 / 120),
             phi0 * (1 - r / 16 - e.x2 / 30)};
  } else {
    out.y = {0, 0, 0, -.4L, -.6L};
  }
  const W F = e.fm + 4 * e.fr / 3;
  if (!(F > 0)) {
    out.status = S::outside_domain;
    return out;
  }
  // The correction is confined to omitted orders of the admitted series.
  out.y[0] = (-e.x2 * out.y[4] - 6 * e.fr * out.y[2]) / (1.5L * F);
  out.constraint = residual(out.y, e);
  W n = std::log(start), end = std::log(target);
  auto rhs = [&](W coordinate, const State &state, State &result) -> S {
    if (!remaining)
      return S::work_limit;
    --remaining;
    ++used;
    const auto current = epoch(b, std::exp(coordinate), k);
    if (current.status != S::ok)
      return current.status;
    result = derivative(state, current);
    for (auto v : result)
      if (!std::isfinite(v))
        return S::overflow;
    return S::ok;
  };
  while (n < end) {
    const W h = std::min(end - n, step / (1 + std::sqrt(e.x2)));
    if (!(h > 0) || n + h == n) {
      out.status = S::conditioning_budget_exceeded;
      return out;
    }
    State k1{}, k2{}, k3{}, k4{}, tmp{};
    S cause = rhs(n, out.y, k1);
    if (cause != S::ok) {
      out.status = cause;
      return out;
    }
    for (unsigned j = 0; j < 5; ++j)
      tmp[j] = out.y[j] + h * k1[j] / 2;
    cause = rhs(n + h / 2, tmp, k2);
    if (cause != S::ok) {
      out.status = cause;
      return out;
    }
    for (unsigned j = 0; j < 5; ++j)
      tmp[j] = out.y[j] + h * k2[j] / 2;
    cause = rhs(n + h / 2, tmp, k3);
    if (cause != S::ok) {
      out.status = cause;
      return out;
    }
    for (unsigned j = 0; j < 5; ++j)
      tmp[j] = out.y[j] + h * k3[j];
    cause = rhs(n + h, tmp, k4);
    if (cause != S::ok) {
      out.status = cause;
      return out;
    }
    for (unsigned j = 0; j < 5; ++j) {
      out.y[j] += h * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]) / 6;
      if (!std::isfinite(out.y[j])) {
        out.status = S::overflow;
        return out;
      }
    }
    n = std::min(end, n + h);
    e = epoch(b, std::exp(n), k);
    if (e.status != S::ok) {
      out.status = e.status;
      return out;
    }
    out.constraint = std::max(out.constraint, residual(out.y, e));
  }
  out.status = S::ok;
  return out;
}
void accept(TransferValue &out, W value, W time, W initial, std::size_t work,
            const TransferPolicy &p) {
  const double rounded = static_cast<double>(value);
  const W arithmetic_error =
      std::abs(value - W(rounded)) +
      128 * std::numeric_limits<W>::epsilon() * W(work + 1) * std::abs(value) +
      64 * std::numeric_limits<double>::epsilon() * std::abs(value);
  const W error = time + initial + arithmetic_error;
  out.time_refinement = static_cast<double>(time);
  out.initial_refinement = static_cast<double>(initial);
  out.absolute_error_estimate = static_cast<double>(error);
  if (!std::isfinite(rounded) || (value != 0 && !std::isnormal(rounded)) ||
      !std::isfinite(out.absolute_error_estimate) ||
      (error > 0 && !std::isnormal(out.absolute_error_estimate))) {
    out.status = S::outside_domain;
    return;
  }
  if (error >
      W(p.absolute_tolerance) + W(p.relative_tolerance) * std::abs(value)) {
    out.status = S::conditioning_budget_exceeded;
    return;
  }
  out.status = S::ok;
  out.value = rounded;
}
W window(W x) {
  if (std::abs(x) < .05L) {
    const W x2 = x * x;
    return 1 - x2 / 10 + x2 * x2 / 280 - x2 * x2 * x2 / 15120 +
           x2 * x2 * x2 * x2 / 1330560;
  }
  return 3 * (std::sin(x) - x * std::cos(x)) / (x * x * x);
}
} // namespace
std::optional<std::size_t> transfer_payload_bound(std::size_t points) noexcept {
  irred::detail::PayloadAccounting b(sizeof(PerfectFluidTransfer) +
                                     sizeof(TransferBatch) + 4096);
  b.add(points, sizeof(TransferRow));
  return b.result();
}
PerfectFluidTransfer prepare_perfect_fluid_transfer(const ThermalBackground &b,
                                                    double initial) {
  PerfectFluidTransfer out;
  if (!arithmetic() || !std::isfinite(initial))
    return out;
  if (b.status() != S::ok) {
    out.status_ = b.status();
    return out;
  }
  const auto &m = b.source();
  if (initial < 1e-20 || initial > 1e-3 || m.omega_b != 0 ||
      !m.species.empty() ||
      m.omega_cdm + W(m.omega_gamma) + m.omega_massless_nonphoton <= 0) {
    out.status_ = S::outside_domain;
    return out;
  }
  out.background_ = b;
  out.initial_ = initial;
  out.status_ = S::ok;
  return out;
}
PerfectFluidTransfer::PerfectFluidTransfer(PerfectFluidTransfer &&o) noexcept {
  *this = std::move(o);
}
PerfectFluidTransfer &
PerfectFluidTransfer::operator=(PerfectFluidTransfer &&o) noexcept {
  if (this != &o) {
    background_ = std::move(o.background_);
    initial_ = o.initial_;
    status_ = o.status_;
    o.background_.reset();
    o.initial_ = 0;
    o.status_ = S::invalid_input;
  }
  return *this;
}
TransferBatch PerfectFluidTransfer::evaluate(std::span<const double> ks,
                                             double a, unsigned outputs,
                                             TransferPolicy p) const {
  TransferBatch out;
  if (status_ != S::ok) {
    out.status = status_;
    return out;
  }
  if (!policy(p) || !outputs || (outputs & ~3u))
    return out;
  const auto bytes = transfer_payload_bound(ks.size());
  if (ks.size() > 65536 || ks.size() > p.maximum_points || !bytes ||
      *bytes > p.maximum_native_bytes ||
      *bytes > std::size_t(1024) * 1024 * 1024) {
    out.status = S::work_limit;
    return out;
  }
  if (!std::isfinite(a)) {
    out.status = S::nonfinite_input;
    return out;
  }
  if (a < initial_ || a > 1) {
    out.status = S::outside_domain;
    return out;
  }
  out.status = S::ok;
  out.scale_factor = a;
  out.requested_outputs = outputs;
  out.rows.reserve(ks.size());
  std::size_t remaining = p.maximum_total_callbacks;
  for (double k : ks) {
    out.rows.push_back({k, {}, {}, 0, 0});
    auto &row = out.rows.back();
    S cause = S::ok;
    if (!std::isfinite(k))
      cause = S::nonfinite_input;
    else if (k <= 0 || k > 1)
      cause = S::outside_domain;
    std::array<Run, 5> runs{};
    std::size_t point_remaining =
        std::min(remaining, p.maximum_callbacks_per_point);
    if (cause == S::ok) {
      for (unsigned i = 0; i < 5; ++i) {
        const W start = W(initial_) / (i == 0 ? 1 : i == 1 ? 2 : 4);
        const W step = W(p.maximum_log_step) / (i < 3 ? 4 : i == 3 ? 2 : 1);
        runs[i] = evolve(*background_, k, start, a, step, point_remaining,
                         row.callbacks);
        if (runs[i].status != S::ok) {
          cause = runs[i].status;
          break;
        }
      }
    }
    remaining -= row.callbacks;
    out.callbacks += row.callbacks;
    if (cause == S::ok) {
      row.maximum_constraint_residual = static_cast<double>(std::max(
          {runs[0].constraint, runs[1].constraint, runs[2].constraint}));
      if (row.maximum_constraint_residual > p.maximum_constraint_residual)
        cause = S::conditioning_budget_exceeded;
    }
    auto output = [&](unsigned mask, unsigned coordinate,
                      TransferValue &value) {
      if (!(outputs & mask))
        return;
      if (cause != S::ok) {
        value.status = cause;
        return;
      }
      const W time =
          8 * std::max(std::abs(runs[2].y[coordinate] - runs[3].y[coordinate]),
                       std::abs(runs[3].y[coordinate] - runs[4].y[coordinate]));
      const W initial =
          8 * std::max(std::abs(runs[0].y[coordinate] - runs[1].y[coordinate]),
                       std::abs(runs[1].y[coordinate] - runs[2].y[coordinate]));
      accept(value, runs[2].y[coordinate], time, initial, row.callbacks, p);
    };
    output(transfer_comoving_cdm, 0, row.comoving_cdm);
    output(transfer_metric, 4, row.metric);
  }
  return out;
}
BandVariance PerfectFluidTransfer::band_variance(const PrimordialBand &law,
                                                 double a, double radius,
                                                 BandVariancePolicy p) const {
  BandVariance out;
  out.primordial = law;
  out.scale_factor = a;
  out.radius_mpc = radius;
  if (status_ != S::ok) {
    out.status = status_;
    return out;
  }
  for (double v : {law.amplitude, law.spectral_index, law.pivot_mpc_inverse,
                   law.minimum_mpc_inverse, law.maximum_mpc_inverse, a, radius,
                   p.absolute_tolerance, p.relative_tolerance})
    if (!std::isfinite(v)) {
      out.status = S::nonfinite_input;
      return out;
    }
  if (!policy(p.transfer) || p.absolute_tolerance < 0 ||
      p.relative_tolerance <= 0 || law.amplitude <= 0 ||
      law.pivot_mpc_inverse <= 0 || law.minimum_mpc_inverse <= 0 ||
      law.maximum_mpc_inverse <= law.minimum_mpc_inverse || radius <= 0 ||
      p.base_panels < 2 || p.base_panels % 2 || p.base_panels > 1024 ||
      background_->source().omega_cdm <= 0)
    return out;
  if (a < initial_ || a > 1) {
    out.status = S::outside_domain;
    return out;
  }
  const std::size_t panels = 4 * p.base_panels;
  // Transfer rows + simultaneous k/integrand/error vectors.
  const auto bound = transfer_payload_bound(panels + 1);
  irred::detail::PayloadAccounting bytes(
      bound.value_or(std::numeric_limits<std::size_t>::max()));
  bytes.add(panels + 1, sizeof(double) + 2 * sizeof(W));
  if (!bound || !bytes.result() ||
      *bytes.result() > p.transfer.maximum_native_bytes) {
    out.status = S::work_limit;
    return out;
  }
  const W lower = std::log(W(law.minimum_mpc_inverse));
  const W span = std::log(W(law.maximum_mpc_inverse) / law.minimum_mpc_inverse);
  std::vector<double> ks(panels + 1);
  for (std::size_t j = 0; j <= panels; ++j)
    ks[j] = static_cast<double>(std::exp(lower + span * W(j) / panels));
  ks.front() = law.minimum_mpc_inverse;
  ks.back() = law.maximum_mpc_inverse;
  for (std::size_t j = 1; j <= panels; ++j)
    if (!(ks[j] > ks[j - 1])) {
      out.status = S::conditioning_budget_exceeded;
      return out;
    }
  // Prevent nested-grid agreement caused by unresolved window/acoustic phases.
  // Each single positive background component supplies an upper eta bound.
  const auto &source = background_->source();
  const W rad = W(source.omega_gamma) + source.omega_massless_nonphoton;
  W eta_bound = std::numeric_limits<W>::infinity();
  if (rad > 0)
    eta_bound = c_km_s * a / (W(source.h0_km_s_mpc) * std::sqrt(rad));
  if (source.omega_cdm > 0)
    eta_bound = std::min(eta_bound, 2 * c_km_s * std::sqrt(W(a)) /
                                        (W(source.h0_km_s_mpc) *
                                         std::sqrt(W(source.omega_cdm))));
  for (std::size_t j = 4; j <= panels; j += 4) {
    const W delta_k = W(ks[j]) - ks[j - 4];
    if (delta_k * radius > .25L ||
        (rad > 0 && delta_k * eta_bound / std::sqrt(3.L) > .5L)) {
      out.status = S::outside_domain;
      return out;
    }
  }
  auto transfer = evaluate(ks, a, transfer_comoving_cdm, p.transfer);
  out.callbacks = transfer.callbacks;
  if (transfer.status != S::ok) {
    out.status = transfer.status;
    return out;
  }
  std::vector<W> integrand(panels + 1), errors(panels + 1);
  for (std::size_t j = 0; j <= panels; ++j) {
    const auto &value = transfer.rows[j].comoving_cdm;
    if (value.status != S::ok) {
      out.status = value.status;
      return out;
    }
    const W t = *value.value, error = value.absolute_error_estimate;
    const W primordial =
        W(law.amplitude) *
        std::pow(W(ks[j]) / law.pivot_mpc_inverse, W(law.spectral_index) - 1);
    if (!(primordial > 0) || !std::isnormal(primordial)) {
      out.status = S::outside_domain;
      return out;
    }
    const W w = window(W(ks[j]) * radius), weight = primordial * w * w;
    integrand[j] = weight * t * t;
    errors[j] = weight * (2 * std::abs(t) * error + error * error);
    if (!std::isfinite(integrand[j]) || !std::isfinite(errors[j])) {
      out.status = S::overflow;
      return out;
    }
  }
  auto simpson = [&](const std::vector<W> &values, std::size_t stride) {
    W sum = values.front() + values.back();
    for (std::size_t j = stride; j < panels; j += stride)
      sum += (j / stride % 2 ? 4 : 2) * values[j];
    return span * stride * sum / (3 * panels);
  };
  const W fine = simpson(integrand, 1), medium = simpson(integrand, 2),
          coarse = simpson(integrand, 4);
  const W quadrature =
      8 * std::max(std::abs(fine - medium), std::abs(medium - coarse));
  const W inherited = simpson(errors, 1),
          total = quadrature + inherited +
                  128 * std::numeric_limits<double>::epsilon() * std::abs(fine);
  out.window_refinement = static_cast<double>(quadrature);
  out.transfer_error_estimate = static_cast<double>(inherited);
  out.absolute_error_estimate = static_cast<double>(total);
  const double value = static_cast<double>(fine),
               sigma = static_cast<double>(std::sqrt(fine));
  if (!(fine > total) || !std::isnormal(value) || !std::isnormal(sigma) ||
      !std::isnormal(out.absolute_error_estimate)) {
    out.status = S::outside_domain;
    return out;
  }
  if (total > W(p.absolute_tolerance) + W(p.relative_tolerance) * fine) {
    out.status = S::conditioning_budget_exceeded;
    return out;
  }
  const W sigma_error = total / (std::sqrt(fine) + std::sqrt(fine - total)) +
                        std::abs(std::sqrt(fine) - W(sigma));
  out.sigma_absolute_error_estimate = static_cast<double>(sigma_error);
  if (!std::isnormal(out.sigma_absolute_error_estimate)) {
    out.status = S::outside_domain;
    return out;
  }
  out.status = S::ok;
  out.variance = value;
  out.sigma = sigma;
  return out;
}
BandVariance PerfectFluidTransfer::sigma8_band(const PrimordialBand &law,
                                               double a,
                                               BandVariancePolicy p) const {
  if (!background_) {
    BandVariance out;
    out.status = status_;
    return out;
  }
  return band_variance(law, a, 800 / background_->source().h0_km_s_mpc, p);
}
} // namespace irred::cosmology
