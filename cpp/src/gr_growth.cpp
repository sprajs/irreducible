#include "irred/gr_growth.hpp"
#include "lcdm_state.hpp"
#include "payload_accounting.hpp"
#include <algorithm>
#include <cfenv>
#include <cmath>
#include <limits>
namespace irred::cosmology {
namespace {
using S = numerics::Status;
bool arithmetic() {
  return std::fegetround() == FE_TONEAREST &&
         std::numeric_limits<long double>::digits >= 64 &&
         std::numeric_limits<long double>::max_exponent >= 16384;
}
struct Context {
  long double q;
};
double integrand(double t, const void *ptr) {
  auto q = static_cast<const Context *>(ptr)->q;
  long double x = t, x2 = x * x;
  return (double)(2 * x2 * x2 / std::pow(1 + q * x2 * x2 * x2, 1.5L));
}
void accept(GrowthValue &out, long double value, long double error,
            double tolerance) {
  const double rounded = (double)value,
               reported = (double)(error + std::abs(value - rounded) +
                                   64 * std::numeric_limits<double>::epsilon() *
                                       std::abs(value));
  if (!std::isfinite(value) || value <= 0 || !std::isnormal(rounded) ||
      !std::isnormal(reported)) {
    out.status = S::outside_domain;
    return;
  }
  if (reported > tolerance * std::abs(rounded)) {
    out.status = S::conditioning_budget_exceeded;
    return;
  }
  out = {S::ok, rounded, reported};
}
} // namespace
std::optional<size_t> growth_payload_bound(size_t points) noexcept {
  irred::detail::PayloadAccounting b(sizeof(GrowthBatch) + 1024);
  b.add(points, sizeof(GrowthRow));
  return b.result();
}
GRGrowth prepare_gr_growth(const Expansion &background) {
  GRGrowth out;
  if (!arithmetic())
    return out;
  if (background.status() != Status::ok ||
      !std::holds_alternative<LCDM>(background.specification())) {
    out.status_ = S::outside_domain;
    return out;
  }
  const auto m = std::get<LCDM>(background.specification()).omega_m;
  if (m < 1e-6 || m > 1) {
    out.status_ = S::outside_domain;
    return out;
  }
  out.background_ = background;
  out.status_ = S::ok;
  return out;
}
GRGrowth::GRGrowth(GRGrowth &&o) noexcept { *this = std::move(o); }
GRGrowth &GRGrowth::operator=(GRGrowth &&o) noexcept {
  if (this != &o) {
    background_ = std::move(o.background_);
    status_ = o.status_;
    o.background_.reset();
    o.status_ = S::invalid_input;
  }
  return *this;
}
GrowthBatch GRGrowth::evaluate(std::span<const double> a, unsigned outputs,
                               GrowthPolicy p) const {
  GrowthBatch out;
  if (status_ != S::ok) {
    out.status = status_;
    return out;
  }
  if (!arithmetic() || !outputs || (outputs & ~3u) ||
      !std::isfinite(p.relative_tolerance) || p.relative_tolerance <= 0 ||
      p.maximum_depth > 60)
    return out;
  const auto bytes = growth_payload_bound(a.size());
  if (a.size() > 65536 || a.size() > p.maximum_points || !bytes ||
      *bytes > p.maximum_native_bytes || *bytes > size_t(1024) * 1024 * 1024) {
    out.status = S::work_limit;
    return out;
  }
  out.rows.reserve(a.size());
  out.requested = outputs;
  out.status = S::ok;
  const auto m = std::get<LCDM>(background_->specification()).omega_m;
  size_t remaining = p.maximum_total_callbacks;
  for (double scale : a) {
    out.rows.push_back({scale, {}, {}, 0});
    auto &r = out.rows.back();
    S cause = S::ok;
    long double D = 0, F = 0, d_error = 0, f_error = 0;
    if (!std::isfinite(scale))
      cause = S::nonfinite_input;
    else if (scale < 1e-8 || scale > 1)
      cause = S::outside_domain;
    else if (m == 1) {
      D = scale;
      F = 1;
    } else {
      const long double aw = scale,
                        q = detail::lcdm_lambda(m) * aw * aw * aw / m;
      // Same retained physical identity supplies the scaled expansion
      // coordinate.
      const long double e2 = detail::lcdm_scaled_expansion(m, aw) / m;
      Context c{q};
      auto integral = numerics::integrate(
          integrand, &c, 0, 1,
          {0, p.relative_tolerance / 32,
           std::min(p.maximum_callbacks_per_point, remaining),
           p.maximum_depth});
      r.callbacks = integral.evaluations;
      out.callbacks += r.callbacks;
      remaining -= r.callbacks;
      cause = integral.status;
      if (cause == S::ok) {
        const long double J = integral.value,
                          delta =
                              integral.error_estimate +
                              64 * std::numeric_limits<double>::epsilon() * J;
        if (!(J > delta)) {
          cause = S::conditioning_budget_exceeded;
        } else {
          D = 2.5L * aw * std::sqrt(e2) * J;
          d_error = D * delta / (J - delta);
          const long double positive = 1 / (std::pow(e2, 1.5L) * J),
                            negative = 1.5L / e2;
          F = positive - negative;
          f_error = positive * delta / (J - delta) +
                    64 * std::numeric_limits<double>::epsilon() *
                        (positive + negative);
        }
      }
    }
    if (outputs & growth_d) {
      if (cause == S::ok)
        accept(r.d, D, d_error, p.relative_tolerance);
      else
        r.d.status = cause;
    }
    if (outputs & growth_f) {
      if (cause == S::ok)
        accept(r.f, F, f_error, p.relative_tolerance);
      else
        r.f.status = cause;
    }
  }
  return out;
}
} // namespace irred::cosmology
