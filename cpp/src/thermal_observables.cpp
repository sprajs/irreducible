#include "irred/thermal_observables.hpp"
#include "flat_geometry.hpp"
#include "payload_accounting.hpp"
#include <algorithm>
#include <cfenv>
#include <cmath>
#include <limits>
#include <utility>
namespace irred::cosmology {
namespace {
using S = numerics::Status;
constexpr unsigned bit(EarlyLateOutput x) { return early_late_mask(x); }
constexpr unsigned ratios = bit(EarlyLateOutput::dm_over_rs) |
                            bit(EarlyLateOutput::dh_over_rs) |
                            bit(EarlyLateOutput::dv_over_rs);
constexpr unsigned integrated =
    bit(EarlyLateOutput::dm_mpc) | bit(EarlyLateOutput::dl_mpc) |
    bit(EarlyLateOutput::dv_mpc) | bit(EarlyLateOutput::dm_over_rs) |
    bit(EarlyLateOutput::dv_over_rs);
constexpr long double arithmetic =
    64.L * std::numeric_limits<double>::epsilon();
bool tolerance(double a, double r) {
  return std::isfinite(a) && std::isfinite(r) && a >= 0 && r >= 0 &&
         (a > 0 || r > 0);
}
bool policy_valid(const ThermalObservablePolicy &p) {
  return std::numeric_limits<long double>::digits >= 64 &&
         std::numeric_limits<long double>::max_exponent >= 16384 &&
         std::fegetround() == FE_TONEAREST && p.maximum_depth <= 60 &&
         tolerance(p.absolute_tolerance_mpc, p.relative_tolerance) &&
         tolerance(p.absolute_tolerance_ratio, p.relative_tolerance_ratio) &&
         tolerance(p.thermal.absolute_tolerance,
                   p.thermal.relative_tolerance) &&
         p.thermal.maximum_depth <= 60;
}
void admit(EarlyLateValue &out, long double value, long double error,
           double absolute, double relative) {
  if (value < 0 || !detail::physical_representable(value)) {
    out.status = S::outside_domain;
    return;
  }
  const double rounded = static_cast<double>(value);
  const long double diagnostic =
      error + arithmetic * value + std::abs(value - rounded);
  // Positive diagnostics are never converted into a fabricated exact zero.
  if (!detail::physical_representable(diagnostic) ||
      diagnostic > absolute + relative * value) {
    out.status = S::conditioning_budget_exceeded;
    return;
  }
  out = {S::ok, rounded, static_cast<double>(diagnostic)};
}
struct CallbackWork {
  std::size_t outer = 0, momentum = 0;
  std::size_t total() const { return outer + momentum; }
};
struct Context {
  const ThermalBackground &background;
  const ThermalObservablePolicy &policy;
  CallbackWork &work;
  std::size_t start;
  long double xmax = 0, a_drag = 0, loading = 0, dependency_relative = 0;
  bool sound = false;
  S status = S::ok;
  ThermalScaledExpansion expansion(long double a) {
    ThermalScaledExpansion out;
    if (status != S::ok) {
      out.status = status;
      return out;
    }
    const auto &p = policy;
    const auto used = work.total() - start;
    if (used > p.maximum_callbacks_per_point ||
        work.total() > p.maximum_total_callbacks ||
        work.momentum > p.thermal.maximum_total_callbacks) {
      status = S::work_limit;
      out.status = status;
      return out;
    }
    auto nested = p.thermal;
    nested.maximum_total_callbacks =
        std::min({p.maximum_callbacks_per_point - used,
                  p.maximum_total_callbacks - work.total(),
                  p.thermal.maximum_total_callbacks - work.momentum});
    out = background.scaled_expansion(a, nested);
    work.momentum += out.callbacks;
    status = out.status;
    if (status == S::ok) {
      // Positive inverse-square-root interval. Maximum sampled dependency is
      // an empirical admission diagnostic, not a continuous certified bound.
      const long double radius =
          std::sqrt(out.a4_e2 / (out.a4_e2 - out.error_estimate)) - 1;
      dependency_relative = std::max(dependency_relative, radius);
    }
    return out;
  }
};
double integrand(double t, const void *ptr) {
  auto &c = *const_cast<Context *>(static_cast<const Context *>(ptr));
  const auto &p = c.policy;
  if (c.status != S::ok || c.work.total() >= p.maximum_total_callbacks ||
      c.work.total() - c.start >= p.maximum_callbacks_per_point) {
    if (c.status == S::ok)
      c.status = S::work_limit;
    return std::numeric_limits<double>::quiet_NaN();
  }
  ++c.work.outer;
  const long double a = c.sound ? c.a_drag * t : std::exp(-c.xmax * t);
  const auto q = c.expansion(a);
  if (q.status != S::ok)
    return std::numeric_limits<double>::quiet_NaN();
  const long double value = c.sound
                                ? 1 / std::sqrt(q.a4_e2 * (1 + c.loading * a))
                                : a / std::sqrt(q.a4_e2);
  const double rounded = static_cast<double>(value);
  if (!detail::physical_representable(value) || !(rounded > 0)) {
    c.status = S::outside_domain;
    return std::numeric_limits<double>::quiet_NaN();
  }
  return rounded;
}
struct Integral {
  S status = S::invalid_input;
  long double value = 0, error = 0;
};
Integral integrate(Context &c, long double scale, long double amplification) {
  Integral out;
  const auto &p = c.policy;
  const auto used = c.work.total() - c.start;
  if (used > p.maximum_callbacks_per_point ||
      c.work.total() > p.maximum_total_callbacks) {
    out.status = S::work_limit;
    return out;
  }
  const auto available = std::min(p.maximum_callbacks_per_point - used,
                                  p.maximum_total_callbacks - c.work.total());
  if (available < 3) {
    out.status = S::work_limit;
    return out;
  }
  const double absolute = static_cast<double>(
      std::min(static_cast<long double>(p.absolute_tolerance_mpc) /
                   (8 * scale * amplification),
               static_cast<long double>(std::numeric_limits<double>::max())));
  const double relative = p.relative_tolerance / 8;
  if (absolute == 0 && relative == 0) {
    out.status = S::conditioning_budget_exceeded;
    return out;
  }
  const auto q = numerics::integrate(
      integrand, &c, 0, 1, {absolute, relative, available, p.maximum_depth});
  out.status = c.status == S::ok ? q.status : c.status;
  if (out.status != S::ok)
    return out;
  out.value = scale * q.value;
  out.error =
      scale * q.error_estimate + out.value * c.dependency_relative +
      scale *
          std::abs(static_cast<long double>(std::nextafter(
                       q.value, std::numeric_limits<double>::infinity())) -
                   q.value) /
          2;
  return out;
}
} // namespace
std::optional<std::size_t>
thermal_observables_payload_bound(std::size_t points, std::size_t species,
                                  std::size_t origins) noexcept {
  irred::detail::PayloadAccounting bytes(sizeof(ThermalObservables) +
                                         sizeof(ThermalObservableBatch) +
                                         sizeof(ThermalPhysicalMapping));
  bytes.add(points, 2 * sizeof(EarlyLateRow));
  bytes.add(species,
            4 * sizeof(ThermalSpecies) + 4 * sizeof(ThermalPhysicalSpecies));
  bytes.add(origins, 4);
  // Conservative bounded nonrecursive contexts/temporary owner headers. Moment
  // quadrature itself has no heap workspace; recursion stack remains excluded.
  bytes.add(1, 2048 + 2 * sizeof(ThermalBackground));
  return bytes.result();
}
ThermalObservables
prepare_thermal_observables(const ThermalObservableRequest &r,
                            ThermalObservablePolicy p) {
  ThermalObservables out;
  if (!policy_valid(p) || r.model.species.size() > 16)
    return out;
  std::size_t origins = 0;
  if (!irred::detail::checked_payload_add(origins, r.source_origin.size(), 1) ||
      !irred::detail::checked_payload_add(origins, r.drag_origin.size(), 1) ||
      !irred::detail::checked_payload_add(origins, 2, 1))
    return out;
  const auto bytes =
      thermal_observables_payload_bound(0, r.model.species.size(), origins);
  if (!bytes)
    return out;
  if (*bytes > p.maximum_native_bytes ||
      r.model.species.size() > p.thermal.maximum_species) {
    out.status_ = S::work_limit;
    return out;
  }
  if (!std::isfinite(r.z_drag)) {
    out.status_ = S::nonfinite_input;
    return out;
  }
  if (r.z_drag < 0 || r.drag_origin.empty() || r.source_origin.empty()) {
    out.status_ = S::outside_domain;
    return out;
  }
  // Retain admitted source provenance even when mapping or physical closure
  // fails. Structurally rejected resource requests retain no payload.
  out.source_ = r;
  auto mapping = map_thermal_physical_model(r.model, p.thermal);
  out.status_ = mapping.status;
  if (mapping.status != S::ok)
    return out;
  auto nested = p.thermal;
  nested.maximum_total_callbacks =
      std::min(nested.maximum_total_callbacks, p.maximum_total_callbacks);
  out.background_ = prepare_thermal_background(*mapping.model, nested);
  out.status_ = out.background_.status();
  return out;
}
ThermalObservables::ThermalObservables(ThermalObservables &&other) noexcept {
  *this = std::move(other);
}
ThermalObservables &
ThermalObservables::operator=(ThermalObservables &&other) noexcept {
  if (this == &other)
    return *this;
  status_ = other.status_;
  source_ = std::move(other.source_);
  background_ = std::move(other.background_);
  other.status_ = S::invalid_input;
  other.source_ = {};
  return *this;
}
ThermalObservableBatch
ThermalObservables::evaluate(std::span<const double> zs, unsigned mask,
                             ThermalObservablePolicy p) const {
  ThermalObservableBatch out;
  if (!policy_valid(p) || !mask || mask >= (thermal_ruler_mask << 1) ||
      zs.size() > 65536)
    return out;
  if (status_ != S::ok) {
    out.status = status_;
    return out;
  }
  // Admit the retained numerical identity before zero-distance shortcuts,
  // output masks or an empty batch can bypass background evaluation.
  if (p.thermal.momentum_method != background_.momentum_method())
    return out;
  std::size_t origins = 0;
  if (!irred::detail::checked_payload_add(origins, source_.source_origin.size(),
                                          1) ||
      !irred::detail::checked_payload_add(origins, source_.drag_origin.size(),
                                          1) ||
      !irred::detail::checked_payload_add(origins, 2, 1))
    return out;
  const auto bytes = thermal_observables_payload_bound(
      zs.size(), source_.model.species.size(), origins);
  if (!bytes)
    return out;
  if (*bytes > p.maximum_native_bytes || zs.size() > p.maximum_points ||
      source_.model.species.size() > p.thermal.maximum_species) {
    out.status = S::work_limit;
    return out;
  }
  out.status = S::ok;
  out.requested_outputs = mask;
  out.rows.reserve(zs.size());
  const long double scale =
      detail::prepare_flat_scale(source_.model.h0_km_s_mpc).distance_mpc;
  CallbackWork work;
  if ((mask & thermal_ruler_mask) || ((mask & ratios) && !zs.empty())) {
    out.ruler = EarlyLateValue{};
    Context c{background_, p, work, work.total()};
    c.sound = true;
    c.a_drag = 1 / (1 + static_cast<long double>(source_.z_drag));
    const auto &mapped = background_.source();
    c.loading = 3 * static_cast<long double>(mapped.omega_b) /
                (4 * static_cast<long double>(mapped.omega_gamma));
    const auto integral = integrate(c, scale * c.a_drag / std::sqrt(3.L), 1);
    out.ruler->status = integral.status;
    if (integral.status == S::ok)
      admit(*out.ruler, integral.value, integral.error,
            p.absolute_tolerance_mpc, p.relative_tolerance);
  }
  for (double z : zs) {
    out.rows.emplace_back();
    auto &row = out.rows.back();
    row.redshift = z;
    const auto start = work.total();
    auto fail = [&](S status, unsigned selected) {
      for (unsigned i = 0; i < early_late_output_count; ++i)
        if (mask & selected & (1u << i))
          row.outputs[i].status = status;
    };
    auto scalar = [&](EarlyLateOutput id, long double value, long double error,
                      double absolute, double relative) {
      if (mask & bit(id))
        admit(row.outputs[static_cast<unsigned>(id)], value, error, absolute,
              relative);
    };
    if (!std::isfinite(z) || z < 0) {
      fail(std::isfinite(z) ? S::outside_domain : S::nonfinite_input, mask);
      continue;
    }
    const long double x = std::log1p(static_cast<long double>(z)),
                      a = std::exp(-x), a2 = a * a;
    Context c{background_, p, work, start};
    c.xmax = x;
    const bool need_local =
        (mask & (~integrated & ((1u << early_late_output_count) - 1))) ||
        (mask & bit(EarlyLateOutput::dv_mpc)) ||
        (mask & bit(EarlyLateOutput::dv_over_rs));
    const auto expansion =
        need_local ? c.expansion(a) : ThermalScaledExpansion{};
    long double e = 0, ee = 0, dh = 0, dhe = 0;
    if (expansion.status == S::ok) {
      e = std::sqrt(expansion.a4_e2) / a2;
      ee = expansion.error_estimate /
           (2 * std::sqrt(expansion.a4_e2 - expansion.error_estimate) * a2);
      dh = scale / e;
      const long double low = scale / (e + ee),
                        high =
                            ee < e
                                ? scale / (e - ee)
                                : std::numeric_limits<long double>::infinity();
      dhe = std::max(dh - low, high - dh);
      if (z == 0) {
        if (mask & bit(EarlyLateOutput::e))
          row.outputs[static_cast<unsigned>(EarlyLateOutput::e)] = {S::ok, 1.,
                                                                    0};
        if (mask & bit(EarlyLateOutput::h_km_s_mpc)) {
          auto &h =
              row.outputs[static_cast<unsigned>(EarlyLateOutput::h_km_s_mpc)];
          if (detail::physical_representable(source_.model.h0_km_s_mpc))
            h = {S::ok, source_.model.h0_km_s_mpc, 0};
          else
            h.status = S::outside_domain;
        }
      } else {
        scalar(EarlyLateOutput::e, e, ee, 0, p.relative_tolerance);
        scalar(EarlyLateOutput::h_km_s_mpc, e * source_.model.h0_km_s_mpc,
               ee * source_.model.h0_km_s_mpc, 0, p.relative_tolerance);
      }
      scalar(EarlyLateOutput::dh_mpc, dh, dhe, p.absolute_tolerance_mpc,
             p.relative_tolerance);
    } else if (need_local)
      fail(expansion.status, mask & ~integrated);
    long double dm = 0, dme = 0, dv = 0, dve = 0;
    S distance_status = S::ok;
    if ((mask & integrated) && z > 0) {
      // Distance-only requests do not require the point's E/H projection; the
      // scaled expansion evaluations within the integral own their
      // dependencies.
      if (c.status != S::ok)
        c.status = S::ok;
      c.dependency_relative = 0;
      const auto integral = integrate(c, scale * x,
                                      (mask & bit(EarlyLateOutput::dl_mpc))
                                          ? 1 + static_cast<long double>(z)
                                          : 1);
      distance_status = integral.status;
      dm = integral.value;
      dme = integral.error;
    }
    if (distance_status == S::ok) {
      scalar(EarlyLateOutput::dm_mpc, dm, dme, p.absolute_tolerance_mpc,
             p.relative_tolerance);
      scalar(EarlyLateOutput::dl_mpc, (1 + static_cast<long double>(z)) * dm,
             (1 + static_cast<long double>(z)) * dme, p.absolute_tolerance_mpc,
             p.relative_tolerance);
      if (expansion.status == S::ok) {
        dv = z == 0 ? 0 : std::cbrt(dm * dm * z * dh);
        const long double lower =
            std::cbrt(std::max(0.L, dm - dme) * std::max(0.L, dm - dme) * z *
                      std::max(0.L, dh - dhe));
        const long double upper =
            std::cbrt((dm + dme) * (dm + dme) * z * (dh + dhe));
        dve = std::max(dv - lower, upper - dv);
        scalar(EarlyLateOutput::dv_mpc, dv, dve, p.absolute_tolerance_mpc,
               p.relative_tolerance);
      } else
        fail(expansion.status,
             bit(EarlyLateOutput::dv_mpc) | bit(EarlyLateOutput::dv_over_rs));
    } else
      fail(distance_status, integrated);
    if (mask & ratios) {
      if (!out.ruler || out.ruler->status != S::ok || !out.ruler->value)
        fail(out.ruler ? out.ruler->status : S::invalid_input, ratios);
      else {
        const long double rs = *out.ruler->value,
                          rse = out.ruler->error_estimate;
        auto ratio = [&](EarlyLateOutput id, long double numerator,
                         long double error) {
          if (!(rs > rse)) {
            fail(S::conditioning_budget_exceeded, bit(id));
            return;
          }
          const long double value = numerator / rs;
          const long double low = std::max(0.L, numerator - error) / (rs + rse),
                            high = (numerator + error) / (rs - rse);
          scalar(id, value, std::max(value - low, high - value),
                 p.absolute_tolerance_ratio, p.relative_tolerance_ratio);
        };
        if (expansion.status == S::ok)
          ratio(EarlyLateOutput::dh_over_rs, dh, dhe);
        if (distance_status == S::ok) {
          ratio(EarlyLateOutput::dm_over_rs, dm, dme);
          if (expansion.status == S::ok)
            ratio(EarlyLateOutput::dv_over_rs, dv, dve);
        }
      }
    }
    row.callbacks = work.total() - start;
  }
  out.outer_callbacks = work.outer;
  out.momentum_callbacks = work.momentum;
  out.callbacks = work.total();
  return out;
}
} // namespace irred::cosmology
