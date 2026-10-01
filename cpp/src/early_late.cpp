#include "irred/early_late.hpp"
#include "early_flat_state.hpp"
#include "flat_geometry.hpp"
#include "payload_accounting.hpp"
#include <algorithm>
#include <cfenv>
#include <cmath>
#include <limits>
namespace irred::cosmology {
namespace {
using S = numerics::Status;
constexpr unsigned index(EarlyLateOutput x) { return static_cast<unsigned>(x); }
constexpr unsigned ratios = early_late_mask(EarlyLateOutput::dm_over_rs) |
                            early_late_mask(EarlyLateOutput::dh_over_rs) |
                            early_late_mask(EarlyLateOutput::dv_over_rs);
constexpr unsigned distances = early_late_mask(EarlyLateOutput::dm_mpc) |
                               early_late_mask(EarlyLateOutput::dl_mpc) |
                               early_late_mask(EarlyLateOutput::dv_mpc) |
                               early_late_mask(EarlyLateOutput::dm_over_rs) |
                               early_late_mask(EarlyLateOutput::dv_over_rs);
struct Context {
  detail::EarlyFlatState state;
  long double xmax;
};
double integrand(double t, const void *ptr) {
  const auto &c = *static_cast<const Context *>(ptr);
  const long double a = std::exp(-c.xmax * t);
  return static_cast<double>(
      a / std::sqrt(detail::early_flat_polynomial(c.state, a)));
}
bool tolerance(double a, double r) {
  return std::isfinite(a) && std::isfinite(r) && a >= 0 && r >= 0 &&
         (a > 0 || r > 0);
}
void admit(EarlyLateValue &out, long double value, long double error,
           double absolute, double relative) {
  if (!detail::physical_representable(value) || value < 0) {
    out.status = S::outside_domain;
    return;
  }
  const double rounded = static_cast<double>(value);
  const long double floor =
      std::abs(value) * (32.L * std::numeric_limits<double>::epsilon() +
                         128.L * std::numeric_limits<long double>::epsilon());
  const long double diagnostic = error + floor + std::abs(value - rounded);
  if (!detail::representable(diagnostic) ||
      diagnostic > absolute + relative * std::abs(value)) {
    out.status = S::conditioning_budget_exceeded;
    return;
  }
  out = {S::ok, rounded, static_cast<double>(diagnostic)};
}
} // namespace
std::optional<std::size_t>
early_late_payload_bound(std::size_t n, std::size_t origins) noexcept {
  irred::detail::PayloadAccounting b(sizeof(EarlyLateBatch));
  b.add(n, 2 * sizeof(EarlyLateRow));
  // source, ruler row and transient sound batch each retain a source string;
  // conservative temporary row and output owner capacities are included.
  b.add(origins, 8);
  b.add(1, 2 * sizeof(SoundHorizonBatch) + 4 * sizeof(SoundHorizonRow));
  return b.result();
}
EarlyLateBatch evaluate_early_late(const SoundHorizonRequest &source,
                                   std::span<const double> zs, unsigned mask,
                                   EarlyLatePolicy p) {
  EarlyLateBatch out;
  if (std::numeric_limits<long double>::digits < 64 ||
      std::numeric_limits<long double>::max_exponent < 16384 ||
      std::fegetround() != FE_TONEAREST || !mask ||
      mask >= (1u << early_late_output_count) ||
      !tolerance(p.absolute_tolerance_mpc, p.relative_tolerance) ||
      !tolerance(p.absolute_tolerance_ratio, p.relative_tolerance_ratio) ||
      p.maximum_depth > 60 || zs.size() > 65536)
    return out;
  if (source.drag_origin.size() == std::numeric_limits<std::size_t>::max())
    return out;
  const auto bytes =
      early_late_payload_bound(zs.size(), source.drag_origin.size() + 1);
  if (!bytes)
    return out;
  if (zs.size() > p.maximum_points || *bytes > p.maximum_native_bytes) {
    out.status = S::work_limit;
    return out;
  }
  // Admission precedes copying origins, reserving rows, and redshift scans.
  out.source = source;
  out.requested_outputs = mask;
  out.rows.reserve(zs.size());
  out.status = S::ok;
  detail::EarlyFlatState state{};
  auto state_status = detail::prepare_early_flat(source.model, state);
  if (state_status == S::ok && !std::isfinite(source.z_drag))
    state_status = S::nonfinite_input;
  if (state_status == S::ok &&
      (source.z_drag < 0 || source.drag_origin.empty()))
    state_status = S::outside_domain;
  if ((mask & ratios) && !zs.empty() && state_status == S::ok) {
    auto sound_policy = p.sound;
    sound_policy.maximum_total_callbacks = std::min(
        sound_policy.maximum_total_callbacks, p.maximum_total_callbacks);
    const auto sound = evaluate_sound_horizon({&source, 1}, sound_policy);
    out.callbacks = sound.callbacks;
    if (!sound.rows.empty())
      out.ruler = sound.rows.front();
    else {
      out.ruler = SoundHorizonRow{};
      out.ruler->status = sound.status;
    }
  }
  const long double scale =
      detail::prepare_flat_scale(source.model.h0_km_s_mpc).distance_mpc;
  for (double z : zs) {
    out.rows.emplace_back();
    auto &row = out.rows.back();
    row.redshift = z;
    auto fail = [&](S status, unsigned selected) {
      for (unsigned i = 0; i < early_late_output_count; ++i)
        if (mask & selected & (1u << i))
          row.outputs[i].status = status;
    };
    if (state_status != S::ok) {
      fail(state_status, mask);
      continue;
    }
    if (!std::isfinite(z)) {
      fail(S::nonfinite_input, mask);
      continue;
    }
    if (z < 0) {
      fail(S::outside_domain, mask);
      continue;
    }
    const long double x = std::log1p((long double)z), a = std::exp(-x);
    const long double e =
        std::sqrt(detail::early_flat_polynomial(state, a)) / (a * a);
    const long double h = source.model.h0_km_s_mpc * e, dh = scale / e;
    auto scalar = [&](EarlyLateOutput id, long double v, long double err,
                      double abs, double rel) {
      if (mask & early_late_mask(id))
        admit(row.outputs[index(id)], v, err, abs, rel);
    };
    scalar(EarlyLateOutput::e, e, 0, 0, p.relative_tolerance);
    scalar(EarlyLateOutput::h_km_s_mpc, h, 0, 0, p.relative_tolerance);
    scalar(EarlyLateOutput::dh_mpc, dh, 0, p.absolute_tolerance_mpc,
           p.relative_tolerance);
    long double dm = 0, dm_error = 0;
    S distance_status = S::ok;
    if ((mask & distances) && z > 0) {
      const std::size_t available =
          std::min(p.maximum_callbacks_per_point,
                   p.maximum_total_callbacks - out.callbacks);
      if (available < 3)
        distance_status = S::work_limit;
      else {
        const long double prefactor = scale * x;
        const double absolute = static_cast<double>(
            std::min((long double)p.absolute_tolerance_mpc / prefactor,
                     (long double)std::numeric_limits<double>::max()));
        // Reserve projection amplification for DL and the ratio consumer.
        long double amplification = 1;
        if (mask & early_late_mask(EarlyLateOutput::dl_mpc))
          amplification = 1.L + z;
        const double requested_absolute =
            static_cast<double>((long double)absolute / amplification);
        if (requested_absolute == 0 && p.relative_tolerance == 0)
          distance_status = S::conditioning_budget_exceeded;
        else {
          Context c{state, x};
          const auto q =
              numerics::integrate(integrand, &c, 0, 1,
                                  {requested_absolute, p.relative_tolerance / 4,
                                   available, p.maximum_depth});
          row.callbacks = q.evaluations;
          out.callbacks += q.evaluations;
          distance_status = q.status;
          if (q.status == S::ok) {
            dm = prefactor * q.value;
            dm_error = prefactor *
                       (q.error_estimate +
                        std::abs((long double)std::nextafter(
                                     q.value,
                                     std::numeric_limits<double>::infinity()) -
                                 q.value) /
                            2);
          }
        }
      }
    }
    if (distance_status == S::ok) {
      scalar(EarlyLateOutput::dm_mpc, dm, dm_error, p.absolute_tolerance_mpc,
             p.relative_tolerance);
      scalar(EarlyLateOutput::dl_mpc, (1.L + z) * dm, (1.L + z) * dm_error,
             p.absolute_tolerance_mpc, p.relative_tolerance);
    } else
      fail(distance_status, distances);
    const long double dv = z == 0 ? 0 : std::cbrt(dm * dm * z * dh);
    // Positive distance uncertainties add deterministically. Linearized DV
    // propagation is guarded below by a finite-radius interval instead.
    long double dv_error = 0;
    if (z > 0 && distance_status == S::ok) {
      const long double dh_floor =
          std::abs(dh) * 64.L * std::numeric_limits<double>::epsilon();
      const long double low = std::cbrt(std::max(0.L, dm - dm_error) *
                                        std::max(0.L, dm - dm_error) * z *
                                        std::max(0.L, dh - dh_floor));
      const long double high =
          std::cbrt((dm + dm_error) * (dm + dm_error) * z * (dh + dh_floor));
      dv_error = std::max(dv - low, high - dv);
    }
    if (distance_status == S::ok)
      scalar(EarlyLateOutput::dv_mpc, dv, dv_error, p.absolute_tolerance_mpc,
             p.relative_tolerance);
    if (mask & ratios) {
      if (!out.ruler || out.ruler->status != S::ok ||
          !out.ruler->sound_horizon_mpc) {
        fail(out.ruler ? out.ruler->status : S::invalid_input, ratios);
        continue;
      }
      const long double rs = *out.ruler->sound_horizon_mpc,
                        rs_error = out.ruler->error_estimate_mpc;
      if (!(rs_error < rs)) {
        fail(S::conditioning_budget_exceeded, ratios);
        continue;
      }
      auto ratio = [&](EarlyLateOutput id, long double numerator,
                       long double err) {
        const long double v = numerator / rs;
        const long double lower =
                              std::max(0.L, numerator - err) / (rs + rs_error),
                          upper = (numerator + err) / (rs - rs_error);
        scalar(id, v, std::max(v - lower, upper - v),
               p.absolute_tolerance_ratio, p.relative_tolerance_ratio);
      };
      ratio(EarlyLateOutput::dh_over_rs, dh,
            std::abs(dh) * 64.L * std::numeric_limits<double>::epsilon());
      if (distance_status == S::ok) {
        ratio(EarlyLateOutput::dm_over_rs, dm, dm_error);
        ratio(EarlyLateOutput::dv_over_rs, dv, dv_error);
      }
    }
  }
  return out;
}
} // namespace irred::cosmology
