#include "irred/growth_amplitude.hpp"
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
S query_status(double a) {
  if (!std::isfinite(a))
    return S::nonfinite_input;
  return a < 1e-8 || a > 1 ? S::outside_domain : S::ok;
}
void fail(GrowthAmplitudeValue &v, S cause) {
  v.availability = Availability::failed;
  v.status = cause;
}
void accept(GrowthAmplitudeValue &out, long double v, long double error,
            const GrowthAmplitudePolicy &p, bool exact = false) {
  const double rounded = static_cast<double>(v);
  error += std::abs(v - rounded);
  if (!exact)
    error += 16 * std::numeric_limits<long double>::epsilon() * std::abs(v);
  const double reported = static_cast<double>(error);
  if (!std::isfinite(v) || !std::isfinite(rounded) || !std::isfinite(error) ||
      error < 0 || !std::isfinite(reported)) {
    fail(out, S::overflow);
    return;
  }
  if ((v != 0 && (rounded == 0 || !std::isnormal(rounded))) ||
      (error != 0 && (reported == 0 || !std::isnormal(reported)))) {
    fail(out, S::outside_domain);
    return;
  }
  out.absolute_error_estimate = reported;
  if (error > p.absolute_tolerance + p.relative_tolerance * std::abs(v)) {
    fail(out, S::conditioning_budget_exceeded);
    return;
  }
  out.availability = Availability::available;
  out.status = S::ok;
  out.value = rounded;
}
void compose(GrowthAmplitudeValue &out, double amplitude, const GrowthValue &d,
             const GrowthValue *f, const GrowthValue &reference,
             const GrowthAmplitudePolicy &p) {
  for (const auto *x : {&d, f, &reference}) {
    if (!x)
      continue;
    if (x->status != S::ok || !x->value) {
      fail(out, x->status == S::ok ? S::invalid_input : x->status);
      return;
    }
    if (!(*x->value > x->absolute_error_estimate)) {
      fail(out, S::conditioning_budget_exceeded);
      return;
    }
  }
  const long double D = *d.value, de = d.absolute_error_estimate,
                    B = *reference.value,
                    be = reference.absolute_error_estimate,
                    F = f ? *f->value : 1,
                    fe = f ? f->absolute_error_estimate : 0, A = amplitude;
  const long double v = A * D * F / B,
                    lower = A * (D - de) * (F - fe) / (B + be),
                    upper = A * (D + de) * (F + fe) / (B - be),
                    error = std::max(std::abs(v - lower), std::abs(upper - v)) +
                            32 * std::numeric_limits<long double>::epsilon() *
                                (std::abs(lower) + std::abs(upper));
  accept(out, v, error, p);
}
} // namespace
S fixed_sigma8_status(const FixedSigma8Source &s) noexcept {
  if (s.convention !=
          Sigma8Convention::
              linear_pressureless_total_matter_top_hat_8_over_h_mpc ||
      s.treatment != AmplitudeTreatment::fixed_supplied ||
      s.amplitude_identity.empty() || s.amplitude_provenance.empty())
    return S::invalid_input;
  if (!std::isfinite(s.sigma8_ref))
    return S::nonfinite_input;
  if (s.sigma8_ref < 0 || (s.sigma8_ref != 0 && !std::isnormal(s.sigma8_ref)))
    return S::outside_domain;
  return query_status(s.reference_scale_factor);
}
std::optional<size_t> growth_amplitude_payload_bound(size_t n) noexcept {
  irred::detail::PayloadAccounting b(sizeof(GrowthAmplitudeBatch));
  b.add(n, sizeof(GrowthAmplitudeRow) + sizeof(double) + sizeof(size_t));
  const auto query = growth_payload_bound(n),
             reference = growth_payload_bound(1);
  if (!query || !reference)
    return {};
  b.add(*query, 1);
  b.add(*reference, 1);
  return b.result();
}
GrowthAmplitudeBatch evaluate_growth_amplitude(const GRGrowth &growth,
                                               const FixedSigma8Source &source,
                                               std::span<const double> a,
                                               unsigned requested,
                                               GrowthAmplitudePolicy p) {
  GrowthAmplitudeBatch out;
  if (!arithmetic() || !requested || (requested & ~3u) ||
      !std::isfinite(p.absolute_tolerance) || p.absolute_tolerance < 0 ||
      !std::isfinite(p.relative_tolerance) || p.relative_tolerance <= 0 ||
      !std::isfinite(p.growth.relative_tolerance) ||
      p.growth.relative_tolerance <= 0 || p.growth.maximum_depth > 60)
    return out;
  irred::detail::PayloadAccounting strings(0);
  strings.add(source.amplitude_identity.size(), 1);
  strings.add(source.amplitude_provenance.size(), 1);
  const auto string_bytes = strings.result(),
             bytes = growth_amplitude_payload_bound(a.size());
  if (!string_bytes || *string_bytes > p.maximum_string_bytes ||
      a.size() > 65536 || a.size() > p.growth.maximum_points || !bytes ||
      *bytes > p.maximum_native_bytes || *bytes > size_t(1024) * 1024 * 1024) {
    out.status = S::work_limit;
    return out;
  }
  out.status = fixed_sigma8_status(source);
  if (out.status != S::ok)
    return out;
  if (growth.status() != S::ok) {
    out.status = growth.status();
    return out;
  }
  out.requested = requested;
  out.rows.reserve(a.size());
  std::vector<double> queries;
  std::vector<size_t> indices;
  queries.reserve(a.size());
  indices.reserve(a.size());
  bool reference_f = false;
  for (size_t j = 0; j < a.size(); ++j) {
    out.rows.push_back({a[j], {}, {}, 0});
    auto &row = out.rows.back();
    auto cause = query_status(a[j]);
    if (cause != S::ok) {
      if (requested & amplitude_sigma8)
        fail(row.sigma8, cause);
      if (requested & amplitude_f_sigma8)
        fail(row.f_sigma8, cause);
    } else if (source.sigma8_ref == 0) {
      if (requested & amplitude_sigma8)
        accept(row.sigma8, 0, 0, p, true);
      if (requested & amplitude_f_sigma8)
        accept(row.f_sigma8, 0, 0, p, true);
    } else if (a[j] == source.reference_scale_factor) {
      if (requested & amplitude_sigma8)
        accept(row.sigma8, source.sigma8_ref, 0, p, true);
      reference_f |= bool(requested & amplitude_f_sigma8);
    } else {
      queries.push_back(a[j]);
      indices.push_back(j);
    }
  }
  const unsigned reference_mask =
      (queries.empty() ? 0u : growth_d) | (reference_f ? growth_f : 0u);
  size_t remaining = p.growth.maximum_total_callbacks;
  if (reference_mask) {
    auto local = p.growth;
    local.maximum_total_callbacks = remaining;
    const double ref[]{source.reference_scale_factor};
    auto computed = growth.evaluate(ref, reference_mask, local);
    out.callbacks = computed.callbacks;
    remaining -= computed.callbacks;
    if (!computed.rows.empty())
      out.reference = computed.rows[0];
    else {
      GrowthRow failed{source.reference_scale_factor, {}, {}, 0};
      failed.d.status = failed.f.status = computed.status;
      out.reference = failed;
    }
    if (reference_f) {
      for (auto &row : out.rows) {
        if (row.scale_factor != source.reference_scale_factor)
          continue;
        const auto &f = out.reference->f;
        if (f.status != S::ok || !f.value)
          fail(row.f_sigma8, f.status);
        else
          accept(row.f_sigma8, (long double)source.sigma8_ref * *f.value,
                 (long double)source.sigma8_ref * f.absolute_error_estimate, p);
      }
    }
  }
  if (!queries.empty()) {
    if (!out.reference->d.value || out.reference->d.status != S::ok) {
      for (auto j : indices) {
        if (requested & amplitude_sigma8)
          fail(out.rows[j].sigma8, out.reference->d.status);
        if (requested & amplitude_f_sigma8)
          fail(out.rows[j].f_sigma8, out.reference->d.status);
      }
    } else {
      auto local = p.growth;
      local.maximum_total_callbacks = remaining;
      auto computed = growth.evaluate(
          queries,
          growth_d | ((requested & amplitude_f_sigma8) ? growth_f : 0u), local);
      out.callbacks += computed.callbacks;
      for (size_t k = 0; k < indices.size(); ++k) {
        auto &row = out.rows[indices[k]];
        if (computed.rows.empty()) {
          if (requested & amplitude_sigma8)
            fail(row.sigma8, computed.status);
          if (requested & amplitude_f_sigma8)
            fail(row.f_sigma8, computed.status);
          continue;
        }
        const auto &v = computed.rows[k];
        row.callbacks = v.callbacks;
        if (requested & amplitude_sigma8)
          compose(row.sigma8, source.sigma8_ref, v.d, nullptr, out.reference->d,
                  p);
        if (requested & amplitude_f_sigma8)
          compose(row.f_sigma8, source.sigma8_ref, v.d, &v.f, out.reference->d,
                  p);
      }
    }
  }
  return out;
}
} // namespace irred::cosmology
