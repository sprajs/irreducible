#include "irred/rsd.hpp"
#include "bao_density_projection.hpp"
#include "payload_accounting.hpp"
#include <algorithm>
#include <cfenv>
#include <cmath>
#include <limits>
#include <utility>
namespace irred::rsd {
namespace {
using S = numerics::Status;
using D = statistics::DensityStatus;
bool arithmetic() {
  return std::fegetround() == FE_TONEAREST &&
         std::numeric_limits<long double>::digits >= 64 &&
         std::numeric_limits<long double>::max_exponent >= 16384;
}
bool known(numerics::Arithmetic a) {
  return a == numerics::Arithmetic::binary64_legacy_v1 ||
         a == numerics::Arithmetic::longdouble_cpu_v1;
}
void state(OutputState &s, S cause) {
  s.availability = cause == S::ok ? cosmology::Availability::available
                                  : cosmology::Availability::failed;
  s.numerical_status = cause;
}
auto source_strings(const DensitySource &s) {
  return std::array{&s.table_identity, &s.covariance_identity,
                    &s.ordering_provenance, &s.calibration_provenance,
                    &s.dependence_provenance};
}
void copied_string(irred::detail::PayloadAccounting &b, const std::string &s,
                   size_t count = 1) {
  if (s.size() == SIZE_MAX)
    b.add(SIZE_MAX, 2);
  else
    b.add(std::max(size_t(32), s.size() + 1), count);
}
statistics::Metadata metadata(const DensitySource &s) {
  statistics::Metadata m;
  m.ordered_ids = s.ordered_ids;
  m.measure = "product of dimensionless synthetic f_sigma8 coordinates";
  m.table_identity = s.table_identity;
  m.uncertainty_identity = s.covariance_identity;
  m.ordering_provenance = s.ordering_provenance;
  m.calibration_provenance = s.calibration_provenance;
  m.dependence_provenance = s.dependence_provenance;
  m.source_semantics = "synthetic linear pressureless total-matter f_sigma8; "
                       "supplied fixed amplitude; comoving top-hat 8 h^-1 Mpc";
  m.input_matrix_convention = "full source covariance in declared exact row "
                              "order; dimensionless f_sigma8 squared";
  return m;
}
} // namespace
std::optional<size_t>
retained_source_payload_bound(const DensitySource &s) noexcept {
  irred::detail::PayloadAccounting b(0);
  b.vector(s.scale_factors);
  b.vector(s.observed);
  b.strings(s.ordered_ids);
  b.strings(s.event_ids);
  b.strings(s.covariance_axis_ids);
  for (const auto *x : source_strings(s))
    b.string(*x);
  return b.result();
}
std::optional<size_t>
preparation_payload_bound(const DensityInput &input,
                          numerics::Arithmetic a) noexcept {
  const size_t n = input.source.scale_factors.size();
  if (!known(a) || (n && n > SIZE_MAX / n))
    return {};
  const auto source = retained_source_payload_bound(input.source);
  if (!source)
    return {};
  irred::detail::PayloadAccounting b(
      sizeof(PreparedDensity) + sizeof(DensityInput) +
      sizeof(statistics::Gaussian) + 2 * sizeof(numerics::Factorization));
  b.add(*source, 1);
  b.vector(input.covariance);
  b.add(n, sizeof(std::string) + 128); // metadata IDs and validation set
  for (const auto &id : input.source.ordered_ids)
    copied_string(b, id, 2);
  for (const auto *x : source_strings(input.source))
    copied_string(b, *x);
  // Retained Gaussian plus factorization-local preparation, matching the
  // existing Gaussian covariance preparation envelope.
  b.add(n * n, 4 * sizeof(double) + 2 * sizeof(long double));
  b.add(n, 16 * sizeof(double) + 16 * sizeof(long double));
  b.add(1, 4096);
  return b.result();
}
std::optional<size_t> density_payload_bound(size_t m, size_t n, size_t strings,
                                            unsigned requested) noexcept {
  if (!requested || (requested & ~7u) || (n && m > SIZE_MAX / n))
    return {};
  irred::detail::PayloadAccounting b(sizeof(DensityBatch));
  b.add(m, sizeof(DensitySlot));
  b.add(strings, 1);
  b.add(m * n, (bool(requested & 2u) + bool(requested & 4u)) *
                   sizeof(cosmology::GrowthAmplitudeValue));
  const auto provider = cosmology::growth_amplitude_payload_bound(n);
  if (!provider)
    return {};
  b.add(*provider, 1);
  b.add(n, 12 * sizeof(double) + 16 * sizeof(long double));
  return b.result();
}
PreparedDensity::PreparedDensity(PreparedDensity &&other) noexcept
    : source_(std::move(other.source_)), gaussian_(std::move(other.gaussian_)),
      status_(std::exchange(other.status_, D::invalid_input)),
      numerical_status_(
          std::exchange(other.numerical_status_, S::invalid_input)) {}
PreparedDensity &PreparedDensity::operator=(PreparedDensity &&other) noexcept {
  if (this != &other) {
    source_ = std::move(other.source_);
    gaussian_ = std::move(other.gaussian_);
    status_ = std::exchange(other.status_, D::invalid_input);
    numerical_status_ =
        std::exchange(other.numerical_status_, S::invalid_input);
  }
  return *this;
}
std::optional<size_t> PreparedDensity::retained_payload_bound() const noexcept {
  irred::detail::PayloadAccounting b(sizeof(*this));
  const auto source = retained_source_payload_bound(source_);
  if (!source)
    return {};
  b.add(*source, 1);
  b.embedded(gaussian_.retained_payload_bound(), sizeof(gaussian_));
  return b.result();
}
PreparedDensity prepare_density(DensityInput input, PreparationPolicy p) {
  PreparedDensity out;
  if (!known(p.arithmetic) || !std::isfinite(p.maximum_forward_sensitivity) ||
      p.maximum_forward_sensitivity <= 0)
    return out;
  if (!arithmetic()) {
    out.status_ = D::unsupported_domain;
    out.numerical_status_ = S::outside_domain;
    return out;
  }
  const auto &s = input.source;
  const size_t n = s.scale_factors.size();
  if (!n || n > p.maximum_queries || n > 65536 ||
      n > p.maximum_matrix_elements / n) {
    out.numerical_status_ = S::work_limit;
    return out;
  }
  if (s.observed.size() != n || s.ordered_ids.size() != n ||
      s.event_ids.size() != n || s.covariance_axis_ids != s.ordered_ids ||
      input.covariance.size() != n * n)
    return out;
  const auto peak = preparation_payload_bound(input, p.arithmetic);
  if (!peak || *peak > p.maximum_native_bytes ||
      *peak > size_t(1024) * 1024 * 1024) {
    out.numerical_status_ = S::work_limit;
    return out;
  }
  irred::detail::PayloadAccounting strings(0);
  for (const auto *x : source_strings(s)) {
    if (x->empty())
      return out;
    strings.add(x->size(), 1);
  }
  for (const auto *ids : {&s.ordered_ids, &s.event_ids, &s.covariance_axis_ids})
    for (const auto &id : *ids) {
      if (id.empty())
        return out;
      strings.add(id.size(), 1);
    }
  const auto chars = strings.result();
  if (!chars || *chars > p.maximum_string_bytes) {
    out.numerical_status_ = S::work_limit;
    return out;
  }
  if (s.role != RowRole::synthetic_control ||
      s.covariance_unit != CovarianceUnit::dimensionless_f_sigma8_squared ||
      s.amplitude_convention !=
          cosmology::Sigma8Convention::
              linear_pressureless_total_matter_top_hat_8_over_h_mpc)
    return out;
  for (size_t j = 0; j < n; ++j) {
    if (!std::isfinite(s.scale_factors[j]) || s.scale_factors[j] < 1e-8 ||
        s.scale_factors[j] > 1 || !std::isfinite(s.observed[j]))
      return out;
    for (size_t k = 0; k < j; ++k)
      if (s.ordered_ids[j] == s.ordered_ids[k])
        return out;
  }
  out.gaussian_ = statistics::prepare_gaussian(
      input.covariance, statistics::MatrixKind::covariance, metadata(s),
      p.maximum_matrix_elements, p.maximum_forward_sensitivity, p.arithmetic);
  out.status_ = out.gaussian_.status();
  out.numerical_status_ = out.gaussian_.numerical_status();
  out.source_ = std::move(input.source);
  // The disposable input matrix dies here; the Gaussian retains its required
  // covariance/factor storage once for all later model evaluations.
  return out;
}
DensityBatch PreparedDensity::evaluate(std::span<const ModelPoint> models,
                                       DensityPolicy p) const {
  DensityBatch out;
  if (status_ != D::finite) {
    out.status = status_;
    out.numerical_status = numerical_status_;
    return out;
  }
  if (!arithmetic()) {
    out.status = D::unsupported_domain;
    out.numerical_status = S::outside_domain;
    return out;
  }
  if (!p.requested || (p.requested & ~7u) ||
      !std::isfinite(p.maximum_forward_sensitivity) ||
      p.maximum_forward_sensitivity <= 0 ||
      !std::isfinite(p.maximum_projection_log_density_error) ||
      p.maximum_projection_log_density_error <= 0)
    return out;
  if (p.arithmetic != gaussian_.arithmetic()) {
    out.status = D::incompatible_metadata;
    return out;
  }
  const size_t n = source_.scale_factors.size();
  if (models.size() > p.maximum_models || n > p.maximum_queries) {
    out.numerical_status = S::work_limit;
    return out;
  }
  irred::detail::PayloadAccounting chars(0), copies(0);
  for (const auto &point : models) {
    chars.add(point.amplitude.amplitude_identity.size(), 1);
    chars.add(point.amplitude.amplitude_provenance.size(), 1);
    copied_string(copies, point.amplitude.amplitude_identity);
    copied_string(copies, point.amplitude.amplitude_provenance);
  }
  const auto strings = chars.result(), copied = copies.result();
  const auto peak =
      copied ? density_payload_bound(models.size(), n, *copied, p.requested)
             : std::nullopt;
  if (!strings || *strings > p.maximum_string_bytes || !peak ||
      *peak > p.maximum_native_bytes || *peak > size_t(1024) * 1024 * 1024) {
    out.numerical_status = S::work_limit;
    return out;
  }
  out.status = D::finite;
  out.numerical_status = S::ok;
  out.slots.reserve(models.size());
  long double inverse_norm = 0;
  if (p.requested & 1u)
    inverse_norm = bao::detail::inverse_covariance_norm_estimate(
        gaussian_.factor_, gaussian_.covariance_, n);
  size_t remaining = p.maximum_total_callbacks;
  for (const auto &point : models) {
    out.slots.emplace_back(point);
    auto &slot = out.slots.back();
    auto local = p.predictions;
    local.growth.maximum_total_callbacks =
        std::min(local.growth.maximum_total_callbacks, remaining);
    auto growth = cosmology::prepare_gr_growth(
        cosmology::prepare(point.expansion, point.geometry));
    auto predicted = cosmology::evaluate_growth_amplitude(
        growth, point.amplitude, source_.scale_factors,
        cosmology::amplitude_f_sigma8, local);
    slot.callbacks = predicted.callbacks;
    out.callbacks += slot.callbacks;
    remaining -= slot.callbacks;
    slot.reference = predicted.reference;
    S cause = predicted.status;
    std::vector<double> mu, residual;
    std::vector<long double> eps;
    if (p.requested & 5u)
      mu.reserve(n);
    if (p.requested & 5u)
      residual.reserve(n);
    if (p.requested & 1u)
      eps.reserve(n);
    if (p.requested & 2u)
      slot.predictions.reserve(n);
    if (p.requested & 4u)
      slot.residuals.reserve(n);
    if (cause == S::ok) {
      for (size_t j = 0; j < n; ++j) {
        const auto &v = predicted.rows[j].f_sigma8;
        if (p.requested & 2u)
          slot.predictions.push_back(v);
        if (v.status != S::ok || !v.value) {
          if (cause == S::ok)
            cause = v.status == S::ok ? S::invalid_input : v.status;
        } else {
          if (p.requested & 5u)
            mu.push_back(*v.value);
          if (p.requested & 1u)
            eps.push_back(v.absolute_error_estimate);
        }
        if (p.requested & 4u) {
          auto r = v;
          r.value.reset();
          if (v.value && v.status == S::ok) {
            const long double wide =
                (long double)source_.observed[j] - *v.value;
            const double rounded = static_cast<double>(wide);
            if (!std::isfinite(rounded) ||
                (wide != 0 && (rounded == 0 || !std::isnormal(rounded)))) {
              r.availability = cosmology::Availability::failed;
              r.status = S::overflow;
            } else {
              r.value = rounded;
              r.absolute_error_estimate += static_cast<double>(
                  std::abs(wide - rounded) +
                  std::numeric_limits<long double>::epsilon() *
                      (std::abs((long double)source_.observed[j]) +
                       std::abs((long double)*v.value)));
            }
          }
          slot.residuals.push_back(r);
        }
      }
    }
    if (p.requested & 2u)
      state(slot.predictions_state, cause);
    if (cause == S::ok && (p.requested & 5u))
      cause = bao::detail::subtract_observations(source_.observed, mu, residual,
                                                 eps);
    if (p.requested & 4u)
      state(slot.residuals_state, cause);
    if (cause == S::ok && (p.requested & 1u))
      cause = bao::detail::project_density(
          gaussian_.factor_, gaussian_, source_.ordered_ids, residual, eps,
          inverse_norm, p.maximum_forward_sensitivity,
          p.maximum_projection_log_density_error,
          slot.projection_log_density_error_estimate, slot.result);
    if (cause != S::ok)
      slot.result.reset();
    if (p.requested & 1u)
      state(slot.density_state, cause);
    slot.numerical_status = cause;
  }
  return out;
}
} // namespace irred::rsd
