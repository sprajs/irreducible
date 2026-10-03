#include "irred/abi.h"
#include "irred/gaussian_posterior.hpp"
#include "payload_accounting.hpp"
#include <algorithm>
#include <cmath>
#include <memory>
#include <new>
#include <optional>
#include <string_view>
#include <utility>
namespace {
namespace s = irred::statistics;
namespace n = irred::numerics;
constexpr size_t elements_cap = 1000000, cases_cap = 65536,
                 string_cap = 1048576, bytes_cap = 268435456,
                 work_cap = 100000000;
template <class T> bool aligned(const T *p) {
  return p && reinterpret_cast<uintptr_t>(p) % alignof(T) == 0;
}
template <class T> bool header(const T *p) {
  return aligned(p) && p->struct_size == sizeof(T) &&
         p->abi_version == IRRED_ABI_VERSION;
}
bool buffer(const irred_f64_buffer &b) {
  return header(&b) && b.element_type == 2 && !b.reserved &&
         b.length <= elements_cap &&
         b.byte_length == b.length * sizeof(double) &&
         (!b.length || aligned(b.data));
}
bool text(const irred_bytes &b, size_t &left) {
  if (!b.length || b.length > 256 || b.length > left || !b.data)
    return false;
  left -= b.length;
  return true;
}
bool strings(const irred_strings &b, size_t &left) {
  if (b.length > cases_cap || b.byte_length != b.length * sizeof(irred_bytes) ||
      (b.length && !aligned(b.data)))
    return false;
  for (size_t i = 0; i < b.length; ++i)
    if (!text(b.data[i], left))
      return false;
  return true;
}
std::string copy(const irred_bytes &b) {
  return {reinterpret_cast<const char *>(b.data),
          static_cast<size_t>(b.length)};
}
std::vector<std::string> copy(const irred_strings &b) {
  std::vector<std::string> out(b.length);
  for (size_t i = 0; i < b.length; ++i)
    out[i] = copy(b.data[i]);
  return out;
}
bool equal(const irred_strings &a, const irred_strings &b) {
  if (a.length != b.length)
    return false;
  for (size_t i = 0; i < a.length; ++i)
    if (a.data[i].length != b.data[i].length ||
        !std::equal(a.data[i].data, a.data[i].data + a.data[i].length,
                    b.data[i].data))
      return false;
  return true;
}
bool unique(const irred_strings &a) {
  std::vector<std::string_view> keys(a.length);
  for (size_t i = 0; i < a.length; ++i)
    keys[i] = {reinterpret_cast<const char *>(a.data[i].data),
               static_cast<size_t>(a.data[i].length)};
  std::sort(keys.begin(), keys.end());
  return std::adjacent_find(keys.begin(), keys.end()) == keys.end();
}
irred_f64_buffer view(std::span<const double> v) {
  return {sizeof(irred_f64_buffer),
          IRRED_ABI_VERSION,
          2,
          0,
          v.data(),
          v.size(),
          v.size_bytes()};
}
irred_bytes bytes(std::string_view s) {
  return {reinterpret_cast<const uint8_t *>(s.data()), s.size()};
}
std::span<const double> span(const irred_f64_buffer &b) {
  return {b.data, static_cast<size_t>(b.length)};
}
} // namespace
struct irred_gaussian_posterior_result {
  // Allocation-free failure constructor: no engaged Gaussian metadata strings.
  std::optional<s::GaussianPosterior> posterior;
  std::vector<s::PosteriorMean> means;
  std::vector<irred_gaussian_posterior_row> rows;
  s::DensityStatus status = s::DensityStatus::invalid_input;
  n::Status numerical_status = n::Status::invalid_input;
  uint32_t phase = 0; // 0 admission, 1 noise, 2 posterior, 3 conditioning
  size_t work = 0, elements = 0, peak = 0, retained = 0;
  const char *arithmetic = "unresolved";
};
extern "C" uint32_t
irred_gaussian_posterior_evaluate(const irred_gaussian_posterior_batch *b,
                                  const irred_gaussian_posterior_policy *q,
                                  irred_gaussian_posterior_result **output) {
  if (!aligned(output))
    return IRRED_INVALID_INPUT;
  *output = nullptr;
  if (!aligned(b) || !aligned(q))
    return IRRED_INVALID_INPUT;
  if (b->abi_version != IRRED_ABI_VERSION ||
      q->abi_version != IRRED_ABI_VERSION)
    return IRRED_ABI_MISMATCH;
  if (!header(b) || !header(q) || q->reserved || q->arithmetic > 1 ||
      q->maximum_elements > elements_cap || q->maximum_cases > cases_cap ||
      q->maximum_string_bytes > string_cap ||
      q->maximum_native_bytes > bytes_cap || q->maximum_work_units > work_cap ||
      !std::isfinite(q->maximum_forward_sensitivity) ||
      q->maximum_forward_sensitivity <= 0 || b->case_count > cases_cap)
    return IRRED_INVALID_INPUT;
  for (const auto *v : {&b->noise_covariance, &b->design, &b->prior_mean,
                        &b->prior_covariance, &b->conditioning_vectors})
    if (!buffer(*v))
      return IRRED_INVALID_INPUT;
  size_t left = string_cap, metadata_slots = 0;
  for (const auto *v :
       {&b->noise_row_ids, &b->noise_event_ids, &b->design_row_ids,
        &b->design_parameter_ids, &b->design_parameter_units, &b->column_units,
        &b->prior_parameter_ids, &b->prior_parameter_units,
        &b->shared_nuisance_ids, &b->conditioning_row_ids,
        &b->conditioning_event_ids, &b->case_ids}) {
    if (!strings(*v, left) ||
        !irred::detail::checked_payload_add(metadata_slots, v->length, 1))
      return IRRED_INVALID_INPUT;
  }
  for (const auto *v :
       {&b->residual_unit, &b->noise_identity, &b->calibration_identity,
        &b->noise_dependence_identity, &b->ordering_provenance,
        &b->design_identity, &b->prior_identity, &b->parameter_measure,
        &b->prior_dependence_identity})
    if (!text(*v, left))
      return IRRED_INVALID_INPUT;
  const size_t nrow = b->noise_row_ids.length,
               p = b->prior_parameter_ids.length, count = b->case_count;
  // Hard element bounds make these products representable even on 32-bit
  // size_t.
  if (!nrow || p < 2 || !count || nrow > elements_cap / nrow ||
      p > elements_cap / p || nrow > elements_cap / p ||
      count > elements_cap / nrow || count > elements_cap / p ||
      b->noise_covariance.length != nrow * nrow ||
      b->design.length != nrow * p || b->prior_mean.length != p ||
      b->prior_covariance.length != p * p ||
      b->conditioning_vectors.length != count * nrow ||
      b->case_ids.length != count)
    return IRRED_INVALID_INPUT;
  if (q->maximum_native_bytes < sizeof(irred_gaussian_posterior_result))
    return IRRED_GAUSSIAN_POSTERIOR_QUOTA_REFUSED;
  try {
    auto result = std::make_unique<irred_gaussian_posterior_result>();
    result->peak = result->retained = sizeof(*result);
    auto refuse = [&](s::DensityStatus status, n::Status numerical) {
      if (status != s::DensityStatus::finite) {
        result->posterior.reset();
        std::vector<s::PosteriorMean>().swap(result->means);
        std::vector<irred_gaussian_posterior_row>().swap(result->rows);
        result->retained = sizeof(*result);
      }
      result->status = status;
      result->numerical_status = numerical;
      *output = result.release();
      return IRRED_OK;
    };
    result->arithmetic =
        q->arithmetic ? "F02/longdouble-cpu/v1" : "F02/binary64-legacy/v1";
    irred::detail::PayloadAccounting numeric(0), work(0);
    numeric.add(nrow * nrow, 1);
    numeric.add(nrow * p, 1);
    numeric.add(p, 1);
    numeric.add(p * p, 2);
    numeric.add(count * nrow, 1);
    numeric.add(count * p, 2);
    work.add(nrow * nrow, nrow);
    work.add(nrow * nrow, p);
    work.add(nrow * p, p);
    work.add(p * p, 8 * p);
    work.add(count, nrow * nrow + nrow * p + p * p);
    const auto units = work.result(), elems = numeric.result();
    if (!units || *units > SIZE_MAX / 32 || !elems)
      return refuse(s::DensityStatus::invalid_input, n::Status::work_limit);
    result->work = *units * 32;
    result->elements = *elems;
    if (count > q->maximum_cases || *elems > q->maximum_elements ||
        result->work > q->maximum_work_units ||
        string_cap - left > q->maximum_string_bytes)
      return refuse(s::DensityStatus::invalid_input, n::Status::work_limit);
    if (b->source_semantics != 1 || b->noise_independence_declared != 1 ||
        b->noise_event_ids.length != nrow ||
        !equal(b->noise_row_ids, b->design_row_ids) ||
        !equal(b->noise_row_ids, b->conditioning_row_ids) ||
        !equal(b->noise_event_ids, b->conditioning_event_ids) ||
        !equal(b->prior_parameter_ids, b->design_parameter_ids) ||
        !equal(b->prior_parameter_units, b->design_parameter_units) ||
        b->prior_parameter_units.length != p || b->column_units.length != p)
      return refuse(s::DensityStatus::incompatible_metadata,
                    n::Status::invalid_input);
    // Bound conversion metadata and setup scratch prospectively before copying.
    // At most two retained/candidate metadata copies; each string's
    // conservative capacity allowance is 256+1 and each original descriptor is
    // charged.
    size_t envelope = sizeof(*result);
    if (!irred::detail::checked_payload_add(envelope, metadata_slots + 9,
                                            2 * (sizeof(std::string) + 257)) ||
        !irred::detail::checked_payload_add(envelope, nrow * nrow,
                                            64 * sizeof(long double)) ||
        !irred::detail::checked_payload_add(envelope, nrow * p,
                                            32 * sizeof(long double)) ||
        !irred::detail::checked_payload_add(envelope, p * p,
                                            64 * sizeof(long double)) ||
        !irred::detail::checked_payload_add(envelope, nrow + p,
                                            64 * sizeof(long double)) ||
        envelope > q->maximum_native_bytes)
      return refuse(s::DensityStatus::invalid_input, n::Status::work_limit);
    result->peak = envelope;
    if (!unique(b->noise_row_ids) || !unique(b->case_ids))
      return refuse(s::DensityStatus::incompatible_metadata,
                    n::Status::invalid_input);
    s::Metadata metadata;
    metadata.ordered_ids = copy(b->noise_row_ids);
    metadata.measure = "product d(" + copy(b->residual_unit) + ")";
    metadata.table_identity = "synthetic supplied conditioning family";
    metadata.uncertainty_identity = copy(b->noise_identity);
    metadata.ordering_provenance = copy(b->ordering_provenance);
    metadata.calibration_provenance = copy(b->calibration_identity);
    metadata.dependence_provenance = copy(b->noise_dependence_identity);
    metadata.source_semantics = "synthetic controls";
    metadata.input_matrix_convention = "covariance";
    s::ParameterPrior prior;
    prior.ordered_parameter_ids = copy(b->prior_parameter_ids);
    prior.parameter_units = copy(b->prior_parameter_units);
    prior.shared_nuisance_ids = copy(b->shared_nuisance_ids);
    prior.mean.assign(b->prior_mean.data, b->prior_mean.data + p);
    prior.covariance.assign(b->prior_covariance.data,
                            b->prior_covariance.data + p * p);
    prior.prior_identity = copy(b->prior_identity);
    prior.design_identity = copy(b->design_identity);
    prior.residual_unit = copy(b->residual_unit);
    prior.parameter_measure = copy(b->parameter_measure);
    prior.dependence_identity = copy(b->prior_dependence_identity);
    prior.noise_independence_declared = true;
    for (size_t j = 0; j < p; ++j)
      if (copy(b->column_units.data[j]) !=
          prior.residual_unit + "/" + prior.parameter_units[j])
        return refuse(s::DensityStatus::incompatible_metadata,
                      n::Status::invalid_input);
    const auto arithmetic = q->arithmetic ? n::Arithmetic::longdouble_cpu_v1
                                          : n::Arithmetic::binary64_legacy_v1;
    const auto source_peak = s::gaussian_preparation_payload_bound(
        nrow, s::MatrixKind::covariance, arithmetic, metadata);
    irred::detail::PayloadAccounting conversion(sizeof(prior));
    conversion.strings(prior.ordered_parameter_ids);
    conversion.strings(prior.parameter_units);
    conversion.strings(prior.shared_nuisance_ids);
    conversion.vector(prior.mean);
    conversion.vector(prior.covariance);
    for (const auto *v :
         {&prior.prior_identity, &prior.design_identity, &prior.residual_unit,
          &prior.parameter_measure, &prior.dependence_identity})
      conversion.string(*v);
    const auto conversion_bytes = conversion.result();
    irred::detail::PayloadAccounting source_phase(sizeof(*result));
    if (source_peak)
      source_phase.add(*source_peak, 1);
    if (conversion_bytes)
      source_phase.add(*conversion_bytes, 1);
    const auto source_phase_peak = source_phase.result();
    if (!source_peak || !conversion_bytes || !source_phase_peak ||
        *source_phase_peak > q->maximum_native_bytes)
      return refuse(s::DensityStatus::invalid_input, n::Status::work_limit);
    result->peak = std::max(result->peak, *source_phase_peak);
    result->phase = 1;
    auto source = s::prepare_gaussian(
        span(b->noise_covariance), s::MatrixKind::covariance,
        std::move(metadata), q->maximum_elements,
        q->maximum_forward_sensitivity, arithmetic);
    if (source.status() != s::DensityStatus::finite)
      return refuse(source.status(), source.numerical_status());
    const auto posterior_peak =
        s::GaussianPosterior::preparation_payload_bound(source, prior);
    if (!posterior_peak ||
        *posterior_peak > q->maximum_native_bytes - sizeof(*result))
      return refuse(s::DensityStatus::invalid_input, n::Status::work_limit);
    result->peak = std::max(result->peak, *posterior_peak + sizeof(*result));
    s::PosteriorPolicy policy{static_cast<size_t>(q->maximum_elements),
                              static_cast<size_t>(q->maximum_native_bytes) -
                                  sizeof(*result),
                              static_cast<size_t>(q->maximum_work_units),
                              q->maximum_forward_sensitivity};
    result->phase = 2;
    result->posterior.emplace(s::GaussianPosterior::prepare(
        std::move(source), span(b->design), source.metadata().ordered_ids,
        std::move(prior), policy));
    if (result->posterior->status() != s::DensityStatus::finite)
      return refuse(result->posterior->status(),
                    result->posterior->numerical_status());
    const auto retained = result->posterior->retained_payload_bound(),
               scratch = result->posterior->evaluation_payload_bound();
    irred::detail::PayloadAccounting evaluation(sizeof(*result));
    evaluation.embedded(retained, sizeof(s::GaussianPosterior));
    evaluation.add(count, sizeof(s::PosteriorMean) +
                              sizeof(irred_gaussian_posterior_row));
    evaluation.add(count * p, 2 * sizeof(double));
    if (scratch)
      evaluation.add(*scratch, 1);
    const auto evaluation_peak = evaluation.result();
    if (!retained || !scratch || !evaluation_peak ||
        *evaluation_peak > q->maximum_native_bytes)
      return refuse(s::DensityStatus::invalid_input, n::Status::work_limit);
    result->peak = std::max(result->peak, *evaluation_peak);
    // Allocate the complete pool after every requested case has been charged.
    result->means.resize(count);
    result->rows.resize(count);
    irred::detail::PayloadAccounting pools(sizeof(*result));
    pools.embedded(retained, sizeof(s::GaussianPosterior));
    pools.vector(result->means);
    pools.vector(result->rows);
    pools.add(count * p, 2 * sizeof(double));
    pools.add(*scratch, 1);
    const auto pool_peak = pools.result();
    if (!pool_peak || *pool_peak > q->maximum_native_bytes)
      return refuse(s::DensityStatus::invalid_input, n::Status::work_limit);
    result->peak = std::max(result->peak, *pool_peak);
    policy.maximum_payload_bytes = *scratch;
    result->phase = 3;
    irred::detail::PayloadAccounting base_owned(sizeof(*result));
    base_owned.embedded(retained, sizeof(s::GaussianPosterior));
    base_owned.vector(result->means);
    base_owned.vector(result->rows);
    auto running_owned = base_owned.result();
    if (!running_owned)
      return refuse(s::DensityStatus::invalid_input, n::Status::work_limit);
    for (size_t i = 0; i < count; ++i) {
      result->means[i] = result->posterior->condition(
          {b->conditioning_vectors.data + i * nrow, nrow},
          result->posterior->source().metadata().ordered_ids, policy);
      const auto &m = result->means[i];
      if (!irred::detail::checked_payload_add(
              *running_owned, m.value.capacity(), sizeof(double)) ||
          !irred::detail::checked_payload_add(
              *running_owned, m.absolute_error_estimates.capacity(),
              sizeof(double)))
        return refuse(s::DensityStatus::invalid_input, n::Status::work_limit);
      size_t remaining_peak = *running_owned;
      if (!irred::detail::checked_payload_add(
              remaining_peak, (count - i - 1) * p, 2 * sizeof(double)) ||
          !irred::detail::checked_payload_add(remaining_peak, *scratch, 1) ||
          remaining_peak > q->maximum_native_bytes)
        return refuse(s::DensityStatus::invalid_input, n::Status::work_limit);
      result->peak = std::max(result->peak, remaining_peak);
      const bool finite = m.status == s::DensityStatus::finite;
      result->rows[i] = {
          sizeof(irred_gaussian_posterior_row),
          IRRED_ABI_VERSION,
          static_cast<uint32_t>(m.status),
          static_cast<uint32_t>(m.numerical_status),
          i,
          view(finite ? std::span<const double>(m.value)
                      : std::span<const double>{}),
          view(finite ? std::span<const double>(m.absolute_error_estimates)
                      : std::span<const double>{}),
          finite ? m.backward_residual : 0,
          finite ? m.estimated_forward_sensitivity : 0};
    }
    irred::detail::PayloadAccounting actual(sizeof(*result));
    actual.embedded(retained, sizeof(s::GaussianPosterior));
    actual.vector(result->means);
    actual.vector(result->rows);
    for (const auto &m : result->means) {
      actual.vector(m.value);
      actual.vector(m.absolute_error_estimates);
    }
    const auto owned = actual.result();
    if (!owned || *owned > q->maximum_native_bytes - *scratch) {
      result->means.clear();
      result->rows.clear();
      return refuse(s::DensityStatus::invalid_input, n::Status::work_limit);
    }
    result->retained = *owned;
    result->peak = std::max(result->peak, *owned + *scratch);
    return refuse(s::DensityStatus::finite, n::Status::ok);
  } catch (const std::bad_alloc &) {
    return IRRED_ALLOCATION_FAILURE;
  } catch (...) {
    return IRRED_EXCEPTION;
  }
}
extern "C" uint32_t
irred_gaussian_posterior_result_view(const irred_gaussian_posterior_result *r,
                                     irred_gaussian_posterior_view *v) {
  if (!aligned(r) || !aligned(v))
    return IRRED_INVALID_INPUT;
  const bool finite = r->status == s::DensityStatus::finite;
  *v = {
      sizeof(*v),
      IRRED_ABI_VERSION,
      static_cast<uint32_t>(r->status),
      static_cast<uint32_t>(r->numerical_status),
      r->phase,
      0,
      r->rows.size(),
      r->rows.data(),
      view(finite ? r->posterior->covariance() : std::span<const double>{}),
      finite ? r->posterior->covariance_relative_error_estimate() : 0,
      sizeof(*r),
      r->work,
      r->elements,
      r->peak,
      r->retained,
      bytes(finite
                ? r->posterior->method_id()
                : "proper-Gaussian-parameter-posterior/whitened-precision/v1"),
      finite ? bytes(r->posterior->source().metadata().arithmetic_id)
             : bytes(r->arithmetic)};
  return IRRED_OK;
}
extern "C" uint32_t
irred_gaussian_posterior_result_destroy(irred_gaussian_posterior_result *r) {
  delete r;
  return IRRED_OK;
}
