#include "irred/abi.h"
#include "irred/gaussian_predictive.hpp"
#include "payload_accounting.hpp"
#include <algorithm>
#include <cmath>
#include <memory>
#include <new>
#include <optional>
#include <string_view>
namespace {
namespace s = irred::statistics;
namespace n = irred::numerics;
using PA = irred::detail::PayloadAccounting;
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
std::string_view text(const irred_bytes &b) {
  return {reinterpret_cast<const char *>(b.data),
          static_cast<size_t>(b.length)};
}
bool utf8(std::string_view v) {
  for (size_t i = 0; i < v.size();) {
    auto c = static_cast<unsigned char>(v[i++]);
    if (c < 0x80)
      continue;
    unsigned count = 0;
    uint32_t code = 0, minimum = 0;
    if ((c & 0xe0) == 0xc0) {
      count = 1;
      code = c & 31;
      minimum = 0x80;
    } else if ((c & 0xf0) == 0xe0) {
      count = 2;
      code = c & 15;
      minimum = 0x800;
    } else if ((c & 0xf8) == 0xf0) {
      count = 3;
      code = c & 7;
      minimum = 0x10000;
    } else
      return false;
    if (count > v.size() - i)
      return false;
    while (count--) {
      auto d = static_cast<unsigned char>(v[i++]);
      if ((d & 0xc0) != 0x80)
        return false;
      code = (code << 6) | (d & 63);
    }
    if (code < minimum || code > 0x10ffff || (code >= 0xd800 && code <= 0xdfff))
      return false;
  }
  return true;
}
bool valid_text(const irred_bytes &b, size_t &left) {
  if (!b.data || !b.length || b.length > 256 || b.length > left ||
      !utf8(text(b)))
    return false;
  left -= b.length;
  return true;
}
bool strings(const irred_strings &b, size_t &left) {
  if (b.length > cases_cap || b.byte_length != b.length * sizeof(irred_bytes) ||
      (b.length && !aligned(b.data)))
    return false;
  for (size_t i = 0; i < b.length; ++i)
    if (!valid_text(b.data[i], left))
      return false;
  return true;
}
bool equal(const irred_strings &a, const irred_strings &b) {
  if (a.length != b.length)
    return false;
  for (size_t i = 0; i < a.length; ++i)
    if (text(a.data[i]) != text(b.data[i]))
      return false;
  return true;
}
bool joined(std::string_view target,
            std::initializer_list<std::string_view> parts) {
  size_t pos = 0;
  for (auto p : parts) {
    if (p.size() > target.size() - pos || target.substr(pos, p.size()) != p)
      return false;
    pos += p.size();
  }
  return pos == target.size();
}
bool unique(const irred_strings &a, std::vector<std::string_view> &scratch) {
  for (size_t i = 0; i < a.length; ++i)
    scratch[i] = text(a.data[i]);
  auto end = scratch.begin() + a.length;
  std::sort(scratch.begin(), end);
  return std::adjacent_find(scratch.begin(), end) == end;
}
bool disjoint(const irred_strings &a, const irred_strings &b,
              std::vector<std::string_view> &scratch) {
  for (size_t i = 0; i < a.length; ++i)
    scratch[i] = text(a.data[i]);
  for (size_t i = 0; i < b.length; ++i)
    scratch[a.length + i] = text(b.data[i]);
  auto end = scratch.begin() + a.length + b.length;
  std::sort(scratch.begin(), end);
  return std::adjacent_find(scratch.begin(), end) == end;
}
std::string copy(const irred_bytes &b) { return std::string(text(b)); }
std::vector<std::string> copy(const irred_strings &b) {
  std::vector<std::string> out(b.length);
  for (size_t i = 0; i < b.length; ++i)
    out[i] = copy(b.data[i]);
  return out;
}
std::span<const double> span(const irred_f64_buffer &b) {
  return {b.data, static_cast<size_t>(b.length)};
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
irred_bytes bytes(std::string_view v) {
  return {reinterpret_cast<const uint8_t *>(v.data()), v.size()};
}
void charge(PA &a, const s::Metadata &m) {
  a.add(1, sizeof(m));
  a.strings(m.ordered_ids);
  for (auto p : {&m.arithmetic_id, &m.measure, &m.table_identity,
                 &m.uncertainty_identity, &m.ordering_provenance,
                 &m.calibration_provenance, &m.dependence_provenance,
                 &m.source_semantics, &m.input_matrix_convention, &m.treatment})
    a.string(*p);
}
void charge(PA &a, const s::ParameterPrior &p) {
  a.add(1, sizeof(p));
  a.strings(p.ordered_parameter_ids);
  a.strings(p.parameter_units);
  a.strings(p.shared_nuisance_ids);
  a.vector(p.mean);
  a.vector(p.covariance);
  for (auto v : {&p.prior_identity, &p.design_identity, &p.residual_unit,
                 &p.parameter_measure, &p.dependence_identity})
    a.string(*v);
}
void charge(PA &a, const s::PredictiveMetadata &m) {
  a.add(1, sizeof(m));
  for (auto v : {&m.ordered_parameter_ids, &m.parameter_units,
                 &m.response_units, &m.training_event_ids, &m.future_event_ids})
    a.strings(*v);
  for (auto v : {&m.future_unit, &m.future_covariance_unit, &m.future_measure,
                 &m.response_identity, &m.conditioning_identity,
                 &m.conditional_noise_identity, &m.dependence_identity})
    a.string(*v);
}
std::optional<size_t> scaled(PA a) {
  auto z = a.result();
  if (!z || *z > SIZE_MAX / 32)
    return {};
  return *z * 32;
}
} // namespace
struct irred_gaussian_predictive_result {
  std::optional<s::GaussianPredictiveConditioning> predictive;
  std::optional<s::PredictiveBatch> batch;
  std::vector<irred_gaussian_predictive_row> rows;
  irred_gaussian_predictive_view view{};
  irred_gaussian_predictive_result() noexcept {
    view.struct_size = sizeof(view);
    view.abi_version = IRRED_ABI_VERSION;
    view.status = static_cast<uint32_t>(s::DensityStatus::invalid_input);
    view.numerical_status = static_cast<uint32_t>(n::Status::invalid_input);
    view.minimum_result_bytes = view.peak_payload_bytes =
        view.retained_payload_bytes = sizeof(*this);
  }
};
extern "C" uint32_t
irred_gaussian_predictive_evaluate(const irred_gaussian_predictive_batch *b,
                                   const irred_gaussian_posterior_policy *q,
                                   irred_gaussian_predictive_result **output) {
  if (!aligned(output))
    return IRRED_INVALID_INPUT;
  *output = nullptr;
  if (!aligned(b) || !aligned(q))
    return IRRED_INVALID_INPUT;
  if (b->abi_version != IRRED_ABI_VERSION ||
      q->abi_version != IRRED_ABI_VERSION ||
      b->training.abi_version != IRRED_ABI_VERSION)
    return IRRED_ABI_MISMATCH;
  if (!header(b) || !header(q) || !header(&b->training) || b->reserved ||
      q->reserved || q->arithmetic > 1 || b->requested_outputs < 1 ||
      b->requested_outputs > 3 || b->future_noise_independence_declared > 1 ||
      b->noise_conditional_on_parameters_declared > 1 ||
      b->training.noise_independence_declared > 1 ||
      q->maximum_elements > elements_cap || q->maximum_cases > cases_cap ||
      q->maximum_string_bytes > string_cap ||
      q->maximum_native_bytes > bytes_cap || q->maximum_work_units > work_cap ||
      !std::isfinite(q->maximum_forward_sensitivity) ||
      q->maximum_forward_sensitivity <= 0 ||
      q->maximum_forward_sensitivity > 1e-10 ||
      b->training.case_count > cases_cap)
    return IRRED_INVALID_INPUT;
  const auto &t = b->training;
  for (auto v :
       {&t.noise_covariance, &t.design, &t.prior_mean, &t.prior_covariance,
        &t.conditioning_vectors, &b->future_noise_covariance,
        &b->future_response, &b->future_vectors})
    if (!buffer(*v))
      return IRRED_INVALID_INPUT;
  size_t left = string_cap, slots = 0;
  for (auto v : {&t.noise_row_ids,
                 &t.noise_event_ids,
                 &t.design_row_ids,
                 &t.design_parameter_ids,
                 &t.design_parameter_units,
                 &t.column_units,
                 &t.prior_parameter_ids,
                 &t.prior_parameter_units,
                 &t.shared_nuisance_ids,
                 &t.conditioning_row_ids,
                 &t.conditioning_event_ids,
                 &t.case_ids,
                 &b->future_noise_row_ids,
                 &b->future_noise_event_ids,
                 &b->future_response_row_ids,
                 &b->future_parameter_ids,
                 &b->future_parameter_units,
                 &b->future_column_units,
                 &b->future_vector_row_ids,
                 &b->future_vector_event_ids})
    if (!strings(*v, left) ||
        !irred::detail::checked_payload_add(slots, v->length, 1))
      return IRRED_INVALID_INPUT;
  for (auto v : {&t.residual_unit, &t.noise_identity, &t.calibration_identity,
                 &t.noise_dependence_identity, &t.ordering_provenance,
                 &t.design_identity, &t.prior_identity, &t.parameter_measure,
                 &t.prior_dependence_identity, &b->future_unit,
                 &b->future_noise_identity, &b->future_calibration_identity,
                 &b->future_noise_dependence_identity,
                 &b->future_ordering_provenance, &b->future_response_identity,
                 &b->future_covariance_unit, &b->future_measure,
                 &b->conditioning_identity, &b->prediction_dependence_identity})
    if (!valid_text(*v, left))
      return IRRED_INVALID_INPUT;
  if (q->maximum_native_bytes < sizeof(irred_gaussian_predictive_result))
    return IRRED_GAUSSIAN_PREDICTIVE_QUOTA_REFUSED;
  std::unique_ptr<irred_gaussian_predictive_result> result;
  auto refuse = [&](s::DensityStatus status, n::Status numerical) {
    result->predictive.reset();
    result->batch.reset();
    std::vector<irred_gaussian_predictive_row>().swap(result->rows);
    auto &v = result->view;
    v.status = static_cast<uint32_t>(status);
    v.numerical_status = static_cast<uint32_t>(numerical);
    v.rows = nullptr;
    v.case_count = 0;
    v.retained_payload_bytes = sizeof(*result);
    *output = result.release();
    return IRRED_OK;
  };
  try {
    result = std::make_unique<irred_gaussian_predictive_result>();
    auto &v = result->view;
    v.requested_outputs = b->requested_outputs;
    const size_t nr = t.noise_row_ids.length, p = t.prior_parameter_ids.length,
                 k = b->future_noise_row_ids.length, count = t.case_count;
    const bool means = b->requested_outputs & 1,
               density = b->requested_outputs & 2;
    if (!nr || p < 2 || !k || !count)
      return refuse(s::DensityStatus::invalid_input, n::Status::invalid_input);
    if (nr > elements_cap / nr || p > elements_cap / p ||
        k > elements_cap / k || nr > elements_cap / p || k > elements_cap / p ||
        count > elements_cap / nr || count > elements_cap / k)
      return refuse(s::DensityStatus::invalid_input, n::Status::work_limit);
    if (t.noise_covariance.length != nr * nr || t.design.length != nr * p ||
        t.prior_mean.length != p || t.prior_covariance.length != p * p ||
        t.conditioning_vectors.length != count * nr ||
        t.case_ids.length != count ||
        b->future_noise_covariance.length != k * k ||
        b->future_response.length != k * p ||
        (density && b->future_vectors.length != count * k) ||
        (!density &&
         (b->future_vectors.length || b->future_vector_row_ids.length ||
          b->future_vector_event_ids.length)))
      return refuse(s::DensityStatus::invalid_input, n::Status::invalid_input);
    PA wc(0), wp(0), wr(0), wa(0), wb(0), el(0), sum(0);
    wc.add(nr * nr, nr);
    wp.add(nr * nr, p);
    wp.add(nr * p, p);
    wp.add(p * p, 8 * p);
    wr.add(k * k, k);
    wa.add(nr * nr, 1);
    wa.add(nr * p, 1);
    wa.add(p * p, 1);
    wa.add(k * p, p);
    wa.add(k * k, p);
    wa.add(k * k, 8 * k);
    wb.add(count, nr * nr + nr * p + p * p + k * p + (density ? k * k : 0));
    const auto c = scaled(wc), post = scaled(wp), r = scaled(wr),
               a = scaled(wa), batch = scaled(wb);
    if (!c || !post || !r || !a || !batch)
      return refuse(s::DensityStatus::invalid_input, n::Status::work_limit);
    for (auto w : {*c, *post, *r, *a, *batch})
      sum.add(w, 1);
    auto total = sum.result();
    if (!total)
      return refuse(s::DensityStatus::invalid_input, n::Status::work_limit);
    v.declared_training_noise_work = *c;
    v.declared_posterior_work = *post;
    v.declared_future_noise_work = *r;
    v.declared_predictive_work = *a;
    v.declared_batch_work = *batch;
    v.declared_total_work = *total;
    for (auto z : {nr * nr, nr * p, p, p * p, k * k, k * p, count * nr})
      el.add(z, 1);
    if (density) {
      el.add(count * k, 1);
      el.add(count, 6);
    }
    if (means)
      el.add(count * k, 2);
    auto elems = el.result();
    if (!elems)
      return refuse(s::DensityStatus::invalid_input, n::Status::work_limit);
    v.pooled_numeric_elements = *elems;
    if (count > q->maximum_cases || *elems > q->maximum_elements ||
        *total > q->maximum_work_units ||
        string_cap - left > q->maximum_string_bytes)
      return refuse(s::DensityStatus::invalid_input, n::Status::work_limit);
    if (t.source_semantics != 1 || t.noise_independence_declared != 1 ||
        b->future_noise_independence_declared != 1 ||
        b->noise_conditional_on_parameters_declared != 1 ||
        t.noise_event_ids.length != nr ||
        b->future_noise_event_ids.length != k ||
        !equal(t.noise_row_ids, t.design_row_ids) ||
        !equal(t.noise_row_ids, t.conditioning_row_ids) ||
        !equal(t.noise_event_ids, t.conditioning_event_ids) ||
        !equal(t.prior_parameter_ids, t.design_parameter_ids) ||
        !equal(t.prior_parameter_ids, b->future_parameter_ids) ||
        !equal(t.prior_parameter_units, t.design_parameter_units) ||
        !equal(t.prior_parameter_units, b->future_parameter_units) ||
        t.prior_parameter_units.length != p || t.column_units.length != p ||
        b->future_column_units.length != p ||
        !equal(b->future_noise_row_ids, b->future_response_row_ids) ||
        (density &&
         (!equal(b->future_noise_row_ids, b->future_vector_row_ids) ||
          !equal(b->future_noise_event_ids, b->future_vector_event_ids))) ||
        text(t.residual_unit) != text(b->future_unit) ||
        !joined(text(b->future_covariance_unit),
                {text(t.residual_unit), "^2"}) ||
        !joined(text(b->future_measure),
                {"product d(", text(t.residual_unit), ")"}) ||
        t.shared_nuisance_ids.length > p)
      return refuse(s::DensityStatus::incompatible_metadata,
                    n::Status::invalid_input);
    for (size_t j = 0; j < p; ++j)
      if (!joined(text(t.column_units.data[j]),
                  {text(t.residual_unit), "/",
                   text(t.prior_parameter_units.data[j])}) ||
          !joined(text(b->future_column_units.data[j]),
                  {text(b->future_unit), "/",
                   text(t.prior_parameter_units.data[j])}))
        return refuse(s::DensityStatus::incompatible_metadata,
                      n::Status::invalid_input);
    PA envelope(sizeof(*result) + 2 * sizeof(s::Metadata) +
                sizeof(s::ParameterPrior) + sizeof(s::PredictiveMetadata));
    const size_t vs = 2 * nr + 2 * k + 5 * p + t.shared_nuisance_ids.length;
    envelope.add(p + p * p, sizeof(double));
    envelope.add(vs, sizeof(std::string));
    envelope.add(vs + 32, 257);
    envelope.add(slots, sizeof(std::string_view));
    auto conv = envelope.result();
    if (!conv || *conv > q->maximum_native_bytes)
      return refuse(s::DensityStatus::invalid_input, n::Status::work_limit);
    v.peak_payload_bytes = std::max<uint64_t>(v.peak_payload_bytes, *conv);
    { // one bounded validation scratch, gone before conversion/native setup
      std::vector<std::string_view> scratch(slots);
      if (scratch.capacity() > slots)
        return refuse(s::DensityStatus::invalid_input, n::Status::work_limit);
      for (auto ids :
           {&t.noise_row_ids, &t.noise_event_ids, &t.prior_parameter_ids,
            &t.shared_nuisance_ids, &t.case_ids, &b->future_noise_row_ids,
            &b->future_noise_event_ids})
        if (!unique(*ids, scratch))
          return refuse(s::DensityStatus::incompatible_metadata,
                        n::Status::invalid_input);
      if (!disjoint(t.noise_row_ids, b->future_noise_row_ids, scratch) ||
          !disjoint(t.noise_event_ids, b->future_noise_event_ids, scratch))
        return refuse(s::DensityStatus::incompatible_metadata,
                      n::Status::invalid_input);
      for (size_t j = 0; j < t.shared_nuisance_ids.length; ++j) {
        bool found = false;
        for (size_t z = 0; z < p; ++z)
          found |= text(t.shared_nuisance_ids.data[j]) ==
                   text(t.prior_parameter_ids.data[z]);
        if (!found)
          return refuse(s::DensityStatus::incompatible_metadata,
                        n::Status::invalid_input);
      }
    }
    s::Metadata cm, rm;
    auto metadata = [&](s::Metadata &m, const irred_strings &ids,
                        const irred_bytes &identity, const irred_bytes &cal,
                        const irred_bytes &dep, const irred_bytes &order) {
      m.ordered_ids = copy(ids);
      m.measure = copy(b->future_measure);
      m.table_identity = "synthetic supplied conditional law";
      m.uncertainty_identity = copy(identity);
      m.calibration_provenance = copy(cal);
      m.dependence_provenance = copy(dep);
      m.ordering_provenance = copy(order);
      m.source_semantics = "synthetic controls";
      m.input_matrix_convention = "covariance";
    };
    metadata(cm, t.noise_row_ids, t.noise_identity, t.calibration_identity,
             t.noise_dependence_identity, t.ordering_provenance);
    metadata(rm, b->future_noise_row_ids, b->future_noise_identity,
             b->future_calibration_identity,
             b->future_noise_dependence_identity,
             b->future_ordering_provenance);
    s::ParameterPrior prior;
    prior.ordered_parameter_ids = copy(t.prior_parameter_ids);
    prior.parameter_units = copy(t.prior_parameter_units);
    prior.shared_nuisance_ids = copy(t.shared_nuisance_ids);
    prior.mean.assign(t.prior_mean.data, t.prior_mean.data + p);
    prior.covariance.assign(t.prior_covariance.data,
                            t.prior_covariance.data + p * p);
    prior.prior_identity = copy(t.prior_identity);
    prior.design_identity = copy(t.design_identity);
    prior.residual_unit = copy(t.residual_unit);
    prior.parameter_measure = copy(t.parameter_measure);
    prior.dependence_identity = copy(t.prior_dependence_identity);
    prior.noise_independence_declared = true;
    s::PredictiveMetadata md;
    md.ordered_parameter_ids = copy(t.prior_parameter_ids);
    md.parameter_units = copy(t.prior_parameter_units);
    md.response_units = copy(b->future_column_units);
    md.training_event_ids = copy(t.noise_event_ids);
    md.future_event_ids = copy(b->future_noise_event_ids);
    md.future_unit = copy(b->future_unit);
    md.future_covariance_unit = copy(b->future_covariance_unit);
    md.future_measure = copy(b->future_measure);
    md.response_identity = copy(b->future_response_identity);
    md.conditioning_identity = copy(b->conditioning_identity);
    md.conditional_noise_identity = copy(b->future_noise_identity);
    md.dependence_identity = copy(b->prediction_dependence_identity);
    md.future_noise_independence_declared =
        md.noise_conditional_on_parameters_declared = true;
    PA converted(sizeof(*result));
    charge(converted, cm);
    charge(converted, rm);
    charge(converted, prior);
    charge(converted, md);
    auto converted_bytes = converted.result();
    if (!converted_bytes || *converted_bytes > *conv)
      return refuse(s::DensityStatus::invalid_input, n::Status::work_limit);
    auto phase = [&](std::optional<size_t> native, PA extra) {
      if (native)
        extra.add(*native, 1);
      auto peak = extra.result();
      if (!native || !peak || *peak > q->maximum_native_bytes)
        return false;
      v.peak_payload_bytes = std::max<uint64_t>(v.peak_payload_bytes, *peak);
      return true;
    };
    const auto arithmetic = q->arithmetic ? n::Arithmetic::longdouble_cpu_v1
                                          : n::Arithmetic::binary64_legacy_v1;
    PA cx(sizeof(*result));
    charge(cx, rm);
    charge(cx, prior);
    charge(cx, md);
    if (!phase(s::gaussian_preparation_payload_bound(
                   nr, s::MatrixKind::covariance, arithmetic, cm),
               cx))
      return refuse(s::DensityStatus::invalid_input, n::Status::work_limit);
    v.phase = 1;
    v.training_noise_attempted = 1;
    v.arithmetic = bytes(q->arithmetic ? "F02/longdouble-cpu/v1"
                                       : "F02/binary64-legacy/v1");
    auto source = s::prepare_gaussian(
        span(t.noise_covariance), s::MatrixKind::covariance, std::move(cm),
        q->maximum_elements, q->maximum_forward_sensitivity, arithmetic);
    if (source.status() != s::DensityStatus::finite)
      return refuse(source.status(), source.numerical_status());
    v.training_noise_prepared = 1;
    PA px(sizeof(*result));
    charge(px, rm);
    charge(px, md);
    if (!phase(s::GaussianPosterior::preparation_payload_bound(source, prior),
               px))
      return refuse(s::DensityStatus::invalid_input, n::Status::work_limit);
    s::PosteriorPolicy pp{static_cast<size_t>(q->maximum_elements),
                          static_cast<size_t>(q->maximum_native_bytes),
                          static_cast<size_t>(q->maximum_work_units - *c),
                          q->maximum_forward_sensitivity};
    v.phase = 2;
    v.posterior_attempted = 1;
    auto posterior = s::GaussianPosterior::prepare(
        std::move(source), span(t.design), source.metadata().ordered_ids,
        std::move(prior), pp);
    if (posterior.status() != s::DensityStatus::finite)
      return refuse(posterior.status(), posterior.numerical_status());
    v.posterior_prepared = 1;
    PA rx(sizeof(*result));
    charge(rx, md);
    rx.add(posterior.retained_payload_bound().value_or(SIZE_MAX), 1);
    if (!phase(s::gaussian_preparation_payload_bound(
                   k, s::MatrixKind::covariance, arithmetic, rm),
               rx))
      return refuse(s::DensityStatus::invalid_input, n::Status::work_limit);
    v.phase = 3;
    v.future_noise_attempted = 1;
    auto future = s::prepare_gaussian(
        span(b->future_noise_covariance), s::MatrixKind::covariance,
        std::move(rm), q->maximum_elements, q->maximum_forward_sensitivity,
        arithmetic);
    if (future.status() != s::DensityStatus::finite)
      return refuse(future.status(), future.numerical_status());
    v.future_noise_prepared = 1;
    PA ax(sizeof(
        *result)); // native bound includes post/R and md once (plus candidates)
    if (!phase(s::GaussianPredictiveConditioning::preparation_payload_bound(
                   posterior, future, md),
               ax))
      return refuse(s::DensityStatus::invalid_input, n::Status::work_limit);
    s::PredictivePolicy ap{
        static_cast<size_t>(q->maximum_elements),
        static_cast<size_t>(q->maximum_native_bytes),
        static_cast<size_t>(q->maximum_work_units - *c - *post - *r),
        q->maximum_forward_sensitivity};
    v.phase = 4;
    v.predictive_attempted = 1;
    result->predictive.emplace(s::GaussianPredictiveConditioning::prepare(
        std::move(posterior), future, span(b->future_response), md, ap));
    if (result->predictive->status() != s::DensityStatus::finite)
      return refuse(result->predictive->status(),
                    result->predictive->numerical_status());
    v.predictive_prepared = 1;
    const s::PredictiveOutputs requested{means, density};
    PA outside(sizeof(*result));
    outside.add(count, sizeof(irred_gaussian_predictive_row));
    charge(outside, md);
    outside.add(future.retained_payload_bound().value_or(SIZE_MAX), 1);
    auto outside_bytes = outside.result(),
         bp = result->predictive->batch_payload_bound(count, requested);
    // Native batch bound includes both embedded owners; deduct their sizeof,
    // already reserved inside the allocation-free result's optional storage.
    PA eval(sizeof(*result));
    eval.add(count, sizeof(irred_gaussian_predictive_row));
    charge(eval, md);
    eval.add(future.retained_payload_bound().value_or(SIZE_MAX), 1);
    eval.embedded(bp, sizeof(s::GaussianPredictiveConditioning) +
                          sizeof(s::PredictiveBatch));
    auto ep = eval.result();
    if (!outside_bytes || !bp || !ep || *ep > q->maximum_native_bytes)
      return refuse(s::DensityStatus::invalid_input, n::Status::work_limit);
    v.peak_payload_bytes = std::max<uint64_t>(v.peak_payload_bytes, *ep);
    result->rows.resize(count);
    if (result->rows.capacity() != count)
      return refuse(s::DensityStatus::invalid_input, n::Status::work_limit);
    ap.maximum_work_units = q->maximum_work_units - *c - *post - *r - *a;
    // Native evaluation includes owner object sizes, so restore their embedded
    // sizes to the available quota after subtracting all outside payload.
    ap.maximum_payload_bytes = q->maximum_native_bytes - *outside_bytes +
                               sizeof(s::GaussianPredictiveConditioning) +
                               sizeof(s::PredictiveBatch);
    v.phase = 5;
    v.batch_called = 1;
    result->batch.emplace(result->predictive->evaluate(
        span(t.conditioning_vectors),
        result->predictive->posterior().source().metadata().ordered_ids,
        span(b->future_vectors), future.metadata().ordered_ids, count,
        requested, ap));
    auto &bat = *result->batch;
    if (bat.status != s::DensityStatus::finite)
      return refuse(bat.status, bat.numerical_status);
    v.batch_admission_completed = 1;
    if (bat.rows.size() != count || bat.work_units != *batch ||
        bat.means.size() != (means ? count * k : 0) ||
        bat.mean_absolute_error_estimates.size() != (means ? count * k : 0) ||
        bat.densities.size() != (density ? count : 0))
      return refuse(s::DensityStatus::numerical_failure,
                    n::Status::invalid_input);
    for (size_t i = 0; i < count; ++i) {
      auto &row = result->rows[i];
      const auto &br = bat.rows[i];
      const bool ok = br.status == s::DensityStatus::finite;
      if (ok && br.numerical_status != n::Status::ok)
        return refuse(s::DensityStatus::numerical_failure,
                      n::Status::invalid_input);
      if (ok && means)
        for (size_t j = 0; j < k; ++j) {
          const auto value = bat.means[i * k + j],
                     error = bat.mean_absolute_error_estimates[i * k + j];
          if (!std::isfinite(value) || !std::isfinite(error) || error < 0 ||
              error / (1 + std::abs(value)) > q->maximum_forward_sensitivity)
            return refuse(s::DensityStatus::numerical_failure,
                          n::Status::conditioning_budget_exceeded);
        }
      row = {};
      row.struct_size = sizeof(row);
      row.abi_version = IRRED_ABI_VERSION;
      row.status = static_cast<uint32_t>(br.status);
      row.numerical_status = static_cast<uint32_t>(br.numerical_status);
      row.case_index = i;
      row.mean = view(ok && means
                          ? std::span<const double>(bat.means).subspan(i * k, k)
                          : std::span<const double>{});
      row.absolute_error_estimates =
          view(ok && means
                   ? std::span<const double>(bat.mean_absolute_error_estimates)
                         .subspan(i * k, k)
                   : std::span<const double>{});
      if (ok && density) {
        const auto &d = bat.densities[i];
        if (d.density.status != s::DensityStatus::finite ||
            d.density.numerical_status != n::Status::ok ||
            !std::isfinite(d.density.log_value) ||
            !std::isfinite(d.quadratic) || !std::isfinite(d.log_determinant) ||
            !std::isfinite(d.normalization) ||
            !std::isfinite(d.backward_residual) ||
            !std::isfinite(d.estimated_forward_sensitivity) ||
            d.quadratic < 0 || d.normalization < 0 || d.backward_residual < 0 ||
            d.estimated_forward_sensitivity < 0 ||
            d.estimated_forward_sensitivity > q->maximum_forward_sensitivity)
          return refuse(s::DensityStatus::numerical_failure,
                        n::Status::invalid_input);
        row.joint_density_available = 1;
        row.log_density = d.density.log_value;
        row.quadratic = d.quadratic;
        row.log_determinant = d.log_determinant;
        row.normalization = d.normalization;
        row.backward_residual = d.backward_residual;
        row.estimated_forward_sensitivity = d.estimated_forward_sensitivity;
      }
    }
    PA owned(sizeof(*result));
    owned.embedded(result->predictive->retained_payload_bound(),
                   sizeof(s::GaussianPredictiveConditioning));
    owned.vector(result->rows);
    owned.vector(bat.rows);
    owned.vector(bat.means);
    owned.vector(bat.mean_absolute_error_estimates);
    owned.vector(bat.densities);
    auto retained = owned.result();
    if (!retained || *retained > q->maximum_native_bytes)
      return refuse(s::DensityStatus::invalid_input, n::Status::work_limit);
    v.retained_payload_bytes = *retained;
    v.peak_payload_bytes = std::max<uint64_t>(v.peak_payload_bytes, *retained);
    v.status = static_cast<uint32_t>(s::DensityStatus::finite);
    v.numerical_status = static_cast<uint32_t>(n::Status::ok);
    v.rows = result->rows.data();
    v.case_count = count;
    *output = result.release();
    return IRRED_OK;
  } catch (const std::bad_alloc &) {
    if (result)
      return refuse(s::DensityStatus::invalid_input, n::Status::work_limit);
    return IRRED_ALLOCATION_FAILURE;
  } catch (const std::length_error &) {
    if (result)
      return refuse(s::DensityStatus::invalid_input, n::Status::work_limit);
    return IRRED_ALLOCATION_FAILURE;
  } catch (...) {
    return IRRED_EXCEPTION;
  }
}
extern "C" uint32_t
irred_gaussian_predictive_result_view(const irred_gaussian_predictive_result *r,
                                      irred_gaussian_predictive_view *v) {
  if (!aligned(r) || !aligned(v))
    return IRRED_INVALID_INPUT;
  *v = r->view;
  v->method = bytes(
      v->predictive_attempted
          ? "proper-Gaussian-repeated-joint-predictive/retained-covariance/v1"
      : v->posterior_attempted && !v->future_noise_attempted
          ? "proper-Gaussian-parameter-posterior/whitened-precision/v1"
      : v->training_noise_attempted
          ? "supplied-Gaussian-covariance-preparation/Cholesky/v1"
          : "unexecuted");
  return IRRED_OK;
}
extern "C" uint32_t
irred_gaussian_predictive_result_destroy(irred_gaussian_predictive_result *r) {
  delete r;
  return IRRED_OK;
}
