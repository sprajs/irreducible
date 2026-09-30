#include "irred/abi.h"
#include "irred/statistics.hpp"
#include "payload_accounting.hpp"
#include "observation_internal.hpp"
#include <cmath>
#include <limits>
#include <memory>
#include <new>
#include <string>
#include <utility>
#include <vector>
namespace {
using irred::statistics::DensityStatus;
static_assert(static_cast<uint32_t>(DensityStatus::finite) ==
              IRRED_GAUSSIAN_STATUS_FINITE);
static_assert(static_cast<uint32_t>(DensityStatus::outside_support) ==
              IRRED_GAUSSIAN_STATUS_OUTSIDE_SUPPORT);
static_assert(static_cast<uint32_t>(DensityStatus::invalid_input) ==
              IRRED_GAUSSIAN_STATUS_INVALID_INPUT);
static_assert(static_cast<uint32_t>(DensityStatus::unsupported_domain) ==
              IRRED_GAUSSIAN_STATUS_UNSUPPORTED_DOMAIN);
static_assert(static_cast<uint32_t>(DensityStatus::numerical_failure) ==
              IRRED_GAUSSIAN_STATUS_NUMERICAL_FAILURE);
static_assert(static_cast<uint32_t>(DensityStatus::incompatible_metadata) ==
              IRRED_GAUSSIAN_STATUS_INCOMPATIBLE_METADATA);
static_assert(
    static_cast<uint32_t>(
        irred::statistics::MatrixValidationScope::full_declared_matrix) ==
    IRRED_GAUSSIAN_MATRIX_VALIDATION_SCOPE_FULL_DECLARED_MATRIX);
static_assert(
    static_cast<uint32_t>(
        irred::statistics::MatrixValidationScope::selected_covariance_only) ==
    IRRED_GAUSSIAN_MATRIX_VALIDATION_SCOPE_SELECTED_COVARIANCE_ONLY);
static_assert(
    static_cast<uint32_t>(irred::statistics::MatrixValidationScope::
                              full_precision_then_marginal) ==
    IRRED_GAUSSIAN_MATRIX_VALIDATION_SCOPE_FULL_PRECISION_THEN_MARGINAL);
template <class T>
bool bounded(const T *data, uint64_t n, uint64_t bytes, uint64_t cap) {
  return n <= cap && n <= std::numeric_limits<size_t>::max() / sizeof(T) &&
         bytes == n * sizeof(T) &&
         (!n || (data && reinterpret_cast<uintptr_t>(data) % alignof(T) == 0));
}
bool valid(const irred_gaussian_policy *p) {
  return p && p->struct_size == sizeof(*p) &&
         p->abi_version == IRRED_ABI_VERSION && !p->reserved && !p->reserved2 &&
         p->maximum_matrix_elements <= std::numeric_limits<size_t>::max() &&
         p->maximum_batch_elements <= IRRED_MAX_BATCH_ELEMENTS &&
         p->maximum_string_bytes <= std::numeric_limits<size_t>::max() &&
         std::isfinite(p->maximum_forward_sensitivity) &&
         p->maximum_forward_sensitivity > 0;
}
bool doubles(const irred_f64_buffer &b, uint64_t cap) {
  return b.struct_size == sizeof(b) && b.abi_version == IRRED_ABI_VERSION &&
         b.element_type == 2 && !b.reserved &&
         bounded(b.data, b.length, b.byte_length, cap);
}
bool bytes(const irred_bytes &b, uint64_t &left) {
  if (!bounded(b.data, b.length, b.length, left))
    return false;
  left -= b.length;
  return true;
}
bool strings(const irred_strings &b, uint64_t expected, uint64_t &left) {
  if (b.length != expected ||
      !bounded(b.data, b.length, b.byte_length, expected))
    return false;
  for (uint64_t i = 0; i < b.length; ++i)
    if (!bytes(b.data[i], left))
      return false;
  return true;
}
std::string copy(const irred_bytes &b) {
  return b.length
             ? std::string(reinterpret_cast<const char *>(b.data), b.length)
             : std::string{};
}
std::vector<std::string> copy(const irred_strings &b) {
  std::vector<std::string> v;
  v.reserve(b.length);
  for (uint64_t i = 0; i < b.length; ++i)
    v.push_back(copy(b.data[i]));
  return v;
}
irred_bytes view(const std::string &s) {
  return {reinterpret_cast<const uint8_t *>(s.data()), s.size()};
}
irred_f64_buffer view(std::span<const double> s) {
  return {sizeof(irred_f64_buffer),
          IRRED_ABI_VERSION,
          2,
          0,
          s.data(),
          s.size(),
          s.size_bytes()};
}
irred_strings view(const std::vector<irred_bytes> &v) {
  return {v.data(), v.size(), v.size() * sizeof(irred_bytes)};
}
} // namespace
struct irred_gaussian {
  irred::statistics::Gaussian native;
  std::vector<irred_bytes> ids;
  std::vector<std::vector<irred_bytes>> prior_ids, kept_ids, complement_ids;
  explicit irred_gaussian(irred::statistics::Gaussian g)
      : native(std::move(g)) {
    ids.reserve(native.metadata().ordered_ids.size());
    kept_ids.reserve(native.selection_history().size());
    complement_ids.reserve(native.selection_history().size());
    prior_ids.reserve(native.priors().size());
    for (const auto &id : native.metadata().ordered_ids)
      ids.push_back(view(id));
    for (const auto &selection : native.selection_history()) {
      kept_ids.emplace_back();
      complement_ids.emplace_back();
      kept_ids.back().reserve(selection.kept_row_ids.size());
      complement_ids.back().reserve(selection.complement_row_ids.size());
      for (const auto &id : selection.kept_row_ids)
        kept_ids.back().push_back(view(id));
      for (const auto &id : selection.complement_row_ids)
        complement_ids.back().push_back(view(id));
    }
    for (const auto &prior : native.priors()) {
      prior_ids.emplace_back();
      prior_ids.back().reserve(prior.applied_row_ids.size());
      for (const auto &id : prior.applied_row_ids)
        prior_ids.back().push_back(view(id));
    }
  }
};
struct irred_gaussian_result {
  std::vector<irred_gaussian_row> rows;
};
extern "C" uint32_t
irred_gaussian_prepare(const irred_prepared *observations, uint32_t selection,
                       const irred_gaussian_policy *p, irred_gaussian **out,
                       uint32_t *semantic, uint32_t *numerical_status) {
  if (out)
    *out = nullptr;
  if (semantic)
    *semantic = IRRED_GAUSSIAN_STATUS_INVALID_INPUT;
  if (numerical_status)
    *numerical_status =
        static_cast<uint32_t>(irred::numerics::Status::invalid_input);
  if (!out || !semantic || !numerical_status || !observations || !p)
    return IRRED_INVALID_INPUT;
  if (p->abi_version != IRRED_ABI_VERSION)
    return IRRED_ABI_MISMATCH;
  if (!valid(p))
    return IRRED_INVALID_INPUT;
  const auto &source = native_observations(observations)->source();
  uint64_t string_left = p->maximum_string_bytes;
  for (const auto *value :
       {&source.table_sha256, &source.uncertainty_sha256,
        &source.calibration_provenance, &source.dependence_provenance,
        &source.ordering_provenance})
    if (!bytes(view(*value), string_left))
      return IRRED_INVALID_INPUT;
  for (const auto &id : source.measurement_ids)
    if (!bytes(view(id), string_left))
      return IRRED_INVALID_INPUT;
  const auto peak = irred::statistics::selected_gaussian_preparation_payload_bound(
      *native_observations(observations), source.values.size(),
      irred::numerics::Arithmetic::binary64_legacy_v1);
  irred::detail::PayloadAccounting combined(peak.value_or(SIZE_MAX));
  combined.add(1, sizeof(irred_gaussian));
  combined.add(source.values.size(), 2 * sizeof(irred_bytes));
  combined.add(2, sizeof(std::vector<irred_bytes>));
  const auto envelope = combined.result();
  if (!peak || !envelope || *envelope > p->maximum_native_bytes)
    return IRRED_INVALID_INPUT;
  try {
    auto result = irred::statistics::prepare_observations(
        *native_observations(observations),
        static_cast<irred::observations::Selection>(selection),
        p->maximum_matrix_elements, p->maximum_forward_sensitivity);
    *semantic = static_cast<uint32_t>(result.status());
    *numerical_status = static_cast<uint32_t>(result.numerical_status());
    if (result.status() == DensityStatus::finite)
      *out = new irred_gaussian(std::move(result));
    return IRRED_OK;
  } catch (const std::bad_alloc &) {
    return IRRED_ALLOCATION_FAILURE;
  } catch (...) {
    return IRRED_EXCEPTION;
  }
}
extern "C" uint32_t irred_gaussian_proper_offset(
    const irred_gaussian *g, const irred_gaussian_prior *prior,
    const irred_gaussian_policy *p, irred_gaussian **out, uint32_t *semantic,
    uint32_t *numerical_status) {
  if (out)
    *out = nullptr;
  if (semantic)
    *semantic = IRRED_GAUSSIAN_STATUS_INVALID_INPUT;
  if (numerical_status)
    *numerical_status =
        static_cast<uint32_t>(irred::numerics::Status::invalid_input);
  if (!g || !prior || !p || !out || !semantic || !numerical_status)
    return IRRED_INVALID_INPUT;
  if (prior->abi_version != IRRED_ABI_VERSION ||
      p->abi_version != IRRED_ABI_VERSION)
    return IRRED_ABI_MISMATCH;
  const auto n = g->ids.size();
  if (n && n > p->maximum_matrix_elements / n)
    return IRRED_INVALID_INPUT;
  uint64_t left = p->maximum_string_bytes;
  if (!valid(p) || prior->struct_size != sizeof(*prior) || prior->reserved ||
      prior->independence_declared > 1 || !doubles(prior->response, n) ||
      prior->response.length != n || !bytes(prior->latent_identity, left) ||
      !strings(prior->ordered_ids, n, left))
    return IRRED_INVALID_INPUT;
  const auto peak = g->native.proper_offset_payload_bound(
      {reinterpret_cast<const char *>(prior->latent_identity.data),
       static_cast<size_t>(prior->latent_identity.length)});
  irred::detail::PayloadAccounting combined(peak.value_or(SIZE_MAX));
  combined.add(2, sizeof(irred_gaussian));
  // Existing views remain live alongside candidate views; reserve exact
  // candidate lengths, while charging existing actual capacities.
  auto existing = [&](const auto &groups) {
    combined.add(groups.capacity(), sizeof(std::vector<irred_bytes>));
    for (const auto &ids : groups) combined.add(ids.capacity(), sizeof(irred_bytes));
  };
  combined.add(g->ids.capacity(), sizeof(irred_bytes));
  existing(g->kept_ids); existing(g->complement_ids); existing(g->prior_ids);
  combined.add(n, 2 * sizeof(irred_bytes) + sizeof(std::string));
  combined.add(g->native.selection_history().size(), 2 * sizeof(std::vector<irred_bytes>));
  combined.add(g->native.priors().size() + 1, sizeof(std::vector<irred_bytes>));
  for (const auto &history : g->native.selection_history()) {
    combined.add(history.kept_row_ids.size(), sizeof(irred_bytes));
    combined.add(history.complement_row_ids.size(), sizeof(irred_bytes));
  }
  for (const auto &prior_item : g->native.priors())
    combined.add(prior_item.applied_row_ids.size(), sizeof(irred_bytes));
  for (size_t i = 0; i < n; ++i)
    combined.add(std::max<size_t>(prior->ordered_ids.data[i].length, std::string{}.capacity()) + 1, 1);
  const auto envelope = combined.result();
  if (!peak || !envelope || *envelope > p->maximum_native_bytes)
    return IRRED_INVALID_INPUT;
  try {
    auto ids = copy(prior->ordered_ids);
    auto result = g->native.proper_offset(
        {prior->response.data, n}, ids, prior->mean, prior->variance,
        copy(prior->latent_identity), prior->independence_declared != 0,
        p->maximum_matrix_elements, p->maximum_forward_sensitivity);
    *semantic = static_cast<uint32_t>(result.status());
    *numerical_status = static_cast<uint32_t>(result.numerical_status());
    if (result.status() == DensityStatus::finite)
      *out = new irred_gaussian(std::move(result));
    return IRRED_OK;
  } catch (const std::bad_alloc &) {
    return IRRED_ALLOCATION_FAILURE;
  } catch (...) {
    return IRRED_EXCEPTION;
  }
}
extern "C" uint32_t irred_gaussian_source_view(const irred_gaussian *g,
                                               irred_gaussian_view *out) {
  if (!out)
    return IRRED_INVALID_INPUT;
  *out = {};
  if (!g)
    return IRRED_INVALID_INPUT;
  const auto &m = g->native.metadata();
  *out = {sizeof(*out),
          IRRED_ABI_VERSION,
          static_cast<uint32_t>(m.matrix_validation_scope),
          view(g->ids),
          view(m.measure),
          view(m.table_identity),
          view(m.uncertainty_identity),
          view(m.ordering_provenance),
          view(m.calibration_provenance),
          view(m.dependence_provenance),
          view(m.source_semantics),
          view(m.input_matrix_convention),
          view(m.treatment),
          view(g->native.mean_shift()),
          g->native.priors().size(),
          g->native.selection_history().size()};
  return IRRED_OK;
}
extern "C" uint32_t irred_gaussian_prior_view(const irred_gaussian *g,
                                              uint64_t index,
                                              irred_gaussian_prior *out) {
  if (!out)
    return IRRED_INVALID_INPUT;
  *out = {};
  if (!g || index >= g->native.priors().size())
    return IRRED_INVALID_INPUT;
  const auto &p = g->native.priors()[index];
  *out = {sizeof(*out),
          IRRED_ABI_VERSION,
          p.independence_assumed ? 1u : 0u,
          0,
          p.mean,
          p.variance,
          view(p.latent_identity),
          view(p.response),
          view(g->prior_ids[index])};
  return IRRED_OK;
}
extern "C" uint32_t irred_gaussian_evaluate(const irred_gaussian *g,
                                            const irred_gaussian_batch *b,
                                            const irred_gaussian_policy *p,
                                            irred_gaussian_result **out) {
  if (out)
    *out = nullptr;
  if (!out || !g || !b || !p)
    return IRRED_INVALID_INPUT;
  if (b->abi_version != IRRED_ABI_VERSION ||
      p->abi_version != IRRED_ABI_VERSION)
    return IRRED_ABI_MISMATCH;
  const auto n = g->ids.size();
  uint64_t left = p->maximum_string_bytes;
  // Reject the total row-major work shape before any copies or output
  // allocation.
  if (!valid(p) || b->struct_size != sizeof(*b) || b->reserved ||
      b->mode > IRRED_GAUSSIAN_MODE_PROFILE_OFFSET_SCORE || !n ||
      b->row_count > p->maximum_batch_elements / n)
    return IRRED_INVALID_INPUT;
  const auto total = b->row_count * n;
  if (!doubles(b->residuals, total) || b->residuals.length != total ||
      !strings(b->ordered_ids, n, left) || !doubles(b->response, n) ||
      b->response.length !=
          (b->mode == IRRED_GAUSSIAN_MODE_PROFILE_OFFSET_SCORE ? n : 0))
    return IRRED_INVALID_INPUT;
  const auto scratch = g->native.evaluation_payload_bound(
      static_cast<size_t>(b->row_count), b->mode == IRRED_GAUSSIAN_MODE_PROFILE_OFFSET_SCORE);
  const auto retained = g->native.retained_payload_bound();
  irred::detail::PayloadAccounting combined(0);
  if (!scratch || !retained) return IRRED_INVALID_INPUT;
  combined.add(*scratch, 1); combined.add(*retained, 1);
  combined.add(1, sizeof(irred_gaussian) - sizeof(irred::statistics::Gaussian));
  combined.add(g->ids.capacity(), sizeof(irred_bytes));
  auto views = [&](const auto &groups) {
    combined.add(groups.capacity(), sizeof(std::vector<irred_bytes>));
    for (const auto &ids : groups) combined.add(ids.capacity(), sizeof(irred_bytes));
  };
  views(g->kept_ids); views(g->complement_ids); views(g->prior_ids);
  combined.add(1, sizeof(irred_gaussian_result));
  combined.add(b->row_count, sizeof(irred_gaussian_row));
  combined.add(n, sizeof(std::string));
  for (size_t i = 0; i < n; ++i)
    combined.add(std::max<size_t>(b->ordered_ids.data[i].length, std::string{}.capacity()) + 1, 1);
  const auto envelope = combined.result();
  if (!envelope || *envelope > p->maximum_native_bytes) return IRRED_INVALID_INPUT;
  try {
    auto ids = copy(b->ordered_ids);
    auto result = std::make_unique<irred_gaussian_result>();
    result->rows.resize(b->row_count);
    for (uint64_t i = 0; i < b->row_count; ++i) {
      auto &row = result->rows[i];
      row.numerical_status = IRRED_NUMERICAL_STATUS_INVALID_INPUT;
      const std::span<const double> residual{b->residuals.data + i * n, n};
      if (b->mode == IRRED_GAUSSIAN_MODE_NORMALIZED_DENSITY) {
        auto r =
            g->native.evaluate(residual, ids, p->maximum_forward_sensitivity);
        row.status = static_cast<uint32_t>(r.density.status);
        row.numerical_status =
            static_cast<uint32_t>(r.density.numerical_status);
        if (r.density.status == DensityStatus::finite) {
          row.log_density = r.density.log_value;
          row.quadratic = r.quadratic;
          row.log_determinant = r.log_determinant;
          row.normalization = r.normalization;
          row.backward_residual = r.backward_residual;
          row.estimated_forward_sensitivity = r.estimated_forward_sensitivity;
        }
      } else {
        if (!g->native.priors().empty()) {
          row.status = IRRED_GAUSSIAN_STATUS_INCOMPATIBLE_METADATA;
          continue;
        }
        auto r = g->native.profile_offset(residual, {b->response.data, n}, ids,
                                          p->maximum_forward_sensitivity);
        row.status = static_cast<uint32_t>(r.status);
        row.numerical_status = static_cast<uint32_t>(r.numerical_status);
        if (r.status == DensityStatus::finite) {
          row.coefficient = r.coefficient;
          row.quadratic = r.quadratic;
          row.backward_residual = r.backward_residual;
          row.estimated_forward_sensitivity = r.estimated_forward_sensitivity;
        }
      }
    }
    *out = result.release();
    return IRRED_OK;
  } catch (const std::bad_alloc &) {
    return IRRED_ALLOCATION_FAILURE;
  } catch (...) {
    return IRRED_EXCEPTION;
  }
}
extern "C" uint32_t irred_gaussian_result_view(const irred_gaussian_result *r,
                                               const irred_gaussian_row **rows,
                                               uint64_t *n) {
  if (!rows || !n)
    return IRRED_INVALID_INPUT;
  *rows = nullptr;
  *n = 0;
  if (!r)
    return IRRED_INVALID_INPUT;
  *rows = r->rows.data();
  *n = r->rows.size();
  return IRRED_OK;
}
extern "C" uint32_t irred_gaussian_destroy(irred_gaussian *g) {
  delete g;
  return IRRED_OK;
}
extern "C" uint32_t irred_gaussian_result_destroy(irred_gaussian_result *r) {
  delete r;
  return IRRED_OK;
}

extern "C" uint32_t
irred_gaussian_selection_view(const irred_gaussian *g, uint64_t index,
                              irred_gaussian_selection *out) {
  if (!out)
    return IRRED_INVALID_INPUT;
  *out = {};
  if (!g || index >= g->native.selection_history().size())
    return IRRED_INVALID_INPUT;
  *out = {sizeof(*out), IRRED_ABI_VERSION,
          view(g->native.selection_history()[index].operation),
          view(g->kept_ids[index]), view(g->complement_ids[index])};
  return IRRED_OK;
}
