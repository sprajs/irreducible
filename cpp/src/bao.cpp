#include "irred/bao.hpp"
#include "piecewise_radial.hpp"
#include <cmath>
#include <limits>
#include <type_traits>
namespace irred::bao {
namespace {
bool normal_or_zero(long double x) {
  if (!std::isfinite(x) || std::abs(x) > std::numeric_limits<double>::max())
    return false;
  const auto y = static_cast<double>(x);
  return x == 0 || (y != 0 && std::fpclassify(y) == FP_NORMAL);
}
bool known(Observable o) {
  return o == Observable::transverse_over_ruler ||
         o == Observable::hubble_over_ruler ||
         o == Observable::volume_over_ruler;
}
} // namespace
namespace {
long double project(long double scale, long double expansion,
                    long double radial, Query query) {
  const auto dm = scale * radial;
  const auto dh = scale / expansion;
  long double value = 0;
  switch (query.observable) {
  case Observable::transverse_over_ruler:
    value = dm;
    break;
  case Observable::hubble_over_ruler:
    value = dh;
    break;
  case Observable::volume_over_ruler:
    value = std::cbrt(static_cast<long double>(query.z) * dm * dm * dh);
    break;
  }
  return value;
}
template <class Background, class EvaluationPolicy>
auto evaluate_observables(const Background &background, Ruler ruler,
                          std::span<const Query> queries,
                          EvaluationPolicy policy) {
  constexpr bool analytic =
      std::is_same_v<Background, cosmology::PiecewiseBackground>;
  using Result = std::conditional_t<analytic, PiecewiseBatch, Batch>;
  using OutputSlot = std::conditional_t<analytic, PiecewiseSlot, Slot>;
  using NativeSlot =
      std::conditional_t<analytic, cosmology::PiecewiseSlot, cosmology::Slot>;
  Result result;
  // Resource checks precede dereferencing or copying any caller query.
  if (queries.size() > policy.maximum_queries ||
      queries.size() > policy.background.maximum_queries ||
      queries.size() > policy.maximum_native_bytes /
                           (sizeof(OutputSlot) + sizeof(cosmology::Query) +
                            sizeof(std::size_t) + sizeof(NativeSlot))) {
    result.status = cosmology::Status::work_limit;
    return result;
  }
  if (!std::isfinite(ruler.h0_rd_km_s))
    return result;
  if (ruler.h0_rd_km_s < 5000 || ruler.h0_rd_km_s > 15000) {
    result.status = cosmology::Status::unsupported_domain;
    return result;
  }
  if (background.status() != cosmology::Status::ok) {
    result.status = background.status();
    return result;
  }
  std::vector<cosmology::Query> admitted;
  std::vector<std::size_t> indices;
  admitted.reserve(queries.size());
  indices.reserve(queries.size());
  result.slots.resize(queries.size());
  for (std::size_t i = 0; i < queries.size(); ++i) {
    result.slots[i].source = queries[i];
    if (!known(queries[i].observable))
      continue;
    admitted.push_back({queries[i].z, queries[i].z,
                        cosmology::Convention::geometric_same_redshift});
    indices.push_back(i);
  }
  // This is one coarse batch, including repeated redshifts as distinct slots.
  auto computed = [&] {
    if constexpr (analytic)
      return cosmology::detail::PiecewiseRadialAccess::radial(
          background, admitted, policy.background);
    else
      return background.evaluate_batch(admitted, policy.background);
  }();
  if (computed.status != cosmology::Status::ok) {
    result.status = computed.status;
    result.slots.clear();
    return result;
  }
  const auto scale = static_cast<long double>(irred::speed_of_light_m_per_s) /
                     1000 / ruler.h0_rd_km_s;
  for (std::size_t j = 0; j < indices.size(); ++j) {
    const auto &source = computed.slots[j];
    auto &slot = result.slots[indices[j]];
    slot.status = source.status;
    slot.numerical_status = source.numerical_status;
    if constexpr (analytic) {
      slot.segment_visits = source.segments_processed;
      result.segment_visits += source.segments_processed;
    } else {
      slot.evaluations = source.evaluations;
      result.evaluations += source.evaluations;
    }
    if (source.status != cosmology::Status::ok)
      continue;
    const auto value = [&] {
      if constexpr (analytic)
        return project(scale, source.expansion_E, source.radial_integral,
                       slot.source);
      else
        return project(scale, source.expansion_E, source.radial_integral,
                       slot.source);
    }();
    if (!normal_or_zero(value) || value < 0 ||
        ((slot.source.z > 0 ||
          slot.source.observable == Observable::hubble_over_ruler) &&
         !(value > 0))) {
      slot.status = cosmology::Status::numerical_failure;
      slot.numerical_status =
          !std::isfinite(value) ||
                  std::abs(value) > std::numeric_limits<double>::max()
              ? numerics::Status::overflow
              : numerics::Status::outside_domain;
      continue;
    }
    slot.dimensionless_value = static_cast<double>(value);
  }
  result.status = cosmology::Status::ok;
  return result;
}
} // namespace
Batch evaluate(const cosmology::Background &background, Ruler ruler,
               std::span<const Query> queries, Policy policy) {
  return evaluate_observables(background, ruler, queries, policy);
}
PiecewiseBatch
evaluate_piecewise(const cosmology::PiecewiseBackground &background,
                   Ruler ruler, std::span<const Query> queries,
                   PiecewisePolicy policy) {
  return evaluate_observables(background, ruler, queries, policy);
}
namespace {
bool valid_policy(const DensityPolicy &p) {
  return std::isfinite(p.maximum_forward_sensitivity) &&
         p.maximum_forward_sensitivity > 0 &&
         (p.arithmetic == numerics::Arithmetic::binary64_legacy_v1 ||
          p.arithmetic == numerics::Arithmetic::longdouble_cpu_v1);
}
bool valid_policy(const PiecewiseDensityPolicy &p) {
  return std::isfinite(p.maximum_forward_sensitivity) &&
         p.maximum_forward_sensitivity > 0 &&
         (p.arithmetic == numerics::Arithmetic::binary64_legacy_v1 ||
          p.arithmetic == numerics::Arithmetic::longdouble_cpu_v1);
}
bool fits(std::size_t count, std::size_t element, std::size_t &remaining) {
  if (count > remaining / element)
    return false;
  remaining -= count * element;
  return true;
}
} // namespace
PreparedDensity prepare_density(const DensityInput &input,
                                DensityPolicy policy) {
  PreparedDensity out;
  const auto n = input.queries.size();
  if (!valid_policy(policy))
    return out;
  if (n > policy.observables.maximum_queries || n == 0 ||
      n > policy.maximum_matrix_elements / n) {
    out.numerical_status_ = numerics::Status::work_limit;
    return out;
  }
  if (input.observed.size() != n || input.ordered_ids.size() != n ||
      input.covariance.size() != n * n)
    return out;
  // Include retained source, Gaussian covariance/factor and bounded preparation
  // vectors, in addition to identity strings; no unbounded source copy first.
  auto bytes = policy.maximum_native_bytes;
  if (!fits(n, sizeof(Query) + sizeof(double) + 2 * sizeof(std::string),
            bytes) ||
      !fits(n * n, 3 * sizeof(double) + sizeof(long double), bytes) ||
      !fits(n, 8 * sizeof(long double) + 4 * sizeof(double), bytes)) {
    out.numerical_status_ = numerics::Status::work_limit;
    return out;
  }
  constexpr std::string_view measure =
      "product of dimensionless BAO released distance-ratio coordinates";
  constexpr std::string_view semantics =
      "released fitted BAO distance summaries; covariance dimensionless ratio "
      "squared";
  constexpr std::string_view matrix =
      "full source covariance in declared exact row order";
  constexpr std::string_view treatment = "normalized Gaussian";
  constexpr std::string_view arithmetic = "F02/longdouble-cpu/v1";
  if (!fits(1, sizeof(PreparedDensity), bytes) ||
      !fits(measure.size() + semantics.size() + matrix.size() +
                treatment.size() + arithmetic.size() + 5,
            1, bytes)) {
    out.numerical_status_ = numerics::Status::work_limit;
    return out;
  }
  auto string_bytes = policy.maximum_string_bytes;
  for (const auto *text :
       {&input.table_identity, &input.covariance_identity,
        &input.ordering_provenance, &input.calibration_provenance,
        &input.dependence_provenance}) {
    if (text->empty())
      return out;
    if (!fits(text->size(), 1, string_bytes) ||
        !fits(text->size() + 1, 2, bytes)) {
      out.numerical_status_ = numerics::Status::work_limit;
      return out;
    }
  }
  if ((input.role != RowRole::released_fitted_distance_summary &&
       input.role != RowRole::synthetic_control) ||
      input.covariance_unit != CovarianceUnit::dimensionless_ratio_squared)
    return out;
  for (std::size_t i = 0; i < n; ++i) {
    if (input.ordered_ids[i].empty() || !known(input.queries[i].observable) ||
        !std::isfinite(input.queries[i].z) || input.queries[i].z < 0 ||
        input.queries[i].z > 5 || !std::isfinite(input.observed[i]))
      return out;
    if (!fits(input.ordered_ids[i].size(), 1, string_bytes) ||
        !fits(input.ordered_ids[i].size() + 1, 2, bytes)) {
      out.numerical_status_ = numerics::Status::work_limit;
      return out;
    }
    for (std::size_t j = 0; j < i; ++j)
      if (input.ordered_ids[i] == input.ordered_ids[j])
        return out;
  }
  statistics::Metadata metadata;
  metadata.ordered_ids = input.ordered_ids;
  metadata.measure =
      "product of dimensionless BAO released distance-ratio coordinates";
  metadata.table_identity = input.table_identity;
  metadata.uncertainty_identity = input.covariance_identity;
  metadata.ordering_provenance = input.ordering_provenance;
  metadata.calibration_provenance = input.calibration_provenance;
  metadata.dependence_provenance = input.dependence_provenance;
  metadata.source_semantics =
      input.role == RowRole::synthetic_control
          ? "synthetic BAO control; covariance dimensionless ratio squared"
          : "released fitted BAO distance summaries; covariance dimensionless "
            "ratio squared";
  metadata.input_matrix_convention =
      "full source covariance in declared exact row order";
  out.gaussian_ = statistics::prepare_gaussian(
      input.covariance, statistics::MatrixKind::covariance, std::move(metadata),
      policy.maximum_matrix_elements, policy.maximum_forward_sensitivity,
      policy.arithmetic);
  out.status_ = out.gaussian_.status();
  out.numerical_status_ = out.gaussian_.numerical_status();
  out.source_ = input;
  return out;
}
DensityBatch PreparedDensity::evaluate(std::span<const ModelQuery> models,
                                       DensityPolicy policy) const {
  return evaluate_common<ModelQuery, DensityBatch>(models, policy);
}
PiecewiseDensityBatch
PreparedDensity::evaluate_piecewise(std::span<const PiecewiseModelPoint> models,
                                    PiecewiseDensityPolicy policy) const {
  return evaluate_common<PiecewiseModelPoint, PiecewiseDensityBatch>(models,
                                                                     policy);
}
template <class Model, class Result, class EvaluationPolicy>
Result PreparedDensity::evaluate_common(std::span<const Model> models,
                                        EvaluationPolicy policy) const {
  constexpr bool analytic = std::is_same_v<Model, PiecewiseModelPoint>;
  using OutputSlot =
      std::conditional_t<analytic, PiecewiseDensitySlot, DensitySlot>;
  using NativeSlot =
      std::conditional_t<analytic, cosmology::PiecewiseSlot, cosmology::Slot>;
  Result batch;
  if (status_ != statistics::DensityStatus::finite) {
    batch.status = status_;
    batch.numerical_status = numerical_status_;
    return batch;
  }
  if (!valid_policy(policy) || policy.arithmetic != gaussian_.arithmetic())
    return batch;
  const auto n = source_.queries.size();
  auto bytes = policy.maximum_native_bytes;
  if (models.size() > policy.maximum_models ||
      !fits(models.size(), sizeof(OutputSlot), bytes) ||
      (models.size() != 0 &&
       (n > std::numeric_limits<std::size_t>::max() / models.size() ||
        !fits(n * models.size(), 2 * sizeof(double), bytes))) ||
      !fits(n,
            sizeof(std::conditional_t<analytic, PiecewiseSlot, Slot>) +
                sizeof(NativeSlot) + sizeof(cosmology::Query) +
                sizeof(std::size_t) + 4 * sizeof(long double) +
                2 * sizeof(double),
            bytes)) {
    batch.numerical_status = numerics::Status::work_limit;
    return batch;
  }
  for (const auto &model : models) {
    const auto id = [&]() -> std::string_view {
      if constexpr (analytic)
        return cosmology::PiecewiseBackground::model_id;
      else
        return model.background.model_id();
    }();
    if (!fits(id.size() + 1, 1, bytes) ||
        !fits(bao::equation_id.size() + 1, 1, bytes) ||
        !fits(gaussian_.metadata().arithmetic_id.size() + 1, 1, bytes)) {
      batch.numerical_status = numerics::Status::work_limit;
      return batch;
    }
  }
  if constexpr (analytic)
    batch.slots.reserve(models.size());
  else
    batch.slots.resize(models.size());
  auto remaining = [&] {
    if constexpr (analytic)
      return policy.observables.background.maximum_segment_visits;
    else
      return policy.observables.background.maximum_total_evaluations;
  }();
  for (std::size_t i = 0; i < models.size(); ++i) {
    if constexpr (analytic)
      batch.slots.emplace_back(models[i]);
    auto &slot = batch.slots[i];
    if constexpr (analytic)
      slot.model_id = cosmology::PiecewiseBackground::model_id;
    else {
      slot.attempted_model = models[i].background.parameters();
      slot.attempted_h0_rd_km_s = models[i].ruler.h0_rd_km_s;
      slot.model_id = models[i].background.model_id();
    }
    slot.equation_id = bao::equation_id;
    slot.arithmetic_id = gaussian_.metadata().arithmetic_id;
    auto local = policy.observables;
    if constexpr (analytic)
      local.background.maximum_segment_visits = remaining;
    else
      local.background.maximum_total_evaluations = remaining;
    auto predicted = [&] {
      if constexpr (analytic) {
        auto background = cosmology::prepare_piecewise_q(
            {piecewise_computational_h0_km_s_mpc, models[i].q});
        return bao::evaluate_piecewise(background, models[i].ruler,
                                       source_.queries, local);
      } else
        return bao::evaluate(models[i].background, models[i].ruler,
                             source_.queries, local);
    }();
    slot.background_status = predicted.status;
    if constexpr (analytic) {
      slot.segment_visits = predicted.segment_visits;
      batch.segment_visits += predicted.segment_visits;
      remaining -= predicted.segment_visits;
    } else {
      slot.evaluations = predicted.evaluations;
      batch.evaluations += predicted.evaluations;
      remaining -= predicted.evaluations;
    }
    if (predicted.status != cosmology::Status::ok) {
      slot.result.density.status = statistics::DensityStatus::numerical_failure;
      slot.numerical_status = predicted.status == cosmology::Status::work_limit
                                  ? numerics::Status::work_limit
                                  : numerics::Status::invalid_input;
      slot.result.density.numerical_status = slot.numerical_status;
      continue;
    }
    slot.predictions.reserve(n);
    slot.residuals.reserve(n);
    bool finite = true;
    for (std::size_t j = 0; j < n; ++j) {
      const auto &row = predicted.slots[j];
      if (row.status != cosmology::Status::ok) {
        slot.background_status = row.status;
        slot.numerical_status = row.numerical_status;
        finite = false;
        break;
      }
      slot.predictions.push_back(row.dimensionless_value);
      const auto residual = static_cast<long double>(source_.observed[j]) -
                            row.dimensionless_value;
      if (!normal_or_zero(residual)) {
        slot.numerical_status =
            std::abs(residual) > std::numeric_limits<double>::max()
                ? numerics::Status::overflow
                : numerics::Status::outside_domain;
        finite = false;
        break;
      }
      slot.residuals.push_back(static_cast<double>(residual));
    }
    if (finite) {
      slot.result = gaussian_.evaluate(slot.residuals, source_.ordered_ids,
                                       policy.maximum_forward_sensitivity);
      slot.numerical_status = slot.result.density.numerical_status;
      finite = slot.result.density.status == statistics::DensityStatus::finite;
    } else {
      slot.result.density.status = statistics::DensityStatus::numerical_failure;
      slot.result.density.numerical_status = slot.numerical_status;
    }
    if (!finite) {
      slot.predictions.clear();
      slot.residuals.clear();
    }
  }
  batch.status = statistics::DensityStatus::finite;
  batch.numerical_status = numerics::Status::ok;
  return batch;
}
} // namespace irred::bao
