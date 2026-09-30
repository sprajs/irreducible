#include "irred/supernova.hpp"
#include <algorithm>
#include <cmath>
#include <type_traits>
#include <utility>
namespace irred::supernova {
Consumer prepare(observations::Prepared observations, Policy policy) {
  const auto &source = observations.source();
  if (source.values.size() > policy.background.maximum_queries) {
    Consumer failed;
    failed.status_ = Status::work_limit;
    return failed;
  }
  std::vector<cosmology::Query> queries;
  if (source.zhd.size() == source.values.size() &&
      source.zhel.size() == source.values.size())
    for (std::size_t i = 0; i < source.values.size(); ++i)
      queries.push_back({source.zhd[i], source.zhel[i],
                         cosmology::Convention::released_zhd_zhel});
  auto ids = source.measurement_ids;
  return Consumer::prepare_common(std::move(observations), ids, queries, false,
                                  policy);
}
Consumer prepare_synthetic(observations::Prepared observations,
                           std::span<const std::string> ids,
                           std::span<const cosmology::Query> queries,
                           Policy policy) {
  return Consumer::prepare_common(std::move(observations), ids, queries, true,
                                  policy);
}
Consumer Consumer::prepare_common(observations::Prepared observations,
                                  std::span<const std::string> ids,
                                  std::span<const cosmology::Query> queries,
                                  bool synthetic, Policy policy) {
  Consumer out;
  if (!(policy.maximum_forward_sensitivity > 0) ||
      !std::isfinite(policy.maximum_forward_sensitivity))
    return out;
  out.synthetic_ = synthetic;
  out.observations_ = std::move(observations);
  const auto &source = out.observations_.source();
  const auto profile = synthetic
                           ? observations::Profile::gaussian_fixture_v1
                           : observations::Profile::pantheon_plus_released_v1;
  const auto role = synthetic ? observations::Role::synthetic_control
                              : observations::Role::released_fitted_summary;
  if (out.observations_.status() != observations::Status::ok ||
      source.profile != profile || source.role != role ||
      source.unit != observations::Unit::magnitude ||
      source.uncertainty != observations::Uncertainty::covariance ||
      ids.size() != source.measurement_ids.size() ||
      !std::equal(ids.begin(), ids.end(), source.measurement_ids.begin()) ||
      queries.size() != source.values.size()) {
    out.status_ = Status::incompatible_metadata;
    return out;
  }
  if (queries.size() > policy.background.maximum_queries) {
    out.status_ = Status::work_limit;
    return out;
  }
  out.source_queries_.assign(queries.begin(), queries.end());
  if (!synthetic) {
    auto selected =
        out.observations_.select(observations::Selection::pantheon_zhd_gt_001);
    if (selected.status != observations::Status::ok)
      return out;
    out.indices_ = std::move(selected.source_indices);
  } else {
    for (std::size_t i = 0; i < queries.size(); ++i) {
      if (!source.source_selection.empty() && !source.source_selection[i])
        continue;
      if (!std::isfinite(queries[i].z_expansion))
        return out;
      if (queries[i].z_expansion > .01)
        out.indices_.push_back(i);
    }
  }
  if (out.indices_.empty())
    return out;
  if (out.indices_.size() > policy.background.maximum_queries) {
    out.status_ = Status::work_limit;
    return out;
  }
  for (auto i : out.indices_) {
    if (queries[i].convention != cosmology::Convention::released_zhd_zhel ||
        !std::isfinite(queries[i].z_observer) ||
        !(queries[i].z_observer > -1)) {
      out.status_ = Status::incompatible_metadata;
      return out;
    }
    out.observed_.push_back(source.values[i]);
    out.ones_.push_back(1);
    out.queries_.push_back(queries[i]);
  }
  // Shared selected-observation kernel: no consumer covariance repair/math.
  auto gaussian = statistics::prepare_selected_observations(
      out.observations_, out.indices_, policy.maximum_matrix_elements,
      policy.maximum_forward_sensitivity, policy.arithmetic);
  out.preparation_status_ = gaussian.status();
  out.preparation_numerical_status_ = gaussian.numerical_status();
  if (gaussian.status() != statistics::DensityStatus::finite) {
    out.status_ = Status::numerical_failure;
    return out;
  }
  out.ids_ = gaussian.metadata().ordered_ids;
  out.profile_ = std::move(gaussian).prepare_offset_profile(
      out.ones_, out.ids_, policy.maximum_forward_sensitivity);
  out.preparation_status_ = out.profile_.status();
  out.preparation_numerical_status_ = out.profile_.numerical_status();
  if (out.profile_.status() != statistics::DensityStatus::finite) {
    out.status_ = Status::numerical_failure;
    return out;
  }
  out.status_ = Status::ok;
  return out;
}
BatchResult Consumer::evaluate_batch(std::span<const ModelPoint> points,
                                     Policy policy) const {
  return evaluate_common<ModelPoint, BatchResult>(points, policy);
}
BatchResultV2 Consumer::evaluate_batch_v2(std::span<const ModelPointV2> points,
                                          Policy policy) const {
  return evaluate_common<ModelPointV2, BatchResultV2>(points, policy);
}
PiecewiseBatchResult
Consumer::evaluate_piecewise_batch(std::span<const PiecewiseModelPoint> points,
                                   PiecewiseEvaluationPolicy policy) const {
  return evaluate_common<PiecewiseModelPoint, PiecewiseBatchResult>(points,
                                                                    policy);
}
template <class Point, class Result, class EvaluationPolicy>
Result Consumer::evaluate_common(std::span<const Point> points,
                                 EvaluationPolicy policy) const {
  Result out;
  if (!(policy.maximum_forward_sensitivity > 0) ||
      !std::isfinite(policy.maximum_forward_sensitivity))
    return out;
  if (status_ != Status::ok) {
    out.status = status_;
    return out;
  }
  if (policy.arithmetic != profile_.arithmetic()) {
    out.status = Status::incompatible_metadata;
    return out;
  }
  if (points.size() > policy.maximum_models) {
    out.status = Status::work_limit;
    return out;
  }
  out.status = Status::ok;
  out.slots.reserve(points.size());
  auto remaining_evaluations = [&]() {
    if constexpr (std::is_same_v<Point, PiecewiseModelPoint>)
      return policy.background.maximum_segment_visits;
    else
      return policy.background.maximum_total_evaluations;
  }();
  for (const auto &point : points) {
    auto row = [&]() {
      if constexpr (std::is_same_v<Point, PiecewiseModelPoint>)
        return PiecewiseSlot{point};
      else if constexpr (std::is_same_v<Point, ModelPoint>)
        return Slot{};
      else
        return SlotV2{point, Slot{}};
    }();
    auto &slot = [&]() -> auto & {
      if constexpr (std::is_same_v<Point, ModelPoint> ||
                    std::is_same_v<Point, PiecewiseModelPoint>)
        return row;
      else
        return row.calculation;
    }();
    if constexpr (!std::is_same_v<Point, PiecewiseModelPoint>)
      slot.source = {point.model, point.omega_m, point.constant_q};
    if constexpr (std::is_same_v<Point, ModelPointV2>)
      if (point.model == cosmology::Model::flat_cpl_late_v1 &&
          point.constant_q != 0) {
        out.slots.push_back(std::move(row));
        continue;
      }
    auto background = [&]() {
      if constexpr (std::is_same_v<Point, PiecewiseModelPoint>) {
        return cosmology::prepare_piecewise_q(
            {computational_h0_km_s_mpc, point.q});
      } else if constexpr (std::is_same_v<Point, ModelPointV2>) {
        if (point.model == cosmology::Model::flat_cpl_late_v1)
          return cosmology::prepare_cpl(
              {computational_h0_km_s_mpc, point.omega_m, point.w0, point.wa});
        // Keep attempted inactive fields: legacy prepare rejects noncanonical
        // values instead of silently replacing them with defaults.
        return cosmology::prepare({point.model, computational_h0_km_s_mpc,
                                   point.omega_m, point.constant_q, point.w0,
                                   point.wa});
      } else {
        // The existing entry remains strict: prepare does not admit model 2.
        return cosmology::prepare({point.model, computational_h0_km_s_mpc,
                                   point.omega_m, point.constant_q});
      }
    }();
    slot.background_status = background.status();
    if (background.status() != cosmology::Status::ok) {
      out.slots.push_back(std::move(row));
      continue;
    }
    if (remaining_evaluations <
        (std::is_same_v<Point, PiecewiseModelPoint> ? 1 : 3)) {
      slot.status = Status::work_limit;
      slot.background_status = cosmology::Status::work_limit;
      slot.numerical_status = numerics::Status::work_limit;
      out.slots.push_back(std::move(row));
      continue;
    }
    auto background_policy = policy.background;
    if constexpr (std::is_same_v<Point, PiecewiseModelPoint>)
      background_policy.maximum_segment_visits = remaining_evaluations;
    else
      background_policy.maximum_total_evaluations = remaining_evaluations;
    auto prediction = background.evaluate_batch(queries_, background_policy);
    slot.background_status = prediction.status;
    if constexpr (std::is_same_v<Point, PiecewiseModelPoint>) {
      slot.segment_visits = prediction.segments_processed;
      out.segment_visits += slot.segment_visits;
      remaining_evaluations -= slot.segment_visits;
    } else {
      for (const auto &value : prediction.slots)
        slot.background_evaluations += value.evaluations;
      remaining_evaluations -= slot.background_evaluations;
    }
    if (prediction.status != cosmology::Status::ok) {
      if constexpr (std::is_same_v<Point, PiecewiseModelPoint>)
        if (prediction.status == cosmology::Status::work_limit)
          slot.numerical_status = numerics::Status::work_limit;
      slot.status = prediction.status == cosmology::Status::work_limit
                        ? Status::work_limit
                        : Status::numerical_failure;
      out.slots.push_back(std::move(row));
      continue;
    }
    bool failed = false;
    for (std::size_t i = 0; i < prediction.slots.size(); ++i) {
      const auto &value = prediction.slots[i];
      slot.numerical_status = value.numerical_status;
      slot.background_status = value.status;
      const double luminosity_shape = [&]() {
        if constexpr (std::is_same_v<Point, PiecewiseModelPoint>)
          return value.geometry.dimensionless_luminosity_shape;
        else
          return value.dimensionless_luminosity_shape;
      }();
      if (value.status != cosmology::Status::ok || !(luminosity_shape > 0)) {
        slot.status = value.status == cosmology::Status::work_limit
                          ? Status::work_limit
                          : Status::numerical_failure;
        failed = true;
        break;
      }
      const auto shape = 5 * std::log10(luminosity_shape);
      slot.shape_magnitudes.push_back(shape);
      slot.base_residuals.push_back(observed_[i] - shape);
    }
    if (!failed) {
      const auto profile = profile_.evaluate(
          slot.base_residuals, ids_, policy.maximum_forward_sensitivity);
      slot.solve_diagnostics = profile;
      slot.profile_status = profile.status;
      slot.numerical_status = profile.numerical_status;
      if (profile.status == statistics::DensityStatus::finite) {
        slot.offset_coefficient = profile.coefficient;
        slot.quadratic = profile.quadratic;
        slot.relative_profile_score = -.5 * profile.quadratic;
        slot.profiled_residuals = profile.adjusted_residuals;
        slot.status = Status::ok;
      } else
        slot.status = Status::numerical_failure;
    }
    if (slot.status != Status::ok) {
      slot.shape_magnitudes.clear();
      slot.base_residuals.clear();
      slot.profiled_residuals.clear();
    }
    out.slots.push_back(std::move(row));
  }
  return out;
}
} // namespace irred::supernova
