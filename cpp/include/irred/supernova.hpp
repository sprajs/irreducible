#pragma once
#include "irred/background.hpp"
#include "irred/observations.hpp"
#include "irred/piecewise_background.hpp"
#include "irred/statistics.hpp"
#include <cstddef>
#include <span>
#include <string>
#include <vector>
namespace irred::supernova {
enum class Status {
  ok,
  invalid_input,
  incompatible_metadata,
  numerical_failure,
  work_limit
};
// H0 is absorbed by the free magnitude offset, not an inferred parameter.
struct ModelPoint {
  cosmology::Model model = cosmology::Model::flat_lcdm_late_v1;
  double omega_m = .3;
  double constant_q = 0;
};
// Version 2 admission is explicit: every active and inactive field is supplied.
// Legacy callers cannot accidentally select CPL through implicit defaults.
struct ModelPointV2 {
  cosmology::Model model;
  double omega_m, constant_q, w0, wa;
  ModelPointV2(cosmology::Model model_id, double matter, double q,
               double present_w, double evolution_w)
      : model(model_id), omega_m(matter), constant_q(q), w0(present_w),
        wa(evolution_w) {}
};
struct Policy {
  cosmology::Policy background;
  numerics::Arithmetic arithmetic = numerics::Arithmetic::binary64_legacy_v1;
  std::size_t maximum_models = 64;
  std::size_t maximum_matrix_elements = 2893401;
  double maximum_forward_sensitivity = 0; // required explicit positive policy
};
struct Slot {
  ModelPoint source;
  Status status = Status::invalid_input;
  cosmology::Status background_status = cosmology::Status::invalid_input;
  numerics::Status numerical_status = numerics::Status::invalid_input;
  statistics::ProfileResult solve_diagnostics;
  statistics::DensityStatus profile_status =
      statistics::DensityStatus::invalid_input;
  double offset_coefficient = 0, quadratic = 0, relative_profile_score = 0;
  // Ordered selected predictions are 5log10((1+zHEL)*I(zHD)); no intercept/H0.
  std::vector<double> shape_magnitudes, base_residuals, profiled_residuals;
  std::size_t background_evaluations = 0;
};
struct BatchResult {
  Status status = Status::invalid_input;
  std::vector<Slot> slots;
};
struct SlotV2 {
  ModelPointV2 source; // Full attempted parameters, including failed rows.
  Slot calculation; // source here is the three-field compatibility projection.
};
struct BatchResultV2 {
  Status status = Status::invalid_input;
  std::vector<SlotV2> slots;
};
// Analytic fixed-five-bin provider; no defaults, priors or legacy model tag.
struct PiecewiseModelPoint {
  std::array<double, 5> q;
  explicit PiecewiseModelPoint(std::array<double, 5> values) : q(values) {}
};
struct PiecewiseEvaluationPolicy {
  cosmology::PiecewisePolicy background;
  numerics::Arithmetic arithmetic =
      static_cast<numerics::Arithmetic>(UINT32_MAX);
  std::size_t maximum_models = 64;
  double maximum_forward_sensitivity = 0;
};
struct PiecewiseSlot {
  explicit PiecewiseSlot(PiecewiseModelPoint point) : source(point) {}
  PiecewiseModelPoint source;
  Status status = Status::invalid_input;
  cosmology::Status background_status = cosmology::Status::invalid_input;
  numerics::Status numerical_status = numerics::Status::invalid_input;
  statistics::ProfileResult solve_diagnostics;
  statistics::DensityStatus profile_status =
      statistics::DensityStatus::invalid_input;
  double offset_coefficient = 0, quadratic = 0, relative_profile_score = 0;
  std::vector<double> shape_magnitudes, base_residuals, profiled_residuals;
  std::size_t segment_visits = 0;
};
struct PiecewiseBatchResult {
  std::size_t segment_visits = 0;
  Status status = Status::invalid_input;
  std::vector<PiecewiseSlot> slots;
};
class Consumer {
public:
  Status status() const noexcept { return status_; }
  numerics::Status preparation_numerical_status() const noexcept {
    return preparation_numerical_status_;
  }
  statistics::DensityStatus preparation_status() const noexcept {
    return preparation_status_;
  }
  const observations::Prepared &observations() const noexcept {
    return observations_;
  }
  std::span<const std::size_t> selected_source_indices() const noexcept {
    return indices_;
  }
  const statistics::Metadata &probability_metadata() const noexcept {
    return profile_.metadata();
  }
  const char *query_provenance() const noexcept {
    return synthetic_
               ? "generated native synthetic queries; explicit source-order IDs"
               : "released table zHD/zHEL source fields";
  }
  const statistics::ProfileOperator &profile_operator() const noexcept {
    return profile_;
  }
  std::span<const cosmology::Query> source_queries() const noexcept {
    return source_queries_;
  }
  std::span<const cosmology::Query> selected_queries() const noexcept {
    return queries_;
  }
  static constexpr const char *prediction_unit = "magnitude";
  static constexpr const char *offset_convention =
      "unconstrained additive magnitude intercept; H0 absorbed";
  static constexpr const char *shape_convention =
      "zHD radial integral with zHEL luminosity prefactor";
  BatchResult evaluate_batch(std::span<const ModelPoint>, Policy) const;
  BatchResultV2 evaluate_batch_v2(std::span<const ModelPointV2>, Policy) const;
  PiecewiseBatchResult
      evaluate_piecewise_batch(std::span<const PiecewiseModelPoint>,
                               PiecewiseEvaluationPolicy) const;
  static constexpr double computational_h0_km_s_mpc = 70;
  static constexpr const char *score_id =
      "W01/intercept-free-profile-relative-score/v1";

private:
  template <class Point, class Result, class EvaluationPolicy>
  Result evaluate_common(std::span<const Point>, EvaluationPolicy) const;
  Status status_ = Status::invalid_input;
  bool synthetic_ = false;
  numerics::Status preparation_numerical_status_ =
      numerics::Status::invalid_input;
  statistics::DensityStatus preparation_status_ =
      statistics::DensityStatus::invalid_input;
  observations::Prepared observations_;
  statistics::ProfileOperator profile_;
  std::vector<std::size_t> indices_;
  std::vector<std::string> ids_;
  std::vector<double> observed_, ones_;
  std::vector<cosmology::Query> queries_, source_queries_;
  friend Consumer prepare(observations::Prepared, Policy);
  friend Consumer prepare_synthetic(observations::Prepared,
                                    std::span<const std::string>,
                                    std::span<const cosmology::Query>, Policy);
  static Consumer prepare_common(observations::Prepared,
                                 std::span<const std::string>,
                                 std::span<const cosmology::Query>, bool,
                                 Policy);
};
// Concrete released Pantheon profile, source mask and zHD>0.01 selection.
// Full raw source remains owned. Selected covariance only is validated; no
// validity claim for the excluded/full probability object. No fitting/sampling,
// normalized-density claim, absolute calibration or H0 identification.
Consumer prepare(observations::Prepared, Policy);
// Native comparison controls: Gaussian synthetic role, explicitly ordered
// source IDs and typed redshifts; same source mask and strict expansion-z>0.01
// rule. Query provenance is generated synthetic input, never a release
// observation.
Consumer prepare_synthetic(observations::Prepared,
                           std::span<const std::string> ordered_source_ids,
                           std::span<const cosmology::Query>, Policy);
} // namespace irred::supernova
