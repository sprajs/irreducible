#pragma once
#include "irred/numerics.hpp"
#include "irred/observations.hpp"
#include <cstdint>
#include <span>
#include <string>
#include <vector>
namespace irred::statistics {
enum class DensityStatus {
  finite,
  outside_support,
  invalid_input,
  unsupported_domain,
  numerical_failure,
  incompatible_metadata
};
struct DensityResult {
  DensityStatus status = DensityStatus::invalid_input;
  double log_value = 0;
  numerics::Status numerical_status = numerics::Status::invalid_input;
};
DensityResult normal_log_density(double x, double mean,
                                 double standard_deviation) noexcept;
DensityResult
poisson_log_mass(std::uint64_t count,
                 double rate) noexcept; // first qualified count domain 0..99
DensityResult selected_standard_normal_positive(
    double x) noexcept; // strict x>0; selection probability exactly 1/2
struct GaussianResult {
  DensityResult density;
  double quadratic = 0, log_determinant = 0, normalization = 0;
  double backward_residual = 0, estimated_forward_sensitivity = 0;
};
struct ProfileResult {
  DensityStatus status = DensityStatus::invalid_input;
  double coefficient = 0, quadratic = 0;
}; // optimized score, NOT a normalized density
enum class MatrixValidationScope {
  full_declared_matrix,
  selected_covariance_only,
  full_precision_then_marginal
};
struct Metadata {
  std::vector<std::string> ordered_ids;
  std::string measure, table_identity, uncertainty_identity,
      ordering_provenance, calibration_provenance, dependence_provenance,
      source_semantics, input_matrix_convention,
      treatment = "normalized Gaussian";
  MatrixValidationScope matrix_validation_scope =
      MatrixValidationScope::full_declared_matrix;
};
// These are original generative priors, not conditional latent posteriors.
struct PriorRecord {
  double mean = 0, variance = 0;
  std::string latent_identity;
  std::vector<double> response;
  std::vector<std::string> applied_row_ids;
  bool independence_assumed = false;
};
struct SelectionRecord {
  std::string operation;
  std::vector<std::string> kept_row_ids, complement_row_ids;
};
enum class MatrixKind { covariance, precision };
class Gaussian {
public:
  DensityStatus status() const noexcept { return status_; }
  const Metadata &metadata() const noexcept { return metadata_; }
  const std::vector<PriorRecord> &priors() const noexcept { return priors_; }
  const std::vector<SelectionRecord> &selection_history() const noexcept {
    return history_;
  }
  std::span<const double> mean_shift() const noexcept { return mean_shift_; }
  std::span<const double> covariance() const noexcept { return covariance_; }
  GaussianResult evaluate(std::span<const double> residual,
                          std::span<const std::string> ordered_ids,
                          double maximum_forward_sensitivity) const;
  Gaussian marginal(std::span<const std::size_t> kept,
                    std::size_t maximum_elements,
                    double maximum_forward_sensitivity) const;
  // Covariance conditional on the complement's CENTERED residual being exactly
  // zero. With proper prior mean shift, base residual on complement equals that
  // shift. Nonzero conditioning needs an explicit conditional mean (not
  // implemented here).
  Gaussian
  conditional_zero_complement(std::span<const std::size_t> kept,
                              std::size_t maximum_elements,
                              double maximum_forward_sensitivity) const;
  ProfileResult profile_offset(std::span<const double> residual,
                               std::span<const double> response,
                               std::span<const std::string> ordered_ids,
                               double maximum_forward_sensitivity) const;
  // Caller explicitly declares D and the proper latent prior as independent.
  // Adds variance*response*response^T exactly once; prior mean explicitly
  // shifts the supplied base residual.
  Gaussian proper_offset(std::span<const double> response,
                         std::span<const std::string> ordered_ids,
                         double prior_mean, double prior_variance,
                         std::string latent_identity,
                         bool independent_prior_declared,
                         std::size_t maximum_elements,
                         double maximum_forward_sensitivity) const;

private:
  DensityStatus status_ = DensityStatus::invalid_input;
  Metadata metadata_;
  std::vector<double> covariance_, mean_shift_;
  std::vector<std::string> latent_ids_;
  std::vector<PriorRecord> priors_;
  std::vector<SelectionRecord> history_;
  numerics::Factorization factor_;
  friend Gaussian prepare_observations(const observations::Prepared &,
                                       observations::Selection, std::size_t,
                                       double);
  friend Gaussian prepare_gaussian(std::span<const double>, MatrixKind,
                                   Metadata, std::size_t, double);
};
Gaussian prepare_gaussian(std::span<const double> matrix, MatrixKind, Metadata,
                          std::size_t maximum_elements,
                          double maximum_forward_sensitivity);
// Selected covariance is a marginal, including when source stores full
// precision. Unknown lineage remains unknown; construction makes no
// independence claim.
Gaussian prepare_observations(const observations::Prepared &,
                              observations::Selection,
                              std::size_t maximum_elements,
                              double maximum_forward_sensitivity);
} // namespace irred::statistics
