#pragma once
#include "irred/statistics.hpp"
namespace irred::statistics {
// This is an explicitly relative quadratic profile, never a density/evidence.
enum class DesignRank {
  unassessed,
  deficient,
  unresolved,
  full_within_conditioning_contract
};
struct DesignMetadata {
  std::vector<std::string> ordered_parameter_ids, parameter_units,
      shared_nuisance_ids;
  std::string residual_unit, design_identity, dependence_identity;
};
struct DesignPolicy {
  std::size_t maximum_elements = 1000000;
  std::size_t maximum_payload_bytes = 256 * 1024 * 1024;
  double maximum_forward_sensitivity = 1e-10;
};
struct DesignResult {
  DensityStatus status = DensityStatus::invalid_input;
  numerics::Status numerical_status = numerics::Status::invalid_input;
  std::vector<double> coefficients, adjusted_residuals;
  double quadratic = 0, relative_log_score = 0;
  double normalized_normal_equation_residual = 0;
  double covariance_solve_backward_residual = 0,
         covariance_solve_forward_sensitivity = 0;
};
class DesignProfile {
public:
  DesignProfile() = default;
  DesignProfile(const DesignProfile &) = delete;
  DesignProfile &operator=(const DesignProfile &) = delete;
  DesignProfile(DesignProfile &&) noexcept;
  DesignProfile &operator=(DesignProfile &&) noexcept;
  // X is row-major; p is the number of explicitly ordered parameters (>=2).
  // Failure leaves source usable; successful preparation consumes it.
  static DesignProfile prepare(Gaussian &&source, std::span<const double> x,
                               std::span<const std::string> ordered_row_ids,
                               DesignMetadata metadata,
                               DesignPolicy policy = {});
  static std::optional<std::size_t>
  preparation_payload_bound(const Gaussian &, std::size_t p,
                            const DesignMetadata &) noexcept;
  std::optional<std::size_t> retained_payload_bound() const noexcept;
  std::optional<std::size_t> evaluation_payload_bound() const noexcept;
  DensityStatus status() const noexcept { return status_; }
  numerics::Status numerical_status() const noexcept {
    return numerical_status_;
  }
  DesignRank rank() const noexcept { return rank_; }
  double equilibrated_gram_condition_inf() const noexcept {
    return gram_condition_;
  }
  const Metadata &metadata() const noexcept { return gaussian_.metadata(); }
  const DesignMetadata &design_metadata() const noexcept {
    return design_metadata_;
  }
  const std::vector<SelectionRecord> &selection_history() const noexcept {
    return gaussian_.selection_history();
  }
  numerics::Arithmetic arithmetic() const noexcept {
    return gaussian_.arithmetic();
  }
  DesignResult evaluate(std::span<const double> residual,
                        std::span<const std::string> ordered_row_ids,
                        DesignPolicy policy = {}) const;

private:
  Gaussian gaussian_;
  DesignMetadata design_metadata_;
  numerics::Factorization gram_factor_;
  std::vector<double> x_, wx_, scales_;
  DensityStatus status_ = DensityStatus::invalid_input;
  numerics::Status numerical_status_ = numerics::Status::invalid_input;
  DesignRank rank_ = DesignRank::unassessed;
  double gram_condition_ = 0, preparation_sensitivity_ = 0;
};
} // namespace irred::statistics
