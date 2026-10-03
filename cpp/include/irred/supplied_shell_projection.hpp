#pragma once
#include "irred/numerics.hpp"
#include <array>
#include <atomic>
#include <memory>
#include <string_view>

namespace irred::projection {
inline constexpr unsigned supplied_shell_signed_projection = 1;
inline constexpr std::string_view supplied_shell_model_id =
    "synthetic/supplied-real-unit-zeta-isotropic-finite-radial-atoms/v1";
inline constexpr std::string_view supplied_shell_method_id =
    "compact-phase-direct-spherical-Bessel-series-alternating-tail/v1";
inline constexpr std::string_view supplied_shell_arithmetic_id =
    "strict-wide-signed-accumulation-separated-argument-term-product-cast-"
    "diagnostics/v1";

enum class SourceMode { scalar_unit_zeta };
enum class GeometryUnits { comoving_mpc_no_h };
enum class FourierConvention { outward_exp_plus_i };
enum class AtomRole { radial_atom, boundary_response };
enum class ProjectionCheck { unavailable, numerical_radius_within_allocation,
                             refused };
enum class Refusal {
  none, source_shape, source_tags, source_strings, source_number,
  resource_policy, allocation, arithmetic_profile, request, phase,
  nonnormal_arithmetic, series_tail, numerical_allocation
};
struct Failure {
  Refusal reason = Refusal::none;
  std::size_t k_index = static_cast<std::size_t>(-1);
  std::size_t shell_index = static_cast<std::size_t>(-1);
  std::size_t ell_index = static_cast<std::size_t>(-1);
};
struct EndpointView {
  std::string_view producer_identity, coordinate_identity, unit_identity;
  double support_lower = 0, support_upper = 0;
  numerics::Status survival_status = numerics::Status::invalid_input;
  std::optional<double> survival, survival_absolute_diagnostic;
};
struct SuppliedShellView {
  std::span<const double> k_mpc_inverse, chi_mpc, amplitudes;
  std::span<const std::string_view> k_ids, shell_ids;
  std::span<const AtomRole> roles;
  std::string_view source_identity, mode_origin, ordering_provenance;
  SourceMode mode = SourceMode::scalar_unit_zeta;
  GeometryUnits units = GeometryUnits::comoving_mpc_no_h;
  FourierConvention convention = FourierConvention::outward_exp_plus_i;
  std::optional<EndpointView> endpoint;
};
struct ProjectionResources {
  std::size_t maximum_total_work = 8000000;
  std::size_t maximum_series_terms = 6000000;
  std::size_t maximum_source_bytes = 65536;
  std::size_t maximum_live_bytes = 4 * 1024 * 1024;
};
struct ProjectionWork {
  std::size_t source_scalar_inspections = 0, source_tag_inspections = 0;
  std::size_t shell_visits = 0, phase_products = 0, Bessel_evaluations = 0;
  std::size_t series_terms = 0, shell_products = 0, signed_additions = 0;
  std::size_t radius_operations = 0, output_casts = 0, total = 0;
  std::size_t source_bytes_inspected = 0, rejected_work_requests = 0;
};
struct ProjectionRequest {
  std::span<const unsigned> ell;
  unsigned outputs = supplied_shell_signed_projection;
  // Tightens the exact rational (1+abs(emitted value))/10^10 by 2^-q.
  unsigned accuracy_reduction_power = 0;
};
struct ProjectionDiagnostics {
  double argument = 0, Bessel_tail = 0, Bessel_arithmetic = 0;
  double shell_products = 0, signed_accumulation = 0;
  double value_cast_loss = 0, radius_assembly_padding = 0;
};
struct ProjectionRow {
  std::size_t k_index = 0, ell_index = 0;
  unsigned ell = 0;
  numerics::Status status = numerics::Status::invalid_input;
  ProjectionCheck check = ProjectionCheck::unavailable;
  std::optional<double> signed_value, absolute_numerical_radius;
  std::optional<ProjectionDiagnostics> diagnostics;
};
namespace detail {
struct IdentifierOffset { std::size_t offset = 0, length = 0; };
struct SourceStorage;
struct SourceMaker;
struct SourceAccess;
struct PreparedState;
struct ResultState;
template <class T> struct TypedDelete {
  std::size_t count = 0;
  void operator()(T *) const noexcept;
};
template <class T> using TypedArray = std::unique_ptr<T[], TypedDelete<T>>;
} // namespace detail
class SuppliedShellSource {
public:
  ~SuppliedShellSource();
  SuppliedShellSource(const SuppliedShellSource &) = delete;
  SuppliedShellSource &operator=(const SuppliedShellSource &) = delete;
  SuppliedShellSource(SuppliedShellSource &&) = delete;
  SuppliedShellSource &operator=(SuppliedShellSource &&) = delete;
  std::span<const double> k_mpc_inverse() const noexcept;
  std::span<const double> chi_mpc() const noexcept;
  std::span<const double> amplitudes() const noexcept;
  std::span<const AtomRole> roles() const noexcept;
  std::string_view k_id(std::size_t) const noexcept;
  std::string_view shell_id(std::size_t) const noexcept;
  std::string_view source_identity() const noexcept;
  std::string_view mode_origin() const noexcept;
  std::string_view ordering_provenance() const noexcept;
  std::optional<EndpointView> endpoint() const noexcept;
  bool has_boundary_response() const noexcept;
  std::size_t source_payload_bytes() const noexcept;
  // Charged requested storage; allocator metadata, caller handles and RSS excluded.
  std::size_t live_payload_bytes() const noexcept;
  std::size_t peak_payload_bytes() const noexcept;

private:
  explicit SuppliedShellSource(detail::SourceStorage &&);
  std::string_view text(detail::IdentifierOffset) const noexcept;
  bool reserve(std::size_t) const noexcept;
  void release(std::size_t) const noexcept;
  std::size_t k_count_ = 0, shell_count_ = 0, byte_count_ = 0;
  std::size_t payload_ = 0, maximum_live_ = 0;
  detail::TypedArray<double> k_, chi_, amplitudes_;
  detail::TypedArray<detail::IdentifierOffset> ids_;
  detail::TypedArray<AtomRole> roles_;
  detail::TypedArray<char> bytes_;
  std::array<detail::IdentifierOffset, 6> labels_{};
  std::optional<EndpointView> endpoint_;
  mutable std::atomic<std::size_t> live_{0}, peak_{0};
  friend struct detail::SourceMaker;
  friend struct detail::SourceAccess;
};
class ProjectionBatch {
public:
  ProjectionBatch();
  ~ProjectionBatch();
  ProjectionBatch(const ProjectionBatch &) = delete;
  ProjectionBatch &operator=(const ProjectionBatch &) = delete;
  ProjectionBatch(ProjectionBatch &&) noexcept;
  ProjectionBatch &operator=(ProjectionBatch &&) noexcept;
  numerics::Status status() const noexcept;
  const SuppliedShellSource *source() const noexcept;
  std::span<const ProjectionRow> rows() const noexcept;
  std::span<const unsigned> requested_ell() const noexcept;
  const ProjectionWork &work_before() const noexcept;
  const ProjectionWork &work_after() const noexcept;
  const ProjectionWork &work_delta() const noexcept;
  Failure failure() const noexcept;
  std::size_t retained_payload_bytes() const noexcept;
  std::size_t peak_payload_bytes() const noexcept;

private:
  void reset() noexcept;
  std::shared_ptr<const SuppliedShellSource> source_;
  std::unique_ptr<detail::ResultState> state_;
  numerics::Status fallback_ = numerics::Status::invalid_input;
  ProjectionWork before_{}, after_{}, delta_{};
  Failure failure_{};
  std::size_t charge_ = 0, peak_ = 0;
  friend class PreparedProjection;
};
class PreparedProjection {
public:
  PreparedProjection();
  ~PreparedProjection();
  PreparedProjection(const PreparedProjection &) = delete;
  PreparedProjection &operator=(const PreparedProjection &) = delete;
  PreparedProjection(PreparedProjection &&) noexcept;
  PreparedProjection &operator=(PreparedProjection &&) noexcept;
  numerics::Status status() const noexcept;
  const SuppliedShellSource *source() const noexcept;
  const ProjectionWork &cumulative_work() const noexcept;
  Failure failure() const noexcept;
  std::size_t peak_payload_bytes() const noexcept;
  ProjectionBatch evaluate(ProjectionRequest);

private:
  void reset() noexcept;
  std::shared_ptr<const SuppliedShellSource> source_;
  std::unique_ptr<detail::PreparedState> state_;
  numerics::Status fallback_ = numerics::Status::invalid_input;
  ProjectionWork fallback_work_{};
  Failure failure_{};
  std::size_t charge_ = 0;
  std::size_t fallback_peak_ = 0;
  friend PreparedProjection acquire_supplied_shell_projection(
      SuppliedShellView, ProjectionResources);
};
PreparedProjection acquire_supplied_shell_projection(
    SuppliedShellView, ProjectionResources = {});
bool supplied_shell_arithmetic_profile() noexcept;
} // namespace irred::projection
