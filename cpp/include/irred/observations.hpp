#pragma once
#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <vector>
namespace irred::observations {
enum class Profile : std::uint32_t { pantheon_plus_released_v1, gaussian_fixture_v1, fits_length_fixture_v1, typed_magnitude_covariance };
enum class Role : std::uint32_t { observed_measurement, released_fitted_summary, synthetic_control, posterior_summary };
enum class Unit : std::uint32_t { magnitude, metre };
enum class Calibration : std::uint32_t { unknown, released_corrected, not_applicable };
enum class Uncertainty : std::uint32_t { none, covariance, precision };
enum class UncertaintyUnit : std::uint32_t { none, magnitude_squared, inverse_magnitude_squared, metre_squared, inverse_metre_squared };
enum class Component : std::uint32_t { unknown, statistical, systematic, total };
enum class Selection : std::uint32_t { all, pantheon_zhd_gt_001 };
enum class Status : std::uint32_t { ok, invalid_shape, invalid_identity, incompatible_semantics, missing_required_value, nonfinite_required_value, resource_limit };
// Source missingness is separate from a computed failure. A masked value has no
// scientific payload; its encoding remains reconstructible in the source asset.
struct Input {
 Profile profile; Role role; Unit unit; Calibration calibration; Uncertainty uncertainty;
 UncertaintyUnit uncertainty_unit=UncertaintyUnit::none; Component component=Component::unknown;
 std::string ordering_provenance; // supplied release-order declaration, not independently inferred from matrix bytes
 std::string calibration_provenance, dependence_provenance; // empty means unknown
 std::string quality_dictionary; // raw flags retained even when dictionary unknown
 std::string table_sha256, uncertainty_sha256;
 std::vector<std::string> measurement_ids, event_ids, uncertainty_axis_ids;
 std::vector<double> values, zhd, zcmb, zhel, uncertainty_matrix;
 std::vector<std::uint8_t> missing, zhd_missing, zcmb_missing, zhel_missing;
 std::vector<std::uint64_t> quality;
 std::vector<std::uint8_t> source_selection;
};
struct Policy { std::size_t maximum_rows; std::size_t maximum_matrix_elements; std::size_t maximum_string_bytes=1048576; };
struct SelectionResult { Status status=Status::invalid_shape; std::vector<std::uint8_t> mask; std::vector<std::size_t> source_indices; };
class Prepared {
public:
 Status status() const noexcept { return status_; }
 const Input& source() const noexcept { return source_; }
 SelectionResult select(Selection) const;
 // Unknown calibration/dependence is retained, never an independence claim.
 bool calibration_state_declared() const noexcept { return source_.calibration!=Calibration::unknown; }
bool calibration_provenance_known() const noexcept { return !source_.calibration_provenance.empty(); }
private:
 Status status_=Status::invalid_shape;
 Input source_{};
 friend Prepared prepare(Input,Policy);
};
// Copies/moves into owned immutable storage. No covariance arithmetic, repair,
// SPD claim, deduplication, or reinterpretation of posterior products.
Prepared prepare(Input,Policy);
}
