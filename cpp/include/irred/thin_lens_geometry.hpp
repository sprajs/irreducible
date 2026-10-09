#pragma once
#include "irred/curved_flrw.hpp"
#include "irred/thermal_observables.hpp"
#include <array>
#include <variant>
namespace irred::cosmology {
struct LensEpochPair { double lens_redshift=0, source_redshift=0; };
enum class LensGeometryOutput : unsigned {
  lens_angular_mpc, source_angular_mpc, lens_source_angular_mpc,
  distance_ratio, critical_surface_density_kg_m2, time_delay_distance_mpc, count
};
inline constexpr unsigned lens_geometry_output_count=
    static_cast<unsigned>(LensGeometryOutput::count);
inline constexpr unsigned lens_geometry_mask(LensGeometryOutput output) {
  const auto n=static_cast<unsigned>(output);
  return n<lens_geometry_output_count?1u<<n:0;
}
struct LensGeometryValue {
  numerics::Status status=numerics::Status::invalid_input;
  std::optional<double> value;
  double error_estimate=0;
};
struct LensGeometryPolicy {
  double absolute_tolerance_mpc=1e-8, relative_tolerance=1e-8;
  double absolute_tolerance_ratio=1e-12;
  double absolute_tolerance_kg_m2=1e-10;
  std::size_t maximum_pairs=2048, maximum_native_bytes=16*1024*1024;
  std::size_t maximum_callbacks_per_pair=20000000, maximum_total_callbacks=100000000;
  unsigned maximum_depth=40;
  // Applies only to the thermal provider, including its retained method check.
  ThermalPolicy thermal=[] { auto p=ThermalPolicy{};
    p.maximum_total_callbacks=100000000;return p; }();
};
struct LensGeometryRow {
  LensEpochPair epochs;
  numerics::Status status=numerics::Status::invalid_input;
  std::array<LensGeometryValue,lens_geometry_output_count> outputs{};
  std::size_t callbacks=0, outer_callbacks=0, momentum_callbacks=0;
};
struct LensGeometryBatch {
  numerics::Status status=numerics::Status::invalid_input;
  unsigned requested_outputs=0;
  std::vector<LensGeometryRow> rows;
  std::size_t callbacks=0, outer_callbacks=0, momentum_callbacks=0;
};
// Owns one actual immutable provider, with no interchangeable input distances.
// Flat FD and nonflat dust/radiation/Lambda retain different source identities.
class ThinLensGeometry {
public:
  ThinLensGeometry()=default;
  explicit ThinLensGeometry(const ThermalObservables &, LensGeometryPolicy={});
  explicit ThinLensGeometry(const CurvedFLRW &, LensGeometryPolicy={});
  ThinLensGeometry(const ThinLensGeometry &)=default;
  ThinLensGeometry &operator=(const ThinLensGeometry &)=default;
  ThinLensGeometry(ThinLensGeometry &&) noexcept;
  ThinLensGeometry &operator=(ThinLensGeometry &&) noexcept;
  numerics::Status status() const noexcept { return status_; }
  std::string_view provider_id() const noexcept;
  const ThermalObservables *thermal_provider() const noexcept;
  const CurvedFLRW *curved_provider() const noexcept;
  std::optional<std::size_t> retained_payload_bound() const noexcept;
  // Includes owner, rows and conservative direct-provider temporary payload;
  // excludes borrowed inputs, allocator metadata, quadrature stack and RSS.
  std::optional<std::size_t> evaluation_payload_bound(std::size_t pairs) const noexcept;
  LensGeometryBatch evaluate(std::span<const LensEpochPair>, unsigned outputs,
                             LensGeometryPolicy={}) const;
private:
  std::variant<std::monostate,ThermalObservables,CurvedFLRW> provider_;
  numerics::Status status_=numerics::Status::invalid_input;
};
inline constexpr std::string_view thin_lens_geometry_id=
    "same-observer-direct-two-epoch-thin-lens-physical-SI-prefactors/v1";
} // namespace irred::cosmology
