#pragma once
#include <cstdint>
#include <span>
#include <string_view>
namespace irred {
inline constexpr std::string_view constant_set_id = "SI-IAU-definitions-v1";
inline constexpr double speed_of_light_m_per_s = 299792458.0;
enum class Dimension : std::uint32_t { dimensionless, length, time, inverse_time };
template<Dimension D> struct Scalar { double value; };
using Length = Scalar<Dimension::length>;
using Duration = Scalar<Dimension::time>;
using ExpansionRate = Scalar<Dimension::inverse_time>;
enum class Unit : std::uint32_t { one=1, metre, kilometre, second, day, inverse_second, parsec, megaparsec, km_per_s_per_mpc };
enum class Role : std::uint32_t { ratio=1, redshift, physical_length, comoving_distance, luminosity_distance, angular_diameter_distance, duration, expansion_rate };
enum class Frame : std::uint32_t { none=0, heliocentric, cmb, model };
enum class LengthConvention : std::uint32_t { none=0, physical, comoving_a0_one };
enum class QuantityStatus : std::uint32_t { ok=0, unknown_unit, dimension_mismatch, role_mismatch, missing_convention, unsupported_transform, nonfinite_input, overflow, underflow, constant_set_mismatch, invalid_domain, invalid_batch };
struct Quantity { double value; Unit unit; Role role; Frame frame=Frame::none; LengthConvention convention=LengthConvention::none; std::string_view constants=constant_set_id; };
struct Target { Unit unit; Role role; Frame frame=Frame::none; LengthConvention convention=LengthConvention::none; std::string_view constants=constant_set_id; };
// Metadata string views borrow caller-owned storage; no retained prepared objects.
struct Conversion { QuantityStatus status=QuantityStatus::invalid_batch; Quantity source{}; Quantity target{}; };
Conversion convert(Quantity source, Target target) noexcept;
QuantityStatus convert_batch(std::span<const Quantity> source, Target target, std::span<Conversion> output) noexcept;
double parsec_in_metres() noexcept;
} // namespace irred
