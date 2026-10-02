#pragma once
namespace irred::cosmology::detail {
struct ThermalCCSample { double rho, pressure; };
// One actual shared node; callers charge one/two scalar integrands before entry.
// Separate translation unit permits forwarding-only linker observation in tests.
ThermalCCSample thermal_cc_sample(double q, long double y, long double scale,
                                 bool need_pressure) noexcept;
} // namespace irred::cosmology::detail
