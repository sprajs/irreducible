#pragma once
#include "irred/background.hpp"
namespace irred::cosmology::detail {
// Same evaluated kernel, pre-cast projection; no second integration or public
// long-double precision promise. BAO combines these before its final F64 cast.
struct RadialAccess {
  static long double expansion(const ExpansionValue &r) noexcept {
    return r.precise_E_;
  }
  static long double expansion(const Radial &r) noexcept {
    return r.precise_E_;
  }
  static long double integral(const Radial &r) noexcept {
    return r.precise_integral_;
  }
};
} // namespace irred::cosmology::detail
