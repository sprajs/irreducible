#pragma once
#include "irred/piecewise_background.hpp"
namespace irred::cosmology::detail {
// Internal precision-preserving projection inputs. Same validation/work as the
// public geometry provider; no public slot/storage or derivative contract
// change.
struct PiecewiseRadialSlot {
  Query source{};
  Status status = Status::invalid_input;
  numerics::Status numerical_status = numerics::Status::invalid_input;
  std::size_t segments_processed = 0;
  long double expansion_E = 0, radial_integral = 0;
};
struct PiecewiseRadialBatch {
  Status status = Status::invalid_input;
  std::vector<PiecewiseRadialSlot> slots;
  std::size_t segments_processed = 0;
};
struct PiecewiseRadialAccess {
  template <bool RadialOnly>
  static auto evaluate(const PiecewiseBackground &, std::span<const Query>,
                       PiecewisePolicy);
  static PiecewiseRadialBatch radial(const PiecewiseBackground &,
                                     std::span<const Query>, PiecewisePolicy);
};
} // namespace irred::cosmology::detail
