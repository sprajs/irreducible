#pragma once
#include "irred/background.hpp"
namespace irred::cosmology {
inline constexpr unsigned growth_d = 1, growth_f = 2;
struct GrowthPolicy {
  double relative_tolerance = 1e-10;
  std::size_t maximum_points = 4096, maximum_native_bytes = 16 * 1024 * 1024;
  std::size_t maximum_callbacks_per_point = 200000,
              maximum_total_callbacks = 4000000;
  unsigned maximum_depth = 40;
};
struct GrowthValue {
  numerics::Status status = numerics::Status::invalid_input;
  std::optional<double> value;
  double absolute_error_estimate = 0;
};
struct GrowthRow {
  double scale_factor = 0;
  GrowthValue d, f;
  size_t callbacks = 0;
};
struct GrowthBatch {
  numerics::Status status = numerics::Status::invalid_input;
  unsigned requested = 0;
  std::vector<GrowthRow> rows;
  size_t callbacks = 0;
};
class GRGrowth {
public:
  GRGrowth() = default;
  GRGrowth(const GRGrowth &) = default;
  GRGrowth &operator=(const GRGrowth &) = default;
  GRGrowth(GRGrowth &&) noexcept;
  GRGrowth &operator=(GRGrowth &&) noexcept;
  numerics::Status status() const noexcept { return status_; }
  const Expansion *background() const noexcept {
    return background_ ? &*background_ : nullptr;
  }
  GrowthBatch evaluate(std::span<const double> scale_factors, unsigned outputs,
                       GrowthPolicy = {}) const;

private:
  std::optional<Expansion> background_;
  numerics::Status status_ = numerics::Status::invalid_input;
  friend GRGrowth prepare_gr_growth(const Expansion &);
};
GRGrowth prepare_gr_growth(const Expansion &);
std::optional<size_t> growth_payload_bound(size_t points) noexcept;
inline constexpr std::string_view gr_growth_id =
    "GR/flat-pressureless-matter-lambda-growing-mode-Heath/v1";
} // namespace irred::cosmology
