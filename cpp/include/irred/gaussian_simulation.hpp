#pragma once
#include "irred/random.hpp"
#include "irred/statistics.hpp"
namespace irred::statistics {
struct GeneratingMean {
  std::vector<double> value;
  std::vector<std::string> ordered_ids, coordinate_units;
  std::string identity, coordinate_measure, generating_law_identity;
};
enum class SimulationOutput : unsigned { values = 1, arithmetic_errors = 2, words = 4 };
struct GaussianSimulationPolicy {
  std::size_t maximum_vectors = 65536, maximum_elements = 1000000;
  std::size_t maximum_payload_bytes = 256 * 1024 * 1024;
  std::size_t maximum_work_units = 100000000;
  double maximum_scaled_arithmetic_error = 1e-10;
  unsigned outputs = 1 | 2;
};
struct GaussianSimulationRow {
  random::Address address;
  numerics::Status status = numerics::Status::invalid_input;
};
struct GaussianSimulationBatch {
  numerics::Status status = numerics::Status::invalid_input;
  std::uint64_t seed = 0;
  std::size_t requested_vectors = 0, dimension = 0, work_units = 0;
  unsigned outputs = 0;
  GeneratingMean mean;
  Metadata source_metadata{"", {}, "", "", "", "", "", "", "", "", ""};
  std::vector<GaussianSimulationRow> rows;
  // Pooled vector-major coordinates. A failed row has no usable numeric slots.
  std::vector<double> values, absolute_error_estimates;
  std::vector<std::array<std::uint32_t, 4>> words;
  const char *method_id() const noexcept { return "Gaussian-generation/retained-colour-halfbin/v1"; }
};
// Borrows an unshifted, prior-free synthetic covariance only during the batch.
// Emitted vectors are a discrete pseudo-Gaussian approximation to N(mean,C).
class GaussianSimulation {
public:
  static std::optional<std::size_t> payload_bound(
      const Gaussian &, const GeneratingMean &, std::size_t count,
      unsigned outputs = 1 | 2) noexcept;
  static std::optional<std::size_t> work_bound(std::size_t dimension,
                                             std::size_t count) noexcept;
  static GaussianSimulationBatch simulate(
      const Gaussian &, const GeneratingMean &,
      std::span<const random::Address>, std::uint64_t seed,
      GaussianSimulationPolicy = {});
};
} // namespace irred::statistics
