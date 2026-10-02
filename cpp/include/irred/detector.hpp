#pragma once
#include "irred/numerics.hpp"
#include "irred/random.hpp"
#include <array>
#include <cstdint>
#include <optional>
#include <span>
#include <string_view>
#include <vector>
namespace irred::detector {
inline constexpr std::string_view model_id =
    "DETECTOR/Poisson-arrivals-fixed-QE-Gaussian-read/v1";
inline constexpr std::string_view generator_id = random::generator_id;
enum class PhotonLaw { unspecified, poisson_arrivals };
struct Input {
  PhotonLaw photon_law = PhotonLaw::unspecified;
  double expected_transmitted_photons = 0, quantum_efficiency = 0;
  double background_expected_electrons = 0, dark_electrons_per_second = 0;
  double observer_exposure_second = 0, read_noise_rms_electrons = 0;
  double gain_electrons_per_adu = 1, bias_adu = 0;
};
struct Moments {
  numerics::Status status = numerics::Status::invalid_input;
  std::optional<double> signal_electrons, poisson_mean_electrons, mean_adu,
      variance_adu;
};
Moments moments(const Input &) noexcept;
using Address = random::Address;
// Counter mapping is injective in (stream,sample) for each caller-owned seed.
// Distinct counter addresses are not a mathematical independence certificate.
std::array<std::uint32_t, 4> random_words(std::uint64_t seed, Address) noexcept;
struct Request {
  Input source;
  Address address;
};
struct Draw {
  Input source;
  Address address;
  numerics::Status status = numerics::Status::invalid_input;
  std::array<std::uint32_t, 4> words{};
  std::optional<std::uint32_t> poisson_electrons;
  std::optional<double> read_electrons, measured_adu;
};
struct Policy {
  std::size_t maximum_rows = 0, maximum_payload_bytes = 0;
};
struct Batch {
  numerics::Status status = numerics::Status::invalid_input;
  std::uint64_t seed = 0;
  std::vector<Draw> rows;
};
std::optional<std::size_t> output_payload_bound(std::size_t) noexcept;
// No mutable RNG state. Duplicate addresses within one batch fail admission.
Batch simulate(std::span<const Request>, std::uint64_t seed, Policy);
} // namespace irred::detector
