#pragma once
#include <array>
#include <cstdint>
#include <string_view>
namespace irred::random {
inline constexpr std::string_view generator_id =
    "RNG/Philox4x32-10-halfbin32-BoxMuller/v1";
struct Address { std::uint64_t stream = 0, sample = 0; };
// Injective counter addresses permit replay, not an independence certificate.
std::array<std::uint32_t, 4> words(std::uint64_t seed, Address) noexcept;
long double halfbin(std::uint32_t) noexcept;
// Same words1/2 cosine expression as the original detector. No sine partner.
long double normal_cosine(const std::array<std::uint32_t, 4> &) noexcept;
} // namespace irred::random
