#include "irred/random.hpp"
#include <cmath>
#include <numbers>
namespace irred::random {
std::array<std::uint32_t, 4> words(std::uint64_t seed, Address a) noexcept {
  // Original transcription of the published Philox4x32 round, not library code.
  std::array<std::uint32_t, 4> c{
      std::uint32_t(a.sample), std::uint32_t(a.sample >> 32),
      std::uint32_t(a.stream), std::uint32_t(a.stream >> 32)};
  std::uint32_t k0 = seed, k1 = seed >> 32;
  for (unsigned r = 0; r < 10; ++r) {
    const std::uint64_t p0 = std::uint64_t(0xD2511F53u) * c[0],
                        p1 = std::uint64_t(0xCD9E8D57u) * c[2];
    c = {std::uint32_t(p1 >> 32) ^ c[1] ^ k0, std::uint32_t(p1),
         std::uint32_t(p0 >> 32) ^ c[3] ^ k1, std::uint32_t(p0)};
    k0 += 0x9E3779B9u;
    k1 += 0xBB67AE85u;
  }
  return c;
}
long double halfbin(std::uint32_t x) noexcept {
  return (static_cast<long double>(x) + .5L) / 4294967296.L;
}
long double normal_cosine(const std::array<std::uint32_t, 4> &w) noexcept {
  return std::sqrt(-2 * std::log(halfbin(w[1]))) *
         std::cos(2 * std::numbers::pi_v<long double> * halfbin(w[2]));
}
} // namespace irred::random
