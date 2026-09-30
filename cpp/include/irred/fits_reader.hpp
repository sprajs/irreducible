#pragma once
#include "irred/observations.hpp"
#include <span>
namespace irred::observations {
enum class DecodeStatus { ok, codec_unavailable, invalid_format, resource_limit };
struct Decoded { DecodeStatus status=DecodeStatus::invalid_format; Input input{}; };
// Narrow synthetic VERIFIER BINTABLE profile only. Decode supplied bytes, never
// re-open a mutable pathname. CFITSIO is format infrastructure, not physics.
Decoded decode_fits_length(std::span<const std::byte>,const std::string& sha256, std::size_t maximum_rows);
}
