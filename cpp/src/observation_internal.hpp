#pragma once
#include "irred/abi.h"
#include "irred/observations.hpp"
#include <memory>
// Internal borrowed view; ownership remains with the C observation handle.
const irred::observations::Prepared *
native_observations(const irred_prepared *) noexcept;

// Retains the same immutable source after the acquisition handle is destroyed.
std::shared_ptr<const irred::observations::Prepared>
shared_native_observations(const irred_prepared *) noexcept;
