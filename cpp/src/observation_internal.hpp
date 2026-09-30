#pragma once
#include "irred/abi.h"
#include "irred/observations.hpp"
// Internal borrowed view; ownership remains with the C observation handle.
const irred::observations::Prepared *
native_observations(const cosmo_prepared *) noexcept;
