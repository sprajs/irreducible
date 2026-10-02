// Same source also compiled against a fresh installed header/archive prefix.
#include "irred/hydrogen_equilibrium.hpp"
#include <array>
#include <cmath>
#include <iostream>
int main() {
  namespace a = irred::atomic;
  std::array states{a::HydrogenState{10000, 6.769876143152199e20}};
  auto result = a::evaluate_hydrogen_equilibrium(states);
  if (result.status != irred::numerics::Status::ok || result.rows.size() != 1 ||
      !result.rows[0].ionized.value || !result.rows[0].neutral.value ||
      std::abs(*result.rows[0].ionized.value - .5) > 1e-12 ||
      std::abs(*result.rows[0].neutral.value - .5) > 1e-12)
    return 1;
  std::cout << "PASS installed hydrogen equilibrium consumer\n";
}
