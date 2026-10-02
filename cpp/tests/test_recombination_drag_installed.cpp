// Compiled separately against a fresh installed header/archive prefix.
#include "irred/recombination_drag.hpp"
#include <array>
#include <cmath>
#include <iostream>
int main() {
  using namespace irred::cosmology;
  PureHydrogenRequest request{{67.4, .02237, .12, 2.7255, 1.7e-5, {}}};
  auto history = prepare_pure_hydrogen_history(request);
  std::array<double, 3> z{1200, 1000, 300};
  auto rows =
      history.evaluate(z, hydrogen_electron_fraction | hydrogen_drag_depth);
  auto root = history.conditional_unit_depth_redshift();
  if (history.status() != irred::numerics::Status::ok ||
      rows.rows.size() != 3 || !rows.rows[0].electron_fraction.value ||
      !rows.rows[1].drag_depth.value || !rows.rows[2].drag_depth.value ||
      !root.value ||
      std::abs(*rows.rows[0].electron_fraction.value - .2908408191239721) >
          3.01e-7 ||
      std::abs(*rows.rows[1].drag_depth.value - .4729377156892689) > 6.73e-7 ||
      *rows.rows[2].drag_depth.value != 0 ||
      std::abs(*root.value - 1053.4095801604812) > .002)
    return 1;
  std::cout << "PASS installed conditional pure-H history consumer\n";
}
