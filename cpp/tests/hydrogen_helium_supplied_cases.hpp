#pragma once
#include "irred/hydrogen_helium_history.hpp"
#include <array>
#include <vector>
namespace hydrogen_helium_supplied_cases {
using namespace irred::cosmology;
inline HydrogenHeliumSuppliedHistoryRequest source(unsigned id = 0) {
  HydrogenHeliumSuppliedHistoryRequest request{
      {{67.4, .02237, .12, 2.7255, 1.7e-5, {}}, .19, .015,
       "synthetic independent nuclei NEXT15-controls/v1#B0", 2700, 300},
      {.9990234375, .5, 7000, "synthetic explicit nonLTE boundary",
       "NEXT15-controls/v1#B0"}};
  if (id == 1) request = {
      {{60, .015, .08, 2.7, 0, {}}, .1, .005,
       "synthetic independent nuclei NEXT15-controls/v1#B1", 2600, 600},
      {.998046875, .25, 6000, "synthetic explicit nonLTE boundary", "NEXT15-controls/v1#B1"}};
  if (id == 2) request = {
      {{80, .03, .15, 2.75, 3e-5, {{0, 1.95, 2}}}, .3, .035,
       "synthetic independent nuclei NEXT15-controls/v1#B2", 2800, 300},
      {.99951171875, .75, 7600, "synthetic explicit nonLTE boundary", "NEXT15-controls/v1#B2"}};
  return request;
}
inline std::vector<double> queries(unsigned id, bool layer = true) {
  const auto request = source(id);
  const double zi = request.history.initial_redshift;
  std::vector<double> z{zi, zi - .01, 2500, 2200, 1900, 1600, 1300, 1100, 800};
  if (id != 1) z.push_back(600);
  z.push_back(request.history.late_redshift);
  if (layer) for (const double offset : std::array<double, 5>{0x1p-20, 0x1p-12, 0x1p-8, 0x1p-4, .5})
    z.push_back(zi - offset);
  return z;
}
} // namespace hydrogen_helium_supplied_cases
