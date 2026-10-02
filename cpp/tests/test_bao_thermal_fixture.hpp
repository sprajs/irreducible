#pragma once
#include "irred/bao_thermal.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
namespace thermal_bao_test {
using namespace irred;
using W = long double;
using S = numerics::Status;
inline unsigned checks = 0;
inline void need(bool ok, const char *s) {
  ++checks;
  if (!ok)
    throw std::runtime_error(s);
}
inline void near(W a, W b, W budget, const char *s) {
  need(std::abs(a - b) <= budget, s);
}
inline cosmology::ThermalObservableRequest massive() {
  return {{67.4,
           .02237,
           .12,
           2.725501088928264,
           1.703989276e-5,
           {{.06, 1.9517599631877276, 2}}},
          1059.95,
          "synthetic supplied drag",
          "matched explicit CLASS constants/FD source"};
}
inline cosmology::ThermalObservableRequest massless() {
  return {{70, .0245, .1225, 2.7255, 1e-5, {}},
          1059.95,
          "synthetic supplied drag",
          "analytic physical density source"};
}
inline bao::ThermalDensityPolicy policy() {
  bao::ThermalDensityPolicy p;
#ifdef IRRED_TEST_NESTED_CC
  p.predictions.thermal.momentum_method=cosmology::ThermalMomentumMethod::nested_clenshaw_curtis;
#endif
  p.maximum_models = 16;
  p.maximum_queries = 16;
  p.maximum_string_bytes = 4096;
  p.maximum_native_bytes = 8 * 1024 * 1024;
  p.maximum_total_callbacks = 200000000;
  p.maximum_forward_sensitivity = 1e-8;
  p.requested = 7;
  return p;
}
inline bao::PreparedDensity prepared(std::span<const W> reference, double z = 1,
                                     bool permute = false) {
  bao::DensityInput d;
  d.queries = {{z, bao::Observable::transverse_over_ruler},
               {z, bao::Observable::hubble_over_ruler},
               {z, bao::Observable::volume_over_ruler}};
  d.observed = {double(reference[0] + .2L), double(reference[1] - .3L),
                double(reference[2] + .1L)};
  d.covariance = {1, .2, .1, .2, 2, .3, .1, .3, 1.5};
  d.ordered_ids = {"DM", "DH", "DV"};
  d.role = bao::RowRole::synthetic_control;
  d.covariance_unit = bao::CovarianceUnit::dimensionless_ratio_squared;
  d.table_identity = "synthetic thermal BAO ratios";
  d.covariance_identity = "synthetic full SPD3";
  d.ordering_provenance = "explicit query order";
  d.calibration_provenance = "conditional supplied drag";
  d.dependence_provenance = "full ordered covariance";
  if (permute) {
    std::swap(d.queries[0], d.queries[2]);
    std::swap(d.observed[0], d.observed[2]);
    std::swap(d.ordered_ids[0], d.ordered_ids[2]);
    auto c = d.covariance;
    constexpr unsigned order[]{2, 1, 0};
    for (unsigned j = 0; j < 3; ++j)
      for (unsigned k = 0; k < 3; ++k)
        d.covariance[j * 3 + k] = c[order[j] * 3 + order[k]];
  }
  return bao::prepare_density(std::move(d),
                              {16, 256, 4096, 8 * 1024 * 1024, 1e-8,
                               numerics::Arithmetic::longdouble_cpu_v1});
}
inline auto evaluate(const bao::PreparedDensity &d,
                     const cosmology::ThermalObservableRequest &r,
                     bao::ThermalDensityPolicy p = policy()) {
  return d.evaluate_thermal(std::span(&r, 1), p);
}
} // namespace thermal_bao_test
