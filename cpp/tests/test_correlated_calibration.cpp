#include "irred/correlated_calibration.hpp"
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>
#include <cmath>
#include <limits>
#include <numbers>
using namespace irred::statistics;
Metadata metadata() {
  Metadata m;
  m.ordered_ids = {"r0", "r1"};
  m.measure = "product d(magnitude)";
  m.table_identity = "synthetic";
  m.ordering_provenance = "explicit synthetic row order";
  return m;
}
CalibrationPrior prior() {
  CalibrationPrior p;
  p.ordered_parameter_ids = {"zero", "colour"};
  p.parameter_units = {"magnitude", "magnitude"};
  p.mean = {.5, -.25};
  p.covariance = {1, .375, .375, .5};
  p.prior_identity = "synthetic proper prior";
  p.response_identity = "fixed X";
  p.residual_unit = "magnitude";
  p.calibration_identity = "synthetic calibration";
  p.dependence_identity = "noise independent of calibration";
  p.measure_identity = "product d(magnitude)";
  p.noise_independence_declared = true;
  return p;
}
Gaussian base() {
  return prepare_gaussian(std::vector<double>{2, .25, .25, 1.5},
                          MatrixKind::covariance, metadata(), 100, 1e-10);
}
void close(double actual, long double expected) {
  assert(std::abs(static_cast<long double>(actual) - expected) <=
         2e-12L * (1 + std::abs(expected)));
}
int main() {
  const std::vector<double> x = {1, .5, -.25, 1}, r = {1.25, -.5};
  const auto ids = metadata().ordered_ids;
  auto g = base();
  auto p = prior();
  auto c = CorrelatedCalibration::prepare(std::move(g), x, ids, p);
  assert(c.status() == DensityStatus::finite &&
         g.status() == DensityStatus::invalid_input);
  assert(c.source().metadata().table_identity == "synthetic" &&
         c.source().priors().empty());
  // Independent explicit 2x2 cofactor arithmetic, no production factor oracle.
  long double a = 2 + 1 + .375L + .125L,
              b = .25L - .25L + .375L - .046875L + .25L,
              d = 1.5L + .0625L - .1875L + .5L;
  long double e0 = r[0] - (.5L - .125L), e1 = r[1] - (-.125L - .25L),
              det = a * d - b * b;
  long double q = (d * e0 * e0 - 2 * b * e0 * e1 + a * e1 * e1) / det;
  auto v = c.evaluate(r, ids);
  assert(v.density.status == DensityStatus::finite);
  close(v.quadratic, q);
  close(v.log_determinant, std::log(det));
  close(v.density.log_value,
        -.5L * (q + std::log(det) +
                2 * std::log(2 * std::numbers::pi_v<long double>)));
  auto moved = std::move(c);
  assert(c.status() == DensityStatus::invalid_input);
  auto &self = moved;
  moved = std::move(self);
  assert(moved.status() == DensityStatus::finite);
  auto zbase = base();
  auto z = CorrelatedCalibration::prepare(
      std::move(zbase), std::vector<double>(4, 0), ids, prior());
  assert(z.status() == DensityStatus::finite);
  close(z.evaluate(r, ids).density.log_value,
        z.source().evaluate(r, ids, 1e-10).density.log_value);
  // Diagonal prior agrees with sequential proper_offset; shared production
  // ancestry.
  auto diag = prior();
  diag.covariance = {1, 0, 0, .5};
  auto dbase = base();
  auto dc = CorrelatedCalibration::prepare(std::move(dbase), x, ids, diag);
  auto scalar = base()
                    .proper_offset(std::vector<double>{1, -.25}, ids, .5, 1,
                                   "a", true, 100, 1e-10)
                    .proper_offset(std::vector<double>{.5, 1}, ids, -.25, .5,
                                   "b", true, 100, 1e-10);
  close(dc.evaluate(r, ids).density.log_value,
        scalar.evaluate(r, ids, 1e-10).density.log_value);
  assert(std::abs(dc.evaluate(r, ids).density.log_value - v.density.log_value) >
         1e-4);
  auto rankbase = base();
  auto rank = CorrelatedCalibration::prepare(
      std::move(rankbase), std::vector<double>{1, 1, 2, 2}, ids, prior());
  assert(rank.status() == DensityStatus::finite);
  for (int mode = 0; mode < 6; ++mode) {
    auto bad = prior();
    if (mode == 0)
      bad.covariance = {1, 1, 1, 1};
    if (mode == 1)
      bad.mean[0] = std::numeric_limits<double>::infinity();
    if (mode == 2)
      bad.noise_independence_declared = false;
    if (mode == 3)
      bad.ordered_parameter_ids[1] = "zero";
    if (mode == 4)
      bad.measure_identity = "other";
    if (mode == 5)
      bad.covariance[0] = std::numeric_limits<double>::denorm_min();
    auto src = base();
    auto f = CorrelatedCalibration::prepare(std::move(src), x, ids, bad);
    assert(f.status() != DensityStatus::finite &&
           src.status() == DensityStatus::finite);
  }
  auto bsrc = base();
  CalibrationPolicy small;
  small.maximum_payload_bytes = 1;
  assert(CorrelatedCalibration::prepare(std::move(bsrc), x, ids, prior(), small)
                 .status() != DensityStatus::finite &&
         bsrc.status() == DensityStatus::finite);
  CalibrationPolicy impossible;
  impossible.maximum_forward_sensitivity = 1e-25;
  assert(moved.evaluate(r, ids, impossible).density.status !=
         DensityStatus::finite);
  assert(
      moved.evaluate(r, std::vector<std::string>{"r1", "r0"}).density.status ==
      DensityStatus::incompatible_metadata);
}
