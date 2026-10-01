// Synthetic linear controls only. Frozen coefficient/residual allocation:
// 2e-12 absolute + 2e-12 relative; quadratic 1e-10 absolute, with 1e-12
// independent small exact controls. Rational references derived by KKT below.
#include "irred/gaussian_design.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
using namespace irred::statistics;
namespace {
int checks = 0;
void check(bool b, const char *s) {
  ++checks;
  if (!b)
    throw std::runtime_error(s);
}
void near(double a, double b, const char *s) {
  check(std::abs(a - b) <= 2e-12 + 2e-12 * std::abs(b), s);
}
Metadata rows(std::size_t n) {
  Metadata m;
  for (std::size_t i = 0; i < n; ++i)
    m.ordered_ids.push_back("row" + std::to_string(i));
  m.measure = "product d(magnitude)";
  m.ordering_provenance = "synthetic explicit order";
  m.calibration_provenance = "synthetic calibration";
  m.dependence_provenance = "known synthetic covariance";
  return m;
}
DesignMetadata columns() {
  return {{"anchor", "host", "zero-point"},
          {"mag", "mag", "mag"},
          {"zero-point"},
          "mag",
          "synthetic two-anchor/host/shared-zero-point v1",
          "single shared nuisance, covariance declared"};
}
Gaussian gaussian(std::size_t n) {
  std::vector<double> c(n * n);
  for (std::size_t i = 0; i < n; ++i)
    c[i * n + i] = 1;
  return prepare_gaussian(c, MatrixKind::covariance, rows(n), n * n, 1e-10);
}
} // namespace
int main() {
  try {
    // Six observations: two anchors, two host pairs, one shared zero point.
    // X columns are disjoint anchor contrast, host contrast, shared zero point.
    // Orthogonal controls r=X*(2,-3,.5)+e with X^T e=0, e=(1,1,-1,-1,0,0).
    const std::vector<double> x = {1, 0,  1, -1, 0, 1, 0, 1,  1,
                                   0, -1, 1, 0,  1, 1, 0, -1, 1};
    const std::vector<double> r = {3.5, -.5, -3.5, 2.5, -2.5, 3.5};
    auto g = gaussian(6);
    auto ids = g.metadata().ordered_ids;
    auto profile = DesignProfile::prepare(std::move(g), x, ids, columns());
    check(profile.status() == DensityStatus::finite,
          "prepare named anchor host shared control");
    check(g.status() == DensityStatus::invalid_input,
          "consume Gaussian source");
    check(profile.rank() == DesignRank::full_within_conditioning_contract,
          "qualified full rank");
    auto v = profile.evaluate(r, ids);
    check(v.status == DensityStatus::finite, "evaluate exact control");
    near(v.coefficients[0], 2, "anchor coefficient");
    near(v.coefficients[1], -3, "host coefficient");
    near(v.coefficients[2], .5, "shared calibration coefficient");
    const std::vector<double> e = {1, 1, -1, -1, 0, 0};
    for (std::size_t i = 0; i < 6; ++i)
      near(v.adjusted_residuals[i], e[i], "adjusted exact residual");
    check(std::abs(v.quadratic - 4) <= 1e-12, "independent quadratic");
    near(v.relative_log_score, -2, "relative score");
    check(v.normalized_normal_equation_residual <= 1e-10,
          "actual output normal defect");
    auto moved = std::move(profile);
    check(profile.status() == DensityStatus::invalid_input,
          "operator move invalidation");
    check(moved.evaluate(r, ids).status == DensityStatus::finite,
          "moved operator works");
    moved = std::move(moved);
    check(moved.status() == DensityStatus::finite, "self move retains");
    auto wrong = ids;
    std::swap(wrong[0], wrong[1]);
    check(moved.evaluate(r, wrong).status ==
              DensityStatus::incompatible_metadata,
          "row order rejects");
    auto nan = r;
    nan[0] = NAN;
    check(moved.evaluate(nan, ids).status == DensityStatus::invalid_input,
          "finite residual contract");
    DesignPolicy tiny;
    tiny.maximum_payload_bytes = 1;
    auto rejected = moved.evaluate(r, ids, tiny);
    check(rejected.status == DensityStatus::numerical_failure &&
              rejected.coefficients.empty() &&
              rejected.adjusted_residuals.empty(),
          "failed evaluation no payload");
    auto gg = gaussian(6);
    auto low = DesignProfile::prepare(std::move(gg), x, ids, columns(), tiny);
    check(low.numerical_status() == irred::numerics::Status::work_limit &&
              gg.status() == DensityStatus::finite,
          "resource failure preserves source");
    auto md = columns();
    md.shared_nuisance_ids.push_back("zero-point");
    check(DesignProfile::prepare(std::move(gg), x, ids, md).status() ==
              DensityStatus::invalid_input,
          "duplicate shared identity");
    md = columns();
    md.shared_nuisance_ids = {"unknown"};
    check(DesignProfile::prepare(std::move(gg), x, ids, md).status() ==
              DensityStatus::invalid_input,
          "unknown shared identity");
    auto zero = x;
    for (std::size_t i = 0; i < 6; ++i)
      zero[i * 3] = 0;
    auto singular = DesignProfile::prepare(std::move(gg), zero, ids, columns());
    check(singular.rank() == DesignRank::deficient &&
              singular.status() != DensityStatus::finite,
          "zero column exact deficient");
    auto duplicate = x;
    for (std::size_t i = 0; i < 6; ++i)
      duplicate[i * 3 + 1] = duplicate[i * 3];
    auto dup = DesignProfile::prepare(std::move(gg), duplicate, ids, columns());
    check(dup.rank() == DesignRank::unresolved &&
              dup.status() != DensityStatus::finite,
          "singular unresolved no dropped columns");
    auto almost = duplicate;
    almost[1] += 1e-8;
    auto ill = DesignProfile::prepare(std::move(gg), almost, ids, columns());
    check(ill.rank() == DesignRank::unresolved &&
              ill.status() != DensityStatus::finite,
          "near collinear unresolved");
    auto scaled = x;
    for (std::size_t i = 0; i < 6; ++i)
      scaled[i * 3] *= 1e100;
    md = columns();
    md.parameter_units[0] = "1e-100 mag";
    auto sp = DesignProfile::prepare(std::move(gg), scaled, ids, md);
    auto sv = sp.evaluate(r, ids);
    check(sv.status == DensityStatus::finite, "unit scaling accepted");
    near(sv.coefficients[0] * 1e100, 2, "unit scaling coefficient");
    near(sv.quadratic, 4, "unit scaling objective");
    auto prior = gaussian(6);
    std::vector<double> response(6, 1);
    prior = prior.proper_offset(response, ids, 0, 1, "proper shared prior",
                                true, 36, 1e-10);
    check(
        DesignProfile::prepare(std::move(prior), x, ids, columns()).status() ==
            DensityStatus::incompatible_metadata,
        "proper prior route explicit rejection");
    check(prior.status() == DensityStatus::finite,
          "prior reject preserves source");
    check(!DesignProfile::preparation_payload_bound(prior, SIZE_MAX, columns()),
          "overflow bound");
    std::printf("gaussian design owner: %d checks passed; named synthetic "
                "controls only\n",
                checks);
    return 0;
  } catch (const std::exception &e) {
    std::fprintf(stderr, "FAIL: %s (%d checks)\n", e.what(), checks);
    return 1;
  }
}
