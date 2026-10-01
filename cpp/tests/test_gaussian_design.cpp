// Synthetic linear controls only. Frozen coefficient/residual allocation:
// 2e-12 absolute + 2e-12 relative; quadratic 1e-10 absolute, with 1e-12
// independent small exact controls. Rational references derived by KKT below.
#include "irred/gaussian_design.hpp"
#include <algorithm>
#include <cfenv>
#include <cmath>
#include <cstdio>
#include <cstring>
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
    check(std::strcmp(moved.method_id(),
                      "retained-whitened-pivoted-householder-qr/v1") == 0,
          "actual QR method identity");
    check(std::strcmp(moved.qr_arithmetic_id(), "longdouble-cpu/v1") == 0,
          "wide QR arithmetic identity");
    check(moved.equilibrated_triangular_condition_inf() >= 1,
          "triangular R conditioning reported");
    const auto owned = moved.retained_payload_bound();
    check(owned && *owned > 6 * 3 * sizeof(double), "owned factor QR metadata");
    const auto eval_bound = moved.evaluation_payload_bound();
    DesignPolicy exact_quota;
    exact_quota.maximum_payload_bytes = *eval_bound;
    check(moved.evaluate(r, ids, exact_quota).status == DensityStatus::finite,
          "evaluation exact byte envelope");
    --exact_quota.maximum_payload_bytes;
    check(moved.evaluate(r, ids, exact_quota).numerical_status ==
              irred::numerics::Status::work_limit,
          "evaluation envelope minus one rejected");
    auto quota_source = gaussian(6);
    const auto prep_bound =
        DesignProfile::preparation_payload_bound(quota_source, 3, columns());
    exact_quota.maximum_payload_bytes = *prep_bound - 1;
    check(DesignProfile::prepare(std::move(quota_source), x, ids, columns(),
                                 exact_quota)
                      .numerical_status() ==
                  irred::numerics::Status::work_limit &&
              quota_source.status() == DensityStatus::finite,
          "preparation envelope minus one retains source");
    exact_quota.maximum_payload_bytes = *prep_bound;
    exact_quota.maximum_elements = x.size();
    auto quota_profile = DesignProfile::prepare(std::move(quota_source), x, ids,
                                                columns(), exact_quota);
    check(quota_profile.status() == DensityStatus::finite,
          "exact preparation envelope and elements");
    auto caller_x = x;
    auto lifetime_source = gaussian(6);
    auto lifetime = DesignProfile::prepare(std::move(lifetime_source), caller_x,
                                           ids, columns());
    caller_x.assign(caller_x.size(), NAN);
    near(lifetime.evaluate(r, ids).quadratic, 4,
         "retains no caller design buffer");
    auto extreme = x;
    for (std::size_t i = 0; i < 6; ++i) {
      extreme[i * 3] = std::ldexp(extreme[i * 3], 900);
      extreme[i * 3 + 1] = std::ldexp(extreme[i * 3 + 1], -900);
    }
    auto extreme_source = gaussian(6);
    auto extreme_profile = DesignProfile::prepare(std::move(extreme_source),
                                                  extreme, ids, columns());
    const auto extreme_fit = extreme_profile.evaluate(r, ids);
    check(extreme_fit.status == DensityStatus::finite,
          "wide coordinate scaling");
    near(std::ldexp(extreme_fit.coefficients[0], 900), 2,
         "large design coordinate");
    near(std::ldexp(extreme_fit.coefficients[1], -900), -3,
         "small design coordinate");
    near(extreme_fit.quadratic, 4, "extreme scaling actual objective");
    std::vector<double> zeros(6, 0);
    const auto zero_fit = moved.evaluate(zeros, ids);
    check(zero_fit.status == DensityStatus::finite && zero_fit.quadratic == 0 &&
              zero_fit.normalized_normal_equation_residual == 0,
          "zero residual absolute defect fallback");
    auto too_small = r;
    for (auto &v : too_small)
      v *= 1e-200;
    const auto underflow = moved.evaluate(too_small, ids);
    check(underflow.status == DensityStatus::numerical_failure &&
              underflow.coefficients.empty() &&
              underflow.adjusted_residuals.empty(),
          "nonzero subnormal objective withheld");
    auto environment_source = gaussian(6);
    const auto previous_round = std::fegetround();
    check(std::fesetround(FE_UPWARD) == 0, "set adversarial rounding");
    const auto environment_profile = DesignProfile::prepare(
        std::move(environment_source), x, ids, columns());
    const auto environment_fit = moved.evaluate(r, ids);
    check(std::fesetround(previous_round) == 0, "restore rounding");
    check(environment_profile.status() == DensityStatus::unsupported_domain &&
              environment_source.status() == DensityStatus::finite &&
              environment_fit.status == DensityStatus::unsupported_domain &&
              environment_fit.coefficients.empty(),
          "unsupported rounding retains source and withholds output");
    // Independent exact triangular witness L=(2;1,3;-1,2,4), z=(1,-2,3).
    const std::vector<double> c = {4, 2, -2, 2, 10, 5, -2, 5, 21};
    const std::vector<double> rhs = {2, -5, 7};
    for (const auto arithmetic :
         {irred::numerics::Arithmetic::binary64_legacy_v1,
          irred::numerics::Arithmetic::longdouble_cpu_v1}) {
      auto factor = irred::numerics::cholesky(c, 3, 9, arithmetic);
      const auto wb = irred::numerics::whitening_payload_bound(3);
      const auto whitened = irred::numerics::whiten(factor, rhs, 3, *wb, 1e-10);
      check(whitened.status == irred::numerics::Status::ok,
            "retained factor triangular whitening");
      check(whitened.value[0] == 1 && whitened.value[1] == -2 &&
                whitened.value[2] == 3,
            "exact independent triangular values");
      check(whitened.backward_residual == 0 &&
                whitened.arithmetic_rounding_estimate > 0,
            "backward vs arithmetic rounding diagnostics");
      const auto short_bytes =
          irred::numerics::whiten(factor, rhs, 3, *wb - 1, 1e-10);
      check(short_bytes.status == irred::numerics::Status::work_limit &&
                short_bytes.value.empty(),
            "whitening byte envelope minus one");
      const auto short_elements =
          irred::numerics::whiten(factor, rhs, 2, *wb, 1e-10);
      check(short_elements.status == irred::numerics::Status::work_limit &&
                short_elements.value.empty(),
            "whitening element quota");
      const auto tight = irred::numerics::whiten(factor, rhs, 3, *wb, 1e-30);
      check(tight.status ==
                    irred::numerics::Status::conditioning_budget_exceeded &&
                tight.value.empty(),
            "whitening diagnostic failure no payload");
      auto poison = rhs;
      poison[0] = INFINITY;
      check(
          irred::numerics::whiten(factor, poison, 3, *wb, 1e-10).value.empty(),
          "nonfinite whitening no payload");
    }
    check(!irred::numerics::whitening_payload_bound(SIZE_MAX),
          "whitening overflow bound");
    std::printf("gaussian design owner: %d checks passed; named synthetic "
                "controls only\n",
                checks);
    return 0;
  } catch (const std::exception &e) {
    std::fprintf(stderr, "FAIL: %s (%d checks)\n", e.what(), checks);
    return 1;
  }
}
