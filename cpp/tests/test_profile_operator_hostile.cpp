// Independent cached-offset contract. C=[[4,1],[1,9]], C^-1*1=(8,3)/35,
// Gram=11/35; r=(2,-3) gives a=7/11 and q=25/11 by exact adjugate.
// Binary64 cannot represent both exact response fractions; stricter later
// budgets must assess the retained response solve, not just a new zero
// residual.
#include "irred/statistics.hpp"
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <new>
#include <stdexcept>
namespace {
bool trace_allocations = false;
std::size_t allocation_count = 0, allocation_bytes = 0;
} // namespace
void *operator new(std::size_t n) {
  if (trace_allocations) {
    ++allocation_count;
    allocation_bytes += n;
  }
  if (auto p = std::malloc(n ? n : 1))
    return p;
  throw std::bad_alloc();
}
void operator delete(void *p) noexcept { std::free(p); }
void operator delete(void *p, std::size_t) noexcept { std::free(p); }
using namespace irred;
namespace {
int checks = 0;
void check(bool ok, const char *why) {
  ++checks;
  if (!ok)
    throw std::runtime_error(why);
}
statistics::Metadata metadata() {
  statistics::Metadata m;
  m.ordered_ids = {"synthetic:A", "synthetic:B"};
  m.measure = "product d(magnitude)";
  m.ordering_provenance = "explicit ordered synthetic fixture IDs";
  return m;
}
observations::Prepared observations_fixture() {
  observations::Input s{};
  s.profile = observations::Profile::gaussian_fixture_v1;
  s.role = observations::Role::synthetic_control;
  s.unit = observations::Unit::magnitude;
  s.calibration = observations::Calibration::not_applicable;
  s.uncertainty = observations::Uncertainty::covariance;
  s.uncertainty_unit = observations::UncertaintyUnit::magnitude_squared;
  s.component = observations::Component::total;
  s.table_sha256 = std::string(64, 'a');
  s.uncertainty_sha256 = std::string(64, 'b');
  s.measurement_ids = {"synthetic:A", "synthetic:B"};
  s.event_ids = {"A", "B"};
  s.uncertainty_axis_ids = s.measurement_ids;
  s.values = {2, -3};
  s.missing = {0, 0};
  s.quality = {0, 0};
  s.source_selection = {1, 1};
  s.ordering_provenance = "generated exact covariance axis order";
  s.uncertainty_matrix = {4, 1, 1, 9};
  return observations::prepare(std::move(s), {2, 4, 4096});
}
} // namespace
int main() {
  auto m = metadata();
  std::array<double, 4> c{4, 1, 1, 9};
  std::array<double, 2> ones{1, 1}, r{2, -3}, zero{0, 0};
  auto g = statistics::prepare_gaussian(c, statistics::MatrixKind::covariance,
                                        m, 4, 1e-10);
  auto op = std::move(g).prepare_offset_profile(ones, m.ordered_ids, 1e-10);
  check(op.status() == statistics::DensityStatus::finite,
        "permissive prepared cache");
  check(g.status() != statistics::DensityStatus::finite,
        "moved source cannot remain scientifically valid");
  check(std::abs(op.cached_response_solution()[0] - 8.L / 35) < 1e-15,
        "independent response adjugate");
  check(std::abs(op.gram() - 11.L / 35) < 1e-15, "independent Gram");
  auto accepted = op.evaluate(r, m.ordered_ids, 1e-10);
  check(accepted.status == statistics::DensityStatus::finite,
        "ordinary cached score");
  check(accepted.numerical_status == numerics::Status::ok,
        "successful profile preserves F02 cause");
  check(std::abs(accepted.coefficient - 7.L / 11) < 1e-10 &&
            std::abs(accepted.quadratic - 25.L / 11) < 1e-10,
        "rational score/offset");
  const auto retained = op.cached_response_forward_sensitivity();
  check(retained > 0 && std::isfinite(retained),
        "nonrepresentable rational response has nonzero diagnostic");
  auto factor = numerics::cholesky(c, 2, 4);
  auto zero_reference = numerics::solve(factor, zero, 1e-10);
  check(zero_reference.status == numerics::Status::ok &&
            retained > zero_reference.estimated_forward_sensitivity,
        "cache diagnostic strictly exceeds zero solve diagnostic");
  const auto narrower =
      (retained + zero_reference.estimated_forward_sensitivity) / 2;
  check(narrower > zero_reference.estimated_forward_sensitivity &&
            narrower < retained,
        "requested budget distinguishes retained solve from new zero solves");
  check(numerics::solve(factor, zero, narrower).status == numerics::Status::ok,
        "new zero solve would pass narrower budget");
  auto tighter = op.evaluate(zero, m.ordered_ids, narrower);
  check(tighter.status == statistics::DensityStatus::numerical_failure,
        "zero new solve cannot hide weaker retained response budget");
  check(tighter.numerical_status ==
            numerics::Status::conditioning_budget_exceeded,
        "conditioning failure preserves F02 cause");
  check(op.evaluate(zero, m.ordered_ids, 0).status !=
            statistics::DensityStatus::finite,
        "zero requested policy invalid");
  auto source = observations_fixture();
  // Real valid span larger than source rows; no invalid span precondition or
  // inaccessible-pointer assumption. Allocation tracing distinguishes bounded
  // error metadata from row/index copies.
  std::array<std::size_t, 3> impossible{0, 1, 0};
  allocation_count = 0;
  allocation_bytes = 0;
  trace_allocations = true;
  statistics::Gaussian empty_error_object;
  trace_allocations = false;
  const auto baseline_count = allocation_count,
             baseline_bytes = allocation_bytes;
  allocation_count = 0;
  allocation_bytes = 0;
  trace_allocations = true;
  auto rejected =
      statistics::prepare_selected_observations(source, impossible, 4, 1e-10);
  trace_allocations = false;
  check(rejected.status() != statistics::DensityStatus::finite &&
            allocation_count == baseline_count &&
            allocation_bytes == baseline_bytes,
        "impossible source index count copies no rows/indices beyond bounded "
        "error object");
  std::array<std::size_t, 2> reversed{1, 0}, duplicated{0, 0}, outside{0, 2},
      valid{0, 1};
  for (auto indices : {reversed, duplicated, outside})
    check(statistics::prepare_selected_observations(source, indices, 4, 1e-10)
                  .status() != statistics::DensityStatus::finite,
          "invalid source-order selection rejects");
  check(statistics::prepare_selected_observations(source, valid, 3, 1e-10)
                .status() != statistics::DensityStatus::finite,
        "selected matrix cap preallocation");
  check(statistics::prepare_selected_observations(source, valid, 4, 1e-10)
                .status() == statistics::DensityStatus::finite,
        "ordinary selected covariance accepts");
  std::printf("{\"suite\":\"independent_profile_operator\",\"checks\":%d,"
              "\"cached_sensitivity\":%.17g,\"passed\":true}\n",
              checks, retained);
}
