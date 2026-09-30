// Independent peer: redshift-segment composite Simpson and exact 2x2
// cofactor Gaussian, versus production analytic segments and Cholesky.
// Observable budget 2e-12+2e-10|reference|; refinement <=10% budget.
// Normalized synthetic density budget1e-8, original24 separately assessed.
#include "irred/bao.hpp"
#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
using namespace irred;
namespace {
unsigned checks = 0;
long double worst = 0, refinement = 0, density_error = 0;
void check(bool b, const char *why) {
  ++checks;
  if (!b)
    throw std::runtime_error(why);
}
void close(double x, long double y) {
  auto budget = 2e-12L + 2e-10L * std::abs(y);
  worst = std::max(worst, std::abs(x - y) / budget);
  check(std::abs(x - y) <= budget, "observable budget");
}
constexpr std::array<long double, 6> edges{0, .1L, .3L, .6L, 1, 2.5L};
long double expansion(long double z, const std::array<double, 5> &q) {
  long double e = 1;
  for (unsigned k = 0; k < 5; ++k) {
    auto end = std::min(z, edges[k + 1]);
    e *= std::pow((1 + end) / (1 + edges[k]), 1 + (long double)q[k]);
    if (z <= edges[k + 1])
      return e;
  }
  throw std::runtime_error("reference domain");
}
long double integral(double z, const std::array<double, 5> &q, unsigned n) {
  long double sum = 0;
  for (unsigned k = 0; k < 5 && z > edges[k]; ++k) {
    auto lo = edges[k], hi = std::min((long double)z, edges[k + 1]),
         h = (hi - lo) / n, s = 0.L;
    for (unsigned j = 0; j <= n; ++j)
      s += (j == 0 || j == n ? 1 : j % 2 ? 4 : 2) / expansion(lo + j * h, q);
    sum += s * h / 3;
  }
  return sum;
}
long double exact(double z, double q) {
  auto u = 1 + (long double)z;
  if (q == -3)
    return (u * u * u - 1) / 3;
  if (q == 2)
    return (1 - 1 / (u * u)) / 2;
  if (q == 0)
    return std::log(u);
  return (1 - std::pow(u, -q)) / q;
}
std::array<std::array<double, 5>, 8> grid{{{0, 0, 0, 0, 0},
                                           {-1, -1, -1, -1, -1},
                                           {.5, .5, .5, .5, .5},
                                           {-.4, -.4, -.2, .1, .3},
                                           {-1, 0, -1, 0, -1},
                                           {-3, 2, -3, 2, -3},
                                           {-3, -3, -3, -3, -3},
                                           {2, 2, 2, 2, 2}}};
bao::DensityInput input() {
  bao::DensityInput x;
  x.queries = {{.51, bao::Observable::transverse_over_ruler},
               {2.33, bao::Observable::hubble_over_ruler}};
  x.observed = {20, 9};
  x.covariance = {4, 1, 1, 9};
  x.ordered_ids = {"DM-control", "DH-control"};
  x.role = bao::RowRole::synthetic_control;
  x.covariance_unit = bao::CovarianceUnit::dimensionless_ratio_squared;
  x.table_identity = "independent generated two-coordinate control";
  x.covariance_identity = "integer covariance determinant35";
  x.ordering_provenance = "explicit generated axis order";
  x.calibration_provenance = "unknown";
  x.dependence_provenance = "unknown";
  return x;
}
bao::DensityPolicy preparation() {
  bao::DensityPolicy p;
  p.maximum_models = 30;
  p.maximum_matrix_elements = 4;
  p.maximum_string_bytes = 10000;
  p.maximum_native_bytes = 1000000;
  p.maximum_forward_sensitivity = 1e-10;
  p.arithmetic = numerics::Arithmetic::longdouble_cpu_v1;
  return p;
}
bao::PiecewiseDensityPolicy evaluation() {
  bao::PiecewiseDensityPolicy p;
  p.maximum_models = 30;
  p.maximum_native_bytes = 1000000;
  p.maximum_forward_sensitivity = 1e-10;
  p.arithmetic = numerics::Arithmetic::longdouble_cpu_v1;
  p.observables.background.maximum_segment_visits = 2000;
  return p;
}
} // namespace
int main() {
  try {
    for (auto q : grid)
      for (double ruler : {5000., 10000., 15000.}) {
        auto bg = cosmology::prepare_piecewise_q({70, q});
        check(bg.status() == cosmology::Status::ok, "prepare reference grid");
        for (double z : {0., .0001, .1, .3, .6, 1., 2.33, 2.5}) {
          std::array<bao::Query, 3> queries{
              {{z, bao::Observable::transverse_over_ruler},
               {z, bao::Observable::hubble_over_ruler},
               {z, bao::Observable::volume_over_ruler}}};
          bao::PiecewisePolicy p;
          auto b = bao::evaluate_piecewise(bg, bao::Ruler(ruler), queries, p);
          check(b.status == cosmology::Status::ok && b.slots.size() == 3,
                "observable batch");
          auto a = integral(z, q, 4096), r = integral(z, q, 8192),
               budget = 2e-12L + 2e-10L * std::abs(r);
          refinement = std::max(refinement, std::abs(a - r) / budget);
          check(std::abs(a - r) <= .1L * budget, "redshift Simpson refinement");
          if (std::all_of(q.begin(), q.end(),
                          [&](double x) { return x == q[0]; })) {
            auto analytic = exact(z, q[0]);
            check(std::abs(r - analytic) <= .1L * budget,
                  "analytic integral control");
            r = analytic;
          }
          const long double scale = 299792.458L / ruler, dm = scale * r,
                            dh = scale / expansion(z, q),
                            dv = std::cbrt(z * dm * dm * dh);
          for (auto &s : b.slots)
            check(s.status == cosmology::Status::ok, "observable finite");
          close(b.slots[0].dimensionless_value, dm);
          close(b.slots[1].dimensionless_value, dh);
          close(b.slots[2].dimensionless_value, dv);
        }
      }
    auto source = input();
    auto edge_bg = cosmology::prepare_piecewise_q({70, grid[5]});
    std::array<bao::Query, 3> malformed_queries{
        {{.3, static_cast<bao::Observable>(99)},
         {std::nextafter(2.5, std::numeric_limits<double>::infinity()),
          bao::Observable::transverse_over_ruler},
         {.3, bao::Observable::transverse_over_ruler}}};
    bao::PiecewisePolicy edge_policy;
    auto edge_result = bao::evaluate_piecewise(edge_bg, bao::Ruler(10000),
                                               malformed_queries, edge_policy);
    check(edge_result.status == cosmology::Status::ok &&
              edge_result.slots[0].status == cosmology::Status::invalid_input &&
              edge_result.slots[1].status ==
                  cosmology::Status::unsupported_domain &&
              edge_result.slots[2].status == cosmology::Status::ok,
          "unknown tag and outside-domain source do not hide valid jump "
          "geometry");
    edge_policy.maximum_queries = 0;
    check(bao::evaluate_piecewise(edge_bg, bao::Ruler(10000), malformed_queries,
                                  edge_policy)
              .slots.empty(),
          "query cap rejects before copies");
    auto dp = preparation();
    auto prepared = bao::prepare_density(source, dp);
    check(prepared.status() == statistics::DensityStatus::finite,
          "density prepared");
    auto policy = evaluation();
    std::vector<bao::PiecewiseModelPoint> points;
    for (auto q : grid)
      for (double ruler : {5000., 10000., 15000.})
        points.emplace_back(q, bao::Ruler(ruler));
    auto out = prepared.evaluate_piecewise(points, policy);
    check(out.slots.size() == 24, "all24 materialized");
    size_t total = 0;
    for (size_t i = 0; i < 24; ++i) {
      auto &s = out.slots[i];
      check(s.result.density.status == statistics::DensityStatus::finite,
            "density finite");
      check(s.source.q == points[i].q &&
                s.source.ruler.h0_rd_km_s == points[i].ruler.h0_rd_km_s,
            "attempt retained");
      auto q = points[i].q;
      long double scale = 299792.458L / points[i].ruler.h0_rd_km_s,
                  r0 = 20 - scale * integral(.51, q, 16384),
                  r1 = 9 - scale / expansion(2.33, q);
      auto quad = (9 * r0 * r0 - 2 * r0 * r1 + 4 * r1 * r1) / 35;
      auto target =
          -.5L * (quad + std::log(35.L) + 2 * std::log(2 * std::acos(-1.L)));
      density_error = std::max(density_error,
                               std::abs(s.result.density.log_value - target));
      check(std::abs(s.result.density.log_value - target) <= 1e-8L,
            "cofactor normalized density");
      check(std::abs(s.result.log_determinant - std::log(35.L)) < 1e-13L,
            "determinant normalization");
      check(s.predictions.size() == 2 && s.residuals.size() == 2,
            "owned payload");
      total += s.segment_visits;
    }
    check(out.segment_visits == total && total <= 2000, "global work count");
    auto limited = policy;
    limited.observables.background.maximum_segment_visits = 8;
    auto partial = prepared.evaluate_piecewise(points, limited);
    check(partial.segment_visits == 8 &&
              partial.slots[0].result.density.status ==
                  statistics::DensityStatus::finite,
          "first model consumes shared analytic budget");
    check(partial.slots[1].segment_visits == 0 &&
              partial.slots[1].numerical_status ==
                  numerics::Status::work_limit &&
              partial.slots[1].predictions.empty(),
          "later model atomic rejection does not reset budget");
    auto tiny = source;
    tiny.queries[0].z = 1e-200;
    auto tiny_owner = bao::prepare_density(tiny, dp);
    check(tiny_owner.status() == statistics::DensityStatus::finite,
          "tiny source structurally retained");
    limited.observables.background.maximum_segment_visits = 6;
    auto failed_work = tiny_owner.evaluate_piecewise(points, limited);
    check(failed_work.segment_visits == 6 &&
              failed_work.slots[0].segment_visits == 6,
          "failed computed geometry consumes analytic work");
    check(failed_work.slots[0].result.density.status !=
                  statistics::DensityStatus::finite &&
              failed_work.slots[0].predictions.empty() &&
              failed_work.slots[0].residuals.empty(),
          "failed geometry cannot publish density payload");
    auto bad = policy;
    bad.maximum_models = 0;
    auto capped = prepared.evaluate_piecewise(points, bad);
    check(capped.slots.empty() &&
              capped.numerical_status == numerics::Status::work_limit,
          "model cap preflight");
    bad = policy;
    bad.maximum_native_bytes = 0;
    capped = prepared.evaluate_piecewise(points, bad);
    check(capped.slots.empty() &&
              capped.numerical_status == numerics::Status::work_limit,
          "byte cap preflight");
    bao::PiecewiseDensityPolicy unset;
    check(prepared.evaluate_piecewise(points, unset).slots.empty(),
          "required precision policy");
    bad = policy;
    bad.arithmetic = numerics::Arithmetic::binary64_legacy_v1;
    check(prepared.evaluate_piecewise(points, bad).slots.empty(),
          "retained precision mismatch");
    auto attempts = points;
    attempts[0].q[0] = 3;
    attempts[1].ruler.h0_rd_km_s = 4999;
    auto failures = prepared.evaluate_piecewise(attempts, policy);
    for (size_t i = 0; i < 2; ++i) {
      check(failures.slots[i].result.density.status !=
                statistics::DensityStatus::finite,
            "domain failure");
      check(failures.slots[i].predictions.empty() &&
                failures.slots[i].residuals.empty(),
            "failed payload omitted");
    }
    check(failures.slots[0].source.q[0] == 3, "invalid source retained");
    auto moved = std::move(prepared);
    check(moved.evaluate_piecewise(points, policy)
                  .slots[0]
                  .result.density.status == statistics::DensityStatus::finite,
          "moved owner");
    auto abandoned = prepared.evaluate_piecewise(points, policy);
    for (const auto &slot : abandoned.slots) {
      check(slot.result.density.status != statistics::DensityStatus::finite,
            "moved-from cannot publish finite density");
      check(slot.predictions.empty() && slot.residuals.empty(),
            "moved-from cannot publish finite arrays");
    }
    bao::PreparedDensity empty;
    check(empty.evaluate_piecewise(points, policy).slots.empty(),
          "default owner");
    std::printf("piecewise BAO peer PASS %u observable_fraction %.18Lg "
                "refinement_fraction %.18Lg density_max_error %.18Lg\n",
                checks, worst, refinement, density_error);
  } catch (const std::exception &e) {
    std::fprintf(stderr, "piecewise BAO peer: %s\n", e.what());
    return 1;
  }
}
