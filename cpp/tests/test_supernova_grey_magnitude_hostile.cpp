// Independent synthetic grey-magnitude contract. Three-dimensional exact
// adjugate covariance, analytic de Sitter geometry; shared libm disclosed.
// Frozen assembled absolute budget 1e-10, magnitude budget 2e-14. No inference.
#include "irred/supernova.hpp"
#include <array>
#include <bit>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
using namespace irred;
namespace {
int checks = 0;
long double worst = 0;
void require(bool v, const char *why) {
  ++checks;
  if (!v)
    throw std::runtime_error(why);
}
void close(long double a, long double b, long double budget = 1e-10L) {
  worst = std::max(worst, std::abs(a - b));
  require(std::abs(a - b) <= budget, "independent analytic budget");
}
auto bits(double x) { return std::bit_cast<std::uint64_t>(x); }
supernova::Policy policy() {
  supernova::Policy p;
  p.arithmetic = numerics::Arithmetic::longdouble_cpu_v1;
  p.maximum_forward_sensitivity = 1e-10;
  p.background.maximum_total_evaluations = 10000;
  return p;
}
constexpr std::array<double, 3> z{.125, .5, 1}, observer{.3, .75, 1.25};
constexpr std::array<long double, 3> noise{.02L, -.03L, .01L};
long double shape(std::size_t i, bool constant) {
  return 5 *
         std::log10((1 + (long double)observer[i]) * (constant ? .5L : z[i]));
}
long double basis(std::size_t i, bool constant) {
  return std::log1p(constant ? .5L : (long double)z[i]) / std::log(2.L);
}
supernova::Consumer make(bool constant = false, double intercept = 17) {
  observations::Input in{};
  in.profile = observations::Profile::gaussian_fixture_v1;
  in.role = observations::Role::synthetic_control;
  in.unit = observations::Unit::magnitude;
  in.calibration = observations::Calibration::not_applicable;
  in.uncertainty = observations::Uncertainty::covariance;
  in.uncertainty_unit = observations::UncertaintyUnit::magnitude_squared;
  in.component = observations::Component::total;
  in.ordering_provenance = "synthetic independent adjugate order";
  in.table_sha256 = std::string(64, 'c');
  in.uncertainty_sha256 = std::string(64, 'd');
  in.measurement_ids = {"short", "peer-long-identity-aaaaaaaaaaaaaaaa",
                        "third"};
  in.event_ids = in.measurement_ids;
  in.uncertainty_axis_ids = in.measurement_ids;
  in.uncertainty_matrix = {4, 1, 0, 1, 3, 1, 0, 1, 2};
  in.missing = {0, 0, 0};
  in.source_selection = {1, 1, 1};
  in.quality = {0, 0, 0};
  std::array<cosmology::Query, 3> qs;
  for (std::size_t i = 0; i < 3; ++i) {
    qs[i] = {constant ? .5 : z[i], observer[i],
             cosmology::Convention::released_zhd_zhel};
    in.values.push_back(static_cast<double>(
        shape(i, constant) + intercept + .13L * basis(i, constant) + noise[i]));
  }
  auto ids = in.measurement_ids;
  auto source = observations::prepare(std::move(in), {3, 9, 4096});
  require(source.status() == observations::Status::ok, "typed synthetic");
  auto c = supernova::prepare_synthetic(std::move(source), ids, qs, policy());
  require(c.status() == supernova::Status::ok, "retained peer consumer");
  return c;
}
std::pair<long double, long double>
exact_profile(std::array<long double, 3> r) {
  // C^-1 = [[5,-2,1],[-2,8,-4],[1,-4,11]]/18;
  // W1=[4,2,8]/18 and 1'W1=14/18. No production solves.
  auto h = 4 * r[0] + 2 * r[1] + 8 * r[2];
  auto q = (5 * r[0] * r[0] + 8 * r[1] * r[1] + 11 * r[2] * r[2] -
            4 * r[0] * r[1] + 2 * r[0] * r[2] - 8 * r[1] * r[2] - h * h / 14) /
           18;
  return {q, h / 14};
}
void absent(const supernova::GreyMagnitudeSlot &s) {
  require(s.calculation.status != supernova::Status::ok, "failed slot");
  require(s.magnitude_shifts.empty() &&
              s.calculation.shape_magnitudes.empty() &&
              s.calculation.base_residuals.empty() &&
              s.calculation.profiled_residuals.empty(),
          "no failed finite arrays");
}
} // namespace
int main() {
  try {
    auto c = make();
    supernova::GreyMagnitudePolicy p{policy(), 1048576};
    const supernova::ModelPointV2 m{cosmology::Model::constant_q_flat_v1, 0, -1,
                                    -1, 0};
    std::array<supernova::ModelPointV2, 1> legacy{m};
    auto old = c.evaluate_batch_v2(legacy, p.evaluation);
    for (double eps : {-.5, -.2, -0., 0., .13, .2, .5}) {
      std::array<supernova::GreyMagnitudePoint, 1> points{{{m, eps}}};
      auto out = c.evaluate_grey_magnitude_batch(points, p);
      require(out.status == supernova::Status::ok && out.slots.size() == 1,
              "owned batch");
      auto &s = out.slots[0];
      require(s.calculation.status == supernova::Status::ok, "finite grey");
      require(bits(s.source.epsilon_mag) == bits(eps), "raw epsilon bits");
      std::array<long double, 3> r;
      for (std::size_t i = 0; i < 3; ++i) {
        r[i] = 17 + (.13L - eps) * basis(i, false) + noise[i];
        close(s.magnitude_shifts[i], eps * basis(i, false), 2e-14L);
        close(s.calculation.shape_magnitudes[i], shape(i, false), 2e-14L);
        close(s.calculation.base_residuals[i], r[i], 2e-14L);
        if (eps != 0)
          require(std::abs(s.magnitude_shifts[i] -
                           eps * std::log1p(observer[i]) / std::log(2.)) > 1e-3,
                  "zhd distinct from observer");
      }
      auto [q, a] = exact_profile(r);
      close(s.calculation.quadratic, q);
      close(s.calculation.relative_profile_score, -q / 2);
      close(s.calculation.offset_coefficient, a);
      if (eps == 0) {
        require(bits(s.calculation.quadratic) ==
                    bits(old.slots[0].calculation.quadratic),
                "zero q exact");
        require(s.calculation.base_residuals ==
                    old.slots[0].calculation.base_residuals,
                "zero residual exact");
        require(s.calculation.profiled_residuals ==
                    old.slots[0].calculation.profiled_residuals,
                "zero canonical exact");
        require(s.calculation.background_evaluations ==
                    old.slots[0].calculation.background_evaluations,
                "zero work exact");
      }
    }
    auto same = make(true), shifted = make(true, 117);
    for (double eps : {-.5, 0., .5}) {
      std::array<supernova::GreyMagnitudePoint, 1> points{{{m, eps}}};
      auto a = same.evaluate_grey_magnitude_batch(points, p),
           b = shifted.evaluate_grey_magnitude_batch(points, p);
      auto expected = exact_profile(noise).first;
      close(a.slots[0].calculation.quadratic, expected);
      close(b.slots[0].calculation.quadratic, expected);
      close(b.slots[0].calculation.offset_coefficient -
                a.slots[0].calculation.offset_coefficient,
            100);
    }
    for (double eps : {std::nextafter(.5, 1.), std::nextafter(-.5, -1.),
                       std::numeric_limits<double>::infinity(),
                       std::numeric_limits<double>::quiet_NaN(),
                       std::numeric_limits<double>::denorm_min()}) {
      std::array<supernova::GreyMagnitudePoint, 1> points{{{m, eps}}};
      auto out = c.evaluate_grey_magnitude_batch(points, p);
      require(out.slots.size() == 1, "attempt materialized");
      absent(out.slots[0]);
      require(bits(out.slots[0].source.epsilon_mag) == bits(eps),
              "failed source bits");
      if (eps == std::numeric_limits<double>::denorm_min())
        require(out.slots[0].calculation.numerical_status ==
                    numerics::Status::outside_domain,
                "tiny shift normal-output policy");
    }
    std::array<supernova::GreyMagnitudePoint, 2> points{{{m, .2}, {m, -.2}}};
    auto no = p;
    no.maximum_native_bytes = 0;
    auto out = c.evaluate_grey_magnitude_batch(points, no);
    require(out.status == supernova::Status::work_limit && out.slots.empty(),
            "preallocation bytecap");
    // Published conservative storage contract, excluding prepared storage,
    // caller inputs and allocator overhead. N=3, P=2 includes suppressed or
    // later failed payload allocation budgets as well as retained vectors.
    const auto bytes = 2 * sizeof(supernova::GreyMagnitudeSlot) +
                       5 * 2 * 3 * sizeof(double) +
                       3 * (sizeof(cosmology::Slot) + 16 * sizeof(double) +
                            8 * sizeof(long double));
    no = p;
    no.maximum_native_bytes = bytes - 1;
    out = c.evaluate_grey_magnitude_batch(points, no);
    require(out.status == supernova::Status::work_limit && out.slots.empty(),
            "declared storage minus one rejected");
    no.maximum_native_bytes = bytes;
    out = c.evaluate_grey_magnitude_batch(points, no);
    require(out.status == supernova::Status::ok && out.slots.size() == 2 &&
                out.slots[0].calculation.status == supernova::Status::ok &&
                out.slots[1].calculation.status == supernova::Status::ok,
            "declared storage exactly admitted");
    no = p;
    no.evaluation.maximum_models = 0;
    out = c.evaluate_grey_magnitude_batch(points, no);
    require(out.status == supernova::Status::work_limit && out.slots.empty(),
            "modelcap");
    no = p;
    no.evaluation.arithmetic = numerics::Arithmetic::binary64_legacy_v1;
    out = c.evaluate_grey_magnitude_batch(points, no);
    require(out.status == supernova::Status::incompatible_metadata &&
                out.slots.empty(),
            "precision mismatch");
    auto first = c.evaluate_grey_magnitude_batch(std::span(points).first(1), p);
    auto work = first.slots[0].calculation.background_evaluations;
    require(work > 0, "actual work");
    no = p;
    no.evaluation.background.maximum_total_evaluations = work;
    out = c.evaluate_grey_magnitude_batch(points, no);
    require(out.slots.size() == 2 &&
                out.slots[0].calculation.status == supernova::Status::ok,
            "first uses shared budget");
    absent(out.slots[1]);
    require(out.slots[1].calculation.status == supernova::Status::work_limit &&
                out.slots[1].calculation.numerical_status ==
                    numerics::Status::work_limit,
            "precise globalwork cause");
    std::printf("grey magnitude independent PASS %d maxerror %.18Lg\n", checks,
                worst);
    return 0;
  } catch (const std::exception &e) {
    std::fprintf(stderr, "grey independent FAIL %s after%d\n", e.what(),
                 checks);
    return 1;
  }
}
