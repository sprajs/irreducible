#include "irred/gaussian_predictive.hpp"
#include <array>
#include <atomic>
#include <cfenv>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <new>
#include <numbers>
#include <stdexcept>
namespace {
std::atomic<bool> armed = false;
std::atomic<size_t> calls = 0, fail_on = 0;
std::atomic<long> live = 0;
} // namespace
#if defined(__GNUC__) || defined(__clang__)
#define NOINLINE __attribute__((noinline))
#else
#define NOINLINE
#endif
NOINLINE void *operator new(size_t n) {
  if (armed && ++calls == fail_on)
    throw std::bad_alloc();
  if (void *p = std::malloc(n ? n : 1)) {
    ++live;
    return p;
  }
  throw std::bad_alloc();
}
NOINLINE void *operator new[](size_t n) { return ::operator new(n); }
NOINLINE void operator delete(void *p) noexcept {
  if (p) {
    --live;
    std::free(p);
  }
}
NOINLINE void operator delete[](void *p) noexcept { ::operator delete(p); }
NOINLINE void operator delete(void *p, size_t) noexcept {
  ::operator delete(p);
}
NOINLINE void operator delete[](void *p, size_t) noexcept {
  ::operator delete(p);
}
namespace {
using namespace irred::statistics;
using S = irred::numerics::Status;
using W = long double;
unsigned checks = 0;
void need(bool ok, const char *what) {
  ++checks;
  if (!ok)
    throw std::runtime_error(what);
}
void near(W actual, W expected, W budget = 2e-12L) {
  need(std::abs(actual - expected) <= budget * (1 + std::abs(expected)),
       "frozen original rational allocation");
}
const std::vector<std::string> rows{"train-cal0", "train-cal1"},
    future{"heldout-cal2", "heldout-cal3"};
const std::array<double, 2> training{1.25, -.5};
const std::array<double, 4> response{1, 2, 1, -1};
Metadata source_metadata(bool heldout = false) {
  Metadata m;
  m.ordered_ids = heldout ? future : rows;
  m.measure = "product d(mag)";
  m.table_identity = heldout ? "distinct synthetic future calibration"
                             : "synthetic training calibration";
  m.uncertainty_identity = "explicit conditional SPD covariance";
  m.ordering_provenance = "declared exact order";
  m.calibration_provenance =
      "shared offset/slope explicitly represented in beta";
  m.dependence_provenance = "conditional full noise covariance";
  m.source_semantics = "synthetic controls";
  return m;
}
ParameterPrior prior() {
  ParameterPrior p;
  p.ordered_parameter_ids = {"offset", "slope"};
  p.parameter_units = {"mag", "mag"};
  p.shared_nuisance_ids = {"offset"};
  p.mean = {.5, -.25};
  p.covariance = {2, .5, .5, 1};
  p.prior_identity = "chosen original synthetic proper prior";
  p.design_identity = "offset+slope fixed t=[0,1]";
  p.residual_unit = "mag";
  p.parameter_measure = "d(offset-mag) d(slope-mag)";
  p.dependence_identity = "proper prior independent of training noise";
  p.noise_independence_declared = true;
  return p;
}
GaussianPosterior posterior(std::array<double, 4> x = {1, 0, 1, 1},
                            ParameterPrior p = prior(),
                            std::array<double, 4> c = {1, .25, .25, 2}) {
  auto g =
      prepare_gaussian(c, MatrixKind::covariance, source_metadata(), 1000,
                       1e-10, irred::numerics::Arithmetic::longdouble_cpu_v1);
  return GaussianPosterior::prepare(std::move(g), x, rows, std::move(p));
}
Gaussian noise(std::array<double, 4> r = {.5, .125, .125, .75},
               Metadata m = source_metadata(true)) {
  return prepare_gaussian(r, MatrixKind::covariance, std::move(m), 1000, 1e-10,
                          irred::numerics::Arithmetic::longdouble_cpu_v1);
}
PredictiveMetadata metadata() {
  PredictiveMetadata m;
  m.ordered_parameter_ids = prior().ordered_parameter_ids;
  m.parameter_units = {"mag", "mag"};
  m.response_units = {"mag/mag", "mag/mag"};
  m.training_event_ids = {"synthetic-train-event0", "synthetic-train-event1"};
  m.future_event_ids = {"synthetic-future-event2", "synthetic-future-event3"};
  m.future_unit = "mag";
  m.future_covariance_unit = "mag^2";
  m.future_measure = "product d(mag)";
  m.response_identity = "offset+slope fixed t=[2,-1]";
  m.conditioning_identity =
      "exact two training values under chosen proper prior";
  m.conditional_noise_identity =
      "R future conditional on offset+slope, excludes marginalized beta";
  m.dependence_identity =
      "future noise independent of original beta prior and training noise";
  m.future_noise_independence_declared = true;
  m.noise_conditional_on_parameters_declared = true;
  return m;
}
GaussianPredictive owner() {
  auto p = posterior();
  auto r = noise();
  return GaussianPredictive::prepare(p, training, rows, r, response,
                                     metadata());
}
void rational_and_retention() {
  auto p = posterior();
  auto r = noise();
  auto m = metadata();
  auto got = GaussianPredictive::prepare(p, training, rows, r, response, m);
  need(got.status() == DensityStatus::finite &&
           p.status() == DensityStatus::finite &&
           r.status() == DensityStatus::finite,
       "borrowing preparation preserves both inputs");
  const W mean[]{-258.L / 668, 855.L / 668},
      cov[]{4444.L / 1336, -773.L / 1336, -773.L / 1336, 2618.L / 1336};
  for (unsigned i = 0; i < 2; ++i) {
    near(got.mean()[i], mean[i]);
    need(got.mean_absolute_error_estimates()[i] > 0,
         "mean diagnostic distinct from variance");
  }
  for (unsigned i = 0; i < 4; ++i) {
    near(got.covariance()[i], cov[i]);
    need(got.covariance_absolute_error_estimates()[i] > 0,
         "covariance arithmetic diagnostic");
  }
  need(got.metadata().future_event_ids == m.future_event_ids &&
           got.training_values()[0] == training[0] &&
           got.prior_identity() == p.prior().prior_identity &&
           got.shared_nuisance_ids()[0] == "offset",
       "exact lineage retained");
  const std::array<double, 2> y{.25, -.5};
  const W det = 11036863.L / 1784896, z0 = y[0] - mean[0], z1 = y[1] - mean[1],
          q = (cov[3] * z0 * z0 - 2 * cov[1] * z0 * z1 + cov[0] * z1 * z1) /
              det;
  const auto density = got.log_density(y, future);
  need(density.density.status == DensityStatus::finite,
       "normalized predictive density admitted");
  near(density.quadratic, q, 2e-11L);
  near(density.log_determinant, std::log(det), 2e-11L);
  near(density.normalization, 2 * std::log(2 * std::numbers::pi_v<W>), 2e-11L);
  near(density.density.log_value,
       -.5L * (q + std::log(det) + 2 * std::log(2 * std::numbers::pi_v<W>)),
       2e-11L);
  const auto *covariance_address = got.covariance().data(),
             *mean_address = got.mean().data();
  auto again = got.log_density(std::array<double, 2>{2, -1}, future);
  need(again.density.status == DensityStatus::finite &&
           got.covariance().data() == covariance_address &&
           got.mean().data() == mean_address,
       "repeated ordered whole-vector evaluation retains covariance and mean");
  p = GaussianPosterior{};
  r = Gaussian{};
  near(got.log_density(y, future).density.log_value, density.density.log_value,
       2e-11L);
  need(got.status() == DensityStatus::finite, "independent of input lifetimes");
  auto moved = std::move(got);
  need(got.status() == DensityStatus::invalid_input && got.mean().empty() &&
           got.covariance().empty() && got.training_values().empty() &&
           got.metadata().future_event_ids.empty(),
       "move invalidates all retained source views");
  auto *self = &moved;
  moved = std::move(*self);
  need(moved.status() == DensityStatus::finite, "self move preserves owner");
  auto other = owner();
  other = std::move(moved);
  need(other.status() == DensityStatus::finite &&
           moved.status() == DensityStatus::invalid_input,
       "move assignment coherent");
  const auto peak = GaussianPredictive::preparation_payload_bound(
      posterior(), noise(), metadata());
  need(peak && other.retained_payload_bound() &&
           other.evaluation_payload_bound() &&
           *peak > *other.retained_payload_bound(),
       "simultaneous payload bound includes borrowed owners/setup");
}
void limits_and_roles() {
  auto p = posterior();
  auto r = noise();
  for (unsigned mode = 0; mode < 14; ++mode) {
    auto m = metadata();
    auto policy = PredictivePolicy{};
    auto rr = source_metadata(true);
    auto a = response;
    if (mode == 0)
      m.future_noise_independence_declared = false;
    if (mode == 1)
      m.noise_conditional_on_parameters_declared = false;
    if (mode == 2)
      m.future_event_ids[0] = m.training_event_ids[1];
    if (mode == 3)
      rr.ordered_ids[0] = rows[0];
    if (mode == 4)
      std::swap(m.ordered_parameter_ids[0], m.ordered_parameter_ids[1]);
    if (mode == 5)
      m.response_units[0] = "unresolved";
    if (mode == 6)
      m.future_measure = "incorrect measure";
    if (mode == 7)
      rr.source_semantics = "posterior summary";
    if (mode == 8)
      policy.maximum_elements = 3;
    if (mode == 9)
      policy.maximum_payload_bytes = 0;
    if (mode == 10)
      policy.maximum_work_units = 0;
    if (mode == 11)
      policy.maximum_forward_sensitivity = 1e-30;
    if (mode == 12)
      a[0] = std::numeric_limits<double>::denorm_min();
    if (mode == 13)
      a[0] = std::numeric_limits<double>::quiet_NaN();
    auto rn = noise({.5, .125, .125, .75}, rr);
    auto bad = GaussianPredictive::prepare(p, training, rows, rn, a,
                                           std::move(m), policy);
    need(bad.status() != DensityStatus::finite && bad.mean().empty() &&
             bad.covariance().empty() && p.status() == DensityStatus::finite &&
             rn.status() == DensityStatus::finite,
         "invalid predictive support preserves borrowed inputs/no payload");
  }
  auto shifted =
      r.proper_offset(std::array<double, 2>{1, 1}, future, 0, .5,
                      "beta already marginalized", true, 1000, 1e-10);
  need(shifted.status() == DensityStatus::finite,
       "prior-bearing source control constructed");
  need(GaussianPredictive::prepare(p, training, rows, shifted, response,
                                   metadata())
               .status() != DensityStatus::finite,
       "same latent prior cannot be applied twice as future noise");
  for (auto matrix : {std::array<double, 4>{1, 1, 1, 1},
                      std::array<double, 4>{1, .1, .2, 1}}) {
    auto invalid = noise(matrix);
    need(GaussianPredictive::prepare(p, training, rows, invalid, response,
                                     metadata())
                 .status() != DensityStatus::finite,
         "invalid future covariance rejected");
  }
  auto got = owner();
  const std::array<double, 2> y{.25, -.5};
  for (unsigned mode = 0; mode < 8; ++mode) {
    auto policy = PredictivePolicy{};
    auto ys = y;
    auto ids = future;
    if (mode == 0)
      policy.maximum_elements = 0;
    if (mode == 1)
      policy.maximum_work_units = 0;
    if (mode == 2)
      policy.maximum_payload_bytes = 0;
    if (mode == 3)
      policy.maximum_forward_sensitivity = 1e-30;
    if (mode == 4)
      ys[0] = std::numeric_limits<double>::denorm_min();
    if (mode == 5)
      ys[0] = std::numeric_limits<double>::infinity();
    if (mode == 6)
      std::swap(ids[0], ids[1]);
    if (mode == 7)
      ys[0] = 1e308;
    need(got.log_density(ys, ids, policy).density.status !=
             DensityStatus::finite,
         "future input/quota/projection refusal");
  }
  const auto round = std::fegetround();
  need(std::fesetround(FE_DOWNWARD) == 0, "set hostile rounding");
  auto unsupported =
      GaussianPredictive::prepare(p, training, rows, r, response, metadata());
  const auto density = got.log_density(y, future);
  need(std::fesetround(round) == 0, "restore rounding");
  need(unsupported.status() == DensityStatus::unsupported_domain &&
           density.density.status == DensityStatus::unsupported_domain,
       "unsupported arithmetic refused for setup and evaluation");
  PredictivePolicy tight;
  tight.maximum_work_units = 127;
  need(got.log_density(y, future, tight).density.numerical_status ==
           S::work_limit,
       "evaluation work ceiling retained");
}
void mathematical_limits() {
  auto p = posterior();
  auto r = noise();
  const std::array<double, 4> zero{};
  auto null =
      GaussianPredictive::prepare(p, training, rows, r, zero, metadata());
  need(null.status() == DensityStatus::finite,
       "null response admitted with SPD future noise");
  for (unsigned i = 0; i < 2; ++i)
    near(null.mean()[i], 0);
  for (unsigned i = 0; i < 4; ++i)
    near(null.covariance()[i], r.covariance()[i]);
  auto no_training = posterior(zero);
  auto prior_predictive = GaussianPredictive::prepare(
      no_training, training, rows, r, response, metadata());
  need(prior_predictive.status() == DensityStatus::finite,
       "null training retains proper prior");
  near(prior_predictive.mean()[0], 0);
  near(prior_predictive.mean()[1], .75);
  near(prior_predictive.covariance()[0], 8.5);
  near(prior_predictive.covariance()[3], 2.75);
  auto diagonal = noise({.5, 0, 0, .75});
  const std::array<double, 4> repeated{1, 0, 1, 0};
  auto shared = GaussianPredictive::prepare(p, training, rows, diagonal,
                                            repeated, metadata());
  need(shared.status() == DensityStatus::finite,
       "repeated responses remain joint law");
  near(shared.covariance()[1], 180.L / 334);
  const std::array<double, 2> y{.25, -.5};
  const auto baseline = owner().log_density(y, future);
  auto pp = prior();
  pp.residual_unit = "scaled-mag";
  auto scaled_post = posterior({2, 0, 2, 2}, pp, {4, 1, 1, 8});
  auto md = source_metadata(true);
  md.measure = "product d(scaled-mag)";
  auto scaled_noise = noise({2, .5, .5, 3}, md);
  auto pm = metadata();
  pm.future_unit = "scaled-mag";
  pm.future_covariance_unit = "scaled-mag^2";
  pm.future_measure = md.measure;
  pm.response_units = {"scaled-mag/mag", "scaled-mag/mag"};
  auto scaled = GaussianPredictive::prepare(
      scaled_post, std::array<double, 2>{2.5, -1}, rows, scaled_noise,
      std::array<double, 4>{2, 4, 2, -2}, pm);
  need(scaled.status() == DensityStatus::finite,
       "explicit future-unit conversion admitted");
  near(scaled.log_density(std::array<double, 2>{.5, -1}, future)
           .density.log_value,
       baseline.density.log_value - 2 * std::log(2.L), 2e-11L);
  auto rank = posterior({1, 1, 2, 2});
  need(
      GaussianPredictive::prepare(rank, training, rows, r, response, metadata())
              .status() == DensityStatus::finite,
      "rank deficient training retains proper predictive law");
  auto one_md = source_metadata(true);
  one_md.ordered_ids.resize(1);
  auto one_noise = prepare_gaussian(
      std::array<double, 1>{.5}, MatrixKind::covariance, one_md, 1000, 1e-10,
      irred::numerics::Arithmetic::longdouble_cpu_v1);
  auto one_pm = metadata();
  one_pm.future_event_ids.resize(1);
  auto one = GaussianPredictive::prepare(p, training, rows, one_noise,
                                         std::array<double, 2>{1, 2}, one_pm);
  need(one.status() == DensityStatus::finite && one.mean().size() == 1,
       "one future coordinate admitted with proper two parameter posterior");
  near(one.mean()[0], -258.L / 668);
  near(one.covariance()[0], 4444.L / 1336);
  const auto one_density =
      one.log_density(std::array<double, 1>{.25}, one_md.ordered_ids);
  need(one_density.density.status == DensityStatus::finite,
       "single future density accepted");
  near(one_density.normalization, std::log(2 * std::numbers::pi_v<W>), 2e-11L);
  auto converted_prior = prior();
  converted_prior.parameter_units = {"changed-mag0", "changed-mag1"};
  const double scale[]{-2, 4};
  for (unsigned i = 0; i < 2; ++i) {
    converted_prior.mean[i] *= scale[i];
    for (unsigned j = 0; j < 2; ++j)
      converted_prior.covariance[i * 2 + j] *= scale[i] * scale[j];
  }
  auto converted = posterior({-.5, 0, -.5, .25}, converted_prior);
  auto converted_metadata = metadata();
  converted_metadata.parameter_units = converted_prior.parameter_units;
  converted_metadata.response_units = {"mag/changed-mag0", "mag/changed-mag1"};
  auto invariant = GaussianPredictive::prepare(
      converted, training, rows, r, std::array<double, 4>{-.5, .5, -.5, -.25},
      converted_metadata);
  need(invariant.status() == DensityStatus::finite,
       "explicit parameter unit conversion admitted");
  near(invariant.log_density(y, future).density.log_value,
       baseline.density.log_value, 2e-11L);
}
void allocations() {
  auto p = posterior();
  auto r = noise();
  size_t preparation_sites = 0;
  calls = 0;
  fail_on = 1;
  armed = true;
  GaussianPredictive empty;
  armed = false;
  need(calls.load() == 0 && empty.status() == DensityStatus::invalid_input,
       "empty predictive construction performs no allocation");
  {
    auto m = metadata();
    calls = 0;
    fail_on = 0;
    armed = true;
    auto measured =
        GaussianPredictive::prepare(p, training, rows, r, response, m);
    armed = false;
    preparation_sites = calls.load();
    need(measured.status() == DensityStatus::finite && preparation_sites > 0,
         "observable preparation allocation count");
  }
  for (size_t point = 1; point <= preparation_sites; ++point) {
    const long before = live.load();
    {
      auto m = metadata();
      calls = 0;
      fail_on = point;
      armed = true;
      auto rejected =
          GaussianPredictive::prepare(p, training, rows, r, response, m);
      armed = false;
      need(rejected.status() != DensityStatus::finite &&
               rejected.numerical_status() == S::work_limit &&
               rejected.mean().empty() && rejected.covariance().empty() &&
               m.future_event_ids.size() == 2,
           "every preparation allocation failure typed/no payload");
    }
    need(live.load() == before && p.status() == DensityStatus::finite &&
             r.status() == DensityStatus::finite,
         "failed preparation leak free and preserves borrowed inputs");
  }
  auto got = owner();
  const std::array<double, 2> y{.25, -.5};
  calls = 0;
  fail_on = 0;
  armed = true;
  const auto measured = got.log_density(y, future);
  armed = false;
  const auto evaluation_sites = calls.load();
  need(measured.density.status == DensityStatus::finite &&
           evaluation_sites > 0 && evaluation_sites < preparation_sites,
       "bounded evaluation scratch allocation sites");
  for (size_t point = 1; point <= evaluation_sites; ++point) {
    const long before = live.load();
    calls = 0;
    fail_on = point;
    armed = true;
    const auto rejected = got.log_density(y, future);
    armed = false;
    need(rejected.density.numerical_status == S::work_limit &&
             live.load() == before && got.status() == DensityStatus::finite,
         "every evaluation allocation failure preserves retained owner");
  }
  need(got.log_density(y, future).density.status == DensityStatus::finite,
       "owner usable after failure injection");
  std::cout << "allocation sites preparation=" << preparation_sites
            << " evaluation=" << evaluation_sites << '\n';
}
} // namespace
int main() {
  try {
    rational_and_retention();
    limits_and_roles();
    mathematical_limits();
    allocations();
    std::cout << "PASS " << checks
              << " proper Gaussian predictive owner controls\n";
  } catch (const std::exception &e) {
    armed = false;
    std::cerr << "FAIL after " << checks << ": " << e.what() << '\n';
    return 1;
  }
}
