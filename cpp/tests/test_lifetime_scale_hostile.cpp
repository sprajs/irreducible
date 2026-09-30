// Ordinary public owner moves, failed causes, and requested final-unit
// controls.
#include "irred/background.hpp"
#include "irred/bao.hpp"
#include "irred/statistics.hpp"
#include "irred/supernova.hpp"
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <new>
#include <stdexcept>
#include <type_traits>
// Reject every ordinary heap allocation during owner move operations. Inputs
// and nonempty assignment destinations are prepared before this switch.
static bool forbid_allocation = false;
// Keep allocator replacement boundaries visible to GCC's allocation-pair
// analysis; inlining free into a standard-container caller gives a false
// mismatched-new-delete warning despite the matching malloc-backed new above.
#if defined(__GNUC__) || defined(__clang__)
#define IRRED_TEST_NOINLINE __attribute__((noinline))
#else
#define IRRED_TEST_NOINLINE
#endif
IRRED_TEST_NOINLINE void *operator new(std::size_t n) {
  if (forbid_allocation)
    throw std::bad_alloc();
  if (void *p = std::malloc(n ? n : 1))
    return p;
  throw std::bad_alloc();
}
IRRED_TEST_NOINLINE void *operator new[](std::size_t n) {
  return ::operator new(n);
}
IRRED_TEST_NOINLINE void operator delete(void *p) noexcept { std::free(p); }
IRRED_TEST_NOINLINE void operator delete[](void *p) noexcept { std::free(p); }
IRRED_TEST_NOINLINE void operator delete(void *p, std::size_t) noexcept {
  std::free(p);
}
IRRED_TEST_NOINLINE void operator delete[](void *p, std::size_t) noexcept {
  std::free(p);
}
#undef IRRED_TEST_NOINLINE
using namespace irred;
namespace {
unsigned checks = 0;
void check(bool x, const char *why) {
  ++checks;
  if (!x)
    throw std::runtime_error(why);
}
statistics::Gaussian gaussian() {
  statistics::Metadata m;
  m.ordered_ids = {"a", "b"};
  m.measure = "magnitude";
  m.ordering_provenance = "synthetic explicit axes";
  return statistics::prepare_gaussian(std::array<double, 4>{4, 1, 1, 9},
                                      statistics::MatrixKind::covariance,
                                      std::move(m), 4, 1e-10);
}
observations::Prepared observation() {
  observations::Input x{};
  x.profile = observations::Profile::typed_magnitude_covariance;
  x.role = observations::Role::synthetic_control;
  x.unit = observations::Unit::magnitude;
  x.calibration = observations::Calibration::not_applicable;
  x.uncertainty = observations::Uncertainty::covariance;
  x.uncertainty_unit = observations::UncertaintyUnit::magnitude_squared;
  x.component = observations::Component::total;
  x.ordering_provenance = "synthetic explicit axes";
  x.table_sha256 = std::string(64, 'a');
  x.uncertainty_sha256 = std::string(64, 'b');
  x.measurement_ids = {"a", "b"};
  x.event_ids = x.measurement_ids;
  x.uncertainty_axis_ids = x.measurement_ids;
  x.values = {12, 15};
  x.uncertainty_matrix = {4, 1, 1, 9};
  x.missing = {0, 0};
  x.quality = {0, 0};
  x.source_selection = {1, 1};
  return observations::prepare(std::move(x), {2, 4, 4096});
}
supernova::Consumer sn_owner() {
  auto source = std::make_shared<const observations::Prepared>(observation());
  return supernova::prepare(
      {source,
       {0, 1},
       {"a", "b"},
       {{.1, {.1, cosmology::Convention::geometric_same_redshift}},
        {.5, {.5, cosmology::Convention::geometric_same_redshift}}}},
      {numerics::Arithmetic::binary64_legacy_v1, 2, 4, 1e-10, 1 << 20});
}
bao::PreparedDensity bao_owner() {
  bao::DensityInput x;
  x.queries = {{.1, bao::Observable::transverse_over_ruler},
               {.5, bao::Observable::hubble_over_ruler}};
  x.observed = {10, 20};
  x.covariance = {4, 1, 1, 9};
  x.ordered_ids = {"a", "b"};
  x.role = bao::RowRole::synthetic_control;
  x.covariance_unit = bao::CovarianceUnit::dimensionless_ratio_squared;
  x.table_identity = "synthetic";
  x.covariance_identity = "synthetic";
  x.ordering_provenance = "explicit axes";
  x.calibration_provenance = "synthetic free ruler; no calibration";
  x.dependence_provenance = "explicit synthetic full covariance";
  return bao::prepare_density(
      std::move(x),
      {2, 4, 4096, 1 << 20, 1e-10, numerics::Arithmetic::binary64_legacy_v1});
}
// Fixed analytic q=-1 screen: 1e-13 relative exceeds ordinary rounding
// in these few unit operations. Shared quantity scale is not an independent
// unit oracle; its separate constant fixtures remain authoritative.
void near(double a, long double b) {
  check(std::isfinite(a) &&
            std::abs((long double)a - b) <= 1e-13L * std::abs(b),
        "analytic final-unit value");
}
} // namespace
int main() {
  try {
    static_assert(std::is_nothrow_move_constructible_v<statistics::Gaussian>);
    static_assert(std::is_nothrow_move_assignable_v<statistics::Gaussian>);
    static_assert(std::is_copy_constructible_v<statistics::Gaussian>);
    const std::array<double, 2> r{2, -3}, response{1, 1};
    const std::array<size_t, 1> kept{0};
    auto obs = observation(), obs_destination = observation();
    auto obs_copy = obs;
    auto sn = sn_owner(), sn_destination = sn_owner();
    auto sn_copy = sn;
    auto ba = bao_owner(), ba_destination = bao_owner();
    auto ba_copy = ba;
    check(obs.status() == observations::Status::ok &&
              sn.status() == supernova::Status::ok &&
              ba.status() == statistics::DensityStatus::finite,
          "adjacent owners prepared");
    forbid_allocation = true;
    auto obs_moved = std::move(obs);
    obs_destination = std::move(obs_moved);
    auto sn_moved = std::move(sn);
    sn_destination = std::move(sn_moved);
    auto ba_moved = std::move(ba);
    ba_destination = std::move(ba_moved);
    forbid_allocation = false;
    check(obs.status() == observations::Status::invalid_shape &&
              obs.select(observations::Selection::all).status !=
                  observations::Status::ok,
          "moved observations cannot select empty success");
    check(sn.status() == supernova::Status::invalid_input &&
              ba.status() == statistics::DensityStatus::invalid_input,
          "moved consumer statuses invalid");
    auto *obs_same = &obs_destination;
    obs_destination = std::move(*obs_same);
    auto *sn_same = &sn_destination;
    sn_destination = std::move(*sn_same);
    auto *ba_same = &ba_destination;
    ba_destination = std::move(*ba_same);
    check(obs_destination.select(observations::Selection::all)
                      .source_indices.size() == 2 &&
              obs_copy.source().values == obs_destination.source().values,
          "observation destination/copy/selfmove");
    supernova::Policy sp;
    sp.arithmetic = numerics::Arithmetic::binary64_legacy_v1;
    sp.maximum_models = 1;
    sp.maximum_native_bytes = 1 << 20;
    sp.maximum_forward_sensitivity = 1e-10;
    sp.requested = 2;
    sp.background.integration =
        numerics::IntegrationPolicy{1e-14, 1e-14, 100000, 30};
    sp.background.maximum_callbacks = 100000;
    sp.background.maximum_native_bytes = 1 << 20;
    sp.background.maximum_queries = 2;
    sp.background.maximum_segment_visits = 100;
    const std::array<supernova::ModelPoint, 1> sm{
        supernova::ModelPoint(cosmology::ConstantQ{-1}, cosmology::FlatFLRW{},
                              supernova::NoMagnitudeEffect{})};
    auto sb = sn_destination.evaluate_batch(sm, sp),
         sc = sn_copy.evaluate_batch(sm, sp), sf = sn.evaluate_batch(sm, sp);
    check(sb.status == supernova::Status::ok && sb.slots.size() == 1 &&
              sb.slots[0].geometry.availability ==
                  cosmology::Availability::available &&
              sb.slots[0].geometric_shape.size() == 2 &&
              sb.slots[0].geometric_shape == sc.slots[0].geometric_shape,
          "SN destination/copy/selfmove retains predictions");
    check(sf.status != supernova::Status::ok && sf.slots.empty(),
          "moved SN no available empty result");
    bao::DensityPolicy bp;
    bp.arithmetic = sp.arithmetic;
    bp.maximum_models = 1;
    bp.maximum_native_bytes = 1 << 20;
    bp.maximum_forward_sensitivity = 1e-10;
    bp.requested = 7;
    bp.observables.maximum_queries = 2;
    bp.observables.maximum_native_bytes = 1 << 20;
    bp.observables.background = sp.background;
    const std::array<bao::ModelPoint, 1> bm{bao::ModelPoint(
        cosmology::ConstantQ{-1}, cosmology::FlatFLRW{}, bao::Ruler(10000))};
    auto bb = ba_destination.evaluate(bm, bp), bc = ba_copy.evaluate(bm, bp);
    check(bb.slots.size() == 1 && bb.slots[0].result &&
              bb.slots[0].result->density.status ==
                  statistics::DensityStatus::finite &&
              bc.slots.size() == 1 && bc.slots[0].result &&
              bb.slots[0].result->density.log_value ==
                  bc.slots[0].result->density.log_value,
          "BAO destination/copy/selfmove retains density");
    auto allocation_g = gaussian(), allocation_destination = gaussian();
    auto allocation_profile_g = gaussian(),
         allocation_profile_destination_g = gaussian();
    auto allocation_profile =
        std::move(allocation_profile_g)
            .prepare_offset_profile(
                response, allocation_profile_g.metadata().ordered_ids, 1e-10);
    auto allocation_profile_destination =
        std::move(allocation_profile_destination_g)
            .prepare_offset_profile(
                response,
                allocation_profile_destination_g.metadata().ordered_ids, 1e-10);
    forbid_allocation = true;
    auto allocation_moved = std::move(allocation_g);
    allocation_destination = std::move(allocation_moved);
    auto allocation_profile_moved = std::move(allocation_profile);
    allocation_profile_destination = std::move(allocation_profile_moved);
    forbid_allocation = false;
    check(allocation_destination.status() ==
                  statistics::DensityStatus::finite &&
              allocation_profile_destination.status() ==
                  statistics::DensityStatus::finite,
          "owner moves allocate no heap payload");
    auto g = gaussian();
    const auto initial = g.evaluate(r, g.metadata().ordered_ids, 1e-10);
    check(initial.density.status == statistics::DensityStatus::finite,
          "initial density finite");
    const auto expected = initial.density.log_value;
    auto copy = g;
    check(copy.status() == statistics::DensityStatus::finite,
          "copy behavior retained");
    auto moved = std::move(g);
    check(g.status() == statistics::DensityStatus::invalid_input,
          "moved Gaussian invalid");
    check(g.marginal(kept, 4, 1e-10).status() !=
              statistics::DensityStatus::finite,
          "moved marginal safe");
    check(g.conditional_zero_complement(kept, 4, 1e-10).status() !=
              statistics::DensityStatus::finite,
          "moved conditional safe");
    check(g.proper_offset(response, {}, 0, 1, "p", true, 4, 1e-10).status() !=
              statistics::DensityStatus::finite,
          "moved prior safe");
    check(
        moved.evaluate(r, moved.metadata().ordered_ids, 1e-10).density.status ==
                statistics::DensityStatus::finite &&
            moved.evaluate(r, moved.metadata().ordered_ids, 1e-10)
                    .density.log_value == expected,
        "move preserves density bits");
    auto *same_g = &moved;
    moved = std::move(*same_g);
    check(moved.status() == statistics::DensityStatus::finite,
          "self move usable");
    copy = std::move(moved);
    check(moved.status() == statistics::DensityStatus::invalid_input,
          "assignment source invalid");
    check(copy.evaluate(r, copy.metadata().ordered_ids, 1e-10).density.status ==
                  statistics::DensityStatus::finite &&
              copy.evaluate(r, copy.metadata().ordered_ids, 1e-10)
                      .density.log_value == expected,
          "assignment density preserved");
    auto p = std::move(copy).prepare_offset_profile(
        response, copy.metadata().ordered_ids, 1e-10);
    // Capture IDs from the resulting owner; the consuming factory invalidates
    // copy.
    check(p.status() == statistics::DensityStatus::finite,
          "profile prepared finite");
    auto pcopy = p;
    const auto initial_profile = p.evaluate(r, p.metadata().ordered_ids, 1e-10);
    check(initial_profile.status == statistics::DensityStatus::finite,
          "initial profile finite");
    const auto profile = initial_profile.quadratic;
    auto pmoved = std::move(p);
    check(p.status() == statistics::DensityStatus::invalid_input,
          "moved profile invalid");
    check(p.evaluate(r, {}, 1e-10).status != statistics::DensityStatus::finite,
          "moved profile safe");
    auto *same_p = &pmoved;
    pmoved = std::move(*same_p);
    const auto moved_profile =
        pmoved.evaluate(r, pmoved.metadata().ordered_ids, 1e-10);
    check(moved_profile.status == statistics::DensityStatus::finite &&
              moved_profile.quadratic == profile,
          "profile self move preserved");
    pcopy = std::move(pmoved);
    check(pcopy.evaluate(r, pcopy.metadata().ordered_ids, 1e-10).status ==
              statistics::DensityStatus::finite,
          "assigned profile usable");
    check(pmoved.status() == statistics::DensityStatus::invalid_input,
          "profile assignment invalidates source");
    auto f = numerics::cholesky(std::array<double, 4>{4, 1, 1, 9}, 2, 4);
    auto fcopy = f;
    auto fmoved = std::move(f);
    check(f.status() == numerics::Status::invalid_input && f.size() == 0,
          "factor source invalid and empty");
    check(numerics::solve(f, r, 1e-10).status ==
              numerics::Status::invalid_input,
          "moved solve rejected");
    auto *same_f = &fmoved;
    fmoved = std::move(*same_f);
    check(numerics::solve(fmoved, r, 1e-10).status == numerics::Status::ok,
          "factor self move");
    fcopy = std::move(fmoved);
    check(numerics::solve(fcopy, r, 1e-10).status == numerics::Status::ok,
          "assigned factor usable");
    check(fmoved.status() == numerics::Status::invalid_input,
          "factor assignment invalidates source");
    statistics::Metadata bad;
    bad.ordered_ids = {"a", "b"};
    bad.measure = "magnitude";
    bad.ordering_provenance = "synthetic";
    auto singular = statistics::prepare_gaussian(
        std::array<double, 4>{1, 1, 1, 1}, statistics::MatrixKind::covariance,
        std::move(bad), 4, 1e-10);
    auto failed_owner = std::move(singular);
    check(singular.status() == statistics::DensityStatus::invalid_input &&
              failed_owner.numerical_status() ==
                  numerics::Status::not_positive_definite,
          "failed owner move preserves cause");
    auto failed = failed_owner.profile_offset(
        r, response, failed_owner.metadata().ordered_ids, 1e-10);
    check(failed.numerical_status == numerics::Status::not_positive_definite,
          "profile retains SPD cause");
    using namespace cosmology;
    auto b = prepare(ConstantQ(-1), FlatFLRW{});
    EvaluationPolicy policy;
    policy.integration = numerics::IntegrationPolicy{1e-14, 1e-12, 10000, 30};
    policy.maximum_queries = 8;
    policy.maximum_callbacks = 10000;
    policy.maximum_native_bytes = 100000;
    std::array<Request, 3> requests{
        Request(.5, 4, {}, PhysicalScale(1e-285)),
        Request(1e-300, 8,
                Observer{1e-300, Convention::geometric_same_redshift},
                PhysicalScale(1e-285)),
        Request(.5, 8, Observer{.5, Convention::geometric_same_redshift},
                PhysicalScale(1e-285))};
    const auto result = b.evaluate(requests, policy);
    check(result.slots[0].clock.value &&
              result.slots[0].clock.value->lookback_seconds.value,
          "clock ignores metre overflow");
    near(*result.slots[0].clock.value->lookback_seconds.value,
         std::log(1.5L) * (megaparsec_in_metres_wide() / 1000) /
             (double)1e-285);
    check(result.slots[1].physical.value.has_value(),
          "finite Mpc outputs survive metre overflow");
    const long double scale = 299792458.L / 1000 / (double)1e-285,
                      dm = scale * (double)1e-300;
    near(result.slots[1].physical.value->radial_mpc, dm);
    near(result.slots[1].physical.value->volume_mpc3_per_sr_per_redshift,
         scale * dm * dm);
    check(!result.slots[2].physical.value &&
              result.slots[2].physical.availability == Availability::failed,
          "actual final volume overflow still fails");
    std::cout << "Lifetime and final-unit controls PASS " << checks << '\n';
  } catch (const std::exception &e) {
    std::cerr << "FAIL " << checks << ": " << e.what() << '\n';
    return 1;
  }
}
