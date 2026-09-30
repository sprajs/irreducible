#include "../src/observation_internal.hpp"
#include "irred/observations.hpp"
#include <cstdlib>
#include <iostream>
#include <new>
#include <stdexcept>
#include <string>
using namespace irred::observations;
namespace {
long countdown = -1, live = 0;
unsigned checks = 0;
void check(bool x, const char *why) {
  ++checks;
  if (!x)
    throw std::runtime_error(why);
}
cosmo_bytes bytes(const std::string &x) {
  return {reinterpret_cast<const uint8_t *>(x.data()), x.size()};
}
cosmo_f64_buffer f64(const double *x, size_t n) {
  return {sizeof(cosmo_f64_buffer), COSMO_ABI_VERSION, 2, 0, x, n,
          n * sizeof(double)};
}
struct Fixture {
  std::string
      a = "x",
      b = std::string(80, 'b'), hash = std::string(64, 'a'),
      covhash = std::string(64, 'c'),
      order =
          "caller declares original covariance axes; not release verification";
  cosmo_bytes ids[2]{bytes(a), bytes(b)}, events[2]{bytes(a), bytes(b)};
  double y[2]{20, 21}, cov[4]{2, .5, .5, 3}, z[2]{.001, .3};
  uint8_t mask[2]{0, 0}, zm[2]{0, 1};
  uint64_t quality[2]{0, 0};
  cosmo_observation_policy policy{2, 4, 4096};
  cosmo_observation_descriptor d{};
  Fixture() {
    d.struct_size = sizeof(d);
    d.abi_version = COSMO_ABI_VERSION;
    d.profile = static_cast<uint32_t>(Profile::typed_magnitude_covariance);
    d.role = static_cast<uint32_t>(Role::observed_measurement);
    d.unit = static_cast<uint32_t>(Unit::magnitude);
    d.calibration = static_cast<uint32_t>(Calibration::unknown);
    d.uncertainty = static_cast<uint32_t>(Uncertainty::covariance);
    d.uncertainty_unit =
        static_cast<uint32_t>(UncertaintyUnit::magnitude_squared);
    d.table_sha256 = bytes(hash);
    d.uncertainty_sha256 = bytes(covhash);
    d.ordering_provenance = bytes(order);
    d.measurement_ids = {ids, 2, sizeof(ids)};
    d.event_ids = {events, 2, sizeof(events)};
    d.uncertainty_axis_ids = d.measurement_ids;
    d.values = f64(y, 2);
    d.uncertainty_matrix = f64(cov, 4);
    d.zhd = f64(z, 2);
    d.zcmb = d.zhel = f64(nullptr, 0);
    d.missing = {mask, 2, 2};
    d.zhd_missing = {zm, 2, 2};
    d.quality = {quality, 2, sizeof(quality)};
  }
};
} // namespace
void *operator new(size_t n) {
  if (countdown == 0)
    throw std::bad_alloc();
  if (countdown > 0)
    --countdown;
  void *p = std::malloc(n ? n : 1);
  if (!p)
    throw std::bad_alloc();
  ++live;
  return p;
}
void operator delete(void *p) noexcept {
  if (p) {
    --live;
    std::free(p);
  }
}
void operator delete(void *p, size_t) noexcept { ::operator delete(p); }
void *operator new[](size_t n) { return ::operator new(n); }
void operator delete[](void *p) noexcept { ::operator delete(p); }
void operator delete[](void *p, size_t) noexcept { ::operator delete(p); }
int main() {
  try {
    Fixture f;
    cosmo_prepared *raw = nullptr;
    uint32_t semantic = 99;
    for (auto role : {Role::observed_measurement, Role::released_fitted_summary,
                      Role::synthetic_control}) {
      f.d.role = static_cast<uint32_t>(role);
      check(cosmo_prepare_observations(&f.d, &f.policy, &raw, &semantic) ==
                    COSMO_OK &&
                raw && semantic == 0,
            "explicit role admission");
      auto owner = shared_native_observations(raw);
      auto second = shared_native_observations(raw);
      const auto *address = native_observations(raw);
      check(owner.get() == address && second.get() == address,
            "shared identical object");
      check(owner->source().role == role &&
                owner->source().calibration == Calibration::unknown &&
                owner->source().calibration_provenance.empty(),
            "no role or calibration upgrade");
      check(owner->source().zhd[0] == .001 &&
                owner->source().zhd_missing[1] == 1,
            "no implicit cut or mask erasure");
      check(cosmo_observation_destroy(raw) == COSMO_OK, "release acquisition");
      raw = nullptr;
      owner.reset();
      check(second.get() == address &&
                second->source().measurement_ids[1] == f.b,
            "final retained owner survives");
    }
    auto reject = [&](uint32_t expected) {
      raw = nullptr;
      semantic = 99;
      check(cosmo_prepare_observations(&f.d, &f.policy, &raw, &semantic) ==
                    COSMO_OK &&
                !raw && semantic == expected,
            "semantic rejection without owner");
    };
    f.d.role = static_cast<uint32_t>(Role::posterior_summary);
    reject(static_cast<uint32_t>(Status::incompatible_semantics));
    f.d.role = static_cast<uint32_t>(Role::synthetic_control);
    f.d.uncertainty = static_cast<uint32_t>(Uncertainty::precision);
    f.d.uncertainty_unit =
        static_cast<uint32_t>(UncertaintyUnit::inverse_magnitude_squared);
    reject(static_cast<uint32_t>(Status::incompatible_semantics));
    f.d.uncertainty = static_cast<uint32_t>(Uncertainty::covariance);
    f.d.uncertainty_unit =
        static_cast<uint32_t>(UncertaintyUnit::magnitude_squared);
    cosmo_bytes reversed[2]{f.ids[1], f.ids[0]};
    f.d.uncertainty_axis_ids = {reversed, 2, sizeof(reversed)};
    reject(static_cast<uint32_t>(Status::invalid_shape));
    f.d.uncertainty_axis_ids = f.d.measurement_ids;
    f.d.zhd_missing = {f.zm, 1, 1};
    reject(static_cast<uint32_t>(Status::incompatible_semantics));
    f.d.zhd_missing = {f.zm, 2, 2};
    bool success = false;
    unsigned failures = 0;
    for (long i = 0; i < 100; ++i) {
      auto before = live;
      countdown = i;
      auto code = cosmo_prepare_observations(&f.d, &f.policy, &raw, &semantic);
      countdown = -1;
      if (code == COSMO_OK) {
        check(raw && semantic == 0, "successful sweep terminal");
        auto keep = shared_native_observations(raw);
        check(cosmo_observation_destroy(raw) == COSMO_OK,
              "sweep release acquisition");
        raw = nullptr;
        f.y[0] = 999;
        f.cov[0] = 999;
        f.a = "mutated";
        check(keep->source().values[0] == 20 &&
                  keep->source().uncertainty_matrix[0] == 2 &&
                  keep->source().measurement_ids[0] == "x",
              "deep immutable input retention");
        keep.reset();
        check(live == before, "last shared owner recovers allocation quota");
        success = true;
        break;
      }
      check(code == COSMO_ALLOCATION_FAILURE && !raw,
            "allocation failure reset");
      check(live == before, "allocation failure cleanup");
      ++failures;
    }
    check(success && failures > 10,
          "all allocation sites swept including shared control block");
    std::cout << "Observation ownership hostile " << checks << " checks, "
              << failures << " allocation failures passed\n";
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
