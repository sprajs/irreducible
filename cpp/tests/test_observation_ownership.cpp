#include "../src/observation_internal.hpp"
#include "irred/abi.h"
#include "irred/observations.hpp"
#include <iostream>
#include <stdexcept>
using namespace irred::observations;
namespace {
void check(bool x) {
  if (!x)
    throw std::runtime_error("observation ownership contract");
}
Input source() {
  Input s{};
  s.profile = Profile::typed_magnitude_covariance;
  s.role = Role::released_fitted_summary;
  s.unit = Unit::magnitude;
  s.calibration = Calibration::unknown;
  s.uncertainty = Uncertainty::covariance;
  s.uncertainty_unit = UncertaintyUnit::magnitude_squared;
  s.table_sha256 = std::string(64, 'a');
  s.uncertainty_sha256 = std::string(64, 'b');
  s.ordering_provenance =
      "declared axes, no independently verified release claim";
  s.measurement_ids = {"row-a", "row-b"};
  s.event_ids = {"event-a", "event-b"};
  s.uncertainty_axis_ids = s.measurement_ids;
  s.values = {21, 23};
  s.uncertainty_matrix = {2, .25, .25, 3};
  s.missing = {0, 0};
  s.quality = {0, 0};
  return s;
}
cosmo_bytes bytes(const std::string &x) {
  return {reinterpret_cast<const uint8_t *>(x.data()), x.size()};
}
cosmo_f64_buffer doubles(const std::vector<double> &x) {
  return {sizeof(cosmo_f64_buffer), COSMO_ABI_VERSION, 2, 0, x.data(), x.size(),
          x.size() * sizeof(double)};
}
} // namespace
int main() {
  try {
    auto s = source();
    auto p = prepare(s, {8, 64, 4096});
    check(p.status() == Status::ok);
    check(p.source().profile == Profile::typed_magnitude_covariance);
    check(p.source().calibration == Calibration::unknown);
    check(p.select(Selection::pantheon_zhd_gt_001).status ==
          Status::incompatible_semantics);
    for (auto role : {Role::observed_measurement, Role::released_fitted_summary,
                      Role::synthetic_control}) {
      auto x = s;
      x.role = role;
      check(prepare(std::move(x), {8, 64, 4096}).status() == Status::ok);
    }
    auto x = s;
    x.role = Role::posterior_summary;
    check(prepare(std::move(x), {8, 64, 4096}).status() ==
          Status::incompatible_semantics);
    x = s;
    x.uncertainty_axis_ids = {"row-b", "row-a"};
    check(prepare(std::move(x), {8, 64, 4096}).status() ==
          Status::invalid_shape);
    x = s;
    x.zhd = {.1, .2};
    x.zhd_missing = {0, 1};
    check(prepare(x, {8, 64, 4096}).status() == Status::ok);
    x.zhd_missing = {0};
    check(prepare(std::move(x), {8, 64, 4096}).status() ==
          Status::incompatible_semantics);
    x = s;
    x.uncertainty_unit = UncertaintyUnit::metre_squared;
    check(prepare(std::move(x), {8, 64, 4096}).status() ==
          Status::incompatible_semantics);
    x = s;
    x.uncertainty = Uncertainty::precision;
    x.uncertainty_unit = UncertaintyUnit::inverse_magnitude_squared;
    check(prepare(std::move(x), {8, 64, 4096}).status() ==
          Status::incompatible_semantics);
    cosmo_observation_descriptor d{};
    d.struct_size = sizeof(d);
    d.abi_version = COSMO_ABI_VERSION;
    d.profile = static_cast<uint32_t>(s.profile);
    d.role = static_cast<uint32_t>(s.role);
    d.unit = static_cast<uint32_t>(s.unit);
    d.calibration = static_cast<uint32_t>(s.calibration);
    d.uncertainty = static_cast<uint32_t>(s.uncertainty);
    d.uncertainty_unit = static_cast<uint32_t>(s.uncertainty_unit);
    d.table_sha256 = bytes(s.table_sha256);
    d.uncertainty_sha256 = bytes(s.uncertainty_sha256);
    d.ordering_provenance = bytes(s.ordering_provenance);
    std::vector<cosmo_bytes> ids, events;
    for (const auto &id : s.measurement_ids)
      ids.push_back(bytes(id));
    for (const auto &id : s.event_ids)
      events.push_back(bytes(id));
    d.measurement_ids = {ids.data(), ids.size(),
                         ids.size() * sizeof(cosmo_bytes)};
    d.event_ids = {events.data(), events.size(),
                   events.size() * sizeof(cosmo_bytes)};
    d.uncertainty_axis_ids = d.measurement_ids;
    d.values = doubles(s.values);
    d.zhd = doubles(s.zhd);
    d.zcmb = doubles(s.zcmb);
    d.zhel = doubles(s.zhel);
    d.uncertainty_matrix = doubles(s.uncertainty_matrix);
    auto mask = [](const auto &v) {
      return cosmo_u8_buffer{v.data(), v.size(), v.size()};
    };
    d.missing = mask(s.missing);
    d.zhd_missing = mask(s.zhd_missing);
    d.zcmb_missing = mask(s.zcmb_missing);
    d.zhel_missing = mask(s.zhel_missing);
    d.source_selection = mask(s.source_selection);
    d.quality = {s.quality.data(), s.quality.size(),
                 s.quality.size() * sizeof(uint64_t)};
    cosmo_observation_policy policy{8, 64, 4096};
    cosmo_prepared *raw = nullptr;
    uint32_t status = UINT32_MAX;
    check(cosmo_prepare_observations(&d, &policy, &raw, &status) == COSMO_OK &&
          status == 0 && raw);
    auto retained = shared_native_observations(raw);
    check(retained.get() == native_observations(raw));
    check(retained.use_count() == 2);
    const auto *address = retained.get();
    check(cosmo_observation_destroy(raw) == COSMO_OK);
    s.values[0] = 999;
    s.uncertainty_matrix.clear();
    check(retained.get() == address && retained.use_count() == 1);
    check(retained->source().values[0] == 21 &&
          retained->source().uncertainty_matrix.size() == 4);
    check(retained->source().measurement_ids[1] == "row-b");
    std::cout << "Shared immutable observation ownership and generic admission "
                 "passed\n";
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
