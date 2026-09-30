#include "../src/payload_accounting.hpp"
#include "irred/bao.hpp"
#include "irred/numerics.hpp"
#include "irred/statistics.hpp"
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace irred;
int main() {
  unsigned count = 0;
  auto check = [&](bool good) {
    ++count;
    if (!good)
      throw std::runtime_error("payload accounting check " +
                               std::to_string(count));
  };
  numerics::Factorization empty;
  check(empty.retained_payload_bound() == sizeof(empty));
  const double C[]{4, 1, 1, 9};
  for (auto arithmetic : {numerics::Arithmetic::binary64_legacy_v1,
                          numerics::Arithmetic::longdouble_cpu_v1}) {
    auto factor = numerics::cholesky(C, 2, 4, arithmetic);
    check(factor.status() == numerics::Status::ok);
    const size_t matrices =
        arithmetic == numerics::Arithmetic::binary64_legacy_v1
            ? 8 * sizeof(double)
            : 4 * (sizeof(double) + sizeof(long double));
    check(factor.retained_payload_bound() == sizeof(factor) + matrices);
    auto moved = std::move(factor);
    check(moved.retained_payload_bound() == sizeof(moved) + matrices);
  }
  detail::PayloadAccounting overflow(1);
  overflow.add(SIZE_MAX, 2);
  check(!overflow.result());
  detail::PayloadAccounting exact(7);
  exact.add(3, 5);
  check(exact.result() == 22);
  bao::DensityInput source;
  source.queries.reserve(7);
  source.observed.reserve(3);
  source.covariance.reserve(19);
  source.ordered_ids.reserve(3);
  source.ordered_ids.emplace_back("a");
  source.ordered_ids[0].reserve(103);
  source.table_identity.reserve(117);
  size_t expected =
      source.queries.capacity() * sizeof(bao::Query) +
      (source.observed.capacity() + source.covariance.capacity()) *
          sizeof(double) +
      source.ordered_ids.capacity() * sizeof(std::string);
  for (const auto &id : source.ordered_ids)
    expected += id.capacity() + 1;
  for (const auto *s :
       {&source.table_identity, &source.covariance_identity,
        &source.ordering_provenance, &source.calibration_provenance,
        &source.dependence_provenance})
    expected += s->capacity() + 1;
  check(bao::retained_source_payload_bound(source) == expected);
  statistics::Metadata metadata;
  metadata.ordered_ids = {"a", "b"};
  metadata.measure = "synthetic";
  metadata.ordering_provenance = "explicit synthetic ordered IDs";
  auto gaussian = statistics::prepare_gaussian(
      C, statistics::MatrixKind::covariance, metadata, 4, 1e-8);
  check(gaussian.status() == statistics::DensityStatus::finite);
  auto before = gaussian.retained_payload_bound();
  check(before && *before > sizeof(gaussian));
  double r[]{1, 2};
  auto value = gaussian.evaluate(r, metadata.ordered_ids, 1e-8);
  check(value.density.status == statistics::DensityStatus::finite &&
        value.quadratic == .6);
  check(gaussian.retained_payload_bound() == before);
  const auto covariance_peak = statistics::gaussian_preparation_payload_bound(
      2, statistics::MatrixKind::covariance,
      numerics::Arithmetic::binary64_legacy_v1, metadata);
  const auto precision_peak = statistics::gaussian_preparation_payload_bound(
      2, statistics::MatrixKind::precision,
      numerics::Arithmetic::longdouble_cpu_v1, metadata);
  check(covariance_peak && *covariance_peak > *before);
  check(precision_peak && *precision_peak > *covariance_peak);
  check(!statistics::gaussian_preparation_payload_bound(
      SIZE_MAX, statistics::MatrixKind::precision,
      numerics::Arithmetic::longdouble_cpu_v1, metadata));
  check(!statistics::gaussian_preparation_payload_bound(
      2, statistics::MatrixKind::covariance,
      static_cast<numerics::Arithmetic>(UINT32_MAX), metadata));
  auto large_metadata = metadata;
  large_metadata.ordered_ids[0].reserve(200000);
  auto larger = statistics::gaussian_preparation_payload_bound(
      2, statistics::MatrixKind::covariance,
      numerics::Arithmetic::binary64_legacy_v1, large_metadata);
  check(larger && *larger > *covariance_peak + 200000);
  const auto proper_peak = gaussian.proper_offset_payload_bound("explicit prior");
  check(proper_peak && *proper_peak > *before + *covariance_peak);
  auto proper = gaussian.proper_offset(
      r, metadata.ordered_ids, 0, 1, "explicit prior", true, 4, 1e-8);
  check(proper.status() == statistics::DensityStatus::finite);
  auto proper_retained = proper.retained_payload_bound();
  check(proper_retained && *proper_peak > *before + *proper_retained);
  check(gaussian.evaluate(r, metadata.ordered_ids, 1e-8).quadratic == value.quadratic);
  std::cout << "Payload accounting " << count << " PASS\n";
}
