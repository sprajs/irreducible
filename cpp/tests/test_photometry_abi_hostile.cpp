// Transport equivalence and lifetime/resource controls, not independent
// physics.
#include "irred/abi.h"
#include "irred/photometry.hpp"
#include <array>
#include <bit>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <new>
#include <stdexcept>
static long fail_at = -1, allocations = 0, live = 0;
#if defined(__GNUC__) || defined(__clang__)
#define NOINLINE __attribute__((noinline))
#else
#define NOINLINE
#endif
NOINLINE void *operator new(std::size_t n) {
  if (fail_at >= 0 && allocations++ == fail_at)
    throw std::bad_alloc();
  if (void *p = std::malloc(n ? n : 1)) {
    ++live;
    return p;
  }
  throw std::bad_alloc();
}
NOINLINE void *operator new[](std::size_t n) { return ::operator new(n); }
NOINLINE void operator delete(void *p) noexcept {
  if (p) {
    --live;
    std::free(p);
  }
}
NOINLINE void operator delete[](void *p) noexcept { ::operator delete(p); }
NOINLINE void operator delete(void *p, std::size_t) noexcept {
  ::operator delete(p);
}
NOINLINE void operator delete[](void *p, std::size_t) noexcept {
  ::operator delete(p);
}
namespace {
unsigned checks = 0;
void check(bool x, const char *w) {
  ++checks;
  if (!x)
    throw std::runtime_error(w);
}
namespace p = irred::photometry;
namespace n = irred::numerics;
irred_photometry_input input() {
  return {sizeof(irred_photometry_input),
          IRRED_ABI_VERSION,
          1e15,
          4e-7,
          8e-7,
          3e-7,
          1.5e-6,
          1e10,
          .5,
          2,
          .75,
          120};
}
p::Input native(const irred_photometry_input &x) {
  return {x.luminosity_watt_per_metre,
          x.rest_lower_metre,
          x.rest_upper_metre,
          x.observed_lower_metre,
          x.observed_upper_metre,
          x.luminosity_distance_metre,
          x.redshift,
          x.collecting_area_square_metre,
          x.optical_transmission,
          x.observer_exposure_second};
}
auto bits(const irred_photometry_input &x) {
  return std::array{x.luminosity_watt_per_metre,
                    x.rest_lower_metre,
                    x.rest_upper_metre,
                    x.observed_lower_metre,
                    x.observed_upper_metre,
                    x.luminosity_distance_metre,
                    x.redshift,
                    x.collecting_area_square_metre,
                    x.optical_transmission,
                    x.observer_exposure_second};
}
irred_photometry_policy policy() {
  return {
      sizeof(irred_photometry_policy), IRRED_ABI_VERSION, 7, 0, 64, 1u << 20};
}
irred_photometry_batch batch(std::span<const irred_photometry_input> r) {
  return {sizeof(irred_photometry_batch), IRRED_ABI_VERSION, r.data(), r.size(),
          r.size_bytes()};
}
irred_photometry_view view(irred_photometry_result *r) {
  irred_photometry_view v{};
  check(irred_photometry_result_view(r, &v) == IRRED_OK, "owned view");
  return v;
}
void scalar(const irred_photometry_scalar &a, const p::Outcome &b) {
  check(a.availability == static_cast<unsigned>(b.availability),
        "availability exact native");
  check(a.numerical_status == static_cast<unsigned>(b.numerical_status),
        "cause exact native");
  if (b.value)
    check(std::bit_cast<std::uint64_t>(a.value) ==
              std::bit_cast<std::uint64_t>(*b.value),
          "value exact native bits");
  else
    check(a.value == 0, "absent initialized wire storage");
}
void reject(irred_photometry_batch b, irred_photometry_policy pol,
            unsigned expected = IRRED_INVALID_INPUT) {
  irred_photometry_result *r = reinterpret_cast<irred_photometry_result *>(1);
  check(irred_photometry_evaluate(&b, &pol, &r) == expected,
        "structural disposition");
  check(r == nullptr, "structural output reset");
}
} // namespace
int main() {
  try {
    std::array rows{input(), input(), input()};
    rows[1].redshift = -1;
    rows[2].optical_transmission = -0.0;
    auto b = batch(rows);
    auto pol = policy();
    irred_photometry_result *r = nullptr;
    check(irred_photometry_evaluate(&b, &pol, &r) == IRRED_OK && r,
          "mixed owned result");
    auto v = view(r);
    check(v.numerical_status == 0 && v.row_count == 3, "mixed outer completed");
    std::array ns{native(rows[0]), native(rows[1]), native(rows[2])};
    auto nb = p::evaluate(ns, {7, 64, 1u << 20});
    for (unsigned i = 0; i < 3; ++i) {
      check(v.rows[i].admission_status ==
                static_cast<unsigned>(nb.rows[i].admission_status),
            "admission exact native");
      const auto a = bits(v.rows[i].source), c = bits(rows[i]);
      for (unsigned j = 0; j < 10; ++j)
        check(std::bit_cast<std::uint64_t>(a[j]) ==
                  std::bit_cast<std::uint64_t>(c[j]),
              "owned source exact bits/order");
      scalar(v.rows[i].incident_band_flux,
             nb.rows[i].flux_watt_per_square_metre);
      scalar(v.rows[i].collected_energy, nb.rows[i].energy_joule);
      scalar(v.rows[i].expected_transmitted_photons,
             nb.rows[i].expected_photons);
    }
    rows[0].redshift = 9;
    check(v.rows[0].source.redshift == .5,
          "caller mutation leaves owned source");
    check(v.model_id.length && v.constants_id.length &&
              v.arithmetic_id.length && v.propagation_id.length,
          "static scientific metadata present");
    check(irred_photometry_result_destroy(r) == IRRED_OK, "destroy result");
    check(irred_photometry_result_destroy(nullptr) == IRRED_OK, "destroy null");
    rows[0] = input();
    b = batch(rows);
    pol.requested_outputs = 2;
    r = nullptr;
    check(irred_photometry_evaluate(&b, &pol, &r) == IRRED_OK,
          "requested energy batch");
    v = view(r);
    check(v.rows[0].incident_band_flux.availability == 0 &&
              v.rows[0].expected_transmitted_photons.availability == 0,
          "unrequested group omitted");
    check(v.rows[0].collected_energy.availability == 1,
          "requested energy available");
    irred_photometry_result_destroy(r);
    pol = policy();
    auto malformed = b;
    malformed.abi_version++;
    reject(malformed, pol, IRRED_ABI_MISMATCH);
    malformed = b;
    malformed.struct_size--;
    reject(malformed, pol);
    auto wrong = pol;
    wrong.reserved = 1;
    reject(b, wrong);
    wrong = pol;
    wrong.requested_outputs = 8;
    reject(b, wrong);
    wrong = pol;
    wrong.requested_outputs = 0;
    reject(b, wrong);
    malformed = b;
    malformed.byte_length--;
    reject(malformed, pol);
    malformed = b;
    malformed.length = std::numeric_limits<std::uint64_t>::max();
    malformed.data = reinterpret_cast<const irred_photometry_input *>(1);
    reject(malformed, pol);
    malformed = b;
    malformed.data = nullptr;
    reject(malformed, pol);
    rows[0].abi_version++;
    reject(b, pol, IRRED_ABI_MISMATCH);
    rows[0] = input();
    alignas(irred_photometry_input)
        std::array<std::byte, sizeof(irred_photometry_input) + 1>
            misaligned{};
    malformed = b;
    malformed.data =
        reinterpret_cast<const irred_photometry_input *>(misaligned.data() + 1);
    reject(malformed, pol);
    pol.maximum_rows = 0;
    r = nullptr;
    check(irred_photometry_evaluate(&b, &pol, &r) == IRRED_OK && r,
          "configured count owned diagnostic");
    v = view(r);
    check(v.numerical_status == static_cast<unsigned>(n::Status::work_limit) &&
              v.row_count == 0,
          "count failure no fake rows");
    irred_photometry_result_destroy(r);
    pol = policy();
    pol.maximum_native_bytes = 0;
    r = nullptr;
    check(irred_photometry_evaluate(&b, &pol, &r) == IRRED_OK && r,
          "zero byte diagnostic owner permitted");
    v = view(r);
    check(v.numerical_status == static_cast<unsigned>(n::Status::work_limit) &&
              v.row_count == 0,
          "byte failure empty diagnostic");
    irred_photometry_result_destroy(r);
    // Find the actual admitted threshold independently of private owner sizeof.
    // Same descriptors and three rows; scalar export choices retain the same
    // wire.
    std::uint64_t lo = 0, hi = 1u << 20;
    while (lo < hi) {
      auto q = policy();
      q.maximum_native_bytes = (lo + hi) / 2;
      r = nullptr;
      check(irred_photometry_evaluate(&b, &q, &r) == IRRED_OK && r,
            "threshold probe transport");
      auto vv = view(r);
      bool ok = vv.numerical_status == 0;
      irred_photometry_result_destroy(r);
      if (ok)
        hi = q.maximum_native_bytes;
      else
        lo = q.maximum_native_bytes + 1;
    }
    auto q = policy();
    q.maximum_native_bytes = lo;
    r = nullptr;
    check(irred_photometry_evaluate(&b, &q, &r) == IRRED_OK && r,
          "exact combined quota admitted");
    check(view(r).row_count == 3, "combined quota actual rows");
    irred_photometry_result_destroy(r);
    q.maximum_native_bytes = lo - 1;
    r = nullptr;
    check(irred_photometry_evaluate(&b, &q, &r) == IRRED_OK && r,
          "one under combined diagnostic");
    check(view(r).row_count == 0, "one under combined no rows");
    irred_photometry_result_destroy(r);
    pol = policy();
    const long baseline = live;
    allocations = 0;
    fail_at = 1000000;
    r = nullptr;
    auto transport = irred_photometry_evaluate(&b, &pol, &r);
    const long sites = allocations;
    fail_at = -1;
    check(transport == IRRED_OK && r, "allocation baseline");
    irred_photometry_result_destroy(r);
    check(live == baseline, "baseline allocations recovered");
    for (long i = 0; i < sites; ++i) {
      allocations = 0;
      fail_at = i;
      r = nullptr;
      transport = irred_photometry_evaluate(&b, &pol, &r);
      fail_at = -1;
      check(transport == IRRED_ALLOCATION_FAILURE && r == nullptr,
            "each allocation failure translated/null");
      check(live == baseline, "each allocation failure cleanup");
    }
    check(sites > 0, "allocation sweep executed");
    b = batch({});
    r = nullptr;
    check(irred_photometry_evaluate(&b, &pol, &r) == IRRED_OK && r,
          "empty admitted owner");
    check(view(r).row_count == 0 && view(r).numerical_status == 0,
          "empty zero output status");
    irred_photometry_result_destroy(r);
    std::cout << "PASS " << checks
              << " photometry ABI controls allocation_sites=" << sites
              << " combined_threshold=" << lo << '\n';
  } catch (const std::exception &e) {
    fail_at = -1;
    std::cerr << "FAIL check " << checks << ": " << e.what() << '\n';
    return 1;
  }
}
