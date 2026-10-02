// Coarse transport/lifetime controls; independent science is in sampled peer
// tests.
#include "irred/abi.h"
#include "irred/sampled_photometry.hpp"
#include <array>
#include <bit>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <new>
#include <stdexcept>
static long fail_at = -1, allocations = 0, live = 0;
static bool nonallocation_exception = false;
#if defined(__GNUC__) || defined(__clang__)
#define NOINLINE __attribute__((noinline))
#else
#define NOINLINE
#endif
NOINLINE void *operator new(std::size_t n) {
  if (fail_at >= 0 && allocations++ == fail_at) {
    if (nonallocation_exception)
      throw 1;
    throw std::bad_alloc();
  }
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
namespace p = irred::photometry;
namespace n = irred::numerics;
unsigned checks = 0;
void check(bool b, const char *what) {
  ++checks;
  if (!b)
    throw std::runtime_error(what);
}
irred_sampled_photometry_curve curve(const char *id, std::span<const double> x,
                                     std::span<const double> y) {
  return {sizeof(irred_sampled_photometry_curve),
          IRRED_ABI_VERSION,
          {reinterpret_cast<const std::uint8_t *>(id),
           std::char_traits<char>::length(id)},
          x.data(),
          y.data(),
          x.size(),
          x.size_bytes()};
}
irred_sampled_photometry_exposure exposure() {
  return {sizeof(irred_sampled_photometry_exposure),
          IRRED_ABI_VERSION,
          0,
          0,
          1e10,
          .5,
          2,
          120};
}
irred_sampled_photometry_policy policy() {
  return {sizeof(irred_sampled_photometry_policy),
          IRRED_ABI_VERSION,
          7,
          0,
          64,
          1u << 20,
          64,
          64,
          64};
}
irred_sampled_photometry_batch
batch(std::span<const irred_sampled_photometry_curve> s,
      std::span<const irred_sampled_photometry_curve> t,
      std::span<const irred_sampled_photometry_exposure> e) {
  return {sizeof(irred_sampled_photometry_batch),
          IRRED_ABI_VERSION,
          s.data(),
          s.size(),
          s.size_bytes(),
          t.data(),
          t.size(),
          t.size_bytes(),
          e.data(),
          e.size(),
          e.size_bytes()};
}
irred_sampled_photometry_view view(irred_sampled_photometry_result *r) {
  irred_sampled_photometry_view v{};
  check(irred_sampled_photometry_result_view(r, &v) == IRRED_OK,
        "view transport");
  return v;
}
void reject(irred_sampled_photometry_batch b, irred_sampled_photometry_policy q,
            unsigned expected = IRRED_INVALID_INPUT) {
  auto *r = reinterpret_cast<irred_sampled_photometry_result *>(1);
  check(irred_sampled_photometry_evaluate(&b, &q, &r) == expected,
        "rejected transport");
  check(r == nullptr, "rejected null output");
}
void scalar(const irred_photometry_scalar &a, const p::Outcome &b) {
  check(a.availability == static_cast<unsigned>(b.availability),
        "exact availability");
  check(a.numerical_status == static_cast<unsigned>(b.numerical_status),
        "exact cause");
  if (b.value)
    check(std::bit_cast<std::uint64_t>(a.value) ==
              std::bit_cast<std::uint64_t>(*b.value),
          "exact native result bits");
  else
    check(a.value == 0, "absent wire payload");
}
} // namespace
int main() {
  try {
    std::array x{4e-7, 6e-7, 8e-7}, y{1e15, 2e15, 4e15}, u{3e-7, 1e-6, 1.5e-6},
        t{0., 1., 0.};
    std::array spectra{curve("spectrum-owned-long-identifier", x, y),
                       curve("second", x, y)};
    std::array passbands{curve("optical", u, t)};
    std::array rows{exposure(), exposure(), exposure(), exposure()};
    rows[1].redshift = -1;
    rows[2].spectrum_index = 99;
    rows[3].spectrum_index = 1;
    rows[3].collecting_area_square_metre = -0.;
    auto b = batch(spectra, passbands, rows);
    auto q = policy();
    irred_sampled_photometry_result *r = nullptr;
    check(irred_sampled_photometry_evaluate(&b, &q, &r) == IRRED_OK && r,
          "owned mixed batch");
    auto v = view(r);
    check(v.numerical_status == 0 && v.row_count == 4 &&
              v.spectrum_count == 2 && v.passband_count == 1,
          "row/pool counts");
    for (unsigned i : {0u, 1u, 3u}) {
      const auto &e = rows[i];
      auto nr = p::evaluate_sampled({{x, y},
                                     {u, t},
                                     e.luminosity_distance_metre,
                                     e.redshift,
                                     e.collecting_area_square_metre,
                                     e.observer_exposure_second},
                                    {7, 64, 64});
      check(v.rows[i].admission_status ==
                static_cast<unsigned>(nr.admission_status),
            "exact admission");
      scalar(v.rows[i].incident_band_flux, nr.flux_watt_per_square_metre);
      scalar(v.rows[i].collected_energy, nr.energy_joule);
      scalar(v.rows[i].expected_transmitted_photons, nr.expected_photons);
    }
    check(v.rows[2].admission_status ==
                  static_cast<unsigned>(n::Status::invalid_input) &&
              v.rows[2].collected_energy.availability == 2,
          "bad reference typed row");
    check(std::bit_cast<std::uint64_t>(
              v.rows[3].source.collecting_area_square_metre) ==
              std::bit_cast<std::uint64_t>(-0.),
          "source signed-zero bits");
    for (unsigned pool = 0; pool < 2; ++pool) {
      const auto *curves = pool ? v.passbands : v.spectra;
      const auto count = pool ? v.passband_count : v.spectrum_count;
      const auto *sources = pool ? passbands.data() : spectra.data();
      for (std::uint64_t i = 0; i < count; ++i) {
        check(curves[i].id.length == sources[i].id.length &&
                  !std::memcmp(curves[i].id.data, sources[i].id.data,
                               curves[i].id.length),
              "owned ID exact bytes/order");
        for (std::uint64_t j = 0; j < curves[i].length; ++j) {
          check(
              std::bit_cast<std::uint64_t>(curves[i].wavelength_metre[j]) ==
                  std::bit_cast<std::uint64_t>(sources[i].wavelength_metre[j]),
              "owned wavelength exact bits/order");
          check(std::bit_cast<std::uint64_t>(curves[i].values[j]) ==
                    std::bit_cast<std::uint64_t>(sources[i].values[j]),
                "owned values exact bits/order");
        }
      }
    }
    const auto old = x[0];
    x[0] = 1.;
    rows[0].redshift = 99;
    check(v.spectra[0].wavelength_metre[0] == old &&
              v.rows[0].source.redshift == .5,
          "immutable sources after mutation");
    check(v.spectra[0].wavelength_metre != x.data() &&
              v.spectra[0].values != y.data() &&
              v.spectra[0].id.data != spectra[0].id.data,
          "owned pool backing");
    check(v.model_id.length && v.method_id.length && v.constants_id.length &&
              v.arithmetic_id.length && v.propagation_id.length,
          "native metadata");
    irred_sampled_photometry_result_destroy(r);
    x[0] = old;
    rows[0] = exposure();
    q.requested_outputs = 2;
    check(irred_sampled_photometry_evaluate(&b, &q, &r) == IRRED_OK,
          "omitted groups");
    v = view(r);
    check(v.rows[0].incident_band_flux.availability == 0 &&
              v.rows[0].expected_transmitted_photons.availability == 0 &&
              v.rows[0].collected_energy.availability == 1,
          "only requested group");
    check(v.rows[2].incident_band_flux.availability == 0 &&
              v.rows[2].collected_energy.availability == 2,
          "invalid reference respects omission");
    irred_sampled_photometry_result_destroy(r);
    q = policy();
    q.maximum_samples = 5;
    check(irred_sampled_photometry_evaluate(&b, &q, &r) == IRRED_OK,
          "row work limit");
    v = view(r);
    check(v.numerical_status == 0 &&
              v.rows[0].admission_status ==
                  static_cast<unsigned>(n::Status::work_limit) &&
              v.rows[0].incident_band_flux.availability == 2,
          "row work limit normalized failed requested");
    irred_sampled_photometry_result_destroy(r);
    q = policy();
    auto original = rows[0];
    rows[0].luminosity_distance_metre = 1e308;
    rows[0].collecting_area_square_metre = 0;
    check(irred_sampled_photometry_evaluate(&b, &q, &r) == IRRED_OK,
          "partial group transport");
    v = view(r);
    check(v.rows[0].admission_status == 0 &&
              v.rows[0].incident_band_flux.availability == 2 &&
              v.rows[0].collected_energy.availability == 1 &&
              v.rows[0].expected_transmitted_photons.availability == 1,
          "underflow flux preserves zero collected groups");
    irred_sampled_photometry_result_destroy(r);
    rows[0] = original;
    q = policy();
    auto bad = b;
    bad.abi_version++;
    reject(bad, q, IRRED_ABI_MISMATCH);
    bad = b;
    bad.struct_size--;
    reject(bad, q);
    bad = b;
    bad.exposure_byte_length--;
    reject(bad, q);
    bad = b;
    bad.spectrum_count = UINT64_MAX;
    bad.spectra = reinterpret_cast<const irred_sampled_photometry_curve *>(1);
    reject(bad, q);
    auto wrong = q;
    wrong.requested_outputs = 8;
    reject(b, wrong);
    wrong = q;
    wrong.reserved = 1;
    reject(b, wrong);
    wrong = q;
    wrong.maximum_total_samples = 65537;
    reject(b, wrong);
    auto saved = spectra[0];
    spectra[0].length = 65537;
    spectra[0].wavelength_metre = reinterpret_cast<const double *>(1);
    reject(b, q);
    spectra[0] = saved;
    spectra[0].byte_length--;
    reject(b, q);
    spectra[0] = saved;
    spectra[0].abi_version++;
    reject(b, q, IRRED_ABI_MISMATCH);
    spectra[0] = saved;
    auto second = spectra[1];
    spectra[1].id = spectra[0].id;
    reject(b, q);
    spectra[1] = second;
    auto rb = rows[0];
    rows[0].struct_size--;
    reject(b, q);
    rows[0] = rb;
    for (unsigned kind = 0; kind < 3; ++kind) {
      q = policy();
      if (kind == 0)
        q.maximum_rows = 0;
      if (kind == 1)
        q.maximum_total_samples = 8;
      if (kind == 2)
        q.maximum_native_bytes = 0;
      check(irred_sampled_photometry_evaluate(&b, &q, &r) == IRRED_OK && r,
            "owned aggregate quota diagnostic");
      v = view(r);
      check(v.numerical_status ==
                    static_cast<unsigned>(n::Status::work_limit) &&
                !v.row_count && !v.spectrum_count && !v.passband_count,
            "empty quota diagnostic");
      irred_sampled_photometry_result_destroy(r);
    }
    std::uint64_t lo = 0, hi = 1u << 20;
    while (lo < hi) {
      q = policy();
      q.maximum_native_bytes = (lo + hi) / 2;
      check(irred_sampled_photometry_evaluate(&b, &q, &r) == IRRED_OK,
            "quota threshold probe");
      bool ok = view(r).numerical_status == 0;
      irred_sampled_photometry_result_destroy(r);
      if (ok)
        hi = q.maximum_native_bytes;
      else
        lo = q.maximum_native_bytes + 1;
    }
    q = policy();
    q.maximum_native_bytes = lo;
    check(irred_sampled_photometry_evaluate(&b, &q, &r) == IRRED_OK &&
              view(r).row_count == 4,
          "exact combined threshold");
    irred_sampled_photometry_result_destroy(r);
    q.maximum_native_bytes--;
    check(irred_sampled_photometry_evaluate(&b, &q, &r) == IRRED_OK &&
              !view(r).row_count,
          "one below threshold");
    irred_sampled_photometry_result_destroy(r);
    q = policy();
    const auto baseline = live;
    allocations = 0;
    fail_at = 1000000;
    auto transport = irred_sampled_photometry_evaluate(&b, &q, &r);
    const auto sites = allocations;
    fail_at = -1;
    check(transport == IRRED_OK && r, "allocation baseline");
    irred_sampled_photometry_result_destroy(r);
    check(live == baseline, "baseline recovered");
    for (long i = 0; i < sites; ++i) {
      allocations = 0;
      fail_at = i;
      transport = irred_sampled_photometry_evaluate(&b, &q, &r);
      fail_at = -1;
      check(transport == IRRED_ALLOCATION_FAILURE && !r,
            "allocation translated");
      check(live == baseline, "allocation cleanup");
    }
    allocations = 0;
    fail_at = sites - 1;
    nonallocation_exception = true;
    transport = irred_sampled_photometry_evaluate(&b, &q, &r);
    fail_at = -1;
    nonallocation_exception = false;
    check(transport == IRRED_EXCEPTION && !r && live == baseline,
          "non-allocation exception containment/cleanup");
    b = batch({}, {}, {});
    check(irred_sampled_photometry_evaluate(&b, &q, &r) == IRRED_OK &&
              view(r).numerical_status == 0,
          "empty admitted batch");
    irred_sampled_photometry_result_destroy(r);
    check(irred_sampled_photometry_result_destroy(nullptr) == IRRED_OK,
          "destroy null");
    std::cout << "PASS " << checks
              << " sampled photometry ABI controls allocation_sites=" << sites
              << " combined_threshold=" << lo << '\n';
  } catch (const std::exception &e) {
    fail_at = -1;
    std::cerr << "FAIL " << checks << ": " << e.what() << '\n';
    return 1;
  }
}
