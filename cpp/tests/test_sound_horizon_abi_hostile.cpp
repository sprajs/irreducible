// Transport equivalence and lifetime/resource controls, not independent
// physics.
#include "irred/abi.h"
#include "irred/sound_horizon.hpp"
#include <array>
#include <bit>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <new>
#include <stdexcept>
#include <string>
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
namespace c = irred::cosmology;
namespace n = irred::numerics;
irred_sound_horizon_policy policy() {
  return {sizeof(irred_sound_horizon_policy),
          IRRED_ABI_VERSION,
          1e-9,
          2e-11,
          100000,
          30,
          32,
          1000000,
          1 << 24};
}
irred_sound_horizon_input input(const std::string &s) {
  return {sizeof(irred_sound_horizon_input),
          IRRED_ABI_VERSION,
          70,
          .3,
          9e-5,
          .05,
          5e-5,
          1059,
          {reinterpret_cast<const uint8_t *>(s.data()), s.size()}};
}
irred_sound_horizon_batch batch(const irred_sound_horizon_input *p, size_t n) {
  return {sizeof(irred_sound_horizon_batch), IRRED_ABI_VERSION, p, n,
          n * sizeof(*p)};
}
irred_sound_horizon_view view(irred_sound_horizon_result *p) {
  irred_sound_horizon_view v{};
  check(irred_sound_horizon_result_view(p, &v) == IRRED_OK, "view transport");
  return v;
}
void bad(irred_sound_horizon_batch b, irred_sound_horizon_policy p,
         uint32_t cause = IRRED_INVALID_INPUT) {
  auto *out = reinterpret_cast<irred_sound_horizon_result *>(uintptr_t(1));
  check(irred_sound_horizon_evaluate(&b, &p, &out) == cause,
        "structural cause");
  check(!out, "structural owner reset");
}
void tests() {
  std::string origin(300, 'x');
  auto x = input(origin);
  auto y = x;
  y.omega_m = 1;
  y.omega_r = std::numeric_limits<double>::denorm_min();
  std::array xs{x, y};
  auto b = batch(xs.data(), xs.size());
  auto p = policy();
  c::SoundHorizonRequest r{
      {x.h0_km_s_mpc, x.omega_m, x.omega_r, x.omega_b, x.omega_gamma},
      x.z_drag,
      origin};
  auto direct = c::evaluate_sound_horizon(
      std::span(&r, 1),
      {p.absolute_tolerance_mpc, p.relative_tolerance,
       p.maximum_callbacks_per_point, unsigned(p.maximum_depth),
       p.maximum_points, p.maximum_total_callbacks, p.maximum_native_bytes});
  check(direct.rows.size() == 1 && direct.rows[0].sound_horizon_mpc.has_value(),
        "native finite");
  irred_sound_horizon_result *out = nullptr;
  check(irred_sound_horizon_evaluate(&b, &p, &out) == IRRED_OK && out,
        "mixed transport");
  // A misaligned offset of a real owner must reject before any dereference.
  // This does not promise validation of arbitrary aligned forged pointers.
  auto *misaligned_owner = reinterpret_cast<irred_sound_horizon_result *>(
      reinterpret_cast<uintptr_t>(out) + 1);
  irred_sound_horizon_view rejected_view{};
  check(irred_sound_horizon_result_view(misaligned_owner, &rejected_view) ==
            IRRED_INVALID_INPUT,
        "misaligned result view rejected");
  check(irred_sound_horizon_result_destroy(misaligned_owner) ==
            IRRED_INVALID_INPUT,
        "misaligned destroy rejected");
  alignas(irred_sound_horizon_view)
      std::array<unsigned char, sizeof(irred_sound_horizon_view) + 1>
          view_storage{};
  check(irred_sound_horizon_result_view(
            out, reinterpret_cast<irred_sound_horizon_view *>(
                     view_storage.data() + 1)) == IRRED_INVALID_INPUT,
        "misaligned view destination rejected");
  check(irred_sound_horizon_result_view(nullptr, &rejected_view) ==
            IRRED_INVALID_INPUT,
        "null owner view rejected");
  check(irred_sound_horizon_result_destroy(nullptr) == IRRED_OK,
        "null destroy remains harmless");
  auto v = view(out);
  check(v.count == 2 && v.numerical_status == 0, "mixed container");
  check(v.rows[0].has_value && v.rows[0].numerical_status == 0, "finite row");
  check(std::bit_cast<uint64_t>(v.rows[0].sound_horizon_mpc) ==
            std::bit_cast<uint64_t>(*direct.rows[0].sound_horizon_mpc),
        "native value bits");
  check(v.rows[0].callbacks == direct.rows[0].callbacks, "native work");
  check(v.rows[1].numerical_status != 0 && !v.rows[1].has_value,
        "failed source no value");
  origin.assign(300, 'y');
  xs[0].h0_km_s_mpc = 99;
  check(v.rows[0].source.h0_km_s_mpc == 70, "owned scalars");
  check(v.rows[0].source.drag_origin.length == 300 &&
            v.rows[0].source.drag_origin.data[0] == 'x',
        "owned long string");
  irred_sound_horizon_result_destroy(out);
  b = batch(&x, 1);
  x = input(origin);
  auto q = p;
  q.maximum_native_bytes = 0;
  check(irred_sound_horizon_evaluate(&b, &q, &out) == IRRED_OK && out,
        "cap diagnostic owner");
  v = view(out);
  check(v.numerical_status == uint32_t(n::Status::work_limit) && !v.count &&
            !v.callbacks,
        "cap empty no work");
  irred_sound_horizon_result_destroy(out);
  uint64_t low = 0, high = 1 << 24;
  while (low + 1 < high) {
    uint64_t mid = (low + high) / 2;
    q = p;
    q.maximum_native_bytes = mid;
    check(irred_sound_horizon_evaluate(&b, &q, &out) == IRRED_OK && out,
          "quota probe");
    v = view(out);
    if (v.count)
      high = mid;
    else
      low = mid;
    irred_sound_horizon_result_destroy(out);
  }
  for (auto cap : {high, high - 1}) {
    q = p;
    q.maximum_native_bytes = cap;
    check(irred_sound_horizon_evaluate(&b, &q, &out) == IRRED_OK,
          "quota boundary");
    v = view(out);
    check((v.count == 1) == (cap == high), "exact/one under");
    irred_sound_horizon_result_destroy(out);
  }
  // Constant compact integrand isolates the arithmetic floor from refinement
  // limits.
  x.omega_m = 0;
  x.omega_r = 1;
  x.omega_b = 0;
  x.omega_gamma = 1;
  x.z_drag = 1000;
  q = p;
  q.absolute_tolerance_mpc = 1e-30;
  q.relative_tolerance = 0;
  check(irred_sound_horizon_evaluate(&b, &q, &out) == IRRED_OK,
        "tight transport");
  v = view(out);
  check(v.count == 1 && !v.rows[0].has_value &&
            v.rows[0].numerical_status ==
                uint32_t(n::Status::conditioning_budget_exceeded),
        "tight precise scientific cause");
  irred_sound_horizon_result_destroy(out);
  auto malformed = b;
  malformed.byte_length++;
  bad(malformed, p);
  malformed = b;
  malformed.abi_version++;
  bad(malformed, p, IRRED_ABI_MISMATCH);
  malformed = b;
  malformed.points = nullptr;
  bad(malformed, p);
  malformed = b;
  malformed.count = UINT64_MAX;
  malformed.points =
      reinterpret_cast<const irred_sound_horizon_input *>(uintptr_t(8));
  bad(malformed, p);
  q = p;
  q.maximum_depth = 61;
  bad(b, q);
  q = p;
  q.relative_tolerance = -1;
  bad(b, q);
  q = p;
  q.maximum_native_bytes = 1073741825ULL;
  bad(b, q);
  x.drag_origin = {nullptr, 1};
  bad(b, p);
  x = input(origin);
  x.abi_version++;
  bad(b, p, IRRED_ABI_MISMATCH);
  x = input(origin);
  const long baseline = live;
  fail_at = 100000;
  allocations = 0;
  check(irred_sound_horizon_evaluate(&b, &p, &out) == IRRED_OK,
        "allocation baseline");
  long sites = allocations;
  fail_at = -1;
  irred_sound_horizon_result_destroy(out);
  check(live == baseline, "baseline cleanup");
  for (long i = 0; i < sites; ++i) {
    fail_at = i;
    allocations = 0;
    auto status = irred_sound_horizon_evaluate(&b, &p, &out);
    fail_at = -1;
    check(status == IRRED_ALLOCATION_FAILURE && !out,
          "allocation precise null");
    check(live == baseline, "allocation cleanup");
  }
  b = batch(nullptr, 0);
  q = p;
  q.maximum_points = 0;
  check(irred_sound_horizon_evaluate(&b, &q, &out) == IRRED_OK,
        "empty transport");
  v = view(out);
  check(!v.count && v.numerical_status == 0, "admitted empty OK");
  irred_sound_horizon_result_destroy(out);
  std::cout << "PASS" << checks
            << " sound ABI controls allocation_sites=" << sites
            << " combined_threshold=" << high << '\n';
}
} // namespace
int main() {
  try {
    tests();
  } catch (const std::exception &e) {
    std::cerr << "FAIL check" << checks << " " << e.what() << '\n';
    return 1;
  }
}
