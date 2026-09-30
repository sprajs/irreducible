// Differential layout regression driver: compile identical source against
// frozen baseline/candidate kernels, compare stdout byte-for-byte. Every public
// binary64 field is encoded by bits; padding is never compared. Test author
// independent of numerical implementation. Dense integer Gram inputs use no
// external assets.
#include "irred/numerics.hpp"
#include <array>
#include <bit>
#include <cfenv>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <new>
#include <stdexcept>
#include <utility>
#include <vector>
using namespace irred::numerics;
namespace {
std::size_t allocated = 0, allocations = 0;
}
void *operator new(std::size_t n) {
  auto p = std::malloc(n ? n : 1);
  if (!p)
    throw std::bad_alloc();
  allocated += n;
  ++allocations;
  return p;
}
void operator delete(void *p) noexcept { std::free(p); }
void operator delete(void *p, std::size_t) noexcept { ::operator delete(p); }
void *operator new[](std::size_t n) { return ::operator new(n); }
void operator delete[](void *p) noexcept { ::operator delete(p); }
void operator delete[](void *p, std::size_t) noexcept { ::operator delete(p); }
namespace {
int checks = 0;
void check(bool b, const char *s) {
  ++checks;
  if (!b)
    throw std::runtime_error(s);
}
void bits(double x) {
  std::printf(" %016llx",
              static_cast<unsigned long long>(std::bit_cast<std::uint64_t>(x)));
}
void emit(const char *label, const Factorization &f, const SolveResult &s,
          std::size_t bytes, std::size_t calls) {
  std::printf("%s %u %zu %u %u %zu %zu", label,
              static_cast<unsigned>(f.status()), f.size(),
              static_cast<unsigned>(f.arithmetic()),
              static_cast<unsigned>(s.status), bytes, calls);
  bits(f.log_determinant());
  bits(f.condition_estimate_inf());
  for (auto x : {s.backward_residual, s.estimated_forward_sensitivity,
                 s.pre_cast_backward_residual, s.pre_cast_sensitivity_estimate,
                 s.output_rounding_error_inf, s.output_rounding_error_relative,
                 s.post_cast_backward_residual})
    bits(x);
  for (auto x : s.value)
    bits(x);
  std::printf("\n");
}
} // namespace
int main() {
  try {
    for (auto arithmetic :
         {Arithmetic::binary64_legacy_v1, Arithmetic::longdouble_cpu_v1})
      for (std::size_t n : {1u, 2u, 5u, 17u}) {
        std::vector<double> A(n * n), C(n * n), r(n), x(n);
        for (std::size_t i = 0; i < n; ++i) {
          A[i * n + i] = 4 + i;
          x[i] = (i % 2 ? -1. : 1.) / (i + 1);
          for (std::size_t j = 0; j < i; ++j)
            A[i * n + j] = static_cast<int>((i + 3 * j) % 5) - 2;
        }
        for (std::size_t i = 0; i < n; ++i)
          for (std::size_t j = 0; j < n; ++j)
            for (std::size_t k = 0; k < n; ++k)
              C[i * n + j] += A[i * n + k] * A[j * n + k];
        for (std::size_t i = 0; i < n; ++i)
          for (std::size_t j = 0; j < n; ++j)
            r[i] += C[i * n + j] * x[j];
        auto initial_bytes = allocated, initial_calls = allocations;
        auto f = cholesky(C, n, n * n, arithmetic);
        auto solved = solve(f, r, 1e-8);
        auto bytes = allocated - initial_bytes,
             calls = allocations - initial_calls;
        check(f.status() == Status::ok && solved.status == Status::ok,
              "dense Gram solves");
        for (std::size_t i = 0; i < n; ++i)
          check(std::abs(solved.value[i] - x[i]) < 1e-11, "known solution");
        emit("Gram", f, solved, bytes, calls);
        auto copy = f;
        auto copied = solve(copy, r, 1e-8);
        check(copied.value == solved.value, "copied factor");
        emit("copy", copy, copied, 0, 0);
        auto moved = std::move(copy);
        auto moved_result = solve(moved, r, 1e-8);
        check(moved_result.value == solved.value, "moved factor");
        emit("moved", moved, moved_result, 0, 0);
        auto empty_move = solve(copy, r, 1e-8);
        check(empty_move.status == Status::invalid_input &&
                  empty_move.value.empty(),
              "moved-from guarded");
        emit("moved-from", copy, empty_move, 0, 0);
        auto tight = solve(f, r, 1e-30);
        check(tight.status == Status::conditioning_budget_exceeded,
              "tight budget");
        emit("tight", f, tight, 0, 0);
        r[0] = std::numeric_limits<double>::quiet_NaN();
        auto invalid = solve(f, r, 1e-8);
        check(invalid.status == Status::nonfinite_input &&
                  invalid.value.empty(),
              "nonfinite RHS");
        emit("nonfinite", f, invalid, 0, 0);
      }
    Factorization empty;
    std::array<double, 1> one{1};
    auto nofactor = solve(empty, one, 1e-8);
    check(nofactor.status == Status::invalid_input, "default guard");
    emit("default", empty, nofactor, 0, 0);
    std::array<double, 4> asymmetric{2, 1, 0, 2}, indefinite{1, 2, 2, 1};
    for (auto arithmetic :
         {Arithmetic::binary64_legacy_v1, Arithmetic::longdouble_cpu_v1}) {
      auto bad = cholesky(asymmetric, 2, 4, arithmetic);
      check(bad.status() == Status::invalid_input, "source symmetry");
      emit("asym", bad, {}, 0, 0);
      auto notspd = cholesky(indefinite, 2, 4, arithmetic);
      check(notspd.status() != Status::ok, "indefinite");
      emit("indefinite", notspd, {}, 0, 0);
      auto cap = cholesky(indefinite, 2, 3, arithmetic);
      check(cap.status() == Status::work_limit, "cap");
      emit("cap", cap, {}, 0, 0);
    }
    auto unknown = cholesky(one, 1, 1, static_cast<Arithmetic>(255));
    check(unknown.status() == Status::invalid_input, "enum guard");
    emit("unknown", unknown, {}, 0, 0);
    auto wide = cholesky(one, 1, 1, Arithmetic::longdouble_cpu_v1);
    auto original_round = std::fegetround();
    check(std::fesetround(FE_DOWNWARD) == 0, "rounding set");
    auto rounding = solve(wide, one, 1e-8);
    check(rounding.status == Status::outside_domain, "changed rounding");
    check(std::fesetround(original_round) == 0, "rounding restore");
    emit("rounding", wide, rounding, 0, 0);
    std::printf("mirror hostile checks %d PASS\n", checks);
  } catch (const std::exception &e) {
    std::fprintf(stderr, "%s\n", e.what());
    return 1;
  }
}
