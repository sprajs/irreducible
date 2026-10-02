#include "irred/recombination_drag.hpp"
#include <array>
#include <cstdlib>
#include <iostream>
#include <new>
// Original allocation-failure instrumentation; no solver/equation reference.
namespace {
bool armed = false;
std::size_t calls = 0, fail_on = 0;
long live = 0;
} // namespace
#if defined(__GNUC__) || defined(__clang__)
#define NOINLINE __attribute__((noinline))
#else
#define NOINLINE
#endif
NOINLINE void *operator new(std::size_t n) {
  if (armed && ++calls == fail_on)
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
int main() {
  using namespace irred::cosmology;
  using S = irred::numerics::Status;
  PureHydrogenRequest source{{67.4, .02237, .12, 2.7255, 1.7e-5, {}}};
  calls = 0;
  armed = true;
  auto original = prepare_pure_hydrogen_history(source);
  armed = false;
  const auto preparation_allocations = calls;
  if (original.status() != S::ok || preparation_allocations == 0)
    return 1;
  std::array<double, 1> z{1000};
  for (std::size_t point = 1; point <= preparation_allocations; ++point) {
    long before = live;
    {
      calls = 0;
      fail_on = point;
      armed = true;
      auto rejected = prepare_pure_hydrogen_history(source);
      armed = false;
      if (rejected.status() != S::work_limit ||
          !rejected.evaluate(z, 3).rows.empty())
        return 2;
    }
    if (live != before)
      return 3;
  }
  long before = live;
  calls = 0;
  fail_on = 1;
  armed = true;
  auto query = original.evaluate(z, 3);
  armed = false;
  if (query.status != S::work_limit || !query.rows.empty() || live != before)
    return 4;
  bool threw = false;
  calls = 0;
  fail_on = 1;
  armed = true;
  try {
    auto copy = original;
    (void)copy;
  } catch (const std::bad_alloc &) {
    threw = true;
  }
  armed = false;
  if (!threw || live != before)
    return 5;
  auto short_source = source;
  short_source.initial_redshift = 1590;
  PureHydrogenPolicy p;
  p.base_intervals = 512;
  auto destination = prepare_pure_hydrogen_history(short_source, p);
  if (destination.status() != S::ok)
    return 6;
  before = live;
  threw = false;
  calls = 0;
  fail_on = 1;
  armed = true;
  try {
    destination = original;
  } catch (const std::bad_alloc &) {
    threw = true;
  }
  armed = false;
  if (!threw || live != before ||
      destination.source()->initial_redshift != 1590) {
    std::cerr << "copy assignment exception changed retained source\n";
    return 7;
  }
  auto unchanged = original.evaluate(z, 3);
  if (!unchanged.rows[0].electron_fraction.value)
    return 8;
  std::cout << "PASS " << preparation_allocations
            << " preparation allocation failures, query, copy and assignment "
               "controls\n";
}
