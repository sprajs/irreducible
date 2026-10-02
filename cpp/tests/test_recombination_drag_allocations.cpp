#include "irred/recombination_drag.hpp"
#include <array>
#include <atomic>
#include <cstdlib>
#include <iostream>
#include <new>
// Original allocation-failure instrumentation; no solver/equation reference.
namespace {
// Observable instrumentation must survive optimization of allocator calls.
std::atomic<bool> armed = false;
std::atomic<std::size_t> calls = 0, fail_on = 0;
std::atomic<long> live = 0;
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
int exercise(irred::cosmology::HydrogenTemperatureModel temperature_model) {
  using namespace irred::cosmology;
  using S = irred::numerics::Status;
  PureHydrogenRequest source{{67.4, .02237, .12, 2.7255, 1.7e-5, {}}};
  source.temperature_model = temperature_model;
  source.model.species.push_back({0, 1.95, 2});
  calls = 0; fail_on = 0; armed = true;
  auto original = prepare_pure_hydrogen_history(source);
  armed = false;
  const auto preparation_allocations = calls.load();
  if (original.status() != S::ok || preparation_allocations == 0) return 1;
  std::array<double, 1> z{1000};
  for (std::size_t point = 1; point <= preparation_allocations; ++point) {
    const long before = live;
    {
      calls = 0; fail_on = point; armed = true;
      auto rejected = prepare_pure_hydrogen_history(source);
      armed = false;
      if (rejected.status() != S::work_limit ||
          !rejected.evaluate(z, 3).rows.empty()) return 2;
    }
    if (live != before) return 3;
  }
  const long before_query = live;
  calls = 0; fail_on = 1; armed = true;
  auto query = original.evaluate(z, 127);
  armed = false;
  if (query.status != S::work_limit || !query.rows.empty() || live != before_query) return 4;
  // Observe the copied owner through the external evaluator. An unused copy
  // can be optimized away, invalidating an allocation-site measurement.
  std::size_t copy_allocations = 0;
  {
    calls = 0; fail_on = 0; armed = true;
    auto copy = original;
    armed = false;
    copy_allocations = calls;
    const auto observed = copy.evaluate(z, 127);
    if (observed.rows.size()!=1 || !observed.rows[0].electron_fraction.value) return 5;
  }
  if (copy_allocations == 0) return 5;
  for (std::size_t point = 1; point <= copy_allocations; ++point) {
    const long before = live;
    bool threw = false;
    calls = 0; fail_on = point; armed = true;
    try {
      auto copy = original;
      armed = false;
      const auto observed = copy.evaluate(z, 127);
      if (observed.rows.size()!=1 || !observed.rows[0].electron_fraction.value) return 5;
    } catch (const std::bad_alloc &) { threw = true; }
    armed = false;
    if (!threw || live != before) return 5;
  }
  auto short_source = source;
  short_source.initial_redshift = 1590;
  PureHydrogenPolicy p;
  // Use a resolved destination: the coupled512 mesh truthfully refuses.
  p.base_intervals = 8192;
  auto destination = prepare_pure_hydrogen_history(short_source, p);
  if (destination.status() != S::ok) return 6;
  const auto old_query = destination.evaluate(z, 127);
  const auto old_work = destination.work();
  const auto old_model = destination.model_identity(), old_method = destination.method_identity();
  const auto old_root = destination.conditional_unit_depth_redshift();
  auto same_value=[](const HydrogenHistoryValue& a,const HydrogenHistoryValue& b) {
    return a.status==b.status && a.value==b.value &&
           a.absolute_error_estimate==b.absolute_error_estimate;
  };
  for (std::size_t point = 1; point <= copy_allocations; ++point) {
    const long before = live;
    bool threw = false;
    calls = 0; fail_on = point; armed = true;
    try { destination = original; }
    catch (const std::bad_alloc &) { threw = true; }
    armed = false;
    if (!threw || live != before || !destination.source() ||
        destination.source()->initial_redshift != 1590 ||
        destination.source()->temperature_model != temperature_model ||
        destination.source()->model.species.size()!=1 ||
        destination.source()->model.species[0].temperature_today_kelvin!=1.95 ||
        destination.model_identity()!=old_model || destination.method_identity()!=old_method ||
        !same_value(destination.conditional_unit_depth_redshift(),old_root)) return 7;
    const auto current_work=destination.work();
    if (current_work.background_evaluations!=old_work.background_evaluations ||
        current_work.momentum_callbacks!=old_work.momentum_callbacks ||
        current_work.rhs_evaluations!=old_work.rhs_evaluations ||
        current_work.equilibrium_solves!=old_work.equilibrium_solves ||
        current_work.temperature_rhs_evaluations!=old_work.temperature_rhs_evaluations) return 7;
    const auto current = destination.evaluate(z, 127);
    if (current.rows.size() != old_query.rows.size() || current.rows.empty() ||
        current.status != old_query.status || current.requested != old_query.requested) return 7;
    const auto& a=current.rows[0]; const auto& b=old_query.rows[0];
    if (a.redshift!=b.redshift || !same_value(a.electron_fraction,b.electron_fraction) ||
        !same_value(a.drag_depth,b.drag_depth) ||
        !same_value(a.matter_temperature_kelvin,b.matter_temperature_kelvin) ||
        !same_value(a.thomson_depth,b.thomson_depth) ||
        !same_value(a.thomson_opacity_per_redshift,b.thomson_opacity_per_redshift) ||
        !same_value(a.visibility_per_redshift,b.visibility_per_redshift) ||
        !same_value(a.finite_endpoint_survival,b.finite_endpoint_survival)) return 7;
  }
  const auto unchanged = original.evaluate(z, 3);
  if (!unchanged.rows[0].electron_fraction.value) return 8;
  std::cout << "PASS " << original.model_identity() << ": " << preparation_allocations
            << " preparation and " << copy_allocations
            << " copy/assignment allocation sites, query and retained-state controls\n";
  return 0;
}
int main() {
  using irred::cosmology::HydrogenTemperatureModel;
  if (const int result=exercise(HydrogenTemperatureModel::prescribed_radiation)) return result;
  return exercise(HydrogenTemperatureModel::evolved_compton_adiabatic);
}
