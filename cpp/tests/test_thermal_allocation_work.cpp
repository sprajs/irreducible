// Linker-only observation of the real linked momentum integrator and allocator.
// The forwarding callback changes no arithmetic and is not a physics reference.
#include "irred/recombination_drag.hpp"
#include "../src/thermal_cc_sample.hpp"
#include <array>
#include <atomic>
#include <cstdlib>
#include <iostream>
#include <new>
#include <stdexcept>
namespace {
std::atomic<bool> armed = false;
std::atomic<std::size_t> allocations = 0, fail_on = 0, callbacks = 0,
                         integration_depth = 0, cc_callbacks = 0, adaptive_callbacks = 0, nested_allocations = 0,
                         post_callback_allocations = 0;
std::atomic<long> live = 0;
unsigned checks = 0;
void require(bool valid, const char *why) {
  ++checks;
  if (!valid)
    throw std::runtime_error(why);
}
void start(std::size_t failure = 0) {
  allocations = 0;
  callbacks = 0;
  cc_callbacks = 0; adaptive_callbacks = 0;
  nested_allocations = 0;
  post_callback_allocations = 0;
  fail_on = failure;
  armed = true;
}
struct Forward {
  irred::numerics::Integrand function;
  const void *context;
};
double observed(double x, const void *context) {
  if (armed) {
    ++callbacks; ++adaptive_callbacks;
  }
  const auto &f = *static_cast<const Forward *>(context);
  return f.function(x, f.context);
}
struct IntegrationScope {
  IntegrationScope() { ++integration_depth; }
  ~IntegrationScope() { --integration_depth; }
};
} // namespace
#if defined(__GNUC__) || defined(__clang__)
#define NOINLINE __attribute__((noinline))
#else
#define NOINLINE
#endif
NOINLINE void *operator new(std::size_t n) {
  if (armed) {
    if (integration_depth)
      ++nested_allocations;
    if (callbacks)
      ++post_callback_allocations;
    if (++allocations == fail_on)
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
// GNU/Clang Unix --wrap routes the unresolved production symbol to this
// observer. The original implementation and all returned diagnostics are used.
extern "C" irred::numerics::ScalarResult
real_integrate(irred::numerics::Integrand, const void *, double, double,
               irred::numerics::IntegrationPolicy)
    asm("__real__ZN5irred8numerics9integrateEPFddPKvES2_ddNS0_17IntegrationPolicyE");
extern "C" irred::numerics::ScalarResult
wrapped_integrate(irred::numerics::Integrand, const void *, double, double,
                  irred::numerics::IntegrationPolicy)
    asm("__wrap__ZN5irred8numerics9integrateEPFddPKvES2_ddNS0_17IntegrationPolicyE");
extern "C" irred::numerics::ScalarResult
wrapped_integrate(irred::numerics::Integrand f, const void *c, double a,
                  double b, irred::numerics::IntegrationPolicy p) {
  const Forward forward{f, c};
  const IntegrationScope scope;
  return real_integrate(observed, &forward, a, b, p);
}
extern "C" irred::cosmology::detail::ThermalCCSample
real_cc_sample(double,long double,long double,bool) noexcept
 asm("__real__ZN5irred9cosmology6detail17thermal_cc_sampleEdeeb");
extern "C" irred::cosmology::detail::ThermalCCSample
wrapped_cc_sample(double,long double,long double,bool) noexcept
 asm("__wrap__ZN5irred9cosmology6detail17thermal_cc_sampleEdeeb");
extern "C" irred::cosmology::detail::ThermalCCSample
wrapped_cc_sample(double q,long double y,long double scale,bool pressure) noexcept {
 if(armed) { callbacks+=pressure?2:1;cc_callbacks+=pressure?2:1; }
 const IntegrationScope scope;
 return real_cc_sample(q,y,scale,pressure);
}
namespace {
using namespace irred::cosmology;
using S = irred::numerics::Status;
void refinement_and_fallback() {
 ThermalPolicy p;p.momentum_method=ThermalMomentumMethod::nested_clenshaw_curtis;
 start();const auto refined=evaluate_thermal_moments(.001,p);armed=false;
 require(refined.status==S::ok&&refined.callbacks==callbacks&&cc_callbacks>1820&&adaptive_callbacks==0,"actual refined-parent node work retained");
 p.maximum_depth=0;start();const auto fallback=evaluate_thermal_moments(.1,p);armed=false;
 require(fallback.status!=S::ok&&fallback.callbacks==callbacks&&cc_callbacks>0&&adaptive_callbacks>0&&!fallback.rho_moment&&!fallback.pressure_moment,"discarded CC parent plus real failed adaptive fallback charged");
 const auto cap=fallback.callbacks-1;p.maximum_callbacks_per_evaluation=cap;
 start();const auto failed=evaluate_thermal_moments(.1,p);armed=false;
 require(failed.status==S::work_limit&&failed.callbacks==callbacks&&callbacks<=cap&&!failed.rho_moment&&!failed.pressure_moment,"failed fallback never resets remaining cap");
 std::cout<<"PASS refined/fallback independently forwarded counters\n";
}
void exercise(unsigned positive, ThermalMomentumMethod method) {
  ThermalPolicy selected;selected.momentum_method=method;
  constexpr double temperature = 1.95;
  constexpr double temperature_ev =
      static_cast<double>(static_cast<long double>(temperature) *
                          (1.380649e-23L / 1.602176634e-19L));
  constexpr double masses[]{.06, .001, .3};
  ThermalFlatModel source{67.4, 5.45e-5, 0, .049, .264, {}};
  for (unsigned i = 0; i < positive; ++i)
    source.species.push_back({masses[i], temperature_ev, 2});
  source.species.push_back({0, temperature_ev, 2});
  start();
  auto original = prepare_thermal_background(source, selected);
  armed = false;
  const auto provider_allocations = allocations.load();
  const auto preparation_work = original.preparation_callbacks();
  require(original.status() == S::ok && preparation_work == callbacks &&
              preparation_work > 0 && provider_allocations > 0,
          "actual successful provider normalization work");
  require(nested_allocations == 0 && post_callback_allocations == 0,
          "all provider storage acquired before momentum callbacks");
  const auto expected_species = original.omega_species_today();
  const auto expected_lambda = original.omega_lambda();
  for (std::size_t point = 1; point <= provider_allocations; ++point) {
    const auto before = live.load();
    bool threw = false;
    start(point);
    try {
      auto failed = prepare_thermal_background(source, selected);
      armed = false;
      require(failed.status() != S::ok, "faulted source acquisition refused");
    } catch (const std::bad_alloc &) {
      armed = false;
      threw = true;
    }
    require(threw && callbacks == 0 && live == before,
            "provider allocation failure precedes every FD callback");
  }
  for (std::size_t cap : {std::size_t(9), preparation_work - 1}) {
    ThermalPolicy p=selected;
    p.maximum_total_callbacks = cap;
    start();
    auto failed = prepare_thermal_background(source, p);
    armed = false;
    require(failed.status() == S::work_limit &&
                failed.preparation_callbacks() == callbacks && callbacks <= cap,
            "failed FD normalization retains actual callbacks within cap");
    require(failed.source().species.empty() &&
                failed.source().h0_km_s_mpc == 0 && nested_allocations == 0,
            "failed provider retains default source semantics");
  }
  auto repeated = prepare_thermal_background(source, selected);
  require(repeated.status() == S::ok &&
              repeated.preparation_callbacks() == preparation_work &&
              repeated.omega_species_today() == expected_species &&
              repeated.omega_lambda() == expected_lambda &&
              source.species.size() == positive + 1,
          "borrowed provider source remains reusable without changed scalars");
  PureHydrogenRequest hydrogen{{67.4, .02237, .12, 2.7255, 0, {}}, 1600, 300,
                               HydrogenTemperatureModel::evolved_compton_adiabatic};
  for (unsigned i = 0; i < positive; ++i)
    hydrogen.model.species.push_back({masses[i], temperature, 2});
  hydrogen.model.species.push_back({0, temperature, 2});
  PureHydrogenPolicy p;
  p.thermal=selected;
  p.base_intervals = 64;
  p.maximum_total_work = preparation_work + 20000;
  start();
  auto bounded = prepare_pure_hydrogen_history(hydrogen, p);
  armed = false;
  const auto hydrogen_allocations = allocations.load();
  require(bounded.status() == S::work_limit && hydrogen_allocations >= 4 &&
              bounded.work().momentum_callbacks == callbacks &&
              bounded.work().total() <= p.maximum_total_work,
          "small cap reaches outer-owned storage then bounds real FD work");
  require(nested_allocations == 0, "nested adaptive integrator allocates no storage");
  std::size_t faults_after_work = 0;
  for (std::size_t point = 1; point <= hydrogen_allocations; ++point) {
    const auto before = live.load();
    start(point);
    {
      auto failed = prepare_pure_hydrogen_history(hydrogen, p);
      armed = false;
      if (callbacks)
        ++faults_after_work;
      require(failed.status() == S::work_limit &&
                  failed.work().momentum_callbacks == callbacks &&
                  failed.work().total() <= p.maximum_total_work &&
                  failed.evaluate(std::array{1000.}, 127).rows.empty(),
              "outer allocation failure cannot erase completed FD work");
    }
    require(live == before, "failed outer acquisition leaks no owned storage");
  }
  require(faults_after_work > 0 && hydrogen.model.species.size() == positive + 1 &&
              hydrogen.model.species.front().mass_ev == .06,
          "faulted post-normalization outer acquisition preserves borrowed source");
  require(PureHydrogenPolicy{}.maximum_total_work == 4000000,
          "combined default cap unchanged");
  std::cout << "PASS positive_species=" << positive
            << " method=" << thermal_momentum_method_id(method)
            << " provider_allocations=" << provider_allocations
            << " normalization_callbacks=" << preparation_work
            << " outer_allocations=" << hydrogen_allocations
            << " post_work_faults=" << faults_after_work << '\n';
}
} // namespace
int main() {
  try {
    refinement_and_fallback();
    for (unsigned positive = 1; positive <= 3; ++positive)
      for(auto method:{ThermalMomentumMethod::direct_adaptive,ThermalMomentumMethod::nested_clenshaw_curtis})
        exercise(positive,method);
    std::cout << "PASS " << checks << " thermal allocation/work controls\n";
  } catch (const std::exception &e) {
    armed = false;
    std::cerr << "FAIL " << checks << ": " << e.what() << '\n';
    return 1;
  }
}
