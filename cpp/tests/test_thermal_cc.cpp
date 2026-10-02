// Opt-in quadrature owner controls. Direct comparisons share physical equations;
// independent numerical references belong to the separately authored peer suite.
#include "irred/thermal_neutrino.hpp"
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <utility>
namespace {
using namespace irred::cosmology;
using S=irred::numerics::Status;
unsigned checks=0;
void require(bool good,const char *why) { ++checks; if(!good) throw std::runtime_error(why); }
void near(double a,double b,long double scale=1) {
  require(std::abs((static_cast<long double>(a)-b)/scale)<=1e-10L+2e-10L*std::abs(b/scale),"matched direct comparison");
}
ThermalPolicy cc() { ThermalPolicy p; p.momentum_method=ThermalMomentumMethod::nested_clenshaw_curtis; return p; }
void moments() {
  auto p=cc();
  for(double y:{0.,1e-12,1e-8,.001,.01,.1,1.,10.,100.,1e4,1e6,1e100}) {
    auto a=evaluate_thermal_moments(y,p),b=evaluate_thermal_moments(y);
    require(a.status==S::ok && b.status==S::ok && a.rho_moment && a.pressure_moment,"both methods admitted");
    const long double scale=std::hypot(1.L,static_cast<long double>(y));
    near(*a.rho_moment,*b.rho_moment,scale);
    near(*a.pressure_moment,*b.pressure_moment,1/scale);
    require(a.rho_error_estimate/scale<=p.absolute_tolerance+p.relative_tolerance*(*a.rho_moment/scale),"rho diagnostic admitted");
    require(a.pressure_error_estimate*scale<=p.absolute_tolerance+p.relative_tolerance*(*a.pressure_moment*scale),"pressure diagnostic admitted");
    if(y==0) {
      require(a.callbacks==0 && b.callbacks==0,"massless has no work");
      require(std::bit_cast<std::uint64_t>(*a.rho_moment)==std::bit_cast<std::uint64_t>(*b.rho_moment),"massless rho bits");
      require(std::bit_cast<std::uint64_t>(*a.pressure_moment)==std::bit_cast<std::uint64_t>(*b.pressure_moment),"massless pressure bits");
    }
  }
  for(std::size_t cap:{0u,1u,2u,3u,7u,31u,129u}) {
    p=cc();p.maximum_callbacks_per_evaluation=cap;
    const auto a=evaluate_thermal_moments(1,p);
    require(a.status==S::work_limit && a.callbacks<=cap && !a.rho_moment && !a.pressure_moment,"failed nodes charged without partial values");
    require(a.callbacks==cap-cap%2,"paired scalar node charges");
  }
  p=cc();p.maximum_depth=0;
  auto fallback=evaluate_thermal_moments(.1,p);
  require(fallback.callbacks>130,"unresolved quadrature attempts retained");
  // A zero-depth policy need not admit the direct fallback either.
  p=cc();p.absolute_tolerance=p.relative_tolerance=1e-30;
  auto bad=evaluate_thermal_moments(1,p);
  require(bad.status==S::conditioning_budget_exceeded && bad.callbacks==0,"unattainable arithmetic allocation refused before nodes");
  p=cc();p.momentum_method=static_cast<ThermalMomentumMethod>(23);
  require(evaluate_thermal_moments(1,p).status==S::invalid_input,"invalid method refused");
  require(thermal_momentum_method_id(p.momentum_method).empty(),"invalid method has no identity");
  p=cc();
  require(evaluate_thermal_moments(-1,p).status==S::outside_domain,"negative y refused");
  require(evaluate_thermal_moments(std::numeric_limits<double>::infinity(),p).status==S::nonfinite_input,"nonfinite y refused");
}
void owner() {
  const ThermalFlatModel source{67.4,5.44e-5,0,.049,.264,{{.06,.000168,2}}};
  auto p=cc();auto a=prepare_thermal_background(source,p),direct=prepare_thermal_background(source);
  require(a.status()==S::ok && direct.status()==S::ok,"both retained owners admitted");
  require(a.momentum_method()==p.momentum_method,"prepared method retained");
  require(a.scaled_expansion(1).status==S::invalid_input,"method mismatch rejected at identity endpoint");
  const std::array<double,4> scales{1,.1,.001,1e-8};
  require(a.evaluate(scales,thermal_e).status==S::invalid_input,"batch mismatch rejected");
  auto x=a.evaluate(scales,thermal_e|thermal_h,p),y=direct.evaluate(scales,thermal_e|thermal_h);
  require(x.status==S::ok && y.status==S::ok,"matched background batches admitted");
  for(std::size_t i=0;i<scales.size();++i) { require(x.rows[i].e.value.has_value() && x.rows[i].h_km_s_mpc.value.has_value(),"requested background outputs"); near(*x.rows[i].e.value,*y.rows[i].e.value); }
  require(x.rows[0].callbacks==0 && *x.rows[0].e.value==1,"normalization identity no work");
  require(a.scaled_expansion(0,p).status==S::ok && a.scaled_expansion(0,p).callbacks==0,"radiation endpoint no nodes");
  auto copy=a;auto moved=std::move(a);
  require(copy.momentum_method()==p.momentum_method && moved.momentum_method()==p.momentum_method,"copy and move retain method");
  require(a.status()==S::invalid_input && a.source().species.empty(),"moved source invalidated");
  auto *self=&moved; moved=std::move(*self);require(moved.status()==S::ok && moved.momentum_method()==p.momentum_method,"self move preserves method");
  p.maximum_total_callbacks=7;
  auto fail=prepare_thermal_background(source,p);
  require(fail.status()==S::work_limit && fail.preparation_callbacks()==7,"density-only preparation exact scalar cap");
  p=cc();p.maximum_native_bytes=0;
  require(prepare_thermal_background(source,p).status()==S::work_limit,"payload admission precedes quadrature");
}
}
int main() { try { moments();owner();std::cout<<"PASS "<<checks<<" thermal CC owner controls\n"; } catch(const std::exception &e) { std::cerr<<"FAIL "<<e.what()<<'\n';return 1; } }
