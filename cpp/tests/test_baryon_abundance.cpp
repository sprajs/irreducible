#include "irred/baryon_abundance.hpp"
#include <array>
#include <bit>
#include <cfenv>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <limits>
#include <new>
#include <stdexcept>
static long fail_at=-1, allocations=0, live=0;
void *operator new(std::size_t n) {
  if (fail_at>=0 && allocations++==fail_at) throw std::bad_alloc();
  if (void *p=std::malloc(n?n:1)) { ++live; return p; }
  throw std::bad_alloc();
}
void *operator new[](std::size_t n) { return ::operator new(n); }
void operator delete(void *p) noexcept { if(p) { --live; std::free(p); } }
void operator delete[](void *p) noexcept { ::operator delete(p); }
void operator delete(void *p,std::size_t) noexcept { ::operator delete(p); }
void operator delete[](void *p,std::size_t) noexcept { ::operator delete(p); }
namespace {
namespace c=irred::cosmology; namespace a=irred::atomic;
using S=irred::numerics::Status; using W=long double;
unsigned checks=0;
void check(bool good,const char *m) { ++checks; if(!good) throw std::runtime_error(m); }
void close(W x,W y,W e,const char *m) { check(std::abs(x/y-1)<=e,m); }
c::BaryonAbundanceSource source() {
  return {.physical_baryon_density=.0224,.helium4_mass_fraction=.245,
    .hydrogen1_effective_mass_kg=1.67353284e-27,.helium4_effective_mass_kg=6.64647907e-27,
    .source_origin=std::string(90,'s'),.mass_origin=std::string(90,'m')};
}
void controls() {
  auto s=source(); auto owner=c::prepare_baryon_abundance(s);
  check(owner.status()==S::ok,"prepare"); s.source_origin="changed";
  check(owner.source()->source_origin.size()==90,"owned source");
  std::array<double,4> grid{1,.5,.001,1}; auto b=owner.evaluate(grid);
  check(b.status==S::ok && b.rows.size()==4,"ordered batch");
  for(auto &r:b.rows) {
    check(r.admission_status==S::ok && r.hydrogen_nuclei.status==S::ok &&
      r.helium_nuclei.status==S::ok,"map admission");
    check(r.hydrogen_nuclei.relative_arithmetic_estimate>0 &&
      r.hydrogen_nuclei.relative_arithmetic_estimate<=1e-15,"diagnostic");
    W nh=*r.hydrogen_nuclei.value, nhe=*r.helium_nuclei.value;
    W mh=owner.source()->hydrogen1_effective_mass_kg,mhe=owner.source()->helium4_effective_mass_kg;
    W y=mhe*nhe/(mh*nh+mhe*nhe);
    close(y,owner.source()->helium4_mass_fraction,5e-16L,"inverse Y");
    close(nhe/nh,W(.245)/(1-W(.245))*mh/mhe,5e-16L,"unequal mass number ratio");
  }
  close(*b.rows[1].hydrogen_nuclei.value,*b.rows[0].hydrogen_nuclei.value*8,2e-16L,"a cubed");
  check(std::bit_cast<uint64_t>(*b.rows[0].hydrogen_nuclei.value)==
    std::bit_cast<uint64_t>(*b.rows[3].hydrogen_nuclei.value),"repeat bits");
  auto copy=owner; auto moved=std::move(copy);
  check(copy.status()==S::invalid_input && !copy.source() && moved.status()==S::ok,"move invalidation");
  auto *self=&moved; moved=std::move(*self); check(moved.status()==S::ok,"self move");
  copy=owner; copy=copy; check(copy.status()==S::ok,"copy self");
  c::BaryonAbundance assigned; assigned=owner;
  check(assigned.source()->source_origin==owner.source()->source_origin,"copy assignment");
  for(double y:{0.,1.,std::nextafter(0.,1.),std::nextafter(1.,0.)}) {
    auto e=source(); e.helium4_mass_fraction=y;
    auto endpoint=c::prepare_baryon_abundance(e); auto r=endpoint.evaluate(std::array{1.}).rows[0];
    check(r.hydrogen_nuclei.status==S::ok,"endpoint H map");
    if(y==std::nextafter(0.,1.))
      check(r.helium_nuclei.status==S::conditioning_budget_exceeded && !r.helium_nuclei.value,"adjacent endpoint positive refusal");
    else check(r.helium_nuclei.status==S::ok,"endpoint He map");
    if(y==0) check(*r.helium_nuclei.value==0 && r.helium_nuclei.relative_arithmetic_estimate==0,"exact He zero");
    if(y==1) check(*r.hydrogen_nuclei.value==0,"exact H zero");
    if(y==0 || y==1) {
      auto l=endpoint.evaluate_equilibrium(std::array{c::BaryonAbundanceLteQuery{.001,8000}},"T supplied");
      check(l.solves==0 && l.rows[0].equilibrium.admission_status==S::outside_domain,"LTE absent refusal");
    }
  }
  auto zero=source(); zero.physical_baryon_density=0;
  auto zr=c::prepare_baryon_abundance(zero).evaluate(std::array{1.}).rows[0];
  check(*zr.hydrogen_nuclei.value==0 && *zr.helium_nuclei.value==0,"zero mass density");
  auto extreme=source(); extreme.physical_baryon_density=std::numeric_limits<double>::max();
  extreme.hydrogen1_effective_mass_kg=std::numeric_limits<double>::min();
  extreme.helium4_effective_mass_kg=std::numeric_limits<double>::max();
  auto er=c::prepare_baryon_abundance(extreme).evaluate(std::array{1.}).rows[0];
  check(er.hydrogen_nuclei.status==S::conditioning_budget_exceeded && !er.hydrogen_nuclei.value &&
    er.helium_nuclei.status==S::ok,"independent representation");
  auto tiny=source(); tiny.physical_baryon_density=std::numeric_limits<double>::denorm_min();
  tiny.hydrogen1_effective_mass_kg=tiny.helium4_effective_mass_kg=std::numeric_limits<double>::max();
  auto tr=c::prepare_baryon_abundance(tiny).evaluate(std::array{1.}).rows[0];
  check(tr.hydrogen_nuclei.status==S::conditioning_budget_exceeded && !tr.hydrogen_nuclei.value,"positive underflow refusal");
  auto invalid=owner.evaluate(std::array{0.,-1.,2.,std::numeric_limits<double>::quiet_NaN()});
  check(invalid.rows[0].admission_status==S::outside_domain &&
    invalid.rows[3].admission_status==S::nonfinite_input,"query failure");
  for(double bad:{-1.,std::numeric_limits<double>::infinity()}) {
    auto t=source();t.physical_baryon_density=bad;
    check(c::prepare_baryon_abundance(t).status()!=S::ok,"source invalid");
  }
  auto badmass=source(); badmass.helium4_effective_mass_kg=std::numeric_limits<double>::denorm_min();
  check(c::prepare_baryon_abundance(badmass).status()==S::outside_domain,"mass normal gate");
  c::BaryonAbundancePolicy p;p.maximum_rows=0;
  check(owner.evaluate(grid,p).status==S::work_limit,"row quota");
  p={};p.maximum_native_bytes=1;
  check(c::prepare_baryon_abundance(source(),p).status()==S::work_limit,"source bytes");
  check(owner.evaluate(grid,p).status==S::work_limit,"map bytes");
  check(owner.evaluate_equilibrium(std::array{c::BaryonAbundanceLteQuery{.001,8000}},"T",{},p).status==S::work_limit,"LTE peak bytes");
  p={};p.maximum_rows=65537;
  check(owner.evaluate(grid,p).status==S::invalid_input,"hard policy cap");
  auto empty=owner.evaluate({});check(empty.status==S::ok && empty.rows.empty(),"empty map");
  std::array<c::BaryonAbundanceLteQuery,3> q{{{.001,8000},{0,8000},{.001,8000}}};
  auto l=owner.evaluate_equilibrium(q,"independently supplied T");
  check(l.status==S::ok && l.rows.size()==3 && l.solves==2,"one batch valid work");
  check(l.rows[1].equilibrium.admission_status==S::outside_domain,"bad a LTE status");
  check(l.rows[0].equilibrium.electron_density.status==S::ok,"shared LTE electron");
  check(l.temperature_origin=="independently supplied T","T lineage");
  auto native=a::evaluate_hydrogen_helium_equilibrium(std::array{l.rows[0].equilibrium.source});
  auto &v=l.rows[0].equilibrium.electron_density;
  check(*v.value==*native.rows[0].electron_density.value &&
    v.relative_arithmetic_estimate>native.rows[0].electron_density.relative_arithmetic_estimate,"inherited diagnostic");
  a::HydrogenHeliumPolicy h;h.requested_outputs=a::electron_number_density; h.maximum_solves=1;
  auto limited=owner.evaluate_equilibrium(q,"T",h);
  check(limited.solves==1 && limited.rows[2].equilibrium.electron_density.status==S::work_limit,"global LTE work");
  check(!limited.rows[0].equilibrium.hydrogen_neutral.value,"mask omission");
  h={};h.maximum_native_bytes=1;
  check(owner.evaluate_equilibrium(q,"T",h).status==S::work_limit,"callee bytes");
  h={};h.requested_outputs=0;
  check(owner.evaluate_equilibrium(q,"T",h).status==S::invalid_input,"mask invalid");
  check(owner.evaluate_equilibrium(q,"").status==S::invalid_input,"missing T origin");
  std::fesetround(FE_DOWNWARD);
  auto rounding=owner.evaluate(std::array{1.}); std::fesetround(FE_TONEAREST);
  check(rounding.rows[0].admission_status==S::invalid_input,"rounding profile");
}
void allocation_controls() {
  auto s=source();auto owner=c::prepare_baryon_abundance(s);
  const std::array<c::BaryonAbundanceLteQuery,1> q{{{.001,8000}}};
  for(int mode=0;mode<3;++mode) {
    bool complete=false; unsigned failures=0;
    for(long i=0;i<30 && !complete;++i) {
      long before=live; allocations=0; fail_at=i;
      try {
        if(mode==0) { auto result=c::prepare_baryon_abundance(s); }
        if(mode==1) { auto result=owner.evaluate(std::array{1.}); }
        if(mode==2) { auto result=owner.evaluate_equilibrium(q,std::string_view(s.source_origin)); }
        complete=true;
      } catch(const std::bad_alloc &) { ++failures; }
      fail_at=-1;
      check(live==before,"allocation cleanup");
    }
    check(complete && failures>0,"every allocation failure reached");
  }
}
}
int main() {
  try { controls(); allocation_controls(); std::cout<<"baryon abundance checks="<<checks<<" failures=0\n"; }
  catch(const std::exception &e) { fail_at=-1; std::cerr<<e.what()<<'\n';return 1; }
}
