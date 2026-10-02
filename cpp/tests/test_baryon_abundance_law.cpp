#include "irred/baryon_abundance_law.hpp"
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
#if defined(__GNUC__) || defined(__clang__)
#define NOINLINE __attribute__((noinline))
#else
#define NOINLINE
#endif
NOINLINE void *operator new(std::size_t n) {
  if (fail_at>=0 && allocations++==fail_at) throw std::bad_alloc();
  if (void *p=std::malloc(n?n:1)) { ++live; return p; }
  throw std::bad_alloc();
}
NOINLINE void *operator new[](std::size_t n) { return ::operator new(n); }
NOINLINE void operator delete(void *p) noexcept { if(p) { --live; std::free(p); } }
NOINLINE void operator delete[](void *p) noexcept { ::operator delete(p); }
NOINLINE void operator delete(void *p,std::size_t) noexcept { ::operator delete(p); }
NOINLINE void operator delete[](void *p,std::size_t) noexcept { ::operator delete(p); }
namespace {
namespace c=irred::cosmology; namespace a=irred::atomic;
using S=irred::numerics::Status; using W=long double;
unsigned checks=0;
void need(bool b,const char *m) { ++checks; if(!b) throw std::runtime_error(m); }
void near(W x,W y,W tolerance,const char *m) {
  need(std::isfinite(x) && std::isfinite(y) && (y==0 ? x==0 : std::abs(x/y-1)<=tolerance),m);
}
std::uint64_t bits(double x) { return std::bit_cast<std::uint64_t>(x); }
struct Fixture {
  std::array<c::BaryonAbundanceLawRow,2> rows{{{"earlier",1e-7},{"later",1.5e-7}}};
  std::array<std::array<double,2>,2> temperature{{{10000,10000},{18000,18000}}};
  std::array<c::BaryonAbundanceLawState,2> states;
  Fixture() {
    c::BaryonAbundanceSource source{.0224,.25,1.67353284e-27,6.64647907e-27,
        std::string(90,'s'),std::string(90,'m')};
    states[0]={"lower",1,source,temperature[0],std::string_view("supplied matter T lower")};
    source.physical_baryon_density*=2;
    states[1]={"upper",3,source,temperature[1],std::string_view("supplied matter T upper")};
  }
  c::BaryonAbundanceLawInput input() const { return {rows,states,"synthetic finite joint law","same state shared across all rows and outputs"}; }
};
void density_controls() {
  Fixture f; auto op=c::prepare_baryon_abundance_law(f.input());
  need(op.status()==S::ok && op.source()->state_preparations==2 && op.source()->normalization_terms==2,"one normalization and one preparation per state");
  near(op.source()->states[0].probability,.25L,0,"rational first weight");
  near(op.source()->states[1].probability,.75L,0,"rational second weight");
  auto result=op.evaluate_density();
  need(result.status==S::ok && result.moments && result.axes.size()==4 && result.attempts.size()==2,"full density joint law");
  need(result.source.get()==op.source().get() && result.work.mapping_rows==4 && result.work.solves==0 && result.work.covariance_products==20,"shared lineage and actual density work");
  const auto &m=*result.moments;
  std::array<W,4> v{};
  for(std::size_t r=0;r<2;++r) {
    v[2*r]=*result.attempts[0].rows[r].hydrogen_nuclei.value;
    v[2*r+1]=*result.attempts[0].rows[r].helium_nuclei.value;
    need(result.axes[2*r].row_index==r && result.axes[2*r].output==c::BaryonAbundanceLawOutput::hydrogen_nuclei &&
        result.axes[2*r+1].output==c::BaryonAbundanceLawOutput::helium_nuclei,"density axis order");
  }
  // Exact rational finite-law algebra, independent of moment accumulation.
  // Mapping against the same owner's values is transport evidence; the peer
  // uses independently derived exact-binary critical-density/charge facts.
  for(std::size_t i=0;i<4;++i) {
    near(m.mean[i],7*v[i]/4,2e-15L,"rational mean7/4");
    need(m.mean_empirical_sensitivity[i]>0 && m.mean_empirical_sensitivity[i]<=1e-12*m.mean[i],"mean empirical allowance");
    for(std::size_t j=0;j<4;++j) {
      near(m.covariance[i*4+j],3*v[i]*v[j]/16,2e-14L,"rational full covariance3/16");
      need(bits(m.covariance[i*4+j])==bits(m.covariance[j*4+i]),"exact stored symmetry");
      need(m.covariance_empirical_sensitivity[i*4+j]>0 &&
          W(m.covariance_empirical_sensitivity[i*4+j])<=1e-8L*std::sqrt(W(m.covariance[i*4+i]))*std::sqrt(W(m.covariance[j*4+j])),"covariance empirical allowance");
    }
  }
  near(m.mean[2]/m.mean[0],std::pow(W(f.rows[0].scale_factor)/f.rows[1].scale_factor,3),2e-15L,"a^-3 mean scaling");
  near(m.covariance[2*4+2]/m.covariance[0],std::pow(W(f.rows[0].scale_factor)/f.rows[1].scale_factor,6),2e-15L,"a^-6 variance scaling");
  auto perm=f.states; std::swap(perm[0],perm[1]); auto in=f.input();in.states=perm;
  auto swapped=c::prepare_baryon_abundance_law(in).evaluate_density();
  need(swapped.status==S::ok,"state permutation admitted");
  for(std::size_t i=0;i<4;++i) near(swapped.moments->mean[i],m.mean[i],2e-15L,"state permutation means");
  for(std::size_t i=0;i<16;++i) near(swapped.moments->covariance[i],m.covariance[i],2e-14L,"state permutation covariance");
  f.states[0].relative_mass*=1024;f.states[1].relative_mass*=1024;
  auto weighted=c::prepare_baryon_abundance_law(f.input()).evaluate_density();
  need(weighted.status==S::ok,"common weight scaling");
  for(std::size_t i=0;i<16;++i) need(bits(weighted.moments->covariance[i])==bits(m.covariance[i]),"power-two weight scaling preserves bits");
  f.states[0].relative_mass=1;f.states[1].relative_mass=3;
  std::array<c::BaryonAbundanceLawState,4> duplicated{f.states[0],f.states[0],f.states[1],f.states[1]};
  const std::array<std::string_view,4> ids{"lower0","lower1","upper0","upper1"};
  for(std::size_t i=0;i<4;++i) { duplicated[i].id=ids[i];duplicated[i].relative_mass/=2; }
  in=f.input();in.states=duplicated;
  auto dup=c::prepare_baryon_abundance_law(in).evaluate_density();
  need(dup.status==S::ok && dup.attempts.size()==4,"duplicated support retains every state");
  for(std::size_t i=0;i<16;++i) near(dup.moments->covariance[i],m.covariance[i],2e-14L,"duplicated support same law");
  f.states[1].abundance=f.states[0].abundance;
  auto constant=c::prepare_baryon_abundance_law(f.input()).evaluate_density();
  need(constant.status==S::ok && constant.work.covariance_products==0,"T-only variation witnesses density constancy");
  for(double x:constant.moments->covariance) need(x==0,"witnessed constant density covariance");
  for(double x:constant.moments->covariance_empirical_sensitivity) need(x==0,"constant dependence cancels shared deterministic diagnostic");
  f.states[0].abundance.helium4_mass_fraction=.25;f.states[1].abundance.helium4_mass_fraction=.75;
  auto anti=c::prepare_baryon_abundance_law(f.input()).evaluate_density();
  need(anti.status==S::ok && anti.moments->covariance[1]<0,"Y induces H/He anticorrelation");
  near(anti.moments->covariance[1],-std::sqrt(W(anti.moments->covariance[0]))*std::sqrt(W(anti.moments->covariance[5])),2e-14L,"rank-one anticorrelation no jitter");
  f.states[0].abundance.helium4_mass_fraction=.25;f.states[1].abundance=f.states[0].abundance;
  f.states[1].abundance.hydrogen1_effective_mass_kg*=2;f.states[1].abundance.helium4_effective_mass_kg*=2;
  auto mass=c::prepare_baryon_abundance_law(f.input()).evaluate_density();
  need(mass.status==S::ok,"supplied mass law");
  near(mass.moments->mean[0],5*v[0]/8,2e-15L,"inverse mass law mean5/8");
  near(mass.moments->covariance[0],3*v[0]*v[0]/64,2e-14L,"inverse mass law variance3/64");
  f.states[0].abundance.helium4_mass_fraction=0;f.states[1].abundance.helium4_mass_fraction=0;
  auto endpoint=c::prepare_baryon_abundance_law(f.input()); auto mapped=endpoint.evaluate_density();
  need(mapped.status==S::ok && mapped.moments->mean[1]==0 && mapped.moments->covariance[5]==0,"exact absent helium density axis");
  auto refused=endpoint.evaluate_equilibrium();
  need(refused.status==S::outside_domain && !refused.moments && refused.attempts.size()==2 && refused.work.solves==0,"absent species LTE preserves all attempts");
  for(const auto &batch:refused.attempts) for(const auto &row:batch.rows)
    need(row.nuclei.helium_nuclei.value && *row.nuclei.helium_nuclei.value==0 &&
        row.equilibrium.electron_density.status==S::outside_domain,"mapping survives LTE absent-species refusal");
}
void lte_controls() {
  Fixture f; auto op=c::prepare_baryon_abundance_law(f.input());
  for(unsigned mask=1;mask<64;++mask) {
    auto result=op.evaluate_equilibrium(mask);
    need(result.status==S::ok && result.moments && result.attempts.size()==2,"all63 LTE masks");
    need(result.axes.size()==2*std::size_t(std::popcount(mask)) && result.work.solves==4 && result.work.mapping_rows==4,"mask axes and one solve per state/row");
    std::size_t iterations=0,evaluations=0;
    for(const auto &batch:result.attempts) {
      iterations+=batch.root_iterations;evaluations+=batch.charge_evaluations;
      for(const auto &row:batch.rows) {
        const auto &e=row.equilibrium;
        std::array<const a::HydrogenHeliumValue*,6> groups{&e.hydrogen_neutral,&e.hydrogen_ionized,&e.helium_neutral,&e.helium_singly_ionized,&e.helium_doubly_ionized,&e.electron_density};
        for(unsigned j=0;j<6;++j) need(mask&(1u<<j) ? groups[j]->value.has_value() : !groups[j]->value && groups[j]->relative_arithmetic_estimate==0,"omitted native group is not zero");
      }
    }
    need(iterations==result.work.root_iterations && evaluations==result.work.charge_evaluations,"exact child counters retained");
  }
  auto in=f.input();in.states=std::span(f.states).first(1);
  auto one=c::prepare_baryon_abundance_law(in).evaluate_equilibrium();
  need(one.status==S::ok,"one-state LTE law");
  for(double x:one.moments->covariance) need(x==0,"one-state exact zero law spread");
  // Original a={1e-7,1.5e-7}, T=1 K failed the native dilute-gas gate
  // before solving (preserved failure-01/probe receipt). Low-density fixed
  // a={1,.5} keeps this control inside that gate without changing any budget.
  f.rows[0].scale_factor=1;f.rows[1].scale_factor=.5;
  f.temperature[0]={1,1};f.temperature[1]={1,1};
  auto cold=c::prepare_baryon_abundance_law(f.input());
  auto neutral=cold.evaluate_equilibrium(a::hydrogen_neutral_fraction|a::helium_neutral_fraction);
  need(neutral.status==S::conditioning_budget_exceeded && !neutral.moments && neutral.attempts.size()==2,"varying inputs equal rounded neutral outputs do not witness zero spread");
  for(const auto &batch:neutral.attempts) for(const auto &row:batch.rows)
    need(row.equilibrium.hydrogen_neutral.value && *row.equilibrium.hydrogen_neutral.value==1 &&
        row.equilibrium.helium_neutral.value && *row.equilibrium.helium_neutral.value==1,"neutral groups survive charged underflow independently");
  in=f.input();in.states=std::span(f.states).first(1);
  auto cold_one=c::prepare_baryon_abundance_law(in);
  auto only_neutral=cold_one.evaluate_equilibrium(a::hydrogen_neutral_fraction|a::helium_neutral_fraction);
  need(only_neutral.status==S::ok,"cold neutral-only requested law available");
  auto charged=cold_one.evaluate_equilibrium();
  need(charged.status==S::conditioning_budget_exceeded && !charged.moments && charged.attempts.size()==1 &&
      charged.attempts[0].rows[0].equilibrium.hydrogen_neutral.value &&
      !charged.attempts[0].rows[0].equilibrium.electron_density.value,"required charged failure withholds aggregate but retains available partners");
  f.temperature[0]={5e5,5e5};f.temperature[1]={18000,18000};
  auto outside=c::prepare_baryon_abundance_law(f.input());
  need(outside.status()==S::ok && outside.evaluate_density().status==S::ok,"normal positive T outside LTE support leaves density useful");
  auto failed=outside.evaluate_equilibrium();
  need(failed.status==S::outside_domain && failed.attempts.size()==2 && failed.attempts[0].rows.size()==2 &&
      failed.attempts[1].rows.size()==2 && failed.work.solves==2,"unsupported T does not omit good later states");
  need(failed.attempts[1].rows[0].equilibrium.electron_density.value.has_value(),"later good state remains available");
}
void resource_controls() {
  Fixture f; auto op=c::prepare_baryon_abundance_law(f.input());
  c::BaryonAbundanceLawPolicy p;p.maximum_solves=1;
  auto solves=op.evaluate_equilibrium(a::electron_number_density,p);
  need(solves.status==S::work_limit && !solves.moments && solves.attempts.size()==2 && solves.work.mapping_rows==4 && solves.work.solves==1,"global solve cap across states");
  need(solves.attempts[0].rows[0].equilibrium.electron_density.value &&
      solves.attempts[1].rows[0].equilibrium.electron_density.status==S::work_limit,"prior success survives remaining zero solve quota");
  p={};p.maximum_charge_evaluations=1;
  auto charges=op.evaluate_equilibrium(a::electron_number_density,p);
  need(charges.status==S::work_limit && charges.attempts.size()==2 && charges.work.solves==4 &&
      charges.work.charge_evaluations==1 && charges.work.root_iterations==0,"failed bracket work charged once across all calls");
  p={};p.maximum_root_iterations=0;
  auto iterations=op.evaluate_equilibrium(a::electron_number_density,p);
  need(iterations.status==S::work_limit && iterations.work.solves==4 && iterations.work.charge_evaluations==0,"zero iteration quota retains attempted solves");
  auto mapped=op.evaluate_density();const auto bound=mapped.work.combined_native_payload_bound;
  p={};p.maximum_native_bytes=bound-1;
  auto bytes=op.evaluate_density(p);
  need(bytes.status==S::work_limit && bytes.attempts.empty() && !bytes.moments,"combined density payload gate before work");
  p.maximum_native_bytes=bound;
  need(op.evaluate_density(p).status==S::ok,"exact density combined payload threshold");
  auto full=op.evaluate_equilibrium();const auto lte_bound=full.work.combined_native_payload_bound;
  p={};p.maximum_native_bytes=lte_bound-1;
  need(op.evaluate_equilibrium(63,p).status==S::work_limit,"combined LTE payload includes child peak");
  p.maximum_native_bytes=lte_bound;
  need(op.evaluate_equilibrium(63,p).status==S::ok,"exact LTE combined payload threshold");
  for(int limit=0;limit<5;++limit) {
    p={};
    if(limit==0)p.maximum_states=1;
    if(limit==1)p.maximum_rows=1;
    if(limit==2)p.maximum_state_rows=3;
    if(limit==3)p.maximum_axes=3;
    if(limit==4)p.maximum_covariance_products=19;
    auto result=op.evaluate_density(p);
    need(result.status==S::work_limit && result.attempts.empty(),"combined dimension/product quota before work");
  }
  p={};p.maximum_states=1025;
  need(op.evaluate_density(p).status==S::invalid_input,"hard policy cap");
  p={};p.mean_relative_sensitivity=2e-12;
  need(op.evaluate_density(p).status==S::invalid_input,"mean allocation cannot be silently weakened");
  p={};p.covariance_relative_sensitivity=2e-8;
  need(op.evaluate_density(p).status==S::invalid_input,"covariance allocation cannot be silently weakened");
  need(op.evaluate_equilibrium(0).status==S::invalid_input && op.evaluate_equilibrium(64).status==S::invalid_input,"invalid masks");
  p={};p.maximum_native_bytes=op.source()->retained_native_payload_bytes-1;
  need(c::prepare_baryon_abundance_law(f.input(),p).status()==S::work_limit,"prepared source payload limit");
  std::fesetround(FE_DOWNWARD);
  auto rounding=op.evaluate_density();auto preparation=c::prepare_baryon_abundance_law(f.input());
  std::fesetround(FE_TONEAREST);
  need(rounding.status==S::invalid_input && rounding.attempts.size()==2 && preparation.status()==S::invalid_input,"strict nearest profile");
}
void domain_and_lifetime_controls() {
  Fixture f; auto in=f.input();
  in.rows={};need(c::prepare_baryon_abundance_law(in).status()==S::invalid_input,"empty rows");
  in=f.input();in.states={};need(c::prepare_baryon_abundance_law(in).status()==S::invalid_input,"empty states");
  in=f.input();in.dependence_origin={};need(c::prepare_baryon_abundance_law(in).status()==S::invalid_input,"missing dependence provenance");
  f.states[1].id=f.states[0].id;need(c::prepare_baryon_abundance_law(f.input()).status()==S::invalid_input,"duplicate state IDs");
  f.states[1].id="upper";f.rows[1].id=f.rows[0].id;
  need(c::prepare_baryon_abundance_law(f.input()).status()==S::invalid_input,"duplicate row IDs");f.rows[1].id="later";
  f.states[0].temperature_kelvin=std::span(f.temperature[0]).first(1);
  need(c::prepare_baryon_abundance_law(f.input()).status()==S::invalid_input,"temperature order shape");
  f.states[0].temperature_kelvin=f.temperature[0];
  for(double w:{0.,-1.,std::numeric_limits<double>::denorm_min(),std::numeric_limits<double>::infinity()}) {
    f.states[0].relative_mass=w;
    need(c::prepare_baryon_abundance_law(f.input()).status()!=S::ok,"illegal probability mass");
  }
  f.states[0].relative_mass=1;
  for(double t:{0.,-1.,std::numeric_limits<double>::denorm_min(),std::numeric_limits<double>::quiet_NaN()}) {
    f.temperature[0][0]=t;
    auto invalid=c::prepare_baryon_abundance_law(f.input());
    need(invalid.status()!=S::ok && invalid.source() && invalid.source()->states.size()==2,"invalid T preserves admitted source lineage");
  }
  f.temperature[0][0]=10000;f.states[1].abundance.helium4_mass_fraction=1.1;
  auto invalid=c::prepare_baryon_abundance_law(f.input());
  need(invalid.status()==S::outside_domain && invalid.source()->states[1].abundance.source() &&
      !invalid.evaluate_density().moments,"invalid abundance state never dropped");
  f.states[1].abundance.helium4_mass_fraction=.25;
  auto owner=c::prepare_baryon_abundance_law(f.input());auto result=owner.evaluate_equilibrium(a::electron_number_density);
  const auto original_omega=owner.source()->states[0].abundance.source()->physical_baryon_density;
  f.temperature[0][0]=999;f.states[0].abundance.source_origin="changed";f.states[0].abundance.physical_baryon_density=100;
  need(owner.source()->states[0].temperature_kelvin[0]==10000 && owner.source()->states[0].abundance.source()->source_origin.size()==90 &&
      owner.source()->states[0].abundance.source()->physical_baryon_density==original_omega,"source owned after caller mutation");
  need(bits(result.attempts[0].rows[0].source.scale_factor)==bits(f.rows[0].scale_factor) &&
      bits(result.attempts[0].rows[0].source.temperature_kelvin)==bits(10000.),"returned source bits and ordering");
  auto copy=owner;auto moved=std::move(copy);
  need(copy.status()==S::invalid_input && !copy.source() && moved.source().get()==owner.source().get(),"move invalidation and immutable shared source");
  auto *self=&moved;moved=std::move(*self);need(moved.status()==S::ok,"self move preserves owner");
  copy=owner;copy=copy;need(copy.status()==S::ok,"copy self preserves owner");
  moved=c::BaryonAbundanceLaw{};owner=c::BaryonAbundanceLaw{};copy=c::BaryonAbundanceLaw{};
  need(result.source->states[0].abundance.source()->source_origin.size()==90 &&
      result.source->rows[0].id=="earlier" && result.attempts[0].rows[0].equilibrium.electron_density.value,"result retains full lineage beyond owner lifetime");
}
void representation_controls() {
  Fixture f;
  f.states[0].relative_mass=std::numeric_limits<double>::min();
  f.states[1].relative_mass=std::numeric_limits<double>::max();
  f.states[0].abundance.physical_baryon_density=f.states[1].abundance.physical_baryon_density=0;
  auto wide=c::prepare_baryon_abundance_law(f.input());
  need(wide.status()==S::ok && wide.source()->states[0].probability>0 &&
      double(wide.source()->states[0].probability)==0,"tiny probability retained beyond binary64 range");
  auto zero=wide.evaluate_density();need(zero.status==S::ok && zero.attempts.size()==2,"tiny-probability state remains an attempt");
  for(double x:zero.moments->mean) need(x==0,"only mathematical absent density is zero");
  f.states[0].relative_mass=1;f.states[1].relative_mass=1e-306;
  f.rows[0].scale_factor=1;
  for(auto &s:f.states) {
    s.temperature_kelvin=s.temperature_kelvin.first(1);
    s.abundance.helium4_mass_fraction=0;
    s.abundance.hydrogen1_effective_mass_kg=1e-22;
  }
  f.states[0].abundance.physical_baryon_density=0;f.states[1].abundance.physical_baryon_density=1;
  auto in=f.input();in.rows=in.rows.first(1);
  auto subnormal=c::prepare_baryon_abundance_law(in).evaluate_density();
  need(subnormal.status==S::ok && subnormal.moments &&
      subnormal.moments->mean[0]>0 && !std::isnormal(subnormal.moments->mean[0]) &&
      subnormal.moments->covariance[0]>0 && !std::isnormal(subnormal.moments->covariance[0]),"positive subnormal mean and covariance admitted with actual cast loss");
  need(subnormal.moments->mean_empirical_sensitivity[0]>0 &&
      W(subnormal.moments->mean_empirical_sensitivity[0])<=1e-12L*subnormal.moments->mean[0] &&
      subnormal.moments->covariance_empirical_sensitivity[0]>0 &&
      W(subnormal.moments->covariance_empirical_sensitivity[0])<=1e-8L*subnormal.moments->covariance[0],"subnormal diagnostics remain nonzero and within allocations");
  f.states[1].abundance.hydrogen1_effective_mass_kg=1e20;
  auto mean_underflow=c::prepare_baryon_abundance_law(in).evaluate_density();
  need(mean_underflow.status==S::conditioning_budget_exceeded && !mean_underflow.moments &&
      mean_underflow.attempts[1].rows[0].hydrogen_nuclei.value,"positive mean underflow refuses rather than zero");
  f.states[0].relative_mass=1;f.states[1].relative_mass=3;
  f.states[0].abundance.physical_baryon_density=1;f.states[1].abundance.physical_baryon_density=2;
  f.states[0].abundance.hydrogen1_effective_mass_kg=f.states[1].abundance.hydrogen1_effective_mass_kg=1e253;
  auto covariance_underflow=c::prepare_baryon_abundance_law(in).evaluate_density();
  need(covariance_underflow.status==S::conditioning_budget_exceeded && !covariance_underflow.moments &&
      covariance_underflow.attempts[0].rows[0].hydrogen_nuclei.value,"positive covariance underflow refuses rather than zero spread");
}
void string_envelope_controls() {
  for(std::size_t length=16;length<=29;++length) {
    Fixture f;
    const std::string origin(length,'p');
    std::array<std::string,2> row_ids{std::string(length,'r'),std::string(length,'r')};
    std::array<std::string,2> state_ids{std::string(length,'s'),std::string(length,'s')};
    row_ids[1].back()='b';state_ids[1].back()='b';
    for(std::size_t i=0;i<2;++i) {
      f.rows[i].id=row_ids[i];f.states[i].id=state_ids[i];
      f.states[i].temperature_origin=origin;
      f.states[i].abundance.source_origin=std::string(length,'a');
      f.states[i].abundance.mass_origin=std::string(length,'m');
    }
    auto in=f.input();in.distribution_origin=origin;in.dependence_origin=origin;
    const std::size_t base=sizeof(c::BaryonAbundanceLaw)+sizeof(c::BaryonAbundanceLawSource)+
        2*sizeof(c::BaryonAbundanceLawRetainedState)+2*sizeof(c::BaryonAbundanceLawOwnedRow)+4*sizeof(double),
        envelope=std::max<std::size_t>(32,length+1), correct=base+12*envelope;
    // Before the fix, this admitted allocation then discovered excess actual
    // capacity. Rejection must precede the FIRST observed allocation.
    c::BaryonAbundanceLawPolicy p;
    p.maximum_native_bytes=base+12*(std::max(length,std::string{}.capacity())+1);
    allocations=0;fail_at=std::numeric_limits<long>::max();
    auto refused=c::prepare_baryon_abundance_law(in,p);
    fail_at=-1;
    need(refused.status()==S::work_limit && !refused.source() && allocations==0,"string preparation envelope refuses before allocation at lengths16..29");
    p.maximum_native_bytes=correct;
    auto op=c::prepare_baryon_abundance_law(in,p);
    need(op.status()==S::ok && op.source()->preparation_native_payload_bound==correct &&
        op.source()->retained_native_payload_bytes<=correct,"conservative admission and actual retained string capacities are distinct");
    auto lte=op.evaluate_equilibrium();need(lte.status==S::ok,"string-envelope LTE control admitted");
    const std::size_t axes=12,cells=axes*axes,values=2*axes,
        peak=op.source()->retained_native_payload_bytes+sizeof(c::BaryonAbundanceLteLawResult)+
        axes*sizeof(c::BaryonAbundanceLawAxis)+2*sizeof(c::BaryonAbundanceLteBatch)+
        4*sizeof(c::BaryonAbundanceLteRow)+2*(axes+cells)*sizeof(double)+
        2*(values+axes+cells)*sizeof(W)+axes*sizeof(unsigned char)+
        2*(sizeof(c::BaryonAbundanceLteQuery)+sizeof(a::HydrogenHeliumState)+sizeof(a::HydrogenHeliumRow))+
        sizeof(a::HydrogenHeliumBatch)+2*envelope;
    need(lte.work.combined_native_payload_bound==peak,"child LTE origin copied-string envelope charged in peak payload");
    p={};p.maximum_native_bytes=peak-1;
    allocations=0;fail_at=std::numeric_limits<long>::max();
    auto rejected=op.evaluate_equilibrium(63,p);
    fail_at=-1;
    need(rejected.status==S::work_limit && rejected.attempts.empty() && allocations==0,"LTE copied-string envelope refuses before first allocation");
    p.maximum_native_bytes=peak;
    need(op.evaluate_equilibrium(63,p).status==S::ok,"exact conservative copied-string LTE threshold");
    for(const auto &batch:lte.attempts) need(batch.temperature_origin.capacity()+1<=envelope,"observed child string capacity fits pre-allocation envelope");
  }
  Fixture f;auto in=f.input();in.distribution_origin=std::string_view("x",std::numeric_limits<std::size_t>::max());
  need(c::prepare_baryon_abundance_law(in).status()==S::work_limit,"string envelope overflow refused without character scan");
}
void allocation_controls() {
  Fixture f;auto op=c::prepare_baryon_abundance_law(f.input());
  for(int mode=0;mode<3;++mode) {
    bool complete=false;unsigned failures=0;
    for(long i=0;i<200 && !complete;++i) {
      const long before=live;allocations=0;fail_at=i;
      try {
        if(mode==0) { auto result=c::prepare_baryon_abundance_law(f.input()); }
        if(mode==1) { auto result=op.evaluate_density(); }
        if(mode==2) { auto result=op.evaluate_equilibrium(); }
        complete=true;
      } catch(const std::bad_alloc &) { ++failures; }
      fail_at=-1;need(live==before,"every allocation failure restores live allocations");
    }
    need(complete && failures>0,"all allocation sites reached through success");
    std::cout<<"allocation_mode="<<mode<<" failures="<<failures<<" live="<<live<<'\n';
  }
}
}
int main() {
  try {
    density_controls();lte_controls();resource_controls();domain_and_lifetime_controls();
    representation_controls();string_envelope_controls();allocation_controls();
    std::cout<<"baryon abundance law checks="<<checks<<" failures=0\n";
  } catch(const std::exception &e) { fail_at=-1;std::fesetround(FE_TONEAREST);std::cerr<<e.what()<<'\n';return 1; }
}
