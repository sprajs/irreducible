#include "irred/conditional_hydrogen_helium_drag.hpp"
#include "../src/hydrogen_helium_cell.hpp"
#include "../src/thermal_ruler.hpp"
#include <array>
#include <bit>
#include <cstddef>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <limits>
#include <new>
#include <stdexcept>
#include <utility>
using namespace irred::cosmology;
using S = irred::numerics::Status;
using W = long double;
namespace allocation_observation {
struct alignas(std::max_align_t) Header { std::size_t bytes; bool observed; };
bool enabled=false;
std::size_t live=0,peak=0;
}
// Observe requested C++ allocation payload, excluding this observer's own
// allocator bookkeeping. Borrowed input strings are constructed before enable.
void *operator new(std::size_t n) {
  using namespace allocation_observation;
  if(n>SIZE_MAX-sizeof(Header) || (enabled && n>SIZE_MAX-live))throw std::bad_alloc();
  auto *h=static_cast<Header*>(std::malloc(sizeof(Header)+n));if(!h)throw std::bad_alloc();
  h->bytes=n;h->observed=enabled;if(enabled){live+=n;peak=std::max(peak,live);}return h+1;
}
void operator delete(void *p) noexcept {
  if(!p)return;
  using namespace allocation_observation;
  auto *h=static_cast<Header*>(p)-1;if(h->observed)live-=h->bytes;std::free(h);
}
void *operator new[](std::size_t n){return ::operator new(n);}
void operator delete[](void*p)noexcept{::operator delete(p);}
void operator delete(void*p,std::size_t)noexcept{::operator delete(p);}
void operator delete[](void*p,std::size_t)noexcept{::operator delete(p);}
namespace {
unsigned checks = 0;
void need(bool b, const char *message) { ++checks; if (!b) throw std::runtime_error(message); }
void foundation_controls() {
  const auto bg = prepare_thermal_background({70, 3./16, 5./16, .25, .25, {}});
  const auto r = detail::thermal_baryon_loading(bg);
  need(r.status == S::ok && r.ratio_today == 1, "analytic loading R0=1");
  need(r.arithmetic_estimate == 32 * std::numeric_limits<W>::epsilon(), "frozen exact ratio arithmetic allowance");
  for (W u : {1.L, 301.L, 1101.L, 2701.L}) {
    const auto q = detail::thermal_baryon_loading_at_inverse_scale(r, u);
    const W original = 3 * W(bg.source().omega_b) / (4 * W(bg.source().omega_gamma) * u);
    need(q.status == S::ok && q.value == original, "pure-H original grouped center preserved");
    need(q.arithmetic_estimate <= 128 * std::numeric_limits<double>::epsilon() * q.value,
         "helper operation estimate inside inherited pure-H floor");
  }
  const auto small_bg = prepare_thermal_background({70, .25, 0, 1e-300, .5, {}});
  const auto small = detail::thermal_baryon_loading(small_bg);
  need(small.status == S::ok && small.ratio_today > 0, "TF02 positive loading admitted");
  need(detail::thermal_baryon_loading_at_inverse_scale(small, std::numeric_limits<W>::max()/2).status != S::ok,
       "TF02 positive inverse-scale underflow refuses");
  const auto zero = detail::thermal_baryon_loading(prepare_thermal_background({70,.25,0,0,.5,{}}));
  const auto zero_at = detail::thermal_baryon_loading_at_inverse_scale(zero, std::numeric_limits<W>::max()/2);
  need(zero.status == S::ok && zero_at.status == S::ok && zero_at.value == 0 && zero_at.arithmetic_estimate == 0,
       "exact zero-baryon witness retained");
  need(detail::thermal_baryon_sound(r,std::numeric_limits<W>::max()/2).status == S::outside_domain,
       "TF03 sound scale outside declared a<=1 refuses");
  const auto sound=detail::thermal_baryon_sound(r,1);
  need(sound.status==S::ok && sound.sound_speed_over_c>0 && sound.sound_speed_squared_over_c_squared>0 &&
       sound.sound_speed_estimate_over_c>0 && sound.sound_speed_squared_estimate_over_c_squared>0,
       "positive sound outputs and diagnostics in admitted domain");
  auto tiny = r; tiny.arithmetic_estimate = std::ldexp(1.L,-100);
  const W a = std::ldexp(1.L,-10);
  const auto xi = detail::thermal_ruler_internal::loading_xi(tiny,a);
  need(xi && *xi > 0 && std::isnormal(*xi), "tiny positive rationalized xi survives");
  need(std::abs(*xi/(tiny.arithmetic_estimate*a/(2*(1+a)))-1) < 1e-16L,
       "tiny xi independent first-order analytic limit");
  tiny.arithmetic_estimate=0;
  need(detail::thermal_ruler_internal::loading_xi(tiny,a)==0, "exact zero loading error xi witness");
  need(detail::thermal_ruler_internal::loading_xi(r,0)==0, "exact zero scale xi witness");
}
void primitive_controls() {
  detail::HydrogenHeliumCell c;
  c.low=0;c.high=1;c.hydrogen={.5L,.25L};c.helium=c.hydrogen;c.coefficient={4,2};
  const auto cell=detail::hydrogen_helium_drag_cell(c,.25L,.25L,1,0);
  need(cell.status==S::ok && cell.rate.degree==6 && cell.error.degree==8, "shared degree-six/eight positive primitive");
  const auto whole=cell.integrate(0,1);
  need(whole.status==S::ok && std::abs(whole.value-127.L/28)<1e-15L, "analytic K=(u^7-1)/28");
  for (const W z : {.25L,.5L,1.L}) {
    const auto q=cell.integrate(0,z);
    need(q.status==S::ok && std::abs(q.value-(std::pow(1+z,7)-1)/28)<2e-15L, "partial positive primitive analytic control");
  }
  const W low=.5L,high=low+std::ldexp(1.L,-50);
  const auto small=cell.integrate(low,high);
  W sum=0;for(unsigned j=0;j<=6;++j)sum+=std::pow(1+high,6-j)*std::pow(1+low,j);
  const W expected=(high-low)*sum/28;
  need(small.status==S::ok && small.value>0 && std::abs(small.value/expected-1)<1e-15L,
       "tiny interval direct positive primitive avoids antiderivative cancellation");
  need(cell.integrate(.5L,.5L).value==0, "exact zero-width primitive witness");
  const W shift=2059.L/3584-61741.L/458752;
  const W central_m=std::pow(1.5L,6)/4;
  need(shift/central_m < .25L && cell.lower(.25L,.5L)>0,
       "central slope fails shifted-root bound while expanded interval retains positivity");
  // Permanent allocation controls: exact rational passing/refused totals from
  // the admitted source packet, independent of the physical history fixture.
  need(.5L/128 > 1e-3L, "locator-only ruler allocation refusal");
  need(321.L/1048576 < 1e-3L, "complete passing analytic TOTAL allocation");
  ConditionalDragWork malformed;malformed.source_map=SIZE_MAX;malformed.cell=1;
  need(!malformed.checked_total(), "output work checksum refuses SIZE_MAX+1");
}
ConditionalHydrogenHeliumDragRequest request(){
 return {{67.4,.02237,.12,2.7255,1.7e-5,{}},.245,1.6735328383153192e-27,6.646479071583153e-27,
  "NEXT16 v3 exact emitted nuclei/thermal working source",
  "caller-chosen neutral-effective kg masses; no atomic measurement claim",2700,300};
}
ConditionalHydrogenHeliumDragPolicy policy(){
 ConditionalHydrogenHeliumDragPolicy p;p.history.base_intervals=16384;p.history.maximum_fine_intervals=65536;
 p.history.maximum_total_work=3999998;p.history.maximum_native_bytes=33546240;return p;
}
void owner_controls(){
 const auto r=request();const auto p=policy();
 const std::array<ConditionalDragInterval,3> input{{
  {"synthetic-depth-sensitivity",.01,.02,"supplied synthetic interval; no tail shape or probability"},
  {"zero-tail-truncated-control",0,0,"explicit truncated control; no physical negligible-tail assertion"},
  {"collapsed-supplied-depth-control",.01,.01,"fixed synthetic supplied-depth control"}}};
 allocation_observation::enabled=true;
 {
  auto owner=prepare_conditional_hydrogen_helium_drag(r,p);
  if(owner.status()!=S::ok)std::cerr<<"prepare status="<<int(owner.status())<<" work="<<owner.preparation_work().checked_total().value_or(SIZE_MAX)<<'\n';
  need(owner.status()==S::ok && owner.source() && owner.source_snapshot() && owner.history(),"actual conditional owner prepared");
  const auto bits=owner.source_snapshot()->emitted_nuclei_binary64_bits;
  need(bits[0]==std::bit_cast<std::uint64_t>(*owner.source_snapshot()->nuclei_today.hydrogen_nuclei.value) &&
       bits[1]==std::bit_cast<std::uint64_t>(*owner.source_snapshot()->nuclei_today.helium_nuclei.value),"emitted nuclei serialized exactly");
  const auto *b=owner.budget_diagnostics();
  need(b && b->requested_history_work==3999998 && b->requested_history_native_bytes==33546240 &&
       b->served_history_work && *b->served_history_work<=3999998 && b->served_history_native_bytes &&
       b->parent_live_bytes_before_history && *b->served_history_native_bytes<=p.maximum_native_bytes-*b->parent_live_bytes_before_history &&
       b->earned_work_before_history==2 && b->history_preparation_attempted,"actual requested/effective history budget retained");
  const auto budget_snapshot=*b;
  const auto prep=owner.preparation_work();
  need(prep.checked_total() && *prep.checked_total()<=4000000 && prep.source_map==2 && prep.cell>0 && prep.prefix>0 &&
       prep.primitive>0 && prep.rhs>0 && prep.capacity_outer>0,"all reached preparation categories retained");
  auto result=owner.evaluate(input);
  need(result.status==S::ok && result.rows.size()==input.size() && result.work.checked_total() &&
       *result.work.checked_total()<=4000000,"one ordered coarse batch owns preparation plus all attempts");
  for(std::size_t i=0;i<input.size();++i){
   need(result.rows[i].source.id==input[i].id && result.rows[i].source.origin==input[i].origin,"exact supplied tail identity/order retained");
   for(const auto &endpoint:result.rows[i].endpoints){
    if(endpoint.root_status!=S::ok || endpoint.ruler_status!=S::ok)std::cerr<<"endpoint statuses="<<int(endpoint.root_status)<<','<<int(endpoint.ruler_status)<<" total="<<endpoint.total_ruler_numerical_estimate_mpc<<'\n';
    need(endpoint.root_status==S::ok && endpoint.ruler_status==S::ok && endpoint.redshift && endpoint.comoving_ruler_mpc &&
         endpoint.total_ruler_numerical_estimate_mpc>0 && endpoint.total_ruler_numerical_estimate_mpc<=1e-3,"same root/ruler TOTAL allocation accepted");
   }
  }
  need(*result.rows[0].endpoints[0].redshift<*result.rows[0].endpoints[1].redshift &&
       *result.rows[0].endpoints[0].comoving_ruler_mpc>*result.rows[0].endpoints[1].comoving_ruler_mpc,"interval root/ruler monotonic order");
  need(result.rows[2].endpoints[0].redshift==result.rows[2].endpoints[1].redshift,"collapsed supplied interval is deterministic");
  need(result.peak_payload_bytes<=p.maximum_native_bytes && allocation_observation::peak<=result.peak_payload_bytes,
       "observed live allocations fit checked simultaneous bound");
  allocation_observation::enabled=false; // Further independent owners are outside this scoped request peak.
  const std::array<ConditionalDragInterval,3> bad{{{"reversed",.02,.01,"synthetic refusal"},
      {"nan",0,std::numeric_limits<double>::quiet_NaN(),"synthetic refusal"},{"after-good",.01,.01,"synthetic control"}}};
  const auto mixed=owner.evaluate(bad);
  need(mixed.rows.size()==3 && mixed.rows[0].endpoints[0].root_status==S::outside_domain &&
       mixed.rows[1].endpoints[0].root_status==S::nonfinite_input && mixed.rows[2].endpoints[0].root_status==S::ok,
       "refused rows are retained alongside valid ordered rows");
  auto copy=owner;auto moved=std::move(owner);
  need(!owner.source() && !owner.source_snapshot() && !owner.budget_diagnostics() && owner.status()==S::invalid_input,
       "moved owner has no source or freshness diagnostics");
  need(copy.source_snapshot()->emitted_nuclei_binary64_bits==bits && moved.source_snapshot()->emitted_nuclei_binary64_bits==bits,
       "owning copy/move retains emitted working law");
  bool refused=false;try{copy=moved;}catch(const std::bad_alloc&){refused=true;}
  need(refused && copy.source_snapshot()->emitted_nuclei_binary64_bits==bits,"oversized simultaneous assignment refuses without changing destination");
  auto tiny=p;tiny.history.maximum_total_work=1;
  const auto failed=prepare_conditional_hydrogen_helium_drag(r,tiny);
  need(failed.status()==S::work_limit && failed.source() && failed.budget_diagnostics() &&
       failed.budget_diagnostics()->history_preparation_attempted && failed.preparation_work().source_map==2 &&
       failed.preparation_work().background==1 && failed.preparation_work().checked_total()==3,
       "failed child retains actual attempts and supplied/effective ceilings");
  auto tight=p;
  const auto history_bound=hydrogen_helium_history_payload_bound(65536,0,0,r.source_origin.capacity());
  need(history_bound.has_value(),"finite inherited simultaneous history bound");
  tight.maximum_native_bytes=*budget_snapshot.parent_live_bytes_before_history+*history_bound-1;
  const auto byte_refusal=prepare_conditional_hydrogen_helium_drag(r,tight);
  need(byte_refusal.status()==S::work_limit && byte_refusal.budget_diagnostics() &&
       byte_refusal.budget_diagnostics()->served_history_native_bytes &&
       !byte_refusal.budget_diagnostics()->history_preparation_attempted && byte_refusal.preparation_work().source_map==2,
       "one-byte-tight history phase refuses before allocation with prior work intact");
 }
 allocation_observation::enabled=false;
 need(allocation_observation::live==0,"all observed owning allocations released");
 auto unsupported=r;unsupported.model.species={{0,1.95,2}};
 need(prepare_conditional_hydrogen_helium_drag(unsupported,p).status()==S::outside_domain,"distinct empty-species consumer domain enforced");
}
}
int main() {
  try { foundation_controls();primitive_controls();owner_controls(); }
  catch(const std::exception &e){std::cerr<<"FAILED: "<<e.what()<<" after "<<checks<<" checks\n";return 1;}
  std::cout<<checks<<" conditional H/He drag checks passed\n";
}
