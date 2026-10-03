#include "irred/gaussian_simulation.hpp"
#include "irred/detector.hpp"
#include "gaussian_recovery_controls.hpp"
#include <sstream>
#include <algorithm>
#include <atomic>
#include <cfenv>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <new>
#include <stdexcept>
namespace {
 std::atomic<bool> armed=false,measure_allocations=false;
 std::atomic<size_t> calls=0,fail_on=0,live_bytes=0,peak_bytes=0;
 struct alignas(std::max_align_t) Allocation {size_t bytes;bool measured;};
}
#if defined(__GNUC__) || defined(__clang__)
#define NOINLINE __attribute__((noinline))
#else
#define NOINLINE
#endif
NOINLINE void *operator new(size_t n) {
 if(armed&&++calls==fail_on)throw std::bad_alloc();
 const size_t amount=n?n:1;if(amount>SIZE_MAX-sizeof(Allocation))throw std::bad_alloc();
 auto* h=static_cast<Allocation*>(std::malloc(amount+sizeof(Allocation)));if(!h)throw std::bad_alloc();
 h->bytes=amount;h->measured=measure_allocations;
 if(h->measured){const size_t now=live_bytes.fetch_add(amount)+amount;peak_bytes.store(std::max(peak_bytes.load(),now));}
 return h+1;
}
NOINLINE void *operator new[](size_t n){return ::operator new(n);}
NOINLINE void operator delete(void*p)noexcept{if(!p)return;auto*h=static_cast<Allocation*>(p)-1;if(h->measured)live_bytes.fetch_sub(h->bytes);std::free(h);}
NOINLINE void operator delete[](void*p)noexcept{::operator delete(p);}
NOINLINE void operator delete(void*p,size_t)noexcept{::operator delete(p);}
NOINLINE void operator delete[](void*p,size_t)noexcept{::operator delete(p);}
namespace {
namespace s=irred::statistics;namespace n=irred::numerics;
unsigned checks=0;void need(bool b,const char*w){++checks;if(!b)throw std::runtime_error(w);}
s::Gaussian source(){s::Metadata m;m.ordered_ids={"a","b"};m.measure="product d(mag)";
 m.source_semantics="synthetic controls";m.table_identity="synthetic covariance";m.ordering_provenance="exact declared source axes";
 return s::prepare_gaussian(std::array<double,4>{4,1,1,1.25},s::MatrixKind::covariance,m,100,1e-10,n::Arithmetic::longdouble_cpu_v1);}
s::GeneratingMean mean(){return {{1,2},{"a","b"},{"mag","mag"},"fixed generating mean","product d(mag)","original synthetic ideal Gaussian"};}
void run(){auto g=source();auto m=mean();const std::array<irred::random::Address,3>a{{{30,5},{30,7},{40,5}}};
 s::GaussianSimulationPolicy p;p.outputs=7;auto b=s::GaussianSimulation::simulate(g,m,a,17,p);
 need(b.status==n::Status::ok&&b.rows.size()==3&&b.values.size()==6&&b.words.size()==6,"coarse pooled vector output");
 for(auto&r:b.rows)need(r.status==n::Status::ok,"all generated rows finite");
 for(size_t i=0;i<3;++i){auto c=s::GaussianSimulation::simulate(g,m,std::span(a.data()+i,1),17,p);
  for(size_t j=0;j<2;++j){need(c.values[j]==b.values[2*i+j]&&c.words[j]==b.words[2*i+j],"exact shard replay");
   need(b.words[2*i+j]==irred::detector::random_words(17,{a[i].stream+j,a[i].sample}),"shared detector words");}}
 const std::array<irred::random::Address,2> overlap{{{30,5},{31,5}}};
 need(s::GaussianSimulation::simulate(g,m,overlap,17,p).rows.empty(),"coordinate interval overlap refused");
 const std::array<irred::random::Address,1> overflow{{{UINT64_MAX,5}}};
 need(s::GaussianSimulation::simulate(g,m,overflow,17,p).rows.empty(),"counter overflow refused");
 // Shape must be rejected before metadata payload traversal or byte admission.
 for(unsigned field=0;field<3;++field)for(size_t size:{0u,1u,3u,1024u}){
  auto malformed=m;if(field==0)malformed.value.resize(size);else if(field==1)malformed.ordered_ids.resize(size);else malformed.coordinate_units.resize(size);
  need(!s::GaussianSimulation::payload_bound(g,malformed,3,7),"malformed mean shape has no public payload bound");
  auto tiny=p;tiny.maximum_payload_bytes=1;armed=true;calls=0;fail_on=1;
  const auto rejected=s::GaussianSimulation::simulate(g,malformed,a,17,tiny);armed=false;
  need(rejected.status==n::Status::invalid_input&&rejected.rows.empty()&&rejected.work_units==0,"structural mean mismatch precedes tiny quota refusal");
  need(calls==0,"malformed mean shape refused without allocation");
 }
 auto bad=m;std::swap(bad.ordered_ids[0],bad.ordered_ids[1]);need(s::GaussianSimulation::simulate(g,bad,a,17,p).rows.empty(),"exact mean axes");
 bad=m;bad.coordinate_measure="wrong";need(s::GaussianSimulation::simulate(g,bad,a,17,p).rows.empty(),"coordinate measure");
 bad=m;bad.value[0]=std::numeric_limits<double>::denorm_min();need(s::GaussianSimulation::simulate(g,bad,a,17,p).rows.empty(),"subnormal mean");
 auto limit=p;limit.maximum_work_units=*s::GaussianSimulation::work_bound(2,3)-1;
 need(s::GaussianSimulation::simulate(g,m,a,17,limit).status==n::Status::work_limit,"combined work boundary");
 limit=p;limit.maximum_payload_bytes=*s::GaussianSimulation::payload_bound(g,m,3,7)-1;
 need(s::GaussianSimulation::simulate(g,m,a,17,limit).status==n::Status::work_limit,"payload byte boundary");
 limit.maximum_payload_bytes+=1;need(s::GaussianSimulation::simulate(g,m,a,17,limit).status==n::Status::ok,"exact payload boundary");
 limit=p;limit.outputs=8;need(s::GaussianSimulation::simulate(g,m,a,17,limit).rows.empty(),"unknown output mask");
 limit=p;limit.outputs=4;auto only=s::GaussianSimulation::simulate(g,m,a,17,limit);
 need(only.values.empty()&&only.absolute_error_estimates.empty()&&only.words==b.words,"request only words");
 limit=p;limit.maximum_scaled_arithmetic_error=1e-30;auto refused=s::GaussianSimulation::simulate(g,m,a,17,limit);
 need(refused.status==n::Status::ok&&refused.rows.size()==3,"row refusals retained");
 for(auto&r:refused.rows)need(r.status==n::Status::conditioning_budget_exceeded,"tight arithmetic request refuses without budget change");
 need(refused.words==b.words,"failed vectors retain exact words");
 const auto rounding=std::fegetround();std::fesetround(FE_UPWARD);
 need(s::GaussianSimulation::simulate(g,m,a,17,p).status==n::Status::outside_domain,"rounding contract");std::fesetround(rounding);
 armed=true;calls=0;fail_on=SIZE_MAX;auto measured=s::GaussianSimulation::simulate(g,m,a,17,p);armed=false;
 const auto sites=calls.load();need(sites>0&&measured.status==n::Status::ok,"observable allocations");
 for(size_t site=1;site<=sites;++site){armed=true;calls=0;fail_on=site;
  auto r=s::GaussianSimulation::simulate(g,m,a,17,p);armed=false;
  need(r.status==n::Status::work_limit || (r.status==n::Status::ok&&r.rows.size()==3),"every allocation refusal typed or retained row");
  need(g.status()==s::DensityStatus::finite,"borrowed covariance survives failure");
  if(r.status==n::Status::ok)need(r.rows[0].status!=n::Status::ok||r.rows[1].status!=n::Status::ok||r.rows[2].status!=n::Status::ok,"failed allocation never silent success");}
 for(size_t length:{16u,29u,43u,400u}){
  s::Metadata md;md.ordered_ids={std::string(length,'a'),std::string(length,'b')};md.measure=std::string(length,'m');
  md.table_identity=std::string(length,'t');md.uncertainty_identity=std::string(length,'u');
  md.ordering_provenance=std::string(length,'o');md.calibration_provenance=std::string(length,'c');
  md.dependence_provenance=std::string(length,'d');md.source_semantics="synthetic controls";
  auto long_source=s::prepare_gaussian(std::array<double,4>{4,1,1,1.25},s::MatrixKind::covariance,md,100,1e-10,n::Arithmetic::longdouble_cpu_v1);
  auto long_mean=m;long_mean.ordered_ids=md.ordered_ids;long_mean.coordinate_units={std::string(length,'u'),std::string(length,'v')};
  long_mean.identity=std::string(length,'i');long_mean.coordinate_measure=md.measure;long_mean.generating_law_identity=std::string(length,'g');
  const auto bound=*s::GaussianSimulation::payload_bound(long_source,long_mean,3,7);
  auto boundary=p;boundary.maximum_payload_bytes=bound-1;
  need(s::GaussianSimulation::simulate(long_source,long_mean,a,17,boundary).status==n::Status::work_limit,"long copied-string lower boundary refused");
  boundary.maximum_payload_bytes=bound;need(live_bytes==0,"allocation instrumentation baseline");peak_bytes=0;measure_allocations=true;
  {auto actual=s::GaussianSimulation::simulate(long_source,long_mean,a,17,boundary);measure_allocations=false;
   need(actual.status==n::Status::ok&&actual.rows.size()==3,"long copied-string exact envelope admitted");
   need(peak_bytes<=bound-*long_source.retained_payload_bound(),"observed requested allocation peak within copied-string envelope");
   std::cout<<"copy_length="<<length<<" observed_peak="<<peak_bytes<<" new_payload_envelope="<<bound-*long_source.retained_payload_bound()<<" arithmetic_length="<<actual.source_metadata.arithmetic_id.size()<<"\n";
   need(actual.source_metadata.arithmetic_id==long_source.metadata().arithmetic_id&&actual.source_metadata.measure==md.measure&&actual.source_metadata.table_identity==md.table_identity&&actual.source_metadata.ordering_provenance==md.ordering_provenance,"long arithmetic measure origin strings preserved");}
  need(live_bytes==0,"tracked outputs and scratch released");
 }
 // Recorder controls are consumers of the same owner used by both campaigns.
 std::ostringstream record;auto*saved=std::cout.rdbuf(record.rdbuf());
 bool stopped=false;try{auto wrong=m;wrong.ordered_ids[0]="wrong";(void)recovery_controls::draw(g,wrong,17,70,9,2);}catch(const std::runtime_error&){stopped=true;}
 auto invalid=recovery_controls::observation({1,0,0,1},{0,0},std::array<double,2>{std::numeric_limits<double>::quiet_NaN(),2},std::array<double,2>{0,0});
 const std::array<double,2>future{std::numeric_limits<double>::infinity(),3};s::GaussianSimulationBatch empty;
 recovery_controls::failed("recorder control",17,9,"training_observation",m.value,invalid,future,empty,empty,empty,0,s::DensityStatus::numerical_failure,recovery_controls::observation_status(invalid));
 std::cout.rdbuf(saved);const auto text=record.str();
 need(stopped&&text.find("batch_admission_refusal")!=std::string::npos&&text.find("\"planned_vectors\":2")!=std::string::npos&&text.find("\"attempted_vectors\":0")!=std::string::npos,"admission refusal recorded before stop with actual counts");
 need(text.find("failed_attempt")!=std::string::npos&&text.find("training_observation")!=std::string::npos&&text.find("\"future\":[{\"ieee754_hex\"")!=std::string::npos&&text.find("\"numerical_status\":"+std::to_string(int(n::Status::overflow)))!=std::string::npos,"original nonfinite observations and concrete statuses retained");
 need(text.find("7ff")!=std::string::npos&&text.find("\"attempted_vectors\":1")!=std::string::npos,"IEEE nonfinite identities and attempted denominator retained");
 recovery_controls::reserve_range(17,700,0,2,2);bool collision=false;
 try{recovery_controls::reserve_range(17,701,1,2,2);}catch(const std::runtime_error&){collision=true;}need(collision,"global campaign range overlap refused");
 bool wrap=false;try{recovery_controls::reserve_range(17,UINT64_MAX,0,1,2);}catch(const std::runtime_error&){wrap=true;}need(wrap,"global campaign counter overflow refused");
 std::cout<<"allocation_sites="<<sites<<"\n";
}
}
int main(){try{run();std::cout<<"PASS "<<checks<<" Gaussian simulation controls\n";}catch(const std::exception&e){armed=false;std::cerr<<"FAIL "<<checks<<": "<<e.what()<<'\n';return 1;}}
