// Bounded one-allocation-at-a-time synthetic control. Execution requires the
// separate root lease; this test never reads a released observation asset.
#include "box_heldout_controls.hpp"
#include "../src/retained_gaussian.hpp"
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <new>
namespace allocation_control {
struct alignas(std::max_align_t) Header {std::size_t bytes;};
std::size_t live=0,initial=0,peak=0,ordinal=0,fail_at=0;
bool armed=false;
constexpr std::size_t maximum_ordinals=512;
void arm(std::size_t at) {initial=live;peak=live;ordinal=0;fail_at=at;armed=true;}
void disarm() noexcept {armed=false;}
std::size_t new_peak() noexcept {return peak>initial?peak-initial:0;}
}
void*operator new(std::size_t n) {
  using namespace allocation_control;
  if(armed&&++ordinal==fail_at)throw std::bad_alloc();
  if(n>std::numeric_limits<std::size_t>::max()-sizeof(Header))throw std::bad_alloc();
  auto*h=static_cast<Header*>(std::malloc(sizeof(Header)+(n?n:1)));
  if(!h)throw std::bad_alloc();h->bytes=n;live+=n;if(armed&&live>peak)peak=live;return h+1;
}
void*operator new[](std::size_t n){return ::operator new(n);}
void operator delete(void*p)noexcept {
  if(!p)return;auto*h=static_cast<allocation_control::Header*>(p)-1;
  allocation_control::live-=h->bytes;std::free(h);
}
void operator delete[](void*p)noexcept{::operator delete(p);}
void operator delete(void*p,std::size_t)noexcept{::operator delete(p);}
void operator delete[](void*p,std::size_t)noexcept{::operator delete(p);}
using namespace irred::statistics;
namespace {
void check(bool value,const char*message){if(!value)throw std::runtime_error(message);}
struct Trace {std::size_t allocations=0,peak_new_requested_bytes=0;bool late_training=false,active_query=false,late_profile=false,late_heldout=false;};
Trace preparation(std::size_t failure) {
  box_heldout_controls::Input in;auto g=in.gaussian();const auto policy=in.policy(g);
  // All by-value argument copies occur before injection. Moves at entry do
  // not allocate, so this control concerns the actual in-body receipt promise.
  auto design=in.design;auto support=in.box;auto metadata=in.heldout;
  GaussianBoxHeldout owner;bool escaped=false;
  allocation_control::arm(failure);
  try {owner=GaussianBoxHeldout::prepare(std::move(g),in.x,in.offsets,in.rows.ordered_ids,
      std::move(design),std::move(support),in.heldout_index,std::move(metadata),policy);}
  catch(const std::bad_alloc&){escaped=true;}
  allocation_control::disarm();
  Trace trace{allocation_control::ordinal,allocation_control::new_peak()};
  check(!escaped,"in-body setup allocation must return a causal refusal");
  const auto&r=owner.preparation();check(r.peak_payload_bound.has_value(),"setup payload receipt retained");
  check(trace.peak_new_requested_bytes<=*r.peak_payload_bound,"setup new requested heap obeys owned payload upper");
  if(owner.status()==DensityStatus::finite)check(g.status()==DensityStatus::invalid_input,"only completed setup consumes source");
  else {
    check(owner.status()==DensityStatus::numerical_failure&&r.numerical_status==irred::numerics::Status::work_limit,"setup allocation cause retained");
    check(g.status()==DensityStatus::finite&&g.input_matrix_kind()==MatrixKind::covariance&&
        std::equal(g.covariance().begin(),g.covariance().end(),in.covariance.begin(),in.covariance.end()),"all failed setup allocation prefixes preserve original covariance");
  }
  check(r.work.whitenings_completed<=r.work.whitening_attempts,"causal setup whitening counters");
  trace.late_heldout=r.stage==BoxHeldoutStage::training_qr&&r.work.training_qr_attempts==1&&
      r.work.training_qr_completed==0&&r.work.whitening_attempts==3&&r.work.whitenings_completed==3;
  return trace;
}
Trace profile(std::size_t failure) {
  box_heldout_controls::Input in;auto g=in.gaussian();auto design=in.design;
  const std::vector<double>x={1,0,0,1,.5,.25};DesignProfile out;bool escaped=false;
  allocation_control::arm(failure);
  try {out=DesignProfile::prepare(std::move(g),x,in.rows.ordered_ids,std::move(design));}
  catch(const std::bad_alloc&){escaped=true;}
  allocation_control::disarm();Trace trace{allocation_control::ordinal,allocation_control::new_peak()};
  // A failure constructing the function-local default owner can precede its
  // receipt. It performs no whitening and still preserves the source.
  if(escaped){check(g.status()==DensityStatus::finite,"pre-receipt profile allocation preserves source");return trace;}
  const auto counts=irred::detail::RetainedQrAccess::preparation_whitenings(out);
  check(counts.second<=counts.first,"causal private profile counters");
  if(out.status()!=DensityStatus::finite) {
    check(g.status()==DensityStatus::finite,"late profile failure preserves source");
    check(out.status()==DensityStatus::numerical_failure&&out.numerical_status()==irred::numerics::Status::work_limit,"profile allocation cause returned");
    trace.late_profile=counts.first==2&&counts.second==2;
  }
  return trace;
}
Trace evaluation(std::size_t failure) {
  box_heldout_controls::Input in;auto g=in.gaussian();const auto policy=in.policy(g,2,2,3);
  auto owner=in.prepare(std::move(g),policy);check(owner.status()==DensityStatus::finite,"fault-control owner prepared");
  const std::vector<double>training={1.25,1.5,1.5,2},values={3,3.25};
  const std::vector<BoxHeldoutRequest>requests={{1,0},{0,1},{1,0}};
  BoxHeldoutBatch out;bool escaped=false;
  allocation_control::arm(failure);
  try {out=owner.evaluate(training,owner.training_row_ids(),values,requests,policy);}
  catch(const std::bad_alloc&){escaped=true;}
  allocation_control::disarm();Trace trace{allocation_control::ordinal,allocation_control::new_peak()};
  check(!escaped,"evaluation allocation must return a causal refusal");
  check(owner.status()==DensityStatus::finite,"evaluation failure retains prepared owner");
  const auto bound=owner.evaluation_payload_bound(2,3);check(bound.has_value(),"evaluation payload bound available");
  check(trace.peak_new_requested_bytes<=*bound,"evaluation new requested heap obeys owned payload upper");
  if(out.status()!=DensityStatus::finite)check(out.status()==DensityStatus::numerical_failure&&out.numerical_status==irred::numerics::Status::work_limit,"outer allocation cause retained");
  if(out.output_layout_available) {
    check(out.requests.size()==3&&out.densities.size()==3&&out.training_normalizations.size()==2&&
        out.training_refusal_witnesses.size()==2&&out.training_offset_subtraction_rounding_estimates.size()==2&&
        out.training_cross_projection_error_estimates.size()==2&&out.training_offset_diagnostics_available.size()==2&&
        out.training_cross_projection_diagnostics_available.size()==2,"earned output layout has matching lengths");
  }
  check(out.work.factor_attempts==0&&out.work.whitenings_completed<=out.work.whitening_attempts,"evaluation retains factor and causal counters");
  for(std::size_t k=0;k<out.training_normalizations.size();++k) {
    const auto&normal=out.training_normalizations[k];
    if(normal.normalization_enclosures_available&&normal.stage==BoxStage::complete&&normal.status==DensityStatus::numerical_failure) {
      check(normal.numerical_status==irred::numerics::Status::work_limit,"late training refusal causal");trace.late_training=true;
      for(std::size_t q=0;q<out.densities.size();++q)if(out.requests[q].training_vector_index==k)
        check(!out.densities[q].density_available,"failed training cache withholds dependent requests");
    }
  }
  for(const auto&result:out.densities)if(result.work.joint_attempts&&result.status!=DensityStatus::finite) {
    check(result.status==DensityStatus::numerical_failure&&result.numerical_status==irred::numerics::Status::work_limit&&
        !result.density_available,"active query allocation retains cause and work delta");trace.active_query=true;
  }
  return trace;
}
}
int main(){try {
  const auto none=std::numeric_limits<std::size_t>::max(),live_before=allocation_control::live;
  const auto prep=preparation(none),prof=profile(none),eval=evaluation(none);
  check(allocation_control::live==live_before,"baseline releases every requested C++ heap byte");
  check(prep.allocations<=allocation_control::maximum_ordinals&&prof.allocations<=allocation_control::maximum_ordinals&&
      eval.allocations<=allocation_control::maximum_ordinals,"frozen allocation-ordinal cap; never enlarge automatically");
  bool late_profile=false,late_training=false,active_query=false,late_heldout=false;
  for(std::size_t k=1;k<=prep.allocations;++k){late_heldout=preparation(k).late_heldout||late_heldout;check(allocation_control::live==live_before,"failed setup has no requested C++ heap leak");}
  for(std::size_t k=1;k<=prof.allocations;++k){late_profile=profile(k).late_profile||late_profile;check(allocation_control::live==live_before,"failed profile has no requested C++ heap leak");}
  for(std::size_t k=1;k<=eval.allocations;++k) {const auto trace=evaluation(k);
    late_training=trace.late_training||late_training;active_query=trace.active_query||active_query;
    check(allocation_control::live==live_before,"failed evaluation has no requested C++ heap leak");}
  check(late_heldout&&late_profile&&late_training&&active_query,"late heldout/private profile/completed-normalizer cache/query allocation paths all challenged");
  // Injection is disabled before stdout or a caller's serializer is used.
  std::cout<<"allocation ordinals setup="<<prep.allocations<<" profile="<<prof.allocations<<" evaluation="<<eval.allocations
      <<" new_requested_heap_peaks="<<prep.peak_new_requested_bytes<<','<<prof.peak_new_requested_bytes<<','<<eval.peak_new_requested_bytes<<'\n';return 0;
}catch(const std::exception&e){allocation_control::disarm();std::cerr<<e.what()<<'\n';return 1;}}
