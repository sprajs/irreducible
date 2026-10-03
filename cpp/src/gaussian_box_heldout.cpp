#include "irred/gaussian_box_heldout.hpp"
#include "box_enclosure.hpp"
#include "retained_gaussian.hpp"
#include "payload_accounting.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <new>
#include <utility>

namespace irred::statistics {
namespace {
using LD=long double;
using Refusal=detail::BoxRefusal;
using Accounting=detail::PayloadAccounting;
// Logical copy/event bytes depend on contents, never spare allocation. The
// same field traversal below still uses capacity-based Accounting for payload.
class LogicalMetadataCopies {
  Accounting total_{0};
public:
  template<class T> void vector(const std::vector<T>&v)noexcept {
    total_.add(v.size(),sizeof(T)); // logical object bytes
    total_.add(v.size(),1);        // element-copy event
    total_.add(1,1);               // vector helper/header event
  }
  void string(const std::string&s)noexcept {
    total_.add(s.size(),1);        // character-copy events
    total_.add(1,8);               // copy/move/helper/terminator prefixes
    total_.add(2,sizeof(std::string)); // inline move buffers; never capacity
  }
  void strings(const std::vector<std::string>&v)noexcept {
    total_.add(1,1);               // string-vector helper entry
    vector(v);for(const auto&s:v)string(s);
  }
  std::optional<std::size_t>result()const noexcept {
    return total_.result();
  }
};
bool same(std::span<const std::string>a,std::span<const std::string>b) {
  return a.size()==b.size()&&std::equal(a.begin(),a.end(),b.begin());
}
bool unique(std::span<const std::string>a) {
  for(std::size_t i=0;i<a.size();++i)
    if(a[i].empty()||std::find(a.begin(),a.begin()+i,a[i])!=a.begin()+i) return false;
  return true;
}
template<class Counter>void charge(Counter &b,const Metadata &m) {
  b.strings(m.ordered_ids);
  for(const auto*s:{&m.arithmetic_id,&m.measure,&m.table_identity,&m.uncertainty_identity,
      &m.ordering_provenance,&m.calibration_provenance,&m.dependence_provenance,
      &m.source_semantics,&m.input_matrix_convention,&m.treatment}) b.string(*s);
}
template<class Counter>void charge(Counter &b,const DesignMetadata &m) {
  b.strings(m.ordered_parameter_ids);b.strings(m.parameter_units);b.strings(m.shared_nuisance_ids);
  b.string(m.residual_unit);b.string(m.design_identity);b.string(m.dependence_identity);
}
template<class Counter>void charge(Counter &b,const BoxSupport&s) {
  b.strings(s.ordered_parameter_ids);b.vector(s.lower);b.vector(s.upper);
  b.string(s.parameter_measure);b.string(s.prior_identity);b.string(s.fixed_coordinate_provenance);
}
template<class Counter>void charge(Counter &b,const BoxHeldoutMetadata&m) {
  for(const auto*s:{&m.source_contract_id,&m.candidate_identity,&m.conditioning_identity,
      &m.original_row_lineage,&m.event_lineage,&m.calibration_dependence_identity,
      &m.heldout_unit,&m.heldout_covariance_unit,&m.heldout_measure,&m.fixed_coordinate_provenance}) b.string(*s);
  b.strings(m.ordered_original_parameter_ids);b.vector(m.active_original_parameter_indices);
}
bool valid(const BoxHeldoutPolicy&p) {
  return detail::box_valid_policy(p.box)&&p.maximum_training_vectors&&p.maximum_candidate_values&&
      p.maximum_requests&&p.maximum_preparation_work_units&&p.maximum_evaluation_work_units&&
      std::isfinite(p.maximum_conditional_log_density_width)&&p.maximum_conditional_log_density_width>0;
}
bool tightens(const BoxHeldoutPolicy&p,const BoxHeldoutPolicy&f) {
  return valid(p)&&p.box.design.maximum_elements<=f.box.design.maximum_elements&&
      p.box.design.maximum_payload_bytes<=f.box.design.maximum_payload_bytes&&
      p.box.design.maximum_forward_sensitivity<=f.box.design.maximum_forward_sensitivity&&
      p.box.maximum_log_probability_width<=f.box.maximum_log_probability_width&&
      p.box.maximum_quantile_width<=f.box.maximum_quantile_width&&
      p.box.cdf_absolute_radius<=f.box.cdf_absolute_radius&&p.box.maximum_cdf_nodes<=f.box.maximum_cdf_nodes&&
      p.box.maximum_bisections<=f.box.maximum_bisections&&p.box.maximum_cdf_evaluations<=f.box.maximum_cdf_evaluations&&
      p.maximum_training_vectors<=f.maximum_training_vectors&&p.maximum_candidate_values<=f.maximum_candidate_values&&
      p.maximum_requests<=f.maximum_requests&&p.maximum_preparation_work_units<=f.maximum_preparation_work_units&&
      p.maximum_evaluation_work_units<=f.maximum_evaluation_work_units&&
      p.maximum_conditional_log_density_width<=f.maximum_conditional_log_density_width;
}
// Upper logical-event envelopes. These include helper entries and failure
// prefixes. Arithmetic overflows refuse before the bounded operation starts.
std::optional<std::size_t> whiten_work(std::size_t n) {
  if(!n)return {};Accounting b(8);b.add(n,n);b.add(n,4);return b.result();
}
std::optional<std::size_t> product(std::size_t a,std::size_t b) {
  std::size_t result=0;if(!detail::checked_payload_add(result,a,b))return {};return result;
}
std::optional<std::size_t> triangle(std::size_t n,bool diagonal=false) {
  const auto nn=product(n,diagonal?n+1:n-1);return nn?std::optional<std::size_t>(*nn/2):std::nullopt;
}
std::optional<std::size_t> factor_work(std::size_t n) {
  if(!n)return {};auto t=triangle(n),u=triangle(n,true),n2=product(n,n);
  auto nm=product(n,n-1),k=nm?product(*nm,n+1):std::nullopt;
  if(!t||!u||!n2||!k)return {};
  // Exact core loop sums include Cholesky dot K, every inverse-column
  // forward/back solve, input/symmetry scans and Gaussian covariance copy.
  Accounting b(32);b.add(*t,n);b.add(*t,n);b.add(*k/6,1);b.add(*n2,8);
  b.add(*t,2);b.add(*u,1);b.add(n,13);return b.result();
}
std::optional<std::size_t> conditions_work(std::size_t p) {
  if(p<2)return {};auto t=triangle(p),u=triangle(p,true),pp=product(p,p);if(!t||!u||!pp)return {};
  Accounting b(2);b.add(*t,p);b.add(*pp,3);b.add(*u,2);b.add(p,7);return b.result();
}
std::optional<std::size_t> qr_work(std::size_t n,std::size_t p) {
  if(p<2||n<p)return {};auto w=whiten_work(n),np=product(n,p),t=triangle(p),u=triangle(p,true),conditions=conditions_work(p);
  auto pp1=product(p,p+1),pprev=product(p-1,p);
  auto cubic_n=pp1?product(*pp1,2*p+1):std::nullopt,cubic_h=pprev?product(*pprev,2*p-1):std::nullopt;
  if(!w||!np||!t||!u||!conditions||!cubic_n||!cubic_h)return {};
  Accounting b(32);b.add(p,*w);b.add(*np,9);b.add(n,1);b.add(p,7);b.add(*u,1);
  b.add(n-p,*u);b.add(n-p,*u);b.add(*cubic_n/6,2);b.add(n,p-1);b.add(n,p-1);
  b.add(*np-*u,1);b.add(*t,1);b.add(n-p,*t);b.add(n-p,*t);b.add(*cubic_h/6,2);
  b.add(p,2);b.add(1,*conditions);return b.result();
}
std::optional<std::size_t> variance_work(std::size_t p) {
  if(p<2)return {};const auto t=triangle(p),u=triangle(p,true);if(!t||!u)return {};
  Accounting b(2);b.add(*t,1);b.add(*u,1);b.add(p,4);return b.result();
}
std::optional<std::size_t> completion_work(std::size_t,std::size_t p) {
  const auto v=variance_work(p);if(!v)return {};
  // Shared completion pV+6p+8 and normalizer22593p+20977, plus one
  // normalization wrapper entry. All bounded private elementary helpers count.
  Accounting b(20986);b.add(p,*v);b.add(p,22599);return b.result();
}
void precharge(BoxHeldoutWork&w,std::optional<std::size_t>cost,std::size_t cap) {
  if(!cost||w.charged_work_units>cap||*cost>cap-w.charged_work_units)
    throw Refusal{numerics::Status::work_limit};
  w.charged_work_units+=*cost;
}
LD checked(LD v) {
  if(!std::isfinite(v)||std::fpclassify(v)==FP_SUBNORMAL) throw Refusal{numerics::Status::outside_domain,v};
  return v;
}
LD rounding(LD x,double reported) {
  return checked(std::abs(x-static_cast<LD>(reported))/(1+std::abs(x)));
}
void withhold(GaussianBoxResult&r,numerics::Status s) {
  if(!r.gaussian_completion_available) {
    const auto step=r.completion_step;const auto index=r.completion_parameter_index;
    r=GaussianBoxResult{};r.completion_step=step;r.completion_parameter_index=index;
  }
  if(!r.endpoint_margins_available) {
    r.standardized_lower_margins.clear();r.standardized_upper_margins.clear();
    r.excluded_mass_upper=0;r.log_prior_volume={};
  }
  if(!r.rectangle_enclosure_available){r.box_probability={};r.log_box_probability={};}
  if(!r.normalization_enclosures_available) {
    r.log_relative_box_integral={};r.log_prior_normalized_relative_evidence={};r.log_observation_normalized_evidence={};
  }
  r.status=DensityStatus::numerical_failure;r.numerical_status=s;
}
void normalize(GaussianBoxResult&r,const detail::RetainedQrView&v,DesignResult&&fit,
    LD q,double logc,const BoxSupport&support,const BoxHeldoutPolicy&policy,BoxHeldoutWork&w) {
  ++w.variance_screen_attempts;detail::qr_complete(r,v,std::move(fit),q,logc,policy.box.design.maximum_forward_sensitivity);
  ++w.variance_screens_completed;++w.tail_bound_attempts;++w.normalization_attempts;
  detail::box_normalize(r,support,policy.box,v.n);
  ++w.tail_bounds_completed;++w.normalizations_completed;r.stage=BoxStage::complete;
  r.status=DensityStatus::finite;r.numerical_status=numerics::Status::ok;
}
} // namespace

struct GaussianBoxHeldout::Impl {
  explicit Impl(const Metadata&source):source_metadata(source){}
  DesignProfile training;
  Metadata source_metadata;
  BoxSupport support;
  BoxHeldoutMetadata metadata;
  BoxHeldoutPolicy frozen;
  std::vector<double> offsets,cross,heldout_design;
  std::vector<LD> h,appended,updated_r,cosine,sine;
  std::size_t heldout=0;
  double ctt=0;
  LD schur=0,sqrt_schur=0;
  double update_condition=0,update_transpose=0,update_sensitivity=0;
  BoxInterval conditional_constant;
  detail::RetainedQrView joint_view()const {
    const auto t=detail::RetainedQrAccess::view(training);
    return {t.n+1,t.p,updated_r,t.scales,t.pivot,update_condition,update_transpose};
  }
};
GaussianBoxHeldout::GaussianBoxHeldout()noexcept=default;
GaussianBoxHeldout::~GaussianBoxHeldout()=default;
GaussianBoxHeldout::GaussianBoxHeldout(GaussianBoxHeldout&&o)noexcept
    :impl_(std::move(o.impl_)),preparation_(std::exchange(o.preparation_,BoxHeldoutPreparation{})){}
GaussianBoxHeldout&GaussianBoxHeldout::operator=(GaussianBoxHeldout&&o)noexcept {
  if(this!=&o){impl_=std::move(o.impl_);preparation_=std::exchange(o.preparation_,BoxHeldoutPreparation{});}return *this;
}
DensityStatus GaussianBoxHeldout::status()const noexcept{return preparation_.status;}
const BoxHeldoutPreparation&GaussianBoxHeldout::preparation()const noexcept{return preparation_;}
std::span<const std::string>GaussianBoxHeldout::original_row_ids()const noexcept {
  return impl_ ? std::span<const std::string>(impl_->source_metadata.ordered_ids) : std::span<const std::string>{};
}
std::span<const std::string>GaussianBoxHeldout::training_row_ids()const noexcept {
  return impl_ ? std::span<const std::string>(impl_->training.metadata().ordered_ids) : std::span<const std::string>{};
}
const BoxHeldoutMetadata&GaussianBoxHeldout::metadata()const noexcept {
  static const BoxHeldoutMetadata absent;return impl_?impl_->metadata:absent;
}
std::optional<std::size_t>GaussianBoxHeldout::preparation_payload_bound(const Gaussian&g,
    const DesignMetadata&d,const BoxSupport&s,const BoxHeldoutMetadata&m)noexcept {
  const auto n=g.metadata().ordered_ids.size(),p=d.ordered_parameter_ids.size();
  if(n<3||p<2||n>SIZE_MAX/n||n>SIZE_MAX/p||p>SIZE_MAX/p)return {};
  const auto mrows=n-1;const auto mm=product(mrows,mrows),mp=product(mrows,p),pp=product(p,p);
  if(!mm||!mp||!pp)return {};
  // Named simultaneous objects, with no heap duplication across moves/returns.
  Accounting common(2*sizeof(GaussianBoxHeldout)+sizeof(Impl)+sizeof(DesignMetadata)+
      sizeof(BoxSupport)+sizeof(BoxHeldoutMetadata)+sizeof(Gaussian)+sizeof(Metadata)+
      3*sizeof(std::vector<double>)+sizeof(std::vector<LD>)+sizeof(numerics::WhiteningResult)+sizeof(std::unique_ptr<Impl>));
  common.embedded(g.retained_payload_bound(),0);
  charge(common,g.metadata());charge(common,g.metadata());
  charge(common,d);charge(common,s);charge(common,m);
  common.add(n,sizeof(double));common.add(mrows,sizeof(double));common.add(p,sizeof(double));
  // Default and emptied named metadata headers still own their inline string
  // capacities. Capacity+1 is counted again consistently with existing bounds.
  const auto empty_string=std::string{}.capacity()+1;
  common.add(4,11*empty_string+sizeof("F02/binary64-legacy/v1")+sizeof("normalized Gaussian"));
  common.add(3,3*empty_string);common.add(1,10*empty_string);
  const auto base=common.result();if(!base)return {};
  auto peak=[&](std::size_t header,auto numeric){Accounting phase(*base);phase.add(1,header);numeric(phase);return phase.result();};
  auto validation=peak(sizeof(Gaussian)+sizeof(Metadata),[&](Accounting&b){
    b.add(*mm,sizeof(double));b.add(*mp,sizeof(double));b.add(mrows,128);
    for(const auto&id:g.metadata().ordered_ids){b.add(id.capacity(),1);b.add(1,1);}});
  auto factor=peak(sizeof(Gaussian)+sizeof(Metadata)+2*sizeof(numerics::Factorization),[&](Accounting&b){
    b.add(*mp,sizeof(double));b.add(*mm,2*sizeof(double)+sizeof(LD));
    Accounting scratch(0);scratch.add(mrows,2*sizeof(LD)+sizeof(double));
    Accounting raw(0);raw.add(*mm,sizeof(double));
    if(!scratch.result()||!raw.result()){b.embedded({},0);return;}b.add(std::max(*scratch.result(),*raw.result()),1);});
  auto cross=peak(2*sizeof(numerics::WhiteningResult),[&](Accounting&b){
    b.add(*mm,2*sizeof(double)+sizeof(LD));b.add(*mp,sizeof(double));b.add(mrows,sizeof(LD));});
  auto training_normalization=peak(2*sizeof(DesignProfile)+sizeof(DesignMetadata)+sizeof(std::vector<LD>)+
      2*sizeof(numerics::WhiteningResult)+sizeof(Gaussian),[&](Accounting&b){
    b.add(*mm,2*sizeof(double)+sizeof(LD));b.add(*mp,2*sizeof(double)+sizeof(LD));
    b.add(mrows,3*sizeof(LD));b.add(p,sizeof(LD));});
  auto training=peak(2*sizeof(DesignProfile)+sizeof(DesignMetadata)+3*sizeof(std::vector<LD>)+
      2*sizeof(numerics::WhiteningResult)+sizeof(Gaussian),[&](Accounting&b){
    b.add(*mm,2*sizeof(double)+sizeof(LD));b.add(*mp,2*sizeof(double)+2*sizeof(LD));
    b.add(mrows,2*sizeof(LD));b.add(p,4*sizeof(LD)+sizeof(std::size_t));});
  auto update=peak(2*sizeof(std::vector<LD>)+2*sizeof(detail::QrConditions),[&](Accounting&b){
    b.add(*mm,sizeof(double)+sizeof(LD));b.add(*mp,sizeof(double)+2*sizeof(LD));
    b.add(mrows,sizeof(LD));b.add(p,2*sizeof(LD)+sizeof(std::size_t));
    b.add(*pp,sizeof(LD));b.add(p,6*sizeof(LD));});
  auto variance=peak(2*sizeof(std::vector<LD>)+3*sizeof(EstimatorVarianceResult),[&](Accounting&b){
    b.add(*mm,sizeof(double)+sizeof(LD));b.add(*mp,sizeof(double)+2*sizeof(LD));
    b.add(mrows,sizeof(LD));b.add(p,2*sizeof(LD)+sizeof(std::size_t));
    b.add(*pp,sizeof(LD));b.add(p,6*sizeof(LD)+sizeof(double));b.add(4,sizeof(LD));b.add(6,empty_string);});
  auto consumption=peak(sizeof(Gaussian),[&](Accounting&b){
    b.add(*mm,sizeof(double)+sizeof(LD));b.add(*mp,sizeof(double)+2*sizeof(LD));
    b.add(mrows,sizeof(LD));b.add(p,2*sizeof(LD)+sizeof(std::size_t));
    b.add(*pp,sizeof(LD));b.add(p,4*sizeof(LD)+sizeof(double));
    b.add(1,11*empty_string);});
  if(!validation||!factor||!cross||!training_normalization||!training||!update||!variance||!consumption)return {};
  return std::max({*validation,*factor,*cross,*training_normalization,*training,*update,*variance,*consumption});
}
std::optional<std::size_t>GaussianBoxHeldout::preparation_work_bound(const Gaussian&g,
    const DesignMetadata&d,const BoxSupport&s,const BoxHeldoutMetadata&m)noexcept {
  const auto n=g.metadata().ordered_ids.size(),p=d.ordered_parameter_ids.size();
  if(n<3||p<2||p>n-1)return {};
  auto payload=preparation_payload_bound(g,d,s,m),factor=factor_work(n-1),qr=qr_work(n-1,p),w=whiten_work(n-1);
  if(!payload||!factor||!qr||!w)return {};
  Accounting names(0);for(const auto*set:{&g.metadata().ordered_ids,&d.ordered_parameter_ids,&d.shared_nuisance_ids,&m.ordered_original_parameter_ids})
    for(const auto&id:*set){names.add(id.size(),1);names.add(1,1);}
  if(!names.result())return {};
  LogicalMetadataCopies copies;charge(copies,g.metadata());charge(copies,d);charge(copies,s);charge(copies,m);
  if(!copies.result())return {};
  Accounting b(256);b.add(*copies.result(),4);b.add(*names.result(),n);b.add(*names.result(),p);b.add(*names.result(),p);b.add(*names.result(),8);
  b.add(1,*factor);b.add(1,*qr);b.add(1,*w);b.add(n,n);b.add(n,n);b.add(n,p);b.add(n,p);b.add(n,p);b.add(n,p);
  b.add(n,p);b.add(n,p);b.add(n,p);b.add(n,p);b.add(n,64);b.add(p,128);
  const auto pp=product(p,p),conditions=conditions_work(p),variance=variance_work(p);if(!pp||!conditions||!variance)return {};
  b.add(*pp,64);b.add(1,*conditions);b.add(p,*variance);
  b.add(1,20764);b.add(1,617);return b.result();
}
std::optional<std::size_t>GaussianBoxHeldout::retained_payload_bound()const noexcept {
  Accounting b(sizeof(*this));if(!impl_)return b.result();b.add(1,sizeof(Impl));
  b.embedded(impl_->training.retained_payload_bound(),sizeof(DesignProfile));charge(b,impl_->source_metadata);
  charge(b,impl_->support);charge(b,impl_->metadata);b.vector(impl_->offsets);b.vector(impl_->cross);
  b.vector(impl_->heldout_design);b.vector(impl_->h);b.vector(impl_->appended);b.vector(impl_->updated_r);
  b.vector(impl_->cosine);b.vector(impl_->sine);return b.result();
}
std::optional<std::size_t>GaussianBoxHeldout::evaluation_payload_bound(std::size_t g,std::size_t r)const noexcept {
  if(!impl_)return {};const auto v=detail::RetainedQrAccess::view(impl_->training);
  const auto gp=product(g,v.p),rp=product(r,v.p);if(!gp||!rp)return {};
  Accounting persistent(2*sizeof(BoxHeldoutBatch)+sizeof(std::vector<std::vector<LD>>));
  persistent.add(g,sizeof(GaussianBoxResult)+2*sizeof(double)+2*sizeof(std::uint8_t)+sizeof(std::optional<LD>)+sizeof(std::vector<LD>));
  persistent.add(*gp,2*sizeof(LD)+2*sizeof(double)+2*sizeof(BoxInterval));persistent.add(g,2*sizeof(LD));
  persistent.add(r,sizeof(BoxHeldoutDensity)+sizeof(BoxHeldoutRequest));persistent.add(*rp,2*sizeof(double)+2*sizeof(BoxInterval));
  const auto fixed=persistent.result();if(!fixed)return {};
  auto phase=[&](std::size_t nwidth,std::size_t pwidth,std::size_t headers){Accounting b(*fixed);
    b.add(v.n,nwidth);b.add(1,nwidth);b.add(v.p,pwidth);b.add(1,headers);b.add(4,sizeof(LD));return b.result();};
  const auto fit_headers=3*sizeof(DesignResult)+3*sizeof(numerics::WhiteningResult)+4*sizeof(std::vector<LD>);
  const auto empty_string=std::string{}.capacity()+1;
  // qr_complete owns weights; qr_variance's three failed/returned headers each
  // contain two empty functional-metadata strings. Training also retains wr.
  const auto training_completion_headers=sizeof(DesignResult)+sizeof(numerics::WhiteningResult)+
      sizeof(std::vector<double>)+3*sizeof(EstimatorVarianceResult)+4*sizeof(std::vector<LD>)+6*empty_string;
  const auto joint_completion_headers=sizeof(DesignResult)+sizeof(std::vector<double>)+
      3*sizeof(EstimatorVarianceResult)+5*sizeof(std::vector<LD>)+6*empty_string;
  // Joint whitening reserve(n) temporarily retains the old m entries beside
  // the newly allocated n entries. Those two buffers coexist at this peak.
  auto training_fit=phase(3*sizeof(LD)+sizeof(double),2*sizeof(LD),fit_headers);
  auto training_completion=phase(2*sizeof(LD)+sizeof(double),3*sizeof(LD)+sizeof(double),training_completion_headers);
  auto joint_fit=phase(3*sizeof(LD)+sizeof(double),3*sizeof(LD),fit_headers);
  auto joint_completion=phase(sizeof(LD)+sizeof(double),4*sizeof(LD)+sizeof(double),joint_completion_headers);
  // Refusal cleanup can construct a fresh result header while the persistent
  // partial output still owns its buffers. No sizeof dominance is assumed.
  auto cleanup=phase(0,0,sizeof(GaussianBoxResult));
  if(!training_fit||!training_completion||!joint_fit||!joint_completion||!cleanup)return {};
  return std::max({*training_fit,*training_completion,*joint_fit,*joint_completion,*cleanup});
}
std::optional<std::size_t>GaussianBoxHeldout::evaluation_work_bound(std::size_t n,std::size_t p,std::size_t g,std::size_t candidates,std::size_t r,std::size_t identity_bytes)noexcept {
  if(n<3||p<2||p>n-1||p>SIZE_MAX/p||n>SIZE_MAX/n)return {};
  auto w=whiten_work(n-1),c=completion_work(n,p);if(!w||!c)return {};
  const auto np=product(n-1,p),nfullp=product(n,p),t=triangle(p);if(!np||!nfullp||!t)return {};
  // Training: two whitenings, shared fit/backsolve/diagnostics, cached QR
  // transform+tail/projection, completion and cache element writes.
  Accounting training(64);training.add(2,*w);training.add(*np,7);training.add(n-1,16);
  training.add(p,13);training.add(*t,1);training.add(1,*c);
  // Joint: one actual-adjusted whitening plus its appended row diagnostics,
  // original residual/cached RHS rotations, shared fit and ratio composition.
  Accounting joint(128);joint.add(1,*w);joint.add(*nfullp,4);joint.add(n,16);
  joint.add(p,16);joint.add(*t,1);joint.add(1,*c);
  if(!training.result()||!joint.result())return {};
  Accounting total(64);total.add(g,*training.result());total.add(r,*joint.result());total.add(g,n-1);
  total.add(r,4);total.add(candidates,2);total.add(n-1,8);total.add(identity_bytes,1);total.add(1,16);return total.result();
}
std::optional<std::size_t>GaussianBoxHeldout::evaluation_work_bound(std::size_t g,std::size_t candidates,std::size_t r)const noexcept {
  if(!impl_)return {};const auto v=detail::RetainedQrAccess::view(impl_->training);
  Accounting identity(0);for(const auto&id:impl_->training.metadata().ordered_ids)identity.add(id.size(),1);
  if(!identity.result())return {};return evaluation_work_bound(v.n+1,v.p,g,candidates,r,*identity.result());
}
GaussianBoxHeldout GaussianBoxHeldout::prepare(Gaussian&&full,std::span<const double>x,
    std::span<const double>offsets,std::span<const std::string>ids,DesignMetadata design,
    BoxSupport support,std::size_t t,BoxHeldoutMetadata metadata,BoxHeldoutPolicy policy) {
  GaussianBoxHeldout out;auto&r=out.preparation_;r.stage=BoxHeldoutStage::source_validation;
  auto fail=[&](DensityStatus ds,numerics::Status ns){r.status=ds;r.numerical_status=ns;return std::move(out);};
  if(!detail::box_supported())return fail(DensityStatus::unsupported_domain,numerics::Status::outside_domain);
  if(full.status()!=DensityStatus::finite)return fail(full.status(),full.numerical_status());
  const auto n=ids.size(),p=design.ordered_parameter_ids.size(),m=n?n-1:0;
  if(!valid(policy)||n<3||p<2||m<p||t>=n||n>SIZE_MAX/p||n>SIZE_MAX/n||
      x.size()!=n*p||offsets.size()!=n||design.parameter_units.size()!=p||
      design.shared_nuisance_ids.size()>p||metadata.ordered_original_parameter_ids.size()!=p+1||
      metadata.active_original_parameter_indices.size()!=p||full.metadata().ordered_ids.size()!=n)return fail(DensityStatus::invalid_input,numerics::Status::invalid_input);
  const auto peak=preparation_payload_bound(full,design,support,metadata);r.peak_payload_bound=peak;
  if(!peak||*peak>policy.box.design.maximum_payload_bytes||n>policy.box.design.maximum_elements/n||
      n>policy.box.design.maximum_elements/p||p>policy.box.design.maximum_elements/p)
    return fail(DensityStatus::numerical_failure,numerics::Status::work_limit);
  try {
    precharge(r.work,preparation_work_bound(full,design,support,metadata),policy.maximum_preparation_work_units);
    if(!same(ids,full.metadata().ordered_ids)||full.input_matrix_kind()!=MatrixKind::covariance||
        full.metadata().matrix_validation_scope!=MatrixValidationScope::full_declared_matrix||
        !full.priors().empty()||!full.mean_shift().empty()||!full.selection_history().empty()||
        full.covariance().size()!=n*n||support.ordered_parameter_ids!=design.ordered_parameter_ids||
        metadata.ordered_original_parameter_ids.size()!=p+1||metadata.active_original_parameter_indices.size()!=p||
        metadata.fixed_original_parameter_index>=p+1||!unique(metadata.ordered_original_parameter_ids)||
        support.fixed_coordinate_provenance.empty()||metadata.fixed_coordinate_provenance!=support.fixed_coordinate_provenance||
        metadata.heldout_unit!=design.residual_unit||metadata.fixed_original_parameter_value!=0)
      return fail(DensityStatus::incompatible_metadata,numerics::Status::invalid_input);
    for(const auto*s:{&metadata.source_contract_id,&metadata.candidate_identity,&metadata.conditioning_identity,
        &metadata.original_row_lineage,&metadata.event_lineage,&metadata.calibration_dependence_identity,
        &metadata.heldout_covariance_unit,&metadata.heldout_measure})
      if(s->empty())return fail(DensityStatus::incompatible_metadata,numerics::Status::invalid_input);
    std::size_t active=0;
    for(std::size_t j=0;j<p+1;++j)if(j!=metadata.fixed_original_parameter_index) {
      if(metadata.active_original_parameter_indices[active]!=j||design.ordered_parameter_ids[active]!=
          metadata.ordered_original_parameter_ids[j])return fail(DensityStatus::incompatible_metadata,numerics::Status::invalid_input);
      ++active;
    }
    if(!std::all_of(x.begin(),x.end(),[](double v){return std::isfinite(v);})||
        !std::all_of(offsets.begin(),offsets.end(),[](double v){return std::isfinite(v);})||
        support.lower.size()!=p||support.upper.size()!=p||support.parameter_measure.empty()||support.prior_identity.empty())
      return fail(DensityStatus::invalid_input,numerics::Status::invalid_input);
    for(std::size_t j=0;j<p;++j)
      if(!std::isfinite(support.lower[j])||!std::isfinite(support.upper[j])||!(support.lower[j]<support.upper[j])||
          (support.lower[j]!=0&&std::fpclassify(support.lower[j])!=FP_NORMAL)||
          (support.upper[j]!=0&&std::fpclassify(support.upper[j])!=FP_NORMAL))
        return fail(DensityStatus::invalid_input,numerics::Status::invalid_input);
    auto candidate=std::make_unique<Impl>(full.metadata());candidate->frozen=policy;candidate->heldout=t;
    candidate->support=std::move(support);candidate->metadata=std::move(metadata);
    candidate->offsets.assign(offsets.begin(),offsets.end());candidate->heldout_design.assign(x.begin()+t*p,x.begin()+(t+1)*p);
    Metadata train_metadata=full.metadata();train_metadata.ordered_ids.erase(train_metadata.ordered_ids.begin()+t);
    std::vector<double>ctt(m*m),xt(m*p);candidate->cross.resize(m);candidate->ctt=full.covariance()[t*n+t];
    r.stage=BoxHeldoutStage::block_extraction;
    for(std::size_t i=0;i<m;++i) {
      const auto oi=i<t?i:i+1;candidate->cross[i]=full.covariance()[oi*n+t];
      for(std::size_t j=0;j<m;++j){const auto oj=j<t?j:j+1;ctt[i*m+j]=full.covariance()[oi*n+oj];}
      for(std::size_t j=0;j<p;++j)xt[i*p+j]=x[oi*p+j];
    }
    r.stage=BoxHeldoutStage::training_factor;++r.work.factor_attempts;
    auto train=prepare_gaussian(ctt,MatrixKind::covariance,std::move(train_metadata),
        policy.box.design.maximum_elements,policy.box.design.maximum_forward_sensitivity,full.arithmetic());
    if(train.status()!=DensityStatus::finite)return fail(train.status(),train.numerical_status());++r.work.factors_completed;
    std::vector<double>().swap(ctt);
    r.stage=BoxHeldoutStage::cross_whitening;++r.work.whitening_attempts;
    auto cross=numerics::whiten(detail::RetainedQrAccess::factor(train),candidate->cross,
        policy.box.design.maximum_elements,policy.box.design.maximum_payload_bytes,policy.box.design.maximum_forward_sensitivity);
    if(cross.status!=numerics::Status::ok)throw Refusal{cross.status};++r.work.whitenings_completed;
    candidate->h=std::move(cross.value);r.cross_whitening_available=true;
    r.cross_whitening_backward_residual=cross.backward_residual;r.cross_whitening_rounding_estimate=cross.arithmetic_rounding_estimate;
    r.stage=BoxHeldoutStage::schur_admission;LD squares=0;
    for(const auto v:candidate->h)squares=checked(squares+checked(v*v));
    candidate->schur=checked(static_cast<LD>(candidate->ctt)-squares);
    if(!(candidate->schur>0))throw Refusal{numerics::Status::not_positive_definite,candidate->schur};
    const auto cancellation=checked((std::abs(static_cast<LD>(candidate->ctt))+squares)/candidate->schur);
    const auto sensitivity=checked(cancellation*(2*cross.arithmetic_rounding_estimate+
        (2*static_cast<LD>(m)+4)*std::numeric_limits<LD>::epsilon()));
    r.schur_variance=detail::box_scalar(candidate->schur);r.schur_cancellation_ratio=detail::box_scalar(cancellation);
    r.schur_sensitivity_estimate=detail::box_scalar(sensitivity);r.schur_available=true;
    if(sensitivity>policy.box.design.maximum_forward_sensitivity)throw Refusal{numerics::Status::conditioning_budget_exceeded,sensitivity};
    r.schur_sqrt_enclosure=detail::box_sqrt(r.schur_variance);r.schur_sqrt_available=true;
    candidate->sqrt_schur=checked(std::sqrt(candidate->schur));
    r.stage=BoxHeldoutStage::training_qr;++r.work.training_qr_attempts;
    const auto train_ids=std::span<const std::string>(train.metadata().ordered_ids);
    candidate->training=DesignProfile::prepare(std::move(train),xt,train_ids,std::move(design),policy.box.design);
    const auto whitening_counts=detail::RetainedQrAccess::preparation_whitenings(candidate->training);
    r.work.whitening_attempts+=whitening_counts.first;r.work.whitenings_completed+=whitening_counts.second;
    if(candidate->training.status()!=DensityStatus::finite)return fail(candidate->training.status(),candidate->training.numerical_status());
    ++r.work.training_qr_completed;r.training_qr_available=true;std::vector<double>().swap(xt);
    const auto view=detail::RetainedQrAccess::view(candidate->training);
    const auto b=detail::RetainedQrAccess::normalized_design(candidate->training);
    r.training_preparation_sensitivity=detail::RetainedQrAccess::preparation_sensitivity(candidate->training);
    r.training_triangular_condition=view.triangular_condition;r.training_transpose_condition=view.transpose_condition;
    candidate->appended.resize(p);LD projection_error=0;
    for(std::size_t j=0;j<p;++j) {
      LD projection=0,absolute=0;
      for(std::size_t i=0;i<m;++i){const auto v=candidate->h[i]*b[i*p+j];projection=checked(projection+v);absolute+=std::abs(v);}
      const auto first=checked(static_cast<LD>(candidate->heldout_design[j])/view.scales[j]);
      candidate->appended[j]=checked((first-projection)/candidate->sqrt_schur);
      projection_error=std::max(projection_error,checked((std::abs(first)+absolute)*(m+4)*
          std::numeric_limits<LD>::epsilon()/(1+std::abs(first-projection))));
    }
    r.appended_projection_rounding_estimate=detail::box_scalar(projection_error);r.appended_response_available=true;
    if(projection_error>policy.box.design.maximum_forward_sensitivity)throw Refusal{numerics::Status::conditioning_budget_exceeded,projection_error};
    r.stage=BoxHeldoutStage::update_qr;++r.work.update_qr_attempts;candidate->updated_r.resize(p*p);
    for(std::size_t i=0;i<p;++i)for(std::size_t j=0;j<p;++j)candidate->updated_r[i*p+j]=j<i?0:view.qr[i*p+j];
    std::vector<LD>row(p);for(std::size_t j=0;j<p;++j)row[j]=candidate->appended[view.pivot[j]];
    candidate->cosine.resize(p);candidate->sine.resize(p);
    for(std::size_t j=0;j<p;++j) {
      const auto a=candidate->updated_r[j*p+j],d=row[j],hyp=checked(std::hypot(a,d));
      if(!(hyp>0))throw Refusal{numerics::Status::singular,hyp};
      const auto c=a/hyp,s=d/hyp;candidate->cosine[j]=c;candidate->sine[j]=s;
      for(std::size_t k=j;k<p;++k){const auto u=candidate->updated_r[j*p+k],v=row[k];
        candidate->updated_r[j*p+k]=checked(c*u+s*v);row[k]=checked(-s*u+c*v);}
    }
    const auto conditions=detail::qr_conditions(candidate->updated_r,p);
    if(conditions.status!=numerics::Status::ok)throw Refusal{conditions.status};
    const auto eta=checked(conditions.triangular*(std::numeric_limits<double>::epsilon()+
        (static_cast<LD>(n)+p)*std::numeric_limits<LD>::epsilon()));
    if(eta>std::min(1e-8,policy.box.design.maximum_forward_sensitivity))throw Refusal{numerics::Status::conditioning_budget_exceeded,eta};
    candidate->update_condition=detail::box_scalar(conditions.triangular);
    candidate->update_transpose=detail::box_scalar(conditions.transpose);candidate->update_sensitivity=detail::box_scalar(eta);
    r.update_triangular_condition=candidate->update_condition;r.update_transpose_condition=candidate->update_transpose;
    ++r.work.variance_screen_attempts;
    std::vector<double>weights(p,0);
    for(std::size_t j=0;j<p;++j){weights[j]=1;
      const auto variance=detail::qr_variance(candidate->joint_view(),weights,policy.box.design.maximum_forward_sensitivity);weights[j]=0;
      if(variance.status!=DensityStatus::finite){r.refusal_index=j;throw Refusal{variance.numerical_status};}}
    ++r.work.variance_screens_completed;
    ++r.work.update_qr_completed;r.update_qr_available=true;
    r.full_source_covariance_log_determinant=detail::RetainedQrAccess::factor(full).log_determinant();
    r.training_covariance_log_determinant=detail::RetainedQrAccess::factor(candidate->training).log_determinant();
    candidate->conditional_constant=detail::box_observation_constant(r.full_source_covariance_log_determinant,r.training_covariance_log_determinant);
    r.covariance_determinants_available=true;
    out.impl_=std::move(candidate);r.retained_payload_bound=out.retained_payload_bound();
    if(!r.retained_payload_bound||*r.retained_payload_bound>*peak||*r.retained_payload_bound>policy.box.design.maximum_payload_bytes) {
      out.impl_.reset();throw Refusal{numerics::Status::work_limit};
    }
    r.retained_frames_available=true;r.status=DensityStatus::finite;r.numerical_status=numerics::Status::ok;r.stage=BoxHeldoutStage::complete;
    // Move construction is noexcept and allocates no default metadata. The
    // final local releases the original full factor only after all gates.
    Gaussian consumed(std::move(full));(void)consumed;return out;
  }catch(const Refusal&f){r.refusal_witness=f.witness;return fail(DensityStatus::numerical_failure,f.status);}
  catch(const std::bad_alloc&){return fail(DensityStatus::numerical_failure,numerics::Status::work_limit);}
}
BoxHeldoutBatch GaussianBoxHeldout::evaluate(std::span<const double>pool,
    std::span<const std::string>ids,std::span<const double>values,
    std::span<const BoxHeldoutRequest>requests,BoxHeldoutPolicy policy)const {
  BoxHeldoutBatch out;
  auto fail=[&](DensityStatus ds,numerics::Status ns){out.status=ds;out.numerical_status=ns;return std::move(out);};
  if(!detail::box_supported())return fail(DensityStatus::unsupported_domain,numerics::Status::outside_domain);
  if(!impl_||status()!=DensityStatus::finite)return fail(status(),preparation_.numerical_status);
  const auto&i=*impl_;const auto view=detail::RetainedQrAccess::view(i.training);const auto m=view.n,p=view.p;
  if(ids.size()!=m)return fail(DensityStatus::incompatible_metadata,numerics::Status::invalid_input);
  if(!tightens(policy,i.frozen)||pool.empty()||pool.size()%m||values.empty()||requests.empty())return out;
  const auto g=pool.size()/m;
  if(m+1>policy.box.design.maximum_elements||p>policy.box.design.maximum_elements)
    return fail(DensityStatus::numerical_failure,numerics::Status::work_limit);
  if(g>policy.maximum_training_vectors||values.size()>policy.maximum_candidate_values||requests.size()>policy.maximum_requests)
    return fail(DensityStatus::numerical_failure,numerics::Status::work_limit);
  const auto payload=evaluation_payload_bound(g,requests.size()),retained=retained_payload_bound();
  out.scratch_output_payload_bound=payload;
  if(!payload||!retained||*retained>policy.box.design.maximum_payload_bytes||
      *payload>policy.box.design.maximum_payload_bytes-*retained)
    return fail(DensityStatus::numerical_failure,numerics::Status::work_limit);
  const auto sensitivity=view.triangular_condition*(std::numeric_limits<double>::epsilon()+
      (static_cast<LD>(m)+p)*std::numeric_limits<LD>::epsilon());
  if(preparation_.training_preparation_sensitivity>policy.box.design.maximum_forward_sensitivity||
      preparation_.schur_sensitivity_estimate>policy.box.design.maximum_forward_sensitivity||
      preparation_.appended_projection_rounding_estimate>policy.box.design.maximum_forward_sensitivity||
      i.update_sensitivity>policy.box.design.maximum_forward_sensitivity||
      sensitivity>std::min(1e-8,policy.box.design.maximum_forward_sensitivity))
    return fail(DensityStatus::numerical_failure,numerics::Status::conditioning_budget_exceeded);
  try {
    Accounting identity(0);for(const auto&id:ids)identity.add(id.size(),1);
    if(!identity.result())throw Refusal{numerics::Status::work_limit};
    precharge(out.work,evaluation_work_bound(m+1,p,g,values.size(),requests.size(),*identity.result()),policy.maximum_evaluation_work_units);
    if(!same(ids,i.training.metadata().ordered_ids))return fail(DensityStatus::incompatible_metadata,numerics::Status::invalid_input);
    for(const auto&r:requests)if(r.training_vector_index>=g||r.candidate_value_index>=values.size())return out;
    if(!std::all_of(pool.begin(),pool.end(),[](double v){return std::isfinite(v);})||
        !std::all_of(values.begin(),values.end(),[](double v){return std::isfinite(v);}))
      return fail(DensityStatus::invalid_input,numerics::Status::nonfinite_input);
    out.requests.assign(requests.begin(),requests.end());out.training_normalizations.resize(g);
    out.training_offset_subtraction_rounding_estimates.resize(g);out.training_refusal_witnesses.resize(g);
    out.training_cross_projection_error_estimates.resize(g);out.training_cross_projection_diagnostics_available.resize(g);
    out.training_offset_diagnostics_available.resize(g);out.densities.resize(requests.size());
    out.output_layout_available=true;
    // Exactly 2p+2 wide scalar entries per training vector: QR head p, actual
    // RHS denominator sums p, orthogonal-tail quadratic and h^T whitened RHS.
    std::vector<std::vector<LD>>cache(g);const auto b=detail::RetainedQrAccess::normalized_design(i.training);
    const auto x=detail::RetainedQrAccess::design(i.training);const auto&factor=detail::RetainedQrAccess::factor(i.training);
    auto whiten=[&](auto rhs,BoxHeldoutWork&w) {
      ++w.whitening_attempts;auto result=numerics::whiten(factor,rhs,policy.box.design.maximum_elements,
          policy.box.design.maximum_payload_bytes,policy.box.design.maximum_forward_sensitivity);
      if(result.status==numerics::Status::ok)++w.whitenings_completed;return result;
    };
    for(std::size_t k=0;k<g;++k) {
      ++out.work.training_attempts;auto&normal=out.training_normalizations[k];
      try {
        std::vector<LD>residual(m);LD error=0;
        for(std::size_t j=0;j<m;++j){const auto original=j<i.heldout?j:j+1;
          residual[j]=checked(static_cast<LD>(pool[k*m+j])-i.offsets[original]);
          double reported=0;if(!detail::qr_output(residual[j],reported))throw Refusal{numerics::Status::outside_domain,residual[j]};
          error=std::max(error,rounding(residual[j],reported));}
        out.training_offset_subtraction_rounding_estimates[k]=detail::box_scalar(error);
        out.training_offset_diagnostics_available[k]=1;
        if(error>policy.box.design.maximum_forward_sensitivity)throw Refusal{numerics::Status::conditioning_budget_exceeded,error};
        normal.completion_step=BoxCompletionStep::source_whitening;
        auto wr=whiten(std::span<const LD>(residual),out.work);if(wr.status!=numerics::Status::ok)throw Refusal{wr.status};
        std::vector<LD>column_scales(p,0);LD cross=0,absolute=0;
        for(std::size_t j=0;j<p;++j)for(std::size_t row=0;row<m;++row)column_scales[j]+=std::abs(b[row*p+j]*wr.value[row]);
        for(std::size_t row=0;row<m;++row){const auto term=checked(i.h[row]*wr.value[row]);
          cross=checked(cross+term);absolute=checked(absolute+std::abs(term));}
        const auto dot_error=checked(absolute*(m+4)*std::numeric_limits<LD>::epsilon());
        // Round upward so reporting cannot understate this empirical estimate.
        out.training_cross_projection_error_estimates[k]=std::nextafter(detail::box_scalar(dot_error),std::numeric_limits<double>::infinity());
        if(dot_error==0)out.training_cross_projection_error_estimates[k]=0;
        out.training_cross_projection_diagnostics_available[k]=1;
        detail::qr_transform(view,detail::RetainedQrAccess::reflectors(i.training),wr.value);
        const auto tail=checked(detail::qr_tail(wr.value,p));
        normal.completion_step=BoxCompletionStep::profile_evaluation;
        auto fit=detail::qr_finish_fit(view,std::span<const LD>(residual),std::span<const LD>(wr.value),
            std::span<const LD>(column_scales),[&](std::size_t row,std::size_t col){return x[row*p+col];},
            [&](std::size_t row,std::size_t col){return b[row*p+col];},
            [&](std::span<const double> adjusted){return whiten(adjusted,out.work);},policy.box.design.maximum_forward_sensitivity);
        if(fit.status!=DensityStatus::finite)throw Refusal{fit.numerical_status};
        normalize(normal,view,std::move(fit),tail,preparation_.training_covariance_log_determinant,i.support,policy,out.work);
        cache[k].resize(2*p+2);std::copy_n(wr.value.begin(),p,cache[k].begin());
        std::copy(column_scales.begin(),column_scales.end(),cache[k].begin()+p);
        cache[k][2*p]=tail;cache[k][2*p+1]=cross;++out.work.training_completed;
      }catch(const Refusal&f){out.training_refusal_witnesses[k]=f.witness;withhold(normal,f.status);}
      catch(const std::bad_alloc&){withhold(normal,numerics::Status::work_limit);}
    }
    for(std::size_t q=0;q<requests.size();++q) {
      const auto request=requests[q];auto&result=out.densities[q];result.stage=BoxHeldoutStage::training_normalization;
      const auto before=out.work;
      const auto&normal=out.training_normalizations[request.training_vector_index];
      if(normal.status!=DensityStatus::finite){result.status=normal.status;result.numerical_status=normal.numerical_status;continue;}
      ++out.work.joint_attempts;result.stage=BoxHeldoutStage::joint_normalization;
      try {
        const auto&cached=cache[request.training_vector_index];std::vector<LD>rhs(cached.begin(),cached.begin()+p),scales(p),residual(m+1);
        LD error=0;
        for(std::size_t row=0;row<m;++row){const auto original=row<i.heldout?row:row+1;
          residual[row]=checked(static_cast<LD>(pool[request.training_vector_index*m+row])-i.offsets[original]);}
        residual[m]=checked(static_cast<LD>(values[request.candidate_value_index])-i.offsets[i.heldout]);
        double reported=0;if(!detail::qr_output(residual[m],reported))throw Refusal{numerics::Status::outside_domain,residual[m]};
        error=std::max(static_cast<LD>(out.training_offset_subtraction_rounding_estimates[request.training_vector_index]),rounding(residual[m],reported));
        result.offset_subtraction_rounding_estimate=detail::box_scalar(error);result.offset_diagnostics_available=true;
        if(error>policy.box.design.maximum_forward_sensitivity)throw Refusal{numerics::Status::conditioning_budget_exceeded,error};
        const auto shifted=checked(residual[m]-cached[2*p+1]),last_original=checked(shifted/i.sqrt_schur);LD last=last_original;
        const auto cross_error=checked((out.training_cross_projection_error_estimates[request.training_vector_index]+
            (std::abs(residual[m])+std::abs(cached[2*p+1]))*4*std::numeric_limits<LD>::epsilon())/(1+std::abs(shifted)));
        result.cross_projection_rounding_estimate=detail::box_scalar(cross_error);result.cross_projection_diagnostics_available=true;
        if(cross_error>policy.box.design.maximum_forward_sensitivity)throw Refusal{numerics::Status::conditioning_budget_exceeded,cross_error};
        LD update_error=0;
        for(std::size_t j=0;j<p;++j){const auto old=rhs[j],c=i.cosine[j],s=i.sine[j];
          rhs[j]=checked(c*old+s*last);const auto next=checked(-s*old+c*last);
          update_error=std::max(update_error,checked((std::abs(old)+std::abs(last))*4*std::numeric_limits<LD>::epsilon()/(1+std::abs(rhs[j])+std::abs(next))));last=next;}
        result.update_rhs_rounding_estimate=detail::box_scalar(update_error);result.update_diagnostics_available=true;
        if(update_error>policy.box.design.maximum_forward_sensitivity)throw Refusal{numerics::Status::conditioning_budget_exceeded,update_error};
        for(std::size_t j=0;j<p;++j)scales[j]=checked(cached[p+j]+std::abs(i.appended[j]*last_original));
        const auto joint=i.joint_view();auto&joint_normal=result.joint_normalization;
        joint_normal.completion_step=BoxCompletionStep::profile_evaluation;
        auto fit=detail::qr_finish_fit(joint,std::span<const LD>(residual),std::span<const LD>(rhs),std::span<const LD>(scales),
            [&](std::size_t row,std::size_t col){return row<m?x[row*p+col]:i.heldout_design[col];},
            [&](std::size_t row,std::size_t col){return row<m?b[row*p+col]:i.appended[col];},
            [&](std::span<const double>adjusted){
              auto wt=whiten(adjusted.first(m),out.work);if(wt.status!=numerics::Status::ok)return wt;
              LD projection=0,absolute=0;for(std::size_t row=0;row<m;++row){const auto term=checked(i.h[row]*wt.value[row]);
                projection=checked(projection+term);absolute=checked(absolute+std::abs(term));}
              const auto last_adjusted=checked((static_cast<LD>(adjusted[m])-projection)/i.sqrt_schur);
              const auto reconstruction=checked(i.sqrt_schur*last_adjusted+projection);
              const auto denominator=checked(std::abs(i.sqrt_schur*last_adjusted)+absolute+std::abs(static_cast<LD>(adjusted[m])));
              const auto backward=denominator==0?0:checked(std::abs(reconstruction-adjusted[m])/denominator);
              wt.backward_residual=detail::box_scalar(std::max(static_cast<LD>(wt.backward_residual),backward));
              wt.value.reserve(m+1);wt.value.push_back(last_adjusted);
              // The appended triangular row's empirical diagnostic includes
              // projection accumulation, division and retained Schur screens.
              const auto adjusted_error=checked((std::abs(static_cast<LD>(adjusted[m]))+absolute)*(m+4)*
                  std::numeric_limits<LD>::epsilon()/(1+std::abs(static_cast<LD>(adjusted[m])-projection)));
              const auto estimate=checked(adjusted_error+preparation_.schur_sensitivity_estimate);
              wt.arithmetic_rounding_estimate=detail::box_scalar(std::max(static_cast<LD>(wt.arithmetic_rounding_estimate),estimate));
              if(estimate>policy.box.design.maximum_forward_sensitivity||wt.backward_residual>policy.box.design.maximum_forward_sensitivity){wt.status=numerics::Status::conditioning_budget_exceeded;wt.value.clear();}
              return wt;
            },policy.box.design.maximum_forward_sensitivity);
        if(fit.status!=DensityStatus::finite)throw Refusal{fit.numerical_status};
        const auto minimum=checked(cached[2*p]+checked(last*last));
        normalize(joint_normal,joint,std::move(fit),minimum,preparation_.full_source_covariance_log_determinant,i.support,policy,out.work);
        ++out.work.joint_completed;result.stage=BoxHeldoutStage::ratio_composition;++out.work.composition_attempts;
        result.scalar_conditional_log_normalization=i.conditional_constant;result.scalar_conditional_constant_available=true;
        result.log_density=detail::box_density(joint_normal.log_relative_box_integral,normal.log_relative_box_integral,
            i.conditional_constant,policy.maximum_conditional_log_density_width);
        result.density_available=true;result.status=DensityStatus::finite;result.numerical_status=numerics::Status::ok;
        result.stage=BoxHeldoutStage::complete;++out.work.compositions_completed;
      }catch(const Refusal&f){result.refusal_witness=f.witness;result.status=DensityStatus::numerical_failure;
        result.numerical_status=f.status;if(result.joint_normalization.status!=DensityStatus::finite)withhold(result.joint_normalization,f.status);}
      catch(const std::bad_alloc&){result.status=DensityStatus::numerical_failure;
        result.numerical_status=numerics::Status::work_limit;withhold(result.joint_normalization,numerics::Status::work_limit);}
      // Per-request actual attempted/completed deltas, separate from the batch
      // precharge. Shared training work is recorded once by the batch.
      result.work=out.work;result.work.charged_work_units=0;
#define IRRED_DELTA(field) result.work.field-=before.field
      IRRED_DELTA(factor_attempts);IRRED_DELTA(factors_completed);IRRED_DELTA(whitening_attempts);IRRED_DELTA(whitenings_completed);
      IRRED_DELTA(training_qr_attempts);IRRED_DELTA(training_qr_completed);IRRED_DELTA(update_qr_attempts);IRRED_DELTA(update_qr_completed);
      IRRED_DELTA(training_attempts);IRRED_DELTA(training_completed);IRRED_DELTA(joint_attempts);IRRED_DELTA(joint_completed);
      IRRED_DELTA(variance_screen_attempts);IRRED_DELTA(variance_screens_completed);IRRED_DELTA(tail_bound_attempts);IRRED_DELTA(tail_bounds_completed);
      IRRED_DELTA(normalization_attempts);IRRED_DELTA(normalizations_completed);IRRED_DELTA(composition_attempts);IRRED_DELTA(compositions_completed);
#undef IRRED_DELTA
    }
    out.status=DensityStatus::finite;out.numerical_status=numerics::Status::ok;return out;
  }catch(const Refusal&f){return fail(DensityStatus::numerical_failure,f.status);}
  catch(const std::bad_alloc&){return fail(DensityStatus::numerical_failure,numerics::Status::work_limit);}
}
} // namespace irred::statistics
