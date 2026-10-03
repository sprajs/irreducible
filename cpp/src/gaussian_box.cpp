#include "irred/gaussian_box.hpp"
#include "payload_accounting.hpp"
#include "box_enclosure.hpp"
#include "retained_gaussian.hpp"
#include <algorithm>
#include <cfenv>
#include <cmath>
#include <limits>
#include <new>
#include <stdexcept>
#include <utility>

namespace irred::statistics {
namespace {
using LD=long double;
using Refusal=detail::BoxRefusal;
void charge(detail::PayloadAccounting &b,const BoxSupport &s) {
  b.strings(s.ordered_parameter_ids);b.vector(s.lower);b.vector(s.upper);
  b.string(s.parameter_measure);b.string(s.prior_identity);b.string(s.fixed_coordinate_provenance);
}
} // namespace

GaussianBox::GaussianBox(GaussianBox &&o) noexcept
    : profile_(std::move(o.profile_)),support_(std::move(o.support_)),
      status_(std::exchange(o.status_,DensityStatus::invalid_input)) {}
GaussianBox &GaussianBox::operator=(GaussianBox &&o) noexcept {
  if (this!=&o) {this->~GaussianBox();new(this) GaussianBox(std::move(o));}
  return *this;
}
std::optional<std::size_t> GaussianBox::preparation_payload_bound(
    const DesignProfile &d,const BoxSupport &s) noexcept {
  detail::PayloadAccounting b(sizeof(GaussianBox));
  b.embedded(d.retained_payload_bound(),sizeof(DesignProfile));charge(b,s);
  return b.result();
}
std::optional<std::size_t> GaussianBox::retained_payload_bound() const noexcept {
  return preparation_payload_bound(profile_,support_);
}
std::optional<std::size_t> GaussianBox::evaluation_payload_bound() const noexcept {
  detail::PayloadAccounting b(sizeof(GaussianBoxResult));
  const auto n=profile_.metadata().ordered_ids.size(),p=support_.lower.size();
  b.embedded(profile_.evaluation_payload_bound(),0);
  b.add(n,4*sizeof(LD));b.add(p,20*sizeof(LD)+8*sizeof(double)+2*sizeof(BoxInterval));
  // Only one variance metadata object is live at once. Arbitrarily long
  // declared IDs/units still count; a fixed per-coordinate charge alone would
  // not bound these copied strings.
  std::size_t largest_id=0,largest_unit=0;
  for (const auto &id:support_.ordered_parameter_ids) largest_id=std::max(largest_id,id.capacity());
  for (const auto &unit:profile_.design_metadata().parameter_units) largest_unit=std::max(largest_unit,unit.capacity());
  b.add(p,sizeof(std::string)+64);b.add(1,largest_id);b.add(1,largest_unit);
  b.add(1,sizeof(LinearFunctionalMetadata)+128);
  return b.result();
}
GaussianBox GaussianBox::prepare(DesignProfile &&d,BoxSupport s,BoxPolicy policy) {
  GaussianBox out;
  if (!detail::box_supported()) {out.status_=DensityStatus::unsupported_domain;return out;}
  if (d.status()!=DensityStatus::finite) {out.status_=d.status();return out;}
  const auto &ids=d.design_metadata().ordered_parameter_ids;
  const auto p=ids.size();
  if (s.ordered_parameter_ids!=ids) {out.status_=DensityStatus::incompatible_metadata;return out;}
  if (!detail::box_valid_policy(policy)||s.lower.size()!=p||s.upper.size()!=p||
      s.parameter_measure.empty()||s.prior_identity.empty()||s.fixed_coordinate_provenance.empty()) return out;
  for (std::size_t j=0;j<p;++j)
    if (!std::isfinite(s.lower[j])||!std::isfinite(s.upper[j])||!(s.lower[j]<s.upper[j])||
        (s.lower[j]!=0&&std::fpclassify(s.lower[j])!=FP_NORMAL)||
        (s.upper[j]!=0&&std::fpclassify(s.upper[j])!=FP_NORMAL)) return out;
  const auto bound=preparation_payload_bound(d,s);
  if (!bound||*bound>policy.design.maximum_payload_bytes||p>policy.design.maximum_elements) {
    out.status_=DensityStatus::numerical_failure;return out;
  }
  out.profile_=std::move(d);out.support_=std::move(s);out.status_=DensityStatus::finite;
  return out;
}
GaussianBoxResult GaussianBox::evaluate(std::span<const double> r,
    std::span<const std::string> ids,BoxMarginalRequest request,BoxPolicy policy) const {
  GaussianBoxResult out;
  detail::BoxCdfCounters work;
  auto fail=[&](DensityStatus status,numerics::Status numerical) {
    if (!out.gaussian_completion_available) {
      const auto step=out.completion_step;
      const auto coordinate=out.completion_parameter_index;
      out=GaussianBoxResult{};
      out.completion_step=step;out.completion_parameter_index=coordinate;
    }
    if (!out.endpoint_margins_available) {
      out.standardized_lower_margins.clear();out.standardized_upper_margins.clear();
      out.excluded_mass_upper=0;out.log_prior_volume={};
    }
    if (!out.rectangle_enclosure_available) {out.box_probability={};out.log_box_probability={};}
    if (!out.normalization_enclosures_available) {
      out.log_relative_box_integral={};out.log_prior_normalized_relative_evidence={};
      out.log_observation_normalized_evidence={};
    }
    if (!out.quantile_enclosure_available) out.requested_quantile={};
    if (!out.endpoint_cdf_enclosures_available) {out.lower_endpoint_cdf={};out.upper_endpoint_cdf={};}
    out.status=status;out.numerical_status=numerical;
    out.cdf_node_evaluations=work.nodes;out.cdf_evaluations=work.evaluations;out.bisections=work.bisections;
    return std::move(out);
  };
  if (!detail::box_supported()) return fail(DensityStatus::unsupported_domain,numerics::Status::outside_domain);
  if (status_!=DensityStatus::finite) return fail(status_,numerics::Status::invalid_input);
  const auto p=support_.lower.size(),n=profile_.metadata().ordered_ids.size();
  if (!detail::box_valid_policy(policy)||request.active_parameter_index>=p||
      !std::isfinite(request.cumulative_probability)||
      !(request.cumulative_probability>0&&request.cumulative_probability<1)) return out;
  const auto bound=evaluation_payload_bound();
  if (!bound||*bound>policy.design.maximum_payload_bytes||n>policy.design.maximum_elements||p>policy.design.maximum_elements)
    return fail(DensityStatus::numerical_failure,numerics::Status::work_limit);
  out.completion_step=BoxCompletionStep::profile_evaluation;
  auto fit=profile_.evaluate(r,ids,policy.design);
  if (fit.status!=DensityStatus::finite) return fail(fit.status,fit.numerical_status);
  try {
    out.completion_step=BoxCompletionStep::source_whitening;
    auto wr=numerics::whiten(profile_.gaussian_.factor_,r,policy.design.maximum_elements,
        policy.design.maximum_payload_bytes,policy.design.maximum_forward_sensitivity);
    if (wr.status!=numerics::Status::ok) throw Refusal{wr.status};
    const auto view=detail::RetainedQrAccess::view(profile_);
    detail::qr_transform(view,detail::RetainedQrAccess::reflectors(profile_),wr.value);
    const auto &factor=detail::RetainedQrAccess::factor(profile_);
    if(factor.status()!=numerics::Status::ok||factor.size()!=n)
      throw Refusal{numerics::Status::invalid_input};
    detail::qr_complete(out,view,std::move(fit),detail::qr_tail(wr.value,p),
        factor.log_determinant(),policy.design.maximum_forward_sensitivity);
    const auto excluded=detail::box_normalize(out,support_,policy,n);
    detail::box_quantile(out,support_,request,policy,work,excluded);
    out.stage=BoxStage::complete;
    out.cdf_node_evaluations=work.nodes;out.cdf_evaluations=work.evaluations;out.bisections=work.bisections;
    out.status=DensityStatus::finite;out.numerical_status=numerics::Status::ok;
    return out;
  } catch (const Refusal &f) { return fail(DensityStatus::numerical_failure,f.status); }
}
} // namespace irred::statistics
