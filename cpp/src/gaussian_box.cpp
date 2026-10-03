#include "irred/gaussian_box.hpp"
#include "payload_accounting.hpp"
#include <algorithm>
#include <cfenv>
#include <cmath>
#include <limits>
#include <new>
#include <stdexcept>
#include <utility>

namespace irred::statistics {
namespace {
struct Refusal { numerics::Status status; };
using LD = long double;
struct I { LD lo, hi; };
constexpr LD infinity = std::numeric_limits<LD>::infinity();
bool supported() {
  return std::numeric_limits<LD>::is_iec559 &&
      std::numeric_limits<LD>::digits >= 64 &&
      std::numeric_limits<LD>::max_exponent >= 16384 &&
      std::fegetround() == FE_TONEAREST;
}
LD checked(LD x) {
  if (!std::isfinite(x) || std::fpclassify(x) == FP_SUBNORMAL)
    throw Refusal{numerics::Status::outside_domain};
  return x;
}
LD dn(LD x) { checked(x); return x == 0 ? 0 : checked(std::nextafter(x, -infinity)); }
LD up(LD x) { checked(x); return x == 0 ? 0 : checked(std::nextafter(x, infinity)); }
I point(LD x) { return {checked(x), checked(x)}; }
I add(I a, I b) {
  const LD l = a.lo+b.lo, h = a.hi+b.hi;
  if ((l == 0 && a.lo != -b.lo) || (h == 0 && a.hi != -b.hi))
    throw Refusal{numerics::Status::outside_domain};
  return {dn(l), up(h)};
}
I neg(I a) { return {-a.hi,-a.lo}; }
I sub(I a, I b) { return add(a,neg(b)); }
LD product(LD a, LD b) {
  const LD x=a*b;
  if (a != 0 && b != 0 && x == 0)
    throw Refusal{numerics::Status::outside_domain};
  return checked(x);
}
I mul(I a, I b) {
  const LD v[] = {product(a.lo,b.lo),product(a.lo,b.hi),
                  product(a.hi,b.lo),product(a.hi,b.hi)};
  return {dn(*std::min_element(v,v+4)),up(*std::max_element(v,v+4))};
}
I reciprocal(I a) {
  if (a.lo <= 0 && a.hi >= 0)
    throw Refusal{numerics::Status::outside_domain};
  LD l=1/a.hi, h=1/a.lo;
  if (l == 0 || h == 0) throw Refusal{numerics::Status::outside_domain};
  return {dn(l),up(h)};
}
I div(I a, I b) { return mul(a,reciprocal(b)); }
I sqrt_point(LD x) {
  if (!(x > 0)) throw Refusal{numerics::Status::outside_domain};
  // libm supplies only a proposal. Squared outward inequalities prove the
  // enclosure, independently of any asserted libm error bound.
  LD l=dn(std::sqrt(x)),h=up(std::sqrt(x));
  for (unsigned k=0;k<8;++k) {
    if (mul(point(l),point(l)).hi <= x &&
        mul(point(h),point(h)).lo >= x) return {l,h};
    l=dn(l);h=up(h);
  }
  throw Refusal{numerics::Status::conditioning_budget_exceeded};
}
I pi() {
  // Exact decimal witnesses: ...2884 < pi < ...2885. One-ulp widening
  // encloses decimal-to-binary parsing as well as this decimal remainder.
  return {dn(3.141592653589793238462643383279502884L),
          up(3.141592653589793238462643383279502885L)};
}
I log_reduced(I x) {
  // x in [1,2], hence 0<=t<=1/3. Every term is positive.
  const I t=div(sub(x,point(1)),add(x,point(1))), t2=mul(t,t);
  I power=t, sum=point(0);
  for (unsigned k=0;k<64;++k) {
    sum=add(sum,div(power,point(2*k+1)));
    power=mul(power,t2);
  }
  const I remainder=div(power,mul(point(129),sub(point(1),t2)));
  return mul(point(2),{sum.lo,add(sum,remainder).hi});
}
I log_point(LD x) {
  if (!(x > 0)) throw Refusal{numerics::Status::outside_domain};
  int exponent=0;
  LD mantissa=std::frexp(x,&exponent)*2;
  --exponent;
  return add(log_reduced(point(mantissa)),
             mul(point(exponent),log_reduced(point(2))));
}
I log_interval(I x) { return {log_point(x.lo).lo,log_point(x.hi).hi}; }
I exp_negative(I x) {
  if (x.lo < 0 || x.hi > 73)
    throw Refusal{numerics::Status::outside_domain};
  // The exact |z|<=12 node has x<=72; its outward mesh enclosure may
  // slightly exceed72. Division256 keeps the certified Taylor domain below
  // 73/256 even with that enclosure, without clamping a node.
  const I v=div(x,point(256));
  I term=point(1),sum=point(1);
  // Degrees 0..31; t32 and geometric ratio v/33 bound the omitted series.
  for (unsigned k=1;k<32;++k) {
    term=div(mul(term,v),point(k));sum=add(sum,term);
  }
  const I next=div(mul(term,v),point(32));
  const I remainder=div(next,sub(point(1),div(v,point(33))));
  I value=reciprocal({sum.lo,add(sum,remainder).hi});
  for (unsigned k=0;k<8;++k) value=mul(value,value);
  return value;
}
LD tail_upper(LD conservative_distance) {
  if (conservative_distance < 0) return 1;
  // Clamp downward before repeated squaring; this avoids positive-tail
  // underflow while preserving a nonzero majorant for arbitrarily far bounds.
  const LD t=std::min(conservative_distance,LD(32));
  I v=add(point(1),div(mul(point(t),point(t)),point(131072)));
  v=reciprocal(v);
  for (unsigned k=0;k<16;++k) v=mul(v,v);
  return mul(point(.5),v).hi;
}
BoxInterval report(I a) {
  checked(a.lo);checked(a.hi);
  double l=static_cast<double>(a.lo),h=static_cast<double>(a.hi);
  if (!std::isfinite(l)||!std::isfinite(h)||
      (l!=0&&std::fpclassify(l)!=FP_NORMAL)||
      (h!=0&&std::fpclassify(h)!=FP_NORMAL)||
      (a.lo!=0&&l==0)||(a.hi!=0&&h==0))
    throw Refusal{numerics::Status::outside_domain};
  if (static_cast<LD>(l)>a.lo) l=std::nextafter(l,-std::numeric_limits<double>::infinity());
  if (static_cast<LD>(h)<a.hi) h=std::nextafter(h,std::numeric_limits<double>::infinity());
  if ((l!=0&&std::fpclassify(l)!=FP_NORMAL)||
      (h!=0&&std::fpclassify(h)!=FP_NORMAL))
    throw Refusal{numerics::Status::outside_domain};
  return {l,h};
}
double scalar(LD x) {
  checked(x);
  const double d=static_cast<double>(x);
  if (!std::isfinite(d)||(x!=0&&(d==0||std::fpclassify(d)!=FP_NORMAL)))
    throw Refusal{numerics::Status::outside_domain};
  return d;
}
bool valid_policy(const BoxPolicy &p) {
  return std::isfinite(p.maximum_log_probability_width)&&p.maximum_log_probability_width>0&&
      std::isfinite(p.maximum_quantile_width)&&p.maximum_quantile_width>0&&
      std::isfinite(p.cdf_absolute_radius)&&p.cdf_absolute_radius>0&&
      p.maximum_cdf_nodes>=2&&p.maximum_bisections>0&&p.maximum_cdf_evaluations>0&&
      std::isfinite(p.design.maximum_forward_sensitivity)&&p.design.maximum_forward_sensitivity>0;
}
void charge(detail::PayloadAccounting &b,const BoxSupport &s) {
  b.strings(s.ordered_parameter_ids);b.vector(s.lower);b.vector(s.upper);
  b.string(s.parameter_measure);b.string(s.prior_identity);b.string(s.fixed_coordinate_provenance);
}
struct CdfWork {
  const BoxPolicy &policy;
  std::size_t nodes=0,evaluations=0,bisections=0;
  I cdf(LD z,LD radius) {
    if (++evaluations>policy.maximum_cdf_evaluations)
      throw Refusal{numerics::Status::work_limit};
    if (z==0) return point(.5);
    const LD b=std::abs(z);
    if (b>12) throw Refusal{numerics::Status::outside_domain};
    const I p=pi();
    const I norm=reciprocal({sqrt_point(2*p.lo).lo,sqrt_point(2*p.hi).hi});
    std::size_t n=2;
    I h,error;
    for (;;) {
      h=div(point(b),point(n));
      error=mul(div(mul(mul(point(12),point(b)),mul(mul(h,h),mul(h,h))),point(180)),norm);
      if (error.hi<=radius/2 && n+1<=policy.maximum_cdf_nodes) break;
      if (n>=policy.maximum_cdf_nodes/2)
        throw Refusal{numerics::Status::work_limit};
      n*=2;
    }
    if (n+1>SIZE_MAX-nodes) throw Refusal{numerics::Status::work_limit};
    nodes+=n+1;
    I sum=point(0);
    for (std::size_t i=0;i<=n;++i) {
      const I x=mul(h,point(i));
      const I v=exp_negative(div(mul(x,x),point(2)));
      sum=add(sum,mul(point(i==0||i==n?1:i%2?4:2),v));
    }
    I integral=mul(div(mul(h,sum),point(3)),norm);
    integral={dn(integral.lo-error.hi),up(integral.hi+error.hi)};
    I result=z>0?add(point(.5),integral):sub(point(.5),integral);
    result.lo=std::max(LD(0),result.lo);result.hi=std::min(LD(1),result.hi);
    if (result.hi-result.lo>2*radius)
      throw Refusal{numerics::Status::conditioning_budget_exceeded};
    return result;
  }
  I inverse(LD probability,LD tolerance) {
    if (!(probability>0&&probability<1))
      throw Refusal{numerics::Status::outside_domain};
    LD lo=-12,hi=12,radius=policy.cdf_absolute_radius;
    // Endpoints must witness the requested probability rather than assuming
    // the finite inverse search domain covers an arbitrarily extreme request.
    const LD endpoint_tail=tail_upper(12);
    if (endpoint_tail>probability||sub(point(1),point(endpoint_tail)).lo<probability)
      throw Refusal{numerics::Status::conditioning_budget_exceeded};
    for (std::size_t iteration=0;iteration<policy.maximum_bisections;++iteration) {
      if (hi-lo<=tolerance) return {lo,hi};
      ++bisections;
      const LD mid=lo+(hi-lo)/2;
      if (mid==lo||mid==hi) throw Refusal{numerics::Status::outside_domain};
      const I value=cdf(mid,radius);
      if (value.hi<probability) lo=mid;
      else if (value.lo>probability) hi=mid;
      else {
        const LD a=std::max(lo,mid-tolerance/2),b=std::min(hi,mid+tolerance/2);
        if (cdf(a,radius).hi<=probability&&cdf(b,radius).lo>=probability)
          return {a,b};
        radius/=4;
        if (!(radius>0)) throw Refusal{numerics::Status::outside_domain};
      }
    }
    throw Refusal{numerics::Status::work_limit};
  }
};
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
  if (!supported()) {out.status_=DensityStatus::unsupported_domain;return out;}
  if (d.status()!=DensityStatus::finite) {out.status_=d.status();return out;}
  const auto &ids=d.design_metadata().ordered_parameter_ids;
  const auto p=ids.size();
  if (s.ordered_parameter_ids!=ids) {out.status_=DensityStatus::incompatible_metadata;return out;}
  if (!valid_policy(policy)||s.lower.size()!=p||s.upper.size()!=p||
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
  CdfWork work{policy};
  auto fail=[&](DensityStatus status,numerics::Status numerical) {
    if (!out.gaussian_completion_available) out=GaussianBoxResult{};
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
  if (!supported()) return fail(DensityStatus::unsupported_domain,numerics::Status::outside_domain);
  if (status_!=DensityStatus::finite) return fail(status_,numerics::Status::invalid_input);
  const auto p=support_.lower.size(),n=profile_.metadata().ordered_ids.size();
  if (!valid_policy(policy)||request.active_parameter_index>=p||
      !std::isfinite(request.cumulative_probability)||
      !(request.cumulative_probability>0&&request.cumulative_probability<1)) return out;
  const auto bound=evaluation_payload_bound();
  if (!bound||*bound>policy.design.maximum_payload_bytes||n>policy.design.maximum_elements||p>policy.design.maximum_elements)
    return fail(DensityStatus::numerical_failure,numerics::Status::work_limit);
  auto fit=profile_.evaluate(r,ids,policy.design);
  if (fit.status!=DensityStatus::finite) return fail(fit.status,fit.numerical_status);
  try {
    auto wr=numerics::whiten(profile_.gaussian_.factor_,r,policy.design.maximum_elements,
        policy.design.maximum_payload_bytes,policy.design.maximum_forward_sensitivity);
    if (wr.status!=numerics::Status::ok) throw Refusal{wr.status};
    // The same retained reflectors, with no second factorization or C readback.
    for (std::size_t k=0;k<p;++k) {
      LD dot=wr.value[k];
      for (std::size_t i=k+1;i<n;++i) dot+=profile_.qr_[i*p+k]*wr.value[i];
      dot*=profile_.tau_[k];wr.value[k]-=dot;
      for (std::size_t i=k+1;i<n;++i) wr.value[i]-=profile_.qr_[i*p+k]*dot;
    }
    LD q=0,logdet=0;
    for (std::size_t i=p;i<n;++i) q+=wr.value[i]*wr.value[i];
    for (std::size_t j=0;j<p;++j)
      logdet+=2*(std::log(profile_.scales_[j])+std::log(std::abs(profile_.qr_[j*p+j])));
    out.minimum_quadratic=scalar(q);
    out.log_design_precision_determinant=scalar(logdet);
    out.reported_profile_quadratic=fit.quadratic;
    out.profile_stationarity=fit.normalized_normal_equation_residual;
    out.unboxed_mean=std::move(fit.coefficients);out.unboxed_variance.resize(p);
    std::vector<double> weights(p,0);
    for (std::size_t j=0;j<p;++j) {
      weights[j]=1;
      LinearFunctionalMetadata fm;fm.functional_identity="box marginal variance/"+support_.ordered_parameter_ids[j];
      fm.output_unit=profile_.design_metadata().parameter_units[j];fm.weight_units.assign(p,"declared coordinate selector");
      const auto v=profile_.estimator_variance(weights,support_.ordered_parameter_ids,std::move(fm),policy.design);
      weights[j]=0;
      if (v.status!=DensityStatus::finite) throw Refusal{v.numerical_status};
      out.unboxed_variance[j]=v.variance;
      out.maximum_variance_sensitivity=std::max(out.maximum_variance_sensitivity,v.estimated_forward_sensitivity);
    }
    std::vector<double> zero(n,0);
    const auto source=profile_.gaussian_.evaluate(zero,ids,policy.design.maximum_forward_sensitivity);
    if (source.density.status!=DensityStatus::finite) throw Refusal{source.density.numerical_status};
    out.source_covariance_log_determinant=source.log_determinant;
    out.gaussian_completion_available=true;out.stage=BoxStage::gaussian_completion;
    I excluded=point(0),volume=point(0);
    out.standardized_lower_margins.resize(p);out.standardized_upper_margins.resize(p);
    for (std::size_t j=0;j<p;++j) {
      const I sigma=sqrt_point(out.unboxed_variance[j]);
      const I lower=div(sub(point(out.unboxed_mean[j]),point(support_.lower[j])),sigma);
      const I upper=div(sub(point(support_.upper[j]),point(out.unboxed_mean[j])),sigma);
      out.standardized_lower_margins[j]=report(lower);out.standardized_upper_margins[j]=report(upper);
      excluded=add(excluded,point(tail_upper(lower.lo)));
      excluded=add(excluded,point(tail_upper(upper.lo)));
      volume=add(volume,log_interval(sub(point(support_.upper[j]),point(support_.lower[j]))));
    }
    const LD e=excluded.hi;
    out.excluded_mass_upper=report({0,e}).upper;
    out.log_prior_volume=report(volume);
    out.endpoint_margins_available=true;out.stage=BoxStage::endpoint_margins;
    if (!(e<1)) throw Refusal{numerics::Status::conditioning_budget_exceeded};
    // -log(1-e)<=e/(1-e): no transcendental error is hidden in the tail gate.
    const I tail_log={-div(point(e),sub(point(1),point(e))).hi,0};
    if (-tail_log.lo>policy.maximum_log_probability_width)
      throw Refusal{numerics::Status::conditioning_budget_exceeded};
    out.box_probability=report({sub(point(1),point(e)).lo,1});
    out.log_box_probability=report(tail_log);
    out.rectangle_enclosure_available=true;out.stage=BoxStage::rectangle_enclosure;
    const I log2pi=log_interval(mul(point(2),pi()));
    const I unboxed=sub(add(point(-out.minimum_quadratic/2),mul(point(static_cast<LD>(p)/2),log2pi)),
                         point(out.log_design_precision_determinant/2));
    const I relative=add(unboxed,tail_log),prior=sub(relative,volume);
    out.log_relative_box_integral=report(relative);
    out.log_prior_normalized_relative_evidence=report(prior);
    out.log_observation_normalized_evidence=report(sub(prior,mul(point(.5),add(
        point(out.source_covariance_log_determinant),mul(point(n),log2pi)))));
    out.normalization_enclosures_available=true;out.stage=BoxStage::normalization_enclosures;
    const auto j=request.active_parameter_index;
    const I sigma=sqrt_point(out.unboxed_variance[j]);
    const I probability=point(request.cumulative_probability);
    const I low_probability=mul(probability,sub(point(1),point(e)));
    const I high_probability=add(probability,mul(sub(point(1),probability),point(e)));
    const LD tolerance=policy.maximum_quantile_width/(4*sigma.hi);
    const I low=work.inverse(low_probability.lo,tolerance),high=work.inverse(high_probability.hi,tolerance);
    I quantile={add(point(out.unboxed_mean[j]),mul(sigma,point(low.lo))).lo,
                add(point(out.unboxed_mean[j]),mul(sigma,point(high.hi))).hi};
    quantile.lo=std::max(quantile.lo,static_cast<LD>(support_.lower[j]));
    quantile.hi=std::min(quantile.hi,static_cast<LD>(support_.upper[j]));
    if (quantile.lo>quantile.hi||quantile.hi-quantile.lo>policy.maximum_quantile_width)
      throw Refusal{numerics::Status::conditioning_budget_exceeded};
    out.requested_quantile=report(quantile);
    if (static_cast<LD>(out.requested_quantile.upper)-out.requested_quantile.lower>
        policy.maximum_quantile_width)
      throw Refusal{numerics::Status::conditioning_budget_exceeded};
    out.quantile_enclosure_available=true;out.stage=BoxStage::quantile_enclosure;
    auto conditional_cdf=[&](LD x) {
      const I z=div(sub(point(x),point(out.unboxed_mean[j])),sigma);
      const I f={work.cdf(z.lo,policy.cdf_absolute_radius).lo,
                 work.cdf(z.hi,policy.cdf_absolute_radius).hi};
      return I{std::max(LD(0),div(sub(f,point(e)),sub(point(1),point(e))).lo),
               std::min(LD(1),div(f,sub(point(1),point(e))).hi)};
    };
    out.lower_endpoint_cdf=report(conditional_cdf(out.requested_quantile.lower));
    out.upper_endpoint_cdf=report(conditional_cdf(out.requested_quantile.upper));
    out.endpoint_cdf_enclosures_available=true;out.stage=BoxStage::complete;
    out.cdf_node_evaluations=work.nodes;out.cdf_evaluations=work.evaluations;out.bisections=work.bisections;
    out.status=DensityStatus::finite;out.numerical_status=numerics::Status::ok;
    return out;
  } catch (const Refusal &f) { return fail(DensityStatus::numerical_failure,f.status); }
}
} // namespace irred::statistics
