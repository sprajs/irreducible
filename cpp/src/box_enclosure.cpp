#include "box_enclosure.hpp"
#include <algorithm>
#include <cfenv>
#include <cmath>
#include <limits>
namespace irred::detail {
using namespace statistics;
namespace {
using Refusal = BoxRefusal;
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
    throw Refusal{numerics::Status::outside_domain,x};
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
bool box_supported() noexcept { return supported(); }
bool box_valid_policy(const BoxPolicy &p) noexcept { return valid_policy(p); }
double box_scalar(long double x) { return scalar(x); }
BoxInterval box_sqrt(double x) { return report(sqrt_point(x)); }
long double box_normalize(GaussianBoxResult &out,const BoxSupport &support,const BoxPolicy &policy,
    std::size_t n) {
  const auto p=support.lower.size();
    I excluded=point(0),volume=point(0);
    out.standardized_lower_margins.resize(p);out.standardized_upper_margins.resize(p);
    for (std::size_t j=0;j<p;++j) {
      const I sigma=sqrt_point(out.unboxed_variance[j]);
      const I lower=div(sub(point(out.unboxed_mean[j]),point(support.lower[j])),sigma);
      const I upper=div(sub(point(support.upper[j]),point(out.unboxed_mean[j])),sigma);
      out.standardized_lower_margins[j]=report(lower);out.standardized_upper_margins[j]=report(upper);
      excluded=add(excluded,point(tail_upper(lower.lo)));
      excluded=add(excluded,point(tail_upper(upper.lo)));
      volume=add(volume,log_interval(sub(point(support.upper[j]),point(support.lower[j]))));
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
  return e;
}
void box_quantile(GaussianBoxResult &out,const BoxSupport &support,BoxMarginalRequest request,
    const BoxPolicy &policy,BoxCdfCounters &counts,long double e) {
  CdfWork work{policy};
  struct Receipt { CdfWork &w; BoxCdfCounters &c; ~Receipt(){c={w.nodes,w.evaluations,w.bisections};} } receipt{work,counts};
    const auto j=request.active_parameter_index;
    const I sigma=sqrt_point(out.unboxed_variance[j]);
    const I probability=point(request.cumulative_probability);
    const I low_probability=mul(probability,sub(point(1),point(e)));
    const I high_probability=add(probability,mul(sub(point(1),probability),point(e)));
    const LD tolerance=policy.maximum_quantile_width/(4*sigma.hi);
    const I low=work.inverse(low_probability.lo,tolerance),high=work.inverse(high_probability.hi,tolerance);
    I quantile={add(point(out.unboxed_mean[j]),mul(sigma,point(low.lo))).lo,
                add(point(out.unboxed_mean[j]),mul(sigma,point(high.hi))).hi};
    quantile.lo=std::max(quantile.lo,static_cast<LD>(support.lower[j]));
    quantile.hi=std::min(quantile.hi,static_cast<LD>(support.upper[j]));
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
}
BoxInterval box_observation_constant(double full,double training) {
  return report(mul(point(-.5),add(sub(point(full),point(training)),
      log_interval(mul(point(2),pi())))));
}
BoxInterval box_density(BoxInterval joint,BoxInterval training,BoxInterval constant,double width) {
  const auto value=report(add(sub({joint.lower,joint.upper},{training.lower,training.upper}),
      {constant.lower,constant.upper}));
  if (static_cast<LD>(value.upper)-value.lower>width)
    throw BoxRefusal{numerics::Status::conditioning_budget_exceeded,{}};
  return value;
}
} // namespace irred::detail
