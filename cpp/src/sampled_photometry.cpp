#include "irred/sampled_photometry.hpp"
#include "photometry_detail.hpp"
#include <algorithm>
#include <array>
#include <cfenv>
namespace irred::photometry {
namespace {
using namespace detail;
Interval sum(Interval a, Interval b) noexcept {
  return {add(a.lower,b.lower).lower,add(a.upper,b.upper).upper};
}
Interval difference(Interval a, Interval b) noexcept {
  return {add(a.lower,-b.upper).lower,add(a.upper,-b.lower).upper};
}
Interval product(Interval a, Interval b) noexcept {
  return {multiply(a.lower,b.lower).lower,multiply(a.upper,b.upper).upper};
}
Interval divide(Interval a, Wide b) noexcept {
  // Outward division; all callers have a nonnegative numerator and positive divisor.
  return {a.lower == 0 ? 0 : std::nextafter(a.lower/b, -INFINITY),
          a.upper == 0 ? 0 : std::nextafter(a.upper/b, INFINITY)};
}
Interval interpolate(Interval x, Interval lo, Interval hi, double y0, double y1) noexcept {
  const auto width = difference(hi,lo), offset = difference(x,lo);
  Interval t{std::max(Wide{0},offset.lower/width.upper),
             std::min(Wide{1},offset.upper/width.lower)};
  t.lower = std::max(Wide{0}, std::nextafter(t.lower,-INFINITY));
  t.upper = std::min(Wide{1}, std::nextafter(t.upper,INFINITY));
  // Convex interpolation keeps all physical coefficients nonnegative;
  // outward subtraction also covers near-zero endpoint weights.
  const Interval reverse{std::max(Wide{0},add(1,-t.upper).lower),
                         std::min(Wide{1},add(1,-t.lower).upper)};
  return sum(scale(reverse,y0),scale(t,y1));
}
Interval knot(const SampledInput &p, std::size_t i, Interval r) noexcept {
  return scale(r,p.spectrum.wavelength_metre[i]);
}
bool close(Interval value) noexcept {
  return value.lower >= 0 && std::isfinite(value.upper) &&
    (value.upper == 0 || (value.lower > 0 &&
     (value.upper-value.lower)/value.lower + 128*std::numeric_limits<Wide>::epsilon() < 1.8e-12L));
}
void checked_store(Outcome &out, Interval value) noexcept {
  if (!close(value)) { fail(out,numerics::Status::conditioning_budget_exceeded); return; }
  store(out,value.lower+(value.upper-value.lower)/2);
}
numerics::Status validate(const SampledInput &p) noexcept {
  for (double v : {p.luminosity_distance_metre,p.redshift,p.collecting_area_square_metre,p.observer_exposure_second})
    if (!std::isfinite(v)) return numerics::Status::nonfinite_input;
  if (!(p.luminosity_distance_metre>0) || p.redshift<0 || p.collecting_area_square_metre<0 || p.observer_exposure_second<0)
    return numerics::Status::outside_domain;
  if (p.spectrum.wavelength_metre.size()<2 || p.passband.wavelength_metre.size()<2 ||
      p.spectrum.wavelength_metre.size()!=p.spectrum.luminosity_watt_per_metre.size() ||
      p.passband.wavelength_metre.size()!=p.passband.optical_transmission.size())
    return numerics::Status::invalid_input;
  for (unsigned axis=0;axis<2;++axis) {
    const auto x=axis?p.passband.wavelength_metre:p.spectrum.wavelength_metre;
    const auto y=axis?p.passband.optical_transmission:p.spectrum.luminosity_watt_per_metre;
    for (std::size_t i=0;i<x.size();++i) {
      if (!std::isfinite(x[i]) || !std::isfinite(y[i])) return numerics::Status::nonfinite_input;
      if (!(x[i]>0) || (i && !(x[i]>x[i-1])) || y[i]<0 || (axis && y[i]>1))
        return numerics::Status::outside_domain;
    }
  }
  return numerics::Status::ok;
}
}
SampledResult evaluate_sampled(const SampledInput &p, SampledPolicy policy) noexcept {
  SampledResult result;
  if (!policy.requested_outputs || (policy.requested_outputs&~7u) ||
      std::numeric_limits<Wide>::digits<64 || std::numeric_limits<Wide>::max_exponent<16384 ||
      std::fegetround()!=FE_TONEAREST) return result;
  const auto ns=p.spectrum.wavelength_metre.size(), nt=p.passband.wavelength_metre.size();
  if (ns>65536 || nt>65536 || ns+nt>65536 || ns+nt>policy.maximum_samples ||
      (ns>=2 && nt>=2 && ns+nt-3>policy.maximum_segments)) {
    result.admission_status=numerics::Status::work_limit;return result;
  }
  result.admission_status=validate(p);
  Outcome *out[]={&result.flux_watt_per_square_metre,&result.energy_joule,&result.expected_photons};
  auto needed=[&](unsigned i){return (policy.requested_outputs&(1u<<i)) && out[i]->availability==Availability::omitted;};
  auto failure=[&](numerics::Status s){for(unsigned i=0;i<3;++i) if(needed(i))fail(*out[i],s);};
  if(result.admission_status!=numerics::Status::ok){failure(result.admission_status);return result;}
  const bool source_zero=std::all_of(p.spectrum.luminosity_watt_per_metre.begin(),p.spectrum.luminosity_watt_per_metre.end(),[](double x){return x==0;});
  const bool transmission_zero=std::all_of(p.passband.optical_transmission.begin(),p.passband.optical_transmission.end(),[](double x){return x==0;});
  if(source_zero)for(unsigned i=0;i<3;++i)if(needed(i))store(*out[i],0);
  if(p.collecting_area_square_metre==0 || p.observer_exposure_second==0 || transmission_zero)
    for(unsigned i=1;i<3;++i)if(needed(i))store(*out[i],0);
  if(!needed(0)&&!needed(1)&&!needed(2))return result;
  const Interval r=add(1,p.redshift);
  std::array<Interval,3> integrals{};
  std::size_t i=0,j=0;
  while(i+1<ns && j+1<nt){
    const auto s0=knot(p,i,r),s1=knot(p,i+1,r);
    const Interval t0{p.passband.wavelength_metre[j],p.passband.wavelength_metre[j]},
                   t1{p.passband.wavelength_metre[j+1],p.passband.wavelength_metre[j+1]};
    const Interval a{std::max(s0.lower,t0.lower),std::max(s0.upper,t0.upper)},
                   b{std::min(s1.lower,t1.lower),std::min(s1.upper,t1.upper)};
    if(b.upper>a.lower){
      const auto width=difference(b,a), source_width=difference(s1,s0);
      if(!(width.lower>0) || !(source_width.lower>0)){
        failure(numerics::Status::conditioning_budget_exceeded);return result;
      }
      const auto l0=interpolate(a,s0,s1,p.spectrum.luminosity_watt_per_metre[i],p.spectrum.luminosity_watt_per_metre[i+1]),
                 l1=interpolate(b,s0,s1,p.spectrum.luminosity_watt_per_metre[i],p.spectrum.luminosity_watt_per_metre[i+1]);
      if(needed(0)) integrals[0]=sum(integrals[0],product(width,divide(sum(l0,l1),2)));
      if(needed(1)||needed(2)){
        const auto q0=interpolate(a,t0,t1,p.passband.optical_transmission[j],p.passband.optical_transmission[j+1]),
                   q1=interpolate(b,t0,t1,p.passband.optical_transmission[j],p.passband.optical_transmission[j+1]);
        const std::array<Interval,3> c{product(l0,q0),divide(sum(product(l0,q1),product(l1,q0)),2),product(l1,q1)};
        if(needed(1))integrals[1]=sum(integrals[1],product(width,divide(sum(sum(c[0],c[1]),c[2]),3)));
        if(needed(2)){
          // Multiply degree-two Bernstein response by linear observed wavelength.
          const auto d0=product(c[0],a),d3=product(c[2],b),
                     d1=divide(sum(product(c[0],b),scale(product(c[1],a),2)),3),
                     d2=divide(sum(scale(product(c[1],b),2),product(c[2],a)),3);
          integrals[2]=sum(integrals[2],product(width,divide(sum(sum(d0,d1),sum(d2,d3)),4)));
        }
      }
    }
    // Never discard a segment whose knot ordering is unresolved. Exact
    // shared knots advance both axes; otherwise outward bounds must order them.
    if(s1.lower==s1.upper && s1.lower==t1.lower){++i;++j;}
    else if(s1.upper<=t1.lower)++i;
    else if(t1.upper<=s1.lower)++j;
    else {failure(numerics::Status::conditioning_budget_exceeded);return result;}
  }
  const Wide propagation=detail::spectral_flux(1,p.luminosity_distance_metre,1+static_cast<Wide>(p.redshift));
  const Interval propagation_bounds{std::nextafter(propagation,-INFINITY),std::nextafter(propagation,INFINITY)};
  const auto collection=multiply(p.collecting_area_square_metre,p.observer_exposure_second);
  for(unsigned k=0;k<3;++k)if(needed(k)){
    auto value=product(integrals[k],propagation_bounds);
    if(k)value=product(value,collection);
    if(k==2)value=divide(value,planck_constant_joule_second*static_cast<Wide>(speed_of_light_m_per_s));
    checked_store(*out[k],value);
  }
  return result;
}
} // namespace irred::photometry
