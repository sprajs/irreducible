#include "irred/thermal_clocks.hpp"
#include "flat_geometry.hpp"
#include "payload_accounting.hpp"
#include "thermal_conformal_epoch.hpp"
#include <algorithm>
#include <cfenv>
#include <cmath>
#include <limits>
#include <utility>
namespace irred::cosmology {
namespace {
using S=numerics::Status;
using W=long double;
constexpr W arithmetic=64.L*std::numeric_limits<double>::epsilon();
constexpr W gyr_seconds=365.25L*86400*1000000000;
constexpr unsigned bit(ThermalClockOutput x) { return thermal_clock_mask(x); }
constexpr unsigned conformal_mask=bit(ThermalClockOutput::comoving_particle_horizon_mpc)|
    bit(ThermalClockOutput::proper_particle_horizon_mpc);
bool tolerance(double a,double r) {
  return std::isfinite(a)&&std::isfinite(r)&&a>=0&&r>=0&&(a>0||r>0);
}
bool valid(const ThermalClockPolicy &p) {
  return std::numeric_limits<W>::digits>=64&&std::numeric_limits<W>::max_exponent>=16384&&
      std::fegetround()==FE_TONEAREST&&p.maximum_depth<=60&&p.maximum_tail_refinements<=512&&
      tolerance(p.absolute_tolerance_gyr,p.relative_tolerance)&&
      tolerance(p.absolute_tolerance_mpc,p.relative_tolerance)&&
      tolerance(p.thermal.absolute_tolerance,p.thermal.relative_tolerance)&&
      p.thermal.maximum_depth<=60&&
      (p.thermal.momentum_method==ThermalMomentumMethod::direct_adaptive||
       p.thermal.momentum_method==ThermalMomentumMethod::nested_clenshaw_curtis);
}
bool valid_source(ThermalStressSourceEstimate p) {
  return std::isfinite(p.relative_component_estimate)&&p.relative_component_estimate>=0&&
      std::isfinite(p.radiation_absolute_estimate)&&p.radiation_absolute_estimate>=0;
}
W source_radius(W d,W a,W lambda,W today,ThermalStressSourceEstimate p) {
  const W a2=a*a,a4=a2*a2;
  return p.relative_component_estimate*(std::max(0.L,d-lambda*a4)+today*a4)+
      p.radiation_absolute_estimate*(1+a4);
}
struct ClockWork {
  std::size_t outer=0,momentum=0;
  std::size_t total() const { return outer+momentum; }
};
struct Context {
  const ThermalBackground &background;
  const ThermalClockPolicy &policy;
  ThermalStressSourceEstimate source;
  W lambda,today;
  ClockWork &work;
  std::size_t start;
  W lower=0,upper=0,dependency=0,cast_fraction=0;
  bool age=false,dust=false,squared=false;
  S status=S::ok;
  bool available() {
    if(work.total()>=policy.maximum_total_callbacks||
       work.total()-start>=policy.maximum_callbacks_per_point||
       work.momentum>policy.thermal.maximum_total_callbacks) {
      status=S::work_limit; return false;
    }
    return true;
  }
  ThermalScaledExpansion expansion(W a) {
    ThermalScaledExpansion result;
    if(status!=S::ok) { result.status=status; return result; }
    if(!available()) { result.status=status; return result; }
    // Tail endpoint and ordinary integrand evaluations are both charged.
    ++work.outer;
    auto nested=policy.thermal;
    nested.maximum_total_callbacks=std::min({policy.maximum_total_callbacks-work.total(),
        policy.maximum_callbacks_per_point-(work.total()-start),
        policy.thermal.maximum_total_callbacks-work.momentum});
    result=background.scaled_expansion(a,nested);work.momentum+=result.callbacks;
    if(result.status!=S::ok) { status=result.status; return result; }
    result.error_estimate+=source_radius(result.a4_e2,a,lambda,today,source);
    if(!(result.a4_e2>result.error_estimate)) {
      result.status=status=S::conditioning_budget_exceeded;return result;
    }
    const W lo=result.a4_e2-result.error_estimate;
    const W radius=result.error_estimate/(std::sqrt(lo)*(std::sqrt(result.a4_e2)+std::sqrt(lo)));
    dependency=std::max(dependency,radius);
    return result;
  }
};
double integrand(double parameter,const void *pointer) {
  auto &c=*const_cast<Context *>(static_cast<const Context *>(pointer));
  if(c.status!=S::ok) return std::numeric_limits<double>::quiet_NaN();
  const W coordinate=c.lower+(c.upper-c.lower)*parameter;
  W value=0;
  if(c.dust) {
    if(!c.available()) return std::numeric_limits<double>::quiet_NaN();
    ++c.work.outer;
    const auto &source=c.background.source();
    const W matter=source.omega_b+static_cast<W>(source.omega_cdm);
    const W u2=coordinate*coordinate,u6=u2*u2*u2;
    const W denominator=matter+c.lambda*u6;
    const W error=arithmetic*denominator+c.source.relative_component_estimate*
        (matter+(matter+c.lambda)*u6);
    if(!(denominator>error)) {
      c.status=S::conditioning_budget_exceeded;return std::numeric_limits<double>::quiet_NaN();
    }
    const W lo=denominator-error;
    c.dependency=std::max(c.dependency,error/(std::sqrt(lo)*(std::sqrt(denominator)+std::sqrt(lo))));
    value=2*(c.age?u2:1)/std::sqrt(denominator);
  } else {
    const W a=c.squared?coordinate*coordinate:coordinate;
    const auto d=c.expansion(a);
    if(d.status!=S::ok) return std::numeric_limits<double>::quiet_NaN();
    value=(c.age?a:1)/std::sqrt(d.a4_e2);
    if(c.squared)value*=2*coordinate;
  }
  value*=c.upper-c.lower;
  const double stored=static_cast<double>(value);
  if(!detail::physical_representable(value)||!std::isfinite(stored)||(value>0&&!(stored>0))) {
    c.status=S::outside_domain;return std::numeric_limits<double>::quiet_NaN();
  }
  if(value>0)c.cast_fraction=std::max(c.cast_fraction,std::abs(value-static_cast<W>(stored))/value);
  return stored;
}
struct Integral { S status=S::invalid_input;W value=0,error=0; };
Integral integrate(Context &c,W absolute,W scale) {
  Integral result;
  if(c.lower==c.upper) { result.status=S::ok;return result; }
  if(!c.available()) { result.status=c.status;return result; }
  const auto remaining=std::min(c.policy.maximum_total_callbacks-c.work.total(),
      c.policy.maximum_callbacks_per_point-(c.work.total()-c.start));
  if(remaining<3) { result.status=S::work_limit;return result; }
  const double abs=static_cast<double>(std::min(absolute/(8*scale),
      static_cast<W>(std::numeric_limits<double>::max())));
  const double rel=c.policy.relative_tolerance/8;
  if(abs==0&&rel==0) { result.status=S::conditioning_budget_exceeded;return result; }
  const auto q=numerics::integrate(integrand,&c,0,1,{abs,rel,remaining,c.policy.maximum_depth});
  result.status=c.status==S::ok?q.status:c.status;
  if(result.status!=S::ok)return result;
  if(!(c.cast_fraction<1)) { result.status=S::conditioning_budget_exceeded;return result; }
  const W response=(std::abs(static_cast<W>(q.value))+q.error_estimate)/(1-c.cast_fraction);
  result.value=q.value;
  result.error=q.error_estimate+response*(c.dependency+c.cast_fraction)+arithmetic*response;
  return result;
}
ThermalClockTail tail(W endpoint,W low,W high,bool age,unsigned refinements) {
  const W numerator=age?endpoint*endpoint/2:endpoint;
  const W lower=numerator/std::sqrt(high),upper=numerator/std::sqrt(low);
  return {endpoint,(lower+upper)/2,(upper-lower)/2+arithmetic*(lower+upper),refinements};
}
void admit(ThermalBackgroundValue &out,Integral value,W scale,double absolute,double relative) {
  out.status=value.status;if(value.status!=S::ok)return;
  const W center=value.value*scale,error=value.error*scale+arithmetic*std::abs(center);
  const double stored=static_cast<double>(center);
  const W estimate=error+std::abs(center-static_cast<W>(stored));
  if(center<0||!detail::physical_representable(center)||!detail::physical_representable(estimate)) {
    out.status=S::outside_domain;return;
  }
  if(estimate>absolute+relative*center) { out.status=S::conditioning_budget_exceeded;return; }
  out={S::ok,stored,static_cast<double>(estimate)};
}
} // namespace
std::optional<std::size_t> thermal_clocks_payload_bound(std::size_t points,std::size_t count) noexcept {
  irred::detail::PayloadAccounting bytes(sizeof(ThermalClocks)+sizeof(ThermalClockBatch));
  bytes.add(points,2*sizeof(ThermalClockRow));bytes.add(count,3*sizeof(ThermalSpecies));
  const auto background=thermal_background_payload_bound(0,count);if(!background)return {};
  bytes.add(1,*background+2048);return bytes.result();
}
ThermalClocks prepare_thermal_clocks(const ThermalBackground &background,
    ThermalClockPolicy policy,ThermalStressSourceEstimate source) {
  ThermalClocks result;
  if(!valid(policy)||!valid_source(source)||policy.thermal.momentum_method!=background.momentum_method())return result;
  if(background.status()!=S::ok) { result.status_=background.status();return result; }
  const auto bytes=thermal_clocks_payload_bound(0,background.source().species.size());
  if(!bytes)return result;
  if(*bytes>policy.maximum_native_bytes||background.source().species.size()>policy.thermal.maximum_species) {
    result.status_=S::work_limit;return result;
  }
  const auto coefficient=detail::ThermalRetainedCoefficientAccess::capture(background);
  if(!coefficient)return result;
  const auto radiation=background.scaled_expansion(0,policy.thermal);
  if(radiation.status!=S::ok) { result.status_=radiation.status;return result; }
  result.background_=background;result.source_estimate_=source;
  result.lambda_=coefficient->lambda_retained;result.omega_species_today_=coefficient->omega_species_today;
  result.radiation_=radiation.a4_e2;
  result.radiation_error_=radiation.error_estimate+
      source_radius(radiation.a4_e2,0,result.lambda_,0,source);
  result.status_=S::ok;return result;
}
ThermalClocks::ThermalClocks(ThermalClocks &&other) noexcept { *this=std::move(other); }
ThermalClocks &ThermalClocks::operator=(ThermalClocks &&other) noexcept {
  if(this==&other)return *this;
  status_=other.status_;background_=std::move(other.background_);source_estimate_=other.source_estimate_;
  radiation_=other.radiation_;radiation_error_=other.radiation_error_;
  lambda_=other.lambda_;omega_species_today_=other.omega_species_today_;
  other.status_=S::invalid_input;other.source_estimate_={};
  other.radiation_=other.radiation_error_=other.lambda_=other.omega_species_today_=0;return *this;
}
ThermalClockBatch ThermalClocks::evaluate(std::span<const double> factors,unsigned requested,
    ThermalClockPolicy policy) const {
  ThermalClockBatch result;
  if(!valid(policy)||requested==0||(requested&~((1u<<thermal_clock_output_count)-1))||
      factors.size()>65536||policy.thermal.momentum_method!=background_.momentum_method())return result;
  if(status_!=S::ok) { result.status=status_;return result; }
  const auto count=background_.source().species.size();const auto bytes=thermal_clocks_payload_bound(factors.size(),count);
  if(!bytes)return result;
  if(factors.size()>policy.maximum_points||count>policy.thermal.maximum_species||*bytes>policy.maximum_native_bytes) {
    result.status=S::work_limit;return result;
  }
  const auto scale=detail::prepare_flat_scale(background_.source().h0_km_s_mpc);
  if(scale.status!=Status::ok)return result;
  const W age_scale=scale.time_seconds/gyr_seconds,horizon_scale=scale.distance_mpc;
  const auto &source=background_.source();
  const W matter=source.omega_b+static_cast<W>(source.omega_cdm);
  const W today=source.omega_gamma+static_cast<W>(source.omega_massless_nonphoton)+matter+omega_species_today_;
  ClockWork work;result.rows.reserve(factors.size());result.requested_outputs=requested;result.status=S::ok;
  for(double a:factors) {
    result.rows.emplace_back();auto &row=result.rows.back();row.scale_factor=a;
    const auto start=work.total(),outer_start=work.outer,momentum_start=work.momentum;
    auto fail=[&](unsigned mask,S status) { for(unsigned i=0;i<thermal_clock_output_count;++i)
      if((requested&mask)&(1u<<i))row.outputs[i].status=status; };
    if(!std::isfinite(a)) { fail(requested,S::nonfinite_input);continue; }
    if(a<0||a>1) { fail(requested,S::outside_domain);continue; }
    Context context{background_,policy,source_estimate_,lambda_,today,work,start};
    const bool need_age=requested&bit(ThermalClockOutput::age_gyr),need_conformal=requested&conformal_mask;
    // If only the proper horizon is requested, allocate its own scale rather
    // than imposing an unrequested comoving-horizon precision.
    const W conformal_scale=horizon_scale*
        ((requested&bit(ThermalClockOutput::comoving_particle_horizon_mpc))?1:a);
    Integral age,conformal;
    if(need_age||need_conformal) {
      if(radiation_==0&&matter==0) {
        age.status=conformal.status=S::outside_domain;
      } else if(a==0) {
        age.status=conformal.status=S::ok;
      } else if(radiation_==0) {
        if(radiation_error_>0||!source.species.empty())age.status=conformal.status=S::conditioning_budget_exceeded;
        else {
          auto dust_integral=[&](bool proper_age,W absolute,W scale) {
            Context interval{background_,policy,source_estimate_,lambda_,today,work,start};
            interval.dust=true;interval.age=proper_age;
            interval.lower=0;interval.upper=std::sqrt(static_cast<W>(a));
            return integrate(interval,absolute,scale);
          };
          if(need_age)age=dust_integral(true,policy.absolute_tolerance_gyr,age_scale);
          if(need_conformal)conformal=dust_integral(false,policy.absolute_tolerance_mpc,conformal_scale);
        }
      } else if(!(radiation_>radiation_error_)) {
        age.status=conformal.status=S::conditioning_budget_exceeded;
      } else {
        W endpoint=a;bool earned=false;
        for(unsigned refinement=0;refinement<=policy.maximum_tail_refinements;++refinement) {
          const auto d=context.expansion(endpoint);
          if(d.status!=S::ok)break;
          const auto age_tail=tail(endpoint,radiation_-radiation_error_,d.a4_e2+d.error_estimate,true,refinement);
          const auto conformal_tail=tail(endpoint,radiation_-radiation_error_,d.a4_e2+d.error_estimate,false,refinement);
          const bool age_ok=!need_age||age_tail.half_width<=policy.absolute_tolerance_gyr/(8*age_scale)+
              policy.relative_tolerance*age_tail.midpoint/8;
          const bool conformal_ok=!need_conformal||conformal_tail.half_width<=policy.absolute_tolerance_mpc/(8*conformal_scale)+
              policy.relative_tolerance*conformal_tail.midpoint/8;
          if(age_ok&&conformal_ok) {
            if(need_age)row.age_tail=age_tail;
            if(need_conformal)row.conformal_tail=conformal_tail;
            earned=true;break;
          }
          endpoint/=2;
          if(!(endpoint>0)) { context.status=S::outside_domain;break; }
        }
        if(!earned)age.status=conformal.status=context.status==S::ok?S::conditioning_budget_exceeded:context.status;
        else {
          auto remaining_integral=[&](bool proper_age,W absolute,W scale) {
            Context interval{background_,policy,source_estimate_,lambda_,today,work,start};
            interval.age=proper_age;interval.squared=true;
            interval.lower=std::sqrt(endpoint);interval.upper=std::sqrt(static_cast<W>(a));
            return integrate(interval,absolute,scale);
          };
          if(need_age) {
            age=remaining_integral(true,policy.absolute_tolerance_gyr,age_scale);
            if(age.status==S::ok) { age.value+=row.age_tail->midpoint;age.error+=row.age_tail->half_width; }
          }
          if(need_conformal) {
            conformal=remaining_integral(false,policy.absolute_tolerance_mpc,conformal_scale);
            if(conformal.status==S::ok) { conformal.value+=row.conformal_tail->midpoint;conformal.error+=row.conformal_tail->half_width; }
          }
        }
      }
      if(need_age)admit(row.outputs[0],age,age_scale,policy.absolute_tolerance_gyr,policy.relative_tolerance);
      if(requested&bit(ThermalClockOutput::comoving_particle_horizon_mpc))
        admit(row.outputs[2],conformal,horizon_scale,policy.absolute_tolerance_mpc,policy.relative_tolerance);
      if(requested&bit(ThermalClockOutput::proper_particle_horizon_mpc))
        admit(row.outputs[3],conformal,horizon_scale*a,policy.absolute_tolerance_mpc,policy.relative_tolerance);
    }
    if(requested&bit(ThermalClockOutput::lookback_gyr)) {
      Integral lookback;
      if(a==1)lookback.status=S::ok;
      else if(radiation_==0&&matter==0&&source.species.empty()&&source_estimate_.radiation_absolute_estimate==0) {
        if(a==0)lookback.status=S::outside_domain;
        else { lookback={S::ok,-std::log(static_cast<W>(a)),0};lookback.error=arithmetic*lookback.value; }
      } else {
        // Separate direct interval context preserves successful finite lookback
        // when an unrelated early-support age/horizon calculation refuses.
        Context interval{background_,policy,source_estimate_,lambda_,today,work,start};
        interval.age=true;
        if(radiation_==0&&matter>0&&radiation_error_==0) {
          interval.dust=true;interval.lower=std::sqrt(static_cast<W>(a));interval.upper=1;
        } else { interval.lower=a;interval.upper=1; }
        lookback=integrate(interval,policy.absolute_tolerance_gyr,age_scale);
      }
      admit(row.outputs[1],lookback,age_scale,policy.absolute_tolerance_gyr,policy.relative_tolerance);
    }
    row.outer_callbacks=work.outer-outer_start;row.momentum_callbacks=work.momentum-momentum_start;
    row.callbacks=work.total()-start;
  }
  result.outer_callbacks=work.outer;result.momentum_callbacks=work.momentum;result.callbacks=work.total();return result;
}
} // namespace irred::cosmology
