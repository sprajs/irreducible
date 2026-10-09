#include "irred/thin_lens_geometry.hpp"
#include "curved_geometry.hpp"
#include "flat_geometry.hpp"
#include "payload_accounting.hpp"
#include "thermal_constants.hpp"
#include <algorithm>
#include <cfenv>
#include <cmath>
#include <limits>
#include <numbers>
#include <utility>
namespace irred::cosmology {
namespace {
using S=numerics::Status;
using W=long double;
constexpr W rounding=128.L*std::numeric_limits<double>::epsilon();
bool valid(const LensGeometryPolicy &p) {
  return std::fegetround()==FE_TONEAREST&&std::numeric_limits<W>::digits>=64&&
      std::isfinite(p.relative_tolerance)&&p.relative_tolerance>0&&p.relative_tolerance<=1e-2&&
      std::isfinite(p.absolute_tolerance_mpc)&&p.absolute_tolerance_mpc>0&&
      std::isfinite(p.absolute_tolerance_ratio)&&p.absolute_tolerance_ratio>0&&
      std::isfinite(p.absolute_tolerance_kg_m2)&&p.absolute_tolerance_kg_m2>0&&
      p.maximum_depth>0&&p.maximum_depth<=60;
}
std::size_t origins(const ThermalObservables &p) {
  std::size_t n=0;
  if(!irred::detail::checked_payload_add(n,p.source().source_origin.capacity(),1)||
     !irred::detail::checked_payload_add(n,p.source().drag_origin.capacity(),1)||
     !irred::detail::checked_payload_add(n,2,1))return std::numeric_limits<std::size_t>::max();
  return n;
}
struct Distance { S status=S::invalid_input;W value=0,error=0; };
Distance angular(Distance radial,W unit,W k,W redshift) {
  if(radial.status!=S::ok)return radial;
  const W x=radial.value/unit,xe=radial.error/unit+rounding*std::abs(x);
  const W center=unit*detail::curved_transverse(x,k)/(1+redshift);
  const W error=unit*detail::curved_transverse_response(x,xe,k)*xe/(1+redshift)+
      rounding*std::abs(center);
  if(!std::isfinite(center)||!std::isfinite(error))return {S::overflow};
  if(!(center>error)||!detail::physical_representable(center))
    return {S::conditioning_budget_exceeded};
  return {S::ok,center,error};
}
LensGeometryValue admit(Distance d,W absolute,W relative) {
  LensGeometryValue out;out.status=d.status;if(out.status!=S::ok)return out;
  const double stored=static_cast<double>(d.value);
  const W error=d.error+rounding*std::abs(d.value)+std::abs(d.value-stored);
  if(!detail::physical_representable(d.value)||!detail::physical_representable(error)||
      !(d.value>error)) { out.status=S::conditioning_budget_exceeded;return out; }
  if(error>absolute+relative*std::abs(d.value)) {
    out.status=S::conditioning_budget_exceeded;return out;
  }
  out.value=stored;out.error_estimate=static_cast<double>(error);return out;
}
// Positive-factor endpoint propagation preserves common-source dependence:
// no variance addition, independence assumption or denominator cancellation.
Distance quotient(Distance a,Distance b,Distance denominator,W scale) {
  for(auto d:{a,b,denominator})if(d.status!=S::ok)return {d.status};
  if(!(a.value>a.error&&b.value>b.error&&denominator.value>denominator.error))
    return {S::conditioning_budget_exceeded};
  const W value=scale*a.value*b.value/denominator.value;
  const W lower=scale*(a.value-a.error)*(b.value-b.error)/(denominator.value+denominator.error);
  const W upper=scale*(a.value+a.error)*(b.value+b.error)/(denominator.value-denominator.error);
  const W error=std::max(value-lower,upper-value)+rounding*std::abs(value);
  if(!std::isfinite(value)||!std::isfinite(error))return {S::overflow};
  return {S::ok,value,error};
}
Distance product(Distance a,Distance b) {
  return quotient(a,b,{S::ok,1,0},1);
}
} // namespace
ThinLensGeometry::ThinLensGeometry(const ThermalObservables &p,LensGeometryPolicy policy) {
  if(!valid(policy))return;
  if(p.status()!=S::ok) { status_=p.status();return; }
  const auto bytes=thermal_observables_payload_bound(0,p.source().model.species.size(),origins(p));
  if(!bytes||*bytes>policy.maximum_native_bytes||policy.maximum_native_bytes<sizeof(*this)) {
    status_=S::work_limit;return;
  }
  provider_.emplace<ThermalObservables>(p);status_=S::ok;
  const auto actual=evaluation_payload_bound(0);
  if(!actual||*actual>policy.maximum_native_bytes) {
    provider_.emplace<std::monostate>();status_=S::work_limit;
  }
}
ThinLensGeometry::ThinLensGeometry(const CurvedFLRW &p,LensGeometryPolicy policy) {
  if(!valid(policy))return;
  if(p.status()!=S::ok) { status_=p.status();return; }
  if(policy.maximum_native_bytes<sizeof(*this)+sizeof(LensGeometryBatch)+sizeof(numerics::ScalarResult)) {
    status_=S::work_limit;return;
  }
  provider_.emplace<CurvedFLRW>(p);status_=S::ok;
}
ThinLensGeometry::ThinLensGeometry(ThinLensGeometry &&p) noexcept { *this=std::move(p); }
ThinLensGeometry &ThinLensGeometry::operator=(ThinLensGeometry &&p) noexcept {
  if(this!=&p) { provider_=std::move(p.provider_);status_=p.status_;
    p.provider_.emplace<std::monostate>();p.status_=S::invalid_input; }
  return *this;
}
const ThermalObservables *ThinLensGeometry::thermal_provider() const noexcept {
  return std::get_if<ThermalObservables>(&provider_);
}
const CurvedFLRW *ThinLensGeometry::curved_provider() const noexcept {
  return std::get_if<CurvedFLRW>(&provider_);
}
std::string_view ThinLensGeometry::provider_id() const noexcept {
  if(thermal_provider())return thermal_observables_equation_id;
  if(curved_provider())return curved_flrw_id;
  return {};
}
std::optional<std::size_t> ThinLensGeometry::retained_payload_bound() const noexcept {
  irred::detail::PayloadAccounting b(sizeof(*this));
  if(const auto p=thermal_provider()) {
    b.vector(p->source().model.species);b.vector(p->background().source().species);
    b.string(p->source().source_origin);b.string(p->source().drag_origin);
  }
  return b.result();
}
std::optional<std::size_t> ThinLensGeometry::evaluation_payload_bound(std::size_t n) const noexcept {
  const auto retained=retained_payload_bound();if(!retained)return {};
  irred::detail::PayloadAccounting b(*retained+sizeof(LensGeometryBatch));
  b.add(n,2*sizeof(LensGeometryRow));
  if(const auto p=thermal_provider()) {
    // Deliberately conservative: inherited bound includes redundant headers
    // and mapping storage, even though direct evaluation never re-maps.
    const auto provider=thermal_observables_payload_bound(0,p->source().model.species.size(),origins(*p));
    if(!provider)return {};
    b.add(*provider,1);
  } else b.add(1,sizeof(numerics::ScalarResult));
  return b.result();
}
LensGeometryBatch ThinLensGeometry::evaluate(std::span<const LensEpochPair> pairs,
    unsigned outputs,LensGeometryPolicy policy) const {
  LensGeometryBatch out;out.requested_outputs=outputs;
  if(status_!=S::ok) { out.status=status_;return out; }
  if(!valid(policy)||pairs.empty()||outputs==0||outputs>=(1u<<lens_geometry_output_count))return out;
  const auto bytes=evaluation_payload_bound(pairs.size());
  if(pairs.size()>65536||pairs.size()>policy.maximum_pairs||!bytes||*bytes>policy.maximum_native_bytes) {
    out.status=S::work_limit;return out;
  }
  const auto thermal=thermal_provider();const auto curved=curved_provider();
  const W h0=thermal?thermal->source().model.h0_km_s_mpc:curved->specification().h0_km_s_mpc;
  const W k=curved?curved->omega_k():0;
  const W unit=(static_cast<W>(speed_of_light_m_per_s)/1000)/h0;
  out.status=S::ok;out.rows.reserve(pairs.size());
  for(const auto pair:pairs) {
    LensGeometryRow row;row.epochs=pair;row.status=S::ok;
    S failure=S::ok;
    if(!std::isfinite(pair.lens_redshift)||!std::isfinite(pair.source_redshift))failure=S::nonfinite_input;
    else if(!(pair.lens_redshift>0&&pair.source_redshift>pair.lens_redshift))failure=S::outside_domain;
    auto radial=[&](double lower,double upper)->Distance {
      const auto used=row.callbacks;
      if(used>=policy.maximum_callbacks_per_pair||out.callbacks>=policy.maximum_total_callbacks||
          (thermal&&out.momentum_callbacks>=policy.thermal.maximum_total_callbacks))return {S::work_limit};
      const auto remaining=std::min(policy.maximum_callbacks_per_pair-used,
                                   policy.maximum_total_callbacks-out.callbacks);
      if(remaining<3)return {S::work_limit};
      Distance d;
      if(thermal) {
        ThermalObservablePolicy p;
        p.absolute_tolerance_mpc=policy.absolute_tolerance_mpc/64;
        p.relative_tolerance=policy.relative_tolerance/64;
        p.maximum_depth=policy.maximum_depth;p.maximum_callbacks_per_point=remaining;
        p.maximum_total_callbacks=remaining;p.maximum_native_bytes=policy.maximum_native_bytes;
        p.thermal=policy.thermal;
        p.thermal.maximum_total_callbacks=std::min(remaining,
            policy.thermal.maximum_total_callbacks-out.momentum_callbacks);
        const auto r=thermal->radial_distance_between(lower,upper,p);
        row.callbacks+=r.callbacks;row.outer_callbacks+=r.outer_callbacks;
        row.momentum_callbacks+=r.momentum_callbacks;
        out.callbacks+=r.callbacks;out.outer_callbacks+=r.outer_callbacks;
        out.momentum_callbacks+=r.momentum_callbacks;
        d.status=r.distance_mpc.status;
        if(r.distance_mpc.value) { d.value=*r.distance_mpc.value;d.error=r.distance_mpc.error_estimate; }
      } else {
        CurvedFLRWPolicy p;
        p.absolute_tolerance_mpc=policy.absolute_tolerance_mpc/64;
        p.relative_tolerance=policy.relative_tolerance/64;p.maximum_depth=policy.maximum_depth;
        p.maximum_callbacks_per_point=remaining;p.maximum_total_callbacks=remaining;
        p.maximum_native_bytes=policy.maximum_native_bytes;
        const auto r=curved->radial_distance_between(lower,upper,p);
        row.callbacks+=r.evaluations;row.outer_callbacks+=r.evaluations;
        out.callbacks+=r.evaluations;out.outer_callbacks+=r.evaluations;
        d={r.status,r.value,r.error_estimate};
      }
      return d;
    };
    std::array<Distance,3> distances{};
    if(failure==S::ok) {
      const auto rl=radial(0,pair.lens_redshift);
      if(rl.status!=S::ok)failure=rl.status;
      Distance rs,interval;
      if(failure==S::ok) { rs=radial(0,pair.source_redshift);if(rs.status!=S::ok)failure=rs.status; }
      auto before_conjugate=[&](Distance r) {
        if(k>=0)return true;
        const W phase=std::sqrt(-k)*(r.value+r.error)/unit;
        return phase+rounding*std::abs(phase)<std::numbers::pi_v<W>;
      };
      if(failure==S::ok&&(!before_conjugate(rl)||!before_conjugate(rs)))failure=S::outside_domain;
      if(failure==S::ok) { interval=radial(pair.lens_redshift,pair.source_redshift);
        if(interval.status!=S::ok)failure=interval.status;
        else if(!before_conjugate(interval))failure=S::outside_domain; }
      if(failure==S::ok)distances={angular(rl,unit,k,pair.lens_redshift),
          angular(rs,unit,k,pair.source_redshift),angular(interval,unit,k,pair.source_redshift)};
    }
    const auto dl=distances[0],ds=distances[1],dls=distances[2];
    for(unsigned n=0;n<lens_geometry_output_count;++n)if(outputs&(1u<<n)) {
      auto &value=row.outputs[n];
      if(failure!=S::ok)value.status=failure;
      else {
        Distance d;W absolute=policy.absolute_tolerance_mpc;
        if(n<3)d=distances[n];
        else if(n==3) { d=quotient(dls,{S::ok,1,0},ds,1);absolute=policy.absolute_tolerance_ratio; }
        else if(n==4) {
          const W c=speed_of_light_m_per_s;
          const W sigma_scale=c*c/(4*std::numbers::pi_v<W>*detail::thermal_gravitational_constant*
                                    megaparsec_in_metres_wide());
          d=quotient(ds,{S::ok,1,0},product(dl,dls),sigma_scale);
          absolute=policy.absolute_tolerance_kg_m2;
        } else d=quotient(dl,ds,dls,1+static_cast<W>(pair.lens_redshift));
        value=admit(d,absolute,policy.relative_tolerance);
      }
      if(value.status!=S::ok&&row.status==S::ok)row.status=value.status;
    }
    if(row.status!=S::ok&&out.status==S::ok)out.status=row.status;
    out.rows.push_back(std::move(row));
  }
  return out;
}
} // namespace irred::cosmology
