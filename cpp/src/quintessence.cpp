#include "irred/quintessence.hpp"
#include <algorithm>
#include <array>
#include <cfenv>
#include <cmath>
#include <limits>
namespace irred::cosmology {
namespace {
using W = long double;
using S = numerics::Status;
using State = std::array<W, 6>; // x,y,Omega_m,Omega_r,ln E, dimensionless chi
constexpr W speed = 299792.458L;
bool arithmetic() {
  return std::fegetround() == FE_TONEAREST &&
         std::numeric_limits<W>::digits >= 64;
}
State add(State a, const State &b, W h) {
  for (size_t i = 0; i < a.size(); ++i) a[i] += h * b[i];
  return a;
}
struct Ode {
  W lambda;
  bool distances;
  size_t cap, calls = 0;
  S status = S::ok;
  State rhs(W n, const State &s) {
    if (calls == cap) { status = S::work_limit; return {}; }
    ++calls;
    for (W v : s) if (!std::isfinite(v)) { status = S::overflow; return {}; }
    W q = 3*s[0]*s[0] + 1.5L*s[2] + 2*s[3], k=std::sqrt(1.5L)*lambda;
    return {-3*s[0]+k*s[1]*s[1]+s[0]*q,
            s[1]*(-k*s[0]+q), s[2]*(-3+2*q), s[3]*(-4+2*q), -q,
            distances ? std::exp(-n-s[4]) : 0};
  }
  State rk(W n, const State &s, W h) {
    auto a=rhs(n,s), b=rhs(n+h/2,add(s,a,h/2)),
         c=rhs(n+h/2,add(s,b,h/2)), d=rhs(n+h,add(s,c,h));
    State out=s;
    for(size_t i=0;i<s.size();++i) out[i]+=h*(a[i]+2*b[i]+2*c[i]+d[i])/6;
    return out;
  }
};
void set(QuintessenceValue &out, W value, W error) {
  double v=static_cast<double>(value);
  error+=std::abs(value-W(v));
  double e=static_cast<double>(error);
  if (W(e)<error) e=std::nextafter(e,std::numeric_limits<double>::infinity());
  if (!std::isfinite(v)||!std::isfinite(e)) {out.status=S::overflow;return;}
  out={S::ok,v,e};
}
void refuse(QuintessenceRow &row, S s, bool distances) {
  row.status=s;
  for(auto *v:{&row.e,&row.h_km_s_mpc,&row.w_phi,&row.omega_phi,&row.omega_m,&row.omega_r}) v->status=s;
  if(distances) row.dm_mpc.status=row.dl_mpc.status=s;
}
}
std::optional<size_t> quintessence_payload_bound(size_t n) noexcept {
  constexpr size_t fixed=sizeof(QuintessenceBatch)+sizeof(State)*12+1024;
  if(n>(std::numeric_limits<size_t>::max()-fixed)/sizeof(QuintessenceRow)) return {};
  return fixed+n*sizeof(QuintessenceRow);
}
QuintessenceBackground prepare_quintessence(ExponentialQuintessence s) {
  QuintessenceBackground out;out.source_=s;
  if(!arithmetic()) return out;
  for(double v:{s.lambda,s.h_anchor_km_s_mpc,s.signed_kinetic_fraction_root,s.potential_fraction,s.radiation_fraction})
    if(!std::isfinite(v)){out.status_=S::nonfinite_input;return out;}
  W x=s.signed_kinetic_fraction_root;
  if(std::abs(s.lambda)>20||s.h_anchor_km_s_mpc<=0||s.h_anchor_km_s_mpc>1e6||
     s.potential_fraction<0||s.radiation_fraction<0||
     x*x+W(s.potential_fraction)+W(s.radiation_fraction)>1) {
    out.status_=S::outside_domain;return out;
  }
  out.status_=S::ok;return out;
}
QuintessenceBatch QuintessenceBackground::evaluate(std::span<const double> a,bool distances,QuintessencePolicy p) const {
  QuintessenceBatch out;
  if(status_!=S::ok){out.status=status_;return out;}
  if(!arithmetic()||!std::isfinite(p.absolute_tolerance)||!std::isfinite(p.relative_tolerance)||
     p.absolute_tolerance<=0||p.relative_tolerance<0||p.maximum_halvings>60) return out;
  auto bytes=quintessence_payload_bound(a.size());
  if(a.size()>65536||a.size()>p.maximum_points||!bytes||*bytes>p.maximum_native_bytes||*bytes>size_t(1024)*1024*1024) {
    out.status=S::work_limit;return out;
  }
  out.rows.reserve(a.size());out.status=S::ok;
  for(double scale:a) {
    out.rows.push_back({});auto &r=out.rows.back();r.scale_factor=scale;
    if(!std::isfinite(scale)) {refuse(r,S::nonfinite_input,distances);out.status=r.status;continue;}
    if(scale<std::exp(-8.)||scale>std::exp(4.)) {refuse(r,S::outside_domain,distances);out.status=r.status;continue;}
    W x=source_.signed_kinetic_fraction_root,v=source_.potential_fraction,rad=source_.radiation_fraction;
    State s{x,std::sqrt(v),1-(x*x+v+rad),rad,0,0},errors{};
    W target=std::log(W(scale)),n=0,h=std::copysign(std::min(.05L,std::abs(target)),target);
    Ode ode{source_.lambda,distances&&scale<=1,
      std::min(p.maximum_callbacks_per_point,p.maximum_total_callbacks-out.callbacks)};
    unsigned halvings=0;
    while(n!=target&&ode.status==S::ok) {
      h=std::copysign(std::min(std::abs(h),std::abs(target-n)),target);
      if(n+h==n||n+h/2==n) {ode.status=S::conditioning_budget_exceeded;break;}
      auto full=ode.rk(n,s,h),mid=ode.rk(n,s,h/2),fine=ode.rk(n+h/2,mid,h/2);
      if(ode.status!=S::ok) break;
      State local{};W ratio=0;
      for(size_t i=0;i<s.size();++i) {
        local[i]=std::abs(fine[i]-full[i])/15+
          64*std::numeric_limits<W>::epsilon()*(std::abs(s[i])+std::abs(fine[i])+std::abs(h));
        W allowed=(p.absolute_tolerance+p.relative_tolerance*std::max(std::abs(s[i]),std::abs(fine[i])))*std::abs(h)/std::abs(target);
        ratio=std::max(ratio,local[i]/allowed);
        if(!std::isfinite(fine[i])) ratio=std::numeric_limits<W>::infinity();
      }
      if(fine[1]<0||fine[2]<0||fine[3]<0) ratio=std::numeric_limits<W>::infinity();
      if(ratio>1) {
        ++r.rejected_steps;
        if(++halvings>p.maximum_halvings){ode.status=S::conditioning_budget_exceeded;break;}
        h/=2;continue;
      }
      halvings=0;++r.accepted_steps;
      for(size_t i=0;i<s.size();++i) errors[i]+=local[i];
      s=fine;n=std::abs(target-n)<=std::abs(h)?target:n+h;
      if(ratio<.03L) h*=2;
    }
    r.callbacks=ode.calls;out.callbacks+=ode.calls;
    if(ode.status!=S::ok){refuse(r,ode.status,distances);out.status=r.status;continue;}
    W kinetic=s[0]*s[0],pot=s[1]*s[1],phi=kinetic+pot;
    W ek=2*std::abs(s[0])*errors[0]+errors[0]*errors[0],
      ep=2*std::abs(s[1])*errors[1]+errors[1]*errors[1], ephi=ek+ep;
    W constraint=std::abs(phi+s[2]+s[3]-1);
    r.constraint_residual=static_cast<double>(constraint);
    if(constraint>p.absolute_tolerance+p.relative_tolerance){refuse(r,S::conditioning_budget_exceeded,distances);out.status=r.status;continue;}
    r.status=S::ok;
    W e=std::exp(s[4]),ee=e*std::expm1(errors[4]);
    set(r.e,e,ee);set(r.h_km_s_mpc,e*source_.h_anchor_km_s_mpc,ee*source_.h_anchor_km_s_mpc);
    if(r.h_km_s_mpc.value && *r.h_km_s_mpc.value<=0) {
      r.h_km_s_mpc.value.reset();r.h_km_s_mpc.status=S::outside_domain;
    }
    set(r.omega_phi,phi,ephi);set(r.omega_m,s[2],errors[2]);set(r.omega_r,s[3],errors[3]);
    if(phi>ephi) set(r.w_phi,(kinetic-pot)/phi,2*ephi/(phi-ephi));
    else r.w_phi.status=S::singular;
    if(distances) {
      if(scale>1) r.dm_mpc.status=r.dl_mpc.status=S::outside_domain;
      else {W dm=-speed/source_.h_anchor_km_s_mpc*s[5],edm=speed/source_.h_anchor_km_s_mpc*errors[5];set(r.dm_mpc,dm,edm);set(r.dl_mpc,dm/scale,edm/scale);}
    }
    for(auto *z:{&r.e,&r.h_km_s_mpc,&r.omega_phi,&r.omega_m,&r.omega_r}) if(z->status!=S::ok) {r.status=z->status;out.status=r.status;}
  }
  return out;
}
} // namespace irred::cosmology
