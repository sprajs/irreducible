#include "irred/curved_flrw.hpp"
#include "curved_geometry.hpp"
#include <algorithm>
#include <cfenv>
#include <cmath>
#include <limits>
#include <numbers>
#include <utility>
namespace irred::cosmology {
namespace {
using S=numerics::Status;
constexpr long double c=299792.458L;
bool arithmetic(){return std::fegetround()==FE_TONEAREST && std::numeric_limits<long double>::digits>=64;}
struct Context {
  CurvedFLRWSpec s; long double k;
  mutable long double dependency=0,cast_fraction=0,coordinate_lower=0,coordinate_width=1;
  mutable S status=S::ok;
};
long double square(const Context& p,long double u){return ((p.s.omega_r*u+p.s.omega_m)*u+p.k)*u*u+p.s.omega_lambda;}
long double scale(const Context& p,long double u){return ((p.s.omega_r*u+p.s.omega_m)*u+std::abs(p.k))*u*u+p.s.omega_lambda;}
bool positive_path(const Context& p,long double end){
  auto good=[&](long double u){return square(p,u)>256*std::numeric_limits<double>::epsilon()*scale(p,u);};
  if(!good(1)||!good(end)) return false;
  if(p.k<0){
    long double root=0;
    if(p.s.omega_r>0) root=-4*p.k/(3*p.s.omega_m+std::sqrt(9.L*p.s.omega_m*p.s.omega_m-32*p.s.omega_r*p.k));
    else if(p.s.omega_m>0) root=-2*p.k/(3*p.s.omega_m);
    if(root>1 && root<end && !good(root)) return false;
  }
  return true;
}
double inverse_e(double z,const void* raw) {
  const auto &p=*static_cast<const Context*>(raw);
  if(p.status!=S::ok)return std::numeric_limits<double>::quiet_NaN();
  const long double physical_z=p.coordinate_lower+p.coordinate_width*z;
  const long double q=square(p,1.L+physical_z);
  const long double error=128*std::numeric_limits<double>::epsilon()*scale(p,1.L+physical_z);
  if(!(q>error)) {
    p.status=S::conditioning_budget_exceeded;return std::numeric_limits<double>::quiet_NaN();
  }
  const long double lower=q-error,value=1/std::sqrt(q);
  p.dependency=std::max(p.dependency,error/(std::sqrt(lower)*(std::sqrt(q)+std::sqrt(lower))));
  const double stored=static_cast<double>(value);
  if(!std::isfinite(stored)||!std::isnormal(stored)) {
    p.status=S::outside_domain;return std::numeric_limits<double>::quiet_NaN();
  }
  p.cast_fraction=std::max(p.cast_fraction,std::abs(value-stored)/value);
  return stored;
}
numerics::ScalarResult integrate_radial(const Context &ctx,double lower,double upper,
                                       const CurvedFLRWPolicy &policy,std::size_t remaining) {
  const long double unit=c/ctx.s.h0_km_s_mpc;
  ctx.status=S::ok;ctx.dependency=0;ctx.cast_fraction=0;
  // Keep the historical observer z mesh. A two-epoch interval instead uses
  // unit support and a wide physical coordinate, including adjacent doubles.
  // No large-distance subtraction or collapsed physical subdivision occurs.
  const bool interval=lower>0;
  ctx.coordinate_lower=interval?lower:0;
  ctx.coordinate_width=interval?static_cast<long double>(upper)-lower:1;
  auto result=numerics::integrate(inverse_e,&ctx,interval?0:lower,interval?1:upper,
      {static_cast<double>(std::min(static_cast<long double>(policy.absolute_tolerance_mpc)/unit/32/
          ctx.coordinate_width,static_cast<long double>(std::numeric_limits<double>::max()))),
       policy.relative_tolerance/32,std::min(remaining,policy.maximum_callbacks_per_point),
       policy.maximum_depth});
  if(ctx.status!=S::ok)result.status=ctx.status;
  if(result.status==S::ok) {
    result.value=static_cast<double>(static_cast<long double>(result.value)*ctx.coordinate_width);
    result.error_estimate=static_cast<double>(static_cast<long double>(result.error_estimate)*ctx.coordinate_width);
    const long double response=(std::abs(static_cast<long double>(result.value))+
        result.error_estimate)/(1-ctx.cast_fraction);
    result.error_estimate=static_cast<double>(result.error_estimate+
        response*(ctx.dependency+ctx.cast_fraction));
  }
  return result;
}


}
CurvedFLRW::CurvedFLRW(CurvedFLRWSpec s):spec_(s){
  if(!arithmetic())return;
  if(!std::isfinite(s.h0_km_s_mpc)||!std::isfinite(s.omega_m)||!std::isfinite(s.omega_r)||!std::isfinite(s.omega_lambda)){status_=S::nonfinite_input;return;}
  if(s.h0_km_s_mpc<1 || s.h0_km_s_mpc>1000 || s.omega_m<0 || s.omega_r<0 || s.omega_lambda<0 || s.omega_m>10 || s.omega_r>10 || s.omega_lambda>10){status_=S::outside_domain;return;}
  omega_k_=1-s.omega_m-s.omega_r-s.omega_lambda;status_=S::ok;
}
CurvedFLRW::CurvedFLRW(CurvedFLRW&& o) noexcept{*this=std::move(o);}
CurvedFLRW& CurvedFLRW::operator=(CurvedFLRW&& o) noexcept{if(this!=&o){spec_=o.spec_;omega_k_=o.omega_k_;status_=o.status_;o.status_=S::invalid_input;}return *this;}
CurvedFLRWBatch CurvedFLRW::evaluate(std::span<const double> zs,CurvedFLRWPolicy policy)const{
  CurvedFLRWBatch out;
  if(status_!=S::ok){out.status=status_;return out;}
  if(!arithmetic()||zs.empty()||!std::isfinite(policy.relative_tolerance)||!std::isfinite(policy.absolute_tolerance_mpc)||policy.relative_tolerance<=0||policy.relative_tolerance>1e-2||policy.absolute_tolerance_mpc<=0||policy.maximum_depth==0||policy.maximum_depth>60||policy.maximum_callbacks_per_point<3||policy.maximum_total_callbacks<3)return out;
  if(zs.size()>65536||zs.size()>policy.maximum_points||policy.maximum_native_bytes<sizeof(out)||zs.size()>(policy.maximum_native_bytes-sizeof(out))/sizeof(CurvedFLRWRow)){out.status=S::work_limit;return out;}
  out.rows.reserve(zs.size());out.status=S::ok;
  Context ctx{spec_,omega_k_};const long double unit=c/spec_.h0_km_s_mpc;
  for(double z:zs){
    CurvedFLRWRow row;row.redshift=z;
    if(!std::isfinite(z)){row.status=S::nonfinite_input;}
    else if(z<0||z>10000||!positive_path(ctx,1.L+z)){row.status=S::outside_domain;}
    else{
      const auto remaining=policy.maximum_total_callbacks-out.callbacks;
      if(z>0 && remaining<3){row.status=S::work_limit;}
      else{
        numerics::ScalarResult integral{S::ok,0,0,0};
        if(z>0) integral=integrate_radial(ctx,0,z,policy,remaining);
        row.callbacks=integral.evaluations;out.callbacks+=row.callbacks;row.status=integral.status;
        if(row.status==S::ok){
          const long double esquare=square(ctx,1.L+z), evalue=std::sqrt(esquare);
          const long double eerror=128*std::numeric_limits<double>::epsilon()*scale(ctx,1.L+z)/evalue;
          const long double chi=integral.value,chierr=integral.error_estimate+64*std::numeric_limits<double>::epsilon()*std::abs(chi);
          const long double q=std::sqrt(std::abs(ctx.k));
          if(ctx.k<0 && q*(chi+chierr)>=std::numbers::pi_v<long double>){row.status=S::outside_domain;}
          else{
            const long double dm=unit*detail::curved_transverse(chi,ctx.k);
            const long double jac=detail::curved_transverse_response(chi,chierr,ctx.k);
            const long double err=unit*jac*chierr+128*std::numeric_limits<double>::epsilon()*std::abs(dm);
            const long double radial=unit*chi, radialerr=unit*chierr;
            const double stored_radial=static_cast<double>(radial),stored_dm=static_cast<double>(dm),stored_da=static_cast<double>(dm/(1.L+z)),stored_dl=static_cast<double>(dm*(1.L+z));
            const double stored_e=static_cast<double>(evalue),stored_h=spec_.h0_km_s_mpc*stored_e;
            if(!std::isfinite(dm)||!std::isfinite(err)||!std::isfinite(stored_radial)||!std::isfinite(stored_dm)||!std::isfinite(stored_da)||!std::isfinite(stored_dl)||!std::isfinite(stored_e)||!std::isfinite(stored_h)){row.status=S::overflow;}
            else if(z>0 && (stored_radial<=0||stored_dm<=0||stored_da<=0||stored_dl<=0)){row.status=S::outside_domain;}
            else if(eerror>policy.relative_tolerance*evalue || err>policy.absolute_tolerance_mpc+policy.relative_tolerance*std::abs(dm)||radialerr>policy.absolute_tolerance_mpc+policy.relative_tolerance*std::abs(radial)){row.status=S::conditioning_budget_exceeded;}
            else{
              row.e=static_cast<double>(evalue);row.e_error_estimate=static_cast<double>(eerror);row.h_error_estimate_km_s_mpc=static_cast<double>(spec_.h0_km_s_mpc*eerror);row.h_km_s_mpc=spec_.h0_km_s_mpc*row.e;
              row.radial_mpc=static_cast<double>(radial);row.transverse_mpc=static_cast<double>(dm);row.angular_mpc=static_cast<double>(dm/(1.L+z));row.luminosity_mpc=static_cast<double>(dm*(1.L+z));
              row.radial_error_estimate_mpc=static_cast<double>(radialerr);row.transverse_error_estimate_mpc=static_cast<double>(err);
            }
          }
        }
      }
    }
    if(row.status!=S::ok && out.status==S::ok)out.status=row.status;
    out.rows.push_back(row);
  }
  return out;
}
numerics::ScalarResult CurvedFLRW::radial_distance_between(double lower,double upper,
    CurvedFLRWPolicy policy) const {
  numerics::ScalarResult out;
  if(status_!=S::ok) { out.status=status_;return out; }
  if(!arithmetic()||!std::isfinite(policy.relative_tolerance)||
      !std::isfinite(policy.absolute_tolerance_mpc)||policy.relative_tolerance<=0||
      policy.relative_tolerance>1e-2||policy.absolute_tolerance_mpc<=0||
      policy.maximum_depth==0||policy.maximum_depth>60||
      policy.maximum_callbacks_per_point<3||policy.maximum_total_callbacks<3)return out;
  if(!std::isfinite(lower)||!std::isfinite(upper)) { out.status=S::nonfinite_input;return out; }
  Context ctx{spec_,omega_k_};
  if(lower<0||upper<lower||upper>10000||!positive_path(ctx,1.L+upper)) {
    out.status=S::outside_domain;return out;
  }
  if(policy.maximum_points<1||policy.maximum_native_bytes<sizeof(out)) { out.status=S::work_limit;return out; }
  if(lower==upper) { out.status=S::ok;return out; }
  const auto result=integrate_radial(ctx,lower,upper,policy,policy.maximum_total_callbacks);
  out.status=result.status;out.evaluations=result.evaluations;if(out.status!=S::ok)return out;
  const long double unit=c/spec_.h0_km_s_mpc,center=unit*result.value;
  const double stored=static_cast<double>(center);
  const long double error=unit*result.error_estimate+
      128*std::numeric_limits<double>::epsilon()*std::abs(center)+std::abs(center-stored);
  if(!std::isfinite(stored)||!(stored>0)||!std::isnormal(stored)||!std::isfinite(error)) {
    out.status=S::outside_domain;return out;
  }
  if(error>policy.absolute_tolerance_mpc+policy.relative_tolerance*center) {
    out.status=S::conditioning_budget_exceeded;return out;
  }
  out.value=stored;out.error_estimate=static_cast<double>(error);return out;
}

}
