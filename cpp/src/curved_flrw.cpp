#include "irred/curved_flrw.hpp"
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
struct Context { CurvedFLRWSpec s; long double k; };
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
double inverse_e(double z,const void* raw){const auto& p=*static_cast<const Context*>(raw); auto q=square(p,1.L+z);return q>0?static_cast<double>(1/std::sqrt(q)):std::numeric_limits<double>::quiet_NaN();}
// Entire even series avoids division by sqrt(|Ok|) and near-flat cancellation.
long double transverse(long double chi,long double k){
  const long double t=k*chi*chi;
  if(std::abs(t)<1e-4L) return chi*(1+t*(1.L/6+t*(1.L/120+t*(1.L/5040+t/362880))));
  const long double q=std::sqrt(std::abs(k));
  return k>0?std::sinh(q*chi)/q:std::sin(q*chi)/q;
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
        if(z>0) integral=numerics::integrate(inverse_e,&ctx,0,z,{static_cast<double>(policy.absolute_tolerance_mpc/unit/32),policy.relative_tolerance/32,std::min(remaining,policy.maximum_callbacks_per_point),policy.maximum_depth});
        row.callbacks=integral.evaluations;out.callbacks+=row.callbacks;row.status=integral.status;
        if(row.status==S::ok){
          const long double esquare=square(ctx,1.L+z), evalue=std::sqrt(esquare);
          const long double eerror=128*std::numeric_limits<double>::epsilon()*scale(ctx,1.L+z)/evalue;
          const long double chi=integral.value,chierr=integral.error_estimate+64*std::numeric_limits<double>::epsilon()*std::abs(chi);
          const long double q=std::sqrt(std::abs(ctx.k));
          if(ctx.k<0 && q*(chi+chierr)>=std::numbers::pi_v<long double>){row.status=S::outside_domain;}
          else{
            const long double dm=unit*transverse(chi,ctx.k);
            const long double jac=ctx.k>0?std::cosh(q*(chi+chierr)):1;
            const long double err=unit*jac*chierr+128*std::numeric_limits<double>::epsilon()*std::abs(dm);
            const long double radial=unit*chi, radialerr=unit*chierr;
            if(!std::isfinite(dm)||!std::isfinite(err)||!std::isfinite(dm*(1.L+z))){row.status=S::overflow;}
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
}
