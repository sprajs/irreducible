// Original external reference driver; no copied external implementation.
// Compile only against pinned CLASS runtime objects. Explicitly match the
// native constants through background physical coefficients, leaving runtime
// source/objects unchanged. This is a matched-FD numerical comparison, not
// unmodified Planck source equality.
extern "C" {
#include "class.h"
}
#include <cmath>
#include <cstdio>
#include <vector>
int main(int argc,char **argv) {
  precision pr{}; background ba{}; thermodynamics th{}; perturbations pt{};
  primordial pm{}; fourier fo{}; transfer tr{}; harmonic hr{}; lensing le{};
  distortions sd{}; output op{}; ErrorMsg err{};
  if(input_init(argc,argv,&pr,&ba,&th,&pt,&tr,&pm,&hr,&fo,&le,&sd,&op,err)==_FAILURE_) {
    std::fprintf(stderr,"input:%s\n",err); return 1;
  }
  const long double pi=acosl(-1),c=299792458.L,ev=1.602176634e-19L,kb=1.380649e-23L/ev;
  const long double hbar=6.62607015e-34L/(2*pi),mpc=648000000000.L*149597870700.L/pi;
  const long double H=67.4L,h2=(H/100)*(H/100),T=2.7255L,tnu=T*.71611L*kb;
  const long double hs=H*1000/mpc,critical=3*hs*hs*c*c/(8*pi*6.67430e-11L)/ev*powl(hbar*c/ev,3);
  const long double og=pi*pi/15*powl(kb*T,4)/critical,r0=cbrtl(4.L/11),nur=3.046L-powl(.71611L/r0,4);
  ba.Omega0_g=og; ba.Omega0_ur=og*(7.L/8)*powl(r0,4)*nur;
  ba.Omega0_b=.02237L/h2; ba.Omega0_cdm=.12L/h2;
  ba.M_ncdm[0]=.06L/tnu;
  const long double coefficient=2*powl(tnu,4)/(2*pi*pi*critical);
  ba.factor_ncdm[0]=coefficient*ba.H0*ba.H0*4*pi*pi*pi;
  double rho=0,pressure=0;
  background_ncdm_momenta(ba.q_ncdm_bg[0],ba.w_ncdm_bg[0],ba.q_size_ncdm_bg[0],
      ba.M_ncdm[0],ba.factor_ncdm[0],0,nullptr,&rho,&pressure,nullptr,nullptr);
  ba.Omega0_ncdm[0]=rho/(ba.H0*ba.H0); ba.Omega0_ncdm_tot=ba.Omega0_ncdm[0];
  ba.Omega0_lambda=1-ba.Omega0_g-ba.Omega0_ur-ba.Omega0_b-ba.Omega0_cdm-ba.Omega0_ncdm_tot;
  if(background_init(&pr,&ba)==_FAILURE_) { std::fprintf(stderr,"background:%s\n",ba.error_message); return 2; }
  std::printf("# matched coefficients %.17g %.17g %.17g %.17g %.17g %.17g\n",ba.Omega0_g,ba.Omega0_ur,ba.Omega0_b,ba.Omega0_cdm,ba.M_ncdm[0],ba.factor_ncdm[0]);
  std::printf("# a age_gyr lookback_gyr comoving_horizon_mpc proper_horizon_mpc\n");
  // CLASS stores proper time in Mpc/c. Convert with this project's exact SI
  // c/IAU Mpc and declared exact Julian year, avoiding CLASS Gyr convention.
  const long double seconds_per_mpc=mpc/c;
  const long double gyr=365.25L*86400*1000000000;
  const long double today_time=ba.background_table[(ba.bt_size-1)*ba.bg_size+ba.index_bg_time];
  std::printf("# source a_initial %.17g CLASS_RD_initial_age H0*t=1/(2H), eta=1/(aH)\n",exp(ba.loga_table[0]));
  std::vector<double> values(ba.bg_size); int last=0;
  for(double a:{.001,.01,.1,.5,1.}) {
    if(background_at_z(&ba,1/a-1,long_info,inter_normal,&last,values.data())==_FAILURE_) {
      std::fprintf(stderr,"at_z:%s\n",ba.error_message);return 3;
    }
    const long double proper_time=values[ba.index_bg_time];
    const long double horizon=ba.conformal_age-values[ba.index_bg_conf_distance];
    std::printf("%.17g %.17Lg %.17Lg %.17Lg %.17Lg\n",a,
        proper_time*seconds_per_mpc/gyr,(today_time-proper_time)*seconds_per_mpc/gyr,horizon,a*horizon);
  }
  return background_free(&ba)==_FAILURE_;
}
