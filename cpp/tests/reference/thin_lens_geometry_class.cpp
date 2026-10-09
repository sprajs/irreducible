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
  if(ba.N_ncdm>0) {
  ba.Omega0_g=og; ba.Omega0_ur=og*(7.L/8)*powl(r0,4)*nur;
  ba.Omega0_b=.02237L/h2; ba.Omega0_cdm=.12L/h2;
  ba.M_ncdm[0]=.06L/tnu;
  const long double coefficient=2*powl(tnu,4)/(2*pi*pi*critical);
  ba.factor_ncdm[0]=coefficient*ba.H0*ba.H0*4*pi*pi*pi;
  double rho=0,pressure=0;
  background_ncdm_momenta(ba.q_ncdm_bg[0],ba.w_ncdm_bg[0],ba.q_size_ncdm_bg[0],
      ba.N_ncdm?ba.M_ncdm[0]:0,ba.N_ncdm?ba.factor_ncdm[0]:0,0,nullptr,&rho,&pressure,nullptr,nullptr);
  ba.Omega0_ncdm[0]=rho/(ba.H0*ba.H0); ba.Omega0_ncdm_tot=ba.Omega0_ncdm[0];
  ba.Omega0_lambda=1-ba.Omega0_g-ba.Omega0_ur-ba.Omega0_b-ba.Omega0_cdm-ba.Omega0_ncdm_tot;
  } else {
    // Matched curved conserved-fluid realization: explicit source fractions.
    // The selected .1 or -.1 curvature comes from the pinned reference ini.
    ba.Omega0_g=.0001;ba.Omega0_ur=0;ba.Omega0_b=.05;ba.Omega0_cdm=.25;
    ba.Omega0_ncdm_tot=0;
    ba.Omega0_lambda=1-ba.Omega0_g-.3-ba.Omega0_k;
  }
  if(background_init(&pr,&ba)==_FAILURE_) { std::fprintf(stderr,"background:%s\n",ba.error_message); return 2; }
  std::printf("# matched coefficients %.17g %.17g %.17g %.17g %.17g %.17g\n",ba.Omega0_g,ba.Omega0_ur,ba.Omega0_b,ba.Omega0_cdm,ba.N_ncdm?ba.M_ncdm[0]:0,ba.N_ncdm?ba.factor_ncdm[0]:0);
  std::printf("# zl zs Dl_Mpc Ds_Mpc Dls_Mpc beta Sigma_physical_kg_m2 Ddt_Mpc\n");
  std::vector<double> lens(ba.bg_size),source(ba.bg_size);int last=0;
  const double pairs[][2]={{.1,.5},{.5,2},{2,10}};
  for(const auto &pair:pairs) {
    if(background_at_z(&ba,pair[0],long_info,inter_normal,&last,lens.data())==_FAILURE_||
       background_at_z(&ba,pair[1],long_info,inter_normal,&last,source.data())==_FAILURE_)return 3;
    const long double dl=lens[ba.index_bg_ang_distance],ds=source[ba.index_bg_ang_distance];
    const long double radial=source[ba.index_bg_conf_distance]-lens[ba.index_bg_conf_distance];
    // Independent background ODE/interpolation witness; moderate separated
    // epochs only. This CLASS table difference is not the production API.
    long double dls=radial/(1+pair[1]);
    if(ba.K>0)dls=sinl(sqrtl(ba.K)*radial)/sqrtl(ba.K)/(1+pair[1]);
    if(ba.K<0)dls=sinhl(sqrtl(-ba.K)*radial)/sqrtl(-ba.K)/(1+pair[1]);
    const long double sigma=c*c/(4*pi*6.67430e-11L)*ds/(dl*dls)/mpc;
    std::printf("%.17g %.17g %.17Lg %.17Lg %.17Lg %.17Lg %.17Lg %.17Lg\n",pair[0],pair[1],dl,ds,dls,dls/ds,sigma,(1+pair[0])*dl*ds/dls);
  }
  return background_free(&ba)==_FAILURE_;
}
