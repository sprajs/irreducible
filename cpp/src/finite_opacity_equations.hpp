#pragma once
#include "irred/finite_opacity_source.hpp"
#include "thermal_conformal_epoch.hpp"
#include <cmath>

namespace irred::cosmology::detail::finite_opacity {
using W = long double;
using Core = FiniteOpacityCore;
using Epoch = ThermalConformalEpoch;
inline constexpr unsigned dg=0, ug=1, sg=2, g0=3, g1=4, g2=5;
inline constexpr unsigned db=6, ub=7, dc=8, uc=9, phi=10;
inline W potential_psi(const Epoch &e, W k, const Core &y) {
  return y[phi] - 6*e.hcal*e.hcal*e.fg*y[sg]/(k*k);
}
inline W potential_prime(const Epoch &e, W k, const Core &y, W psi) {
  return -e.hcal*psi + 1.5L*e.hcal*e.hcal/k *
      (e.fc*y[uc]+e.fb*y[ub]+4*e.fg*y[ug]/3);
}
// Independently authored MB equations. All coordinates are dimensionless;
// derivatives carry Mpc^-1. F2=2*sigma, not a second evolved core coordinate.
inline Core core_derivative(const Epoch &e, W k, W opacity,
                            const Core &y, W f3, W p3) {
  Core d; // every entry below is assigned once
  const W psi=potential_psi(e,k,y), pp=potential_prime(e,k,y,psi);
  const W pi=2*y[sg]+y[g0]+y[g2], ratio=4*e.fg/(3*e.fb);
  d[dg]=-4*k*y[ug]/3+4*pp;
  d[ug]=k*(y[dg]/4-y[sg]+psi)+opacity*(y[ub]-y[ug]);
  d[sg]=4*k*y[ug]/15-3*k*f3/10-opacity*y[sg]+opacity*pi/20;
  d[g0]=-k*y[g1]+opacity*(-y[g0]+pi/2);
  d[g1]=k*(y[g0]-2*y[g2])/3-opacity*y[g1];
  d[g2]=k*(2*y[g1]-3*p3)/5+opacity*(-y[g2]+pi/10);
  d[db]=-k*y[ub]+3*pp;
  d[ub]=-e.hcal*y[ub]+k*psi+ratio*opacity*(y[ug]-y[ub]);
  d[dc]=-k*y[uc]+3*pp;
  d[uc]=-e.hcal*y[uc]+k*psi;
  d[phi]=pp;
  return d;
}
inline W hamiltonian(const Epoch &e, W k, const Core &y,
                     W &absolute_term_scale) {
  const W psi=potential_psi(e,k,y), pp=potential_prime(e,k,y,psi);
  const W a=k*k*y[phi], b=3*e.hcal*(pp+e.hcal*psi);
  const W c=1.5L*e.hcal*e.hcal*(e.fc*y[dc]+e.fb*y[db]+e.fg*y[dg]);
  absolute_term_scale=std::abs(a)+std::abs(b)+std::abs(c);
  return a+b+c;
}
inline Core finite_seed(const Epoch &e, W k) {
  Core y{};
  const W f=e.fb+e.fc+4*e.fg/3, x2=e.x2;
  y[phi]=-1/(1.5L+2*x2/(9*f));
  const W v=y[phi]/2, density=-3*v-2*x2*y[phi]/(3*f);
  y[db]=y[dc]=density; y[dg]=4*density/3;
  y[ug]=y[ub]=y[uc]=k*v/e.hcal;
  return y;
}
inline std::array<W,4> raw_source(const Epoch &e, W k, W opacity,
                                W survival, const Core &y) {
  const W psi=potential_psi(e,k,y), pp=potential_prime(e,k,y,psi);
  const W visibility=opacity*survival, pi=2*y[sg]+y[g0]+y[g2];
  return {survival*pp+visibility*y[dg]/4,
          survival*k*psi+visibility*y[ub], visibility*pi/8,
          std::sqrt(6.L)*visibility*pi/8};
}
} // namespace irred::cosmology::detail::finite_opacity
