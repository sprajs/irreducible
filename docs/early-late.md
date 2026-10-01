# One early and late expansion identity

The native C++ `irred/early_late.hpp` API uses the same `EarlyFlatModel` and
supplied-drag `SoundHorizonRequest` as the sound-horizon calculation. It supplies
conditional distances and BAO ratio predictions only. It supplies no BAO data
likelihood, fit, drag prediction, thermal history, growth or CMB qualification.

The shared compiled state owns the safe nonnegative flat closure and
P(a)=Omega_r+Omega_m*a+Omega_Lambda*a^4, with E(a)=sqrt(P(a))/a^2.
Its domain is exactly the sound model domain: finite H0>0 in km/s/Mpc,
Omega_r>0, Omega_gamma>0, 0<=Omega_b<=Omega_m,
Omega_gamma<=Omega_r, and Omega_Lambda=1-Omega_m-Omega_r>=0.
All queries have finite z>=0. The request retains finite z_drag>=0 and a
nonempty drag_origin even when only distances are requested. No unrelated
radiation-free model or free H0*r_d parameter enters this calculation.

For the declared flat FLRW, same-redshift, photon-conserving geometry,
H=H0*E, D_H=c/H, D_M=(c/H0)*integral_0^z dz'/E(z'),
D_L=(1+z)*D_M and D_V=cbrt(D_M^2*z*D_H). The speed of light is
exactly 299792.458 km/s and all distances and the ruler are in Mpc.
Production quadrature uses x=log(1+z), t=x/x_max, with integrand
exp(-x)/sqrt(P(exp(-x))) on t in [0,1]. This resolves tiny redshifts
without subtracting nearly equal scale factors. There is no high-z cutoff.
At z=0, D_M,D_L,D_V and their ruler ratios are exactly zero; E=1 and
D_H=c/H0. A requested positive output that cannot be represented as a
normal binary64 value fails honestly.

A mask requests only needed outputs. One batch owns its model/drag source,
redshifts, per-output status and optional scalar/empirical error diagnostic.
Only requested distance integrals run; the ruler is evaluated once per batch
only if a ratio is requested. Distance success survives a failed ruler.
D_M/r_s,D_H/r_s,D_V/r_s use that same model and supplied drag identity.
Their diagnostics propagate the ruler and numerator intervals additively;
no independent-random-error assumption is made. These are numerical
admission diagnostics, not rigorous bounds or observational qualification.

Policies separately declare Mpc and dimensionless ratio tolerances, distance
callback/depth limits, total callbacks including sound, point count and a
conservative simultaneous owned-payload allowance. Admission precedes origin
copies, row allocation and redshift scans. The bounded sound policy also
applies. Payload excludes allocator metadata, recursion stack and RSS.
Unrequested entries remain absent with invalid_input status. Batch success
means processing was admitted; every requested scalar retains its own status.
Arithmetic requires FE_TONEAREST and long double with at least 64 mantissa
bits and exponent range 16384. No automatic jitter, dropped query or weakened
budget is used.

The frozen comparison allocation is 1e-9 Mpc+2e-11*|reference| for distances,
2e-11 relative for E,H,D_H, and 1e-11+5e-11*|reference| for ratios.
Independent-reference refinement must occupy at most 5% of each allocation.
Analytic radiation-only and Lambda=0 matter+radiation controls, H0 scaling,
zero/tiny/high redshift, invalid closure/subsets, quotas and existing sound
regressions are required. Independent direct-z quadrature has different
coordinates from production; shared polynomial helpers are not independent
evidence. Named tested controls do not qualify every arbitrary request.

The native [conditional BAO consumer](bao-conditional.md) composes these ratios with a retained full ordered Gaussian observation object and checks downstream projection sensitivity separately. The supplied drag origin and shared parameter/unit identity remain explicit; physical validity of a released compression is a separate condition.
