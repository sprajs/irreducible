# Same-state cosmic clocks and particle horizon

Implementation contract for the stacked native clock package, based on
verified Irreducible main `79e71ab9edee2b15c9d01830d2c86434748eabce`
(PR #72 head `5ef0a1627983f1b15199d768bc45b4328702e5ff`). The compiled
interface and tests in this branch are preparation work; no build, numerical
qualification or completed consumer is claimed before its allocated validation.

A retained `ThermalBackground` supplies the same photon, dust, explicitly
massless radiation, zero-chemical-potential FD relic and nonnegative flat Lambda
state as the existing stress and observable owners. The source equations and
constant conventions remain in [thermal neutrino](thermal-neutrino.md) and
[effective thermal neutrino](effective-neutrino.md). This package adds proper
age, direct lookback and comoving/proper particle horizon. No BBN, recombination,
changing early relativistic degrees of freedom, perturbations, event horizon,
observed stellar age or full standard cosmology follows.

## Equations, time convention and domain

For D(a)=a^4 E(a)^2 and conformal time d eta=dt/a,

    H0*t(a) = integral_0^a a'/sqrt(D(a')) da'
    H0*eta(a) = integral_0^a 1/sqrt(D(a')) da'
    H0*lookback(a) = integral_a^1 a'/sqrt(D(a')) da'.

Age and lookback are Gyr. One Julian year is exactly 365.25*86400 seconds;
one Gyr is 1e9 such years. The H0 inverse uses the existing exact SI/IAU Mpc
conversion. Comoving particle horizon is c*eta in Mpc, with
c=299792.458 km/s; proper particle horizon is a*c*eta. Neither is a future
event horizon. Lookback is evaluated directly on [a,1], rather than subtracting
large ages near a=1.

The declared model is extrapolated formally to a=0 with its existing constant
photon/FD species content. This does not reproduce electron/positron annihilation,
inflation or a changing early thermal history. Finite 0<=a<=1 is admitted,
with each requested output retaining its own status. Radiation and matter
models have exact age/horizon zero at a=0, and lookback zero at a=1. A pure
Lambda state has divergent age and particle horizon from a=0 and refuses
finite big-bang outputs, including a fabricated zero at its nominal endpoint.
Its direct lookback at positive a is a finite-interval control. Invalid and
failed rows remain in supplied order.

## Early support and numerical estimates

For D0=D(0)>0, nonnegative dust/Lambda and the increasing FD density moment
make D(a) nondecreasing. For 0<=a<=ae, the effective source therefore gives

    ae^2/[2sqrt(D(ae))] <= H0*t(ae) <= ae^2/[2sqrt(D0)]
    ae/sqrt(D(ae)) <= H0*eta(ae) <= ae/sqrt(D0).

Endpoint numerical estimates enter as D(ae)+eD(ae) for the lower endpoint and
D0-eD0 for the upper endpoint; require positive denominator intervals. Retain
the midpoint and half-width of each early contribution. Halve ae until each
requested tail consumes at most one eighth of its own allocation, retaining
all endpoint/refinement work and failures. Integrate the remaining interval;
no early support is dropped or renormalized. The estimates inherited from D
are empirical, so these numerical brackets are not continuous certified bounds.

A radiation-free source with no relics and positive dust uses u=sqrt(a) to
remove the integrable endpoint singularity: age integrand
2*u^2/sqrt(Omega_m+Omega_Lambda*u^6), and conformal integrand
2/sqrt(Omega_m+Omega_Lambda*u^6). Pure Lambda is handled as the distinct
zero-dust divergent case. The source's retained Lambda coefficient and its
normalization diagnostic remain the owner; no new closure is fitted.

Outer quadrature, early half-width, sampled inverse-square-root dependency,
reporting casts and arithmetic estimates are added. Optional deterministic
source-conversion estimates use the effective-neutrino diagnostic convention,
including the present-normalization direction into Lambda; they are not physical
input uncertainty. Positive unrepresentable values/diagnostics refuse instead
of disappearing. Default clock admission is **1e-8 Gyr + 2e-10 relative**;
horizon admission is **1e-8 Mpc + 2e-10 relative**. Independent reference
refinement/uncertainty must occupy <=5% of these frozen allocations. Near-zero
lookback uses its absolute allocation. Requested outputs alone incur their
integrals. Proper/comoving horizon share one conformal integral.

Policy bounds points, whole-batch/per-point combined outer and FD work,
whole-batch momentum work, outer/momentum depth, tail refinements and owned
payload before allocation. All failed endpoints, callbacks and direct fallback
work count against the original limits; no implicit larger cap is installed.
Copies own independently retained state, moves invalidate the old owner, and
method/rounding admission remains consistent with the retained thermal owner.

## Planned independent controls and acceptance

Pure radiation, Einstein--de Sitter and flat dust+Lambda proper age provide
analytic controls. The latter uses

    H0*t(a)=2/(3sqrt(OmegaLambda))
             *asinh(sqrt(OmegaLambda/Omega_m)*a^(3/2)),

with the exact zero-Lambda limit handled independently. Tests will retain
endpoint/divergence distinctions, signed input/resource/copy/move refusals,
requested-output omission, near-one direct lookback and tail/callback accounting.
Separate high-precision/refined FD clock integrals and matched-constant CLASS
age/conformal-time witnesses are required before qualification. Differential
clock identities use a separate finite-difference allocation and refinement;
they do not replace the independent clock references. The relation
chi(a)=c*(eta(1)-eta(a)) against existing distance outputs is a valuable
shared-state consistency check, not independent source validation.

Reference execution was approved by the focused Astra contract review before
any clock build or reference run. Independent FD clock controls will use
fixed high-precision Gauss--Legendre momentum quadrature on [0,80], explicit
exponential-tail estimates, and a full-support squared-coordinate outer
Gauss--Legendre rule at independently doubled orders/precision. This differs
from the production adaptive rule and retained early brackets. Matched CLASS
controls retain the earlier source/coefficient match and separately record
CLASS finite-start age/conformal corrections. No unmodified CLASS source or
observational qualification is inherited.

Five-point log-a derivative diagnostics have separate empirical allocations:
**2e-6 Gyr + 2e-7 relative** for proper-age derivatives and
**2e-4 Mpc + 2e-7 relative** for comoving-horizon derivatives. Step-halving
refinement must occupy <=5% of that derivative allocation. These diagnostics
are consistency checks, not the independent 1e-8 clock/horizon references.
