# Retained two-epoch thin-lens geometry

Implementation contract frozen before reference execution, 2026-10-09.
The compiled public interface is `irred/thin_lens_geometry.hpp`; source review and
numerical acceptance remain distinct.

Existing owners: CurvedFLRW supplies observer radial/transverse/DA/DL, but its
guide excludes two-epoch lens-source distance. Its private transverse function
is the sole existing sin/sinh/near-flat-series law. ThermalObservables retains
an explicit flat photon/dust/FD/Lambda state and distances; it supplies no
lens prefactors. SISThinLens already has flat massless EarlyLate distances and
synthetic SIS/sheet/Fermat propagation. NFW requires supplied Sigma_crit and
D_l; Plummer lens uses supplied physical kg/m² critical density. Late flat
volume and sampled/temporal spectral flux already exist. Their relabeling is
outside this new package.

Choose one new compiled consumer: same-observer lens geometry from a typed
retained ThermalObservables OR CurvedFLRW provider. These retain separate
physical identities. No thermal+curvature model is created. Factories accept
these actual immutable provider types, not arbitrary distance values/string
identity; a bounded owned provider variant is a compiled two-owner choice,
not a runtime plugin or recipe. No background re-fit/re-normalization or drag
prediction. Thermal provider drag provenance remains present but unqueried.

For finite 0<zl<zs within the selected provider's redshift/path domain:

  X(z)=integral_0^z dz/E(z), DH=c/H0, K=Omega_k retained by provider
  DM_l=DH*S_K(X_l), DM_s=DH*S_K(X_s)
  DeltaX=integral_zl^zs dz/E(z)  (DIRECT, not large-distance subtraction)
  Dl=DM_l/(1+zl), Ds=DM_s/(1+zs), Dls=DH*S_K(DeltaX)/(1+zs)
  beta=Dls/Ds
  Sigma_crit=c²/(4*pi*G)*Ds_Mpc/(Dl_Mpc*Dls_Mpc)/Mpc_metres
                                                PHYSICAL kg/m²
  D_time_delay=(1+zl)*Dl*Ds/Dls             Mpc

S_K(x)=sinh(sqrt(K)*x)/sqrt(K) for K>0, x for K=0, and
sin(sqrt(-K)*x)/sqrt(-K) for K<0. Extract/factor the existing CurvedFLRW
transverse law into ONE shared owner; both old and new callers use it. Factor
existing retained radial integration into bounded direct interval helpers;
no duplicated E² or FD/photon evolution equations. Thermal has K=0. The
curved owner must expose its actual retained K and complete positive-path
admission, rather than recomputing closure with different rounding.

Sigma uses exact SI c and existing exact SI/IAU Mpc. Adopt the SAME existing
thermal G=6.67430e-11 m³ kg^-1 s^-2 value from one constants owner; it is a
fixed declared CODATA central value, not an exact SI definition or physical
uncertainty prediction. Physical versus comoving Sigma are distinct: a
comoving surface convention would divide by(1+zl)² and is omitted. Do not
invent an exact Msun mass: NFW composition requires an explicitly supplied
mass-unit conversion; physical-SI Plummer composition is directly compatible.
D_time_delay is a geometric coefficient; actual delay also requires a Fermat
potential difference. No lens mass/light/kinematic law, shear catalogue,
photometric-redshift population, mass inference or observed lens is qualified.

Retain ordered pair rows including invalid/refused rows, requested-output
masks, source/provider identity and owned input pair lineage. Dl/Ds/Dls,
beta, Sigma and D_time_delay are independently status-bearing outputs.
Foreground/equal-redshift pairs are outside this bounded finite-prefactor
slice; never silently turn an unresolved positive interval into a zero-efficiency
source or drop that row. A foreground source population law is separate.

Require positive resolved Dl,Ds,Dls lower intervals for finite prefactors.
Closed geometry must stay strictly before the first conjugate point for all
requested observer paths, including numerical phase intervals. Keep the
permitted branch beyond the transverse maximum: cos can be negative for
pi/2<sqrt(-K)*X_s<pi. Never recover that cosine via an unsigned square root
or apply a flat distance subtraction. beta may exceed1 in curved geometry;
no artificial upper cutoff. Direct DeltaX resolves close pairs; refine under
original quotas or refuse if its distance interval includes0.

Frozen empirical admission allocations (approved before implementation):
  distances and D_time_delay: 1e-8 Mpc + 1e-8 relative
  beta: 1e-12 + 1e-8 relative
  Sigma_crit: 1e-10 kg/m² + 1e-8 relative
Interval-style monotone positive-factor propagation includes inherited provider
errors, radial-to-S_K response, finite series/truncation/arithmetic and casts,
unit conversion and positive quotient denominators. Common-source correlations
can make this conservative; never invent independence. Existing provider
budgets/physical identity are unchanged. Empirical inherited errors are not
continuous certified bounds. Native bytes, points/pairs, per-pair and total
outer/FD callbacks and depth remain explicit bounded policies; failed nodes
and work are retained under original quotas. No implicit larger resource cap.

Independent controls before acceptance: analytic flat EdS and open Milne
Dls=DH/(1+zs)*sinh(log[(1+zs)/(1+zl)]); independent high-precision direct
FLRW/FD integrals and matched-source nonflat CLASS distance/clock controls;
near-flat +/-K joins; close epochs; closed permitted branch beyond pi/2;
closed first-conjugate approach/refusal; beta>1; inverse/unit/physical-comoving
checks and SI profile composition. Reference refinement and propagated
uncertainty must occupy <=5% of EACH final output allocation, including Sigma
and D_time_delay, rather than merely qualifying distances. Hostile inputs,
mask omissions, lifetime/move/copy, arithmetic and strict callback/payload
refusals remain required. Actual measurement substitution requires its own
source/mass/line-of-sight/calibration/selection and likelihood gates.

Source verification note: Prospector independently read/rendered exact Suyu
arXiv:0910.2773v2, sections2.1--2.2/Eqs1--7, supporting physical Sigma and
D_time_delay plus Fermat and mass-sheet limitations. Hogg9905116v4 Eq19 is
explicitly restricted to Omega_k>=0 and is not used as a closed-space shortcut.
Pinned CLASS angular_distance_from_to directly evaluates sin/sinh of the
radial difference and source-redshift denominator. Exact hashes/coverage are
owned in the separate thin-lens-geometry-source-locators.json receipt.
The explicit Mpc-to-metres divisor in Sigma is a unit clarification before
implementation/reference execution; proposed numerical allocations are unchanged.


`ThinLensGeometry` copies one actual provider after a bounded payload check.
`thermal_provider()` and `curved_provider()` expose its immutable source; the
inapplicable getter returns null. Copies own independent source storage; moves
invalidate their source, while self move preserves it. `evaluate(pairs, mask,
policy)` retains pair order, all refused rows and requested-output statuses.
An invalid mask or structural resource refusal returns no rows. Each admitted
pair computes three radial intervals for the same observer geometry; unrequested
final prefactors stay absent. There is no ruler request or photon/FD remapping.
The thermal momentum-method policy must agree with the retained provider.

The default policy bounds 2048 pairs, 16 MiB of declared native payload, 20M
callbacks per pair and 100M callbacks overall; nested thermal work also obeys
its independent policy quota. All three queries debit the original pair and
whole-batch quotas, including failed callbacks. `retained_payload_bound()`
charges actual visible vector/string capacities; `evaluation_payload_bound(n)`
adds conservative row and inherited provider scratch. These exclude borrowed
inputs, allocator metadata, quadrature recursion stack and RSS.

```cpp
using namespace irred::cosmology;
ThinLensGeometry geometry(CurvedFLRW({70, .3, .0001, .6999}));
const LensEpochPair pair{.5, 2};
const auto result = geometry.evaluate({&pair, 1},
    lens_geometry_mask(LensGeometryOutput::critical_surface_density_kg_m2) |
    lens_geometry_mask(LensGeometryOutput::time_delay_distance_mpc));
// Inspect status and optional values before composing a lens/source law.
```

For the initial reference execution, the independent script uses fixed
Gauss--Legendre quadrature in log(1+z) at 60/80 decimal digits, radial orders
64/96 and FD momentum orders160/240 on q=0..80. The omitted positive FD tail,
including today's normalization contribution, is explicitly propagated.
Refinement and that tail must fit the unchanged 5% final-output allocation.
The CLASS witness uses separately refined background ODE/table controls and
matched physical coefficients; table subtraction is used only for its moderate
separated-epoch witness. Production always integrates the interval directly.
The nextafter-close analytic/reference witness never uses table subtraction.
These sources are original validation drivers, not a Python production API.


The permanent `thin_lens_geometry_reference_contract` retains 204 comparisons:
150 high-precision outputs across 25 EdS/Milne/open/closed/closed-dust/FD pairs,
and 54 matched CLASS outputs across nine separated FD/open/closed pairs.
The independent reference refinement occupied at most3.238e-10 of each final
allocation; native versus high precision at most2.010e-6. Refined CLASS
uncertainty occupied at most0.0001245 (below0.05), and CLASS versus high
precision with refinement at most0.007830 (below1). These are finite empirical
controls, not whole-domain certification or observational qualification.

The first native attempt refused adjacent-double curved epochs because the
original double-z quadrature could not subdivide their support; the failure,
callbacks and original witness are preserved. The direct two-epoch helper now
uses t=0..1, wide z=zl+(zs-zl)t, and multiplies both integral and error by
zs-zl. Observer integrations keep their existing z mesh. Thermal intervals
use log1p((zs-zl)/(1+zl)) rather than subtracting large logarithms. No budgets
or quotas were weakened. The shared curvature series remainder for |K X²|<1e-4
is below the retained128-binary64-epsilon arithmetic allocation.

The existing thermal provider owns the mapped physical state and its empirical
numerical estimates. This consumer inherits that exact emitted state; it does
not add physical-input, G, redshift or source-conversion uncertainty, or certify
an upstream convention against a released lens system. Such uncertainties need
a separately declared source law. The independent physical-source calculations
also test mapping/cast differences under the final numerical allocations.
