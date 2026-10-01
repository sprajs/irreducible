# One thermal background for distances and a supplied-drag ruler

The standalone C++20 `irred/thermal_observables.hpp` API prepares one flat
photon, baryon, CDM, explicitly supplied other massless radiation, collisionless
thermal-relic and Lambda state. That same retained state supplies E, H,
D_H, D_M, D_L, D_V, the conditional sound horizon and three distance/ruler
ratios. It has no CLI or C ABI route. It predicts no recombination, ionization,
drag epoch, growth, perturbation spectrum, fit or posterior.

The physical source is `ThermalPhysicalModel`. It requires finite positive
H0 in km/s/Mpc and Tcmb in Kelvin, separately supplied nonnegative physical
baryon/CDM/other-massless densities, and at most sixteen explicit relic species.
Here a physical density is `omega_i = Omega_i h^2`, with `h=H0/100`; the mapped
`ThermalFlatModel::omega_*` fields instead hold today's critical-density
fractions Omega_i. There is no total matter-density input. Relics are not
included in the supplied baryon or CDM densities. The other massless component
contains no photons and no explicitly supplied species.

Each species declares its mass in eV, positive present momentum temperature
in Kelvin and positive statistical weight g. Its frozen distribution is the
zero-chemical-potential Fermi–Dirac law described in
[thermal neutrino](thermal-neutrino.md). g counts populated states; g=2
can describe one neutrino/antineutrino pair with the same distribution. The
API has no implicit temperature ratio, effective neutrino number, hierarchy,
mass sum, fitted 93.14-eV conversion or external-engine default. A caller that
starts from a temperature ratio must first record the ratio and derive the
explicit species Kelvin temperature from its declared Tcmb.

The mapping uses the exact SI kB value 1.380649e-23 J/K and exact
1.602176634e-19 J/eV, with the same existing h, c, IAU Mpc and fixed CODATA2018
G convention as the thermal background. It sets T_eV=kB*T_K/eV_J and
rho_gamma,0=(pi^2/15)*T_eV^4, then divides by the shared critical energy density.
It divides the other physical densities by h^2 and converts species temperatures
once. Positive values that cannot be stored as normal binary64 fractions or
temperatures fail; zero physical densities remain valid controls. The mapping
identity is `thermal_physical_mapping_id`. This is a declared conversion and
source convention, not full Planck base LCDM.

A `ThermalObservableRequest` retains this physical source, finite z_drag>=0,
nonempty `drag_origin` and nonempty `source_origin`. Origins are caller-supplied
provenance labels and do not certify their contents. Structural and payload
admission precede source copies. An admitted source remains available even if
mapping or physical closure fails. Successful preparation retains its own
physical source and one independently usable `ThermalBackground`, including the
present relic normalization and its closure diagnostic. Copies own independent
storage. Moves transfer source and background together and invalidate the old
owner; self move preserves the owner.

The thermal kernel alone owns

\[
P(a)=a^4E(a)^2=\Omega_\gamma+\Omega_{\mathrm{other,massless}}
 +(\Omega_b+\Omega_{\mathrm{cdm}})a+\Omega_\Lambda a^4
 +\sum_i\frac{g_iT_{i0}^4}{2\pi^2\rho_{\mathrm{crit},0}}
 I_\rho(m_i a/T_{i0}).
\]

`scaled_expansion` evaluates this wide numerical coordinate on 0<=a<=1.
At a=0 the exact massless FD moment gives its finite radiation limit; the
low-level radiation-free owner returns its exact finite P(0)=0 coordinate.
Nonzero scaled coordinates and their diagnostics must remain normal in the
wide arithmetic profile; positive-a underflow is refused.
The inverse-distance/ruler consumer instead requires its positive supplied
photon component. No finite E(0) is fabricated; public E/H queries continue to
require a>0. At a=1
the retained model normalization gives E=1 and H=H0 exactly without momentum
work. A nonnegative flat Lambda closure is required and never clipped.

For finite z>=0, a=1/(1+z), the declared flat FLRW and photon-conserving geometry
uses H=H0*E, D_H=c/H, D_M=(c/H0)*integral_0^z dz'/E(z'),
D_L=(1+z)*D_M and D_V=cbrt(D_M^2*z*D_H). Distances are Mpc and
c=299792.458 km/s. The distance quadrature uses x=log1p(z), t in [0,1]
and integrand exp(-x*t)/sqrt(P(exp(-x*t))), with prefactor (c/H0)*x.
This avoids cancellation at tiny positive z and has no high-redshift cutoff.
At z=0 the transverse, luminosity and volume distances and their ratios are
exactly zero.

With a_drag=1/(1+z_drag), the conditional tight-coupling ruler is

\[
r_s=\frac{c}{H_0\sqrt3}\int_0^{a_\mathrm{drag}}
 \frac{da}{\sqrt{P(a)}\sqrt{1+3\Omega_ba/(4\Omega_\gamma)}}.
\]

The same photon and baryon fractions define the baryon loading. The finite
scaled radiation endpoint enters this integral directly. The supplied drag
redshift is not inferred by the engine. z_drag=0 is a mathematical control of
this prescribed approximation; it does not claim late baryon/photon tight
coupling. The ratios D_M/r_s, D_H/r_s and D_V/r_s use this same source identity.

The nine `EarlyLateOutput` coordinates and masks retain their existing units
and meanings. `thermal_ruler_mask` requests the single optional batch ruler
directly; any ratio also requests it once. Ruler-only batches perform no
per-redshift momentum work, while preserving admitted row order. A distance-only
request runs no separate point E evaluation. D_H-only runs no distance
quadrature. Unrequested scalars stay absent. Batch success means processing
admission; each requested scalar retains its own status and optional value.
A failed ruler does not relabel a separately successful distance. Positive
physical values and positive diagnostics cannot silently underflow to zero.

`ThermalObservablePolicy` bounds points, species, per-point combined work,
whole-batch combined work, momentum-only whole-batch work, outer depth and
nested momentum depth/callback work. The batch reports `outer_callbacks`,
`momentum_callbacks` and their sum `callbacks`; each row's callbacks are its
combined work, excluding the once-per-batch ruler. Preparation callback counts
remain available from the retained background. The conservative payload bound
includes retained physical/mapped sources, preparation temporaries, origin
copies, row capacities and bounded nested nonrecursive scratch; allocator
metadata, recursion stack and RSS are excluded. No unbounded cache or table
is built.

Default external acceptance allocations are 2e-10 relative for E/H,
1e-8 Mpc + 2e-10 relative for distances/r_s, and 1e-10 + 5e-10 relative
for ratios. Production outer integration reserves one eighth of the requested
absolute/relative distance allocation, including luminosity amplification.
The nested normalized momentum policy remains 1e-12 + 2e-12 relative, with
an explicit default 100-million aggregate callback ceiling for this consumer
(20 million combined per point). The standalone background policy retains its
4-million aggregate default. Its normalization and species diagnostics propagate
additively into P. Each positive
integral retains its maximum sampled inverse-sqrt dependency radius, added to
outer quadrature and arithmetic/storage diagnostics. Derived D_H, D_V and ratios
use positive intervals. These are empirical numerical admission diagnostics,
not certified continuous error bounds or independent random uncertainties.
Unattainable or exceeded user budgets fail without changing tolerance.

Owner regressions retain mapping/units, finite radiation endpoints, geometric
identities, output omission, independent partial failures, lifetime rules,
invalid domains, positive representability and nested resource limits. The
peer suite uses originally written direct-redshift/sqrt(a) fixed-panel controls,
independent high-precision FD evaluations and a matched external CLASS source
with explicit constants/species/temperature/statistical-amplitude conventions.
Reference refinement must occupy at most 5% of each frozen comparison allocation;
failed reference refinements are retained as failures, not expected values.
No external implementation code or assets enter the library. These named
comparisons remain conditional on the supplied physical source and drag epoch;
they do not establish recombination closure, arbitrary-domain accuracy,
observational qualification or a posterior.
