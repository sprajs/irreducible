# Collisionless thermal relic background

The standalone C++20 `irred/thermal_neutrino.hpp` API calculates the energy
density and pressure of explicit decoupled zero-chemical-potential Fermi–Dirac
species, then prepares a flat photon, massless-radiation, baryon, CDM, relic and
Lambda background. It supplies E(a) and H(a), with no CLI or C ABI route.
This remains the background prerequisite. The separate
[thermal observables](thermal-observables.md) consumer now uses this same retained
state for distances and a supplied-drag ruler. Neither provider supplies
recombination, ionization, drag prediction, perturbations, growth or CMB spectra.

A species supplies finite mass `mass_ev >= 0`, present momentum temperature
`temperature_today_ev > 0` and explicit positive `statistical_weight`. The
weight counts populated states: g=2 describes a neutrino plus antineutrino with
the same distribution. It is not a mass sum or an effective neutrino number.
There is no implicit Tcmb, Neff, temperature ratio, hierarchy, degeneracy or
93.14-eV mapping. The explicit `ThermalPhysicalModel` mapping supplies Kelvin
temperatures and physical baryon/CDM/other-massless densities omega_i=Omega_i h^2;
it converts with exact SI kB=1.380649e-23 J/K and derives the photon fraction
from the supplied blackbody Tcmb. See the thermal-observables guide for this
source convention. The frozen distribution is f(q)=1/(exp(q)+1), q=p a/Tnu0;
its momentum temperature scales as 1/a even after the species becomes
nonrelativistic. No interactions or production history are evolved.

In natural units hbar=c=kB=1, define y=m a/Tnu0 and

\[
I_\rho(y)=\int_0^\infty\frac{q^2\sqrt{q^2+y^2}}{e^q+1}\,dq,\qquad
I_P(y)=\frac13\int_0^\infty\frac{q^4}{\sqrt{q^2+y^2}(e^q+1)}\,dq.
\]

The species state is

\[
\rho(a)=\frac{g}{2\pi^2}\left(\frac{T_{\nu0}}a\right)^4I_\rho(y),\qquad
P(a)=\frac{g}{2\pi^2}\left(\frac{T_{\nu0}}a\right)^4I_P(y).
\]

Returned physical values are energy density and pressure in eV^4, rather than
eV/m^3 or mass density. The exact massless controls are
I_rho(0)=7 pi^4/120 and I_P(0)=I_rho(0)/3. The nonrelativistic leading limits
are I_rho/y -> 3 zeta(3)/2 and y I_P -> 15 zeta(5)/2; the production routine
integrates the distribution rather than replacing it with these limits.
Differentiating the frozen distribution gives
`d rho / d ln a = -3 (rho + P)`. Continuity and both limits are permanent
numerical controls.

A prepared `ThermalBackground` copies its admitted source once and retains
its present relic normalization. It distinguishes four nonnegative supplied
fractions `omega_gamma`, `omega_massless_nonphoton`, `omega_b` and `omega_cdm`.
Photons are separate from the supplied nonphoton massless component; baryons
and CDM contain no relic contribution. Each fraction refers to today's critical
energy density. There is deliberately no total Omega_m input.

\[
\rho_{\mathrm{crit},0}=\frac{3H_0^2c^2}{8\pi G},\quad
\Omega_\Lambda=1-\Omega_\gamma-\Omega_{\mathrm{massless,nonphoton}}
-\Omega_b-\Omega_{\mathrm{cdm}}
-\sum_i\frac{\rho_i(1)}{\rho_{\mathrm{crit},0}},
\]

\[
E^2(a)=\frac{\Omega_\gamma+\Omega_{\mathrm{massless,nonphoton}}}{a^4}
+\frac{\Omega_b+\Omega_{\mathrm{cdm}}}{a^3}
+\sum_i\frac{\rho_i(a)}{\rho_{\mathrm{crit},0}}
+\Omega_\Lambda,\qquad H(a)=H_0 E(a).
\]

H0 is supplied in km/s/Mpc. The implementation reuses the engine's exact
SI speed of light, Planck constant and IAU parsec convention. It fixes the
SI elementary charge conversion to 1.602176634e-19 J/eV and the
[CODATA 2018 G value](https://pml.nist.gov/cuu/pdf/RevModPhys.93.025010.pdf)
to 6.67430e-11 m^3 kg^-1 s^-2. Converting natural energy density to SI uses
(eV_J)^4/(hbar c)^3. G is a measured constant: treating this central value as
fixed is an explicit model/conversion convention. Its physical uncertainty is
not covered by numerical quadrature estimates. The current constants identity is
`SI2019-exact-h-c-kB-eV-IAU2012-AU-CODATA2018-G-fixed`, including the explicit
kB conversion now used by the physical mapping. The existing eV-source
background retains the same numerical constants; historical receipts bearing
`SI2019-exact-h-c-eV-IAU2012-AU-CODATA2018-G-fixed` remain immutable.
These constants are identified by `thermal_neutrino_constants_id`; the compiled
numerical method and arithmetic profile have separate `thermal_neutrino_method_id` and
`thermal_neutrino_arithmetic_id` identities.

All physical species-state and E/H background queries require finite 0<a<=1, finite
positive H0 and a nonnegative flat Lambda closure. Zero supplied fractions are
valid analytic controls. At most sixteen explicit species are supported.
Preparation subtracts the largest nonnegative contributions first, rejecting
a strictly positive closure excess even when a naive sum would round to one.
A closure within the present-density numerical diagnostic is refused rather
than clipped or silently replaced. The prepared normalization gives E(1)=1
and H(1)=H0 with no new momentum work. The same owner exposes the wide
scaled coordinate P(a)=a^4 E(a)^2 on 0<=a<=1, using the exact FD massless
moment at the finite radiation endpoint a=0. This does not make E(0) finite.
E/H and the separate distance/ruler consumer reuse that one scaled equation.

A prepared owner supports independent copies. A move transfers its species,
normalization and numerical state together, leaving the source invalid with
no usable fractions or evaluation rows. Self move preserves the owner. These
lifetime rules prevent a lost species vector from being evaluated with a
stale successful closure.

`evaluate_thermal_moments` takes finite binary64 y>=0. Production integrates
I_rho/s and s I_P, where s=hypot(1,y), with binary64 callbacks and wide scaled
arithmetic. Direct momentum intervals of width four resolve the known smooth
FD support before adaptive quadrature. A bounded search selects a cutoff up to
256. The omitted tail uses f(q)<=exp(-q), sqrt(q^2+y^2)<=q+y and the
upper envelopes for inverse energy, integrated with integer incomplete-gamma
polynomials. This is an analytic envelope for the omitted mathematical tail;
the quadrature estimate and floating-point diagnostics are empirical, not
universal error bounds.

The policy applies absolute plus relative tolerance separately to each
normalized moment. The total diagnostic adds quadrature, omitted-tail and
arithmetic/cast terms. An arithmetic allowance refuses unattainable budgets.
No user tolerance is weakened. The policy bounds callbacks per momentum
evaluation and across preparation or a complete background batch, depth,
points, species and a conservative simultaneously owned payload. Stack,
allocator metadata and RSS are outside that payload allowance. The arithmetic
profile requires FE_TONEAREST and long double with at least 64 mantissa bits
and exponent range 16384. A requested positive public value or diagnostic
that cannot be stored as a normal binary64 number fails; positive zero is not
a substitute for an unrepresentable result.

Background masks request E, H or both. Each requested scalar retains its
status and optional value; an unrequested scalar is absent. Dimensional H
failure leaves representable E usable. Batch success describes processing
admission, with individual query failures retained. Species dependence and
present closure diagnostics propagate additively to E and H, with no
independent-random-error assumption. The exact a=1 identities carry zero
numerical diagnostics because they reuse the model's defined normalization.

Named owner controls test analytic limits, continuity, state weight,
component closure, callback/payload limits, invalid inputs and dimensional
failures. Independent peer references use an originally written
high-precision momentum calculation, distinct from production quadrature;
reference refinement must occupy at most 5% of the frozen comparison
allocation, 1e-10 + 2e-10 times the normalized reference moment and 2e-10
relative for E and H. Permanent peer controls also freeze a matched
[CLASS v3.3.0 source](https://github.com/lesgourg/class_public/tree/0ceb7a9a4c1e444ef5d5d56a8328a0640be91b18)
background at z=0, 0.1, 1, 10, 100, 1000 and 1e6, with source-refinement
change at most 5.972e-12 relative. The fixture declares supplied photon and
massless fractions, baryon/CDM physical densities, one 0.06-eV species,
explicit momentum temperature and g=2. An independently derived constant
conversion adjusts the external statistical amplitude to match SI h, G, eV
and Mpc conventions. Its original inputs and refinement receipts remain
separate from production. This is a named background comparison, with no
implicit Planck or external-default mapping. Species, temperature, weight,
massless/photon accounting and constants must match explicitly.
The [CLASS IV background treatment](https://cds.cern.ch/record/1345134)
is a scientific source; the independently executed CLASS reference has its
own numerical approximations and ancestry. No external implementation code
or assets are copied into this library. Passing named checks does not establish
arbitrary-domain accuracy, full Planck base LCDM, a posterior or an
observational interpretation. See [scientific contracts](scientific-contracts.md)
and [testing](testing.md) for the distinction.
