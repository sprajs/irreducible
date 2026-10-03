# A declared effective-fluid background and matched finite ruler

The standalone C++20 `irred/effective_fluid.hpp` API implements an independently
authored compensated density shape. Its identity is
`independent-compensated-density-nine-halves-thermal-FD-background/v1`.
It is a background approximation inspired by a limited EDE source review;
the original scalar normalization, initial field coordinate, modified CLASS,
thermal history, perturbations and full-fit mapping remain unresolved.

Map an explicit `ThermalPhysicalModel` once with the existing
`map_thermal_physical_model`, then supply its complete emitted `ThermalFlatModel`
in `EffectiveFluidRequest`. Read [thermal observables](thermal-observables.md)
for physical versus fractional densities, fixed constants, Kelvin/eV conversion
and explicit collisionless zero-chemical-potential FD species. Retain every
species and its order, mass, temperature and populated-state weight. The other
massless fraction excludes photons and explicit species. The caller retains
the original physical input and mapping witnesses; this nominal emitted model
does not propagate their arithmetic metadata as physical uncertainty.

For a0=1, X=rho_fluid/rho_crit,ref,0, Xc>=0 and 0<ac<1, the compiled law is

    X(a) = 2 Xc / (1 + (a/ac)^(9/2))
    w(a) = -1 + (3/2) (a/ac)^(9/2) / (1 + (a/ac)^(9/2))
    dX/dln(a) = -(9/2) X(a) (a/ac)^(9/2) / (1 + (a/ac)^(9/2))
    E_new(a)^2 = E_ref(a)^2 + X(a) - X(1)
    Omega_Lambda,new = Omega_Lambda,ref - X(1).

The conservation equation `dX/dln(a)+3(1+w)X=0` includes zero amplitude without
division by density. X(ac)=Xc and w(ac)=-1/4. The early density tends to 2Xc;
the late density dilutes as a^(-9/2). No scalar perturbation sound speed follows
from this background equation of state. The actual retained wide reference
Lambda is captured once, and negative or numerically unresolved compensated
closure is refused. H0 and all standard-sector coordinates stay fixed.

`scaled_expansion` supplies P=a^4 E_new^2 on 0<=a<=1 using the same retained
thermal background. E/H queries require 0<a<=1. The zero-amplitude path preserves
the reference scaled values, diagnostics and work exactly. The radiation
endpoint and today P=1 also reuse the reference. Positive fluid contributions
use bounded ratio powers and a cancellation-resistant near-today expression.
Unresolved positive products are refused. Diagnostics remain empirical
operation/libm and sampled quadrature estimates, with no certified continuous
background or fluid-averaging error claim.

Supply a fixed `ruler_endpoint` with 0<a_end<=1 and nonempty source/endpoint
origins. The existing shared loading/ruler helper integrates

    r_s = c/(H0 sqrt(3)) integral_0^a_end da /
          [sqrt(P(a)) sqrt(1+3 Omega_b a/(4 Omega_gamma))].

The ruler is in Mpc. The same complete species, baryon/photon loading and endpoint
define a reference comparison. Because X(a)>=X(1) on this domain, the compensated
expansion cannot decrease and the matched finite ruler cannot increase. This
comparison does not predict a drag or last-scattering endpoint, a thermal
history, CMB spectra or an observational fit.

`evaluate` takes an ordered coarse scale-factor batch. Masks `thermal_e`,
`thermal_h` and `effective_fluid_ruler` request E, H in km/s/Mpc and one fixed
batch ruler. Omitted values stay absent. Individual invalid rows and scalar
refusals retain their order; a ruler refusal preserves independent successful
scalars. Preparation retains admitted source and performed work on later
closure failure. Copies own storage, moves invalidate the old complete owner,
and self move preserves it.

The unchanged `ThermalObservablePolicy` controls tolerances and resources:
default E/H allocation 2e-10 relative, ruler 1e-8 Mpc+2e-10 relative, 20 million
combined callbacks per point/ruler, 100 million per batch and momentum work,
depth30, 4096 points and 16MiB owned payload. Preparation counts are separate;
all failed evaluation attempts count. The standalone thermal default stays
4 million. Momentum selection must match preparation even for empty/today
batches. Payload accounting includes source copies, retained normalization,
rows, labels and fixed scratch, while excluding borrowed inputs, allocator
metadata, stack and RSS. No table, server or runtime equation framework is used.
The payload bound applies to preparation/evaluation and move-return of this
retained owner with its operation temporaries. Additional caller-owned copies
and copy assignment's arbitrary old-target/replacement coexistence are outside
that policy envelope; callers own those lifetimes. Copy semantics carry no
resource-policy parameter.

The permanent source controls cover exact zero/today/transition/conservation,
ordered duplication, domains, closure, resource refusal, masks and lifetime.
The peer uses an independently written log-density and direct-a fixed Simpson
integration of an exactly dyadic empty-species polynomial reference, with each
4096/8192 refinement, positive dependency propagation, sum/scale arithmetic and
measured reference reporting cast loss together occupy at most 5% of each frozen
ruler allocation. It also requires the discrepancy to fit the actual native
reported diagnostic plus that complete empirical reference estimate, retaining
the original whole-allocation comparison. The log/exp/sqrt operation estimates
are conditional empirical allowances; continuous libm qualification stays open.
The public SDK
consumer maps one explicit synthetic positive-mass relic source and checks a
matched zero/positive-amplitude background and finite ruler, and preserves the
earned positive-FD normalization work on a later negative-closure refusal. Shared thermal/FD
ancestry is not an independent massive-relic certificate. These new controls
are written but have not yet been compiled or executed; runtime qualification
requires the frozen source review and a separately allocated compute job.
