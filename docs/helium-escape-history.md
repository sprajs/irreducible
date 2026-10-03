# Conditional helium escape history

This source candidate evolves a supplied HeII/He population history with HeI
584 Å resonance escape, hydrogen continuum absorption and 591 Å intercombination.
It has not yet been compiled or numerically accepted. Build consumers are
registered in source; installed execution, independent trajectory and same-run
CLASS comparison are pending.
It is a distinct supplied-driver approximation; the existing
[singlet H/He history](hydrogen-helium-history.md) and its defaults keep their identity.

## Inputs and the useful result

`prepare_helium_escape_history` borrows a `HeliumEscapeHistoryRequest` and owns one
immutable driver in the returned `HeliumEscapeHistory`. Driver knots retain their
original order and binary64 scalars: decreasing redshift, H in s^-1, physical nH
in m^-3, radiation temperature Tr in K and xHI=HI/H. Supply f=nHe/nH, initial
q=HeII/He and one nonempty producer identity. Declare the role as supplied
cosmological history or synthetic prescribed bath. That role is not a certificate.
No H0, abundance, nuclei mass or neutrino composition is inferred.

The driver is geometric interpolation of H,nH,Tr,xHI in ell=ln(a)=-ln(1+z), with
no extrapolation. The emitted finite driver defines the conditional IVP. Refining
an external producer grid changes that supplied input law and needs separate
identity/error evidence. Rates use Tm=Tr. Actual producer Tmat departures belong
in its scope receipt; matter temperature is not evolved here.

Evaluate an ordered coarse redshift span, choosing `helium_escape_fraction` (1),
`helium_escape_electron_density` (2) and `helium_escape_thomson_opacity` (4).
Each requested output has its own status, optional value and absolute numerical
diagnostic. Invalid query rows retain valid neighbours. A refused or unrequested
value is absent, never a zero or a dropped row. A failed preparation supplies no
trajectory but retains lawful acquired source metadata, work and partial witnesses.

The source pointer is const. Copies acquire independent storage with simultaneous
owner payload checks; copy resource refusal throws `length_error`, allocation
failure throws `bad_alloc`. Copy assignment retains its destination until a
complete replacement is acquired. Moves invalidate the original. Ordinary copies
do not rerun the physical preparation or report new scientific evidence.

## Population and rate ownership

Let p=1-xHI, y=fq and xHeI=f(1-q). One shared native charge equation forms
ne=nH*p+nHe*q, nHe=nH*f; the rate uses this same ne/nH. Fixed alpha/me ratios are
one. Continuum-equilibrated excited populations and detailed balance give the
net helium ground-state formation per H per second

\[
 F_{\rm He}=Y_{\downarrow}(x_e y-x_{\rm HeI}s),\qquad
 \frac{dq}{d\ln a}=\frac{Y_{\downarrow}[(1-q)s-qx_e]}{H},\qquad
 \frac{dq}{dz}=-\frac{1}{1+z}\frac{dq}{d\ln a}.
\]

Proper time satisfies dln(a)/dt=H. The same shared charge emits ne and
qT=c*sigmaT*ne/[H(1+z)]. There is no separate hydrogen derivative in this helium
IVP; an electron derivative would be -dxHI/dt+f*dq/dt. This does not feed changed
helium charge back into the supplied H-neutral or thermal driver.

`helium_escape_rates.hpp` owns the equation and its analytic q derivative. Its
finite scalar facts retain CLASS 3.3.0 commit 0ceb7a9 / HyRec2020 helium.c SHA
cc113cdc4b2d4ffa0762ab82d5485816bbbb3bbf63cd3441be5206c4c3bf6587.
The existing quantum-density owner supplies the T^(3/2) shape, normalized to the
source 2.414194e21 m^-3 at 1 K. This does not replace that source coefficient with a
new mass convention. Atomic energies in Kelvin, widths and escape coefficients
are fixed source facts without a qualified joint uncertainty law.

With s0=4Q/nH, s=s0 exp(-285325/Tr), y2s=exp(46090/Tr)/s0 and
y2p=3 exp(39101/Tr)/s0, the downward rate is
Ydown=50.94*y2s+1.7989e9*y2p*P_eff in s^-1. The 50.94/s two-photon asset differs
from the old singlet 51.3/s; two-photon decay is not a newly added mechanism.
The newly included continuum/intercombination escape uses

\[
 I=\eta_c^{-1}=9.15776\,10^{28}\frac{H}{n_Hx_{\rm HI}},\quad
 \tau=4.277\,10^{-14}\frac{n_Hf(1-q)}H,\quad
 \tau_c=\frac{\Gamma_{\rm inc}\tau}{4\pi^2 I}.
\]

I has units s^-1; tau and tau_c are dimensionless. The incoherent width includes
2P→2S coupling with its stimulated term, then 3S/3D/4S/4D/5S/5D contributions.
The complete seven-term width and finite constants are in the shared rate owner;
no external tables or runtime formulas are input.

The enhancement is E=sqrt(1+pi²tau_c)+7.74tau_c/(1+70tau_c). For
B=1-exp(-1.023e-7tau), phi=6.14e13/I and R=.964525exp(2947/Tr),

\[
 P_{\rm eff}=\frac{E(1-B e^{-\phi})+BR}{\tau}.
\]

The 591 line removes a fraction of the 584 contribution while adding its own
escape. The combined numerator is evaluated positively without clipping the
signed intercombination increment. P_eff is an effective multiplier, not an
independent probability capped at one.

## Domain and physical limits

The bounded profile is 1900<=z<=2800, first knot 2600–2800, last 1900–2200;
4300<=Tr<=7800K, 5e8<=nH<=8e9m^-3, 1e-14<=H<=1e-12s^-1 and 0.04<=f<=0.12.
Require strict 0<xHI,q<1, tau>=100 and f*q>=1e-6. The latter matches the current
CLASS active-stage comparison threshold; this model refuses an interval crossing
it instead of silently freezing, clipping or extending its domain. These ranges
are engineering/source restrictions, not a certified physical validity region.

Every attempted rate, including rejected numerical trials, contributes to work
and retains its temperature/tau range. Witness endpoints are outward binary64
projections of actual long-double coordinates. The excluded HeIII activity
Q exp(-631462.7/Tr)/ne is checked at every retained node, with maximum log<=log(1e-12).
It is a scope witness, not a full kinetic population bound between nodes. Query
scope checks use the actual queried state as well. No HeIII, triplet populations,
H kinetics, heat feedback, reionization, integrated optical depth, drag epoch,
survival, visibility or CMB spectrum is supplied.

The acquired arXiv1011.3758v2 helium section supports Tm=Tr, steady radiation,
continuum-equilibrated excited levels, damping wings, neglected resonant scattering
and a blackbody continuum source. Its Eq.104 assumes HI Saha/xp approximately one;
current source explicitly accepts actual xHI since 2012. The selected model owns
that conditional generalization, not a fully coupled nonequilibrium HI emissivity
or photon/heat solution. Omitted higher levels, collisions and atomic/model errors
remain separate from numerical error.

The primary states the enhancement fit within 0.8% for positive tau_c in its
simplified transfer problem. Its historical maximum 0.3% full-history comparison
is attributed to its own staging/cosmology/reference. Neither is a numerical
budget or physical qualification of this consumer. The acquired Eq.99/121
hydrogen-sign/convention tension is preserved in source receipts; the owned
helium population derivation omits that ambiguous total-electron term.

## Numerical diagnostics, resources and acceptance

The scalar backward-Euler solve in ln(a) uses the shared analytic derivative and
physical-box safeguarded Newton corrections. All supplied knots are boundaries
on three nested meshes; Richardson2 and the three-mesh difference estimate
refinement. Root corrections, finite-driver/libm/arithmetic diagnostics and
query curvature are retained. These are empirical estimates, not certified global
bounds. Failure/refinement refusals must be preserved. The output cast is charged.
Shared ne/opacity inherit raw q errors even when q's separately requested group
refuses its own allocation.

Allocations are Eq=1e-8+2e-6|q|, Ene=nHe*1e-8+2e-6|ne| and
EqT=2e-9+3e-6|qT|. For EACH comparison the complete numerical reference error,
including refinement, arithmetic, casts, boundary, interpolation and projection,
must be<=5% of that allocation. The discrepancy must fit native diagnostic plus
that combined reference error. Separate components do not each earn5% slack.
Signed scalar-rate controls use the declared absolute positive reaction scale,
not a relative error to a stationary zero. Atomic/physical/source-law uncertainty
is not hidden in these numerical allocations.

Defaults remain 4096 base intervals, hard finest 65536, at most 4097 source knots,
4096 query rows, 4M charged source/driver/rate evaluations and 32 MiB requested owned
payload, with hard 1 GiB. The owner reports preparation work; each query batch reports
its actual driver evaluations and checks preparation plus that batch against the
same 4M limit. Repeated independent batches have no mutable lifetime counter;
a consumer controller must report their complete aggregate work. Work includes
failed attempts; no fresh background or ODE evaluation occurs at queries.
Remaining work is checked after source import before numerical setup. Requested payload is preflighted;
actual string/vector capacities and simultaneous preparation/output/replacement
owners are checked. This excludes allocator overhead/RSS; complete allocation-
fault and measured peak validation is still pending. Rejected trajectory/output
storage is released, including its capacity. A captured source survives a later
refusal when its retained capacity fits the original cap. If acquisition itself
exceeds that cap, rejected source reservations are released; lawful metadata,
work and attempted witnesses remain, while source_capture_complete is false.
This refusal owner is not a whole-process proof for an arbitrarily tiny cap.
Nearest strict FP,
binary64 inputs and binary long double with at least 64 bits of precision and
maximum exponent at least16384 are required; unsupported arithmetic is a refusal
at preparation and evaluation.

Native controls declare stationary analytic and nonstationary synthetic histories,
independent source-shaped scalar algebra and an original adaptive order3 Radau-IIA
reference with separate driver/interpolation/root code. They share physical
coefficients/libm ancestry, not a separate physical theory. The reference keeps
partial failures and all attempted work; every requested q/ne/opacity output
must pass its combined error gate. Source controls and installed consumer remain
unexecuted until a root compute lease.

## Same-run reference and source notices

Repro owns a proposed same-run CLASS stage exporter: raw HII/H, HeII/He, H-neutral,
fHe/nH0/mass convention, H/s, Tr/Tmat, boundary, phase, active/frozen branch and
rate diagnostics. Exporter implementation and execution remain pending. Required
locations are thermodynamics.c:2593–2678/2063/2114,
wrap_hyrec.c:174/193–200 and helium.c:118–179. Raw states must be distinguished from
rejected trials and phase-smoothed table output. Capture the current workspace
immediately after the successful derivative call at 2947–2949, before the previous
stage evaluation at 2973 overwrites it; commit only after the complete source
callback succeeds. Retain derivative signs before the minus-z inversion and keep
repeated writes/generations, partial failure and whole-run completion separate.
Standard total xe cannot supply those ion fields. Keep CLASS fHe=YHe/[3.9715(1-YHe)] and its nuclei mass convention;
do not substitute native neutral masses or the value four. Input/build/source/precision
and full mapped comparison identity must be frozen before acceptance.

CLASS and native Thomson assets differ; report raw opacity and any explicitly
declared sigma-native/sigma-CLASS scalar projection separately. Producer grid,
phase, solver and projection numerical errors enter that comparison's complete
reference allowance. No shared coefficient ancestry establishes independent
physical evidence or full recombination/CMB qualification.

Primary: Ali-Haimoud/Hirata, DOI10.1103/PhysRevD.83.043513, acquired
arXiv1011.3758v2 HTML SHA cfc5a7b8f940968d5168e5c3aa422cbbe8f035af3feab957bf788eb6e32b7e38.
The continuum delegate is DOI10.1103/PhysRevD.77.083006; .083008 is the historical
full-history benchmark. Atomic source bodies/individual uncertainty generation,
journal/preprint equivalence and figure pixels remain unreviewed.
The matched upstream README at 09e8243 requests 2020/2011 companion citations but
supplies no explicit redistribution permission. This source is independently
authored equations; no HYREC code/assets or external table has been copied and
no permission is inferred from citation instructions. Full scientific scopes,
physical validity, inference and observations remain open.
