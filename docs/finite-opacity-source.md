# Supplied finite-opacity scalar source

This C++ consumer evolves photons, cold baryons, CDM and smooth Lambda in a
flat GR background. It owns an explicitly supplied finite conformal opacity
history and one finite initial mode. The first implementation returns computed
trajectories and raw temperature/polarization source channels, with numerical
admission withheld until the missing integrated error and reference gates are
earned. It does not supply a full neutrino CLASS model, regular asymptotic
initial conditions, a recombination opacity law, a primordial spectrum or C_l.

`irred/finite_opacity_source.hpp` is independently usable. Move a
`FiniteOpacityRequest` into `prepare_finite_opacity_source`, inspect preparation
status and its immutable identity, then call `produce`. Both the producer and
result are move-only. Results share the original const identity and boundary;
they may outlive preparation. Producer moves invalidate prior getters, and
self move preserves them. A refused evaluation retains the original inputs,
mapped state, boundary, work and reached prefix. No partial grid is exported as
a complete `ContinuousCmbSource`.
The full reached F/G tails keep the matching per-k dimensionless core, eta and
actual a, including refusals between exported source endpoints.

The physical inputs are H0, physical omega_b/omega_cdm and Tcmb. Additional
massless density and explicit species must be absent. The existing mapper runs
once, the existing thermal background runs once, and its retained coefficient
receipt is captured once. Each actual stage P(a) query is charged before the
shared conformal conversion and its opt-in no-species diagnostic. Conditional
coefficient diagnostics do not establish their integrated source response.

The immutable binary64 opacity knots define linear K(eta), in Mpc^-1:
K=a n_e sigma_T times the metres per Mpc conversion. There is no c factor,
extrapolation, knot deletion or opacity rescaling. Require strictly increasing
positive eta and k, K in [0,10^12], total optical depth at most 512,
2^-512<=a_i<a_observer<=1 and k(eta_observer-eta_i)<=64. The supplied clock
endpoints are conditionally checked against the same retained thermal age;
the source integrates its actual a(eta) and retains observer-a mismatch. The
endpoint check and age refinement witnesses do not earn the missing clock
propagation error. Failed positive survival or narrowing is refused rather than
replaced by zero. K=0 is an exact supplied no-scattering law.

Conventions are conformal Newtonian
ds^2=a^2[-(1+2 psi)deta^2+(1-2 phi)dx^2], Fourier exp(+ik.x), theta=ik.v,
and MB F/G moments sum_l (-i)^l(2l+1)F_l P_l. F0=delta_gamma,
F1=4theta_gamma/(3k), F2=2sigma_gamma and Pi=2sigma_gamma+G0+G2.
The eleven core coordinates store theta/k, so state and full-stage residual
scales are dimensionless. Photons retain anisotropic stress:
psi=phi-6 Hcal^2 fg sigma_gamma/k^2; the shared momentum constraint supplies
phi'. Photon/baryon collision forces cancel in their enthalpy-weighted sum.
The compiled hierarchy retains the polarized collision terms and cold baryon
Euler law; it makes no tight-coupling or streaming switch.

The finite seed has zero relative entropy/slip, local zeta=-phi+delta_b/3=1,
and exactly zero F_l>=2 and every G_l. It satisfies the real finite Hamiltonian
algebra before emission; seed cast losses and emitted residuals are retained.
It can contain finite-start transients and has its own mode ID. All nine
attempts reuse the same emitted seed, with L=48/96/192 and h/h2/h4. The finite
tail uses k X_(L-1)-(L+1)X_L/eta-K X_L. A two-stage order-three Radau IIA solve
eliminates the two 2x2 tail chains into a pivoted 22x22 core, recovers both
tails and checks the original full-stage equations. There is no dense
hierarchy fallback. Clock stages use a separately bounded fixed-point solve.

For w=exp(-tau), g=Kw and tau'=-K, the complete raw split channels are
T0=w phi'+g delta_gamma/4, T1=w k psi+g theta_b/k, T2=g Pi/8 and
Psrc=sqrt(6)g Pi/8, all in Mpc^-1 per the same signed mode. The outward-sky
source T0+i mu T1-P2(mu)T2 uses the existing continuous projection and positive
scalar E convention. No integration by parts drops finite endpoints.
Quarter-cell binary64 output times preserve every original knot and are the
actual trajectory endpoint law. Complete outputs use eta-major/k-minor order.

The shared boundary retains tau_i and strictly positive w_i, monopole
delta_gamma_i/4 and dipole theta_gamma_i/k. Its temperature contribution is
w_i[monopole j_l(x_i)+dipole j_l'(x_i)]. Initial E is exactly zero because of
the explicit seed, rather than a generic adiabatic assumption. The continuous
operator computes only the forced support integral. This producer does not
encode a delta cell or silently add/renormalize the boundary. It retains the
positive omitted-boundary bound w_i(|monopole|+|dipole|). Complete transfer
composition needs a producer-owned radial/error receipt.

Hard caps are two k, 1024 original times, 4096 fine times, 100000 aggregate
attempted steps, 1000000 background/clock calls, 200000 **whole coupled-stage
solve attempts**, two billion charged numeric destination updates and 64 MiB
simultaneous owned payload. Actual 2x2 inversions and core factorizations have
separate counters. Mandatory refusal-prefix writes reserve ledger space before
a solve and are recorded separately until executed. Failed preparation work
is never reset by evaluation. Capacity, strings with terminators, original
owners, retained attempts, source arrays and model scratch are counted;
requested replacement buffers are checked before reserve and actual capacity
after reserve. Allocator metadata, RSS and arbitrary caller clones are excluded.
Begun step counts increment locally only after their aggregate reservation;
denied requests have a separate counter and perform no next-step computation.
The aggregate begun/denied counts equal their sums over all attempt receipts.
The logical record ledger charges a node's 33 initialized fields, 33 population
updates and 33 owned-copy fields separately. Each channel diagnostic charges
12 initialized fields (including three absent-error discriminators), nine
population updates and 12 owned-copy fields. Reservations precede construction,
assignments and publication. These are declared buffer allowances, rather than
compiler-store counts or physical-input uncertainties. Private prefix controls
cover denials before local population and before owned publication.
Both core and tail stage checks reject nonfinite derivatives, terms, residuals,
scales and normalized values before maximum aggregation; the prior complete
state is retained on refusal.

Complete raw output currently returns `conditioning_budget_exceeded` with
`source_numerically_admitted=false`. Time/hierarchy differences and measured
cast losses are witnesses. Integrated correlated background/clock,
arithmetic/linear forward response and source-grid errors remain absent.
The original allocations are 1e-8+3e-5|X| for dimensionless state and
1e-11+3e-5|S| Mpc^-1 for channels; each of the five required contributions
must occupy at most one fifth. A complete independently authored angular/time/
precision reference must occupy at most five percent of **each** comparison
allocation. Endpoint agreement, shared ancestry, a build or CI cannot earn
those gates. Neither the downstream projection's missing radial/source-grid
errors nor the positive boundary are assigned zero.

The native controls separately assemble a small dense stage algorithm, check
the finite Hamiltonian/zeta seed, exchange cancellation, zero collision energy,
polarized -0.3K fast decay, stationary leading 16/45 shear, masks, lifetime,
one-move export and refusal accounting. These are engineering/analytic
controls; the dense small hierarchy is not the independent angular reference.
The installed public-header/archive consumer also produces one bounded positive
finite source and checks its unadmitted discriminator and surviving boundary.

At frozen commit `04332db2239122ea2a2685e2375ebb15e98686ea`, the first bounded
trial passed both native contract binaries and the fresh installed consumer.
Its seven serial configure/build/run/install steps all passed in 67.410 seconds,
including fresh setup and compilation. Actual strict C++20 flags, installed
archive/header identity and all 89 unchanged source/runtime input pins were
checked; an independent receipt readback also passed. This passage covers the
named narrow trial. Full repository and integrated-main CI checks are separate.
The compiled controls checked the work/payload caps and the one-begun,
one-denied between-node refusal. Full-run internal ledger values were not
printed, so this passage supplies no measured internal work-count claim.
Numerical admission, the independent angular/precision reference and physical
qualification remain unearned; complete raw output retains the unadmitted
discriminator above.

The equations were independently authored from retained MB, CLASS and
line-of-sight source readings. No CLASS code/assets or paper bodies are copied
into production. [The source pin manifest](finite-opacity-source-provenance.json)
retains all 49 immutable proposal members, exact baseline and the original
pre-code cap finding. Original paper/code body copies and licensing receipts
remain in ignored evidence. Redistribution gaps would need resolution before
copying those bodies; no such copying is required by this implementation.
