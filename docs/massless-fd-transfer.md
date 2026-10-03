# Pure massless FD collisionless transfer

This branch contains a source implementation of a distinct native C++ consumer:
flat GR, positive pressureless CDM, nonnegative Lambda and explicitly supplied
zero-mass Fermi–Dirac species. Strict GNU Release compilation and four narrow
native/shared-consumer tests passed at checkpoint `4e52cab4`. The current peer
repairs, full numerical comparison and actual installed caller remain pending.
This page states the bounded candidate contract; it does not award a numerical
or observational qualification.

The physical identity is
`GR/flat-pure-explicit-massless-FD-CDM-lambda-unit-zeta/v1`. It excludes physical
photons, baryons, additional massless density and positive-mass relics. It does
not provide ordinary full-species matter transfer, a primordial spectrum,
full-support sigma8, recombination, visibility or a primary CMB spectrum. The
existing [perfect-fluid transfer](linear-transfer.md) has a separate identity
and radiation closure. NEXT-14 and NEXT-19 remain open in the sole [roadmap](roadmap.md).

## State and source conventions

`prepare_massless_fd_transfer(background, initial_a)` retains the actual
[thermal background](thermal-neutrino.md), source species order, temperatures,
statistical weights, momentum method and densities. The source requires
`omega_gamma = omega_massless_nonphoton = omega_b = 0`, `omega_cdm > 0` and
one to sixteen explicit species with exactly zero mass, positive temperature
and positive weight. Equal massless brightness modes may be aggregated because
they share the declared adiabatic initial mode; their physical source entries
are retained. A photon-temperature mapper must not substitute a different
background for this source.

The selected original equations are Ma & Bertschinger's conformal Newtonian
metric and Einstein equations (23), massless collisionless transport (48),
hierarchy (50), approximate outgoing endpoint (51) and leading mode (98), from
[the original paper](https://arxiv.org/abs/astro-ph/9506072). The implementation
and centered reference are independently derived algebra, not copied reference
engine code. The source review binds those selected claims; it does not promote
all equations or physics in the paper to implemented capability.

The conventions are
`ds² = a²[-(1+2 psi)d eta² + (1-2 phi)dx²]`, spatial `exp(+i k x)` and
`F(mu) = sum_l (-i)^l (2l+1) F_l P_l(mu)`.
`delta_r = F0`, `theta_r = 3 k F1/4`, `sigma_r = F2/2`.
Conformal age is in Mpc, k and conformal H are in inverse Mpc; the retained H0
and speed of light supply the conversion. Outputs are signed fields per unit
asymptotic curvature zeta. The leading pure-radiation mode has
`psi = -10/19`, `phi = -14/19`, `delta_c = 15/19` and
`sigma_r = psi (k eta)²/15`. Newtonian `delta_c` differs from the requested
comoving `Delta_c = delta_c + 3 Vc`, with `Vc = Hcal theta_c/k²`.

The native state uses `Delta_c`, relative entropy `S = delta_c - 3 delta_r/4`,
`dV = Vr - Vc`, `Vc`, `phi`, scaled shear `sigma_r/(k/Hcal)²`, conformal age
and F3 through FL. The same retained `P = a⁴ E²` supplies H and fractions;
fractions are not renormalized to close a rounded background. Closure and the
background identity defect are reported separately. Their inherited precision
and the existing binary64 Lambda accessor are shared ancestry, not independent
background evidence.

The initializer is the explicitly projected leading approximation. It keeps
leading phi and scaled shear, sets S and dV to zero and projects Delta and Vc
onto the Hamiltonian constraint. In pure radiation its leading comoving
coefficient is `Delta_c = (7/19) zeta (k eta_i)²`. The regular growing series
has coefficient `1/4` and a nonzero higher-order dV; this implementation does
not silently substitute that series or renormalize zeta at finite start.
Leading and projected primitive states, their corrections and constraint
residuals remain in each attempt witness.

## Bounded native request

The ordered k batch is in `[1e-7, 1e-2] Mpc^-1`, common output a in
`[1e-4, 1]`, and the upper conformal phase must not exceed 20. Initial a is in
`[4e-20, 1e-6]`; all three native starts require matter/radiation ratio at most
`1e-8`, initial phase at most `1e-4`, Lambda fraction at most `1e-12` and
initial a below output a/16. The early age uses the radiation/CDM integral and
an explicit Lambda tail, radiation and arithmetic allowance. Independent
adaptive Simpson integration of the same retained P checks endpoint age.

`evaluate(k, a, mask, policy)` requests any nonempty subset of
`massless_fd_comoving_cdm`, `massless_fd_spatial_potential` and
`massless_fd_lapse_potential`. Order, duplicate modes and invalid original rows
are retained. Unsupported input, exhausted work/storage or an unearned error
share withholds the requested value. A batch status does not erase per-row
refusals. Copy/move ownership remains ordinary native library ownership; no
new CLI, ABI, runtime expression or cache service is introduced.

Each field has epsilon `1e-7 + 1e-4 |X|`. Five native empirical allocations
(time, finite start, hierarchy, background/age and arithmetic/cast) each need
at most epsilon/6 and together at most 5 epsilon/6. They are refinement and
sensitivity estimates for the bounded request, not rigorous global forward
error certificates. Accuracy can be tightened, but original tolerances and
work caps cannot be widened through the policy.

The eight native anchors are finest `(a_i/4,L128,h/4)`, starts a_i and a_i/2
at finest settings, h and h/2 at finest start/L, L32 and L64 at finest
start/time, and joint coarse `(a_i,L32,h)`. All levels use the same L128-based
complete mesh. The conservative mesh uses `min(.04,.04 exp(-.04)/x,Hcal eta_lower/130)`;
every attempted RK stage separately checks `129 h <= Hcal eta_lower` and
`h x <= .04`. Each combined step also checks its actual conformal-age
increment plus wide arithmetic allowance against `.04`. No
unrecorded retry resets a failed stage. The eight numerical-input trials are
plus/minus radiation normalization, the shared P envelope, the Lambda getter
cast and initial age, all at the matched finest settings. The final arithmetic
control rounds only the six core perturbation components to binary64 at RK
stage/final assembly; high hierarchy moments and age retain wide arithmetic.
Every attempted active RHS and stage assignment counts. Reserved State padding
and witness copies are not hierarchy updates.

The original limits are 16 k rows, one million RHS attempts and 128 million
active scalar updates per row, four million RHS and 512 million updates per
batch, four million background queries, 200,000 age quadrature evaluations
and 16 MiB native payload. Failed stages/controls consume their completed
work. The first cap refusal is evidence; a later method amendment cannot erase
it. The native 17-attempt witness retains status, initializer, state/metric
epochs, fields, counters, phase maxima and constraint maximum. A background
query refusal withholds the requested result. Its diagnostic lapse witness
remains available precisely when the retained successful metric epoch matches
the actual state epoch; otherwise `lapse_available=false` withholds psi.
An extra query failure may leave a valid same-epoch cached witness. Delta, phi,
scaled shear and age survive independently of lapse availability.

## Independent reference and actual caller

The test-only reference transports a centered angular representation of the
original MB48 brightness, uses actual stored Gauss–Legendre moments and Gram
translations, and evolves the metric through the original Einstein TRACE
rather than the production momentum equation. It uses conformal-time adaptive
DP54, 64/128/256 angular nodes and earlier starts. Its numerical off-constraint
TRACE law adds `-kappa E`, with kappa 2, plus 1.5/3 controls. Exact constrained
solutions are unchanged. The homogeneous constraint stability proof covers
Einstein E/W modes, not unwanted physical modes or every angular discretization
error. Original E, W, TRACE, correction, slip/conservation/background forcing
and recentering histories remain recorded.

Source review found incomplete forcing-to-output accounting in the first peer
implementation. Its normalized-conservation integral and endpoint output scale
cannot establish an output error allowance. Complete state/event response and
retained-background derivative accounting remain blocking method work. The
source records preserve that original proposal and the subsequent repairs;
their presence does not qualify a peer result.

Time/kappa, angle, start, background/age and arithmetic/recentering reference
shares each need epsilon/30, total epsilon/6. Its independent controls include
closed characteristic/Bessel transport, actual-moment identities, constraint
propagation/stability, Gram translations, radiation and EdS mathematical
limits, species splitting/order and adversarial refusals. Shared source and
thermal background ancestry are explicitly retained. Every attempted/rejected
stage record is byte-guarded: 256 MiB per k case, 2 GiB campaign and 8 KiB per
record, alongside two million RHS, 512 million updates and 16 MiB payload.
No raw long-double padding is serialized. CTest's 300-second timeout applies
to the whole peer test binary, including all its cases; runtime limits are not
performance qualification.

`cpp/tests/test_installed_massless_fd.cpp` is the actual small SDK caller:
H0 70 km/s/Mpc, fractional CDM .3, zero photon/additional/baryon densities,
one species `(mass=0 eV, T_today=.0002 eV, weight=2)`, initial a `1e-14`,
output a `1e-4`, ordered k `[1e-7,.01]`, and all three fields. It compiles
against independently installed public headers and `libirred_core.a` through
`tools/check_install.py`, the repository's current standalone SDK route.
No CMake package or `find_package` route is implied. Native GCC/Clang controls
also compile that caller. A low-k Delta below the absolute allocation has no
claimed relative accuracy. None of these synthetic controls is an observation
or supplies a posterior or physical primordial amplitude.
