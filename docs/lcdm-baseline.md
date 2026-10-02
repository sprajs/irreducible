# A bounded LambdaCDM baseline

Use an established cosmology to test the stack before interpreting a new theory.
The reference is a declared physical model, data release and calculation; there
is no single parameter point that must reproduce every published result.
[Planck 2018 VI](https://arxiv.org/abs/1807.06209) reports parameters inferred
under its model and likelihood. Those posterior summaries are useful reference
results, not additional independent observations.

The [Reproducible baseline packet](https://github.com/sprajs/reproducible/tree/main/experiments/lcdm-baseline)
consumes the compiled library. [Prospector](https://github.com/sprajs/prospector)
owns reviewed papers and its versioned candidate. This guide describes the
engine boundary; the [roadmap](roadmap.md) remains the sole development plan.

## Physical identities

The full reference is spatially flat GR with cold dark matter, baryons, photons,
neutrinos, a cosmological constant, specified primordial modes and thermal
history. The Planck base model fixes the neutrino mass sum to 0.06 eV. Its
massive-neutrino energy density evolves between relativistic and nonrelativistic
regimes. Neither a constant radiation fraction nor pressureless matter at every
epoch implements that transition. A complete parameter-table reference also
names the likelihood/data combination and distinguishes physical densities
omega_i = Omega_i h^2 from fractional densities Omega_i.

The runnable native variant has pressureless matter, massless radiation and
Lambda, with

    E(z)^2 = Omega_r (1+z)^4 + Omega_m (1+z)^3 + 1-Omega_m-Omega_r.

It supplies photon and baryon subsets and a drag redshift with an origin. Flat
FLRW, a common observer/expansion redshift and photon conservation give
D_H=c/H, D_M=(c/H0) integral dz/E, D_L=(1+z)D_M and
D_V=(z D_M^2 D_H)^(1/3). The ruler is the tight-coupling sound-horizon integral
at that supplied drag redshift. [Early/late](early-late.md) and
[sound horizon](sound-horizon.md) own the equations, numerical domains and units.

This is an explicitly changed model when compared with full Planck base LCDM.
A supplied Planck drag value does not restore its missing neutrino or thermal
physics. Conversely, the CLI `background.evaluate` LCDM model omits radiation;
it must not silently replace the native early/late state. H0 scaling at fixed
fractions tests units and identifiability: distances and ruler scale as 1/H0,
while their ratios remain invariant. Holding physical densities fixed is a
different parameter variation.

## What the bounded test can establish

| Behavior | Available calculation | Evidence and boundary |
| --- | --- | --- |
| Expansion and flat distances | Native shared early/late state; radiation-free CLI background separately | Analytic limits, independent coordinates and named Astropy background controls; explicit thermal relic E/H and the native thermal distance/ruler consumer now share a tested state with supplied physical densities/temperatures; full-reference likelihood and perturbation closure remain missing |
| Supplied-drag ruler and BAO ratios | Native shared massless or thermal state and conditional BAO density | All released ratio rows and their full ordered covariance; fixed-point numerical comparison, without predicting drag or reproducing a posterior |
| Relative SN distance shape | CLI/native free-offset profile | Named released-input controls; no absolute H0 information or newly combined SN/BAO likelihood |
| Spectral projection and synthetic measurement | Native sampled/temporal photometry, finite shared passband calibration and bounded detector/censoring | Analytic and independent time/frequency and characteristic-function controls; declared source/distance/calibration/arrival law, without measured-instrument or source-population qualification |
| Absolute calibration | Native synthetic ladder, conditional estimator variance/sampling law and proper correlated calibration density | Named synthetic recovery, exact noise-response and held-out controls; released physical-axis mapping, coverage and observational reconstruction remain separate |
| Pressureless growth | Native radiation-free GR D/f | EdS and independent ODE/refinement controls; no radiation/relic transfer, RSD or shear likelihood |
| Hydrogen thermal prerequisites | Native supplied-state H/He LTE and conditional pure-H history with bounded explicit FD relics | Rational/polynomial/Decimal, resolved RK4 and direct-momentum/coupled Radau/refinement controls; distinct prescribed/evolved matter-temperature identities, finite-endpoint Thomson visibility and truncated drag root, without full thermal/physical-epoch qualification |
| Synthetic thin-lens images and delays | Native SIS point source, explicit mass sheet and Gaussian PSF pixels | Independent potential/Jacobian and pixel integration controls; compatible massless distances and synthetic degeneracies, without measured-system mass or instrument qualification |
| CMB and physical drag prediction | Conditional background/ruler and bounded hydrogen pieces | Helium/multilevel/full-endpoint visibility closure, radiation/metric/massive-neutrino perturbations, primordial modes and line-of-sight predictions remain missing |

The DESI DR2 mean/covariance products are
[released fitted compressions](https://arxiv.org/abs/2503.14738), with estimator,
reconstruction and covariance assumptions. Keep their exact bytes, source
revision, redshifts, observable tags, axes and covariance together. A numerical
reference dataset is "gold" only for those declared identities and comparisons;
it is neither a statement that each observed value equals the model prediction
nor permission to reuse that compression under every theory. Unknown dependence
with another probe remains unknown.

The limited comparison allocates 1e-9 Mpc + 2e-11 relative to distances,
2e-11 relative to dimensionless E, 1e-11 + 5e-11 relative to ratios, and absolute 1e-8 separately to the quadratic,
log determinant, Gaussian normalization and log density. Independent-reference
refinement occupies at most 5% of each allocation. Record stricter producer
policies separately. Preserve original failed attempts before changing a
reference, policy or implementation. These named numerical checks provide no
universal accuracy guarantee or inference qualification.

The first real-data attempt can fail because a requested producer allowance is
below combined dependency diagnostics, even when the underlying quadratures
finish. Ratio admission includes numerator/ruler interval propagation and
binary64 arithmetic floors; it needs its own feasible allowance. Diagnose the
reported status, preserve that attempt and justify a revised producer policy
against the unchanged downstream allocation. Do not relax the comparison,
projection gate or physical covariance to manufacture a pass.

## Recorded diagnostic

The October 2026 baseline uses chosen controls H0=67.4 km/s/Mpc,
Omega_m=0.315, Omega_r=9.2e-5, Omega_b=0.049, Omega_gamma=5.45e-5 and supplied
z_drag=1059. These fractions are independent inputs near familiar reference
values; they are not a Planck posterior draw or a temperature-derived parameter
mapping. The second point doubles H0 at fixed fractions and drag.

The initial producer ratio allowance of 1e-14 absolute + 2e-14 relative refused
the combined diagnostics with `conditioning_budget_exceeded`; those attempts
were preserved. After source review, the ratio producer allowance became
1e-14 absolute + 2e-13 relative, while sound/distance producer policies and all
external comparison, refinement and projection gates stayed unchanged. The
bounded massless variant passed 64 comparisons on all 13 DESI DR2 rows: the
quadratic was 28.070790034811825, with independent discrepancy about 4.2e-14,
and the projection diagnostic was about 1.51e-10 against its 1e-8 allowance.
The packet owns exact source/build/run identities and retained failed attempts.
This fixed-point diagnostic is neither a fitted cosmological result nor a claim
that these observations must equal the prediction.

Later [released-ladder](https://github.com/sprajs/reproducible/tree/main/experiments/released-ladder)
and [SN observer/passband](https://github.com/sprajs/reproducible/tree/main/experiments/sn-observer-passband)
packets retain separate SDK/input identities. They test a fixed-coordinate
linear target and source-defined observer/optical controls, respectively. A
chosen radiation-free LCDM point with H0=70 and Omega_m=0.3 is a diagnostic,
not a fitted Planck result. Calibrated data need not equal a deterministic
prediction, and numerical acceptance does not close frame, selection,
calibration or joint-dependence qualification. Historical receipts stay immutable.

Independent source-to-consumer review also found that the experiment's twin-H0
expansion check reused the distance check's additive Mpc allowance. A retained
negative witness showed that the old harness admitted a dimensionless E shift
outside the declared expansion allocation. The corrected check applies the
original E allocation, with a permanent negative regression and a fresh
committed-source run. This repairs the experiment gate; it changes no engine
equation, physical covariance or declared comparison allocation.

The clean engine at db4765838fc404a489a0115ee69713f2ea6cd2f4 also passed 12
selected existing native owner/peer executables for expansion, conditional BAO,
sampled/calibrated photometry, synthetic ladder and proper correlated
calibration. Those tests retain their named mathematical scope; full Planck,
CMB, growth, lensing and observational ladder results remain outside this run.

## Cross-project handoff

1. **Prospector** pins the primary paper version, source hash, coverage and
   source-defined claim. Its candidate names the full reference and any changed
   runnable variant, required assets, lost scope and blockers.
2. **Reproducible** reviews that candidate, pins its revision/path/hash and exact
   input products, declares parameter mapping and acceptance before execution,
   then runs a small consumer of the installed Irreducible library. It owns
   experiment orchestration, plots, immutable attempts and concise findings.
3. **Irreducible** owns compiled equations, numerics and shared likelihoods.
   A failed experiment becomes a source-linked discrepancy witness or explicit
   closure blocker. Implement the smallest useful repair with independent native
   controls and review affected consumers; do not put replacement physics in the
   experiment controller.
4. Re-run the pinned experiment against a clean rebuilt engine after a repair.
   Publish source changes through each repository's reviewed PR/CI workflow;
   record actual source/build/executable/input identities in each attempt.
   A changed build does not rewrite earlier receipts or promote a former failure.

Assign one integration owner per repository and an explicit cross-project
coordinator. Inspect live status, worktrees, open PRs and built discovery before
assigning work. Give each worker exclusive paths, frozen interfaces, acceptance,
evidence destination and resources. Use the requested exact model; record a
capacity failure and reuse correctly configured workers instead of substituting.
The four-job local compute allowance is shared across chats, agents and
worktrees; allocate jobs explicitly and set Rust test threads explicitly.

The native early/late, thermal-observable and conditional-density interfaces have no CLI
or C ABI route. A bounded experiment-specific installed SDK consumer is a real
execution path, not an implemented general paper runner. Verify its archive,
headers, build manifest, compiler, consumer source and executable; bind the
actual input and output statuses. CLI discovery continues to describe its ten
implemented operations.

Public packets and source tests stay small. Bulk third-party inputs, papers and
full run records stay in their designated ignored storage, with terms reviewed
before redistribution. A published scientific result needs a durable archive;
local receipts and expiring CI artifacts alone do not supply one. Separate
execution, numerical comparison, inference and interpretation in the account.
