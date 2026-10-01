# Roadmap

Irreducible is intended to be a simulation-first, agent-operated physics engine for cosmology and the astronomical observations that constrain it. One owned implementation of an equation should serve synthetic experiments, specific objects, observable prediction and later fitting. The existing engine is a useful numerical starting point; it is not yet that complete forward model. [Capabilities](capabilities.md) describe implemented calculations, and [gaps](gaps.md) separate present data and software from proposals.

This is the active development plan. Historical calculations and comparison records are evidence, not competing instructions. Advance a capability when its own physical, data and numerical prerequisites are ready; do not require completion of whole fields or a general inference ecosystem.

## Shared physics, explicit observational operators

The forward direction is:

**Theory and initial/boundary conditions → physical state → source/object population → propagation → instrument and selection → predicted observables → probability model.**

The acquisition direction preserves original observations and reductions so that a likelihood can compare corresponding predicted and measured objects. Images, spectra, light curves, catalogue statistics and timing are different observables of shared physical states. They should reuse geometry, radiation, gravity, units and numerical kernels where their assumptions agree. They must not silently share incompatible closures, frames or calibration conventions.

A theory implementation declares its equations, matter content, geometry, observers, degrees of freedom, closure relations, domain, initial conditions and supported observables. A background H(z) does not supply density perturbations, metric potentials, recombination, source populations or instrument response. For example, the [CLASS perturbation module](https://github.com/lesgourg/class_public/blob/master/source/perturbations.c) depends on background and thermodynamic state as distinct inputs; that architecture is a reference, not an automatic qualification of this engine.

A reusable computation needs a concrete scientific owner, not a universal physics language. Agents add ordinary compiled source, rebuild and test. There is no runtime expression interpreter, plugin framework or server prerequisite. Empirical relations and physical simulations can coexist if each is labelled: an empirical Cepheid period–luminosity–metallicity relation is a useful conditional model without claiming to solve stellar evolution.

A concrete ownership proposal is a cosmology module for theory-specific background, thermal and perturbation state; reusable spectra and image operations for source radiation, convolution and measurement; and observational-analysis consumers that compose those with instrument, selection and probability contracts. Common units, frames, numerics and instrument response have one owner. Add compiled types/functions only when a real consumer exists: this is an ownership boundary, not a promise of empty namespaces or folders.

The common numerical/data layer should support units and frames, stable linear algebra, quadrature, differential equations, interpolation with stated errors, probability, spectra, images and uncertainty propagation. Add each operation for a real consumer and verify its contract. Rendering and interactive visualization belong in external tools and are the lowest implementation priority.

## Dependency graph and first executable verticals

The forward model advances through small complete calculations. [Deterministic photometry](photometry.md), its native sampled extension and [the supplied-drag sound horizon](sound-horizon.md) implement bounded first pieces; their named numerical checks do not complete the verticals below. Calibrated-response uncertainty, noise and recovery, a joint ladder, and predicted thermal/drag history remain proposals, with no claim that those comparisons or new data acquisition have run. A scoped joint fit can start as soon as the forward models and dependence contracts for that subset qualify; it need not wait for CMB or nonlinear simulations.

### 1. Synthetic photometric observation and recovery

**Bounded first increment:** deterministic radiometry for a finite constant rest-frame wavelength spectrum, a rectangular observed passband and optical transmission. `photometry.predict` supplies incident flux, collected energy and expected transmitted photons under a supplied distance and explicit propagation convention; see [photometry](photometry.md). This is not a noisy experiment or recovery model.

**Native sampled increment:** the standalone C++ photometry API integrates finite piecewise-linear wavelength spectra and optical passbands with explicit support, normalization and fixed-calibration assumptions; see [photometry](photometry.md). Its constant/linear controls and independent frequency-coordinate refinement concern the declared sampled model. CLI ingestion and calibration-uncertainty propagation remain next prerequisites. Noise, recorded random streams, selection and recovery need their own scope and checks.

**Inputs:** an explicitly identified source spectrum or empirical time-dependent SN spectral template, redshift/time convention, distance or flux normalization, extinction hypothesis, measured passband and detector response, observation times, noise and selection parameters. Begin with one passband and a simple analytic source before a multi-band SN.

**Calculation:** transform source radiation into observer-frame spectral flux, integrate the declared photon- or energy-weighted passband response, then predict counts or calibrated flux. A later image variant adds PSF convolution, pixel integration, detector effects and sky/background. Generate noise with a recorded random stream; model detection/censoring separately from the incident signal. Spectrum, attenuation, geometry and instrument effects remain distinct.

**Outputs:** noiseless predictions, simulated measured values, uncertainties, detection decisions and an immutable truth/assumption record; a reduction/recovery result refers to that same experiment. Use signed flux rather than replacing non-detections with positive magnitudes.

**Gate:** analytic constant-spectrum/passband limits, dimensional and energy/photon-count checks, redshift/time mapping, independent integration, resolution refinement and recovery of known synthetic parameters. Challenge selection and non-detection handling. A trained template comparison establishes template implementation agreement, not a physical explanation of SN luminosity. Released SN training and calibration are inseparable from their declared assumptions; see the [Pantheon+ calibration analysis](https://arxiv.org/abs/2112.03864).

This supplies a practical simulation-first foundation for spectra, photometry and eventual raw-image comparisons without waiting for a full stellar or hydrodynamic model.

### 2. Conditional absolute calibration and H₀

**Native executable increment:** the [synthetic calibration ladder](calibration-ladder.md) now predicts and jointly recovers host moduli, empirical Cepheid coefficients, SN luminosity, one shared zero point and H₀ under an explicitly supplied fixed reference distance shape. Two anchors and a distinct calibration measurement give an identifiable synthetic control with ordered covariance. Independent KKT, direct/compressed relative-score and remove/refit held-out anchor/Cepheid tests establish named numerical recovery. Released-data fitting, uncertainty qualification and a normalized parameter posterior remain separate work.

**Linear prerequisite:** the standalone C++ [Gaussian design profile](gaussian-design.md) retains an ordered multi-column design and covariance factor, declares shared nuisance identities and uses conservative rank/conditioning admission. The native ladder consumes this prerequisite; the scalar CLI offset remains a separate relative shape calculation. Normalized multi-column proper-prior marginalization, released-data fitting and observational qualification remain missing contracts.

**Inputs:** named anchor likelihoods, Cepheid periods, fluxes/colours, metallicities, host and instrument identities, covariance/calibration responses, selection/crowding/extinction declarations, calibrator SN observations and Hubble-flow SN observations. Every duplicated object or calibrator has one identity.

**Calculation:** an explicit empirical Cepheid relation, such as absolute magnitude = intercept + period term + metallicity term, connects anchor distances to host distances. Host SN calibration connects those distances to a named SN luminosity model; the Hubble-flow forward prediction then depends on the specified geometry and H₀ definition. Wesenheit or extinction choices, relation breaks, parallax offsets, peculiar velocities and population/selection terms are model components, not preprocessing truths. The [SH0ES measurement paper](https://arxiv.org/abs/2112.04510) gives a concrete released ladder and sensitivity variants to reconstruct conditionally.

**Outputs:** predicted anchor/Cepheid/SN observations, correlated calibration parameters, residuals and a joint likelihood or explicitly relative profile. Keep empirical calibration separate from a physical stellar model. The existing free-offset SN shape profile has no absolute H₀ information; inserting H₀ into a background query does not measure it.

**Gate:** synthetic ladder recovery, exact linear Gaussian limits, covariance and shared-nuisance controls, direct versus compressed likelihood comparison, source/host linkage and held-out anchor/host prediction. A historical correlated host-distance compression may be a reference fixture only under its released selection, relation and nuisance measure. Do not multiply it into a likelihood that already used those calibration observations.

Once these prerequisites pass, their scoped joint calibration/H₀ fit is useful immediately. Additional independent routes—TRGB, masers, surface-brightness fluctuations, sirens and time-delay lenses—join only with their own forward contracts and actual data.

### 3. Early-time state and a bounded sound-horizon calculation

**Bounded first increment:** `cosmology.sound_horizon` implements an early-time expansion and photon–baryon sound-speed state, followed by the sound-horizon integral at an explicitly supplied drag redshift; see [the contract](sound-horizon.md). This is a conditional ruler calculation, not a prediction of the drag epoch. Full linear perturbations and CMB spectra are later consumers, not prerequisites for this integral.

**Inputs:** a defined baseline expansion/matter/radiation content, photon and baryon densities, units and early-time domain; a declared drag redshift and its origin. Use the tightly coupled photon–baryon sound-speed convention and specify how the high-redshift tail is bounded. The input drag epoch is part of scientific identity and cannot inherit qualification from an unrelated cosmology.

**Calculation:** compute r_s(z_d) = integral from z_d to infinity of c_s(z)/H(z) dz, with c_s = c/sqrt(3(1+R)) and R = 3ρ_b/(4ρ_γ) under that approximation. The first implementation supplies z_d explicitly. A subsequent increment predicts the drag epoch using an identified ionization/recombination and baryon-drag optical-depth contract, with atomic assets and a thermal history. Distinguish the last-scattering sound horizon from the drag ruler and from an empirically free BAO scale. The [CLASS thermodynamics implementation](https://github.com/lesgourg/class_public/blob/master/source/thermodynamics.c) is a primary reference for this dependency, not an imported authority or production requirement.

**Outputs:** the conditional ruler, integrand/state diagnostics, declared compactification/tail treatment, empirical numerical error estimates and independently allocated checks, and the origin/status of z_d. Later thermal work returns its predicted drag epoch separately. The current free H₀r_d parameter does not identify H₀ and r_d independently.

**Gate:** dimensions and scaling, analytic simplified expansion/sound-speed limits, independent quadrature, tail treatment and refinement. A downstream BAO observable/log-density allocation is a separate future gate. When z_d becomes predicted, add independent thermal/drag comparisons and propagate its uncertainty. Neither this integral nor a matching ruler qualifies growth, lensing or CMB predictions.

Before using this ruler in a physical-ruler BAO likelihood, connect early and late expansion through one consistent matter/radiation, geometry, unit and parameter identity. The current free ruler remains a separate conditional calculation. Predicting the drag epoch then requires its thermal/ionization contract.

Background plus thermal/ionization history can subsequently support a declared perturbation closure, primordial modes, transfer functions and line-of-sight projection. Each later observable needs its own reference and likelihood allocation.

These three verticals share foundations and can progress independently where inputs are available. Their outcomes determine the next useful consumer rather than creating empty modules in advance.

## Cosmological and object-specific branches

| Branch | Required forward physics and observational inputs | Verification before inference |
| --- | --- | --- |
| SN spectra/light curves and populations | Empirical or physical source evolution, radiation/dust, passbands, calibration, host populations, redshift/time mapping and detection selection | Synthetic recovery, spectral/temporal refinement, independent flux projection and held-out events; distinguish fitted distances from calibrated flux or pixels |
| BAO, clustering and RSD | Early ruler or explicit free-ruler hypothesis; initial spectrum, transfer/growth, velocities, bias, reconstruction, survey window/mask, scale cuts and nonlinear prescription | Analytic/linear limits, estimator and window controls, mock recovery and validity of compression under the tested theory |
| CMB temperature, polarization and lensing | Thermal/recombination/reionization and perturbation closure, primordial modes, line-of-sight transfer, lensing, foreground/instrument model and covariance | Independent spectra, precision/refinement at likelihood level, estimator response and shared primary/lensing dependence |
| A specific strong lens | Lens/source/observer geometry, mass and line-of-sight structure, source brightness/variability, image PSF/noise, time delays and stellar kinematics | Synthetic image/time-delay recovery, degeneracy controls, kinematic/anisotropy sensitivity and actual system data; cosmological distances do not determine the lens mass |
| Nonlinear structure and light cones | Theory-consistent initial conditions, gravity/field constraints, particle/mesh or fluid evolution, boundary conditions, baryonic/radiative closures and ray tracing | Conservation, force/solver controls, resolution/convergence, reproducible initial phases and projection errors; subgrid agreement is not physical truth |

DESI's [released clustering products](https://data.desi.lbl.gov/doc/releases/dr1/vac/full-shape-bao-clustering/) contain window and covariance matrices, mocks and likelihood scale-cut/systematic processing. Those are conditional reduced products, not theory-neutral galaxy positions. An alternative theory may require redoing coordinate mapping, templates, reconstruction or a lower-level likelihood. For strong lenses, [mass-sheet transformations](https://arxiv.org/abs/2011.06002) make mass, kinematic and line-of-sight assumptions essential to time-delay distance inference.

Weak/cosmic shear needs source-redshift distributions, galaxy shape and PSF measurement, intrinsic alignments, matter/metric potentials and projection, with survey selection and covariance shared with LSS. It is a separate observational operator from a specific strong lens or CMB lensing. Redshift drift and peculiar velocities can become small named consumers once their observer/time and velocity conventions are available.

Propagation contracts declare photon conservation, reciprocity/distance duality and redshift/time mapping before using a flux law. FLRW luminosity-distance relations must not be applied automatically to another theory. Detector projection must conserve the declared flux/count quantity and apply effective area and photon-energy conversion exactly once; foreground and source-frame attenuation remain distinct.

Relativity must be represented at the level required by the observable: metric geometry and null propagation for light, gravitational and observer clock conventions for timing, and field/constraint evolution for dynamical spacetime where needed. Weak-field or multi-plane approximations are useful named models, not substitutes silently applied to all theories.

Additional modalities enter the same dependency graph:

- **TRGB and surface-brightness fluctuations:** population/colour/extinction, crowding/completeness, luminosity calibration and actual imaging/catalogue assets.
- **Masers:** geometric/kinematic disk model, warp and accelerations, angular/velocity measurements and calibration covariance.
- **Standard sirens:** waveform and detector calibration/noise, orientation/distance, host/redshift or galaxy-catalog selection and detection probability. GW and electromagnetic propagation may differ in an alternative theory. [GW170817](https://arxiv.org/abs/1710.05835) is a concrete primary reference, not an assumption-free distance.
- **Chronometers:** stellar population/SFH/metallicity modelling, differential ages and progenitor selection; an age-labelled table is not a direct H(z) measurement.
- **BBN:** early expansion and thermal state, nuclear reaction network/rates and uncertainty, abundance extraction and astrophysical systematics.
- **21-cm observations:** neutral fraction, spin/kinetic temperatures, source/radiative backgrounds, foregrounds and instrument frequency/beam response.
- **Clusters:** halo abundance/growth, gas/feedback, mass–observable relation and selection, with lensing or other mass-calibration dependence.

Each branch needs identified available products; a literature topic or file extension is not a dataset. [Gaps](gaps.md) records what has actually been located.

## Alternative theories and agent-driven paper tests

A paper implementation begins with faithful transcription of its original equations and definitions. Any repaired closure, changed convention or extension is a separately identified variant. Check dimensions, conservation, limiting cases, stability and closure before fitting. State unsupported observables rather than filling them with a GR approximation without a distinct model identity.

Timescape requires its averaged geometry, clock lapse, bare/dressed quantities and observer/redshift/distance mapping. It is not a new H(z) passed into flat-FLRW projection. Its original [clock/averaging formulation](https://arxiv.org/abs/gr-qc/0702082) and [observable predictions](https://arxiv.org/abs/0909.0749) are specifications to inspect; author viability claims are not validation here. A missing perturbation or nonlinear closure limits the supported test. H₀ labels must state which observer-defined quantity they mean.

Modified gravity declares action/field equations or the limited phenomenological closure, couplings/frame, degrees of freedom, metric slip/forces, initial modes and screening where applicable. The [Bellini–Sawicki scalar-tensor formulation](https://arxiv.org/abs/1404.3713) illustrates additional linear functions beyond background expansion within its stated class; it is not a universal replacement for arbitrary theories.

An agent may systematically test papers, but it must:

1. Register the exact model, supported predictions and preprocessing dependencies.
2. Declare compatible data, nuisance assumptions, priors and overlap/cross-covariance before comparison.
3. Verify numerical and synthetic cases, then held-out or posterior-predictive cases appropriate to the claim.
4. Separate incomplete, not computable, unsupported, numerically failed, data-incompatible and validated conditional results.
5. Report complexity, prior sensitivity and search multiplicity. The lowest χ² among many tried models is not evidence of truth; searching many papers creates a selection problem.

Inference declares its target measure and estimand: profile scores, relative targets, proper-prior marginal densities and evidence are different products. MAP depends on coordinates; an improper offset integral cannot supply absolute evidence. Check identifiability, domain extent, consequential tails and support overlap before using reweighting or reporting moments. A nonexistent moment or missed mode cannot be repaired by a convergence diagnostic. Numerical failure remains distinct from mathematically established zero probability.

Do not automate a paper's credibility from its popularity or outsider status. Evaluate computable claims and documented assumptions using the same standards for this engine and references.

## Execution and extension design

A future bounded run recipe would declare pinned data sources, compiled models, typed steps and resources in one strict JSON graph. Start with structural resolution and current requests, then narrow acquisition and retained-data reuse. Validate identities and dependencies before scientific execution. Recipes record changed variants and conditional reproduction status; they do not certify papers. The [run-recipe design](run-recipes.md) explains the proposed contracts and format decision. No recipe runner is implemented.

Rust owns bounded acquisition, structural decoding, source identities and execution records. The current retained stream is serial; parallel acquisition and larger streamed I/O are proposed execution work, not current performance claims. C++ owns the shared numerical/physical inner loops. Reuse an immutable source and retained preparation across batches. Assess established format libraries with representative data, correctness, allocation and throughput measurements before writing a custom parser. Parallel acquisition must not overwhelm memory or compete invisibly with compute.

Design the data ownership for GPU and distributed workloads now: explicit host/device transfers, bounded memory arenas, batch/shard ownership, halos or reductions, checkpoint schema, RNG streams, restart semantics and numerical precision. Implement CUDA or MPI/distributed paths when an actual qualified workload and measured bottleneck justify them. Preserve a portable CPU baseline and compare total setup, transfer, compute and memory at matched observable/likelihood accuracy. Thread count and vector width alone are not acceleration evidence.

Randomness belongs to the experiment: record generator/stream assignment and seed or counter, including sharding and restart. Distinguish deterministic reproducibility from acceptable statistical equivalence. Checkpoint publication must distinguish an atomic manifest from power-loss durability, and reject incomplete shard generations. Stable work identities determine RNG streams rather than arrival order; a changed rank/device topology needs a tested restart contract. Checkpoints bind physical state, assets, model and numerical configuration; restarting a different model cannot inherit the old identity. The current four-job development coordination budget is a local working rule, not a product-scale limit.

A milestone is complete only when a useful path is implemented, its physical/data scope and error allocation are explicit, independent gates pass, records describe actual work, and the interface/docs agree. A common kernel defect can affect every modality consistently; shared code reduces duplication, not the need for independent evidence.

## Permanent regression and comparison lifecycle

Error allocation belongs to the downstream observable, likelihood and estimand. Deterministic biases are not independent random errors to combine in quadrature; a two-setting agreement alone does not establish asymptotic convergence. Required outputs are declared before execution and cannot become optional after failure.

Before changing a calculation, declare its equation, domain, consumer error budget and reference ancestry. Qualification combines analytic limits, genuinely independent algorithms or derivations, refinement/high precision and adversarial failure controls. Preserve the first discrepant input/result before changing an expectation; explain whether the cause is a defect, input, convention or intentional model change.

Migrate the smallest sufficient licensed comparison facts and reconstructible derivations into permanent C++ mathematical/scientific tests and Rust boundary/record tests. Routine CI must not require the historical Python tree, an external reference engine or a temporary comparison environment. If an optional oracle is needed for a stronger certificate, its dependency, input identities and unavailable status remain explicit; it does not replace portable regressions.

After the migrated tests pass and their owners release the artifacts, delete only inventoried disposable comparison scripts, environments and duplicate outputs. Preserve originals, source identities, reference uncertainty/derivations, reconstruction facts and immutable failed/accepted scientific receipts. Test migration and cleanup are part of completing a capability, not a reason to erase inconvenient history.
