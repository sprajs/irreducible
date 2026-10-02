# Capabilities

Use the executable's `describe --json` output for its exact build, ABI revision and request schema. Irreducible has one current interface; obsolete model-specific routes are removed together with their callers. Historical records retain their original identities and meanings.

| Operation | Calculation | Scope |
| --- | --- | --- |
| `fixture.checked_i64_add` | Overflow-checked integer addition | Exact integer contract |
| `quantity.convert` | Typed physical conversions | Explicit units, roles, frames and conventions |
| `numerics.scalar_batch` | Compiled scalar numerical methods | Method-specific arithmetic and domain |
| `cosmology.sound_horizon` | Conditional comoving sound horizon | Supplied drag redshift; flat pressureless matter + massless radiation + Lambda |
| `background.evaluate` | Requested flat-FLRW expansion and projections | LCDM, constant q, CPL or fixed five-bin q |
| `observations.prepare` | Immutable typed source preparation | Structural checks; no probability or lineage upgrade |
| `statistics.gaussian` | Normalized Gaussian density, offset profile or proper latent prior | Explicit ordered residuals and nuisance assumptions |
| `supernova.profile` | Conditional single-offset magnitude profile | Same expansion types, with no source effect or an explicit grey magnitude effect |
| `bao.density` | Conditional normalized free-ruler Gaussian density | Same expansion types and explicit H0rd |
| `photometry.predict` | Deterministic rectangular or sampled spectrum/passband prediction | Supplied distance/redshift; incident flux, collected energy and transmitted photon expectation |

Background requests select expansion, radial, luminosity_shape, clock, physical or kinematics groups. Dimensionless E and clock integrals remain usable when a requested dimensional H or lookback-time projection fails. Observer and physical-scale failures do not erase independent outputs. E-only and DH-only calculations avoid radial integration. Exact redshift bits share nodes within a model batch; repeated model specifications across model rows currently recompute their geometry under the same global work budget.

Observation preparation distinguishes measured quantities, released fitted summaries and synthetic controls. The generic magnitude-covariance profile does not confer Pantheon release identity, calibration, standardization or independence. Source hashes, axis declarations and selection remain explicit. A retained process owns the same immutable observation object through its dependent consumers, even after the parent handle is released; factors are prepared once and reused across evaluations.

Gaussian covariance validation applies to the selected principal block. Precision input requires the full precision factor and inversion before marginal selection: slicing precision would mean conditioning. Ordered residual and response IDs must match. A proper named latent prior retains mean, variance, response and independence declaration. An offset profile is a relative score, not normalized density or evidence.

Native supernova comparisons cover the named original seven-point, CPL four-point, fixed-q six-point and conditional grey twelve-point cases on the pinned 1590-row sample. The grey cases do not reproduce a historical 1820-row analysis. Stable/full historical fixed-q comparisons pass their allocations; five compressed approximations fail the reference allocation and remain withheld. These facts do not qualify a posterior, an arbitrary supplied sample or every point in the model domain.

Native BAO comparisons include the named eleven-point and fixed-q twenty-four-point cases on pinned thirteen-row inputs. The latter use an optional test-only directed GMP/MPFR interval oracle. The unchanged log-density budget is 1e-8, allocated as reference 2e-9, factor/solve 3e-9 and background 5e-9. Earlier conservative failures remain historical evidence; the directed certificate establishes the named cases, with a narrow margin at the highest-residual endpoint. Production has no GMP/MPFR dependency. SN and BAO remain separate conditional consumers; no joint likelihood or source independence is inferred.

Required numerical checks that pass satisfy default `numerical_contract` assurance and exit 0. Explicit `qualified` assurance exits 6 when applicable named evidence is absent. Completed scientific failures exit 2. Interpretation remains unqualified; no arbitrary request receives automatic scientific qualification.

Measured mirrored dense preparation improved one pinned fixture on one compiler/machine; that historical measurement does not establish performance of the consolidated interface. Optimization follows matched-quality measurements, not the presence of a wider type or more threads. Samplers, priors over cosmological parameters, smoothing campaigns and whole-domain certificates remain outside these capabilities.

The [conditional sound horizon](sound-horizon.md) uses an exact compact scale-factor interval, with explicitly supplied photon/baryon fractions and drag-redshift provenance. It predicts neither thermal history nor drag epoch; no existing BAO likelihood qualification is inherited.

See [photometry](photometry.md) for the bounded deterministic projection.

The standalone C++ library also provides a retained multi-column Gaussian design profile. Preparation consumes a Gaussian owner and retains its covariance factor and design calculations across residual evaluations. Ordered row/parameter identities, units and declared shared nuisance coordinates are explicit. This native API has no CLI request or C ABI exposure yet; its relative score is not a normalized density or evidence. See [Gaussian design](gaussian-design.md) for its retained covariance-whitened pivoted QR, conditioning admission and limits.

The standalone C++ library also provides `photometry::evaluate_sampled` for finite piecewise-linear rest-wavelength luminosity and observed optical transmission. It integrates the declared interpolation model, with supplied distance and redshift, to incident band flux, transmitted energy and expected photons. A coarse ABI2 batch and the same `photometry.predict` CLI operation now own pooled curves once and evaluate indexed exposures; source roles/provenance remain explicit. Sample interpolation and fixed calibration are assumptions; their uncertainty is not a numerical error estimate. See the sampled contract in [photometry](photometry.md).

The standalone C++ [synthetic calibration ladder](calibration-ladder.md) predicts anchor moduli, Cepheid period/metallicity relations, calibrator SNe and Hubble-flow SNe with a supplied reference distance shape and one shared calibration coordinate. A retained ordered covariance/design supports a relative joint fit and a checked H0 projection. Exact Gaussian/model provenance, conventions, identities and rank admission are required. Named noiseless, correlated, KKT, direct-versus-compressed and held-out anchor/Cepheid controls establish numerical recovery within their allocations. This native API supplies no CLI/C ABI route, observational H0 measurement, normalized posterior or model evidence.

The standalone C++ [early/late calculation](early-late.md) uses the same flat pressureless-matter, massless-radiation and Lambda parameter identity as the conditional sound horizon. It predicts requested E, H, D_H, D_M, D_L, D_V and conditional D_M/r_s, D_H/r_s and D_V/r_s, evaluating the supplied-drag ruler once per batch. Each requested output retains its own numerical status; successful distances survive ruler failure. Named analytic, independent-coordinate quadrature and Astropy background comparisons are numerical evidence. This provider supplies no CLI/C ABI route, drag prediction or CMB qualification. A distinct [conditional BAO density](bao-conditional.md) now composes its ratios with an ordered retained Gaussian observation object.

The standalone C++ [correlated calibration operator](correlated-calibration.md) computes a normalized observed-residual Gaussian density after integrating a declared proper correlated calibration prior. It retains conditional source and effective covariance factors, explicit parameter/row order, measures and provenance. Independent cofactor and direct two-dimensional prior integration controls test named synthetic cases, including null and rank-deficient responses. This native API adds no CLI/C ABI route, posterior, model evidence or automatic ladder/photometry uncertainty qualification.

The native [finite passband calibration law](photometry-calibration.md) propagates explicitly shared, physically valid optical-transmission states through the sampled operator. It returns multi-band energy/photon expectations and full calibration covariance, with numerical sensitivity reported separately. Required state failures or unresolved spread withhold aggregate moments. This adds no CLI/C ABI route, photon shot noise, source/distance uncertainty or qualification of a measured calibration distribution.

The native [conditional BAO density](bao-conditional.md) reuses the retained full covariance factor and one shared early/late physical state for supplied-drag ruler ratios. Projection sensitivity has a separate admission gate; predictions survive a refused density. Ratios and density are H0 invariant at fixed fractions and drag. This adds no CLI/C ABI route, drag prediction, observational H0 inference, or automatic validity of a released distance compression.

For a source-backed standard-cosmology comparison, use the [bounded LambdaCDM baseline](lcdm-baseline.md). Its native consumer and blocked full-model sectors are separate from CLI discovery and runtime qualification.

The native [Gaussian design](gaussian-design.md) also returns a requested
linear-estimator contrast variance using its retained QR. The
[synthetic ladder](calibration-ladder.md) consumes that variance for an explicit
conditional H0 sampling law, including its biased expectation and ordered
quantiles. These assume fixed design, supplied Gaussian observation covariance
and a declared generating mean. They add no CLI/C ABI operation, parameter
posterior, demonstrated confidence coverage or observational H0 result.

The standalone C++ [thermal-neutrino background](thermal-neutrino.md) evolves
explicit collisionless zero-chemical-potential Fermi–Dirac species between
relativistic and nonrelativistic regimes. It supplies density and pressure in
eV^4 and retained flat E/H batches, separating photons, other massless radiation,
baryons, CDM and relics. Independent high-precision, continuity and matched CLASS
controls test named cases. It has no CLI/C ABI, automatic Neff/temperature/mass
mapping, thermal/ionization history or perturbations. The distinct native
[thermal observables](thermal-observables.md) consumer explicitly maps physical
baryon/CDM/other-massless densities and Kelvin temperatures, then uses this same
retained state for E/H, flat distances and a supplied-drag ruler/ratios. It adds no
CLI/C ABI operation, automatic Neff/mass hierarchy or drag prediction. The native
[thermal conditional BAO density](bao-thermal.md) now composes this same state
with the retained ordered observation covariance and an unchanged 1e-8 density
projection gate. Independent cofactor, polynomial/refinement and matched massive
reference controls supply named numerical evidence. Fixed physical densities
do not imply H0 cancellation. Released compression validity and observational
qualification remain separate; the massless conditional model keeps its identity.

The native [proper Gaussian parameter posterior](gaussian-posterior.md) retains
fixed observation/design factors and an explicitly independent full proper
Gaussian prior. It returns conditional mean, covariance and normalized density
in declared parameter units/order/measure. Analytic rational, independent
cofactor and direct joint integration/refinement controls cover named synthetic
cases, including rank deficiency. This adds no CLI/C ABI route, nonlinear
cosmological posterior, released-box prior or observational qualification.

The native [GR growth](gr-growth.md) consumer retains the radiation-free flat LCDM background for a pressureless growing mode D(a) and f(a), normalized by D/a→1. Exact Einstein–de Sitter and independent ODE/refinement controls test its bounded domain. It adds no CLI/ABI operation, scale-dependent transfer, radiation/relic perturbations, RSD likelihood or observational qualification.

The native [ground-state hydrogen equilibrium](hydrogen-equilibrium.md) operator
solves the supplied-temperature/physical-density Saha model with separately
represented positive ionized and neutral fractions. Exact SI definitions and
sourced atomic central values have distinct identities. Rational, Decimal and
independent logit/refinement controls test the declared dilute pure-H model;
LTE, partition truncation and atomic uncertainty remain physical assumptions.

The native [hydrogen–helium equilibrium](hydrogen-helium-equilibrium.md)
operator closes ground-state H/4He stage fractions with one charge-neutral
electron density. Temperature and both physical nuclei densities are supplied.
Original polynomial and Decimal110/150 controls, trace-species limits and tiny
positive output/refusal tests support its numerical contract. This supplies an
LTE state, not a helium kinetic history or abundance/cosmology mapping.

The retained native [pure-H history](recombination-drag.md) composes that initial
condition with the shared thermal background and source-defined
effective three-level rates. It predicts a non-equilibrium electron fraction
under either prescribed radiation temperature or coupled Compton/adiabatic
matter-temperature evolution, a late-endpoint drag depth and a separate
conditional unit-depth root. The coupled variant also predicts Thomson depth,
scattering rate and per-redshift finite-endpoint visibility with boundary survival.
Independent direct-SI RK4 and coupled Radau/refinement controls test the distinct
solver identities. It adds no CLI/ABI route, helium, multilevel correction,
reionization, physical drag-epoch qualification or perturbation/CMB prediction.
Positive-mass profiles now admit up to three explicit collisionless FD species
with separate model IDs and actual nested momentum work. Original direct-SI
momentum/Radau controls test the bounded masses/temperatures; massless scalar
arithmetic and identities are retained. The four-million default work cap still
refuses expensive profiles; a larger explicit caller cap is needed.
The supplied-drag ruler remains separate.

The native [synthetic SIS thin lens](sis-thin-lens.md) retains a compatible
massless-radiation early/late source and computes two distances once, without
evaluating its supplied-drag ruler. An axial point source, SIS template and
explicit nonnegative mass sheet give image positions, signed magnifications,
fluxes and relative observer delays. A circular Gaussian PSF integrates image
flux into rectangular pixels, including an explicit half-open point limit.
Inherited distance diagnostics, topology margins, separate output masks and
bounded pixel work retain typed failures. Independent potential/Jacobian and
Gaussian-density integration controls test named synthetic recovery and
mass-sheet/H0 degeneracies. The tangent-plane, thin-lens and fixed source/PSF
assumptions remain explicit. This SDK adds no CLI/C ABI route, measured-system
mass or kinematics qualification, lens-mass inference, posterior, arbitrary
multi-plane propagation or CMB lensing prediction.

The native [temporal photometry](temporal-photometry.md) owner retains finite bilinear rest-time/wavelength luminosity grids and fixed optical passbands. Coarse observer exposures return full-interval mean flux, collected energy and expected photons with explicit source epoch and time-coverage status. Analytic and independent time/frequency integration controls test the declared model. It adds no template reader, source/calibration uncertainty, noise, selection, inference or CLI/ABI operation.

The native [detector and censoring](detector-selection.md) operator separately composes expected photons with declared Poisson arrivals, fixed quantum efficiency/background/dark current and Gaussian read noise. Caller-owned seed/stream/sample addresses support exact integer replay and bounded simulation. The same measurement law scores discrete or continuous detections and threshold-only nondetections, with an explicit selected-only measure. Analytic, characteristic-function/refinement and synthetic censored-recovery controls test this bounded model. Source/calibration uncertainty, measured-instrument qualification, saturation/digitization, population recovery and CLI/ABI ingestion remain open.

The native [Gaussian predictive](gaussian-predictive.md) operator conditions a
fixed linear proper-prior model once and retains `A*mu` and one factor of
`R*+A*V*A^T` for a joint future vector. Internally correlated future noise is
supplied and declared independent of training noise/prior; shared parameter
uncertainty induces covariance across predictions. Exact covariance-route
controls, direct prior/noise integration and normalized native mass quadrature
test the synthetic scope. Disjoint ordered row/event IDs, units and measures are
explicit declarations. This adds no observational qualification, nonlinear
prediction, train/future cross-noise covariance or released finite-box inference.
