# Capability and data gaps

This assessment separates the current numerical engine from the proposed [roadmap](roadmap.md). It does not qualify a cosmological interpretation, claim that historical analyses were rerun, or treat a file's presence as a valid likelihood. Use [capabilities](capabilities.md) and actual executable discovery for the implemented interface.

## What exists

The current engine exposes ten operations: exact integer addition, typed physical conversion, scalar numerical methods, immutable observation preparation, requested flat-FLRW background outputs, Gaussian calculations, conditional SN magnitude profiles and conditional free-ruler BAO densities, deterministic rectangular or sampled-passband photometry, and a conditional supplied-drag sound horizon. The late-time backgrounds are LCDM, constant q, CPL and fixed five-bin q. They are bounded CPU calculations with explicit source, arithmetic, numerical and execution contracts. Named wide-arithmetic gates use Linux/GCC long double with at least 64 mantissa bits; hosts where long double equals double are not covered, and Ubuntu CI is a separate engineering check. The early ruler separately assumes flat pressureless matter, massless radiation and Lambda; photon/baryon fractions and drag redshift are supplied, not inferred, and no physical-ruler BAO qualification follows. Gaussian nuisance elimination through the CLI supports a scalar offset. The standalone C++ [multi-column design profile](gaussian-design.md) supplies a bounded linear prerequisite; its relative score is separate from the native synthetic ladder and proper correlated calibration density described below.

No current runtime request receives a named scientific qualification. Passing the default numerical_contract means that required numerical checks passed; it is not inference or interpretation qualification. Named comparison evidence applies to its recorded inputs/build/domain, not automatically to a new request.

Shared source ownership and retained factors support repeated evaluations. Named original-input SN and BAO comparisons and analytic/adversarial controls provide evidence for stated cases. They do not establish all-domain accuracy, physical completeness, a joint fit or a new H₀ measurement. The current SN free-offset profile does not retain absolute calibration information; free H₀r_d BAO does not identify H₀ or r_d independently. A supplied dimensional H₀ in a background projection is an input, not a measured result.

The main missing capability is a coherent forward physical state with qualified observational consumers: source populations/templates and instrument/selection laws; an observationally qualified absolute distance-ladder calibration model; complete thermal/perturbation evolution; measured lens systems; nonlinear initial conditions/dynamics/light cones; and qualified GPU/distributed execution. Bounded sampled/temporal, detector and thermal operators below do not complete those consumers. Current historical performance evidence does not benchmark the consolidated interface or future workloads. A new benchmark must measure representative phases, setup/memory and matched scientific quality.

Engineering remains serial CPU execution with bounded dense matrices and narrow source-specific ASCII adapters. General survey image/spectrum ingestion is absent; the optional FITS codec has test-only synthetic BINTABLE coverage and is unavailable in the product. The native detector has an addressed generator and bounded replay controls; full simulation-workload RNG, distributed restart/sharding and device execution remain unqualified. There is no consolidated-interface throughput benchmark. Native SDK consumers can compose compiled models; adding a CLI route also requires coordinated Rust descriptors, schema, ABI, discovery and records. Centralize structural descriptors where useful; keep equations in compiled scientific owners.

## Standard-model diagnostic gaps

The [bounded LambdaCDM baseline](lcdm-baseline.md) distinguishes full Planck base
LCDM from the current supported approximation. The native shared early/late
state treats radiation as massless and matter as pressureless at every epoch.
The separate native [thermal-neutrino provider](thermal-neutrino.md) now
supplies explicit relic density/pressure and flat E/H with a mass transition,
independent high-precision and matched CLASS controls. The native
[thermal distance/ruler consumer](thermal-observables.md) now maps supplied
physical densities and explicit temperatures/species and retains one state for
distances and the conditional supplied-drag ruler. It preserves the massless
model as a distinct physical identity. The native [thermal BAO density](bao-thermal.md)
now composes these predictions with retained ordered covariance and density-level
projection checks. Released compression validity remains unqualified. A native
conditional pure-H history is described below; full thermal/ionization and
physical drag prediction remain open. There is no implicit Neff/mass hierarchy. The CLI late LCDM background omits radiation and has a
different physical identity; native early/late and conditional BAO currently
require a C++ SDK consumer rather than a CLI/C ABI request.

The historical released-ladder admission failure at source 99ce262 remains
preserved. The native retained whitened pivoted QR at source f844080 now admits
the exact same 3492-row/47-column compact products with both covariance arithmetic
profiles and unchanged budgets. All coefficients and the objective agree with
independent LAPACK QR/SVD references; q=3552.759330295523. This closes the named
numerical admission blocker. The complete physical parameter dictionary, released
calibration assumptions and observational H0 reconstruction
still need separate qualification. Native QR contrast variance and a synthetic
ladder H0 sampling law now supply conditional uncertainty under fixed Gaussian
observation noise and an explicitly assumed generating mean. Named synthetic
checks do not qualify released-data parameter mapping or coverage. The relative
profile is not a normalized posterior. A native [proper Gaussian parameter
posterior](gaussian-posterior.md) now supplies fixed-linear conditioning under
an explicitly independent SPD Gaussian prior; it does not normalize the
released finite-box target or reproduce an observational posterior. Source review now identifies original
column 46 as 5 log10(H0 in km/s/Mpc); its formal contrast variance on the same
3492x47 products agrees with independent SVD/QR. A later [Reproducible source audit](https://github.com/sprajs/reproducible/tree/main/experiments/released-ladder)
now joins all 2150 initial Cepheid rows to the pinned primary table and identifies
all 37 host columns. Period/metallicity and the H0 contrast are corroborated;
six later source-supported anchor/nuisance coordinates are now identified, while
original axis 44 has an unresolved physical label and literal fixed zero support.
The audit preserves the paper/release N1365 count discrepancy and four original
table joins that remain unresolved. Reproducible now pins the 46-dimensional
released box target with axis 44 fixed, preserving all 3492 rows and full covariance;
its unboxed fixed-coordinate profile is distinct from box-supported inference.
Three source constraint-mean sensitivity pairs at fixed design/covariance are
recorded. They do not normalize the box posterior, supply a full systematic
uncertainty or qualify observational held-out predictions.
The native [Gaussian predictive](gaussian-predictive.md) prerequisite now retains
a normalized joint future distribution for fixed synthetic linear responses
under the proper Gaussian parameter prior and declared independent future noise.
It propagates shared parameter covariance without treating the predictive law
as another measurement. Native [addressed Gaussian generation and bounded
recovery campaigns](gaussian-simulation.md) now test distinct fixed-truth and
prior-predictive laws, retaining covariance, random identities and every refusal.
Seeded agreement does not prove RNG independence; predictive preparation still
rebuilds its invariant covariance for each new training vector. Observational held-out qualification still requires
resolved object/calibration/selection dependence, and the released finite-box
target remains a separate measure and normalization problem.
See [Gaussian design](gaussian-design.md).

The native [hydrogen equilibrium](hydrogen-equilibrium.md) model now supplies
separate ground-state pure-H fractions under supplied temperature/density.
The [H/He mixture equilibrium](hydrogen-helium-equilibrium.md) now supplies
shared-electron ground-state LTE fractions for independently supplied nuclei
densities. [Supplied baryon abundance](baryon-abundance.md) now maps explicit
physical density, He4 mass fraction and neutral effective masses into that LTE
consumer. The native [finite joint abundance law](baryon-abundance-law.md)
propagates supplied correlated density/mass/temperature states into full
population moments, with separate numerical diagnostics and truthful required
state refusals. The law itself needs source qualification; G and atomic-asset
uncertainty remain excluded. Abundance prediction, helium kinetic rates and history remain open;
LTE cannot replace the missing helium kinetics. The
[conditional history](recombination-drag.md) advances a bounded effective
three-level model with the shared thermal background and explicit choice of
prescribed Tm=Tr or coupled Compton/adiabatic matter temperature. The coupled
variant supplies Thomson depth, per-redshift scattering rate and unnormalized
finite-endpoint visibility with boundary survival. Both integrate drag optical
depth from an explicit late endpoint. Its unit-depth
root depends on that endpoint and physical approximation; it is not a qualified
cosmological drag epoch. Bounded positive-mass relics now feed the same retained
H source, with direct momentum costs charged. An explicit nested CC momentum route
now has the same consumer allocations and independent comparisons; direct
remains the default, and matched workload evidence governs promotion. No hierarchy/Neff mapping or helium
kinetics is added. Helium, multilevel rates, atomic/model uncertainty,
qualified full optical-depth endpoints and reionization remain prerequisites
before physical-ruler or CMB coupling. Finite-endpoint visibility does not
predict the present-day last-scattering distribution or CMB spectra.

The native [GR growth](gr-growth.md) operator supplies a bounded radiation-free
pressureless growing mode. The [conditional amplitude/RSD consumer](growth-rsd.md)
now adds supplied-reference sigma8/f*sigma8 and a retained full-covariance
synthetic density. Source amplitude and survey-estimator/AP/window validity
remain explicit prerequisites for measured RSD. Matter transfer, radiation/relic perturbations, general lensing potentials and propagation, CMB spectra, stellar/source populations, measured detector/selection laws and survey recovery remain separate gaps. A known cosmology makes these missing
sectors concrete; an H(z) or distance match does not supply them. Use the
[roadmap](roadmap.md) to choose the next qualified prerequisite after the
experiment identifies a blocker.

The [Reproducible observer/passband packet](https://github.com/sprajs/reproducible/tree/main/experiments/sn-observer-passband)
now passes 6093 conditional comparisons on 321 released SDSS redshift pairs/twins
and a pinned historical 910-knot optical response. Coasting and chosen
radiation-free LCDM geometry use independent analytic/refined references.
Unresolved sky/frame/velocity lineage and calibrated PHOT units/time/error law
block an observational SN fit. A fitted 102-coordinate zero-point covariance and
nine systematic template variants do not provide a measured optical-state law;
the 102-versus-105 calibration mapping and joint probability masses remain open.

The native [synthetic SIS slice](sis-thin-lens.md) now supplies one bounded
point-source image/delay and Gaussian-PSF pixel calculation using compatible
massless-radiation distances. Its synthetic mass-sheet/H0 controls demonstrate
remaining degeneracies. Actual system images, PSF/noise, source variability,
mass/environment constraints and kinematics have not been qualified; neither
these synthetic controls nor background distances determine a lens mass.
General or multi-plane propagation and cosmological/CMB lensing remain separate
physics consumers.

## Data products are not interchangeable

The historical research collection contains several levels of information:

| Role | Located examples | Assumptions that remain attached |
| --- | --- | --- |
| Detector data | Selected HST raw dark frames | Detector calibration, reference validity, noise and observing configuration; these are not a complete raw SN survey |
| Calibrated images/flux | HST FLTs and flat/PAM assets; DES flux/error/time/filter/flag tables; RAISIN photometry | Instrument response, zero points, extraction, resampling, masking and selection |
| Fitted summaries | Pantheon corrected magnitudes/distance columns and covariance; RAISIN systematic distance variants; host ages; compressed BAO means/covariance | Light-curve training, bias corrections, peculiar velocity/extinction/population assumptions, covariance construction and compression validity |
| Empirical/calibration assets | SALT3 templates and colour dispersion, passbands/SEDs, calibration operators | Training sample, conventions, uncertainty and range of validity |
| Synthetic controls | SNANA/BBC simulation and FITRES products, mocks, operator experiments | Chosen source population, selection, nuisance model, random streams and simulator closure |
| Reference results | Historical chains, posterior summaries, model predictions and validation records | Original theory, priors, assets, preprocessing and numerical settings; they are not independent observations |

The [Pantheon+ release paper](https://arxiv.org/abs/2112.03863) distinguishes light-curve release from subsequent distance/bias and cosmological inference, while its [calibration analysis](https://arxiv.org/abs/2112.03864) describes retraining and calibration covariance. DESI's [clustering products](https://data.desi.lbl.gov/doc/releases/dr1/vac/full-shape-bao-clustering/) include windows, mocks, covariance and processed likelihoods. Neither release supplies universally assumption-free reduced data.

A theory change must identify which reduction assumptions remain valid. It may require lower-level observations or re-running affected calibration, selection, coordinate mapping, reconstruction or estimator response. Keep original measurements and all transformations; do not silently reinterpret fitted quantities as raw physical observations.

## Concrete availability and missing inputs

The following is a targeted inventory, not a claim that every listed product was decoded or scientifically validated. ACT version wording comes from its sampled README. Planck product labels come from filesystem paths and provenance/configuration records; this audit did not validate their binary contents, completeness, likelihood version or scientific usability.

| Modality | Products identified in the historical collection | Next data/physics requirement |
| --- | --- | --- |
| Optical/NIR SN | Pantheon table/covariance; DES calibrated fluxes and fitted predictions/masks; RAISIN photometry/distances; Dovekie reference declarations | Complete event/calibration/selection linkage for a chosen forward experiment; actual templates, uncertainties and independent calibration data |
| Images/instruments | Selected HST science FLTs, raw/calibrated darks, flats/PAM and operator arrays | Specific raw-science exposure coverage, PSF/background/extraction model, valid reference files and pixel-level verification |
| Cepheids/anchors | Released ladder design, data and covariance FITS, period/metallicity table, overlap records; correlated 37-host compression | Anchor/Cepheid/SN identities, shared calibration covariance, extinction/crowding/selection contract and absolute forward model |
| BAO/LSS | Thirteen-row compressed BAO inputs; external BAO packages and DESI reference configurations/chains; Lyα/full-shape lineage/material | Actual catalogue/estimator/window/reconstruction identities and theory compatibility; physical ruler or explicitly free ruler, perturbations/bias/RSD |
| CMB primary/lensing | ACT DR6 v1.2 README and likelihood-data paths; Planck-labelled primary/native-likelihood files and ACT/Planck/SPT/MUSE-related reference material | Complete pinned likelihood assets, thermal/perturbation predictions, foreground/instrument/estimator assumptions and shared covariance |
| Specific strong lenses | No complete named-system image/time-delay/kinematics dataset verified by this inventory | Lens/source/line-of-sight boundaries, images/PSF, delays, kinematics and degeneracy-breaking observations |
| TRGB, masers, sirens | No dedicated raw dataset verified by this inventory | Identified source products and instrument/calibration/selection covariance; appropriate stellar, geometric or waveform forward models |
| Chronometers, BBN, 21 cm, clusters | Host-age summaries exist, but no complete dedicated measurement chain for these modalities was verified | Actual measurements and population/nuclear/radiation/gas/selection models; age-labelled rows are not direct chronometers |

“Not verified” is not proof of absence: binary containers and many archives remain unread. Filename searches also have false positives. A future acquisition task must identify the exact scientifically usable product, version, role, units, selection and uncertainty, not merely add a folder.

### A useful conditional ladder reference

A historical study records an SN-free Gaussian ladder compression with 3,138 rows, 45 parameters, 37 correlated host-distance parameters and eight nuisance coordinates. It removes 354 SN rows and the SN absolute-magnitude/H₀ columns. The full host covariance is retained rather than replacing it with independent host errors.

This is a possible reconstruction/compression fixture, not raw photometry, an independent H₀ prior or normalized evidence under improper nuisance priors. It remains conditional on released Cepheid selection, Wesenheit law, metallicities, photometric corrections, anchors and covariance. Its host labels and overlap must be checked before reuse, and the same calibration information cannot appear twice in a joint likelihood. The [primary ladder measurement](https://arxiv.org/abs/2112.04510) provides the relevant measurement and sensitivity context.

## Theory closure and inference gaps

Background expansion is only one part of a theory. Growth/RSD need matter and velocity perturbations; lensing needs metric potentials and propagation; CMB needs thermal/ionization and radiation/metric perturbations, primordial modes and projection. Nonlinear simulations need theory-consistent initial conditions, force/field evolution, constraints, boundaries and declared baryonic/subgrid closure. The [upstream CLASS source](https://github.com/lesgourg/class_public/blob/master/source/perturbations.c) illustrates distinct dependencies; copying a reference engine would not resolve scientific assumptions automatically.

Timescape's averaged geometry and clock/observer mapping cannot be replaced by flat-FLRW distances with a new H(z). Its [original formulation](https://arxiv.org/abs/gr-qc/0702082) and [observable paper](https://arxiv.org/abs/0909.0749) require bare/dressed definitions and supported observer quantities. Missing perturbation or other closure must stay explicit. Modified gravity likewise needs more than a background function; the [Bellini–Sawicki linear scalar-tensor formulation](https://arxiv.org/abs/1404.3713) is one scoped example, not a universal closure.

A particular lens requires actual mass/source/environment and instrument constraints; [mass-sheet degeneracy](https://arxiv.org/abs/2011.06002) limits time-delay inference. A siren distance depends on waveform, detector, orientation and selection; [GW170817](https://arxiv.org/abs/1710.05835) does not remove those assumptions. Alternative gravity may change waveform or propagation, so inherited GR likelihoods need an explicit compatibility decision.

Joint fitting is currently missing, but it is not postponed until every modality is complete. A qualified Cepheid/anchor/SN subset can support its own joint calibration fit. Every subset needs shared nuisance parameters, repeated-object/calibrator identities, survey overlap and relevant cross-covariance. Unknown dependence is not independence. Tests of multiple theories need complexity/prior sensitivity, held-out or posterior-predictive checks and search-multiplicity controls; lowest χ² is not a theory verdict.

## Historical code and licensing

Historical studies cover expansion, selection/censoring, population inference, host ages, dust, spectral/light-curve fitting, passbands, detector geometry, calibration, infrared observations and joint cosmology. These can supply source-linked mathematical controls, discrepancy witnesses and independent comparison algorithms. They are not automatic production implementations or qualified posterior results. Some historical joint analyses explicitly record higher-accuracy likelihood sensitivity beyond their allocation; their posterior summaries must not be adopted as unquestioned expected truth.

Provenance manifests record source origins, extracted symbols, pinned identities and adaptations. Review those records before selecting a small reusable fixture. Preserve original equations and separately identify repaired variants. Keep legacy Python/reference engines outside production and isolate external comparisons with their exact code, assets, flags and licence conditions.

The sampled licence manifest declares MIT for MC-Age and BSD-2-Clause for ACT lenslike, but leaves many Pantheon/Dovekie/DES/SNANA/BAO/SFD entries undeclared and Cobaya labelled Other/NOASSERTION. These are historical declarations, not a current legal determination. Inspect exact pinned per-file notices and dataset terms before copying or redistribution. Public availability does not imply public-domain status, and a code licence does not necessarily cover observational assets.

## Inventory coverage and limits

A read-only filesystem traversal of the historical research root recorded 245,904 file entries and about 108.3 GB of logical file sizes. About 102.8 GB lies in its work area, including simulations, dependencies, copied reference engines, chains and real scientific assets. The curated data directory is about 765 MB. Counts include duplicates and may include symlinks; volume is not independent information or scientific coverage.

The provenance manifests enumerate 391 curated inputs, 4,211 study-source records and 2,926 restoration inputs. Their declarations were structurally inspected, not all rehashed or revalidated. Initial targeted content inspection covered 17 named text/source/manifest files and 189 first-line or FITS-header samples. A subsequent full read covered 11 key authored shared/ladder modules and 29 foundation, scientific, data, execution, agent-interface, reproduction and validation owner plans. These counts describe content coverage, not validation of every equation or the 4,211 source records. No whole-archive binary decoding, posterior rerun or 102-GB scientific reanalysis was performed. Large images/spectra, chains, most source bodies and many containers remain unread. Dependency/cache code was not counted as authored research. The detailed inventory remains a local review artifact rather than a large public path catalogue.

The analytic one-passband experiment, native synthetic shared-calibration ladder and conditional sound-horizon integral with supplied drag epoch now supply bounded executable foundations in the roadmap. Predicting that epoch adds an explicit thermal/ionization/drag contract; full perturbation and CMB work are separate later dependencies.

The next useful work is to choose an actual vertical from the [roadmap](roadmap.md), pin its inputs and forward assumptions, and qualify its own observable and downstream error. Neither broad inventory completion nor a successful numerical call closes these physical and data gaps.

Deterministic photometry predicts three radiometric outputs for a finite constant rest spectrum and supplied distance through the CLI. The standalone C++ sampled operator extends this to declared piecewise-linear wavelength spectra and optical passbands. Pooled sampled CLI/ABI ingestion is implemented. The native [temporal photometry](temporal-photometry.md) consumer adds declared spectral-time grids and observer exposures. Native [detector and censoring](detector-selection.md) now adds a bounded declared Poisson/Gaussian measurement law, addressed simulation and a threshold likelihood with synthetic recovery. Measured instrument laws, source populations and astrophysical recovery remain open. The native [finite calibration law](photometry-calibration.md) propagates a declared shared passband ensemble. Its [joint detector composition](optical-detector.md) now marginalizes every state's conditional count/censoring law and retains the shared all-band selection denominator; moments alone cannot replace that law. Acquiring and qualifying an actual calibration distribution remains separate.

The native [synthetic ladder](calibration-ladder.md) implements an empirical supplied-shape joint relative fit with anchors, Cepheids, calibrator/Hubble-flow SNe and one shared calibration coordinate. That relative fit has no posterior measure. The separate [proper-prior consumer](calibration-predictive.md) supplies normalized Gaussian conditioning and joint future prediction for declared synthetic inputs. Actual source/selection/calibration reconstruction, the released box target and an observed H0 result remain separate gates.

The native early/late operator now shares the sound-horizon model's matter/radiation and flat geometry identity with distance and conditional BAO-ratio predictions. The native [conditional BAO density](bao-conditional.md) composes those predictions with a retained ordered ratio covariance. Its supplied drag epoch, fixed matter content and named numerical controls do not establish a predicted thermal history, physical validity of a released compression or perturbation/CMB closure; see [early and late expansion](early-late.md).

The native [correlated proper calibration](correlated-calibration.md) calculation supplies normalized observed-residual densities for fixed response and an explicitly independent proper latent prior. It closes this statistical prerequisite for named synthetic controls; response construction, calibrated-data dependencies, posterior inference and propagation through a physical observation model remain separate work.

Synthetic proper-prior ladder conditioning and joint future prediction are now
implemented in the [native consumer](calibration-predictive.md). The exact
linear controls distinguish fixed-truth sampling coverage from a conditional
prior-predictive law. [Seeded synthetic recovery](gaussian-simulation.md) now
tests both ensembles with full error covariance and separately targeted coverage,
preserving every attempt and the finite emitted-law approximation. A released
box-supported target, observational held-out event/calibration dependence and
source-supported relation/selection variants remain separate work.
