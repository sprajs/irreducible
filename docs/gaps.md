# Capability and data gaps

This assessment separates the current numerical engine from the proposed [roadmap](roadmap.md). It does not qualify a cosmological interpretation, claim that historical analyses were rerun, or treat a file's presence as a valid likelihood. Use [capabilities](capabilities.md) and actual executable discovery for the implemented interface.

## What exists

The current engine exposes ten operations: exact integer addition, typed physical conversion, scalar numerical methods, immutable observation preparation, requested flat-FLRW background outputs, Gaussian calculations, conditional SN magnitude profiles and conditional free-ruler BAO densities, deterministic rectangular-passband photometry, and a conditional supplied-drag sound horizon. The late-time backgrounds are LCDM, constant q, CPL and fixed five-bin q. They are bounded CPU calculations with explicit source, arithmetic, numerical and execution contracts. Named wide-arithmetic gates use Linux/GCC long double with at least 64 mantissa bits; hosts where long double equals double are not covered, and Ubuntu CI is a separate engineering check. The early ruler separately assumes flat pressureless matter, massless radiation and Lambda; photon/baryon fractions and drag redshift are supplied, not inferred, and no physical-ruler BAO qualification follows. Gaussian nuisance elimination through the CLI supports a scalar offset. The standalone C++ [multi-column design profile](gaussian-design.md) supplies a bounded linear prerequisite; its relative score is separate from the native synthetic ladder and proper correlated calibration density described below.

No current runtime request receives a named scientific qualification. Passing the default numerical_contract means that required numerical checks passed; it is not inference or interpretation qualification. Named comparison evidence applies to its recorded inputs/build/domain, not automatically to a new request.

Shared source ownership and retained factors support repeated evaluations. Named original-input SN and BAO comparisons and analytic/adversarial controls provide evidence for stated cases. They do not establish all-domain accuracy, physical completeness, a joint fit or a new H₀ measurement. The current SN free-offset profile does not retain absolute calibration information; free H₀r_d BAO does not identify H₀ or r_d independently. A supplied dimensional H₀ in a background projection is an input, not a measured result.

The main missing capability is a coherent forward physical state and its observational operators: source spectra/light curves, instrument/pixel response and selection; an observationally qualified absolute distance-ladder calibration model; early thermal/perturbation evolution; specific lenses; nonlinear initial conditions/dynamics/light cones; and qualified GPU/distributed execution. Current historical performance evidence does not benchmark the consolidated interface or these future workloads. A new benchmark must measure representative phases, setup/memory and matched scientific quality.

Engineering remains serial CPU execution with bounded dense matrices and narrow source-specific ASCII adapters. General survey image/spectrum ingestion is absent; the optional FITS codec has test-only synthetic BINTABLE coverage and is unavailable in the product. There is no qualified distributed restart/sharding, RNG or device execution, nor a consolidated-interface throughput benchmark. Adding a model currently requires coordinated native variants/validation and Rust descriptors, schema, ABI, discovery and records. Centralize structural descriptors where useful; keep equations in compiled scientific owners.

## Standard-model diagnostic gaps

The [bounded LambdaCDM baseline](lcdm-baseline.md) distinguishes full Planck base
LCDM from the current supported approximation. The native shared early/late
state treats radiation as massless and matter as pressureless at every epoch.
The separate native [thermal-neutrino provider](thermal-neutrino.md) now
supplies explicit relic density/pressure and flat E/H with a mass transition,
independent high-precision and matched CLASS controls. It is not yet connected
to the early/late distance or supplied-drag ruler consumer. Temperature/species
and physical-density mapping still need an explicit source contract; no operator
predicts thermal/ionization history. Drag redshift is supplied. The CLI late LCDM background omits radiation and has a
different physical identity; native early/late and conditional BAO currently
require a C++ SDK consumer rather than a CLI/C ABI request.

The historical released-ladder admission failure at source 99ce262 remains
preserved. The native retained whitened pivoted QR at source f844080 now admits
the exact same 3492-row/47-column compact products with both covariance arithmetic
profiles and unchanged budgets. All coefficients and the objective agree with
independent LAPACK QR/SVD references; q=3552.759330295523. This closes the named
numerical admission blocker. Physical parameter labels, the released calibration
assumptions, the released estimator contrast and observational H0 reconstruction
still need separate qualification. Native QR contrast variance and a synthetic
ladder H0 sampling law now supply conditional uncertainty under fixed Gaussian
observation noise and an explicitly assumed generating mean. Named synthetic
checks do not qualify released-data parameter mapping or coverage. The relative
profile is not a normalized posterior.

No current operator predicts matter transfer/growth, lensing potentials or
images, recombination/drag, CMB spectra, stellar/source populations, detector
noise, detection or survey recovery. A known cosmology makes these missing
sectors concrete; an H(z) or distance match does not supply them. Use the
[roadmap](roadmap.md) to choose the next qualified prerequisite after the
experiment identifies a blocker.

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

Deterministic photometry predicts three radiometric outputs for a finite constant rest spectrum and supplied distance through the CLI. The standalone C++ sampled operator extends this to declared piecewise-linear wavelength spectra and optical passbands. Sampled CLI ingestion, time dependence, source populations, noise, selection and recovery remain proposed. The native [finite calibration law](photometry-calibration.md) now propagates a declared shared passband ensemble; acquiring and qualifying an actual calibration distribution remains separate.

The native [synthetic ladder](calibration-ladder.md) now implements an empirical supplied-shape joint relative fit with anchors, Cepheids, calibrator/Hubble-flow SNe and one shared calibration coordinate. Named independent recovery and held-out controls do not establish fitting of the released ladder, posterior uncertainty or an observed H0 result. Actual source/selection/calibration reconstruction and any proper-prior normalization remain separate gates.

The native early/late operator now shares the sound-horizon model's matter/radiation and flat geometry identity with distance and conditional BAO-ratio predictions. The native [conditional BAO density](bao-conditional.md) composes those predictions with a retained ordered ratio covariance. Its supplied drag epoch, fixed matter content and named numerical controls do not establish a predicted thermal history, physical validity of a released compression or perturbation/CMB closure; see [early and late expansion](early-late.md).

The native [correlated proper calibration](correlated-calibration.md) calculation supplies normalized observed-residual densities for fixed response and an explicitly independent proper latent prior. It closes this statistical prerequisite for named synthetic controls; response construction, calibrated-data dependencies, posterior inference and propagation through a physical observation model remain separate work.
