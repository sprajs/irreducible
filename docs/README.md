# Documentation

Irreducible is being built to follow a scientific idea all the way to what an instrument would observe. Shared physics should connect simulations of sources and populations to light propagation, spectra, images and catalogue predictions, then support joint fitting when the necessary models and data are ready.

Today there is a smaller compiled engine you can use through a non-interactive CLI or the C++ library. These guides separate those available calculations from the planned engine. Agents can add a model or reader in source, test it and rebuild; external fitting and plotting tools remain useful consumers.

Read [AGENTS.md](../AGENTS.md) first, then choose a guide:

| Task | Guide |
| --- | --- |
| Run a bounded standard-cosmology comparison across projects | [LambdaCDM baseline](lcdm-baseline.md) |
| Development direction and missing prerequisites | [Roadmap](roadmap.md) and [gaps](gaps.md) |
| Build and run | [Getting started](getting-started.md) |
| See what is implemented | [Capabilities](capabilities.md), then `irred describe --json` |
| Plan a repeatable paper workflow (proposed) | [Run recipes](run-recipes.md) |
| Predict a synthetic rectangular or sampled-passband observation | [Photometry](photometry.md) |
| Integrate a declared spectral-time source over observer exposures in C++ | [Temporal photometry](temporal-photometry.md) |
| Simulate bounded detector counts and score threshold censoring in C++ | [Detector and censoring](detector-selection.md) |
| Propagate a finite shared passband calibration law in C++ | [Photometry calibration](photometry-calibration.md) |
| Marginalize shared optical states in a joint detector/censoring law | [Optical-state detector](optical-detector.md) |
| Make a request | [CLI contract](cli.md) |
| Interpret an output | [Run records](run-records.md) |
| Add a model, reader or calculation | [Development](development.md) and [architecture](architecture.md) |
| Evaluate a supplied-drag early-time ruler | [Conditional sound horizon](sound-horizon.md) |
| Use one early/late state for distances and conditional ruler ratios in C++ | [Early and late expansion](early-late.md) |
| Evaluate explicit thermal relic density, pressure and flat E/H in C++ | [Thermal neutrino background](thermal-neutrino.md) |
| Map physical densities/temperatures and evaluate thermal distances/ruler ratios | [Thermal observables](thermal-observables.md) |
| Evaluate a supplied-drag BAO density with a retained ordered covariance in C++ | [Massless conditional BAO](bao-conditional.md) or [thermal conditional BAO](bao-thermal.md) |
| Calculate a pressureless GR growing mode in C++ | [GR growth](gr-growth.md) |
| Predict bounded perfect-fluid CDM transfer and a declared primordial band variance | [Perfect-fluid transfer](linear-transfer.md) |
| Inspect the open primary CMB projection prerequisites | [Primary CMB contract](primary-cmb-projection.md) |
| Predict supplied-amplitude growth and a synthetic full-covariance RSD density | [Conditional growth/RSD](growth-rsd.md) |
| Evaluate supplied-temperature/density hydrogen ionization in C++ | [Hydrogen equilibrium](hydrogen-equilibrium.md) |
| Evaluate supplied-temperature/density H/He equilibrium with shared electrons | [Hydrogen–helium equilibrium](hydrogen-helium-equilibrium.md) |
| Map supplied neutral-mass baryon/He abundance and compose shared-electron LTE | [Baryon abundance](baryon-abundance.md) |
| Propagate a supplied joint abundance/mass/temperature law | [Finite abundance law](baryon-abundance-law.md) |
| Prepare pure-H ionization/temperature and finite-endpoint optical histories | [Recombination and drag](recombination-drag.md) |
| Predict a synthetic SIS point-source image, delay and Gaussian PSF pixels in C++ | [SIS thin lens](sis-thin-lens.md) |
| Check scientific assumptions | [Scientific contracts](scientific-contracts.md) |
| Profile an ordered multi-column Gaussian model in C++ | [Gaussian design](gaussian-design.md) |
| Normalize a broad finite uniform-box Gaussian target and bound one coordinate quantile | [Gaussian box](gaussian-box.md) |
| Calculate conditional estimator variance and a synthetic H0 sampling law | [Gaussian design](gaussian-design.md) and [calibration ladder](calibration-ladder.md) |
| Marginalize a correlated proper calibration prior in C++ | [Correlated calibration](correlated-calibration.md) |
| Condition a proper Gaussian parameter prior in a fixed linear model | [Gaussian parameter posterior](gaussian-posterior.md) |
| Predict a normalized joint future vector with proper prior and independent future noise | [Gaussian predictive](gaussian-predictive.md) |
| Condition many original training vectors with one retained fixed future covariance | [Gaussian predictive](gaussian-predictive.md) and [ladder predictive](calibration-predictive.md) |
| Condition a proper synthetic ladder prior and predict joint held-out rows | [Ladder predictive](calibration-predictive.md) |
| Generate addressed full-covariance Gaussian vectors and run synthetic recovery campaigns | [Gaussian simulation](gaussian-simulation.md) |
| Predict and recover a synthetic anchor/Cepheid/SN ladder in C++ | [Calibration ladder](calibration-ladder.md) |
| Run tests or add a fixture | [Testing](testing.md) and [fixture provenance](../cpp/tests/fixtures/README.md) |
| Maintain these docs | [Maintenance](maintenance.md) |

These docs travel with the source. Read the version for your commit and check the actual build's capabilities and qualification state before using a result.
