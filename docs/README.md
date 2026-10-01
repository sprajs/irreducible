# Documentation

Irreducible is being built to follow a scientific idea all the way to what an instrument would observe. Shared physics should connect simulations of sources and populations to light propagation, spectra, images and catalogue predictions, then support joint fitting when the necessary models and data are ready.

Today there is a smaller compiled engine you can use through a non-interactive CLI or the C++ library. These guides separate those available calculations from the planned engine. Agents can add a model or reader in source, test it and rebuild; external fitting and plotting tools remain useful consumers.

Read [AGENTS.md](../AGENTS.md) first, then choose a guide:

| Task | Guide |
| --- | --- |
| Development direction and missing prerequisites | [Roadmap](roadmap.md) and [gaps](gaps.md) |
| Build and run | [Getting started](getting-started.md) |
| See what is implemented | [Capabilities](capabilities.md), then `irred describe --json` |
| Plan a repeatable paper workflow (proposed) | [Run recipes](run-recipes.md) |
| Predict a synthetic rectangular-passband observation | [Photometry](photometry.md) |
| Make a request | [CLI contract](cli.md) |
| Interpret an output | [Run records](run-records.md) |
| Add a model, reader or calculation | [Development](development.md) and [architecture](architecture.md) |
| Evaluate a supplied-drag early-time ruler | [Conditional sound horizon](sound-horizon.md) |
| Use one early/late state for distances and conditional ruler ratios in C++ | [Early and late expansion](early-late.md) |
| Check scientific assumptions | [Scientific contracts](scientific-contracts.md) |
| Profile an ordered multi-column Gaussian model in C++ | [Gaussian design](gaussian-design.md) |
| Marginalize a correlated proper calibration prior in C++ | [Correlated calibration](correlated-calibration.md) |

| Predict and recover a synthetic anchor/Cepheid/SN ladder in C++ | [Calibration ladder](calibration-ladder.md) |
| Run tests or add a fixture | [Testing](testing.md) and [fixture provenance](../cpp/tests/fixtures/README.md) |
| Maintain these docs | [Maintenance](maintenance.md) |

These docs travel with the source. Read the version for your commit and check the actual build's capabilities and qualification state before using a result.
