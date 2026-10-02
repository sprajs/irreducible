# Irreducible

What should we see in the sky if an idea about the Universe is right?

Irreducible is being built to help turn a paper's equations into predictions we can compare with observations. A different expansion history might change how bright distant supernovae appear. A theory of gravity might change a lens image or the way galaxies cluster. The aim is to follow those ideas through shared physics to spectra, images, light curves and catalogues, rather than build a separate calculation for every test.

We want to explore the expansion rate H₀, Cepheid and supernova distances, lensing, large-scale structure and the cosmic microwave background. Today the engine offers a smaller set of tools: flat-FLRW backgrounds, conditional supernova and BAO calculations, and Gaussian calculations. You can drive them through a non-interactive CLI or use the C++ library; fitting and plotting can live in your own tools.

The native [H/He equilibrium](https://github.com/sprajs/irreducible/blob/main/docs/hydrogen-helium-equilibrium.md) calculation supplies ground-state fractions with shared electrons for supplied temperature and nuclei densities. The [supplied abundance consumer](https://github.com/sprajs/irreducible/blob/main/docs/baryon-abundance.md) maps explicit physical density, He4 mass fraction and neutral effective masses into that state. A predicted abundance and helium kinetic history remain separate prerequisites.

The [proper-prior synthetic ladder](https://github.com/sprajs/irreducible/blob/main/docs/calibration-predictive.md) shares anchor/Cepheid/SN equations with relative recovery and retains a joint law for explicitly disjoint future rows. It supplies no observational H0 result. Thermal consumers can select [nested FD integration](https://github.com/sprajs/irreducible/blob/main/docs/thermal-neutrino.md) at unchanged numerical allocations; direct integration remains the default.

A future recipe format would put the data sources, model choices and calculation steps in one file, so another person or agent could repeat a paper test and see its assumptions. There is no recipe runner yet. New equations belong in compiled code with tests and review.

## Try, explore or contribute

- [Build and run an example](https://github.com/sprajs/irreducible/blob/main/docs/getting-started.md), then browse [what works today](https://github.com/sprajs/irreducible/blob/main/docs/capabilities.md) and the [CLI guide](https://github.com/sprajs/irreducible/blob/main/docs/cli.md).
- Test the [bounded LambdaCDM baseline](https://github.com/sprajs/irreducible/blob/main/docs/lcdm-baseline.md) through the reviewed Prospector → Reproducible → Irreducible workflow; full Planck closure remains distinct from the supported massless variant.
- Read the [roadmap](https://github.com/sprajs/irreducible/blob/main/docs/roadmap.md) for the next experiments and [gaps](https://github.com/sprajs/irreducible/blob/main/docs/gaps.md) for the missing physics and data.
- Have a model, reader or useful test to add? [Contributions](https://github.com/sprajs/irreducible/blob/main/CONTRIBUTING.md) are welcome, including agent-written pull requests. Start with [agent instructions](https://github.com/sprajs/irreducible/blob/main/AGENTS.md) and the [development guide](https://github.com/sprajs/irreducible/blob/main/docs/development.md).

The [documentation index](https://github.com/sprajs/irreducible/blob/main/docs/README.md) has the full guides. They live alongside the code so that the instructions and calculations can be reviewed together.

The [conditional growth/RSD consumer](https://github.com/sprajs/irreducible/blob/main/docs/growth-rsd.md) predicts supplied-reference sigma8 and f*sigma8 with a full-covariance synthetic density. Radiation/relic transfer and the validity of a measured survey compression remain separate prerequisites.
