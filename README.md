# Irreducible

What would the sky look like if a cosmological idea were true?

Irreducible is being built to turn that question into calculations you can run. Imagine reading a paper that proposes a different history of cosmic expansion. You could implement its equations, predict how bright distant supernovae should appear, and compare that prediction with observations. As the engine grows, the same idea should face other tests too: galaxy clustering, gravitational lenses and the cosmic microwave background.

The goal is to share the physics behind those predictions, from simulated sources to the spectra, images and catalogues we measure. We want to explore questions such as how fast the Universe is expanding, how Cepheids calibrate supernova distances, and what bent light tells us about gravity. Agents can help write and test models; researchers can use the calculations with their own fitting and plotting tools.

Today you can run flat-FLRW background calculations, conditional supernova and BAO comparisons, and Gaussian calculations through the **`irred`** CLI or C++ library. These are the starting tools, not yet a complete simulator or joint cosmological analysis. [Current capabilities](docs/capabilities.md) explain their scope; the [roadmap](docs/roadmap.md) shows what comes next.

The native [Gaussian predictive](docs/gaussian-predictive.md) calculation retains
a normalized joint future distribution under a proper Gaussian prior, fixed
linear responses and declared independent future noise. Shared calibration
uncertainty induces covariance across those synthetic predictions. The
[proper-prior synthetic ladder](docs/calibration-predictive.md) now composes
that law with shared anchor/Cepheid/SN equations and explicit held-out lineage.
Its posterior H0 projection remains distinct from fixed-truth sampling coverage
and observational calibration.

The native [pressureless GR growing mode](docs/gr-growth.md) predicts D(a) and f(a) under an explicit radiation-free LCDM approximation, with independent differential-equation controls.

The native [hydrogen equilibrium](docs/hydrogen-equilibrium.md) calculation uses
supplied temperature/density and fixed atomic assets. The native
[H/He equilibrium](docs/hydrogen-helium-equilibrium.md) extension solves shared
electron neutrality for independently supplied nuclei densities. The
[supplied baryon abundance](docs/baryon-abundance.md) consumer maps explicit
physical density, He4 mass fraction and neutral effective masses into that state.
This supplies no abundance prediction or helium kinetic history. A separate
[pure-hydrogen history](docs/recombination-drag.md) solves a bounded effective
three-level model with either prescribed radiation temperature or coupled
Compton/adiabatic matter-temperature evolution. The coupled state supplies
Thomson optical depth and finite-endpoint visibility with explicit survival mass.
Bounded explicit massive relics can feed its retained thermal expansion state,
with distinct identities and measured momentum work. The default work cap stays
fixed; expensive profiles require an explicit larger allowance. An explicit
[nested FD integration](docs/thermal-neutrino.md) option accelerates tested
thermal consumers with the same numerical allocations; direct integration
remains the default.
Its truncated drag depth and conditional unit-depth root remain distinct from
a physical drag epoch and the existing supplied-drag ruler.

The native [synthetic SIS thin lens](docs/sis-thin-lens.md) composes compatible
massless-radiation distances with an axial point source, an explicit mass sheet
and Gaussian PSF pixels. It predicts images and relative arrival delays;
synthetic recovery and degeneracy controls do not qualify a measured lens system.

The native [spectral-time photometry](docs/temporal-photometry.md) operator integrates finite source grids over declared observer exposures, with explicit clipping and source epoch.
The native [detector and censoring](docs/detector-selection.md) calculation composes expected photons with a declared Poisson/read-noise law, addressed random draws and a threshold likelihood retaining nondetections.

## Test an established baseline

The [bounded LambdaCDM baseline](docs/lcdm-baseline.md) connects reviewed papers in
[Prospector](https://github.com/sprajs/prospector) to pinned experiments in
[Reproducible](https://github.com/sprajs/reproducible) using this compiled engine.
It checks supported predictions and released-data likelihood calculations while
recording missing physics. The full Planck base model and the current native
massless-radiation approximation have distinct identities; numerical agreement
is scoped evidence, not a guarantee across every cosmological domain.

The native library now also has an explicit thermal-neutrino
[density/pressure and E/H provider](docs/thermal-neutrino.md), tested against
matched CLASS controls. The native [thermal distance/ruler consumer](docs/thermal-observables.md)
now maps explicitly supplied physical densities and temperatures into one retained
state for distances and a conditional supplied-drag ruler. The native
[thermal BAO density](docs/bao-thermal.md) reuses that state and the retained
ordered covariance, with explicit density-level numerical admission.
`photometry.predict` accepts finite sampled spectra and optical passbands in a
bounded pooled batch, preserving declared input roles and calibration limits.
[Conditional estimator variance](docs/gaussian-design.md) and the
[synthetic ladder's H₀ sampling law](docs/calibration-ladder.md) test a separate
calibration prerequisite. A [proper Gaussian parameter posterior](docs/gaussian-posterior.md)
adds normalized fixed-linear conditioning under a declared Gaussian prior. These native calculations are available through the
C++ SDK; their guides distinguish tested controls from observational inference.
The sole active [roadmap](docs/roadmap.md) coordinates the reference baseline, joint
fits, native closures and measured verticals; [gaps](docs/gaps.md) records current blockers.

## Repeat a paper’s calculation

A planned recipe format would keep a calculation's data sources, model choices and steps in one file, so another person or agent could repeat it and inspect its assumptions. The [run-recipe proposal](docs/run-recipes.md) compares JSON, TOML and YAML as a design option governed by the active roadmap; the runner does not exist yet. New physics will still be added as tested compiled code, reviewed through a pull request—not as equations executed from a configuration file.

## Try it

Follow the [build setup](docs/getting-started.md), then inspect your build and run a small example:

```sh
python3 tools/build.py
target/debug/irred describe --json
target/debug/irred run tests/fixtures/background-late.json runs/example
```

This example calculates the expansion rate and physical distances at redshift 1 for a flat LCDM model, and saves the results with a run record. The [CLI guide](docs/cli.md) explains scientific requests, and [gaps](docs/gaps.md) describes the physics and data still needed. Agents working on the code should read [AGENTS.md](AGENTS.md); the [documentation index](docs/README.md) has the full guides.

## Contributing

PRs are welcome, including agent-written ones. Bring a coherent change with the checks and explanation it needs; useful contributions do not have to fit a fixed size. See [contributing](CONTRIBUTING.md).

The [CI workflow](.github/workflows/repository.yml) checks native, Rust and CLI paths separately and produces a scoped Ubuntu build artifact. [Testing and version identity](docs/testing.md#ci-checks-and-build-artifacts) explains the deliberate source version, per-run identity, notices and artifact limits.

## License and citation

[BSD 3-Clause](LICENSE).

If you use Irreducible in research, please cite it and include the commit or release you used:

```bibtex
@software{irreducible,
  author = {Prajs, Szymon},
  title = {Irreducible},
  year = {2026},
  url = {https://github.com/sprajs/irreducible}
}
```

S3 is the persistent store for useful bulk data, private evidence and unfinished
cross-run handoffs. Local ignored directories are temporary scratch/cache.
Use the [common storage layout](docs/shared-data.md), publish and pin before
ending work, and evict local copies only after exact-version byte verification.
Read the [cloud startup guide](https://github.com/sprajs/reproducible/blob/main/docs/cloud-startup.md)
for the bounded shared-catalog readiness check; Reproducible owns the startup script.
