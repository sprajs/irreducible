# Irreducible

What would the sky look like if a cosmological idea were true?

Irreducible is being built to turn that question into calculations you can run. Imagine reading a paper that proposes a different history of cosmic expansion. You could implement its equations, predict how bright distant supernovae should appear, and compare that prediction with observations. As the engine grows, the same idea should face other tests too: galaxy clustering, gravitational lenses and the cosmic microwave background.

The goal is to share the physics behind those predictions, from simulated sources to the spectra, images and catalogues we measure. We want to explore questions such as how fast the Universe is expanding, how Cepheids calibrate supernova distances, and what bent light tells us about gravity. Agents can help write and test models; researchers can use the calculations with their own fitting and plotting tools.

Today you can run flat-FLRW background calculations, conditional supernova and BAO comparisons, and Gaussian calculations through the **`irred`** CLI or C++ library. These are the starting tools, not yet a complete simulator or joint cosmological analysis. [Current capabilities](docs/capabilities.md) explain their scope; the [roadmap](docs/roadmap.md) shows what comes next.

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
matched CLASS controls. Its distance/ruler coupling remains a next step.
[Conditional estimator variance](docs/gaussian-design.md) and the
[synthetic ladder's H₀ sampling law](docs/calibration-ladder.md) test a separate
calibration prerequisite. These native calculations are available through the
C++ SDK; their guides distinguish tested controls from observational inference.

## Repeat a paper’s calculation

A planned recipe format would keep a calculation's data sources, model choices and steps in one file, so another person or agent could repeat it and inspect its assumptions. The [run-recipe proposal](docs/run-recipes.md) compares JSON, TOML and YAML and defines the intended first step; the runner does not exist yet. New physics will still be added as tested compiled code, reviewed through a pull request—not as equations executed from a configuration file.

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
