# Irreducible

Irreducible is a simulation-first physics engine for agents investigating cosmology and astronomical observations. The goal is one shared implementation of theories and physical models that predicts across supernovae, Cepheids and H₀, lensing, large-scale structure, BAO, CMB and other modalities.

Theory and initial conditions determine a physical state. Source populations, propagation, instruments and selection turn that state into predicted observations. Images, spectra, light curves and catalogues should reuse shared physics where their assumptions agree. Integrated joint fitting is a later consumer of those forward models; external fitters and samplers are also valid.

The interface is non-interactive: **`irred`** accepts explicit requests and returns machine-readable results. Agents add alternative theories or readers in compiled source, test, rebuild and use the same engine. Rust owns acquisition and control; C++ owns the equations and numerical work. Portable CPU execution is implemented; CUDA and distributed execution require concrete workloads, ownership contracts and measured benefit. Visualization is a low priority and can remain external.

It's early. The current implementation has numerical foundations, quantity conversions, retained observation preparation, flat-FLRW expansion, conditional SN profiles, free-ruler BAO densities and run records. A complete physical forward model, observational simulation and integrated joint fitting remain planned. See [what exists today](docs/capabilities.md).

The [roadmap](docs/roadmap.md) is the active development plan; [gaps](docs/gaps.md) separates present capabilities and data from proposals.

## Using it

**Agents: read [AGENTS.md](AGENTS.md), then [docs/README.md](docs/README.md).** The docs describe the available commands, inputs, tests and limitations.

After the [build setup](docs/getting-started.md):

```sh
python3 tools/build.py
target/debug/irred describe --json
target/debug/irred run tests/fixtures/exact-add.json runs/example
```

That last command is a small infrastructure example. For scientific work, check the operation and qualification status reported by your build.

## Contributing

PRs are welcome, including agent-written ones. Thanks for helping. [Contributing](CONTRIBUTING.md).

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
