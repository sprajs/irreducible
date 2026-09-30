# Irreducible

I want to set my agents to work on cosmology. I'm writing Irreducible to give them a fast, shared implementation of the physics and physical models they need.

I think of it as something an agent can reach for in the way it reaches for a computer algebra system when doing mathematics. Here the job is to read observations, evaluate physical models and work out what follows from different assumptions. The agent decides what to investigate; Irreducible provides the calculations.

It's built for non-interactive use through **`irred`**, a CLI with explicit requests and machine-readable results. Long commands and detailed request files are fine. The interface should let an agent say exactly what it wants to calculate, with which data, model and configuration.

I want one shared implementation of each equation and convention, with different physical models built on top. When an agent needs a new data reader or wants to try different physics, it should be able to add it, test it, rebuild and use it in the same way. Rust handles the interface and data; C++ handles the physics and numerical work. Repeated calculations should stay in compiled code and run in batches.

This doesn't need to contain every part of an analysis. An agent can use an external fitter or sampler and call Irreducible for model predictions. Some calculations may need fitting internally, and that's fine when there's a reason for it. The aim is to build the parts that make the science possible and make them fast.

It's early. The current implementation has numerical foundations, quantity conversions, observation preparation and run records; the physical models are still being built. See [what exists today](docs/capabilities.md).

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
