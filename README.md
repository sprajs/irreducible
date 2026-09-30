# Irreducible

**Ask a harder question. Keep the calculation honest.**

Cosmology begins with a deceptively simple question: what does the universe actually tell us? Answering it can mean navigating a stack of scripts, conventions, fitted summaries and inherited assumptions before reaching the calculation you wanted to inspect.

Irreducible is being built to bring that work into one coherent scientific tool. An agent can discover its capabilities, prepare a bounded request and run compiled calculations through a single executable, **`irred`**. A researcher can follow the same path back through the equations, inputs and tests. Changing a physical hypothesis means changing an explicit model, rebuilding it and testing what follows.

The name is a design ambition: reduce the machinery until what remains earns its place. Shared equations, reusable numerical kernels, a small command interface and as little repeated work as possible. Speed matters because a question worth asking once is often worth asking a million times—with different parameters, assumptions or data.

> **Agents: start with [AGENTS.md](AGENTS.md), then [docs/README.md](docs/README.md).**
> `docs/` is the canonical, agent-first documentation: task routes, executable commands, scientific contracts and known limits. Discover your actual build with `irred describe --json` before attempting an analysis.

## What you might do with it

Imagine taking a published result and asking which part came from observations and which part came from the model. Then changing one assumption, keeping the rest fixed, and rerunning the calculation with a record of exactly what changed.

That is the workflow Irreducible is working toward:

- **Investigate an assumption.** Change a compiled physical model and compare its consequences against a stated numerical budget.
- **Reproduce a calculation.** Identify original inputs, conventions and dependencies, and say plainly when a reproduction is exact, approximate, conditional or blocked.
- **Explore efficiently.** Reuse prepared data and send batches into native code instead of launching a process for each likelihood evaluation.
- **Work with agents.** Let an agent handle orchestration while explicit schemas, native tests and qualification rules constrain what its results can mean.

These are the intended uses, not a claim that every step is available today. Irreducible is in early development. The current public executable provides strict JSON requests, immutable local run records, checked integer fixtures, tagged quantity conversions and a small scalar numerical interface. The C++ library contains additional numerical building blocks. End-to-end cosmological models, likelihood workflows and inference are still being developed. See [capabilities and limits](docs/capabilities.md).

## Try the foundation

The current build profile targets a Linux CPU toolchain with Rust, a C++20 compiler, Python and CMake. Follow [getting started](docs/getting-started.md) for the complete setup, including dependency preparation. Once prepared:

```sh
python3 tools/build.py
target/debug/irred describe --json
target/debug/irred run tests/fixtures/exact-add.json runs/first-run
```

The last command runs an exact integer fixture and writes a local result and receipt. It is a quick way to inspect the interface before doing scientific work. Scientific requests currently return **exit 6** when a calculation completes but the required numerical qualification has not been registered. Read the receipt as well as the process status; a finite answer alone is not an accepted scientific result.

## One tool, explicit responsibilities

Rust handles commands, request parsing, data movement and run records. C++20 owns physical definitions, numerical kernels and the developing scientific engine. The boundary carries coarse batches; the C++ library can also be used independently. Physical hypotheses live in small compiled modules with explicit identities.

Performance work follows measurement at matched numerical quality. Portable CPU execution is the baseline; wider SIMD, GPU and distributed execution are future work. There is no claim of a demonstrated speed advantage yet.

The project is open to questioning a model—including an established one. That freedom comes with a practical obligation: make the assumptions visible, preserve failed comparisons, and test the replacement. Agreement with another package is useful evidence; it is not a substitute for understanding the calculation.

## Join in

PRs, questions, counterexamples and better tests are welcome. **All AI tools are welcome**, including agent-written contributions. Contributors remain responsible for what they submit. Read [CONTRIBUTING.md](CONTRIBUTING.md) and our short [code of conduct](CODE_OF_CONDUCT.md): be nice to each other, and argue about the work.

- [Documentation](docs/README.md) — the source of truth, written for agents first and readable by humans.
- [GitHub Wiki](https://github.com/sprajs/irreducible/wiki) — an entry point to the documentation.
- [Issues](https://github.com/sprajs/irreducible/issues) — bugs, ideas and scientific counterexamples.

## License and citation

Irreducible is available under the permissive [BSD 3-Clause license](LICENSE). You can use, modify and redistribute it, including commercially, under those terms. If it contributes to research, please cite the software and the exact commit or release you used; [CITATION.cff](CITATION.cff) provides citation metadata. Scientific citation is requested separately from the license conditions. Data and third-party assets retain their own terms.
