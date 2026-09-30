# Working in Irreducible

**Irreducible is a tool for agents doing physics, starting with cosmology.** Build a fast, shared implementation of equations and physical models that agents can call, inspect and extend. The primary interface is a **non-interactive CLI** with explicit inputs and machine-readable outputs. Long commands, detailed request files and bounded batches are welcome; interactive prompts and a GUI are not the intended workflow.

Agents should be able to add data readers, process new observations, introduce physical models, change configurations and build other tools on top. Keep one implementation of each shared equation and convention; alternative physical models can use those same foundations. Let concrete scientific questions determine what gets built.

**A complete fitting or inference framework is not a project requirement.** External fitters, samplers and orchestration are valid consumers. Add internal fitting only when a specific calculation needs it. Do not delay a useful physics capability to build unrelated analysis machinery, or treat a broad roadmap as a feature checklist.

## Start here

Read [docs/README.md](docs/README.md), [capabilities](docs/capabilities.md) and the relevant guide. Inspect existing code and discover the actual executable before claiming a capability exists. Identify the equation/model, data semantics, numerical domain and affected consumers. Preserve unrelated work and coordinate shared edits and builds.

## Implementation

Rust owns CLI/configuration, structural parsing, acquisition and run records. C++20 owns physical definitions, models and numerical calculations, including statistical kernels needed by their consumers. Keep the C++ library independently usable. Use coarse batches, one ABI schema source, explicit buffer ownership and no cross-language unwinding or per-row FFI.

Add hypotheses in source as small compiled modules with explicit model IDs, then rebuild and test. Do not execute expressions from data files or introduce a runtime physics plugin framework. Changed models need new scientific identities and cannot inherit stale results or qualification.

State equations, assumptions and the consumer error budget before implementation or optimization. Reuse justified infrastructure; keep production physical definitions in the shared core. External scientific software can provide independent comparisons or drive an analysis. Review licenses and attribution before copying code or assets.

Measure performance at matched numerical quality. Start with algorithms, layout and batching; use vectorization, multicore or accelerators when measurements justify them. Retain a portable baseline and explicit precision/ISA contracts. Rust and C++ share a resource budget.

## Scientific checks

Tests must challenge the contract: analytic limits, independent algorithms, high precision, refinement and adversarial inputs. Shared ancestry is not independent evidence. Keep durable C++ tests and Rust boundary/CLI tests in the real build, with concise fixture sources, derivations, conventions and tolerance rationale.

Preserve a failed comparison and explain its cause before changing an expected value. Never hide a discrepancy with jitter, dropped rows or weaker tolerances. Add a regression and review affected consumers at observable, likelihood or estimand level as applicable. A failed gate withholds the affected claim; independent work can continue.

Keep observations, fitted summaries, calibration/training assets, assumptions and synthetic truth distinct. Preserve source identities, units, frames, covariance ordering, overlap, shared calibrators and selection. Unknown dependence remains unknown. A posterior summary is not an independent measurement.

Record source/build/model/data identities, resolved requests, resource/RNG policy, numerical settings and diagnostics. Report execution, numerical, inference and interpretation status separately where applicable. A completed calculation is not automatically accepted science. Label reproductions exact, approximate, conditional or blocked, with reasons.

## Interface and documentation

Use `irred describe --json` to discover implemented operations and [docs/cli.md](docs/cli.md) for their contract. Return structured results and explicit failures; never silently substitute inputs or change the requested physical meaning. If a capability is missing, implement and validate the required slice before using it.

Public documentation lives in **`docs/` and is written for agents first**. Update it and executable discovery when behavior changes. Keep README as a short introduction and the Wiki as navigation. Distinguish current behavior from intended design.

Generic C++ code uses `irred`, `cpp/include/irred/` and `irred_core`. Actual cosmological models may use `irred::cosmology`. Existing C ABI symbols and scientific/schema IDs keep their identities unless their contracts change.

## Validation and publication

Follow [testing](docs/testing.md). Shared builds use at most four jobs. Run the source/flag mutation test only in an isolated checkout, directly through its compiled binary, with no competing build there.

Commit source, tests and documentation; respect `.gitignore`. Do not force-add ignored material or delete original research and immutable records as routine cleanup. Keep small coherent commits, stage explicit paths and preserve unfinished work.

The integration owner coordinates workers, commits validated milestones and pushes `origin main`, verifying the remote SHA. This is standing user authorization for this repository. Workers do not independently commit or push. Force-pushes require explicit authorization for the specific rewrite.

Report what changed, checks actually run, limitations and the next useful step. Passing tests or documentation work alone establishes no new scientific result or speedup.
