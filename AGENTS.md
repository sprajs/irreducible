# Working in Irreducible

**Read [docs/README.md](docs/README.md) at the start of every task.** It is the canonical agent-first entry point. Public documentation must be sufficient to build and work on a fresh clone. Discover actual executable capabilities before making scientific claims.

## Start every task

1. Read the documentation index, [capabilities](docs/capabilities.md), and the relevant task guide. Identify development versus analysis, the equation/model, input product, numerical domain and affected consumers. Read existing code and local instructions before editing.
2. If a local `Plan/` directory exists, read its `STATUS.md`, `README.md` and `BUILD-ORDER.md`. Its build order owns local implementation sequencing; topic plans own detailed contracts. Historical reviews and reference snapshots are evidence, not new requirements. Public contributors do not need this private planning package.
3. Preserve unrelated work and original inputs. Coordinate with the integration owner before builds or shared edits. Use explicit paths when staging. Scientific statuses begin unqualified; implementation and passing tests alone do not change that.

## Development

Implement the smallest coherent slice whose prerequisites have qualified. Follow [development](docs/development.md) and [scientific contracts](docs/scientific-contracts.md). Do not jump to headline fits before the mathematical and statistical contracts exist.

Rust owns CLI/configuration, structural parsing, acquisition and run records. C++20 owns physical definitions/conversions, models, numerical kernels, likelihoods and inference loops. CUDA is a future optional C++ backend. Use coarse batches, one ABI schema source, explicit ownership, no cross-language unwinding and no per-row FFI. Keep the C++ library independently usable.

Add hypotheses directly in source as small compiled modules with explicit model IDs. Rebuild and rerun. Do not introduce a runtime physics plugin framework or duplicate shared equations across analyses. A changed model changes scientific identity; it must not inherit stale cached results or qualification. Preserve immutable historical results locally.

Identify equations and assumptions, build the comparison case and justify the consumer error budget before optimizing. Own physical/statistical definitions. Reuse justified general infrastructure; do not import an astronomy or inference package as the production scientific engine. Review actual licenses and attribution before copying code or assets.

Measure hotspots and total wall time/memory at matched numerical quality. Improve algorithms, layout and batching before selecting vectorization, multicore, SIMD or accelerator work. Maintain a portable baseline, explicit precision/ISA contracts and one resource budget for Rust and C++. More threads or wider vectors are not automatically faster.

Tests must challenge the contract using analytic limits, independent algorithms, high precision, refinement and invalid cases. Shared ancestry is not independent evidence. Before changing an expected result, preserve the old result locally and explain whether the cause is a defect, convention, input or intentional theory change. Never hide discrepancies with jitter, dropped rows or weakened tolerances.

Keep durable C++ tests and Rust ABI/wrapper/CLI tests wired into the real build. Keep small necessary fixture provenance with the tests: sources, derivations, input conventions, uncertainty and tolerance rationale. Developer-only oracle environments, bulk outputs and receipts belong locally. After migration, native tests must run without external comparison software or local evidence folders. See [testing](docs/testing.md).

Review shared kernels against affected consumers at observable, log-likelihood and estimand levels. A failed gate withholds the affected claim. Preserve its cause and reconstructibility locally, add a durable discrepancy regression, and continue independent work. Neither an established package nor this repository is automatically right.

## Analysis

Use `irred describe --json` from the actual build to discover operations and schemas. Read [the CLI contract](docs/cli.md); proposed commands are not executable capabilities. Check source identities, units, frames, calibrations and required qualifications before a bounded run. Prefer prepared data and coarse batches over per-evaluation processes.

Keep observations, fitted summaries, empirical training/model assets, assumptions and synthetic truth distinct. Track repeated sources, survey overlap, covariance ordering, shared calibrators and selection. Posterior summaries are not independent measurements. Unknown dependence or unavailable original assets limits the inference; it is not silently assumed absent.

Record exact source/build/model/data identities, resolved configuration, resource allocation, RNG policy, numerical settings, output hashes and diagnostics. Report execution, numerical, inference and interpretation statuses separately. Completion is not acceptance. Do not infer scientific certainty from code reuse or reference agreement.

If an equation/model is missing, switch to development, qualify the new compiled capability and rebuild before analysis. Never execute agent-generated expressions from data files. Source adapters must preserve scientific object semantics.

Label reproduction exact, approximate, conditional or blocked, with reasons. Do not claim historical research was rerun here without actual runs. Check the estimand's error, support and sensitivity before posterior claims; weak reweighting or an unresolved numerical screen is insufficient.

## Repository and publication

- Product: **Irreducible**. Executable: **`irred`**. Generic C++ code uses the `irred` namespace, `cpp/include/irred/` and the `irred_core` library. A future `irred::cosmology` module is reserved for actual cosmological models. Legacy `cosmo_*` C ABI names, `COSMO_*` tags and scientific/schema IDs remain wire identifiers; product naming does not rewrite scientific identity.
- **`docs/` is canonical and agent-first.** Update it with every public command, contract or workflow change. Keep README and executable discovery consistent. The Wiki is a navigation layer linking to `docs/`.
- **Never commit `Plan/`, `evidence/`, run stores, acquired datasets, local environments or scratch output.** They stay local and ignored. Do not force-add them. Retain concise provenance and small intentional fixtures under the public tests/docs instead.
- Do not erase original research, local planning reviews or immutable accepted records as routine cleanup. Ignore/untrack local material without deleting it. Git is for the tool, tests and documentation; users retain scientific run records separately.
- The integration owner makes small coherent commits at validated milestones and before handoff, coordinates stable snapshots, stages explicit related paths, then pushes `origin main` and verifies the remote SHA. This is standing user authorization for this repository. Workers do not independently commit or push. Never force-push without explicit user authorization for that specific rewrite.

## Validation and handoff

Use the commands in [getting started](docs/getting-started.md) and [testing](docs/testing.md). Shared builds use at most four jobs. The source/flag mutation test is exclusive: use an isolated checkout and run its compiled binary directly while no other build operates there.

Update public capability/contract documentation when behavior changes; retain detailed receipts locally. A handoff records implemented scope, commands actually run, comparisons, failures, unresolved decisions and the next dependency. If local plans exist, update their status and link their local receipts.

Report changed files, actual validation and material limits. Documentation alone establishes no implementation, speedup, scientific agreement or exhaustive coverage.
