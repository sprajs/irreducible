# Working in Irreducible

This repository is both a software-development project and, once implemented, an analysis tool. Choose the appropriate track below. The repository contains a reviewed planning package and a tested minimal infrastructure executable. Discover the actual built capabilities; scientific calculations remain unqualified unless their own evidence says otherwise. Do not turn proposed commands into claims that they have run.

## Start every task

1. Read [Plan/STATUS.md](Plan/STATUS.md), [the main plan](Plan/README.md) and [BUILD-ORDER.md](Plan/BUILD-ORDER.md). The build order owns sequencing; topic plans own their contracts; reference snapshots and review history are evidence rather than active requirements.
2. Identify the requested capability, equation/model, data product, numerical domain and relevant owner plans. Read existing code and local instructions before changing it. Preserve unrelated work and original inputs.
3. Distinguish development from analysis. Record what exists, what has been tested, and what remains planned. All scientific capability statuses start unqualified until actual evidence exists.

## Development track

Implement the smallest coherent slice whose prerequisite capabilities have qualified. Start with S00 and [FIRST-TASKS.md](Plan/08-delivery/FIRST-TASKS.md); do not jump to headline cosmological fits before the mathematical and statistical contracts exist.

Rust owns CLI/configuration, structural parsing, acquisition and run records. C++20 owns physical definitions/conversions, models, numerical kernels, likelihoods and inference loops. CUDA kernels are optional C++ backends. Follow [E05](Plan/04-execution/E05-rust-cpp-boundary.md): coarse batches, one schema source, explicit ownership, no cross-language unwinding and no per-row FFI. The C++ library remains usable independently of the CLI.

Add or modify physical hypotheses directly in source, rebuild and rerun. Use small compiled modules and explicit model IDs; do not add a runtime physics plugin framework or duplicate shared equations in separate analyses. Model changes are legitimate scientific work, but they change scientific identity and must not inherit stale cache entries or qualification. Keep immutable historical outputs and their model identity.

For a change, identify equations/assumptions, implement the comparison case, justify the consumer error budget, then write and optimize code. Own domain physics/statistical definitions. Reuse justified infrastructure under [DEPENDENCY-POLICY.md](Plan/08-delivery/DEPENDENCY-POLICY.md); do not smuggle an astronomy or inference package into the production calculation. Upstream algorithms are references; copied code/assets require actual license and attribution review.

Optimize repeated hotspots after measurement: algorithm/layout/batching, compiler vectorization, multicore, runtime-dispatched SIMD, selected intrinsics or small assembly, and GPU/cluster execution where useful. Maintain a portable baseline and explicit precision/ISA contracts. Compare total wall time and memory at matched numerical quality. A wider vector or more threads is not automatically faster. Keep Rust I/O and C++ compute within one resource budget.

Tests must challenge the stated mathematical/scientific contract, not mirror implementation logic. Use analytic limits, independent algorithms, high precision, refinement and adversarial invalid cases. Shared ancestry is not independent evidence. Before updating an expected result, preserve the old result and explain whether a defect, changed convention, changed input or intentional theory change caused the difference. Never hide a discrepancy by adding jitter, dropping rows, weakening tolerance or replacing a failed result silently.

Deliver durable C++ scientific tests and Rust ABI/safe-wrapper/CLI tests wired into the real build and actually run. Follow V01 for independent fixtures and pinned provenance. Python/IDL/astronomy reference scripts and environments are developer-only comparison tooling: inventory them, migrate accepted comparisons and discrepancy regressions into native tests, then clean temporary outputs/tooling/environments after concise immutable evidence preserves reconstructibility. Do not delete original historical research, planning-review receipts or immutable accepted records.

Review changes to shared kernels against their affected consumers. Compare at observable, log-likelihood and estimand level. On a failed gate follow [V02](Plan/06-validation/V02-discrepancy-protocol.md), withhold the affected capability/claim and continue useful independent work. Do not reflexively defer to an established package or to this codebase: follow the equations and evidence.

## Analysis track

Use discovery from the actual built executable to learn capabilities/schemas. The broad CLI described in [A01](Plan/05-agent-interface/A01-cli-contract.md) is proposed; only commands listed below are currently implemented. Once it exists, resolve the full request, check input identities/units/frames/calibrations and required qualifications, then run a bounded batch or compiled analysis. Prefer retained prepared data to one process per likelihood evaluation.

Keep observations, fitted summaries, empirical training/model assets, assumptions and synthetic truth distinct. Track repeated sources, cross-survey overlap, covariance ordering, shared calibrators and selection. Do not substitute a posterior summary for an independent measurement. Unknown dependence or unavailable original assets limits the corresponding inference; it is not silently assumed absent.

Record exact source/build/model/data identities, resolved configuration, resource allocation, RNG policy, numerical settings, output hashes and diagnostics. Report execution, numerical, inference and interpretation statuses separately. A process that completed is not necessarily an accepted result. Do not infer scientific certainty from code reuse or reference agreement.

If the requested equation or model is missing, switch to development: add the compiled capability, compare it, update qualification and rebuild before analysis. Do not execute agent-generated expressions from data files. New readers are narrow source adapters and must preserve the same scientific object semantics.

Reproduction must be labelled exact, approximate, blocked or conditional with reasons. Use [R01](Plan/07-reproduction/R01-migration-map.md) and [R02](Plan/07-reproduction/R02-acceptance-cases.md) for old-research cases. Do not claim the original repo's historical validation was rerun in Irreducible. For posterior claims check the specific estimand's error, support and sensitivity; do not promote a weakly supported reweighting or an unresolved numerical screen.

## Commands that exist now

From the repository root:

```bash
python3 Plan/tools/check_plan.py
```

This validates planning structure, reference snapshot hashes, catalogue routing and stage dependencies. It does not compile Rust/C++, test equations or rerun science. No dependencies need installing to read or structurally check this plan beyond Python 3.9+.

Implemented and verified S00 commands ([receipt](Plan/08-delivery/implementation/verifier-s00-gate.json), [setup](Plan/08-delivery/implementation/BUILD.md)):

```bash
python3 tools/build.py
.build-tools/bin/ctest --test-dir build/native --output-on-failure
cargo test --locked --offline -j4
target/debug/irred describe --json
target/debug/irred run tests/fixtures/exact-add.json runs/example
```

The exclusive native source/flag mutation suite is compiled with `cargo test --locked --offline -j4 --test build_identity --no-run`, then its emitted binary is run directly with `--ignored --exact source_receipt_flags_and_cache_identity` while no other build is running. These commands verify infrastructure and bounded native tests; they do not establish cosmological inference. The implemented `quantity.convert.v1` request in `tests/fixtures/quantity-length.json` also produces ordered tagged conversion results with preserved source metadata. Until matching numerical evidence is registered, its completed finite run returns exit 6 and `accepted=false`; no execution success implies numerical qualification. The build uses four jobs; coordinate shared builds before running it. Add further production commands only after actual implementation and verification.

## Handoff and completion

Update status/manifests when evidence changes, linking the receipt rather than merely marking a checkbox. A handoff records the exact stage, implemented slice, actual commands, comparisons, failures, unresolved decisions and next dependency. Keep one authoritative main plan and build order; put review findings in the review ledger and integrate accepted corrections into their owner documents.

Report changed files, validation actually performed and material limitations. No claim of implementation, acceleration, scientific agreement or exhaustive software coverage follows from completing documentation alone.

## Product identity and validated publication

The product is **Irreducible** and its executable is **irred**. The checkout path remains `cosmology`; C++ scientific namespaces, schema/scientific IDs and immutable historical receipts keep their original identities. Public repository: https://github.com/sprajs/irreducible ; remote `origin`.

The integration owner makes small coherent commits at validated milestones and before handoff, stages explicit related paths and coordinates a stable snapshot with workers before committing. Preserve unrelated and unfinished work. Push validated milestone commits with `git push origin main` once initial publication ownership is released; never force-push. Report actual build/test and qualification state, including failures and pending gates. Workers do not independently commit or push.
