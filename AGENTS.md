# Working in Irreducible

Irreducible is a personal, agent-first physics provider. Rust owns non-interactive configuration, acquisition, structural parsing and records; C++20 owns shared equations, physical models, numerical kernels and concrete statistical consumers. External fitting and orchestration are valid consumers. A complete inference ecosystem is not required. Reusable supernova calculations are the first substantial target; visualization belongs in external tools.

Use [docs/README.md](docs/README.md), [capabilities](docs/capabilities.md) and the relevant guide when their context helps; inspect actual built discovery/code before claiming implementation. Local planning/history is context, not a feature checklist. Preserve unrelated work and original inputs.

## Lasting design rules

Keep one implementation of each shared equation/convention. Add models in compiled source, not expressions from data or a runtime plugin framework. Physical models, geometry, observed data and source effects are distinct. Build the smallest useful consumer of qualified prerequisites.

Maintain one coherent current interface. Do not create compatibility branches for hypothetical users. Breaking changes may remove obsolete routes, layouts, admission rules and redundant parameters together with their callers/tests/docs. Preserve scientific identities, immutable historical records and meaningful mathematical tests; do not preserve obsolete live plumbing without an actual consumer need.

Use coarse batches and explicit ownership. Acquire/validate an immutable observation object once, pass it through preparation, retain factors across evaluations and avoid copies/readbacks merely to cross a layer. The C++ library stays independently usable. Request only outputs needed by the consumer; no server/cache framework is required to retain a bounded calculation.

## Evidence and records

State equations, assumptions, domains and justified error budgets before coding or optimizing. Test analytic limits, independent algorithms, high precision/refinement and adversarial invalid inputs. Shared ancestry is not independent evidence. Preserve a failed result before changing an expectation; never hide discrepancies with jitter, dropped rows or weaker budgets.

Keep observations, fitted summaries, assumptions and synthetic controls distinct. Track units, frames, calibration, source/axis order, selection and dependence; unknown overlap is not independence. Numerical acceptance does not establish inference or interpretation. Operation code owns output IDs, method/arithmetic and check status; the common recorder must not infer them from operation names.

Measure performance at matched quality and resource limits, including setup/memory. Keep a portable baseline and explicit arithmetic/ISA contracts. Review affected consumers whenever shared kernels change. Review licensing before copying code/assets.

## Repository traps and publication

[schema/abi.json](schema/abi.json) owns the shared C/Rust ABI; regenerate with `tools/generate_abi.py`, never hand-edit bindings. No cross-language exception unwinding or per-row FFI. Public documentation is in `docs/`; keep executable discovery and docs consistent with actual gates.

Coordinate shared files/builds; total compiler/compute budget is four jobs. Run source/flag mutation tests only in an isolated checkout, directly through the compiled test binary, without a competing build. Fresh native directories avoid stale objects when snapshot copies preserve mtimes.

Stage explicit coherent source/test/docs paths, respecting `.gitignore`. Never force-add local planning/evidence or delete originals/immutable receipts as cleanup. Snapshot unfinished work before a migration; remove only inventoried disposable artifacts released by their owners.

The integration owner has standing authorization to commit validated milestones, push `origin main` and verify the remote SHA. Workers do not publish independently. Force pushes require specific user authorization. Report actual commands/results, limits, failures and the next useful dependency.
