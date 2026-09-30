# Working in Irreducible

Irreducible is building a shared physics engine that takes a scientific idea through to predicted observations. Start with the question and the equations it needs. Keep available calculations distinct from the larger simulation and joint-fitting goals.

Rust handles non-interactive configuration, data acquisition, structural parsing and records. C++20 handles physical models, shared equations, numerical kernels and statistical calculations. The library remains useful independently of the CLI. External fitting, orchestration and plotting tools are welcome consumers.

Use the single active [roadmap](docs/roadmap.md), [gaps](docs/gaps.md), [docs/README.md](docs/README.md), [capabilities](docs/capabilities.md) and the relevant guide when their context helps; inspect actual built discovery/code before claiming implementation. Local planning/history is context, not a feature checklist. Preserve unrelated work and original inputs.

## Lasting design rules

Keep each shared equation and convention in one place. Add models in compiled source, not expressions from data or a runtime plugin framework. Physical models, geometry, observed data and source effects are distinct. Build the smallest useful consumer of qualified prerequisites.

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

Work on a `codex/` branch and open a coherent, reviewable pull request. Scope follows the change: a coupled interface migration can be large, while unrelated work belongs in another PR. Use validated intermediate commits and explain dependencies. Commits, branch pushes and PR creation are part of the normal workflow; direct pushes to `main` and merges require an explicit user request. Coordinate one integration owner when agents share files. Review staged paths and the PR diff, report actual checks and limitations, and inspect CI for the latest integrated PR candidate. Do not force-push without specific authorization.

After an authorized merge, verify the PR head/merge SHA and wait for CI on that exact `main` commit, fetch/prune, update `main` fast-forward-only, and prove the branch tip is an ancestor of fetched `origin/main` before using `git branch -d` from a clean checkout. Preserve extra commits and branches used by active worktrees; never force-delete them. Remote merged branches are automatically deleted; that does not clean local clones. Follow the concrete procedure in [development](docs/development.md), and report deferred cleanup rather than discarding work. If `main` fails, stop new merges, preserve evidence and prioritize a repair/revert PR without weakening checks or bypassing explicit merge authorization.

Treat wiki pages as reviewed documentation: edit their repository mirror in `docs/wiki/` with the PR, then sync the approved pages after merge. Publishing the separate wiki before merge requires an explicit user request. Never present a proposed recipe or command as an implemented capability.
