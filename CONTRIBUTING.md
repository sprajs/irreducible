# Contributing

PRs are welcome. Good contributions include code, documentation, counterexamples, clearer contracts, faster algorithms and tests that expose a wrong answer.

**All AI tools are allowed.** Human-written, AI-assisted and agent-written changes are welcome on the same terms. You are responsible for understanding the change, checking it, and having the right to contribute it. Do not submit generated output you cannot explain. Mention material limits or missing validation; no special AI disclosure ritual is required.

## Start small

Read [AGENTS.md](AGENTS.md) and [docs/README.md](docs/README.md). For a large model, new dependency or architectural change, open an issue to discuss the contract and intended scope. Small fixes can go straight to a PR.

1. Make a focused change on a branch. Keep unrelated cleanup separate.
2. Explain the problem and resulting behavior. Scientific changes need equations, assumptions, domains and an error-budget rationale.
3. Run the relevant [tests](docs/testing.md), including native and boundary tests for changed scientific behavior. State the exact checks and any that you could not run.
4. Update public docs and executable discovery when their contract changes.
5. Review the diff for temporary files, data rights and copied material before opening the PR.

Tests should challenge the contract rather than duplicate implementation logic. Preserve failed comparisons locally and add a small durable regression. Never fix a failing test merely by widening a tolerance or regenerating its expected result.

## What belongs in a PR

Commit source, schemas, documentation, native tests and small intentional fixtures with concise provenance. Keep `Plan/`, `evidence/`, acquired datasets, run stores, oracle environments and bulk logs local. Public tests must work without those folders or external comparison software. A useful fixture note states the derivation/source, conventions, uncertainty and tolerance rationale.

Third-party code and assets require compatible terms and proper attribution. Research data is not automatically covered by this repository's license. Production physics/statistics stays in the owned C++ core; reference tools are for independent comparison, not hidden runtime dependencies.

## Review and licensing

Be specific, kind and willing to change your mind. Review the work, not the person. See the short [code of conduct](CODE_OF_CONDUCT.md).

By submitting a contribution, you agree that your original contribution may be distributed under the project's [BSD 3-Clause license](LICENSE). Preserve applicable third-party notices. There is no separate contributor agreement.
