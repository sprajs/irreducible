# Development workflow

Start with [AGENTS.md](../AGENTS.md), [architecture](architecture.md) and [scientific contracts](scientific-contracts.md). Use [getting started](getting-started.md) to prepare the toolchain.

## Build for a scientific use

Start with a calculation an agent needs. Add a physical model, reader or numerical operation with explicit inputs and machine-readable results. Keep repeated work in compiled batches. External fitters and samplers can call the tool; implement internal fitting only for a named consumer that needs it. A broad roadmap does not require implementing every analysis method.

1. Identify the consumer and qualified prerequisites. State the equation, semantics, domain and error budget.
2. Design an independent comparison and adversarial cases. Check reference ancestry and terms.
3. Implement the smallest coherent change, with explicit model/operation identity where appropriate.
4. Add native scientific tests and boundary/CLI tests. Run the affected tests and investigate discrepancies.
5. Update the public documentation and discovery metadata. Report remaining limits separately from implemented behavior.

`schema/abi.json` owns shared ABI definitions; regenerate through the build driver. Keep changes to generated bindings alongside their schema changes. Use coarse calls and explicit buffer ownership. C++ exceptions and Rust panics must not unwind across the ABI.

## Build discipline

`python3 tools/build.py` is the integrated entry point. It configures CMake, generates bindings/manifests and builds Rust against the native library. The default profile uses Debug. `--profile release` uses a separate native build directory and optimized Rust/C++ binaries, with conservative floating-point flags and four jobs in both profiles. Changing flags, compiler, architecture or backend requires recording identity and revalidating the affected contract.

Do not mutate source while another build is reading it. Coordinate shared CPU/memory budgets. The exclusive build-identity test belongs in an isolated checkout and must run directly as described in [testing](testing.md).

Use `irred` for generic C++ namespaces and `irred_core` for the native library. Public headers live under `cpp/include/irred/`. Expansion models use `irred::cosmology`; supernova and BAO consumers use `irred::supernova` and `irred::bao`. Cosmology is not the root for generic numerical or observation infrastructure. Keep one current public interface; coherent changes may replace obsolete C ABI names and callers together. Preserve immutable scientific/method/source identities rather than iteration plumbing.

## Publishing

Use coherent commits with explicit paths, respect `.gitignore`, and review `git diff --cached`. Keep the docs current with changes to public behavior.

Choose a coherent PR scope. A shared-interface migration may need code, callers, tests and documentation to land together, even when the diff is large. Keep useful validated intermediate commits, document dependencies and separate unrelated work. A dependent PR must say which base change it needs; validate the final integrated candidate rather than treating earlier component checks as its gate.

Agents create a `codex/` branch; human contributors can choose their branch name. Test locally, commit coherent changes and push the branch. Open a PR when the candidate is ready for review. Include the checks actually run, scientific or engineering limits, and any unresolved discrepancy. Inspect CI for the PR commit. Coordinate one integration owner for shared agent work; this does not grant permission to push `main` or merge. Those actions, and force pushes, require an explicit user request.

CI starts when a PR targeting `main` is opened, reopened or receives new commits, and checks the merged candidate. Ordinary branch pushes before a PR do not start CI. Each push to `main` runs the post-merge checks; manual runs are also available. Superseded runs for the same PR are cancelled, while `main` runs are never cancelled by this policy. Update a dependent PR's integration base and rerun when needed. The required `public-tree` check must pass against an up-to-date base, and review conversations must be resolved. This is a personal project: no additional reviewer count is required.

Merges use merge commits to preserve the structured commits and their ancestry; squash, rebase and automatic merging are disabled. Merge authorization still comes from the user, not from a green check.

If the integrated candidate fails, preserve the failure and fix or explicitly block the affected change; passing component tests do not override a failed integrated gate. If `main` turns red after a merge, stop new merges and withhold the affected claims, preserve the failing commit and logs, and prioritize a focused repair or revert PR. A revert or repair still needs explicit authorization to merge; do not bypass the failure with a direct push or weakened checks.

### Clean up after a verified merge

After an authorized merge, verify the PR state, its exact head commit and merge commit. Wait for the post-merge CI result on that exact `main` commit. GitHub automatically deletes the merged remote branch; local cleanup remains the agent's responsibility and cannot be done safely by a GitHub runner.

For local cleanup, first check `git status --short` and `git worktree list`. Preserve a dirty checkout, ongoing work, extra unmerged commits or a branch used by another worktree. Record the verified PR head and merge SHA from `gh pr view <number> --json state,headRefOid,mergeCommit` (or GitHub). The state must be `MERGED`, the merge SHA must have passing post-merge CI, and the local branch tip must equal that PR head. If the branch is already absent, report cleanup complete without attempting deletion. Then use this procedure, substituting the actual branch name and verified commit IDs:

```sh
(
set -eu
task_branch=codex/your-change
pr_head=REPLACE_WITH_VERIFIED_HEAD_SHA
merge_sha=REPLACE_WITH_VERIFIED_MERGE_SHA
test -z "$(git status --porcelain)"
git fetch origin --prune
test "$(git rev-parse "$task_branch")" = "$pr_head"
git merge-base --is-ancestor "$merge_sha" origin/main
git merge-base --is-ancestor main origin/main
git switch main
git merge --ff-only origin/main
git merge-base --is-ancestor "$task_branch" origin/main
git branch -d "$task_branch"
)
```

The subshell stops on the first failed command without exiting your interactive shell. The ancestry check must pass before deletion: `git branch -d` alone may compare against the feature branch upstream, which is not proof of integration into `main`. Verify the fetched remote SHA against the merged commit and any later known merges. The local `main` ancestry check also rejects an ahead or diverged local `main`; a fast-forward-only merge alone would report an ahead branch as already up to date. If the local tip contains extra work or any check fails, preserve the branch and report why cleanup is deferred. Never use `git branch -D`, reset away work or force a fast-forward. There is no cleanup daemon, global hook or background helper.

Review documentation and wiki mirror changes in the same PR when they explain the capability. [Maintenance](maintenance.md) describes the separate wiki publication step.

`python3 tools/build.py --profile release --jobs 3` selects three build jobs while preserving the Release compiler settings. `--jobs` accepts 1 through 4 and defaults to 4; the selected value is recorded in the build manifest. Installation checks use the explicit matching profile: `python3 tools/check_install.py --profile release` (default debug).

Before fitting, resolve the actual covariance/precision, selection, source roles and nuisance dependencies. A bounded controller preserves requested ordering and all failed items rather than fabricating success by dropping rows. A source change invalidates derived scientific identities. Regenerated bytes need their own identity even when a decoded-content comparison establishes equivalence. Label reconstruction exact, approximate, conditional or blocked separately from algorithm agreement and physical validation.
