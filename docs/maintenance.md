# Maintaining the docs

The README introduces the project; [the documentation index](README.md) helps readers choose a task. Keep commands, inputs, result semantics and build instructions with the code they describe. The [roadmap](roadmap.md) owns future direction; discovery from the built executable owns what can run today.

Make documentation changes on a branch and review them in a pull request; agents use the `codex/` prefix. Check examples against the current interface. Label future designs as proposals and avoid commands that suggest an unimplemented feature is available. Run `python3 tools/check_docs.py` for tracked files and relative Markdown links; code changes also need their affected native and Rust tests.

## Wiki publication

`docs/wiki/Home.md` and `docs/wiki/_Sidebar.md` are the reviewed mirror of the GitHub Wiki. The wiki has its own Git history at `https://github.com/sprajs/irreducible.wiki.git` and does not use the repository's PR flow. Normally, merge the reviewed mirror first, then copy those exact pages into a wiki checkout, commit and push. Publishing before merge needs an explicit user request.

Before syncing, compare the live wiki with the mirror and preserve unrelated edits. Verify the pushed wiki commit and read back the pages. Links to new documentation must resolve: use merged pages, or an immutable commit link clearly labelled as a proposed design. Never imply that an open PR is already on `main`.

Keep citation information at the bottom of the README; add a release or paper citation when one actually exists.

## Repository checks

[The workflow](../.github/workflows/repository.yml) checks PRs targeting `main` when they open, reopen or receive commits, and checks `main` again after a push. Ordinary feature-branch pushes do not start a duplicate run. Manual dispatch is available once this workflow is on the default branch. Only a newer run for the same PR cancels its predecessor; main and manual runs stay independent.

Keep the six required check contexts aligned with the `CI` workflow and its GitHub Actions app: C++ tests (gcc), C++ tests (clang), Rust unit tests, CLI integration tests, Release application and package, and Documentation and generated contracts. When renaming checks, first obtain passing new contexts on the actual PR candidate, then atomically replace the old required set and verify protection; never leave an unprotected gap. Main protection requires an up-to-date checked PR and resolved conversations, including for administrators. No path filters or draft skips omit the required check. Repository settings keep merge commits, automatic deletion of merged remote branches and branch-update support enabled; auto-merge stays disabled. These settings do not replace independent review and actual green required CI; agents perform routine qualified merges under the user’s standing authorization.
