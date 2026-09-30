# Maintaining the docs

The README introduces the project; [the documentation index](README.md) helps readers choose a task. Keep commands, inputs, result semantics and build instructions with the code they describe. The [roadmap](roadmap.md) owns future direction; discovery from the built executable owns what can run today.

Make documentation changes on a branch and review them in a pull request; agents use the `codex/` prefix. Check examples against the current interface. Label future designs as proposals and avoid commands that suggest an unimplemented feature is available. Run `python3 tools/check_docs.py` for tracked files and relative Markdown links; code changes also need their affected native and Rust tests.

## Wiki publication

`docs/wiki/Home.md` and `docs/wiki/_Sidebar.md` are the reviewed mirror of the GitHub Wiki. The wiki has its own Git history at `https://github.com/sprajs/irreducible.wiki.git` and does not use the repository's PR flow. Normally, merge the reviewed mirror first, then copy those exact pages into a wiki checkout, commit and push. Publishing before merge needs an explicit user request.

Before syncing, compare the live wiki with the mirror and preserve unrelated edits. Verify the pushed wiki commit and read back the pages. Links to new documentation must resolve: use merged pages, or an immutable commit link clearly labelled as a proposed design. Never imply that an open PR is already on `main`.

Keep citation information at the bottom of the README; add a release or paper citation when one actually exists.
