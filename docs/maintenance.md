# Repository maintenance

The public repository contains the tool, schemas, tests, small intentional fixtures and documentation. Keep private planning, bulk evidence, run stores, acquired datasets and environments local. `.gitignore` prevents accidental additions; `tools/check_docs.py` also rejects those roots if someone force-adds them.

## Documentation ownership

`docs/README.md` is the agent-first entry point. Update relevant guides with public behavior changes, especially operation IDs, exit codes, qualification status and build requirements. README tells the story and links to the guides. AGENTS provides working rules. Avoid maintaining a second detailed command reference in the Wiki.

Wiki entry pages live under `docs/wiki/`. After a reviewed docs change, copy those Markdown files to a checkout of `https://github.com/sprajs/irreducible.wiki.git`, commit and push. The Wiki's first page must be initialized through GitHub before its Git repository can be cloned. Its pages link to canonical docs in the main repository.

## Checks and releases

The repository hygiene workflow checks the tracked public tree and relative file links on pushes and PRs. It is not a scientific qualification workflow. Run the native and Rust tests for code changes and record the tested scope in the PR or handoff.

Before a release, build/test the exact revision, document toolchain/domain limits and verify the distributable from a clean checkout. Add an actual release version and date to `CITATION.cff` only when that release exists. A package manifest's development version is not evidence of a published release. Consider a research archive/DOI when a citable release is ready; do not invent one.

Preserve third-party notices and record data licenses separately. The BSD license covers original repository contributions; it does not grant rights to every dataset a user might supply.
