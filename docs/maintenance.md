# Maintaining the docs

`docs/README.md` is the agent entry point. Keep commands, inputs, output semantics and build instructions current with the code. README introduces the project; detailed contracts belong in the guides.

Wiki navigation is maintained in `docs/wiki/`. Copy those files into a checkout of `https://github.com/sprajs/irreducible.wiki.git`, commit and push when the links change.

Run `python3 tools/check_docs.py` to check tracked files and relative Markdown links. For code changes, run the relevant native and Rust tests too. Keep citation information at the bottom of README; add a release or paper citation when one actually exists.
