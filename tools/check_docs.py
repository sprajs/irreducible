#!/usr/bin/env python3
"""Check tracked publication hygiene and relative Markdown file links, offline."""
import json
from pathlib import Path, PurePosixPath
import re
import subprocess
import sys
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[1]
tracked = set(subprocess.check_output(
    ["git", "ls-files", "-z"], cwd=ROOT).decode().rstrip("\0").split("\0"))
errors = []
private_roots = {"Plan", "evidence", "runs", "build", "target", ".build-tools",
                 ".venv", "venv", "data", "datasets", "scratch", "tmp"}
for name in sorted(tracked):
    path = PurePosixPath(name)
    if path.parts[0] in private_roots:
        errors.append(f"local-only path is tracked: {name}")
    if path.name == ".env" or (path.name.startswith(".env.") and path.name != ".env.example"):
        errors.append(f"environment file is tracked: {name}")
    if name == "tests/verification/reference_oracles.py":
        errors.append(f"temporary oracle tooling is tracked: {name}")

required = {"README.md", "AGENTS.md", "LICENSE", "CITATION.cff", "CONTRIBUTING.md",
            "CODE_OF_CONDUCT.md", "SECURITY.md", "docs/README.md"}
for name in sorted(required - tracked):
    errors.append(f"missing tracked repository file: {name}")

for name in sorted(tracked):
    if not name.endswith(".md"):
        continue
    body = (ROOT / name).read_text()
    # Only check normal inline Markdown file links; skip fenced examples.
    prose = re.sub(r"^```[^\n]*\n.*?^```\s*$", "", body, flags=re.M | re.S)
    for dest in re.findall(r"\[[^\]\n]*\]\(([^\s)]+)(?:\s+\"[^\"]*\")?\)", prose):
        dest = dest.strip("<>")
        url = urlsplit(dest)
        if url.scheme or url.netloc or not url.path:
            continue
        target = ((ROOT / name).parent / unquote(url.path)).resolve()
        try:
            rel = target.relative_to(ROOT).as_posix()
        except ValueError:
            errors.append(f"{name}: link leaves repository: {dest}")
            continue
        exists = rel in tracked or any(p.startswith(rel.rstrip('/') + '/') for p in tracked)
        if not exists:
            errors.append(f"{name}: link is not in the public tree: {dest}")
    for block in re.findall(r"^```json\s*\n(.*?)^```", body, flags=re.M | re.S):
        try:
            json.loads(block)
        except json.JSONDecodeError as error:
            errors.append(f"{name}: invalid JSON example: {error}")

if errors:
    print("\n".join(errors), file=sys.stderr)
    sys.exit(1)
print(f"Public-tree and documentation checks passed ({len(tracked)} tracked files).")
