"""Source version and CI provenance; execution labels never replace build_id."""
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import tomllib


def source_version(root):
    return tomllib.loads((Path(root) / "Cargo.toml").read_text())["package"]["version"]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def ci_metadata(root, manifest, discovery, environment, event, checked_sha):
    version = source_version(root)
    if manifest.get("source_version") != version or discovery.get("version") != version:
        raise ValueError("source/manifest/executable version mismatch")
    if discovery["build"]["build_id"] != manifest["build_id"]:
        raise ValueError("executable build identity mismatch")
    if manifest.get("git_head") != checked_sha or environment["GITHUB_SHA"] != checked_sha:
        raise ValueError("checked-out commit mismatch")
    if not re.fullmatch(r"[0-9a-f]{40}", checked_sha):
        raise ValueError("invalid checked-out SHA")
    run = environment["GITHUB_RUN_ID"]
    attempt = environment["GITHUB_RUN_ATTEMPT"]
    if not run.isdecimal() or not attempt.isdecimal():
        raise ValueError("invalid CI run identity")
    event_name = environment["GITHUB_EVENT_NAME"]
    pr = event.get("pull_request")
    if event_name not in ("pull_request", "push", "workflow_dispatch") or (event_name == "pull_request") != bool(pr):
        raise ValueError("CI event/payload mismatch")
    if pr and any(not re.fullmatch(r"[0-9a-f]{40}", pr[k]["sha"]) for k in ("head", "base")):
        raise ValueError("invalid PR source SHA")
    return {"schema_version": 1, "source_version": version,
            "ci_identity": f"run-{run}.attempt-{attempt}.sha-{checked_sha}",
            "run_id": run, "run_attempt": attempt,
            "event_name": environment["GITHUB_EVENT_NAME"],
            "checked_out_sha": checked_sha,
            "pr_head_sha": pr["head"]["sha"] if pr else None,
            "pr_base_sha": pr["base"]["sha"] if pr else None,
            "checkout_role": "synthetic_pr_merge" if pr else "source_commit",
            "build_id": manifest["build_id"],
            "platform_scope": "Ubuntu 24.04 x86_64; runtime-library compatibility required"}


def write_metadata(root, destination):
    root = Path(root)
    manifest = json.loads((root / "build/build-manifest-release.json").read_text())
    discovery = json.loads(subprocess.check_output([str(root / "target/release/irred"),
                                                  "describe", "--json"], text=True))
    event = json.loads(Path(os.environ["GITHUB_EVENT_PATH"]).read_text())
    sha = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip()
    value = ci_metadata(root, manifest, discovery, os.environ, event, sha)
    Path(destination).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    return value
