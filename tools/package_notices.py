"""Preserve upstream notices for the resolved runtime graph and Rust library."""
import json
from pathlib import Path
import shutil
import subprocess
from ci_identity import digest


def collect(root, destination):
    root, destination = Path(root), Path(destination)
    version = subprocess.check_output(["rustc", "-vV"], text=True)
    host = next(line.split(": ", 1)[1] for line in version.splitlines() if line.startswith("host:"))
    metadata = json.loads(subprocess.check_output(
        ["cargo", "metadata", "--locked", "--offline", "--format-version", "1", "--filter-platform", host],
        cwd=root, text=True))
    packages = {p["id"]: p for p in metadata["packages"]}
    nodes = {p["id"]: p for p in metadata["resolve"]["nodes"]}
    seen, pending = set(), [metadata["resolve"]["root"]]
    while pending:
        current = pending.pop()
        if current in seen:
            continue
        seen.add(current)
        for dep in nodes[current]["deps"]:
            package = packages[dep["pkg"]]
            if any("proc-macro" in t["kind"] for t in package["targets"]):
                continue
            if any(k["kind"] is None for k in dep["dep_kinds"]):
                pending.append(dep["pkg"])
    inventory = []
    for key in sorted(seen - {metadata["resolve"]["root"]}):
        package = packages[key]
        directory = Path(package["manifest_path"]).parent
        candidates = [p for p in directory.iterdir() if p.is_file() and
                      any(label in p.name.upper() for label in ("LICENSE", "COPYRIGHT", "NOTICE", "UNLICENSE"))]
        if package.get("license_file"):
            candidates.append(directory / package["license_file"])
        if not candidates:
            raise ValueError(f"missing upstream notices: {package['name']}")
        target = destination / f"{package['name']}-{package['version']}"
        target.mkdir(parents=True)
        copied = {}
        for source in sorted(set(candidates)):
            shutil.copyfile(source, target / source.name)
            copied[source.name] = digest(source)
        inventory.append({"name": package["name"], "version": package["version"],
                          "declared_license": package.get("license"),
                          "source": package.get("source"), "notice_sha256": copied})
    sysroot = Path(subprocess.check_output(["rustc", "--print", "sysroot"], text=True).strip())
    docs = sysroot / "share/doc/rust"
    library_notice = docs / "COPYRIGHT-library.html"
    if not library_notice.is_file():
        raise ValueError("Rust standard-library notice missing; install the pinned rust-docs component")
    target = destination / "rust-standard-library"
    target.mkdir()
    shutil.copyfile(library_notice, target / library_notice.name)
    licenses = docs / "licenses"
    if not licenses.is_dir():
        raise ValueError("Rust notice licence texts unavailable")
    shutil.copytree(licenses, target / "licenses")
    std_hashes = {str(p.relative_to(target)): digest(p) for p in target.rglob("*") if p.is_file()}
    result = {"runtime_crates": inventory, "target": host, "rustc": version,
              "rust_library_notice_sha256": std_hashes,
              "scope": "Resolved normal runtime edges; build/dev/proc-macro-only crates excluded. Faithful installed Rust library notice covers potential library components, not a symbol-level linkage audit. Dynamic system libraries are not bundled."}
    (destination / "inventory.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result
