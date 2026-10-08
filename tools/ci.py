#!/usr/bin/env python3
"""Execute separated CI suites and construct a scoped, verified build artifact."""
import argparse
import json
from pathlib import Path
import shutil
import platform
import subprocess
import tarfile
from ci_identity import digest, source_version, write_metadata
from package_notices import collect
from ci_suites import native_build_targets, classify_native, native_regex, cli_targets

ROOT = Path(__file__).resolve().parents[1]


def run(*args):
    subprocess.run(args, cwd=ROOT, check=True)


def targets():
    metadata = json.loads(subprocess.check_output(
        ["cargo", "metadata", "--format-version", "1", "--no-deps", "--locked", "--offline"],
        cwd=ROOT, text=True))
    return next(p["targets"] for p in metadata["packages"]
                if Path(p["manifest_path"]).resolve() == ROOT / "Cargo.toml")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=["native", "rust-unit", "cli", "package"])
    parser.add_argument("--profile", choices=["debug", "release"], default="debug")
    parser.add_argument("--compiler", choices=["gcc", "clang"], default="gcc")
    parser.add_argument("--suite", choices=["fast", "full"], default="fast",
                        help="fast engineering checks; full includes scientific regressions")
    parser.add_argument("--jobs", type=int, choices=range(1, 5), default=2)
    args = parser.parse_args()
    cmake = str(ROOT / ".build-tools/bin/cmake")
    if args.action == "native":
        compiler = shutil.which("g++" if args.compiler == "gcc" else "clang++")
        if not compiler:
            raise SystemExit("declared native compiler unavailable")
        build = f"build/ci-native-{args.compiler}"
        run("python3", "tools/generate_abi.py")
        run(cmake, "-S", "cpp", "-B", build, "-DCMAKE_BUILD_TYPE=Debug",
            f"-DCMAKE_CXX_COMPILER={compiler}", "-DCMAKE_CXX_FLAGS=",
            "-DIRRED_TEST_CFITSIO=OFF", "-DIRRED_REFERENCE_GMP_MPFR=OFF")
        inventory = json.loads(subprocess.check_output(
            [str(ROOT / ".build-tools/bin/ctest"), "--test-dir", build, "--show-only=json-v1"], text=True))
        print("Native suite inventory:", json.dumps(classify_native(inventory), sort_keys=True), flush=True)
        selected = native_build_targets(args.suite)
        run(cmake, "--build", build, "--parallel", str(args.jobs),
            *(["--target", *selected] if selected else []))
        run(str(ROOT / ".build-tools/bin/ctest"), "--test-dir", build,
            "--output-on-failure", "--no-tests=error", "-j1",
            *(["-R", native_regex()] if args.suite == "fast" else []))
    elif args.action == "rust-unit":
        unit_flags = [item for t in targets() if "bin" in t["kind"] for item in ("--bin", t["name"])]
        has_library = any("lib" in t["kind"] for t in targets())
        if has_library:
            unit_flags.append("--lib")
        if not unit_flags:
            raise SystemExit("no unit-test targets found")
        run("cargo", "test", "--locked", "--offline", f"-j{args.jobs}",
            *(["--release"] if args.profile == "release" else []), *unit_flags,
            "--", "--test-threads=1")
        if has_library:
            run("cargo", "test", "--locked", "--offline", f"-j{args.jobs}",
                *(["--release"] if args.profile == "release" else []), "--doc",
                "--", "--test-threads=1")
    elif args.action == "cli":
        integration_targets = sorted(t["name"] for t in targets() if "test" in t["kind"])
        if not integration_targets:
            raise SystemExit("no integration targets found")
        integration_targets = cli_targets(integration_targets, args.suite)
        print(f"Integration suite {args.suite} targets:", ", ".join(integration_targets), flush=True)
        run("cargo", "test", "--locked", "--offline", f"-j{args.jobs}",
            *(["--release"] if args.profile == "release" else []),
            *[item for name in integration_targets for item in ("--test", name)],
            "--", "--test-threads=1")
    else:
        system = platform.freedesktop_os_release()
        if system.get("ID") != "ubuntu" or system.get("VERSION_ID") != "24.04" or platform.machine() != "x86_64":
            raise SystemExit("package requires actual Ubuntu 24.04 x86_64 host")
        destination = ROOT / "build/ci-package"
        if destination.exists():
            shutil.rmtree(destination)
        destination.mkdir(parents=True)
        prefix = destination / "irred"
        run(cmake, "--install", "build/native-release", "--prefix", str(prefix))
        (prefix / "bin").mkdir()
        shutil.copyfile(ROOT / "target/release/irred", prefix / "bin/irred")
        (prefix / "bin/irred").chmod(0o755)
        shutil.copyfile(ROOT / "LICENSE", prefix / "LICENSE")
        shutil.copyfile(ROOT / "build/build-manifest-release.json", prefix / "build-manifest.json")
        metadata = write_metadata(ROOT, prefix / "ci-metadata.json")
        collect(ROOT, prefix / "third-party-notices")
        (prefix / "README.txt").write_text(
            "Irreducible CI build artifact, not a tagged release or universal installer.\n"
            "Built on Ubuntu 24.04 x86_64; compatible dynamic system libraries are required.\n"
            "Run bin/irred describe --json or bin/irred version --json.\n"
            "Cargo source version is deliberate; CI run/attempt/SHA identify this execution; build_id identifies the recorded source/tool/flags content.\n"
            "Verify SHA256SUMS before use. lib/ and include/ form the standalone native SDK.\n"
            "Upstream licence declarations/notices and hashes are in third-party-notices/.\n")
        for name in ("lib/libirred_core.a", "include/irred/abi.h", "bin/irred"):
            if not (prefix / name).is_file() or not (prefix / name).stat().st_size:
                raise SystemExit(f"missing package payload: {name}")
        files = sorted(p for p in prefix.rglob("*") if p.is_file())
        (prefix / "SHA256SUMS").write_text("".join(
            f"{digest(p)}  {p.relative_to(prefix)}\n" for p in files))
        name = f"irred-{source_version(ROOT)}-ubuntu24.04-x86_64-run{metadata['run_id']}-attempt{metadata['run_attempt']}.tar.gz"
        archive = destination / name
        with tarfile.open(archive, "w:gz") as package:
            package.add(prefix, arcname="irred")
        (destination / (name + ".sha256")).write_text(f"{digest(archive)}  {name}\n")
        print(archive)


if __name__ == "__main__":
    main()
