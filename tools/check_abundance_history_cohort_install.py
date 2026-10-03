#!/usr/bin/env python3
"""Install and check the fixed two-state public SDK proof caller.

This invokes compilation and scientific execution. Its caller must first obtain
the coordinated compute lease. Complete reference qualification is withheld.
"""
import argparse
import hashlib
import json
import pathlib
import shlex
import shutil
import subprocess
import tempfile

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--profile", choices=["debug", "release"], default="debug")
parser.add_argument("--native", type=pathlib.Path, help="Explicit fresh matching native build directory")
parser.add_argument("--require-library-no-elision", action="store_true")
parser.add_argument("--output-receipt", type=pathlib.Path)
args = parser.parse_args()
root = pathlib.Path(__file__).resolve().parents[1]
native = args.native.resolve() if args.native else root / ("build/native" if args.profile == "debug" else "build/native-release")
cache_file = native / "CMakeCache.txt"
cache = {}
if cache_file.is_file():
    for line in cache_file.read_text().splitlines():
        if line and not line.startswith(("#", "//")) and ":" in line and "=" in line:
            key, value = line.split("=", 1)
            cache[key.split(":", 1)[0]] = value
compiler = cache.get("CMAKE_CXX_COMPILER") or shutil.which("c++")
cmake = root / ".build-tools/bin/cmake"
if not compiler or not cmake.is_file():
    raise SystemExit("Configured native compiler/CMake is unavailable")

# Checking metadata alone does not certify the actual compiled object. The exact
# build/source/runtime receipt must also be reviewed by the integration owner.
commands_file = native / "compile_commands.json"
core_commands = []
if commands_file.is_file():
    for entry in json.loads(commands_file.read_text()):
        command = entry.get("arguments") or shlex.split(entry.get("command", ""))
        if any("irred_core.dir" in token for token in command):
            core_commands.append(command)
def requests_no_elision(command):
    flags = [token for token in command if token in ("-fno-elide-constructors", "-felide-constructors")]
    return bool(flags) and flags[-1] == "-fno-elide-constructors"


library_no_elision = bool(core_commands) and all(requests_no_elision(c) for c in core_commands)
if args.require_library_no_elision and not library_no_elision:
    raise SystemExit("Library no-elision compile metadata is missing; complete gate withheld")

members = [root / "cpp/tests/test_installed_abundance_history_cohort.cpp",
           root / "cpp/tests/abundance_history_cohort.hpp",
           root / "cpp/tests/abundance_history_cohort_source.hpp"]
receipt = {
    "role": "fixed_synthetic_two_half_exact_emitted_working_input",
    "profile": args.profile,
    "native": str(native),
    "native_source_directory": cache.get("CMAKE_HOME_DIRECTORY"),
    "native_build_type": cache.get("CMAKE_BUILD_TYPE"),
    "caller_no_elision": True,
    "library_no_elision_compile_metadata": library_no_elision,
    "complete_reference_qualified": False,
    "physical_joint_law_qualified": False,
    "source_sha256": {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest() for p in members},
    "run_status": "not_run",
    "stages": [],
}


def run_stage(name, command):
    result = subprocess.run(command, text=True, capture_output=True)
    receipt["stages"].append({"name": name, "command": command, "exit_code": result.returncode,
                              "stdout": result.stdout, "stderr": result.stderr})
    print(result.stdout, end="")
    if result.stderr:
        print(result.stderr, end="")
    if result.returncode:
        raise RuntimeError(f"{name} refused with exit code {result.returncode}")
    return result


try:
    if cache.get("CMAKE_HOME_DIRECTORY") and pathlib.Path(cache["CMAKE_HOME_DIRECTORY"]).resolve() != root / "cpp":
        raise RuntimeError("Native build is from another source checkout; exact source gate withheld")
    identity = subprocess.run(["git", "rev-parse", "HEAD"], cwd=root, text=True, capture_output=True, check=True)
    receipt["source_head"] = identity.stdout.strip()
    with tempfile.TemporaryDirectory(prefix="irred-abundance-history-install-") as temporary:
        scratch = pathlib.Path(temporary)
        prefix = scratch / "prefix"
        run_stage("install", [str(cmake), "--install", str(native), "--prefix", str(prefix)])
        installed_header = prefix / "include/irred/hydrogen_helium_history.hpp"
        if not installed_header.is_file():
            raise RuntimeError("Missing installed public history header")
        receipt["installed_history_header_sha256"] = hashlib.sha256(installed_header.read_bytes()).hexdigest()
        output = scratch / "cohort-consumer"
        command = [compiler, "-std=c++20", "-Wall", "-Wextra", "-Wpedantic", "-fno-fast-math",
                   "-ffp-contract=off", "-fno-elide-constructors", str(members[0]),
                   "-I", str(prefix / "include"), str(prefix / "lib/libirred_core.a"), "-o", str(output)]
        run_stage("caller_compile", command)
        receipt["installed_library_sha256"] = hashlib.sha256((prefix / "lib/libirred_core.a").read_bytes()).hexdigest()
        receipt["caller_sha256"] = hashlib.sha256(output.read_bytes()).hexdigest()
        result = run_stage("caller_execution", [str(output)])
        receipt["run_status"] = "passed"
        receipt["exit_code"] = result.returncode
except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as error:
    receipt["run_status"] = "failed"
    receipt["failure"] = str(error)
    raise SystemExit(str(error)) from error
finally:
    if args.output_receipt:
        args.output_receipt.parent.mkdir(parents=True, exist_ok=True)
        args.output_receipt.write_text(json.dumps(receipt, indent=2) + "\n")
print("Fresh installed cohort caller passed; complete reference and physical law remain withheld")
