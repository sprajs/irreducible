# Getting started

Irreducible is currently built from source. The executable is `irred`; there is no packaged release or supported installer yet.

## Arithmetic on Apple Silicon

Apple Silicon uses binary64 for both `double` and `long double`. The existing
`F02/longdouble-cpu/v1` owner requires at least 64 mantissa bits and a wider
exponent range, so that owner refuses on this representation. The distinct
binary64 owner remains available; changing arithmetic is a new calculation
identity, not an automatic substitute for the qualified wide calculation.
The draft ideal acoustic model also requires wide arithmetic.

The Apple Silicon CI check records its actual Apple Clang version and
architecture, compiles the production numerical owner, tests an exact analytic
SPD solve with binary64, and verifies that unsupported wide preparation refuses.
It is an engineering arithmetic check, not a full macOS CLI/SDK build, a
MacBook execution receipt, or scientific qualification of binary64 consumers.
The same test runs in the ordinary Linux native suite.

## Prepare a Linux checkout

You need Git, Python 3.11 or newer with `venv`, Rust/Cargo with Rust 2024 edition support and the `rustfmt` component, and a C++20 compiler/linker. The locally tested development environment uses Rust/Cargo 1.98.0, GCC 16.2.1 and CMake 4.1.3 on 64-bit Linux. The qualified wide-arithmetic suites require long double with at least 64 mantissa bits; targets where long double equals double are not qualified by these gates. Ubuntu CI is a separate engineering portability check, not inherited scientific validation. Other toolchains require their own validation; these versions describe the tested host, not a proven minimum-version matrix.

```sh
git clone https://github.com/sprajs/irreducible.git
cd irreducible
rustup component add rustfmt
python3 -m venv .build-tools
.build-tools/bin/python -m pip install cmake==4.1.3
cargo fetch --locked
python3 tools/build.py
```

Initial dependency preparation uses the network. Subsequent Cargo build/test commands use the lockfile and offline cache. If offline resolution fails, run `cargo fetch --locked` with network access; do not silently change dependency versions.

The build driver expects `.build-tools/bin/cmake`. It generates the ABI bindings, configures the native Debug build, records build identity, builds C++ with four jobs, and builds Rust against that library. It deliberately rejects undeclared compiler flags, target/profile overrides and Cargo config overrides in the checkout, ancestor directories or Cargo home. Start from a clean shell configuration if an override is rejected. See [development](development.md) before changing the build profile.

## Inspect and run

```sh
target/debug/irred version --json
target/debug/irred describe --json
target/debug/irred run tests/fixtures/exact-add.json runs/first-run
```

The exact fixture returns `[0, 2, -4, 9223372036854775807]`, with a receipt. This verifies basic parsing, the Rust/C++ boundary and local record writing. It does not test cosmology.

```sh
target/debug/irred run tests/fixtures/quantity-length.json runs/quantity-example
```

The quantity example preserves source metadata and converts metre values to kilometres. A result whose required numerical checks pass exits **0** by default, with acceptance scoped to `numerical_contract`. Explicit `--assurance qualified` exits **6** when applicable evidence is absent; interpretation remains unqualified. See [CLI status semantics](cli.md).

## Verify the build

```sh
.build-tools/bin/ctest --test-dir build/native --output-on-failure
cargo test --locked --offline -j4
python3 tools/check_docs.py
```

The source/flag mutation test is ignored by ordinary Cargo test runs and has a separate exclusive procedure in [testing](testing.md). Coordinate shared builds; do not run a second build or a source-mutation test in an active build tree.

`build/`, `target/`, `.build-tools/` and `runs/` are generated locally.

For optimized binaries, use the separately identified [Release profile](testing.md#explicit-release-profile). Debug remains the default; build both profiles before comparing performance at matched numerical quality.
