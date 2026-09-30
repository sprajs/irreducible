# Testing

The native tests run from a source checkout using the fixtures included in the repository. External comparison software is not required.

## Optional directed interval reference

The test-only `IRRED_REFERENCE_GMP_MPFR` profile requires installed dynamic GMP/GMPXX and MPFR libraries. The integrated product build forces it OFF; these libraries are not production dependencies. Missing headers or dynamic libraries make this explicit profile unavailable. No libraries or original datasets are bundled.

The following separate reference build and no-assets controls have been executed with GMP/GMPXX 6.3.0 and MPFR 4.2.2:

```sh
.build-tools/bin/cmake -S cpp -B /tmp/irred-reference-profile-check -DCMAKE_BUILD_TYPE=Release -DIRRED_REFERENCE_GMP_MPFR=ON
.build-tools/bin/cmake --build /tmp/irred-reference-profile-check --target test_bao_piecewise_interval --parallel 1
.build-tools/bin/ctest --test-dir /tmp/irred-reference-profile-check -R '^bao_piecewise_interval_controls$' --output-on-failure
```

The executable target is excluded from the default build, so build it explicitly before running its control test. Controls challenge exact rational and directed interval arithmetic; they do not establish an original-data scientific gate. Original-data mode requires an external SHA guard and separate component-budget review. Pin actual library/header versions and hashes, compiler flags, reference source and raw input identities for that comparison. The current host's MPFR headers declare LGPL 3 or later; GMPXX headers declare LGPL 3 or later or GPL 2 or later. Review the installed notices for other enabled environments.

## Ordinary checks

After [dependency preparation](getting-started.md), run:

```sh
python3 tools/build.py
.build-tools/bin/ctest --test-dir build/native --output-on-failure
cargo test --locked --offline -j4
python3 tools/check_docs.py
python3 tools/check_install.py
```

C++ tests exercise scientific/numerical contracts and hostile ABI cases. Rust tests cover safe wrappers, request parsing, records and the CLI. `check_docs.py` checks tracked publication paths and relative Markdown file links; it does not validate equations, external websites or every Markdown feature.

## Exclusive build-identity test

This test changes source comments and CMake flags, then rebuilds and restores them. Run it in an isolated checkout, with no other build using that tree:

```sh
cargo test --locked --offline -j4 --test build_identity --no-run
```

Cargo prints the compiled test executable path. Run **that exact path**, followed by:

```text
--ignored --exact source_receipt_flags_and_cache_identity
```

Do not launch it through a nested Cargo test invocation. Ordinary Cargo testing intentionally ignores it. It checks that source changes affect build identity, mutable local files do not, unsupported overrides fail, and poisoned CMake flags are reset. Preserve the original checkout if the process is forcibly terminated before restoration.

## Fixture ancestry

[Native fixture provenance](../cpp/tests/fixtures/README.md) records the essential constants, derivations and budgets. Keep small explicit expected values and independent algorithms with the tests. Preserve the distinction between exact algebra, high-precision comparisons and shared-library checks.

Full external-tool output and development receipts stay local. Once a comparison becomes a durable regression, routine tests should not require rerunning the external software. Do not delete original historical research or accepted records merely to make the public tree smaller.

A passing test is scoped evidence. It does not qualify arbitrary inputs, every supported toolchain, cosmological inference or a performance claim. Quadrature missing a narrow unsampled feature is a known limitation; a finite estimate is not a proof of accuracy.

## Observation and optional codec tests

Ordinary native tests include immutable observation ownership, exact masks/row order, uncertainty asymmetry retention, hostile descriptors and scoped allocation-failure cleanup. Rust ingestion/CLI tests retain raw source bytes, challenge same-path mutations, semantic failures and nonfinite/missing distinctions. They use small synthetic inputs and require no acquired survey files.

The product build deliberately reports FITS unavailable. A separate optional native codec profile requires installed CFITSIO and can be tested with `-DCOSMOLOGY_TEST_CFITSIO=ON` in a separate CMake build directory. It covers one handcrafted synthetic BINTABLE profile, not universal FITS support. The wrapper resets this test-only option to OFF for its supported product build. CFITSIO is file-format infrastructure; its installed version/library identity and NASA permissive notice must be reviewed for any enabled deployment.

`check_install.py` installs to a fresh temporary prefix, rejects the old generic include root, and compiles/runs the durable `test_installed_consumer.cpp` against the installed `libirred_core.a`. It checks independent library usability and cleans its temporary installation.

## Explicit Release profile

Build and test the optimized profile with the same native and Rust cases:

```sh
python3 tools/build.py --profile release
.build-tools/bin/ctest --test-dir build/native-release --output-on-failure
cargo test --release --locked --offline -j4
target/release/irred describe --json
```

Debug remains the default and uses `build/native`, `target/debug` and `build/build-manifest.json`. Release uses `build/native-release`, `target/release` and `build/build-manifest-release.json`. Cargo selects the matching native archive and embeds the matching manifest; do not manually cross-link profiles. Both retain strict floating-point options and panic containment. Release adds C++ `-O3 -DNDEBUG` and the pinned standard Cargo release settings; its source/build/profile identity differs from Debug. Native checks remain active under NDEBUG. These build checks establish no workload speedup.

After building both profiles in an isolated checkout, run the durable identity/negative witness:

```sh
cargo test --locked --offline -j4 --test profile_identity -- --ignored
```

This checks each executable's describe/version manifest, identical source digests with distinct profile/build/native-archive identities, and a deliberately wrong native fixture compiled with NDEBUG. That fixture must fail, proving its required checks were executed. It compiles one temporary native test, cleans it after success, and never mutates production source.

## Current native and public regression suites

Ordinary CTest includes the current expansion, retained supernova, BAO and Gaussian boundary suites, independent analytic/quadrature/cofactor controls, ownership/allocation failures and payload accounting. Ordinary Cargo tests include current one-shot and retained-stream requests, source identities, quota failures and immutable run records. Request schema is 2; physical operation identifiers can retain their `.v1` suffix. Scientific qualification remains bounded to named tested cases.

`test_w01_reference` is an optional independent LDLT condition/reconstruction control, not a retired consumer interface. Existing native fixtures retain historical direct/stable and compressed-route discrepancies separately. No expected values are regenerated during builds.
