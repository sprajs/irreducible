# Testing

The native tests run from a source checkout using the fixtures included in the repository. External comparison software is not required.

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

## Optional original-input supernova comparison

The native consumer and independent synthetic/LDLT checks run in ordinary tests. The original released-input comparison is explicit and potentially expensive; acquire the exact assets pinned by `cpp/tests/fixtures/w01_historical.hpp` before invoking it. No dataset is bundled. Build its optional target in the chosen profile:

```sh
.build-tools/bin/cmake --build build/native --target test_w01_scores --parallel 4
```

Set `IRRED_W01_TABLE`, `IRRED_W01_COVARIANCE` and `IRRED_W01_NATIVE_HARNESS` to the original files and built `test_w01_scores` executable, then run:

```sh
cargo test --locked --offline -j4 --test w01_reference -- --ignored --nocapture
```

The Rust wrapper verifies raw hashes before running C++; it contains no physical equations. The native harness compares the frozen direct/stable historical scores separately from the approximate compressed route. This is a named fixed-input comparison, not parameter fitting or cosmological inference. The implemented supernova interface has separate native and CLI transport gates.

The retained supernova interface has a durable native allocation/ownership suite (`test_supernova_abi_hostile`) and Rust CLI suite (`cargo test --locked --offline -j4 --test supernova_cli`). Its generated transport controls use analytic de Sitter/profile limits and intentionally invalid unselected covariance entries. They are separate from the optional SHA-verified original seven-point regression. Native Debug and Release suites exercise the same checks; the interface also has actual Release CLI coverage. A passing transport regression does not qualify arbitrary source/model requests.

For the optional native BAO eleven-point comparison, explicitly acquire the assets pinned by `cpp/tests/fixtures/bao_reference.hpp`. The guard rejects changed bytes and checks them again after execution:

```bash
IRRED_BAO_MEAN=/path/to/desi_gaussian_bao_ALL_GCcomb_mean.txt \
IRRED_BAO_COVARIANCE=/path/to/desi_gaussian_bao_ALL_GCcomb_cov.txt \
IRRED_BAO_NATIVE_HARNESS="$PWD/build/native-release/test_bao" \
cargo test --release --locked --offline -j4 --test bao_reference -- --ignored --nocapture
```

This is native fixed-case numerical coverage; it supplies no BAO CLI, posterior or joint-probe qualification.

Explicit mixed-model v2 transport has native `test_cpl_v2_abi_hostile` and Rust `cpl_v2_cli` suites. For the optional original four-point CPL CLI comparison, set the same exact `IRRED_W01_TABLE` and `IRRED_W01_COVARIANCE` assets and run:

```sh
cargo test --release --locked --offline -j4 --test w01_reference released_assets_and_cli_cpl_four_points -- --ignored --nocapture
```

This verifies source hashes before/after execution, full attempted parameter identities and immutable result records against the pinned native relative-profile fixture. It contains no Rust physical equations and does not qualify arbitrary requests.

## Original BAO CLI transport check

With the exact assets pinned in `cpp/tests/fixtures/bao_reference.hpp`, set `IRRED_BAO_MEAN` and `IRRED_BAO_COVARIANCE`, then run:

```sh
cargo test --release --locked --offline -j4 --test bao_reference released_assets_and_cli_eleven_points -- --ignored --nocapture
```

The guard verifies raw hashes before and after execution and compares eleven ordered outputs to the independent GL8/long-double LDLT fixture. Exact direct native/ABI parity is a separate ordinary hostile test. Optional `IRRED_BAO_CLI_RECORD_DIRECTORY` captures the executable and discovery before execution plus the request, output and immutable store. No dataset is bundled or required by ordinary CI.
