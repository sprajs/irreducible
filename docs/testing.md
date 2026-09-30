# Testing

Tests are the durable verification interface. They must run from a public source checkout without local plans, comparison receipts or external astronomy/oracle environments.

## Ordinary checks

After [dependency preparation](getting-started.md), run:

```sh
python3 tools/build.py
.build-tools/bin/ctest --test-dir build/native --output-on-failure
cargo test --locked --offline -j4
python3 tools/check_docs.py
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
