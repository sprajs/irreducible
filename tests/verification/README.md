# Independent native verification

Accepted S00 cases live in `cpp/tests/test_abi_hostile.cpp`, `tests/cli.rs` and
`tests/build_identity.rs`. Temporary Python bridge/record/build-mutation drivers
were removed after their native equivalents passed. The exact rational covariance
oracle generator remains pending S04 migration; it imports no production code.
Historical supernova assets are read-only; W01 byte identities are recorded in
`Plan/08-delivery/implementation/verifier-w01-inventory.json`.

## Actual qualification commands

From the repository root, using no other concurrent build/mutation worker:

```sh
python3 tools/build.py
.build-tools/bin/ctest --test-dir build/native --output-on-failure
cargo test --locked --offline -j4
cargo test --locked --offline -j4 --test build_identity --no-run
```

The last command prints the compiled `tests/build_identity.rs` executable. Run
that exact executable directly with `--ignored --exact
source_receipt_flags_and_cache_identity`. This avoids nesting Cargo while testing
rebuilds. The actual 2026-09-30 command was:

```sh
target/debug/deps/build_identity-77c4b34e31ec65b9 --ignored --exact source_receipt_flags_and_cache_identity
```

Ordinary `cargo test` explicitly ignores that exclusive mutation test and cannot
establish its gate. The test temporarily edits Rust/C++ source comments, restores
original bytes on unwind, verifies rebuilt identity, excludes mutable receipts,
rejects declared unsupported overrides and poisons then resets CMake flags.

The S00 qualification receipt pins the exact build, commands, binary digests,
accepted scope and remaining limits. Failed identity findings and original/retest
receipts are preserved independently. No cosmological capability follows from
these engineering checks. The current layout fixture applies to 64-bit Linux;
other targets require their own layout and execution evidence.
