# Irreducible

Irreducible is an agent-operated scientific engine with a Rust control/data layer and an owned C++20 numerical and scientific library. Its executable is `irred`.

The implemented CPU-only build has passed the S00 infrastructure gate and the bounded S01 unit-scale/metadata/ABI gate. A narrow S02 numerical library is implemented and under staged validation; its scalar/reduction CLI bridge is pending independent qualification. Scientific CLI calculations currently remain unaccepted with exit 6 until applicable exact-build numerical evidence is registered. No cosmological inference or historical research reproduction is claimed. CUDA and distributed execution remain planned.

Build and inspect the actual executable using the [tested build instructions](Plan/08-delivery/implementation/BUILD.md):

```bash
python3 tools/build.py
.build-tools/bin/ctest --test-dir build/native --output-on-failure
cargo test --locked --offline -j4
target/debug/irred describe --json
target/debug/irred run tests/fixtures/exact-add.json runs/example
```

- [Main plan](Plan/README.md) explains the design and topic owners.
- [Build order](Plan/BUILD-ORDER.md) defines implementation dependencies and gates.
- [Status](Plan/STATUS.md) records current evidence and limitations.
- [Agent instructions](AGENTS.md) cover development, analysis and validated publication.

The earlier supernova research supplies cases, identified assets and known limitations. Its outputs are not automatically independent validation or results of Irreducible. Public source: [sprajs/irreducible](https://github.com/sprajs/irreducible).
