# Architecture

Irreducible has one executable and an independently usable native scientific library.

| Layer | Responsibility | Source |
| --- | --- | --- |
| Rust control/data layer | Commands, structural parsing, safe wrappers, data movement, run records | `src/` |
| Versioned C ABI | Coarse batches, explicit buffers/ownership, status propagation | [schema/abi.json](../schema/abi.json), generated bindings |
| C++20 core | Physical definitions, conversions, numerical methods and developing models | `cpp/include/cosmology/`, `cpp/src/` |
| Verification | Mathematical fixtures, independent algorithms, domain and discrepancy regressions | `cpp/tests/`, `tests/` |

`tools/generate_abi.py` generates bindings from one shared schema. Change the schema and regenerate rather than hand-maintaining competing enum/layout definitions. Catch exceptions on the C++ side and translate status explicitly. Rust wrappers enforce lifetime/buffer rules. No exception may unwind across the boundary.

Equations belong in compiled source with explicit model/equation identities. Agents propose changes through normal development; data files cannot inject executable physics. Shared definitions avoid subtly different copies of the same equation across analyses.

`tools/build.py` owns the integrated build. Cargo links the native static library; CMake builds C++ and its tests. Manifests make toolchain and input identities inspectable. The conservative Debug profile is an engineering baseline, not a release-speed benchmark.

A future accelerator implements an explicit numerical contract. It must preserve a portable baseline and be compared at matched quality, including setup/transfer costs and downstream error. A backend label is not qualification.
