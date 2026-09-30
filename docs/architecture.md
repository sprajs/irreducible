# Architecture

Irreducible supplies compiled physics and model calculations through a non-interactive CLI and an independently usable C++ library. Agents handle the investigation: they can compose requests, drive external fitting tools and extend the source with new readers or models.

The core shares equations and conventions across models. It needs statistical or fitting kernels only when a concrete calculation requires them; providing every analysis method is not an architectural goal.

| Layer | Responsibility | Source |
| --- | --- | --- |
| Rust control/data layer | Commands, structural parsing, safe wrappers, data movement, run records | `src/` |
| Versioned C ABI | Coarse batches, explicit buffers/ownership, status propagation | [schema/abi.json](../schema/abi.json), generated bindings |
| C++20 core | Physical definitions, conversions, numerical methods and developing models | `cpp/include/irred/`, `cpp/src/` |
| Verification | Mathematical fixtures, independent algorithms, domain and discrepancy regressions | `cpp/tests/`, `tests/` |

`tools/generate_abi.py` generates bindings from one shared schema. Change the schema and regenerate rather than hand-maintaining competing enum/layout definitions. Catch exceptions on the C++ side and translate status explicitly. Rust wrappers enforce lifetime/buffer rules. No exception may unwind across the boundary.

Equations belong in compiled source with explicit model/equation identities. Agents propose changes through normal development; data files cannot inject executable physics. Shared definitions avoid subtly different copies of the same equation across analyses.

`tools/build.py` owns the integrated build. Cargo links the native static library; CMake builds C++ and its tests. Manifests make toolchain and input identities inspectable. The conservative Debug profile is an engineering baseline, not a release-speed benchmark.

A future accelerator implements an explicit numerical contract. It must preserve a portable baseline and be compared at matched quality, including setup/transfer costs and downstream error. A backend label is not qualification.

Generic C++ modules use `irred` (`irred::numerics`, `irred::observations`) and install under `include/irred/`; the standalone archive is `libirred_core.a`. Reserve `irred::cosmology` for an actual physical-model module. C ABI names are structural transport identities, not a promise of indefinite compatibility. Coherent breaking migrations may replace routes/layouts and their callers; immutable scientific records retain the exact contracts they used.

Operation wrappers return a typed outcome with output IDs, actual method/arithmetic, numerical checks and resource metadata. The common recorder records that outcome without classifying operation-name families. Retained native factors already serve batches; duplicated acquisition/preparation and request-level ownership are separate migration targets.

The native observation handle stores one immutable shared observation object. Internal retained-source access keeps that same object alive after the acquisition handle is released; source storage remains alive until its last owner releases it. The generic `typed_magnitude_covariance` source profile admits explicitly declared measured values, fitted summaries or synthetic controls with magnitude covariance axes. It does not apply a Pantheon cut, verify a release identity or upgrade calibration/dependence declarations. The narrow ASCII adapter uses explicit `ROW_ID EVENT_ID MAG` columns and supplied covariance ordering.
