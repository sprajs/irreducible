# Capabilities and limits

This page describes the public command surface at this documentation revision. Development moves quickly: `describe --json` and source at the exact commit take precedence. A native C++ function is not automatically a CLI operation.

| Surface | Available now | Limit |
| --- | --- | --- |
| `describe --json`, `version --json` | Product, build, ABI/schema, commands and capability metadata | No automatic scientific qualification |
| `fixture.checked_i64_add.v1` | Checked integer batches and infrastructure fault tests | Engineering fixture, not science |
| `quantity.convert.v1` | Tagged scale conversions with role/frame/convention metadata | No physical model or frame transformation; finite results currently unaccepted |
| `numerics.scalar_batch.v1` | Compensated sum, log-sum-exp, log1p, expm1, positive log-gamma | Bounded interface; finite results currently unaccepted |
| C++ numerical library | Integration and dense SPD/Gaussian building blocks with native tests | Specific sampled contracts, not unrestricted domain qualification |
| Local run store | Objects, resolved specifications, attempt receipts, hashes | Replay/cache/resume service not implemented |

The scientific qualification list in CLI discovery is currently empty. Native comparisons and regression tests are useful development evidence, but cannot be promoted into blanket acceptance of a scientific request. The integer fixture may be accepted as an engineering run while its discovery metadata still says `unqualified`.

## Direction of travel

The intended progression is physical definitions and numerical contracts; observation/data semantics and source adapters; probability and likelihood building blocks; physical models; inference; and reproducible end-to-end scientific analyses. Each consumer needs qualified prerequisites and its own error budget. Work in progress is not a promise that its CLI already exists.

GPU/distributed execution, broad cosmological analyses and a demonstrated performance advantage remain future work. Historical research results are not established merely by importing their ideas or documenting them.

## Known numerical limits

Adaptive quadrature can miss a narrow unsampled feature. Small SPD examples do not qualify arbitrary large or ill-conditioned matrices. Agreement through a shared system library is not independent numerical evidence. Domain rejection, including supported quantity roles and underflow policy, is part of the contract. See [testing](testing.md) and [fixture provenance](../cpp/tests/fixtures/README.md).
