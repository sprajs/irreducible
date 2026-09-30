# CLI contract

```sh
target/debug/irred describe --json
target/debug/irred version --json
target/debug/irred run tests/fixtures/exact-add.json runs/example
```

`describe` and `version` currently return the same metadata document. `run` takes exactly a request-file path and a local store-directory path. Requests are non-interactive and explicit. An agent can construct them in a script or an external fitting loop; use batches when available to avoid one process per evaluation. A request file configures compiled capabilities and cannot inject executable expressions.

## Requests

Requests use `schema_version: 1` and an explicit `operation`. Operation-specific parsing rejects unknown fields. Do not add fields from a proposed schema or coerce an unsupported scientific object into an existing operation.

An exact infrastructure request:

```json
{"schema_version":1,"operation":"fixture.checked_i64_add.v1","a":[1,2],"b":[3,4]}
```

Arrays must have compatible lengths and sums must fit signed 64-bit integers. The optional `fault` field and `COSMOLOGY_TEST_*` environment controls are for tests; leave them unset in normal runs.

For tagged conversions, use the complete [quantity request](../tests/fixtures/quantity-length.json). Metadata includes unit, physical role, frame, convention and constant-set identity. Enum definitions come from [schema/abi.json](../schema/abi.json). Unit conversion does not imply a physical-role, epoch or reference-frame transformation.

A scalar/reduction request:

```json
{"schema_version":1,"operation":"numerics.scalar_batch.v1","method":"compensated_sum","values":[10000000000000000,1,-10000000000000000]}
```

Methods are `compensated_sum`, `log_sum_exp`, `log1p`, `expm1` and `log_gamma_positive`. The first two reduce a batch; the others evaluate ordered elements. Domain/size restrictions are enforced by the implementation. Diagnostics are labelled as empirical or arithmetic diagnostics, not certified bounds.

## Prepared ASCII observations

Use [the synthetic observation request](../tests/fixtures/observations-gaussian.json) for a bounded local example. `metadata` declares profile, role, unit, calibration, uncertainty type/unit and component. `selection` is `all` or the compiled `pantheon_zhd_gt_001` selector. Explicit resources limit asset bytes, rows, matrix elements and metadata strings. Unknown fields/tags and incompatible scientific semantics fail closed.

`table` and optional `uncertainty` are local paths resolved from the process working directory. Each asset is read once into bounded bytes, hashed and retained. The scientific specification uses those content hashes, so changing bytes at the same path changes scientific identity. The unlabeled matrix's axis linkage is a supplied ordering convention, not something independently verified from its bytes; declare `ordering_provenance` explicitly.

Output preserves source rows and full uncertainty, and separately reports the selected mask and ordered source indices. The full matrix is an immutable object reference, not a large stdout array. ASCII readers do not decode quality flags: zeros are explicitly labelled structural placeholders, not all-clear flags. Missingness and nonfinite source values have separate masks; required selected missing/nonfinite values fail. There is no covariance projection, numerical repair, event deduplication or inference in this operation.

## Outputs and exit codes

For an evaluated request, stdout contains `result` and `receipt`. Errors also produce a structured failure on stderr. Some early failures occur before a complete receipt can be written. See [run records](run-records.md).

| Code | Meaning | Agent action |
| --- | --- | --- |
| `0` | Command succeeded; the exact fixture can be accepted | Inspect operation and receipt; discovery is not analysis |
| `2` | Usage, parsing, unsupported request or evaluation failure | Read stderr and any result/receipt; fix the cause |
| `6` | Scientific calculation completed but numerical qualification is missing | Preserve the output as unaccepted; do not report qualified science |

Aborts/signals are separate and may leave an incomplete attempt. Do not treat every nonzero exit as “no output”, or every finite value as acceptance. Current scientific operations return `accepted=false` and `numerical: "not_assessed"` even after finite execution. Future qualification behavior must be documented and tested before an agent relies on it.

`statistics.gaussian_batch.v1` is a retained Gaussian batch calculation. Its bounded native/CLI interface checks passed; completed calculations remain unqualified outside an applicable numerical qualification. See `tests/fixtures/gaussian-density.json` for a bounded synthetic request. The nested `observations` object uses the existing acquisition request. Supply exact `ordered_ids`, row-major `residuals`, a `residual_unit` matching the observations, `response_unit: "one"`, and explicit total-element and forward-sensitivity budgets.

`mode: "normalized_density"` returns normalized log densities and their quadratic, determinant and normalization terms. An optional `proper_prior` declares a named original latent prior, its mean/variance, dimensionless response, exact row IDs and independence assumption. These are retained in output. `mode: "profile_offset_score"` instead requires a response and returns a fitted coefficient and quadratic score; density and normalization are explicitly not applicable. A profile response and proper prior are mutually exclusive. Profile coefficients and prior means use the residual unit; prior variance uses its square. No optimizer or sampler is involved.

Outputs retain prepared distribution metadata separately from evaluation treatment, ordered IDs, matrix validation scope, kept/complement selection history and prior provenance. The full source uncertainty artifact remains unchanged. Failed batch rows omit finite numerical payloads. Source uncertainty and support assumptions remain the caller's explicit declarations.

## Background parameter/query batches

`background.parameter_query_batch.v1` accepts explicit `parameters`, shared `queries` and a resource/numerical `policy`. See `tests/fixtures/background-late.json`. Compiled model IDs are `flat_lcdm_late_v1` and `constant_q_flat_v1`; inactive parameter fields must be zero. Queries distinguish `geometric_same_redshift` from `released_zhd_zhel`. The structural schema is exposed as `abi_schema` in discovery; C++ owns physical domain checks.

Rows are ordered by parameter, then query. Each retains source parameters/query, equation and constant identities, numerical cause and callback count. Successful rows report distances in Mpc, lookback time in seconds and volume per steradian per redshift. Failed rows omit physical values. `maximum_total_evaluations` applies across every parameter/query pair. `maximum_native_output_bytes` bounds native owned slot storage; the CLI separately caps JSON output at 65,536 slots. Count, product and allocation arithmetic are checked before copying. The operation returns its calculation directly in `result`, without an observation-acquisition wrapper.

Bounded native/CLI interface checks passed. Completed calculations remain unqualified with exit 6. Adaptive error estimates are empirical, and no curvature, CPL, age or sound-horizon calculation is provided.
