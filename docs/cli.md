# CLI contract

```sh
target/debug/irred describe --json
target/debug/irred version --json
target/debug/irred run tests/fixtures/exact-add.json runs/example
```

`describe` and `version` currently return the same metadata document. `run` takes a request-file path, a local store-directory path and optional `--assurance numerical_contract|qualified` (default `numerical_contract`). Requests are non-interactive and explicit. An agent can construct them in a script or an external fitting loop; use batches when available to avoid one process per evaluation. A request file configures compiled capabilities and cannot inject executable expressions.

## Requests

Requests use `schema_version: 1` and an explicit `operation`. Operation-specific parsing rejects unknown fields. JSON decimal numbers use the pinned binary64 roundtrip parser; resolved specifications retain the decoded source bits. Do not add fields from a proposed schema or coerce an unsupported scientific object into an existing operation.

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
| `6` | Numerical contract passed but explicit qualified assurance is unsatisfied | Preserve output; do not claim request-level scientific qualification |

Aborts/signals may leave an incomplete attempt. Finite output alone is not acceptance: receipts identify requested assurance, its satisfaction and accepted scope. Default numerical success can be accepted while interpretation remains unqualified; top-level numerical checks summarize actual required output contracts. Explicit qualified assurance requires applicable evidence and does not invent it.

`statistics.gaussian_batch.v1` is a retained Gaussian batch calculation. Its bounded native/CLI interface checks passed; completed calculations remain unqualified outside an applicable numerical qualification. See `tests/fixtures/gaussian-density.json` for a bounded synthetic request. The nested `observations` object uses the existing acquisition request. Supply exact `ordered_ids`, row-major `residuals`, a `residual_unit` matching the observations, `response_unit: "one"`, and explicit total-element and forward-sensitivity budgets.

`mode: "normalized_density"` returns normalized log densities and their quadratic, determinant and normalization terms. An optional `proper_prior` declares a named original latent prior, its mean/variance, dimensionless response, exact row IDs and independence assumption. These are retained in output. `mode: "profile_offset_score"` instead requires a response and returns a fitted coefficient and quadratic score; density and normalization are explicitly not applicable. A profile response and proper prior are mutually exclusive. Profile coefficients and prior means use the residual unit; prior variance uses its square. No optimizer or sampler is involved.

Outputs retain prepared distribution metadata separately from evaluation treatment, ordered IDs, matrix validation scope, kept/complement selection history and prior provenance. The full source uncertainty artifact remains unchanged. Failed batch rows omit finite numerical payloads. Source uncertainty and support assumptions remain the caller's explicit declarations.

## Background parameter/query batches

`background.parameter_query_batch.v1` accepts explicit `parameters`, shared `queries` and a resource/numerical `policy`. See `tests/fixtures/background-late.json`. Compiled model IDs are `flat_lcdm_late_v1` and `constant_q_flat_v1`; inactive parameter fields must be zero. Queries distinguish `geometric_same_redshift` from `released_zhd_zhel`. The structural schema is exposed as `abi_schema` in discovery; C++ owns physical domain checks.

Rows are ordered by parameter, then query. Each retains source parameters/query, equation and constant identities, numerical cause and callback count. Successful rows report distances in Mpc, lookback time in seconds and volume per steradian per redshift. Failed rows omit physical values. `maximum_total_evaluations` applies across every parameter/query pair. `maximum_native_output_bytes` bounds native owned slot storage; the CLI separately caps JSON output at 65,536 slots. Count, product and allocation arithmetic are checked before copying. The operation returns its calculation directly in `result`, without an observation-acquisition wrapper.

Bounded native/CLI interface checks passed. Completed calculations remain scientifically unqualified; numerical success exits 0 by default or 6 under explicit qualified assurance. Adaptive error estimates are empirical, and v1 provides no curvature, CPL, age or sound-horizon calculation.

## Retained supernova profile batches

`supernova.profile_batch.v1` reuses a nested observation request with profile `pantheon_plus_released_v1` and selection `pantheon_zhd_gt_001`. Provide a `models` array of `{model, omega_m, constant_q}` using the two v1 background IDs. There is no H0, CPL, latent-prior or optimizer field: this target profiles one additive magnitude offset and returns a relative score, not a normalized density or evidence.

The explicit `policy` declares `arithmetic` (`binary64_legacy_v1` or `longdouble_cpu_v1`), `include_residual_arrays`, positive `maximum_forward_sensitivity`, quadrature tolerances/depth/callback budgets and model/source/matrix/array/native-output caps. See the generated [run schema](../schema/run.schema.json) for required fields. One callback budget covers the entire native parameter batch. The array cap counts all four native per-model arrays, including adjusted diagnostic workspace, even when three exported arrays are omitted. Native output bytes do not include allocator overhead or prepared source/factor storage; those have separate source/matrix limits.

`result.observations` retains acquired rows and immutable full matrix artifacts; `result.calculation` retains selected row/source indices, zHD/zHEL conventions, original asset identities, arithmetic and equation IDs, profile scores and numerical causes. Optional shape/base/profiled residual arrays are empty when omitted. Failed rows omit finite payloads. Successful profile rows explicitly mark density and normalization not applicable.

Execution completion, numerical checks and interpretation are separate: successful native checks report `checks_passed`, while interpretation remains unqualified and numerical acceptance is scoped explicitly. Structural transport failures report no invented numerical cause. Named original-input seven-point validation coverage is inspectable development evidence and does not automatically qualify an arbitrary source/model request. Completed numerically checked batches return exit 0 by default; explicit qualified assurance returns exit 6 when applicable evidence is absent; numerical/domain/resource failures return exit 2 with their structured cause.

## Explicit CPL versioned batches

`background.parameter_query_batch.v2` and `supernova.profile_batch.v2` admit the additional `flat_cpl_late_v1` model. Every row requires explicit `model`, `omega_m`, `constant_q`, `w0` and `wa`; background rows also require `h0_km_s_mpc`. LCDM/CPL require `constant_q: 0`, constant-q requires `omega_m: 0`, and legacy models require `w0: -1, wa: 0`. Missing or unknown fields fail; inactive values are never silently ignored. These are the current compiled row/model restrictions; they are not an indefinite compatibility promise.

The v2 operations use the same prepared observations, conventions, resource policy and shared native physics as v1. A single callback budget covers mixed-model rows in source order. Outputs retain all attempted parameters even when a row fails, while failed rows omit finite scientific payloads. Named native comparisons and hostile transport tests cover explicit CPL; completed requests remain unqualified. The optional four-point original-input CLI regression exercises transport of the separately validated relative-profile calculation, without normalized density, evidence or inference claims.

## Retained BAO density batches

`bao.gaussian_batch.v1` accepts `source`, `prepare_policy`, `models` and `evaluation_policy`. Each model requires the five explicit v2 cosmological fields plus `h0_rd_km_s`. The adapter records fixed computational H0=70 km/s/Mpc; H0rd is a free empirical ruler parameter, without early-universe sound-horizon physics. See the generated run schema for all resource fields.

Source profiles are `desi_dr2_all_gccomb_13_v1` (exactly thirteen mean rows and a 13×13 covariance ASCII matrix) and `synthetic_inline_v1`. Literal source tags are `DM_over_rs`, `DH_over_rs` and `DV_over_rs`; unknown tags fail. Effective-redshift, free-ruler and computational-H0 conventions, calibration/dependence provenance and declared covariance ordering are required. Format-family declaration is distinct from exact asset identity; unlabeled matrix axes remain a supplied declaration. Original bytes and the full matrix are retained as immutable objects.

Separate explicit preparation and evaluation policies retain arithmetic, positive sensitivity guard, callback and allocation caps. A single callback budget covers mixed models in source order. Owned outputs retain IDs, attempted parameters, normalized log density, quadratic, log determinant and normalization; failed rows omit finite payloads. Optional predictions/residuals remain bounded even when export is disabled. Execution and numerical status remain separate from scientific acceptance: numerically checked calculations return exit 0 by default and remain unqualified for interpretation; explicit qualified assurance returns exit 6 when applicable evidence is absent. Named original-input eleven-point native comparisons do not establish arbitrary-request validity or cross-probe independence.

## Analytic fixed-bin background batches

`background.piecewise_query_batch.v1` requires `parameters` rows `{h0_km_s_mpc, q: [five values]}`, shared background `queries` and an analytic `policy`. It uses the compiled five-bin provider, with no new v1/v2 model tag. Policy fields are `maximum_parameters`, `maximum_queries`, `maximum_slots`, `maximum_native_output_bytes` and `maximum_total_segment_visits`. One analytic segment budget covers every row, including admitted work on failed queries. A query lacking enough remaining budget is rejected atomically; a later cheaper query may still run. There are no quadrature tolerances or callback fields. The native byte cap conservatively includes retained source/output and native scratch, excluding caller input and allocator overhead; JSON output is separately limited to 65,536 slots.

Results retain attempted H0/q/query values and source order. `geometry` exists only for successful slots. `derivatives` declares readable q convention, jerk availability and boolean availability flags; unavailable values are absent. Unequal internal jumps have right-limit q and no ordinary jerk, while equal adjacent bins retain ordinary jerk. Failed slots retain `not_assessed` derivative semantics. A finite batch container does not mean all rows passed. Structural errors have an error ID without invented numerical payload. Numerically checked requests return exit 0 by default and remain unqualified for interpretation; explicit qualified assurance returns exit 6 when applicable evidence is absent; failed calculations return structured causes.

`supernova.piecewise_profile_batch.v1` accepts the same nested retained `observations` request as the other supernova operations, `models: [{"q": [q0, q1, q2, q3, q4]}]`, and a required analytic policy. It requires explicit arithmetic and positive forward-sensitivity budget, with source-row/full-matrix preparation caps and separate selected-query, model, four-array, aggregate native-byte and total analytic-segment evaluation caps. See the generated schema for exact fields. No quadrature tolerances or callback settings are supplied. Failed materialized rows retain attempted q values without a finite score; configured evaluation quota exhaustion yields an owned empty work-limit result. Preparation caps retain structural rejection semantics.

`bao.piecewise_gaussian_batch.v1` reuses the released thirteen-row or truthful synthetic BAO source adapter. Each required model is `{"q": [q0, q1, q2, q3, q4], "h0_rd_km_s": value}`. Preparation and evaluation policies are separate flat analytic packets with no quadrature fields. Preparation bounds source rows, full matrix elements, strings and native preparation bytes. Evaluation bounds models, selected queries, two retained arrays per model, total analytic segments and aggregate native output/scratch bytes even when exports are suppressed. Arithmetic and positive sensitivity budgets apply to both phases; valid precision mismatch delegates an owned empty native semantic failure. The computational H0=70 convention is recorded and does not add an independent H0 inference parameter. Distinct typed native handles preserve the existing BAO owner/layout and policies.

Numerical-contract assurance exits 0 only when every required scientific output reports `checks_passed` (genuinely nonnumerical outputs may be `not_applicable`). Completed scientific failures exit 2 and keep execution `completed`; malformed requests, transport failures and internal faults have a separate failed execution. `--assurance qualified` exits 6 only when numerical checks pass but applicable request-level qualification is missing. `accepted_scope` and requested/satisfied assurance make this scope explicit; a numerical success never claims inference or interpretation. Assurance is part of execution identity, not a change to the physical equations.
