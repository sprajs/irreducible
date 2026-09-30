# CLI contract

```sh
target/debug/irred describe --json
target/debug/irred version --json
target/debug/irred run tests/fixtures/exact-add.json runs/example
```

`describe` and `version` currently return the same metadata document. `run` takes exactly a request-file path and a local store-directory path. There is no implicit interactive session, expression evaluator, inference command or automatic replay command.

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

## Outputs and exit codes

For an evaluated request, stdout contains `result` and `receipt`. Errors also produce a structured failure on stderr. Some early failures occur before a complete receipt can be written. See [run records](run-records.md).

| Code | Meaning | Agent action |
| --- | --- | --- |
| `0` | Command succeeded; the exact fixture can be accepted | Inspect operation and receipt; discovery is not analysis |
| `2` | Usage, parsing, unsupported request or evaluation failure | Read stderr and any result/receipt; fix the cause |
| `6` | Scientific calculation completed but numerical qualification is missing | Preserve the output as unaccepted; do not report qualified science |

Aborts/signals are separate and may leave an incomplete attempt. Do not treat every nonzero exit as “no output”, or every finite value as acceptance. Current scientific operations return `accepted=false` and `numerical: "not_assessed"` even after finite execution. Future qualification behavior must be documented and tested before an agent relies on it.
