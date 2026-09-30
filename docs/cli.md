# CLI contract

The executable is non-interactive. Inspect `irred describe --json` for its build and current request schema. The C ABI uses one current revision, with size, tag, count and version checks; obsolete layouts are not supported alongside it.

```sh
irred describe --json
irred version --json
irred run REQUEST.json STORE [--assurance numerical_contract|qualified]
irred stream STORE [LIMITS.json]
```

`run` evaluates one complete request and writes immutable input, resolved specification, output and execution records. It uses the same preparation/evaluation implementation as retained sessions. Rust parses and acquires sources; C++ owns equations and numerical calculations. No expression from a request is executable physics.

Current model/observation requests use schema version 2 and operation `background.evaluate`, `observations.prepare`, `statistics.gaussian`, `supernova.profile` or `bao.density`. Integer, quantity and scalar numerical fixtures also require request schema 2, retaining their stable physical operation identifiers. Use the generated schema for exact required fields, caps and enums; unknown fields are rejected.

Expansion descriptors contain only active parameters: LCDM omega_m, constant q, CPL omega_m/w0/wa, or five ordered q coefficients. FlatFLRW is explicit. Observer convention and physical scale are supplied for projections that need them. Request only the required output groups; partial failures preserve independent successful values. Source effects are separate from expansion, with none or explicit grey log1p magnitude shift in the SN consumer.

SN requests contain observations, an explicit selected index/coordinate descriptor or the named Pantheon selection, preparation policy, models, requested outputs and numerical policy. BAO requests contain a typed source, preparation policy, models with free H0rd, requested outputs and numerical policy. Gaussian requests contain current observations, explicit selection, ordered residual IDs and units, mode, optional response or proper prior, numerical guard, batch cap and preparation-byte cap. Raw ASCII adapters retain supplied source/axis/calibration/dependence declarations without independently verifying unlabeled matrix axes or published release identity.

## Retained sessions

`stream` accepts one JSON object per line. Actions are `prepare_observations`, `prepare_supernova`, `prepare_bao`, `background_evaluate`, `supernova_profile`, `bao_density` and `release`. Preparation returns a process-local monotonic handle and a scientific preparation digest. Evaluations bind the actual source, selection and preparation digest; temporary handle IDs and unrelated retained occupancy do not alter scientific identity.

Limits JSON explicitly supplies `sources`, `consumers`, `retained_bytes`, `commands`, `input_line_bytes` and `output_line_bytes`. Defaults are 16 source and consumer handles, 1 GiB retained payload, 4096 commands, 16 MiB input lines and 64 MiB output lines. Output limits must be at least 1024 bytes. Limits cover documented payload scopes, excluding caller input, allocator overhead and RSS. Native preparation peak allowance is reserved before allocations; shared sources stay charged until their last dependent owner is released.

Each admitted frame receives exactly one complete `{receipt,result}` reply. Receipts bind build, executable, raw input and reply objects; result retains the typed disposition and output checks. Oversized lines are discarded to their newline with bounded memory, then rejected. EOF closes the session. After the Nth admitted reply the process closes without reading an N+1 frame or emitting an extra request reply. A preparation reply that cannot fit rolls back the unpublished handle and charge while consuming its ID. Release is encoded before mutation. JSON is never truncated into a partial reply.

## Assurance and failures

Default `numerical_contract` returns exit 0 only when all required numerical checks pass. Structural source preparation uses `structural_source_contract`; checked integer addition uses `exact_integer_contract`. Acceptance names its scope and does not establish interpretation.

A completed scientific failure returns failed checks and exit 2. Malformed requests, resource/transport failures and internal exceptions remain execution failures with separate causes. Explicit `--assurance qualified` returns exit 6 only after numerical success when applicable qualification evidence is unavailable. Historical receipts are immutable and retain their original acceptance policy.

Preparation structural/resource rejection publishes no handle or invented empty source. Admitted numerical preparation failures preserve source and cause. Evaluations can return owned empty or partial diagnostics for configured exhaustion or incompatible preparation; failed groups never expose finite scientific payload. Suppressed exports still count toward retained native payload and workspace quotas.
