# Independent S00 verification contract

Owner: verifier. Contract specified before implementation inspection, 2026-09-30.
This file implements no capability and records no passing scientific gate.

## T00 bridge gate

The arithmetic fixture uses dyadic integers/fractions with exactly representable
products/sums; expected values are derived using Python Fraction, never the core.
Direct and CLI outputs must agree with that oracle and preserve batch ordering.
Their scientific operation/input identity must agree; binary and attempt identities
must remain distinct. Empty input must have an explicit successful empty output.

Construct descriptors independently from the public C header. Exercise ABI mismatch,
undersized struct, nonzero reserved fields if forbidden, wrong type/address space,
misalignment, nonzero count with null pointer, short byte length, overflowed count,
shape/stride mismatch and null output destination. Only pointers to owned live test
allocations are used; arbitrary invalid addresses are not a memory-safety test.
Failure must be stable integer status and null output handle. A valid pointer with
zero elements follows the documented empty contract. Verify layout and enum values
against Rust bindings mechanically rather than assuming repr(C) suffices.

Inject std::exception, unknown exception and bad_alloc separately. Inject diagnostic
allocation failure and cleanup failure independently; observe stable original error,
null optional diagnostics and no unwind. Scope/double-close safe-wrapper tests must
challenge ownership. Raw double free is not promised safe unless explicitly stated.
Use sanitizer execution when supported; absence of findings is scoped evidence.

Check Rust panic policy as configured, both controlled Rust errors and forced process
abort. An abort is tested in a child process and requires an already durable incomplete
attempt; no catch_unwind promise is inferred. No cross-language callbacks are needed.

## T01 record gate

Reconstruct an output from manifest bytes, source/build identity, exact input digests,
resolved request, equation ID and output digest. Rewrite input bytes at the same path
while preserving size to challenge path-only and size-only identity. Repeat a request;
retain separate attempt IDs or explicitly documented immutable idempotent semantics.
Preexisting result files must never be overwritten. Change Rust and C++ source
individually in isolated build copies and require build identity change. Verify source
identity includes dirty content and required untracked source, not just Git HEAD.

Kill after attempt creation and before finalization. The incomplete attempt must be
recognizable, excluded from completed results, and preserve old results. Corrupt a
stored digest object, omit a dependency, and present stale model/schema/build/input
identity: refuse reuse rather than silently recompute or accept stale results.
Atomic rename alone is not evidence of durable fsync; record the implemented policy.

Execution completed/failed/cancelled/incomplete and numerical/inference/interpretation
statuses remain distinct. S00 arithmetic cannot qualify cosmological science. Tagged
outside_support and numerical_failure must be distinguishable if claimed; otherwise
the capability descriptor must mark that future interface slice unavailable.

## Evidence and promotion

Each actual check records command, executable digest, input identity, expected and
observed results, exit code and status. Unimplemented injection surfaces are missing
gates, never assumed pass. Direct/CLI agreement has shared core ancestry and qualifies
only the bridge. Failed expectations freeze both observations under V02. S00 advances
only for its explicitly tested scope, with deferred broader E05 interfaces listed.

Build identity mutation matrix: compiled Rust/C++ sources, C header/schema generator,
model source, build scripts/configuration and dependency locks must affect scientific
build identity. Mutable receipts/output directories must not affect it. Add/remove an
output receipt and rebuild unchanged source: identity must stay stable, preventing a
receipt/rebuild feedback loop. Mutation experiments use isolated copies when possible;
otherwise preserve original bytes and restore in finally. Keep each discrepancy and
retest in a separate receipt; never replace the first failed observation.
