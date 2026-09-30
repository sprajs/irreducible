# Run records

A run store is the local directory supplied to `irred run`. It contains SHA-256-addressed `objects/` and `attempts/` records. It belongs outside version control.

The request bytes are retained as an object. A run records an incomplete attempt before calculation, then stores outputs and, when available, a resolved scientific specification and complete attempt. The incomplete record remains as part of the history. A terminated process may have an incomplete record and some objects without a complete record.

Records identify the source revision/status captured by the build, build identity, executable digest, runtime libraries, input/output digests, resolved specification and execution identity. Current deterministic operations record RNG as not applicable and declare precision/backend and resources. Source changes after building do not change an existing executable's embedded metadata; rebuild when source changes.

Object publication checks existing content rather than trusting a digest-shaped filename. Preserve immutable outputs. Never edit a receipt to make a failed run appear accepted. Objects on disk do not imply an implemented cache or resume facility.

| Status | Question |
| --- | --- |
| Execution | Did the calculation complete? |
| Numerical | Does it meet the applicable numerical contract and budget? |
| Inference | Are the estimator, support, diagnostics and uncertainty adequate? |
| Interpretation | What follows from the result and assumptions? |

The scientific CLI can complete while numerical qualification is not assessed. Inference may be not applicable. Keep these states explicit in agent summaries.

Full receipts, raw comparisons and acquired assets belong in local or purpose-built archival storage. For a published scientific result, preserve its exact run record and data identities in a suitable research archive and cite that separately. The source repository carries native tests and concise fixture provenance, not every development run. Removing material from Git does not mean deleting local historical originals.

Prepared observation runs retain table/matrix source bytes as independent objects. The resolved specification records their digests, semantic metadata, selection, ordering declarations and resource policy. A separate full-matrix JSON object stores all native retained values, axis IDs, uncertainty kind/unit, ordering provenance and source-asset digest; stdout references its digest and element count. Original columns/fields and source missing/nonfinite masks prevent unused values from being silently reinterpreted. No ordering, calibration or independence claim follows merely from retaining those records.

Operation outcomes declare their actual output IDs, methods/arithmetic and required numerical checks. A completed owned scientific failure remains execution `completed` with failed checks; request/transport/internal failures are distinct. The top-level numerical field summarizes required checks and never stands in for interpretation. Receipts record `assurance.requested`, `assurance.satisfied`, `accepted_scope` and `accepted`; default numerical-contract acceptance can be true while interpretation remains unqualified. Explicit qualified assurance returns exit 6 only after numerical success with missing applicable evidence. The assurance request belongs to execution identity and does not change the scientific specification. Existing immutable receipts retain their historical policy.

Mixed numerical failures preserve useful rows and return exit 2. Named regression coverage is bounded and is not automatically inherited by arbitrary source/model requests.

Versioned CPL batches retain every attempted `w0`/`wa` and inactive parameter value in ordered source rows and the resolved scientific specification, including failed rows. `supernova.profile_batch.v2` uses the same separation of completed execution, actual numerical checks and unqualified interpretation as v1. Changing an explicit parameter changes scientific identity; the v1 contract never supplies implicit CPL defaults.

Analytic supernova records retain all five q coefficients, the actual cached arithmetic, analytic segment policy, selected source identities and admitted segment counts. Source-row and full-matrix limits apply to preparation; evaluation byte accounting includes native scratch and all four internal per-model arrays even when exports are suppressed. Caller input and allocator overhead are outside this byte estimate. Numerical failure remains distinct from missing scientific acceptance: completed failures return exit 2; numerically successful unqualified calculations exit 0 by default or 6 under explicit qualified assurance.

Analytic BAO records retain required five-q/ruler model identities and separate actual preparation/evaluation analytic policies. Full source rows, declared covariance axes, raw assets and the complete matrix artifact remain immutable. Failed materialized rows retain attempted parameters without finite density fields; configured evaluation quota exhaustion can return an owned empty numerical work-limit result. Named native twenty-four-point coverage is recorded separately from request applicability, which remains unqualified.
