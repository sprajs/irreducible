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

For `supernova.profile_batch.v1`, a completed process and its native checks are reported separately from scientific qualification. The resolved specification includes actual source content identities, selected ordering, model batch and explicit arithmetic/resource policy. The receipt retains that policy and precision ID. Its per-output `numerical` field reports the calculation's checks (`checks_passed`, `failed`, or `not_evaluated` for a structural transport failure); the top-level numerical qualification remains `not_assessed`, and interpretation/acceptance remain unqualified. Mixed numerical failures can leave a completed execution record and useful successful rows while the command returns exit 2. Exit 6 means finite execution with missing applicable acceptance, not a numerical failure. Original-input regression coverage is named and bounded; it is not silently inherited by a changed source or model request.

Versioned CPL batches retain every attempted `w0`/`wa` and inactive parameter value in ordered source rows and the resolved scientific specification, including failed rows. `supernova.profile_batch.v2` uses the same separation of completed execution, actual numerical checks and unqualified interpretation as v1. Changing an explicit parameter changes scientific identity; the v1 contract never supplies implicit CPL defaults.

Analytic supernova records retain all five q coefficients, the actual cached arithmetic, analytic segment policy, selected source identities and admitted segment counts. Source-row and full-matrix limits apply to preparation; evaluation byte accounting includes native scratch and all four internal per-model arrays even when exports are suppressed. Caller input and allocator overhead are outside this byte estimate. Numerical failure remains distinct from missing scientific acceptance: completed failures return exit 2, finite unqualified calculations exit 6.
