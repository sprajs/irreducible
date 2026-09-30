# Run recipes (proposed)

A run recipe would preconfigure a bounded sequence of data acquisition, preparation and compiled calculations. A person or agent could review it, rerun the same paper comparison, or create an explicitly changed variant without reconstructing a conversation. This is proposed work: the current [CLI](cli.md) accepts individual requests and retained-session commands, not recipe files. No recipe schema or runner is implemented.

Use **strict JSON as the single canonical recipe format**. It matches the existing JSON requests, records and [request schema](../schema/run.schema.json), and fits schema-constrained agent tools. This is a project integration decision, not a claim that JSON is universally easiest, fastest or most accurate for language models. Keep explanatory `description` fields and surrounding Markdown for human commentary.

## Format decision

| Format | Useful properties | Contract needed here |
| --- | --- | --- |
| JSON | Small value model; ordered arrays; direct fit to current requests and JSON Schema | Reject duplicate decoded member names before lossy parsing; reject unknown fields; require finite, representable typed numbers; use arrays for meaningful order; no comments |
| TOML | Comments, explicit scalar types, readable configuration tables and arrays of tables | Define conversion of dates and nested step tables into operation types; prohibit special nonfinite floats in scientific parameters; duplicate keys are invalid |
| YAML | Comments, readable mappings/sequences, tags and aliases for richer graphs | Pin version and scalar-resolution schema, string keys and allowed types; bound or forbid aliases/cycles and custom tags; enforce unique keys; preserve source bytes before conversion |

[RFC 8259](https://www.rfc-editor.org/rfc/rfc8259.html) recommends unique JSON names but describes inconsistent receiver behavior for duplicates. JSON itself excludes NaN/Infinity; it does not guarantee that a legal exponent fits binary64. Interoperable binary64 integers are restricted to the exact range around ±2^53; larger counters require an explicit typed encoding. Preserve signed zero and validated scientific float bits in resolved records. Object member order is not a scientific ordering contract.

[TOML 1.0](https://toml.io/en/v1.0.0) prohibits repeated keys and permits comments, date/time values and `inf`/`nan`. [YAML 1.2.2](https://yaml.org/spec/1.2.2/) requires unique mapping keys and supports aliases, tags and different scalar-resolution schemas. These are legitimate features, not proof that YAML is unsafe; they enlarge the profile a scientific runner would have to specify. Supporting three runtime formats would add conversion and identity tests without a demonstrated first-use benefit.

[OpenAI Structured Outputs](https://developers.openai.com/api/docs/guides/structured-outputs) and [function calling](https://developers.openai.com/api/docs/guides/function-calling) support JSON Schema-constrained output/tool arguments, with supported-subset restrictions. [Google's structured-output interface](https://ai.google.dev/gemini-api/docs/structured-output) also supports a JSON Schema subset. Neither establishes scientific correctness. Editing a JSON file does not activate API constrained generation: those guarantees apply only when an integration actually uses the supported API mode. Recipe validation remains provider-independent and introduces no OpenAI runtime dependency. Generate any tool-compatible schema projection from the authoritative recipe contract and test it for agreement; do not maintain a second competing definition.

Format studies exist: [StructEval](https://arxiv.org/abs/2505.20139) evaluates structured generation and conversion across formats, while [structured context engineering](https://arxiv.org/abs/2602.05447) uses SQL generation to study input context. Reading context, converting formats and authoring valid scientific recipes are different tasks. Their results are bounded by the models, tasks and evaluation settings studied, not a universal ranking of formats. This proposal has not benchmarked agent recipe-writing accuracy, token use or speed.

## What a recipe would contain

The following is a compact **illustration of proposed structure**, not an accepted request, a complete schema example or a command to execute. Placeholder hashes, URLs and declarations must be replaced and verified. The referenced request file would contain the complete current typed BAO request; a future resolver would bind its declared local mean-file input to the acquired artifact. Its covariance and other required inputs must also be pinned in a real recipe.

```json
{
  "kind": "run_recipe",
  "description": "Conditional paper comparison; not an exact reproduction claim",
  "build": {"revision_hint": "<commit>", "build_id": "<exact build identity>",
            "capabilities": ["bao.density"]},
  "limits": {"steps": 2, "download_bytes": 1048576, "retained_bytes": 67108864},
  "sources": [{
    "id": "means", "url": "https://example.invalid/release/means.txt",
    "release": "<release>", "sha256": "<64 hexadecimal digits>",
    "license": "<verified terms or unresolved restriction>",
    "role": "released_fitted_summary", "credential_reference": null
  }],
  "steps": [
    {"id": "fetch", "kind": "acquire", "source": "means", "output": "mean_file"},
    {"id": "compare", "kind": "compiled_request", "operation": "bao.density",
     "request_file": {"path": "bao-request.json", "sha256": "<64 hexadecimal digits>"},
     "inputs": {"mean_file": {"artifact": "fetch.mean_file"}},
     "assurance": "numerical_contract"}
  ]
}
```

References are typed names, not interpolation expressions. The complete operation request owns active model parameters, geometry/source effects, numerical policies, requested outputs and source metadata. The `request_file` notation packages that request; it does not introduce a template, scientific-parameter override or merge language. Artifact binding resolves only declared typed file inputs. The revision hint is descriptive; the exact build identity precondition must cover source state, compiler flags and backend rather than assuming a commit alone identifies the executable. A recipe names compiled capabilities already present in the requested build. It cannot define new equations, execute shell/Python, evaluate expressions, install packages or invent missing physics. An agent develops a missing model in source, verifies it and rebuilds before creating the executable recipe.

Immutable file artifacts and process-local prepared owners are different reference types. A saved `observations.prepare` output is not a serialized native observation object. A runner may keep a source or factor alive using the existing retained-session semantics, with bounded lifetime and final release. Handles/pointers cannot be stored as reusable cross-process scientific inputs. A new process must reconstruct an owner from identified assets and preparation policy; restartable factor serialization is outside the first slice.

## Admission and execution rules

Use one authoritative future schema with closed typed step variants and typed outputs. Reuse current operation definitions rather than copying scientific fields into a second model language. JSON Schema provides annotations and validation structure, but graph compatibility requires additional checks. Its [`default` annotation](https://json-schema.org/understanding-json-schema/reference/annotations) does not insert values; any application defaults must be explicit in the resolved record.

Before acquisition, validate the entire structural graph: unique IDs, existing references, no cycles, output/input types, explicit ordered axes, legal operation/model combinations, requested build/capabilities, step/resource limits and credential names. Preserve array order; never sort covariance axes or model vectors to make them match. After bounded acquisition, validate actual hashes, decoded shapes, units/frames, selection, masks, calibration, artifact roles and dependence before scientific execution. Remote data cannot be fully admitted from a URL alone.

Acquisition needs explicit allowed schemes/hosts, redirect count and cross-origin credential policy, timeouts, concurrency, download and decoded-size bounds, and license/access disposition. Credentials are references to an external credential provider, never secret values in recipes or receipts. Record actual origin and redirects without credential-bearing URLs. An unresolved license may block redistribution even if local access is allowed. An unknown hash permits discovery only under an explicitly separate admission policy; it cannot silently become a pinned repeatable input.

The initial executor should be a serial local bounded runner over existing typed operations and retained owners. No scheduler, server or distributed checkpoint service is required. Define cancellation at safe boundaries, release all owners, preserve an incomplete attempt, and publish no false completed result. Failed required outputs fail their step; downstream dependencies are blocked. Independent successful outputs remain inspectable as partial results, never silently accepted as the full recipe. No retry may change assets, tolerances, rows or models without a new identified variant.

Reuse within a run is explicit ownership reuse. Omit automatic cross-run cache/resume initially. Any later reuse must verify asset hashes, model/source-effect code, selection/covariance, numerical contract and applicable build/backend identity; a matching filename or old handle is insufficient. Recipe steps do not turn current immutable run stores into an implemented cache.

## Records and scientific meaning

Retain exact recipe/request/asset bytes, their hashes, the fully resolved typed graph and defaults, actual executable/build/library identities, per-step input/output digests and dependency identities. Separate stable scientific identity from execution details such as handles, occupancy, credentials references, assurance request and effective quotas. Record requested versus actual resource/backend configuration; an unexpected build must reject a build precondition rather than quietly substitute it.

Keep execution, numerical, inference and interpretation states distinct, following [run records](run-records.md). A recipe's success is not automatic qualification. For stochastic future steps, record generator, seed/counter and stream assignment; deterministic current operations use not-applicable RNG. A seed alone cannot promise identical GPU/distributed reductions, topology or library behavior. Define bit reproducibility or statistical equivalence per backend and consumer budget before claiming it.

For a paper workflow, identify paper/version, equations, data release, preprocessing, nuisance measure and the exact tested claim. Preserve faithful transcription separately from repaired or extended models. Declare exact, approximate, conditional or blocked reproduction with reasons. Model-dependent calibration, reconstruction and selection must be recomputed or justified under the new theory. A joint comparison requires repeated-object/shared-calibrator identities and overlap/cross-covariance; unknown dependence is not independence. A list of separate likelihood steps is not a qualified joint fit.

## Small first implementation and its gate

First implement structural validation and resolution for a small acyclic recipe, pinned local artifacts and existing compiled requests. Add a bounded acquisition step for one narrow supported adapter, then in-process preparation reuse. Keep all scientific calculations in their current owners. Multi-paper search, automatic fitting, remote scheduling and distributed restart are later concrete consumers, not prerequisites.

Freeze acceptance cases and consumer error allocations before candidate outputs: duplicate/unknown fields, cycles, wrong reference types/order, unsupported capabilities, hash/build mismatch, missing assets, acquisition quota/redirect failures, cancellation, partial failures and owner cleanup. Assert one-shot versus recipe numeric/source identity at the same policy, and scientific identity stability across ephemeral handles. Challenge stale reuse explicitly. Existing mathematical oracles and tolerances remain unchanged; orchestration agreement is transport evidence, not an independent scientific algorithm. Test each new stochastic/backend claim with its own reproducibility and downstream allocation.

This document is a design recommendation. No recipe runner, recipe validation gate, acquisition performance result or paper reproduction has been executed by writing it.
