//! Current requested-output expansion operation, shared by one-shot/stream.
use crate::{
    bridge::expansion::{self, Policy},
    model_spec::{BackgroundRequest, ExpansionSpec, Geometry},
    outcome::{NumericalCheck, Outcome, OutputCheck},
};
use serde::Deserialize;
use serde_json::json;

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Request {
    schema_version: u32,
    operation: String,
    geometry: Geometry,
    models: Vec<ExpansionSpec>,
    queries: Vec<BackgroundRequest>,
    numerical_policy: Policy,
}
pub(crate) fn execute(input: &[u8]) -> Result<Outcome, String> {
    let request: Request = serde_json::from_slice(input).map_err(|e| e.to_string())?;
    if request.schema_version != 2 || request.operation != "background.evaluate" {
        return Err("UNSUPPORTED_SPECIFICATION".into());
    }
    let calculation = expansion::evaluate(
        &request.models,
        request.geometry,
        &request.queries,
        &request.numerical_policy,
    )?;
    let requested:Vec<_>=calculation.checks.iter().map(|(id,_)|json!({"id":id,"required":true,"numerical_gate":"required","inference_gate":"not_applicable"})).collect();
    let spec = json!({"schema_version":2,"operation":"background.evaluate","geometry":request.geometry,
        "models":request.models,"queries":request.queries,"numerical_policy":request.numerical_policy,
        "requested_outputs":requested,"constants_id":"SI-IAU-definitions-v1",
        "node_reuse_scope":"exact expansion-z bits within one model batch; no cross-model cache"});
    let resources = json!({"maximum_models":request.numerical_policy.maximum_models,
        "maximum_queries":request.numerical_policy.maximum_queries,"maximum_slots":request.numerical_policy.maximum_slots,
        "maximum_callbacks":request.numerical_policy.maximum_callbacks,"maximum_segment_visits":request.numerical_policy.maximum_segment_visits,
        "maximum_native_bytes":request.numerical_policy.maximum_native_bytes,
        "byte_scope":"combined wrapper/native owned payload and sequential workspace; excludes borrowed caller input, allocator metadata and RSS"});
    let mut outcome = Outcome::scientific(
        spec,
        calculation.output,
        "expansion",
        json!("compiled_flat_flrw_requested_projections"),
        json!("binary64_storage_longdouble_intermediate"),
        resources,
        "named native reference cases; arbitrary input applicability and interpretation remain unqualified",
    );
    outcome.outputs=calculation.checks.into_iter().map(|(id,passed)|OutputCheck{
        check_kind:"numerical_contract",id,required:true,
        numerical:if passed{NumericalCheck::ChecksPassed}else{NumericalCheck::Failed},inference:"not_applicable",
        interpretation:"unqualified",evidence:vec![],validation_coverage:"compiled requested projection contract; named fixture coverage separate from request qualification"
    }).collect();
    Ok(outcome)
}
