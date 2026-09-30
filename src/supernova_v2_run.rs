//! Resolve immutable observations and explicit native policy; no equations here.
use crate::{
    bridge::{SupernovaModelV2Request, SupernovaPolicyV2Request, supernova_v2_evaluate},
    observation_run,
};
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};
use std::path::Path;
#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
struct Request {
    schema_version: u32,
    operation: String,
    observations: Value,
    models: Vec<SupernovaModelV2Request>,
    policy: SupernovaPolicyV2Request,
}
pub(crate) fn execute(input: &[u8], store: &Path) -> Result<(Value, Value), String> {
    let request: Request = serde_json::from_slice(input).map_err(|e| e.to_string())?;
    if request.schema_version != 1 || request.operation != "supernova.profile_batch.v2" {
        return Err("UNSUPPORTED_SPECIFICATION".into());
    }
    let nested = serde_json::to_vec(&request.observations).map_err(|e| e.to_string())?;
    let (observations, output) =
        observation_run::execute_calculation(&nested, store, |source, selection, policy| {
            if std::env::var_os("COSMOLOGY_TEST_PANIC").is_some() {
                panic!("injected Rust panic");
            }
            supernova_v2_evaluate(source, selection, policy, &request.models, &request.policy)
                .map(Some)
        })?;
    let spec = json!({"schema_version":1,"operation":request.operation,"equation_id":"W01/released-profile-score/v1",
        "observations":observations,"models":request.models,"policy":request.policy,
        "target":"single-offset relative profile score; normalized density and evidence not applicable",
        "units":{"shape":"mag","offset":"mag","quadratic":"one","relative_profile_score":"one"},
        "resource_interpretation":{"maximum_native_output_bytes":"owned native model result slots and arrays; allocator overhead and source/factor storage separately bounded",
        "maximum_array_elements":"all four native per-model arrays, including adjusted diagnostic workspace, even if JSON residual arrays omitted",
        "maximum_total_evaluations":"one native callback budget across all models and selected observations"},
        "requested_outputs":[{"id":"profile_scores","required":true,"numerical_gate":"required","inference_gate":"not_applicable"}]});
    Ok((spec, output))
}
