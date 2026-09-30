//! Resolve immutable observations and explicit native policy; no equations here.
use crate::{
    bridge::{
        SupernovaPiecewiseModelRequest, SupernovaPiecewisePolicyRequest,
        supernova_piecewise_evaluate,
    },
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
    models: Vec<SupernovaPiecewiseModelRequest>,
    policy: SupernovaPiecewisePolicyRequest,
}
fn calculation(input: &[u8], store: &Path) -> Result<(Value, Value), String> {
    let request: Request = serde_json::from_slice(input).map_err(|e| e.to_string())?;
    if request.schema_version != 1 || request.operation != "supernova.piecewise_profile_batch.v1" {
        return Err("UNSUPPORTED_SPECIFICATION".into());
    }
    let nested = serde_json::to_vec(&request.observations).map_err(|e| e.to_string())?;
    let (observations, output) =
        observation_run::execute_calculation(&nested, store, |source, selection, policy| {
            if std::env::var_os("COSMOLOGY_TEST_PANIC").is_some() {
                panic!("injected Rust panic");
            }
            supernova_piecewise_evaluate(
                source,
                selection,
                policy,
                &request.models,
                &request.policy,
            )
            .map(Some)
        })?;
    let spec = json!({"schema_version":1,"operation":request.operation,"equation_id":"W01/released-profile-score/v1",
        "observations":observations,"models":request.models,"policy":request.policy,
        "target":"single-offset relative profile score; normalized density and evidence not applicable",
        "units":{"shape":"mag","offset":"mag","quadratic":"one","relative_profile_score":"one"},
        "resource_interpretation":{"maximum_native_output_bytes":"conservative owned native slots/arrays/source identities and query scratch; excludes caller input/allocator overhead; preparation separately source-row/full-matrix bound",
        "maximum_array_elements":"all four native per-model arrays, including adjusted diagnostic workspace, even if JSON residual arrays omitted",
        "maximum_total_segment_visits":"one analytic admitted segment budget across all models and selected observations; includes failed admitted work"},
        "requested_outputs":[{"id":"profile_scores","required":true,"numerical_gate":"required","inference_gate":"not_applicable"}]});
    Ok((spec, output))
}

pub(crate) fn execute(input: &[u8], store: &Path) -> Result<crate::outcome::Outcome, String> {
    let (spec, output) = calculation(input, store)?;
    let arithmetic = spec["policy"]["arithmetic"].clone();
    let resources = spec["policy"].clone();
    Ok(crate::outcome::Outcome::scientific(
        spec,
        output,
        "profile_scores",
        json!("retained_single_offset_profile_analytic_segments"),
        arithmetic,
        resources,
        "named original-input six-point piecewise regression; request applicability not established",
    ))
}
