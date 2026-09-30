//! Explicit request resolution; no physical equations or unit conversions here.
use crate::bridge::{
    BackgroundParametersRequest, BackgroundPolicyRequest, BackgroundQueryRequest,
    background_evaluate,
};
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};
use std::path::Path;
#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
struct Request {
    schema_version: u32,
    operation: String,
    parameters: Vec<BackgroundParametersRequest>,
    queries: Vec<BackgroundQueryRequest>,
    policy: BackgroundPolicyRequest,
}
fn calculation(input: &[u8], _store: &Path) -> Result<(Value, Value), String> {
    let request: Request = serde_json::from_slice(input).map_err(|e| e.to_string())?;
    if request.schema_version != 1 || request.operation != "background.parameter_query_batch.v1" {
        return Err("UNSUPPORTED_SPECIFICATION".into());
    }
    let spec = json!({"schema_version":1,"operation":request.operation,"equation_id":"P01/compiled-late-background-batch/v1","parameters":request.parameters,"queries":request.queries,"policy":request.policy,"units":{"h0":"km/s/Mpc","redshift":"one","distance":"Mpc","lookback":"s","volume":"Mpc^3/sr/redshift"},"ordering":"parameter major, query minor","resource_interpretation":{"maximum_native_output_bytes":"native owned slot storage; JSON is separately bounded to 65536 slots","maximum_total_evaluations":"one integration-callback budget across all parameters and queries"},"requested_outputs":[{"id":"background_slots","required":true,"numerical_gate":"required","inference_gate":"not_applicable"}]});
    if std::env::var_os("COSMOLOGY_TEST_PANIC").is_some() {
        panic!("injected Rust panic");
    }
    let result = background_evaluate(&request.parameters, &request.queries, &request.policy)?;
    Ok((spec, result))
}

pub(crate) fn execute(input: &[u8], store: &Path) -> Result<crate::outcome::Outcome, String> {
    let (spec, output) = calculation(input, store)?;
    let arithmetic = json!("binary64_storage_longdouble_intermediate");
    let resources = spec["policy"].clone();
    Ok(crate::outcome::Outcome::scientific(
        spec,
        output,
        "background_slots",
        json!("adaptive_simpson"),
        arithmetic,
        resources,
        "bounded native model/reference and interface tests; request applicability not established",
    ))
}
