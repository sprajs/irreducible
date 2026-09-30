//! Resolve a fixed compiled analytic provider request without physics in Rust.
use crate::bridge::piecewise;
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};
use std::path::Path;
#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
struct Request {
    schema_version: u32,
    operation: String,
    parameters: Vec<piecewise::Parameters>,
    queries: Vec<piecewise::Query>,
    policy: piecewise::Policy,
}
fn calculation(input: &[u8], _store: &Path) -> Result<(Value, Value), String> {
    let request: Request = serde_json::from_slice(input).map_err(|e| e.to_string())?;
    if request.schema_version != 1 || request.operation != "background.piecewise_query_batch.v1" {
        return Err("UNSUPPORTED_SPECIFICATION".into());
    }
    let spec = json!({"schema_version":1,"operation":request.operation,"equation_id":"P01/fixed-five-bin-q-flat-kinematic/v1","parameters":request.parameters,"queries":request.queries,"policy":request.policy,"units":{"h0":"km/s/Mpc","q":"one","redshift":"one","distance":"Mpc","lookback":"s","volume":"Mpc^3/sr/redshift"},"ordering":"parameter major, query minor","resource_interpretation":{"maximum_native_output_bytes":"conservative retained input/output plus native scratch estimate; excludes caller input and allocator overhead","maximum_total_segment_visits":"one analytic admitted segment budget across successful and failed rows; no quadrature callbacks"},"requested_outputs":[{"id":"piecewise_background_slots","required":true,"numerical_gate":"required","inference_gate":"not_applicable"}]});
    if std::env::var_os("COSMOLOGY_TEST_PANIC").is_some() {
        panic!("injected Rust panic");
    }
    let result = piecewise::evaluate(&request.parameters, &request.queries, &request.policy)?;
    Ok((spec, result))
}

pub(crate) fn execute(input: &[u8], store: &Path) -> Result<crate::outcome::Outcome, String> {
    let (spec, output) = calculation(input, store)?;
    let arithmetic = json!("binary64_storage_longdouble_analytic_intermediate");
    let resources = spec["policy"].clone();
    Ok(crate::outcome::Outcome::scientific(
        spec,
        output,
        "piecewise_background_slots",
        json!("analytic_piecewise_segments"),
        arithmetic,
        resources,
        "bounded native fixed-five-bin comparisons; request applicability not established",
    ))
}
