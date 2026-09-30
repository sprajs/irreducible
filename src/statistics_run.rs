//! Request resolution only; C++ owns the retained Gaussian and batch equations.
use crate::{
    bridge::{ProperPrior, gaussian_batch},
    observation_run,
};
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};
use std::path::Path;
#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
struct Prior {
    response: Vec<f64>,
    ordered_ids: Vec<String>,
    mean: f64,
    variance: f64,
    latent_identity: String,
    independence_declared: bool,
}
#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
struct Request {
    schema_version: u32,
    operation: String,
    observations: Value,
    mode: String,
    residual_unit: String,
    response_unit: String,
    ordered_ids: Vec<String>,
    residuals: Vec<Vec<f64>>,
    response: Option<Vec<f64>>,
    proper_prior: Option<Prior>,
    selection: String,
    maximum_preparation_bytes: u64,
    maximum_evaluation_bytes: u64,
    maximum_batch_elements: u64,
    maximum_forward_sensitivity: f64,
}
fn calculation(input: &[u8], store: &Path) -> Result<(Value, Value, Value), String> {
    let r: Request = serde_json::from_slice(input).map_err(|e| e.to_string())?;
    if r.schema_version != 2 || r.operation != "statistics.gaussian" {
        return Err("UNSUPPORTED_SPECIFICATION".into());
    }
    let profile = match r.mode.as_str() {
        "normalized_density" => false,
        "profile_offset_score" => true,
        _ => return Err("UNKNOWN_GAUSSIAN_MODE".into()),
    };
    if profile != r.response.is_some() || (profile && r.proper_prior.is_some()) {
        return Err("INCOMPATIBLE_NUISANCE_TREATMENT".into());
    }
    if r.response_unit != "one" || r.observations["metadata"]["unit"] != r.residual_unit {
        return Err("INCOMPATIBLE_RESIDUAL_UNITS".into());
    }
    if !r.maximum_forward_sensitivity.is_finite() || r.maximum_forward_sensitivity <= 0.0 {
        return Err("INVALID_NUMERICAL_BUDGET".into());
    }
    let nested = serde_json::to_vec(&r.observations).map_err(|e| e.to_string())?;
    let acquired = observation_run::acquire(&nested, store, 1 << 30)?;
    let prior = r.proper_prior.as_ref().map(|p| ProperPrior {
        response: &p.response,
        ordered_ids: &p.ordered_ids,
        mean: p.mean,
        variance: p.variance,
        latent_identity: &p.latent_identity,
        independence_declared: p.independence_declared,
    });
    let resources = &r.observations["resources"];
    let policy = [
        "maximum_rows",
        "maximum_matrix_elements",
        "maximum_string_bytes",
    ]
    .map(|field| {
        resources[field]
            .as_u64()
            .ok_or("INVALID_OBSERVATION_RESOURCES")
    })
    .into_iter()
    .collect::<Result<Vec<_>, _>>()?;
    let effective_preparation_bytes = r.maximum_preparation_bytes.min(
        (1usize << 30)
            .checked_sub(acquired.retained_bytes)
            .ok_or("RETAINED_BYTE_LIMIT")? as u64,
    );
    let output = gaussian_batch(
        &acquired.owner,
        &r.selection,
        &r.residual_unit,
        [policy[0], policy[1], policy[2]],
        &r.residuals,
        &r.ordered_ids,
        r.response.as_deref(),
        prior,
        r.maximum_batch_elements,
        r.maximum_forward_sensitivity,
        effective_preparation_bytes,
        r.maximum_evaluation_bytes.min(
            (1usize << 30)
                .checked_sub(acquired.retained_bytes)
                .ok_or("RETAINED_BYTE_LIMIT")? as u64,
        ),
    )?;
    let observations = acquired.specification;
    let spec = json!({"schema_version":2,"operation":r.operation,"equation_id":"F03/prepared-Gaussian/v1","observations":observations,"preparation_digest":acquired.preparation_digest,"selection":r.selection,"maximum_preparation_bytes":r.maximum_preparation_bytes,"maximum_evaluation_bytes":r.maximum_evaluation_bytes,"mode":r.mode,"residual_unit":r.residual_unit,"response_unit":r.response_unit,"ordered_ids":r.ordered_ids,"residuals":r.residuals,"response":r.response,"proper_prior":r.proper_prior,"maximum_batch_elements":r.maximum_batch_elements,"maximum_forward_sensitivity":r.maximum_forward_sensitivity,"requested_outputs":[{"id":"gaussian_batch","required":true,"numerical_gate":"required","inference_gate":"not_applicable"}]});
    let actual_resources = json!({"aggregate_payload_ceiling":1usize << 30,"shared_source_retained_bytes":acquired.retained_bytes,"effective_preparation_bytes":effective_preparation_bytes,"effective_evaluation_bytes":r.maximum_evaluation_bytes.min(((1usize<<30)-acquired.retained_bytes) as u64)});
    Ok((spec, output, actual_resources))
}

pub(crate) fn execute(input: &[u8], store: &Path) -> Result<crate::outcome::Outcome, String> {
    let (spec, output, actual_resources) = calculation(input, store)?;
    let arithmetic = json!("binary64_legacy_v1");
    let resources = json!({"actual":actual_resources,"maximum_preparation_bytes":spec["maximum_preparation_bytes"],"maximum_evaluation_bytes":spec["maximum_evaluation_bytes"],"maximum_batch_elements":spec["maximum_batch_elements"],"maximum_forward_sensitivity":spec["maximum_forward_sensitivity"]});
    let method = spec["mode"].clone();
    Ok(crate::outcome::Outcome::scientific(
        spec,
        output,
        "gaussian_batch",
        method,
        arithmetic,
        resources,
        "bounded native Gaussian and interface tests; request applicability not established",
    ))
}
