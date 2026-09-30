//! Request resolution only; C++ owns the retained Gaussian and batch equations.
use crate::{
    bridge::{gaussian_batch, ProperPrior},
    observation_run,
};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
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
    maximum_batch_elements: u64,
    maximum_forward_sensitivity: f64,
}
pub(crate) fn execute(input: &[u8], store: &Path) -> Result<(Value, Value), String> {
    let r: Request = serde_json::from_slice(input).map_err(|e| e.to_string())?;
    if r.schema_version != 1 || r.operation != "statistics.gaussian_batch.v1" {
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
    let (observations, output) =
        observation_run::execute_calculation(&nested, store, |source, selection, policy| {
            let prior = r.proper_prior.as_ref().map(|p| ProperPrior {
                response: &p.response,
                ordered_ids: &p.ordered_ids,
                mean: p.mean,
                variance: p.variance,
                latent_identity: &p.latent_identity,
                independence_declared: p.independence_declared,
            });
            gaussian_batch(
                source,
                selection,
                policy,
                &r.residuals,
                &r.ordered_ids,
                r.response.as_deref(),
                prior,
                r.maximum_batch_elements,
                r.maximum_forward_sensitivity,
            )
            .map(Some)
        })?;
    let spec = json!({"schema_version":1,"operation":r.operation,"equation_id":"F03/prepared-Gaussian/v1","observations":observations,"mode":r.mode,"residual_unit":r.residual_unit,"response_unit":r.response_unit,"ordered_ids":r.ordered_ids,"residuals":r.residuals,"response":r.response,"proper_prior":r.proper_prior,"maximum_batch_elements":r.maximum_batch_elements,"maximum_forward_sensitivity":r.maximum_forward_sensitivity,"requested_outputs":[{"id":"gaussian_batch","required":true,"numerical_gate":"required","inference_gate":"not_applicable"}]});
    Ok((spec, output))
}
