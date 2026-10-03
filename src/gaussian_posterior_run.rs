//! Strict inline structural resolution; all posterior equations remain native.
use crate::gaussian_input::{Conditioning, Design, Noise, Policy, Prior};
use crate::{
    bridge::gaussian_posterior,
    outcome::{NumericalCheck, Outcome, OutputCheck},
};
use serde::{Deserialize, Serialize};
use serde_json::json;
#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct Request {
    pub schema_version: u32,
    pub operation: String,
    pub source_semantics: String,
    pub noise: Noise,
    pub design: Design,
    pub parameter_prior: Prior,
    pub conditioning: Conditioning,
    pub resource_policy: Policy,
}
pub(crate) fn execute(bytes: &[u8]) -> Result<Outcome, String> {
    if bytes.len() > 16 * 1024 * 1024 {
        return Err("INPUT_LIMIT".into());
    }
    crate::strict_json::validate(bytes).map_err(|_| "INVALID_JSON")?;
    let r: Request = serde_json::from_slice(bytes).map_err(|e| e.to_string())?;
    if r.schema_version != 2
        || r.operation != "statistics.gaussian_posterior"
        || r.source_semantics != "synthetic_controls"
    {
        return Err("UNSUPPORTED_SPECIFICATION".into());
    }
    if r.conditioning
        .vectors
        .lengths
        .iter()
        .any(|&n| n != r.noise.ordered_row_ids.len())
    {
        return Err("INVALID_BATCH_SHAPE".into());
    }
    let c = gaussian_posterior::evaluate(&r)?;
    let mut spec = serde_json::to_value(&r).map_err(|_| "RECORD_ENCODING")?;
    spec["requested_outputs"] = json!([{"id":"posterior_covariance","required":true,"numerical_gate":"required","inference_gate":"not_applicable"},{"id":"posterior_means","required":true,"numerical_gate":"required","inference_gate":"not_applicable"}]);
    spec["batch_meaning"] = json!(
        "family of conditionals under one fixed noise/design/independent proper prior; not a joint parameter posterior or IID generated draws"
    );
    let mut outcome = Outcome::scientific(
        spec,
        c.output,
        "posterior_means",
        c.method,
        c.arithmetic,
        c.resources,
        "named native proper-Gaussian posterior controls; supplied synthetic prior/data applicability unqualified",
    );
    outcome.outputs=[("posterior_covariance",c.covariance_passed),("posterior_means",c.means_passed)].into_iter().map(|(id,pass)|OutputCheck{check_kind:"numerical_contract",id,required:true,numerical:if pass{NumericalCheck::ChecksPassed}else{NumericalCheck::Failed},inference:"not_applicable",interpretation:"unqualified",evidence:vec![],validation_coverage:"original two-dimensional rational/cofactor and native/ABI/CLI parity; inherited scientific ancestry"}).collect();
    Ok(outcome)
}
