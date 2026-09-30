//! Operation-owned radiometric outputs and records, without duplicate equations.
use crate::{
    bridge::photometry::{self, Input, Output, Policy},
    outcome::{NumericalCheck, Outcome, OutputCheck},
};
use serde::Deserialize;
use serde_json::json;
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Request {
    schema_version: u32,
    operation: String,
    source_model: String,
    propagation: String,
    constants_id: String,
    inputs: Vec<Input>,
    requested_outputs: Vec<Output>,
    resource_policy: Policy,
}
pub(crate) fn execute(bytes: &[u8]) -> Result<Outcome, String> {
    let request: Request = serde_json::from_slice(bytes).map_err(|e| e.to_string())?;
    if request.schema_version != 2
        || request.operation != "photometry.predict"
        || request.source_model != "constant_rest_luminosity_rectangular_band"
        || request.propagation != "isotropic_luminosity_distance_standard_redshift"
        || request.constants_id != "si_2019_radiometric_definitions"
    {
        return Err("UNSUPPORTED_SPECIFICATION".into());
    }
    let calculation = photometry::evaluate(
        &request.inputs,
        &request.requested_outputs,
        &request.resource_policy,
    )?;
    let outputs:Vec<_>=calculation.checks.iter().map(|(id,_)|json!({"id":id,"required":true,"numerical_gate":"required","inference_gate":"not_applicable"})).collect();
    let spec = json!({"schema_version":2,"operation":request.operation,"source_model":calculation.model,"propagation":calculation.propagation,"constants_id":calculation.constants,
        "inputs":request.inputs,"requested_outputs":outputs,"resource_policy":request.resource_policy,"source_role":"analytic_synthetic_control","time_convention":"observer exposure; steady rest spectrum; standard redshift time mapping","transmission_role":"constant optical transmission, not quantum efficiency or electronic gain"});
    let resources = json!({"maximum_rows":request.resource_policy.maximum_rows,"maximum_native_bytes":request.resource_policy.maximum_native_bytes,"compute_threads":1,"byte_scope":"successful combined native/wrapper allocated input/result payload; excludes fixed bounded empty work-limit diagnostic owner (also allocated at quota zero), borrowed caller storage, scalar call frames, allocator overhead and RSS"});
    let mut outcome = Outcome::scientific(
        spec,
        calculation.output,
        "photometry",
        json!("analytic_rectangular_band_constant_rest_luminosity"),
        json!(calculation.arithmetic),
        resources,
        "named deterministic numerical controls; no stochastic recovery or request qualification",
    );
    outcome.outputs=calculation.checks.into_iter().map(|(id,passed)|OutputCheck{check_kind:"numerical_contract",id,required:true,numerical:if passed{NumericalCheck::ChecksPassed}else{NumericalCheck::Failed},inference:"not_applicable",interpretation:"unqualified",evidence:vec![],validation_coverage:"named analytic/frequency-quadrature controls; source and detector assumptions conditional"}).collect();
    Ok(outcome)
}
