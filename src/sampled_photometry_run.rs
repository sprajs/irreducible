//! Sampled photometry owns source roles, method identity and group checks.
use crate::{
    bridge::{
        photometry::Output,
        sampled_photometry::{self, Bounded, Exposure, Passband, Policy, Spectrum},
    },
    outcome::{NumericalCheck, Outcome, OutputCheck},
};
use serde::{Deserialize, Serialize};
use serde_json::json;
#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
struct Request {
    schema_version: u32,
    operation: String,
    source_model: String,
    propagation: String,
    constants_id: String,
    spectra: Bounded<Spectrum, 4096>,
    passbands: Bounded<Passband, 4096>,
    exposures: Bounded<Exposure, 65536>,
    requested_outputs: Bounded<Output, 3>,
    resource_policy: Policy,
}
pub(crate) fn execute(bytes: &[u8]) -> Result<Outcome, String> {
    if bytes.len() > 16 * 1024 * 1024 {
        return Err("INPUT_LIMIT".into());
    }
    crate::strict_json::validate(bytes).map_err(|_| "INVALID_JSON")?;
    let r: Request = serde_json::from_slice(bytes).map_err(|e| e.to_string())?;
    if r.schema_version != 2
        || r.operation != "photometry.predict"
        || r.source_model != "piecewise_linear_rest_luminosity_observed_optical_passband"
        || r.propagation != "isotropic_luminosity_distance_standard_redshift"
        || r.constants_id != "si_2019_radiometric_definitions"
    {
        return Err("UNSUPPORTED_SPECIFICATION".into());
    }
    let sampled = sampled_photometry::evaluate(
        &r.spectra,
        &r.passbands,
        &r.exposures,
        &r.requested_outputs,
        &r.resource_policy,
    )?;
    let c = sampled.calculation;
    let mut spec = serde_json::to_value(&r).map_err(|_| "RECORD_ENCODING")?;
    spec["source_model"] = json!(c.model);
    spec["constants_id"] = json!(c.constants);
    spec["propagation"] = json!(c.propagation);
    spec["requested_outputs"]=json!(c.checks.iter().map(|(id,_)|json!({"id":id,"required":true,"numerical_gate":"required","inference_gate":"not_applicable"})).collect::<Vec<_>>());
    spec["interpolation"]=json!("exact declared piecewise linear in wavelength; zero outside finite support; samples do not qualify the true source or response");
    spec["time_convention"] =
        json!("observer exposure; steady rest spectrum; standard redshift time mapping");
    spec["transmission_role"]=json!("observed optical transmission; not detector quantum efficiency or electronic gain; uncertainty excluded");
    let mut resources = serde_json::to_value(&r.resource_policy).map_err(|_| "RECORD_ENCODING")?;
    resources["compute_threads"] = json!(1);
    resources["byte_scope"]=json!("native owner, pooled curve owners/ID bytes/two double arrays, exported curve descriptors, result rows; excludes fixed empty work-limit diagnostic owner, borrowed Rust input and wire descriptors, scalar call frames, allocator overhead and RSS");
    let mut outcome=Outcome::scientific(spec,c.output,"photometry",json!(sampled.method),json!(c.arithmetic),resources,"named deterministic sampled numerical controls; declared source roles/provenance preserved but not qualified; no calibration covariance, stochastic recovery or inference");
    outcome.outputs=c.checks.into_iter().map(|(id,passed)|OutputCheck{check_kind:"numerical_contract",id,required:true,numerical:if passed{NumericalCheck::ChecksPassed}else{NumericalCheck::Failed},inference:"not_applicable",interpretation:"unqualified",evidence:vec![],validation_coverage:"named native constant/linear polynomial and independent frequency-coordinate controls; transport exact parity and ownership controls"}).collect();
    Ok(outcome)
}
