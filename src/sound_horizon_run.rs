use crate::{
    bridge::sound_horizon::{self, Point, Policy},
    outcome::{NumericalCheck, Outcome, OutputCheck},
};
use serde::Deserialize;
use serde_json::json;
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Request {
    schema_version: u32,
    operation: String,
    points: Vec<Point>,
    numerical_policy: Policy,
}
pub(crate) fn execute(input: &[u8]) -> Result<Outcome, String> {
    let r: Request = serde_json::from_slice(input).map_err(|e| e.to_string())?;
    if r.schema_version != 2 || r.operation != "cosmology.sound_horizon" {
        return Err("UNSUPPORTED_SPECIFICATION".into());
    }
    let (output, passed) = sound_horizon::evaluate(&r.points, &r.numerical_policy)?;
    let specification = json!({"schema_version":2,"operation":r.operation,"points":r.points,"numerical_policy":r.numerical_policy,
 "model_id":output["metadata"]["model_id"],"equation_id":output["metadata"]["equation_id"],
 "constants_id":output["metadata"]["constants_id"],"coordinate":output["metadata"]["coordinate"],"requested_outputs":["sound_horizon_mpc"],
 "scope":"supplied drag redshift; massless radiation and pressureless matter; no predicted thermal history or BAO likelihood qualification"});
    let resources = json!({"maximum_points":r.numerical_policy.maximum_points,"maximum_total_callbacks":r.numerical_policy.maximum_total_callbacks,"maximum_native_bytes":r.numerical_policy.maximum_native_bytes,
 "byte_scope":"successful native/C-boundary calculation admission payload; fixed bounded empty work-limit diagnostic owner permitted outside this quota; excludes borrowed caller buffers, stack, allocator overhead and RSS"});
    let arithmetic = output["metadata"]["arithmetic_id"].clone();
    let mut outcome = Outcome::scientific(
        specification,
        output,
        "sound_horizon",
        json!("compact_scale_factor_adaptive_simpson"),
        arithmetic,
        resources,
        "named analytic/refinement controls; arbitrary requests remain unqualified",
    );
    outcome.outputs=vec![OutputCheck{check_kind:"numerical_contract",id:"sound_horizon_mpc",required:true,numerical:if passed{NumericalCheck::ChecksPassed}else{NumericalCheck::Failed},inference:"not_applicable",interpretation:"unqualified",evidence:vec![],validation_coverage:"conditional supplied-drag approximation; independent named native budgets separate from runtime qualification"}];
    Ok(outcome)
}
