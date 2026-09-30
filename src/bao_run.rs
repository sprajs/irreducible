//! Acquire bounded source bytes and retain specifications; no Rust physics.
use crate::{
    bao_ingestion::Source,
    bridge::{BaoModelRequest, BaoPolicyRequest, bao_evaluate},
    records::{hash, publish},
};
use serde::Deserialize;
use serde_json::{Value, json};
use std::path::Path;
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Request {
    schema_version: u32,
    operation: String,
    source: Source,
    prepare_policy: BaoPolicyRequest,
    models: Vec<BaoModelRequest>,
    evaluation_policy: BaoPolicyRequest,
}
pub(crate) fn execute(input: &[u8], store: &Path) -> Result<(Value, Value), String> {
    let r: Request = serde_json::from_slice(input).map_err(|e| e.to_string())?;
    if r.schema_version != 1 || r.operation != "bao.gaussian_batch.v1" {
        return Err("UNSUPPORTED_SPECIFICATION".into());
    }
    if r.prepare_policy.maximum_rows > 4096
        || r.prepare_policy.maximum_matrix_elements > 16777216
        || r.models.len() > 64
    {
        return Err("RESOURCE_LIMIT".into());
    }
    let source = r.source.decode(
        r.prepare_policy.maximum_rows as usize,
        r.prepare_policy.maximum_matrix_elements as usize,
    )?;
    let objects = store.join("objects");
    for asset in [&source.table_asset, &source.covariance_asset] {
        if hash(asset.bytes()) != asset.sha256() {
            return Err("ASSET_IDENTITY_CHANGED".into());
        }
        publish(&objects.join(asset.sha256()), asset.bytes())?;
    }
    let matrix = json!({"kind":"full_declared_bao_covariance","dimension":source.rows.len(),"axis_ids":source.covariance_axis_ids,"ordering_provenance":source.ordering_provenance,"ordering_assessment":"supplied declaration; not independently verified from unlabeled bytes","unit":"ratio_squared","source_asset_sha256":source.covariance_asset.sha256(),"values":source.covariance,"nonfinite":source.covariance.iter().map(|x|!x.is_finite()).collect::<Vec<_>>()});
    let matrix_bytes = serde_json::to_vec(&matrix).map_err(|e| e.to_string())?;
    let matrix_digest = hash(&matrix_bytes);
    publish(&objects.join(&matrix_digest), &matrix_bytes)?;
    let resolved = json!({"profile":match r.source {Source::Released{..}=>"desi_dr2_all_gccomb_13_v1",Source::Synthetic{..}=>"synthetic_inline_v1"},"role":source.role,"unit":"one","covariance_unit":"ratio_squared","rows":source.rows,"source_z_nonfinite":source.rows.iter().map(|r|!r.z.is_finite()).collect::<Vec<_>>(),"source_value_nonfinite":source.rows.iter().map(|r|!r.value.is_finite()).collect::<Vec<_>>(),"original_fields":source.original_fields,"table_asset_sha256":source.table_asset.sha256(),"covariance_asset_sha256":source.covariance_asset.sha256(),"full_covariance":{"object_digest":matrix_digest,"elements":source.covariance.len(),"axis_ids":source.covariance_axis_ids},"ordering_provenance":source.ordering_provenance,"axis_assessment":"supplied declaration; not independently verified","calibration_provenance":source.calibration_provenance,"dependence_provenance":source.dependence_provenance,"redshift_convention":source.redshift_convention,"ruler_convention":source.ruler_convention,"computational_h0_convention":source.computational_h0_convention,"computational_h0_km_s_mpc":70});
    if std::env::var_os("COSMOLOGY_TEST_PANIC").is_some() {
        panic!("injected Rust panic");
    }
    let calculation = bao_evaluate(&source, &r.models, &r.prepare_policy, &r.evaluation_policy)?;
    let success = calculation["numerical_status"] == "checks_passed";
    let spec = json!({"schema_version":1,"operation":r.operation,"source":resolved,"models":r.models,"prepare_policy":r.prepare_policy,"evaluation_policy":r.evaluation_policy,"target":"normalized Gaussian on ordered free-ruler fitted distance-ratio coordinates; no joint-probe inference","requested_outputs":[{"id":"bao_densities","required":true,"numerical_gate":"required","inference_gate":"not_applicable"}]});
    Ok((
        spec,
        json!({"kind":if success{"finite"}else{"failure"},"error_id":if success{None}else{Some("BAO_NUMERICAL_EVALUATION_FAILURE")},"source":resolved,"calculation":calculation}),
    ))
}
