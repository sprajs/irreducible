//! Local structural adapters and retained acquisition; scientific policy is C++.
use crate::{
    bridge::{ObservationInput, ObservationMetadata, prepare_observations},
    ingestion::{gaussian_fixture, pantheon_covariance, pantheon_plus, read_asset},
    records::{hash, publish},
};
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};
use std::path::{Path, PathBuf};
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Resources {
    maximum_asset_bytes: usize,
    maximum_rows: u64,
    maximum_matrix_elements: u64,
    maximum_string_bytes: u64,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Request {
    schema_version: u32,
    operation: String,
    table: PathBuf,
    uncertainty: Option<PathBuf>,
    metadata: ObservationMetadata,
    selection: String,
    resources: Resources,
    #[serde(default)]
    calibration_provenance: String,
    #[serde(default)]
    dependence_provenance: String,
    #[serde(default)]
    quality_dictionary: String,
    #[serde(default)]
    ordering_provenance: String,
    #[serde(default)]
    source_selection: Vec<u8>,
}
pub(crate) fn execute(input: &[u8], store: &Path) -> Result<crate::outcome::Outcome, String> {
    let (spec, output) = execute_calculation(input, store, |_, _, _| Ok(None))?;
    let resources = spec["resource_budget"].clone();
    Ok(crate::outcome::Outcome::structural(
        spec,
        output,
        "prepared_observations",
        json!("typed_observation_preparation"),
        json!("binary64_source_storage_no_scientific_transformation"),
        resources,
        "structural typed observations and source identity; no inference acceptance",
    ))
}
// Reuse exactly the same structural acquisition and immutable source retention.
pub(crate) fn execute_calculation(
    input: &[u8],
    store: &Path,
    calculation: impl FnOnce(&ObservationInput, &str, [u64; 3]) -> Result<Option<Value>, String>,
) -> Result<(Value, Value), String> {
    let r: Request = serde_json::from_slice(input).map_err(|e| e.to_string())?;
    if r.schema_version != 1 || r.operation != "observations.prepare.v1" {
        return Err("UNSUPPORTED_SPECIFICATION".into());
    }
    let limits = &r.resources;
    if limits.maximum_asset_bytes > 256 * 1024 * 1024
        || limits.maximum_rows > 10000
        || limits.maximum_matrix_elements > 25000000
        || limits.maximum_string_bytes > 16 * 1024 * 1024
    {
        return Err("RESOURCE_LIMIT".into());
    }
    let asset = read_asset(&r.table, limits.maximum_asset_bytes)?;
    publish(&store.join("objects").join(asset.sha256()), asset.bytes())?;
    let table = match r.metadata.profile.as_str() {
        "pantheon_plus_released_v1" => pantheon_plus(asset, limits.maximum_rows as usize)?,
        "gaussian_fixture_v1" => gaussian_fixture(asset, limits.maximum_rows as usize)?,
        _ => return Err("UNAVAILABLE_SOURCE_PROFILE".into()),
    };
    let mut matrix_values = vec![];
    let mut axes = vec![];
    let mut matrix_digest = String::new();
    let mut reader_ordering_provenance = String::new();
    if let Some(path) = &r.uncertainty {
        let asset = read_asset(path, limits.maximum_asset_bytes)?;
        publish(&store.join("objects").join(asset.sha256()), asset.bytes())?;
        let matrix = pantheon_covariance(
            asset,
            table.measurement_ids.clone(),
            limits.maximum_matrix_elements as usize,
        )?;
        if matrix.dimension != table.values.len() {
            return Err("AXIS_SHAPE".into());
        }
        reader_ordering_provenance = matrix.ordering_provenance;
        matrix_digest = matrix.asset.sha256().to_owned();
        matrix_values = matrix.values;
        axes = matrix.axis_ids;
    }
    let n = table.values.len();
    let nz = table.zhd.len();
    let source = ObservationInput {
        metadata: r.metadata,
        table_sha256: table.asset.sha256().to_owned(),
        uncertainty_sha256: matrix_digest,
        calibration_provenance: r.calibration_provenance,
        dependence_provenance: r.dependence_provenance,
        quality_dictionary: format!(
            "ASCII adapter: no decoded quality column; zero placeholders are not all-clear flags; supplied dictionary declaration: {}",
            r.quality_dictionary
        ),
        ordering_provenance: r.ordering_provenance,
        measurement_ids: table.measurement_ids,
        event_ids: table.event_ids,
        uncertainty_axis_ids: axes,
        values: table.values,
        zhd: table.zhd,
        zcmb: table.zcmb,
        zhel: table.zhel,
        uncertainty_matrix: matrix_values,
        missing: vec![0; n],
        zhd_missing: vec![0; nz],
        zcmb_missing: vec![0; nz],
        zhel_missing: vec![0; nz],
        source_selection: r.source_selection,
        quality: table.quality,
    };
    // Canonical identity uses acquired immutable bytes, not mutable pathname identity.
    let spec = json!({"schema_version":1,"operation":r.operation,"equation_id":"D01/prepared-observations/v1","table_digest":source.table_sha256,"uncertainty_digest":source.uncertainty_sha256,"metadata":source.metadata,"selection":r.selection,
 "calibration_provenance":source.calibration_provenance,"dependence_provenance":source.dependence_provenance,"quality_dictionary":source.quality_dictionary,"ordering_provenance":source.ordering_provenance,"reader_ordering_provenance":reader_ordering_provenance,"source_selection":source.source_selection,
 "resource_budget":{"maximum_asset_bytes":limits.maximum_asset_bytes,"maximum_rows":limits.maximum_rows,"maximum_matrix_elements":limits.maximum_matrix_elements,"maximum_string_bytes":limits.maximum_string_bytes},
 "requested_outputs":[{"id":"prepared_observations","required":true,"numerical_gate":"required","inference_gate":"not_applicable"}]});
    if std::env::var_os("COSMOLOGY_TEST_PANIC").is_some() {
        panic!("injected Rust panic");
    }
    let prepared = match prepare_observations(
        &source,
        &r.selection,
        [
            limits.maximum_rows,
            limits.maximum_matrix_elements,
            limits.maximum_string_bytes,
        ],
    ) {
        Ok(value) => value,
        Err(error) => return Err(error),
    };
    #[derive(Serialize)]
    struct RetainedMatrix<'a> {
        schema_version: u32,
        values: &'a [f64],
        axis_ids: &'a [String],
        uncertainty: &'a str,
        unit: &'a str,
        ordering_provenance: &'a str,
        source_asset_digest: &'a str,
    }
    let retained = RetainedMatrix {
        schema_version: 1,
        values: &prepared.retained_matrix,
        axis_ids: &source.uncertainty_axis_ids,
        uncertainty: &source.metadata.uncertainty,
        unit: &source.metadata.uncertainty_unit,
        ordering_provenance: &source.ordering_provenance,
        source_asset_digest: &source.uncertainty_sha256,
    };
    let matrix_bytes = serde_json::to_vec(&retained).map_err(|e| e.to_string())?;
    let matrix_artifact_digest = hash(&matrix_bytes);
    publish(
        &store.join("objects").join(&matrix_artifact_digest),
        &matrix_bytes,
    )?;
    // No row deletion, uncertainty projection, conversion, or dependence assumption.
    let output = json!({"kind":"finite","metadata":source.metadata,"table_digest":source.table_sha256,"uncertainty_digest":source.uncertainty_sha256,
 "source_missing":{"values":source.missing,"zhd":source.zhd_missing,"zcmb":source.zcmb_missing,"zhel":source.zhel_missing},"source_nonfinite":{"values":source.values.iter().map(|v|!v.is_finite()).collect::<Vec<_>>(),"zhd":source.zhd.iter().map(|v|!v.is_finite()).collect::<Vec<_>>(),"zcmb":source.zcmb.iter().map(|v|!v.is_finite()).collect::<Vec<_>>(),"zhel":source.zhel.iter().map(|v|!v.is_finite()).collect::<Vec<_>>()},"source_value_encoding":"nonfinite binary64 slots encode JSON null and have separate explicit nonfinite masks; raw original fields retained","measurement_ids":source.measurement_ids,"event_ids":source.event_ids,"quality":source.quality,"values":prepared.source_values,"zhd":source.zhd,"zcmb":source.zcmb,"zhel":source.zhel,
 "uncertainty_axis_ids":source.uncertainty_axis_ids,"full_uncertainty_matrix":{"object_digest":matrix_artifact_digest,"elements":prepared.retained_matrix.len(),"axis_ids":source.uncertainty_axis_ids,"unit":source.metadata.uncertainty_unit,"retention":"full matrix unchanged; no projection or repair"},"selected_mask":prepared.mask,"selected_source_indices":prepared.source_indices,
 "calibration_provenance":source.calibration_provenance,"dependence_provenance":source.dependence_provenance,"quality_dictionary":source.quality_dictionary,"ordering_provenance":source.ordering_provenance,
 "quality_flags_available":false,"quality_placeholder_semantics":"no decoded quality column; zeros are structural placeholders, not all-clear flags","original_columns":table.original_columns,"original_fields":table.original_fields,"raw_assets_retained":true,"reader_ordering_provenance":reader_ordering_provenance,"ordering_verified_from_unlabeled_matrix":false,"inference_independence":"not_asserted"});
    let calculation = calculation(
        &source,
        &r.selection,
        [
            limits.maximum_rows,
            limits.maximum_matrix_elements,
            limits.maximum_string_bytes,
        ],
    )?;
    let output = if let Some(result) = calculation {
        json!({"kind":result["kind"],"observations":output,"calculation":result})
    } else {
        output
    };
    Ok((spec, output))
}
