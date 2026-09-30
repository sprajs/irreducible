//! One structural acquisition and immutable native source; no selection/factor.
use crate::{
    bridge::{
        ObservationInput, ObservationMetadata,
        observations::{Prepared, prepare_shared_bounded},
    },
    ingestion::{magnitude_table, pantheon_covariance, pantheon_plus, read_asset},
    records::{hash, publish},
};
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};
use sha2::{Digest, Sha256};
use std::path::{Path, PathBuf};

#[derive(Clone, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct Resources {
    maximum_asset_bytes: usize,
    maximum_rows: u64,
    maximum_matrix_elements: u64,
    maximum_string_bytes: u64,
    maximum_preparation_bytes: usize,
}
#[derive(Deserialize, Serialize, PartialEq)]
#[serde(rename_all = "snake_case")]
enum Export {
    SourceValues,
    Covariance,
    OriginalFields,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Request {
    schema_version: u32,
    operation: String,
    table: PathBuf,
    uncertainty: Option<PathBuf>,
    metadata: ObservationMetadata,
    resources: Resources,
    exports: Vec<Export>,
    calibration_provenance: String,
    dependence_provenance: String,
    quality_dictionary: String,
    ordering_provenance: String,
    source_selection: Vec<u8>,
}
pub(crate) struct Acquired {
    pub owner: Prepared,
    pub specification: Value,
    pub output: Value,
    pub preparation_digest: String,
    pub retained_bytes: usize,
}
pub(crate) fn execute(input: &[u8], store: &Path) -> Result<crate::outcome::Outcome, String> {
    let acquired = acquire(input, store, 1 << 30)?;
    Ok(crate::outcome::Outcome::structural(
        acquired.specification,
        acquired.output,
        "prepared_observations",
        json!("immutable_typed_source_preparation"),
        json!("binary64_source_storage"),
        json!({"retained_payload_bytes":acquired.retained_bytes}),
        "structural source contract; no covariance SPD or scientific lineage qualification",
    ))
}
fn bits_digest(values: &[f64]) -> String {
    let mut digest = Sha256::new();
    for value in values {
        digest.update(value.to_le_bytes());
    }
    format!("{:x}", digest.finalize())
}
// Includes bounded Rust ASCII/token storage and native construction together.
// Conservative supported container envelope; caller request frame, allocator
// metadata and RSS are excluded. Checked before reading either source asset.
fn acquisition_bounds(input: &[u8]) -> Result<(usize, usize, usize), String> {
    let r: Request = serde_json::from_slice(input).map_err(|e| e.to_string())?;
    if r.schema_version != 2 || r.operation != "observations.prepare" {
        return Err("UNSUPPORTED_SPECIFICATION".into());
    }
    if r.resources.maximum_asset_bytes > 256 * 1024 * 1024
        || r.resources.maximum_rows > 10000
        || r.resources.maximum_matrix_elements > 25000000
        || r.resources.maximum_string_bytes > 16 * 1024 * 1024
    {
        return Err("HARD_RESOURCE_LIMIT".into());
    }
    for (index, export) in r.exports.iter().enumerate() {
        if r.exports[..index].contains(export) {
            return Err("DUPLICATE_EXPORT".into());
        }
    }
    let table = usize::try_from(std::fs::metadata(&r.table).map_err(|_| "INPUT_READ")?.len())
        .map_err(|_| "RESOURCE_LIMIT")?;
    let covariance = match &r.uncertainty {
        Some(path) => usize::try_from(std::fs::metadata(path).map_err(|_| "INPUT_READ")?.len())
            .map_err(|_| "RESOURCE_LIMIT")?,
        None => 0,
    };
    if table > r.resources.maximum_asset_bytes || covariance > r.resources.maximum_asset_bytes {
        return Err("ASSET_BYTE_LIMIT".into());
    }
    let terms = [
        (table, 64usize),
        (covariance, 2),
        (
            usize::try_from(r.resources.maximum_matrix_elements).map_err(|_| "RESOURCE_LIMIT")?,
            32,
        ),
        (
            usize::try_from(r.resources.maximum_rows).map_err(|_| "RESOURCE_LIMIT")?,
            2048,
        ),
        (
            usize::try_from(r.resources.maximum_string_bytes).map_err(|_| "RESOURCE_LIMIT")?,
            64,
        ),
    ];
    let mut total = 4096usize;
    for (count, width) in terms {
        total = total
            .checked_add(count.checked_mul(width).ok_or("RESOURCE_LIMIT")?)
            .ok_or("RESOURCE_LIMIT")?;
    }
    if total > r.resources.maximum_preparation_bytes {
        return Err("PREPARATION_BYTE_LIMIT".into());
    }
    Ok((total, table, covariance))
}
pub(crate) fn peak_bound(input: &[u8]) -> Result<usize, String> {
    acquisition_bounds(input).map(|bounds| bounds.0)
}
pub(crate) fn acquire(
    input: &[u8],
    store: &Path,
    remaining_bytes: usize,
) -> Result<Acquired, String> {
    let (peak, table_limit, covariance_limit) = acquisition_bounds(input)?;
    if peak > remaining_bytes {
        return Err("PREPARATION_BYTE_LIMIT".into());
    }
    let r: Request = serde_json::from_slice(input).map_err(|e| e.to_string())?;
    if r.schema_version != 2 || r.operation != "observations.prepare" {
        return Err("UNSUPPORTED_SPECIFICATION".into());
    }
    if r.resources.maximum_asset_bytes > 256 * 1024 * 1024
        || r.resources.maximum_rows > 10000
        || r.resources.maximum_matrix_elements > 25000000
        || r.resources.maximum_string_bytes > 16 * 1024 * 1024
    {
        return Err("HARD_RESOURCE_LIMIT".into());
    }
    let raw = read_asset(&r.table, table_limit)?;
    publish(&store.join("objects").join(raw.sha256()), raw.bytes())?;
    let table = match r.metadata.profile.as_str() {
        "pantheon_plus_released_v1" => pantheon_plus(raw, r.resources.maximum_rows as usize)?,
        "typed_magnitude_covariance" | "gaussian_fixture_v1" => {
            magnitude_table(raw, r.resources.maximum_rows as usize)?
        }
        _ => return Err("UNAVAILABLE_SOURCE_PROFILE".into()),
    };
    let mut covariance = Vec::new();
    let mut axes = Vec::new();
    let mut covariance_identity = String::new();
    let mut reader_ordering = String::new();
    if let Some(path) = &r.uncertainty {
        let raw = read_asset(path, covariance_limit)?;
        publish(&store.join("objects").join(raw.sha256()), raw.bytes())?;
        let matrix = pantheon_covariance(
            raw,
            table.measurement_ids.clone(),
            r.resources.maximum_matrix_elements as usize,
        )?;
        if matrix.dimension != table.values.len() {
            return Err("AXIS_SHAPE".into());
        }
        covariance = matrix.values;
        axes = matrix.axis_ids;
        covariance_identity = matrix.asset.sha256().to_owned();
        reader_ordering = matrix.ordering_provenance;
    }
    let n = table.values.len();
    let nz = table.zhd.len();
    let source = ObservationInput {
        metadata: r.metadata,
        table_sha256: table.asset.sha256().to_owned(),
        uncertainty_sha256: covariance_identity,
        calibration_provenance: r.calibration_provenance,
        dependence_provenance: r.dependence_provenance,
        quality_dictionary: format!(
            "ASCII adapter has no decoded quality flags; zero placeholders are not all-clear. Supplied declaration: {}",
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
        uncertainty_matrix: covariance,
        missing: vec![0; n],
        zhd_missing: vec![0; nz],
        zcmb_missing: vec![0; nz],
        zhel_missing: vec![0; nz],
        source_selection: r.source_selection,
        quality: table.quality,
    };
    let owner = prepare_shared_bounded(
        &source,
        [
            r.resources.maximum_rows,
            r.resources.maximum_matrix_elements,
            r.resources.maximum_string_bytes,
        ],
        peak,
    )?;
    let view = owner.source_view()?;
    let values_identity = bits_digest(view.values()?);
    let covariance_bits_identity = bits_digest(view.covariance()?);
    let specification = json!({"schema_version":2,"operation":"observations.prepare","source_type":source.metadata,
        "table_identity":source.table_sha256,"covariance_identity":source.uncertainty_sha256,
        "values_f64le_digest":values_identity,"covariance_f64le_digest":covariance_bits_identity,
        "ordered_ids":source.measurement_ids,"covariance_axis_ids":source.uncertainty_axis_ids,
        "calibration_provenance":source.calibration_provenance,"dependence_provenance":source.dependence_provenance,
        "ordering_provenance":source.ordering_provenance,"quality_dictionary":source.quality_dictionary,"source_selection":source.source_selection,
        "resources":r.resources,"requested_exports":r.exports,"requested_outputs":[{"id":"prepared_observations","required":true,"numerical_gate":"required","inference_gate":"not_applicable"}]});
    let mut output = json!({"kind":"finite","source_row_count":n,"source_type":source.metadata,"table_digest":source.table_sha256,
        "uncertainty_digest":source.uncertainty_sha256,"values_f64le_digest":values_identity,"covariance_f64le_digest":covariance_bits_identity,
        "measurement_ids":source.measurement_ids,"event_ids":source.event_ids,"covariance_axis_ids":source.uncertainty_axis_ids,
        "calibration_provenance":source.calibration_provenance,"dependence_provenance":source.dependence_provenance,
        "ordering_provenance":source.ordering_provenance,"reader_ordering_provenance":reader_ordering,
        "ordering_verified_from_unlabeled_matrix":false,"quality_flags_available":false,"raw_assets_retained":true,"inference_independence":"not_asserted"});
    if r.exports.contains(&Export::SourceValues) {
        output["values"] = json!(view.values()?);
        output["source_nonfinite"] = json!({"values":source.values.iter().map(|x|u8::from(!x.is_finite())).collect::<Vec<_>>(),"zhd":source.zhd.iter().map(|x|u8::from(!x.is_finite())).collect::<Vec<_>>(),"zcmb":source.zcmb.iter().map(|x|u8::from(!x.is_finite())).collect::<Vec<_>>(),"zhel":source.zhel.iter().map(|x|u8::from(!x.is_finite())).collect::<Vec<_>>()});
        output["source_missing"] = json!({"values":source.missing,"zhd":source.zhd_missing,"zcmb":source.zcmb_missing,"zhel":source.zhel_missing});
    }
    if r.exports.contains(&Export::Covariance) {
        let artifact = json!({"schema_version":2,"values":view.covariance()?,"axis_ids":source.uncertainty_axis_ids,
            "unit":source.metadata.uncertainty_unit,"source_asset_digest":source.uncertainty_sha256,"ordering_provenance":source.ordering_provenance});
        let bytes = serde_json::to_vec(&artifact).map_err(|e| e.to_string())?;
        let digest = hash(&bytes);
        publish(&store.join("objects").join(&digest), &bytes)?;
        output["full_uncertainty_matrix"] = json!({"object_digest":digest,"elements":view.covariance()?.len(),"axis_ids":source.uncertainty_axis_ids});
    }
    if r.exports.contains(&Export::OriginalFields) {
        output["original_columns"] = json!(table.original_columns);
        output["original_fields"] = json!(table.original_fields);
    }
    drop(view);
    let preparation_digest = hash(&serde_json::to_vec(&specification).map_err(|e| e.to_string())?);
    let metadata_bytes = serde_json::to_vec(&specification)
        .map_err(|e| e.to_string())?
        .len()
        .checked_add(
            serde_json::to_vec(&output)
                .map_err(|e| e.to_string())?
                .len(),
        )
        .ok_or("RESOURCE_LIMIT")?;
    let retained_bytes = owner
        .retained_bytes()?
        .checked_add(metadata_bytes.checked_mul(64).ok_or("RESOURCE_LIMIT")?)
        .and_then(|n| n.checked_add(4096))
        .ok_or("RESOURCE_LIMIT")?;
    if retained_bytes > peak {
        return Err("PREPARATION_BYTE_LIMIT".into());
    }
    Ok(Acquired {
        owner,
        specification,
        output,
        preparation_digest,
        retained_bytes,
    })
}
