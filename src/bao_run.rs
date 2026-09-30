//! Retained BAO acquisition. The structural reader is shared with one-shot use.
use crate::{
    bao_ingestion::Source,
    bridge::bao::{self, PreparationPolicy, Prepared},
    records::{hash, publish},
};
use serde_json::{Value, json};
use std::path::Path;
pub(crate) struct Acquired {
    pub owner: Prepared,
    pub specification: Value,
    pub preparation_digest: String,
    pub retained_bytes: usize,
    pub status: u32,
    pub numerical_status: u32,
    pub effective_policy: PreparationPolicy,
}
// Supported ASCII/container envelope. Caller frame and allocator/RSS excluded.
pub(crate) fn acquisition_bound(
    source: &Source,
    policy: &PreparationPolicy,
) -> Result<(usize, usize), String> {
    if policy.maximum_queries > 4096
        || policy.maximum_matrix_elements > 16777216
        || policy.maximum_string_bytes > 16 * 1024 * 1024
    {
        return Err("HARD_RESOURCE_LIMIT".into());
    }
    let (rows, matrix_elements) = match source {
        Source::Released { .. } => (13usize, 169usize),
        Source::Synthetic {
            rows, covariance, ..
        } => (rows.len(), covariance.len()),
    };
    if rows > usize::try_from(policy.maximum_queries).map_err(|_| "RESOURCE_LIMIT")?
        || matrix_elements
            > usize::try_from(policy.maximum_matrix_elements).map_err(|_| "RESOURCE_LIMIT")?
    {
        return Err("RESOURCE_LIMIT".into());
    }
    let metadata = serde_json::to_vec(source)
        .map_err(|_| "RECORD_ENCODING")?
        .len();
    let asset = match source {
        Source::Released {
            mean,
            covariance,
            maximum_asset_bytes,
            ..
        } => {
            if *maximum_asset_bytes > 16 * 1024 * 1024 {
                return Err("HARD_RESOURCE_LIMIT".into());
            }
            let a = usize::try_from(std::fs::metadata(mean).map_err(|_| "INPUT_READ")?.len())
                .map_err(|_| "RESOURCE_LIMIT")?;
            let b = usize::try_from(
                std::fs::metadata(covariance)
                    .map_err(|_| "INPUT_READ")?
                    .len(),
            )
            .map_err(|_| "RESOURCE_LIMIT")?;
            let size = a.max(b);
            if size > *maximum_asset_bytes {
                return Err("ASSET_BYTE_LIMIT".into());
            }
            size
        }
        Source::Synthetic { .. } => serde_json::to_vec(source)
            .map_err(|_| "RECORD_ENCODING")?
            .len(),
    };
    let bytes = asset
        .checked_mul(128)
        .and_then(|n| n.checked_add(matrix_elements.checked_mul(64)?))
        .and_then(|n| n.checked_add(rows.checked_mul(2048)?))
        .and_then(|n| n.checked_add(metadata.checked_mul(128)?))
        .and_then(|n| n.checked_add(4096))
        .ok_or("RESOURCE_LIMIT")?;
    Ok((bytes, asset))
}
pub(crate) fn acquire(
    mut source: Source,
    requested_policy: &PreparationPolicy,
    mut effective_policy: PreparationPolicy,
    store: &Path,
    remaining: usize,
) -> Result<Acquired, String> {
    let (rust_peak, asset_limit) = acquisition_bound(&source, &effective_policy)?;
    let native = remaining
        .checked_sub(rust_peak)
        .ok_or("PREPARATION_BYTE_LIMIT")?;
    effective_policy.maximum_native_bytes = u64::try_from(native).map_err(|_| "RESOURCE_LIMIT")?;
    if let Source::Released {
        maximum_asset_bytes,
        ..
    } = &mut source
    {
        *maximum_asset_bytes = asset_limit;
    }
    let decoded = source.decode(
        usize::try_from(effective_policy.maximum_queries).map_err(|_| "RESOURCE_LIMIT")?,
        usize::try_from(effective_policy.maximum_matrix_elements).map_err(|_| "RESOURCE_LIMIT")?,
    )?;
    publish(
        &store.join("objects").join(decoded.table_asset.sha256()),
        decoded.table_asset.bytes(),
    )?;
    publish(
        &store
            .join("objects")
            .join(decoded.covariance_asset.sha256()),
        decoded.covariance_asset.bytes(),
    )?;
    let prepared = bao::prepare(&decoded, &effective_policy)?;
    let specification = json!({"operation":"bao.density","source":{"role":decoded.role,"unit":"one","covariance_unit":"ratio_squared",
        "original_fields":decoded.original_fields,"rows":decoded.rows,"covariance_axis_ids":decoded.covariance_axis_ids,"table_identity":decoded.table_asset.sha256(),"covariance_identity":decoded.covariance_asset.sha256(),
        "ordering_provenance":decoded.ordering_provenance,"ordering_assessment":"supplied declaration; not independently verified",
        "calibration_provenance":decoded.calibration_provenance,"dependence_provenance":decoded.dependence_provenance,
        "redshift_convention":decoded.redshift_convention,"ruler_convention":decoded.ruler_convention},"preparation_policy":requested_policy,"native_metadata":prepared.metadata});
    let encoded = serde_json::to_vec(&specification).map_err(|_| "RECORD_ENCODING")?;
    let digest = hash(&encoded);
    publish(&store.join("objects").join(&digest), &encoded)?;
    let charge = prepared
        .retained_bytes
        .checked_add(encoded.len().checked_mul(64).ok_or("RESOURCE_LIMIT")?)
        .and_then(|n| n.checked_add(4096))
        .ok_or("RESOURCE_LIMIT")?;
    if charge > remaining {
        return Err("RETAINED_BYTE_LIMIT".into());
    }
    Ok(Acquired {
        owner: prepared.owner,
        specification,
        preparation_digest: digest,
        retained_bytes: charge,
        status: prepared.status,
        numerical_status: prepared.numerical_status,
        effective_policy,
    })
}
