//! Explicit analytic BAO transport; all calculations stay in the shared core.
use super::generated::*;
use crate::bao_ingestion::Decoded;
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};
use std::{ffi::c_void, ptr};
#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct Model {
    pub q: [f64; 5],
    pub h0_rd_km_s: f64,
}
#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct Policy {
    pub arithmetic: String,
    pub include_predictions: bool,
    pub include_residuals: bool,
    pub maximum_models: u64,
    pub maximum_rows: u64,
    pub maximum_matrix_elements: u64,
    pub maximum_string_bytes: u64,
    pub maximum_native_bytes: u64,
    pub maximum_array_elements: u64,
    pub maximum_native_output_bytes: u64,
    pub maximum_queries: u64,
    pub maximum_total_segment_visits: u64,
    pub maximum_forward_sensitivity: f64,
}
fn policy(p: &Policy) -> Result<BaoPiecewisePolicy, String> {
    Ok(BaoPiecewisePolicy {
        struct_size: std::mem::size_of::<BaoPiecewisePolicy>() as u32,
        abi_version: ABI_VERSION,
        arithmetic: supernova_tag_id("arithmetic", &p.arithmetic)
            .ok_or("UNSUPPORTED_ARITHMETIC")?,
        include_predictions: p.include_predictions as u32,
        include_residuals: p.include_residuals as u32,
        reserved: 0,
        maximum_models: p.maximum_models,
        maximum_rows: p.maximum_rows,
        maximum_matrix_elements: p.maximum_matrix_elements,
        maximum_string_bytes: p.maximum_string_bytes,
        maximum_native_bytes: p.maximum_native_bytes,
        maximum_array_elements: p.maximum_array_elements,
        maximum_native_output_bytes: p.maximum_native_output_bytes,
        maximum_queries: p.maximum_queries,
        maximum_total_segment_visits: p.maximum_total_segment_visits,
        maximum_forward_sensitivity: p.maximum_forward_sensitivity,
    })
}
use super::bao::{buffer, bytes, doubles, status, strings, text};
struct Prepared(*mut c_void);
impl Drop for Prepared {
    fn drop(&mut self) {
        unsafe {
            cosmo_bao_piecewise_destroy(self.0);
        }
    }
}
struct ResultOwner(*mut c_void);
impl Drop for ResultOwner {
    fn drop(&mut self) {
        unsafe {
            cosmo_bao_piecewise_result_destroy(self.0);
        }
    }
}
pub(crate) fn evaluate(
    source: &Decoded,
    models: &[Model],
    preparation: &Policy,
    evaluation: &Policy,
) -> Result<Value, String> {
    let count = source.rows.len();
    if count > 4096 || models.len() > 64 || count as u64 > preparation.maximum_rows {
        return Err("RESOURCE_LIMIT".into());
    }
    let pp = policy(preparation)?;
    let ep = policy(evaluation)?;
    let queries = source
        .rows
        .iter()
        .map(|r| {
            Ok(BaoQuery {
                z: r.z,
                observable: bao_tag_id("observable", &r.observable)
                    .ok_or("UNSUPPORTED_OBSERVABLE")?,
                reserved: 0,
            })
        })
        .collect::<Result<Vec<_>, String>>()?;
    let observed = source.rows.iter().map(|r| r.value).collect::<Vec<_>>();
    let ids = source.rows.iter().map(|r| bytes(&r.id)).collect::<Vec<_>>();
    let axes = source
        .covariance_axis_ids
        .iter()
        .map(|id| bytes(id))
        .collect::<Vec<_>>();
    let d = BaoDescriptor {
        struct_size: std::mem::size_of::<BaoDescriptor>() as u32,
        abi_version: ABI_VERSION,
        role: bao_tag_id("role", source.role).ok_or("UNSUPPORTED_ROLE")?,
        covariance_unit: bao_tag_id("covariance_unit", "ratio_squared").unwrap(),
        queries: queries.as_ptr(),
        query_count: queries.len() as u64,
        query_byte_length: std::mem::size_of_val(queries.as_slice()) as u64,
        observed: buffer(&observed),
        covariance: buffer(&source.covariance),
        ordered_ids: strings(&ids),
        covariance_axis_ids: strings(&axes),
        table_identity: bytes(source.table_asset.sha256()),
        covariance_identity: bytes(source.covariance_asset.sha256()),
        ordering_provenance: bytes(&source.ordering_provenance),
        calibration_provenance: bytes(&source.calibration_provenance),
        dependence_provenance: bytes(&source.dependence_provenance),
        redshift_convention: bytes(&source.redshift_convention),
        ruler_convention: bytes(&source.ruler_convention),
        computational_h0_convention: bytes(&source.computational_h0_convention),
    };
    let mut raw = ptr::null_mut();
    status(unsafe { cosmo_bao_piecewise_prepare(&d, &pp, &mut raw) })?;
    if raw.is_null() {
        return Err("NULL_PREPARED".into());
    }
    let prepared = Prepared(raw);
    let mut v: BaoPiecewiseSourceView = unsafe { std::mem::zeroed() };
    status(unsafe { cosmo_bao_piecewise_source_view(prepared.0, &mut v) })?;
    if v.struct_size != std::mem::size_of::<BaoPiecewiseSourceView>() as u32
        || v.abi_version != ABI_VERSION
    {
        return Err("INVALID_NATIVE_SOURCE_VIEW".into());
    }
    let attempted = models
        .iter()
        .map(|m| BaoPiecewiseModel {
            q: m.q,
            h0_rd_km_s: m.h0_rd_km_s,
        })
        .collect::<Vec<_>>();
    let batch = BaoPiecewiseBatch {
        struct_size: std::mem::size_of::<BaoPiecewiseBatch>() as u32,
        abi_version: ABI_VERSION,
        models: attempted.as_ptr(),
        model_count: attempted.len() as u64,
        model_byte_length: std::mem::size_of_val(attempted.as_slice()) as u64,
    };
    let mut raw = ptr::null_mut();
    status(unsafe { cosmo_bao_piecewise_evaluate(prepared.0, &batch, &ep, &mut raw) })?;
    if raw.is_null() {
        return Err("NULL_RESULT".into());
    }
    let result = ResultOwner(raw);
    let mut p = ptr::null();
    let mut length = 0;
    let mut batch_status = 0;
    let mut numerical = 0;
    let mut segments = 0;
    status(unsafe {
        cosmo_bao_piecewise_result_view(
            result.0,
            &mut p,
            &mut length,
            &mut batch_status,
            &mut numerical,
            &mut segments,
        )
    })?;
    if length > models.len() as u64
        || (length > 0
            && (p.is_null() || (p as usize) % std::mem::align_of::<BaoPiecewiseRow>() != 0))
    {
        return Err("INVALID_NATIVE_ROWS".into());
    }
    let rows = if length == 0 {
        &[][..]
    } else {
        unsafe { std::slice::from_raw_parts(p, length as usize) }
    };
    let mut output = Vec::with_capacity(rows.len());
    let mut all_finite = batch_status == 0 && rows.len() == models.len();
    for (i, row) in rows.iter().enumerate() {
        let src = &row.source_parameters;
        let expected = &attempted[i];
        if row.struct_size != std::mem::size_of::<BaoPiecewiseRow>() as u32
            || row.abi_version != ABI_VERSION
            || row.reserved != 0
            || row.model_index != i as u64
            || src
                .q
                .iter()
                .zip(expected.q.iter())
                .any(|(a, b)| a.to_bits() != b.to_bits())
            || src.h0_rd_km_s.to_bits() != expected.h0_rd_km_s.to_bits()
        {
            return Err("NATIVE_SOURCE_MISMATCH".into());
        }
        let finite = row.status == 0;
        all_finite &= finite;
        let values = [
            row.log_density,
            row.quadratic,
            row.log_determinant,
            row.normalization,
            row.backward_residual,
            row.estimated_forward_sensitivity,
        ];
        if finite && values.iter().any(|v| !v.is_finite()) {
            return Err("NONFINITE_NATIVE_PAYLOAD".into());
        }
        let predictions = doubles(&row.predictions, count)?;
        let residuals = doubles(&row.residuals, count)?;
        if (!finite && (!predictions.is_empty() || !residuals.is_empty()))
            || (finite
                && ((evaluation.include_predictions && predictions.len() != count)
                    || (evaluation.include_residuals && residuals.len() != count)))
        {
            return Err("INVALID_NATIVE_ARRAY_SHAPE".into());
        }
        output.push(json!({"identity":{"model_index":i,"source_parameters":models[i],"model_id":text(&row.model_id)?,"arithmetic_id":text(&row.arithmetic_id)?,"equation_id":text(&row.equation_id)?,"ruler_convention":text(&row.ruler_convention_id)?,"computational_h0_convention":text(&row.computational_h0_convention_id)?},"status":gaussian_status_name(row.status).ok_or("UNKNOWN_NATIVE_STATUS")?,"background_status":background_status_name(row.background_status).ok_or("UNKNOWN_NATIVE_CAUSE")?,"numerical_status":numerical_status_name(row.numerical_status).ok_or("UNKNOWN_NATIVE_CAUSE")?,"segment_visits":row.segment_visits,"result":if finite {json!({"kind":"finite","log_density":row.log_density,"quadratic":row.quadratic,"log_determinant":row.log_determinant,"normalization":row.normalization,"backward_residual":row.backward_residual,"estimated_forward_sensitivity":row.estimated_forward_sensitivity,"predictions":predictions,"residuals":residuals})}else{json!({"kind":"failure"})}}));
    }
    Ok(
        json!({"kind":"bao_piecewise_normalized_gaussian_batch","execution_status":"completed","numerical_status":if all_finite{"checks_passed"}else{"failed"},"interpretation_status":"unqualified","preparation":{"status":gaussian_status_name(v.status).ok_or("UNKNOWN_NATIVE_STATUS")?,"numerical_status":numerical_status_name(v.numerical_status).ok_or("UNKNOWN_NATIVE_CAUSE")?,"arithmetic_id":text(&v.arithmetic_id)?,"equation_id":text(&v.equation_id)?,"constants_id":text(&v.constants_id)?,"policy":preparation},"batch_status":gaussian_status_name(batch_status).ok_or("UNKNOWN_NATIVE_STATUS")?,"batch_numerical_status":numerical_status_name(numerical).ok_or("UNKNOWN_NATIVE_CAUSE")?,"evaluation_policy":evaluation,"segment_visits":segments,"rows":output}),
    )
}
