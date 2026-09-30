//! Coarse retained native BAO transport. All scientific equations stay in C++.
use super::generated::*;
use crate::bao_ingestion::Decoded;
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};
use std::{ffi::c_void, ptr};
#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct Model {
    pub model: String,
    pub omega_m: f64,
    pub constant_q: f64,
    pub w0: f64,
    pub wa: f64,
    pub h0_rd_km_s: f64,
}
#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct Policy {
    pub background: super::background_v2::Policy,
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
    pub maximum_forward_sensitivity: f64,
}
fn policy(p: &Policy) -> Result<BaoPolicy, String> {
    let q = &p.background;
    Ok(BaoPolicy {
        struct_size: std::mem::size_of::<BaoPolicy>() as u32,
        abi_version: ABI_VERSION,
        background: BackgroundPolicy {
            struct_size: std::mem::size_of::<BackgroundPolicy>() as u32,
            abi_version: ABI_VERSION,
            maximum_depth: q.maximum_depth,
            reserved: 0,
            maximum_parameters: q.maximum_parameters,
            maximum_queries: q.maximum_queries,
            maximum_slots: q.maximum_slots,
            maximum_native_output_bytes: q.maximum_native_output_bytes,
            maximum_total_evaluations: q.maximum_total_evaluations,
            maximum_evaluations_per_integral: q.maximum_evaluations_per_integral,
            absolute_tolerance: q.absolute_tolerance,
            relative_tolerance: q.relative_tolerance,
        },
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
        maximum_forward_sensitivity: p.maximum_forward_sensitivity,
    })
}
pub(super) fn bytes(s: &str) -> Bytes {
    Bytes {
        data: s.as_ptr(),
        length: s.len() as u64,
    }
}
pub(super) fn buffer(s: &[f64]) -> F64Buffer {
    F64Buffer {
        struct_size: std::mem::size_of::<F64Buffer>() as u32,
        abi_version: ABI_VERSION,
        element_type: 2,
        reserved: 0,
        data: s.as_ptr(),
        length: s.len() as u64,
        byte_length: std::mem::size_of_val(s) as u64,
    }
}
pub(super) fn strings(s: &[Bytes]) -> Strings {
    Strings {
        data: s.as_ptr(),
        length: s.len() as u64,
        byte_length: std::mem::size_of_val(s) as u64,
    }
}
pub(super) fn text(s: &Bytes) -> Result<String, String> {
    if s.length > 16777216 || (s.length > 0 && s.data.is_null()) {
        return Err("INVALID_NATIVE_TEXT".into());
    }
    if s.length == 0 {
        return Ok(String::new());
    }
    std::str::from_utf8(unsafe { std::slice::from_raw_parts(s.data, s.length as usize) })
        .map(str::to_owned)
        .map_err(|_| "INVALID_NATIVE_UTF8".into())
}
pub(super) fn doubles(s: &F64Buffer, cap: usize) -> Result<Vec<f64>, String> {
    if s.struct_size == 0
        && s.abi_version == 0
        && s.element_type == 0
        && s.reserved == 0
        && s.data.is_null()
        && s.length == 0
        && s.byte_length == 0
    {
        return Ok(vec![]);
    }
    if s.struct_size != std::mem::size_of::<F64Buffer>() as u32
        || s.abi_version != ABI_VERSION
        || s.element_type != 2
        || s.reserved != 0
        || s.length > cap as u64
        || s.length.checked_mul(8) != Some(s.byte_length)
        || (s.length > 0 && (s.data.is_null() || (s.data as usize) % 8 != 0))
    {
        return Err("INVALID_NATIVE_ARRAY".into());
    }
    Ok(if s.length == 0 {
        vec![]
    } else {
        unsafe { std::slice::from_raw_parts(s.data, s.length as usize) }.to_vec()
    })
}
struct Prepared(*mut c_void);
impl Drop for Prepared {
    fn drop(&mut self) {
        unsafe {
            cosmo_bao_destroy(self.0);
        }
    }
}
struct ResultOwner(*mut c_void);
impl Drop for ResultOwner {
    fn drop(&mut self) {
        unsafe {
            cosmo_bao_result_destroy(self.0);
        }
    }
}
pub(super) fn status(s: u32) -> Result<(), String> {
    if s == OK {
        Ok(())
    } else {
        Err(format!("CORE_STATUS_{s}"))
    }
}
pub(crate) fn evaluate(
    source: &Decoded,
    models: &[Model],
    preparation: &Policy,
    evaluation: &Policy,
) -> Result<Value, String> {
    let count = source.rows.len();
    if count > 4096
        || models.len() > 64
        || count as u64 > preparation.maximum_rows
        || models.len() as u64 > evaluation.maximum_models
        || count
            .checked_mul(models.len())
            .and_then(|x| x.checked_mul(2))
            .ok_or("RESOURCE_LIMIT")? as u64
            > evaluation.maximum_array_elements
    {
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
    status(unsafe { cosmo_bao_prepare(&d, &pp, &mut raw) })?;
    if raw.is_null() {
        return Err("NULL_PREPARED".into());
    }
    let prepared = Prepared(raw);
    let mut v: BaoSourceView = unsafe { std::mem::zeroed() };
    status(unsafe { cosmo_bao_source_view(prepared.0, &mut v) })?;
    if v.struct_size != std::mem::size_of::<BaoSourceView>() as u32 || v.abi_version != ABI_VERSION
    {
        return Err("INVALID_NATIVE_SOURCE_VIEW".into());
    }
    let attempted = models
        .iter()
        .map(|m| {
            Ok(BaoModel {
                parameters: SupernovaModelV2 {
                    model: background_v2_tag_id("model", &m.model).ok_or("UNSUPPORTED_MODEL")?,
                    reserved: 0,
                    omega_m: m.omega_m,
                    constant_q: m.constant_q,
                    w0: m.w0,
                    wa: m.wa,
                },
                h0_rd_km_s: m.h0_rd_km_s,
            })
        })
        .collect::<Result<Vec<_>, String>>()?;
    let batch = BaoBatch {
        struct_size: std::mem::size_of::<BaoBatch>() as u32,
        abi_version: ABI_VERSION,
        models: attempted.as_ptr(),
        model_count: attempted.len() as u64,
        model_byte_length: std::mem::size_of_val(attempted.as_slice()) as u64,
    };
    let mut raw = ptr::null_mut();
    status(unsafe { cosmo_bao_evaluate(prepared.0, &batch, &ep, &mut raw) })?;
    if raw.is_null() {
        return Err("NULL_RESULT".into());
    }
    let result = ResultOwner(raw);
    let mut p = ptr::null();
    let mut length = 0;
    let mut batch_status = 0;
    let mut numerical = 0;
    let mut callbacks = 0;
    status(unsafe {
        cosmo_bao_result_view(
            result.0,
            &mut p,
            &mut length,
            &mut batch_status,
            &mut numerical,
            &mut callbacks,
        )
    })?;
    if length > models.len() as u64
        || (length > 0 && (p.is_null() || (p as usize) % std::mem::align_of::<BaoRow>() != 0))
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
        if row.model_index != i as u64
            || src.parameters.model != expected.parameters.model
            || src.parameters.reserved != 0
            || [
                src.parameters.omega_m,
                src.parameters.constant_q,
                src.parameters.w0,
                src.parameters.wa,
                src.h0_rd_km_s,
            ]
            .into_iter()
            .zip([
                expected.parameters.omega_m,
                expected.parameters.constant_q,
                expected.parameters.w0,
                expected.parameters.wa,
                expected.h0_rd_km_s,
            ])
            .any(|(a, b)| a.to_bits() != b.to_bits())
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
        output.push(json!({"identity":{"model_index":i,"source_parameters":models[i],"model_id":text(&row.model_id)?,"arithmetic_id":text(&row.arithmetic_id)?,"equation_id":text(&row.equation_id)?,"ruler_convention":text(&row.ruler_convention_id)?,"computational_h0_convention":text(&row.computational_h0_convention_id)?},"status":gaussian_status_name(row.status).ok_or("UNKNOWN_NATIVE_STATUS")?,"background_status":background_status_name(row.background_status).ok_or("UNKNOWN_NATIVE_CAUSE")?,"numerical_status":numerical_status_name(row.numerical_status).ok_or("UNKNOWN_NATIVE_CAUSE")?,"evaluations":row.evaluations,"result":if finite {json!({"kind":"finite","log_density":row.log_density,"quadratic":row.quadratic,"log_determinant":row.log_determinant,"normalization":row.normalization,"backward_residual":row.backward_residual,"estimated_forward_sensitivity":row.estimated_forward_sensitivity,"predictions":predictions,"residuals":residuals})}else{json!({"kind":"failure"})}}));
    }
    Ok(
        json!({"kind":"bao_normalized_gaussian_batch","execution_status":"completed","numerical_status":if all_finite{"checks_passed"}else{"failed"},"interpretation_status":"unqualified","preparation":{"status":gaussian_status_name(v.status).ok_or("UNKNOWN_NATIVE_STATUS")?,"numerical_status":numerical_status_name(v.numerical_status).ok_or("UNKNOWN_NATIVE_CAUSE")?,"arithmetic_id":text(&v.arithmetic_id)?,"equation_id":text(&v.equation_id)?,"constants_id":text(&v.constants_id)?,"policy":preparation},"batch_status":gaussian_status_name(batch_status).ok_or("UNKNOWN_NATIVE_STATUS")?,"batch_numerical_status":numerical_status_name(numerical).ok_or("UNKNOWN_NATIVE_CAUSE")?,"evaluation_policy":evaluation,"evaluations":callbacks,"rows":output}),
    )
}
