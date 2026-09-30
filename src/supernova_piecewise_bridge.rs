//! Structural retained analytic profile transport; C++ owns all equations.
use super::generated::*;
type WireModel = SupernovaPiecewiseModel;
type WireBatch = SupernovaPiecewiseBatch;
type WireSlot = SupernovaPiecewiseSlot;
use super::generated::{
    cosmo_supernova_piecewise_evaluate as evaluate_ffi,
    cosmo_supernova_piecewise_result_destroy as destroy_ffi,
    cosmo_supernova_piecewise_result_view as view_ffi,
};

use super::observations::{ObservationInput, copied, prepare_retained};
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};
use std::{ffi::c_void, ptr};

#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct Model {
    pub q: [f64; 5],
}
#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct Policy {
    pub arithmetic: String,
    pub include_residual_arrays: bool,
    pub maximum_models: u64,
    pub maximum_source_rows: u64,
    pub maximum_matrix_elements: u64,
    pub maximum_queries: u64,
    pub maximum_array_elements: u64,
    pub maximum_native_output_bytes: u64,
    pub maximum_total_segment_visits: u64,
    pub maximum_forward_sensitivity: f64,
}
struct Consumer(*mut c_void);
impl Drop for Consumer {
    fn drop(&mut self) {
        unsafe {
            cosmo_supernova_destroy(self.0);
        }
    }
}
struct ResultOwner(*mut c_void);
impl Drop for ResultOwner {
    fn drop(&mut self) {
        unsafe {
            destroy_ffi(self.0);
        }
    }
}
fn text(v: &Bytes) -> Result<String, String> {
    let data = copied(v.data, v.length, 1_048_576)?;
    String::from_utf8(data).map_err(|_| "INVALID_SOURCE_ENCODING".into())
}
fn doubles(v: &F64Buffer, cap: usize) -> Result<Vec<f64>, String> {
    // Optional exports use an entirely zero descriptor to denote absence.
    if v.struct_size == 0
        && v.abi_version == 0
        && v.element_type == 0
        && v.reserved == 0
        && v.data.is_null()
        && v.length == 0
        && v.byte_length == 0
    {
        return Ok(Vec::new());
    }
    if v.struct_size != std::mem::size_of::<F64Buffer>() as u32
        || v.abi_version != ABI_VERSION
        || v.element_type != 2
        || v.reserved != 0
        || v.length.checked_mul(8) != Some(v.byte_length)
    {
        return Err("INVALID_BUFFER_VIEW".into());
    }
    let values = copied(v.data, v.length, cap)?;
    if values.iter().any(|x| !x.is_finite()) {
        return Err("INVALID_FINITE_PAYLOAD".into());
    }
    Ok(values)
}
fn ids(v: &Strings, cap: usize) -> Result<Vec<String>, String> {
    let n = usize::try_from(v.length).map_err(|_| "INVALID_SOURCE_LENGTH")?;
    if n > cap
        || v.length.checked_mul(std::mem::size_of::<Bytes>() as u64) != Some(v.byte_length)
        || (n > 0 && (v.data.is_null() || (v.data as usize) % std::mem::align_of::<Bytes>() != 0))
    {
        return Err("INVALID_SOURCE_VIEW".into());
    }
    let entries = if n == 0 {
        &[][..]
    } else {
        unsafe { std::slice::from_raw_parts(v.data, n) }
    };
    entries.iter().map(text).collect()
}
pub(crate) fn evaluate(
    source: &ObservationInput,
    selection: &str,
    observation_policy: [u64; 3],
    models: &[Model],
    p: &Policy,
) -> Result<Value, String> {
    if selection != "pantheon_zhd_gt_001" || source.metadata.profile != "pantheon_plus_released_v1"
    {
        return Err("INCOMPATIBLE_SUPERNOVA_SOURCE_SELECTION".into());
    }
    if models.len() > 64
        || source.values.len() > 4096
        || source.values.len() as u64 > p.maximum_source_rows
        || source.uncertainty_matrix.len() as u64 > p.maximum_matrix_elements
    {
        return Err("RESOURCE_LIMIT".into());
    }
    let (observations, _) = prepare_retained(source, selection, observation_policy)?;
    let policy = SupernovaPiecewisePolicy {
        struct_size: std::mem::size_of::<SupernovaPiecewisePolicy>() as u32,
        abi_version: ABI_VERSION,
        arithmetic: supernova_tag_id("arithmetic", &p.arithmetic)
            .ok_or("UNKNOWN_ARITHMETIC_POLICY")?,
        include_residual_arrays: u32::from(p.include_residual_arrays),
        maximum_models: p.maximum_models,
        maximum_source_rows: p.maximum_source_rows,
        maximum_matrix_elements: p.maximum_matrix_elements,
        maximum_array_elements: p.maximum_array_elements,
        maximum_native_output_bytes: p.maximum_native_output_bytes,
        maximum_queries: p.maximum_queries,
        maximum_total_segment_visits: p.maximum_total_segment_visits,
        maximum_forward_sensitivity: p.maximum_forward_sensitivity,
    };
    let mut raw = ptr::null_mut();
    let status = unsafe { cosmo_supernova_piecewise_prepare(observations.0, &policy, &mut raw) };
    let consumer = Consumer(raw);
    if status != OK {
        return Err(format!("CORE_STATUS_{status}"));
    }
    if raw.is_null() {
        return Err("NULL_SUPERNOVA_OWNER".into());
    }
    // Every field is an integer, float, pointer or C descriptor; zero is valid initialization.
    let mut view: SupernovaView = unsafe { std::mem::zeroed() };
    view.struct_size = std::mem::size_of::<SupernovaView>() as u32;
    view.abi_version = ABI_VERSION;
    if unsafe { cosmo_supernova_source_view(consumer.0, &mut view) } != OK {
        return Err("INVALID_SOURCE_VIEW".into());
    }
    let ordered_ids = ids(&view.ordered_ids, 4096)?;
    if view.selected_source_indices.length.checked_mul(8)
        != Some(view.selected_source_indices.byte_length)
    {
        return Err("INVALID_SOURCE_VIEW".into());
    }
    let indices = copied(
        view.selected_source_indices.data,
        view.selected_source_indices.length,
        4096,
    )?;
    let expansion = doubles(&view.selected_z_expansion, 4096)?;
    let observer = doubles(&view.selected_z_observer, 4096)?;
    if view.struct_size != std::mem::size_of::<SupernovaView>() as u32
        || view.abi_version != ABI_VERSION
        || view.arithmetic != policy.arithmetic
        || view.source_row_count != source.values.len() as u64
        || indices.len() != ordered_ids.len()
        || expansion.len() != ordered_ids.len()
        || observer.len() != ordered_ids.len()
        || view.matrix_validation_assessed > 1
        || view.matrix_validation_scope > 2
    {
        return Err("INVALID_SOURCE_VIEW".into());
    }
    for (i, index) in indices.iter().enumerate() {
        let index = usize::try_from(*index).map_err(|_| "INVALID_SOURCE_INDEX")?;
        if index >= source.values.len()
            || ordered_ids[i] != source.measurement_ids[index]
            || expansion[i].to_bits() != source.zhd[index].to_bits()
            || observer[i].to_bits() != source.zhel[index].to_bits()
        {
            return Err("INVALID_SOURCE_IDENTITY".into());
        }
    }
    let matrix_scope = match view.matrix_validation_scope {
        GAUSSIAN_MATRIX_VALIDATION_SCOPE_FULL_DECLARED_MATRIX => "full_declared_matrix",
        GAUSSIAN_MATRIX_VALIDATION_SCOPE_SELECTED_COVARIANCE_ONLY => "selected_covariance_only",
        GAUSSIAN_MATRIX_VALIDATION_SCOPE_FULL_PRECISION_THEN_MARGINAL => {
            "full_precision_then_marginal"
        }
        _ => return Err("INVALID_MATRIX_SCOPE".into()),
    };
    let metadata = json!({"status":supernova_status_name(view.status).ok_or("INVALID_SCIENTIFIC_STATUS")?,
                "preparation_status":gaussian_status_name(view.preparation_status).ok_or("INVALID_SCIENTIFIC_STATUS")?,
                "preparation_numerical_status":numerical_status_name(view.preparation_numerical_status).ok_or("INVALID_NUMERICAL_STATUS")?,
                "arithmetic_id":text(&view.arithmetic_id)?, "score_id":text(&view.score_id)?, "constant_set_id":text(&view.constant_set_id)?,
                "radial_equation_id":text(&view.radial_equation_id)?, "shape_convention":text(&view.shape_convention)?,
                "offset_convention":text(&view.offset_convention)?, "query_provenance":text(&view.query_provenance)?,
                "table_identity":text(&view.table_identity)?, "uncertainty_identity":text(&view.uncertainty_identity)?,
                "ordering_provenance":text(&view.ordering_provenance)?, "calibration_provenance":text(&view.calibration_provenance)?,
                "dependence_provenance":text(&view.dependence_provenance)?, "source_semantics":text(&view.source_semantics)?,
                "matrix_validation_scope":if view.matrix_validation_assessed == 1 {Some(matrix_scope)}else{None}, "matrix_validation_assessed":view.matrix_validation_assessed == 1,
                "ordered_ids":ordered_ids, "selected_source_indices":indices, "source_row_count":view.source_row_count,
                "z_expansion":expansion, "z_observer":observer});
    let parameters: Vec<_> = models.iter().map(|m| WireModel { q: m.q }).collect();
    let descriptor = WireBatch {
        struct_size: std::mem::size_of::<WireBatch>() as u32,
        abi_version: ABI_VERSION,
        models: parameters.as_ptr(),
        model_count: parameters.len() as u64,
        model_bytes: std::mem::size_of_val(parameters.as_slice()) as u64,
    };
    let mut raw = ptr::null_mut();
    let status = unsafe { evaluate_ffi(consumer.0, &descriptor, &policy, &mut raw) };
    let result = ResultOwner(raw);
    if status != OK {
        return Err(format!("CORE_STATUS_{status}"));
    }
    if raw.is_null() {
        return Err("NULL_SUPERNOVA_RESULT".into());
    }
    let mut slots = ptr::null();
    let mut count = 0;
    let mut batch_status = 0;
    let mut native_segments = 0;
    if unsafe {
        view_ffi(
            result.0,
            &mut slots,
            &mut count,
            &mut batch_status,
            &mut native_segments,
        )
    } != OK
        || count > models.len() as u64
        || (count > 0
            && (slots.is_null() || (slots as usize) % std::mem::align_of::<WireSlot>() != 0))
    {
        return Err("INVALID_RESULT_VIEW".into());
    }
    if view.status == 0 && batch_status == 0 && count != models.len() as u64 {
        return Err("INVALID_RESULT_LENGTH".into());
    }
    let slots = if count == 0 {
        &[][..]
    } else {
        unsafe { std::slice::from_raw_parts(slots, count as usize) }
    };
    let mut rows = Vec::with_capacity(slots.len());
    let mut evaluations = 0u64;
    let mut failed = view.status != 0 || batch_status != 0;
    for (i, s) in slots.iter().enumerate() {
        let m = &parameters[i];
        if s.struct_size as usize != std::mem::size_of::<WireSlot>()
            || s.abi_version != ABI_VERSION
            || s.model_index != i as u64
            || s.source_parameters
                .q
                .iter()
                .zip(m.q)
                .any(|(a, b)| a.to_bits() != b.to_bits())
            || s.reserved != 0
            || s.has_profile_payload > 1
        {
            return Err("INVALID_RESULT_ORDER".into());
        }
        let result_ids = ids(&s.ordered_ids, 4096)?;
        if s.selected_count.checked_mul(8) != Some(s.selected_index_bytes) {
            return Err("INVALID_RESULT_SOURCE".into());
        }
        let result_indices = copied(s.selected_source_indices, s.selected_count, 4096)?;
        if result_ids != ordered_ids
            || result_indices != indices
            || doubles(&s.expansion_z, 4096)?
                .iter()
                .map(|x| x.to_bits())
                .collect::<Vec<_>>()
                != expansion.iter().map(|x| x.to_bits()).collect::<Vec<_>>()
            || doubles(&s.observer_z, 4096)?
                .iter()
                .map(|x| x.to_bits())
                .collect::<Vec<_>>()
                != observer.iter().map(|x| x.to_bits()).collect::<Vec<_>>()
        {
            return Err("INVALID_RESULT_SOURCE".into());
        }
        evaluations = evaluations
            .checked_add(s.segment_visits)
            .ok_or("INVALID_WORK_ACCOUNTING")?;
        let mut identity = json!({"model_index":s.model_index,"model_id":text(&s.model_id)?,"arithmetic_id":text(&s.arithmetic_id)?,"score_id":text(&s.score_id)?,"radial_equation_id":text(&s.radial_equation_id)?});
        identity["source_parameters"] =
            serde_json::to_value(&models[i]).map_err(|e| e.to_string())?;
        let cause = numerical_status_name(s.numerical_status).ok_or("INVALID_NUMERICAL_STATUS")?;
        let diagnostics = json!({"numerical_status":cause,"background_status":background_status_name(s.background_status).ok_or("INVALID_SCIENTIFIC_STATUS")?,
                    "profile_status":gaussian_status_name(s.profile_status).ok_or("INVALID_SCIENTIFIC_STATUS")?,"segment_visits":s.segment_visits});
        if s.status != 0 {
            failed = true;
            if s.has_profile_payload != 0
                || s.shape_magnitudes.length != 0
                || s.base_residuals.length != 0
                || s.profiled_residuals.length != 0
            {
                return Err("FAILED_FINITE_PAYLOAD".into());
            }
            rows.push(json!({"identity":identity,"result":{"kind":"failure","status":supernova_status_name(s.status).ok_or("INVALID_SCIENTIFIC_STATUS")?},"diagnostics":diagnostics}));
            continue;
        }
        let numbers = [
            s.offset_coefficient,
            s.quadratic,
            s.relative_profile_score,
            s.backward_residual,
            s.estimated_forward_sensitivity,
            s.coefficient_solve_backward_residual,
            s.coefficient_solve_forward_sensitivity,
            s.residual_l1,
            s.solution_norm_inf,
            s.adjusted_residual_l1,
            s.adjusted_solution_norm_inf,
        ];
        if s.has_profile_payload != 1
            || s.numerical_status != 0
            || numbers.iter().any(|x| !x.is_finite())
        {
            return Err("INVALID_FINITE_PAYLOAD".into());
        }
        let arrays = json!({"shape_magnitudes":doubles(&s.shape_magnitudes,4096)?,"base_residuals":doubles(&s.base_residuals,4096)?,"profiled_residuals":doubles(&s.profiled_residuals,4096)?});
        let expected = if p.include_residual_arrays {
            ordered_ids.len() as u64
        } else {
            0
        };
        if s.shape_magnitudes.length != expected
            || s.base_residuals.length != expected
            || s.profiled_residuals.length != expected
        {
            return Err("INVALID_RESULT_LENGTH".into());
        }
        rows.push(json!({"identity":identity,"result":{"kind":"finite","offset_coefficient":s.offset_coefficient,"quadratic":s.quadratic,
                    "relative_profile_score":s.relative_profile_score,"density_status":"not_applicable","normalization_status":"not_applicable","arrays":arrays},
                    "diagnostics":{"numerical_status":cause,"segment_visits":s.segment_visits,"backward_residual":s.backward_residual,
                    "estimated_forward_sensitivity":s.estimated_forward_sensitivity,"coefficient_solve_backward_residual":s.coefficient_solve_backward_residual,
                    "coefficient_solve_forward_sensitivity":s.coefficient_solve_forward_sensitivity,"residual_l1":s.residual_l1,"solution_norm_inf":s.solution_norm_inf,
                    "adjusted_residual_l1":s.adjusted_residual_l1,"adjusted_solution_norm_inf":s.adjusted_solution_norm_inf}}));
    }
    if evaluations != native_segments || evaluations > p.maximum_total_segment_visits {
        return Err("INVALID_GLOBAL_WORK_ACCOUNTING".into());
    }
    Ok(
        json!({"kind":if failed {"failure"}else{"finite"},"error_id":if failed {Some("SUPERNOVA_EVALUATION_FAILURE")}else{None},
                "execution_status":"completed","numerical_status":if failed {"failed"}else{"checks_passed"},"interpretation_status":"unqualified",
                "validation_coverage":{"profile":"original_inputs_piecewise_six_points_v1","applicability":"named regression coverage only; this request is not automatically qualified"},
                "source":metadata,"batch_status":supernova_status_name(batch_status).ok_or("INVALID_BATCH_STATUS")?,"segment_visits":evaluations,"rows":rows}),
    )
}
