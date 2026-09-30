//! Coarse retained Gaussian calculation; all probability algebra lives in C++.
use super::{
    generated::*,
    observations::{Prepared, bytes, copied, doubles, strings},
};
use serde_json::{Value, json};
use std::{ffi::c_void, ptr};
struct Gaussian(*mut c_void);
impl Drop for Gaussian {
    fn drop(&mut self) {
        unsafe {
            cosmo_gaussian_destroy(self.0);
        }
    }
}
struct Batch(*mut c_void);
impl Drop for Batch {
    fn drop(&mut self) {
        unsafe {
            cosmo_gaussian_result_destroy(self.0);
        }
    }
}
pub(crate) struct ProperPrior<'a> {
    pub response: &'a [f64],
    pub ordered_ids: &'a [String],
    pub mean: f64,
    pub variance: f64,
    pub latent_identity: &'a str,
    pub independence_declared: bool,
}
fn semantic(transport: u32, status: u32, raw: *mut c_void) -> Result<Gaussian, String> {
    let owner = Gaussian(raw);
    if transport != OK {
        return Err(format!("CORE_STATUS_{transport}"));
    }
    if status != GAUSSIAN_STATUS_FINITE {
        return Err(format!(
            "GAUSSIAN_{}",
            gaussian_status_name(status).ok_or("INVALID_SCIENTIFIC_STATUS")?
        ));
    }
    if raw.is_null() {
        return Err("NULL_GAUSSIAN".into());
    }
    Ok(owner)
}
pub(crate) fn gaussian_batch(
    source: &Prepared,
    selection: &str,
    residual_unit: &str,
    observation_policy: [u64; 3],
    residuals: &[Vec<f64>],
    ordered_ids: &[String],
    response: Option<&[f64]>,
    prior: Option<ProperPrior<'_>>,
    maximum_batch_elements: u64,
    maximum_forward_sensitivity: f64,
    maximum_preparation_bytes: u64,
    maximum_evaluation_bytes: u64,
) -> Result<Value, String> {
    if response.is_some() && prior.is_some() {
        return Err("PROFILE_AND_PROPER_PRIOR_MUTUALLY_EXCLUSIVE".into());
    }
    let n = ordered_ids.len();
    let total = residuals.len().checked_mul(n).ok_or("RESOURCE_LIMIT")?;
    if total as u64 > maximum_batch_elements
        || maximum_batch_elements > MAX_BATCH_ELEMENTS
        || residuals.iter().any(|r| r.len() != n)
    {
        return Err("INVALID_BATCH_SHAPE".into());
    }
    let source_rows = source.source_view()?.values()?.len();
    let selected = selected_indices(source, selection, source_rows)?;
    let mut policy = GaussianPolicy {
        struct_size: std::mem::size_of::<GaussianPolicy>() as u32,
        abi_version: ABI_VERSION,
        reserved: 0,
        reserved2: 0,
        maximum_matrix_elements: observation_policy[1],
        maximum_batch_elements,
        maximum_string_bytes: observation_policy[2],
        maximum_forward_sensitivity,
        maximum_native_bytes: maximum_preparation_bytes,
    };
    let mut raw = ptr::null_mut();
    let mut status = u32::MAX;
    let mut numerical_status = u32::MAX;
    let tag = observation_tag_id("selection", selection).ok_or("UNKNOWN_SELECTION")?;
    let transport = unsafe {
        cosmo_gaussian_prepare(
            source.0,
            tag,
            &policy,
            &mut raw,
            &mut status,
            &mut numerical_status,
        )
    };
    if transport == OK && status != GAUSSIAN_STATUS_FINITE {
        let _drop = Gaussian(raw);
        return Ok(
            json!({"kind":"failure","status":gaussian_status_name(status).ok_or("INVALID_SCIENTIFIC_STATUS")?,"numerical_status":numerical_status_name(numerical_status).ok_or("INVALID_NUMERICAL_STATUS")?,"phase":"preparation","evaluations":[],"selected_source_indices":selected}),
        );
    }
    let mut owner = semantic(transport, status, raw)?;
    if let Some(prior) = prior {
        let ids: Vec<_> = prior.ordered_ids.iter().map(|x| bytes(x)).collect();
        let descriptor = GaussianPrior {
            struct_size: std::mem::size_of::<GaussianPrior>() as u32,
            abi_version: ABI_VERSION,
            independence_declared: prior.independence_declared as u32,
            reserved: 0,
            mean: prior.mean,
            variance: prior.variance,
            latent_identity: bytes(prior.latent_identity),
            response: doubles(prior.response),
            ordered_ids: strings(&ids),
        };
        raw = ptr::null_mut();
        status = u32::MAX;
        let transport = unsafe {
            cosmo_gaussian_proper_offset(
                owner.0,
                &descriptor,
                &policy,
                &mut raw,
                &mut status,
                &mut numerical_status,
            )
        };
        if transport == OK && status != GAUSSIAN_STATUS_FINITE {
            let _drop = Gaussian(raw);
            return Ok(
                json!({"kind":"failure","status":gaussian_status_name(status).ok_or("INVALID_SCIENTIFIC_STATUS")?,"numerical_status":numerical_status_name(numerical_status).ok_or("INVALID_NUMERICAL_STATUS")?,"phase":"proper_prior_preparation","evaluations":[],"selected_source_indices":selected}),
            );
        }
        owner = semantic(transport, status, raw)?;
    }
    policy.maximum_native_bytes = maximum_evaluation_bytes;
    let ids: Vec<_> = ordered_ids.iter().map(|x| bytes(x)).collect();
    let flat: Vec<_> = residuals.iter().flatten().copied().collect();
    let descriptor = GaussianBatch {
        struct_size: std::mem::size_of::<GaussianBatch>() as u32,
        abi_version: ABI_VERSION,
        mode: if response.is_some() {
            GAUSSIAN_MODE_PROFILE_OFFSET_SCORE
        } else {
            GAUSSIAN_MODE_NORMALIZED_DENSITY
        },
        reserved: 0,
        row_count: residuals.len() as u64,
        residuals: doubles(&flat),
        ordered_ids: strings(&ids),
        response: doubles(response.unwrap_or(&[])),
    };
    raw = ptr::null_mut();
    let transport = unsafe { cosmo_gaussian_evaluate(owner.0, &descriptor, &policy, &mut raw) };
    let batch = Batch(raw);
    if transport != OK {
        return Err(format!("CORE_STATUS_{transport}"));
    }
    if raw.is_null() {
        return Err("NULL_GAUSSIAN_BATCH".into());
    }
    let mut rows = ptr::null();
    let mut length = 0;
    if unsafe { cosmo_gaussian_result_view(batch.0, &mut rows, &mut length) } != OK
        || length != residuals.len() as u64
    {
        return Err("INVALID_RESULT_VIEW".into());
    }
    // Generated structs are C-layout; copy only after checking the bounded owner view.
    let values = if length == 0 {
        &[][..]
    } else {
        if rows.is_null() || (rows as usize) % std::mem::align_of::<GaussianRow>() != 0 {
            return Err("INVALID_RESULT_VIEW".into());
        }
        unsafe { std::slice::from_raw_parts(rows, length as usize) }
    };
    let output: Result<Vec<_>, String> = values.iter().map(|row| {
        let label = gaussian_status_name(row.status).ok_or("INVALID_SCIENTIFIC_STATUS")?;
        let numerical_label = numerical_status_name(row.numerical_status).ok_or("INVALID_NUMERICAL_STATUS")?;
        if row.status != GAUSSIAN_STATUS_FINITE {
            return Ok(json!({"kind":if row.status==GAUSSIAN_STATUS_OUTSIDE_SUPPORT {"outside_support"} else {"failure"}, "status":label,"numerical_status":numerical_label}));
        }
        if response.is_some() {
            if row.numerical_status != 0 || [row.coefficient,row.quadratic,row.backward_residual,row.estimated_forward_sensitivity].iter().any(|x|!x.is_finite()) {return Err("INVALID_FINITE_PAYLOAD".into());}
            return Ok(json!({"kind":"finite","profile_coefficient":row.coefficient,"quadratic":row.quadratic,"numerical_status":numerical_label,"backward_residual":row.backward_residual,"estimated_forward_sensitivity":row.estimated_forward_sensitivity,"diagnostic_scope":"adjusted-residual solve; empirical sensitivity, not certified score error","density":null,"normalization":"not_applicable","meaning":"optimized offset score; not normalized density"}));
        }
        if row.numerical_status!=0 || [row.log_density,row.quadratic,row.log_determinant,row.normalization,row.backward_residual,row.estimated_forward_sensitivity].iter().any(|x|!x.is_finite()) {
            return Err("INVALID_FINITE_PAYLOAD".into());
        }
        Ok(json!({"kind":"finite","log_density":row.log_density,"quadratic":row.quadratic,"log_determinant":row.log_determinant,"normalization":row.normalization,"backward_residual":row.backward_residual,"estimated_forward_sensitivity":row.estimated_forward_sensitivity}))
    }).collect();
    let output = output?;
    let failed = output.iter().any(|row| row["kind"] == "failure");
    let mut view = std::mem::MaybeUninit::<GaussianView>::uninit();
    if unsafe { cosmo_gaussian_source_view(owner.0, view.as_mut_ptr()) } != OK {
        return Err("INVALID_SOURCE_VIEW".into());
    }
    let view = unsafe { view.assume_init() };
    if view.struct_size != std::mem::size_of::<GaussianView>() as u32
        || view.abi_version != ABI_VERSION
    {
        return Err("INVALID_SOURCE_VIEW".into());
    }
    let text = |b: &Bytes| -> Result<String, String> {
        String::from_utf8(copied(b.data, b.length, observation_policy[2] as usize)?)
            .map_err(|_| "INVALID_METADATA_ENCODING".into())
    };
    let decode_ids = |v: &Strings| -> Result<Vec<String>, String> {
        let length = usize::try_from(v.length).map_err(|_| "INVALID_RESULT_LENGTH")?;
        if length > source_rows
            || v.byte_length
                != v.length
                    .checked_mul(std::mem::size_of::<Bytes>() as u64)
                    .ok_or("INVALID_RESULT_LENGTH")?
            || (length > 0
                && (v.data.is_null() || (v.data as usize) % std::mem::align_of::<Bytes>() != 0))
        {
            return Err("INVALID_ID_VIEW".into());
        }
        let values = if length == 0 {
            &[][..]
        } else {
            unsafe { std::slice::from_raw_parts(v.data, length) }
        };
        values.iter().map(&text).collect()
    };
    let native_ids = decode_ids(&view.ordered_ids)?;
    let matrix_scope = match view.matrix_validation_scope {
        GAUSSIAN_MATRIX_VALIDATION_SCOPE_FULL_DECLARED_MATRIX => "full_declared_matrix",
        GAUSSIAN_MATRIX_VALIDATION_SCOPE_SELECTED_COVARIANCE_ONLY => "selected_covariance_only",
        GAUSSIAN_MATRIX_VALIDATION_SCOPE_FULL_PRECISION_THEN_MARGINAL => {
            "full_precision_then_marginal"
        }
        _ => return Err("INVALID_MATRIX_SCOPE".into()),
    };
    if view.prior_count > 1 {
        return Err("UNSUPPORTED_PRIOR_VIEW".into());
    }
    let mut priors = Vec::new();
    for index in 0..view.prior_count {
        let mut prior = std::mem::MaybeUninit::<GaussianPrior>::uninit();
        if unsafe { cosmo_gaussian_prior_view(owner.0, index, prior.as_mut_ptr()) } != OK {
            return Err("INVALID_PRIOR_VIEW".into());
        }
        let prior = unsafe { prior.assume_init() };
        if prior.struct_size != std::mem::size_of::<GaussianPrior>() as u32
            || prior.abi_version != ABI_VERSION
            || prior.reserved != 0
            || prior.independence_declared > 1
            || !prior.mean.is_finite()
            || !prior.variance.is_finite()
        {
            return Err("INVALID_PRIOR_VIEW".into());
        }
        priors.push(json!({"mean":prior.mean,"variance":prior.variance,"latent_identity":text(&prior.latent_identity)?,"response":copied(prior.response.data,prior.response.length,n)?,"applied_row_ids":decode_ids(&prior.ordered_ids)?,"independence_declared":prior.independence_declared == 1,"meaning":"original generative proper latent prior; not conditional posterior"}));
    }
    if view.selection_count > 1 {
        return Err("UNSUPPORTED_SELECTION_VIEW".into());
    }
    let mut selection_history = Vec::new();
    for index in 0..view.selection_count {
        let mut record = std::mem::MaybeUninit::<GaussianSelection>::uninit();
        if unsafe { cosmo_gaussian_selection_view(owner.0, index, record.as_mut_ptr()) } != OK {
            return Err("INVALID_SELECTION_HISTORY".into());
        }
        let record = unsafe { record.assume_init() };
        if record.struct_size != std::mem::size_of::<GaussianSelection>() as u32
            || record.abi_version != ABI_VERSION
        {
            return Err("INVALID_SELECTION_HISTORY".into());
        }
        selection_history.push(json!({"operation":text(&record.operation)?,"kept_row_ids":decode_ids(&record.kept_row_ids)?,"complement_row_ids":decode_ids(&record.complement_row_ids)?}));
    }
    let mean_shift = copied(view.mean_shift.data, view.mean_shift.length, n)?;
    Ok(
        json!({"kind":if failed {"failure"} else {"finite"},"error_id":if failed {Some("GAUSSIAN_EVALUATION_FAILURE")}else{None},"mode":if response.is_some(){"profile_offset_score"}else{"normalized_density"},"ordered_ids":native_ids,"selected_source_indices":selected,"measure":text(&view.measure)?,"input_matrix_convention":text(&view.input_matrix_convention)?,"matrix_validation_scope":matrix_scope,"source_semantics":text(&view.source_semantics)?,"table_identity":text(&view.table_identity)?,"uncertainty_identity":text(&view.uncertainty_identity)?,"ordering_provenance":text(&view.ordering_provenance)?,"calibration_provenance":text(&view.calibration_provenance)?,"dependence_provenance":text(&view.dependence_provenance)?,"prepared_distribution_treatment":text(&view.treatment)?,"residual_unit":residual_unit,"response_unit":"one","evaluation_treatment":if response.is_some(){"profile score; normalized density not applicable"}else{"normalized Gaussian density"},"mean_shift":mean_shift,"mean_shift_encoding":"empty means zero shift for each ordered row","priors":priors,"selection_history":selection_history,"rows":output}),
    )
}

fn selected_indices(
    source: &Prepared,
    selection: &str,
    source_rows: usize,
) -> Result<Vec<u64>, String> {
    let tag = observation_tag_id("selection", selection).ok_or("UNKNOWN_SELECTION")?;
    let (mut raw, mut status) = (ptr::null_mut(), 0);
    let code = unsafe { cosmo_observation_select(source.0, tag, &mut raw, &mut status) };
    let owner = super::Owned(raw);
    if code != OK {
        return Err(format!("CORE_STATUS_{code}"));
    }
    if status != OBSERVATION_STATUS_OK || raw.is_null() {
        return Err("INVALID_SELECTION_VIEW".into());
    }
    let (mut mask, mut n, mut indices, mut k) = (ptr::null(), 0, ptr::null(), 0);
    let code =
        unsafe { cosmo_result_selection_view(owner.0, &mut mask, &mut n, &mut indices, &mut k) };
    if code != OK || n != source_rows as u64 {
        return Err("INVALID_SELECTION_VIEW".into());
    }
    copied(indices, k, source_rows)
}
