//! One pooled borrowed-input call, one owned posterior result; no Rust equations.
use super::{
    generated::*,
    observations::{bytes, copied, doubles, strings},
};
use crate::gaussian_posterior_run::Request;
use serde_json::{Value, json};
use std::{mem::size_of, ptr};
struct Owned(*mut std::ffi::c_void);
impl Drop for Owned {
    fn drop(&mut self) {
        unsafe {
            irred_gaussian_posterior_result_destroy(self.0);
        }
    }
}
pub(crate) struct Calculation {
    pub output: Value,
    pub method: Value,
    pub arithmetic: Value,
    pub resources: Value,
    pub covariance_passed: bool,
    pub means_passed: bool,
}
const METHOD: &str = "proper-Gaussian-parameter-posterior/whitened-precision/v1";
pub(crate) fn evaluate(r: &Request) -> Result<Calculation, String> {
    let q = &r.resource_policy;
    let arithmetic = match q.arithmetic.as_str() {
        "binary64" => 0,
        "wide" => 1,
        _ => return Err("UNKNOWN_ARITHMETIC".into()),
    };
    let policy = GaussianPosteriorPolicy {
        struct_size: size_of::<GaussianPosteriorPolicy>() as u32,
        abi_version: ABI_VERSION,
        arithmetic,
        reserved: 0,
        maximum_elements: q.maximum_elements,
        maximum_cases: q.maximum_cases,
        maximum_string_bytes: q.maximum_string_bytes,
        maximum_native_bytes: q.maximum_native_bytes,
        maximum_work_units: q.maximum_work_units,
        maximum_forward_sensitivity: q.maximum_forward_sensitivity,
    };
    let noise_rows: Vec<_> = r.noise.ordered_row_ids.iter().map(|x| bytes(x)).collect();
    let noise_events: Vec<_> = r.noise.event_ids.iter().map(|x| bytes(x)).collect();
    let design_rows: Vec<_> = r.design.ordered_row_ids.iter().map(|x| bytes(x)).collect();
    let design_ids: Vec<_> = r
        .design
        .ordered_parameter_ids
        .iter()
        .map(|x| bytes(x))
        .collect();
    let design_units: Vec<_> = r.design.parameter_units.iter().map(|x| bytes(x)).collect();
    let column_units: Vec<_> = r.design.column_units.iter().map(|x| bytes(x)).collect();
    let prior_ids: Vec<_> = r
        .parameter_prior
        .ordered_parameter_ids
        .iter()
        .map(|x| bytes(x))
        .collect();
    let prior_units: Vec<_> = r
        .parameter_prior
        .parameter_units
        .iter()
        .map(|x| bytes(x))
        .collect();
    let shared: Vec<_> = r
        .parameter_prior
        .shared_nuisance_ids
        .iter()
        .map(|x| bytes(x))
        .collect();
    let condition_rows: Vec<_> = r
        .conditioning
        .ordered_row_ids
        .iter()
        .map(|x| bytes(x))
        .collect();
    let condition_events: Vec<_> = r.conditioning.event_ids.iter().map(|x| bytes(x)).collect();
    let cases: Vec<_> = r.conditioning.case_ids.iter().map(|x| bytes(x)).collect();
    let b = GaussianPosteriorBatch {
        struct_size: size_of::<GaussianPosteriorBatch>() as u32,
        abi_version: ABI_VERSION,
        source_semantics: 1,
        noise_independence_declared: r.parameter_prior.noise_independence_declared as u32,
        case_count: r.conditioning.vectors.lengths.len() as u64,
        noise_covariance: doubles(&r.noise.covariance_row_major),
        design: doubles(&r.design.values_row_major),
        prior_mean: doubles(&r.parameter_prior.mean),
        prior_covariance: doubles(&r.parameter_prior.covariance_row_major),
        conditioning_vectors: doubles(&r.conditioning.vectors.values),
        noise_row_ids: strings(&noise_rows),
        noise_event_ids: strings(&noise_events),
        design_row_ids: strings(&design_rows),
        design_parameter_ids: strings(&design_ids),
        design_parameter_units: strings(&design_units),
        column_units: strings(&column_units),
        prior_parameter_ids: strings(&prior_ids),
        prior_parameter_units: strings(&prior_units),
        shared_nuisance_ids: strings(&shared),
        conditioning_row_ids: strings(&condition_rows),
        conditioning_event_ids: strings(&condition_events),
        case_ids: strings(&cases),
        residual_unit: bytes(&r.noise.residual_unit),
        noise_identity: bytes(&r.noise.noise_identity),
        calibration_identity: bytes(&r.noise.calibration_identity),
        noise_dependence_identity: bytes(&r.noise.dependence_identity),
        ordering_provenance: bytes(&r.noise.ordering_provenance),
        design_identity: bytes(&r.design.design_identity),
        prior_identity: bytes(&r.parameter_prior.prior_identity),
        parameter_measure: bytes(&r.parameter_prior.parameter_measure),
        prior_dependence_identity: bytes(&r.parameter_prior.dependence_identity),
    };
    let mut raw = ptr::null_mut();
    let transport = unsafe { irred_gaussian_posterior_evaluate(&b, &policy, &mut raw) };
    let owner = Owned(raw);
    let requested_method = json!({"requested":METHOD,"actual":null,"executed":false});
    let requested_arithmetic = json!({"requested_source":q.arithmetic,"actual_source":null,"posterior_products":null,"output_storage":null});
    let mut resources = serde_json::to_value(q).map_err(|_| "RECORD_ENCODING")?;
    resources["compute_threads"] = json!(1);
    resources["native_byte_scope"] = json!(
        "charged result owner, native posterior/source ownership, conversion metadata, phase scratch and pooled mean/error/case views; excludes borrowed Rust request/wire buffers, JSON serialization, stack, allocator overhead and RSS"
    );
    if transport == GAUSSIAN_POSTERIOR_QUOTA_REFUSED {
        if !raw.is_null() {
            return Err("INVALID_QUOTA_OWNER".into());
        }
        return Ok(Calculation {
            output: json!({"kind":"failure","phase":"admission","status":"invalid_input","numerical_status":"work_limit","native_payload_present":false,"requested_case_count":b.case_count,"covariance":null,"rows":[],"meaning":"completed original-cap resource refusal; no native calculation payload"}),
            method: requested_method,
            arithmetic: requested_arithmetic,
            resources,
            covariance_passed: false,
            means_passed: false,
        });
    }
    if transport != OK {
        return Err(format!("CORE_STATUS_{transport}"));
    }
    if raw.is_null() {
        return Err("NULL_POSTERIOR_RESULT".into());
    }
    let mut v = std::mem::MaybeUninit::<GaussianPosteriorView>::uninit();
    if unsafe { irred_gaussian_posterior_result_view(owner.0, v.as_mut_ptr()) } != OK {
        return Err("INVALID_POSTERIOR_VIEW".into());
    }
    let v = unsafe { v.assume_init() };
    if v.struct_size != size_of::<GaussianPosteriorView>() as u32
        || v.abi_version != ABI_VERSION
        || v.reserved != 0
        || v.phase > 3
    {
        return Err("INVALID_POSTERIOR_VIEW".into());
    }
    let phase = [
        "admission",
        "noise_preparation",
        "posterior_preparation",
        "conditioning",
    ][v.phase as usize];
    let label = gaussian_status_name(v.status).ok_or("INVALID_SCIENTIFIC_STATUS")?;
    let numerical = numerical_status_name(v.numerical_status).ok_or("INVALID_NUMERICAL_STATUS")?;
    resources["actual"] = json!({"declared_work_units":v.declared_work_units,"pooled_numeric_elements":v.pooled_numeric_elements,"peak_payload_bytes":v.peak_payload_bytes,"retained_payload_bytes":v.retained_payload_bytes,"minimum_result_bytes":v.minimum_result_bytes,"phase":phase});
    let decode = |b: &Bytes| -> Result<String, String> {
        String::from_utf8(copied(b.data, b.length, 256)?)
            .map_err(|_| "INVALID_METADATA_ENCODING".into())
    };
    let finite = v.status == GAUSSIAN_STATUS_FINITE;
    if !finite {
        if v.case_count != 0 || v.covariance.length != 0 {
            return Err("FAILED_POSTERIOR_PAYLOAD".into());
        }
        return Ok(Calculation {
            output: json!({"kind":"failure","phase":phase,"status":label,"numerical_status":numerical,"native_payload_present":true,"requested_case_count":b.case_count,"covariance":null,"rows":[]}),
            method: requested_method,
            arithmetic: requested_arithmetic,
            resources,
            covariance_passed: false,
            means_passed: false,
        });
    }
    let p = r.parameter_prior.mean.len();
    let count = r.conditioning.vectors.lengths.len();
    if v.numerical_status != 0
        || v.case_count != count as u64
        || v.covariance.length != p.checked_mul(p).ok_or("RESOURCE_LIMIT")? as u64
        || !v.covariance_relative_error_estimate.is_finite()
        || v.covariance_relative_error_estimate < 0.0
    {
        return Err("INVALID_FINITE_POSTERIOR_VIEW".into());
    }
    fn values(b: &F64Buffer, expected: usize) -> Result<Vec<f64>, String> {
        if b.struct_size != size_of::<F64Buffer>() as u32
            || b.abi_version != ABI_VERSION
            || b.element_type != 2
            || b.reserved != 0
            || b.length != expected as u64
            || b.byte_length != b.length.checked_mul(8).ok_or("INVALID_RESULT_LENGTH")?
        {
            return Err("INVALID_POSTERIOR_BUFFER".into());
        }
        let v = copied(b.data, b.length, 1000000)?;
        if v.iter().any(|x| !x.is_finite()) {
            return Err("INVALID_FINITE_PAYLOAD".into());
        }
        Ok(v)
    }
    if count > 0
        && (v.rows.is_null()
            || (v.rows as usize) % std::mem::align_of::<GaussianPosteriorRow>() != 0)
    {
        return Err("INVALID_POSTERIOR_ROWS".into());
    }
    let mut out = Vec::with_capacity(count);
    let mut means_passed = true;
    for (i, row) in unsafe { std::slice::from_raw_parts(v.rows, count) }
        .iter()
        .enumerate()
    {
        if row.struct_size != size_of::<GaussianPosteriorRow>() as u32
            || row.abi_version != ABI_VERSION
            || row.case_index != i as u64
        {
            return Err("INVALID_POSTERIOR_ROW".into());
        }
        let status = gaussian_status_name(row.status).ok_or("INVALID_SCIENTIFIC_STATUS")?;
        let num = numerical_status_name(row.numerical_status).ok_or("INVALID_NUMERICAL_STATUS")?;
        if row.status != GAUSSIAN_STATUS_FINITE {
            means_passed = false;
            if row.mean.length != 0 || row.absolute_error_estimates.length != 0 {
                return Err("FAILED_POSTERIOR_PAYLOAD".into());
            }
            out.push(json!({"kind":"failure","case_index":i,"case_id":r.conditioning.case_ids[i],"status":status,"numerical_status":num}));
            continue;
        }
        if row.numerical_status != 0
            || !row.backward_residual.is_finite()
            || !row.estimated_forward_sensitivity.is_finite()
        {
            return Err("INVALID_FINITE_PAYLOAD".into());
        }
        let mean = values(&row.mean, p)?;
        let errors = values(&row.absolute_error_estimates, p)?;
        if errors.iter().any(|&x| x < 0.0) {
            return Err("INVALID_NUMERICAL_ERROR".into());
        }
        out.push(json!({"kind":"finite","case_index":i,"case_id":r.conditioning.case_ids[i],"status":status,"numerical_status":num,"mean":mean,"absolute_error_estimates":errors,"backward_residual":row.backward_residual,"estimated_forward_sensitivity":row.estimated_forward_sensitivity}));
    }
    let method = decode(&v.method)?;
    if method != METHOD {
        return Err("INVALID_POSTERIOR_METHOD".into());
    }
    Ok(Calculation {
        output: json!({"kind":if means_passed{"finite"}else{"failure"},"phase":phase,"status":label,"numerical_status":numerical,"native_payload_present":true,"ordered_row_ids":r.noise.ordered_row_ids,"event_ids":r.noise.event_ids,"ordered_parameter_ids":r.parameter_prior.ordered_parameter_ids,"parameter_units":r.parameter_prior.parameter_units,"parameter_measure":r.parameter_prior.parameter_measure,"covariance":values(&v.covariance,p*p)?,"covariance_relative_error_estimate":v.covariance_relative_error_estimate,"numerical_error_meaning":"empirical native arithmetic/conditioning diagnostics; distinct from physical posterior variance","rows":out,"batch_meaning":"family of conditionals under one fixed noise/design/prior; no IID or joint-across-cases law"}),
        method: json!({"requested":METHOD,"actual":method,"executed":true}),
        arithmetic: json!({"requested_source":q.arithmetic,"actual_source":decode(&v.arithmetic)?,"posterior_products":"longdouble-cpu/v1","output_storage":"binary64"}),
        resources,
        covariance_passed: true,
        means_passed,
    })
}
