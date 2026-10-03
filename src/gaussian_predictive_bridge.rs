//! One borrowed coarse call; result views remain alive through JSON encoding.
use super::{
    generated::*,
    observations::{bytes, copied, doubles, strings},
};
use crate::gaussian_predictive_run::Request;
use serde_json::{Value, json};
use std::{
    mem::{align_of, size_of},
    ptr,
};
struct Owned(*mut std::ffi::c_void);
impl Drop for Owned {
    fn drop(&mut self) {
        unsafe {
            irred_gaussian_predictive_result_destroy(self.0);
        }
    }
}
pub(crate) struct Calculation {
    pub output: Value,
    pub method: Value,
    pub arithmetic: Value,
    pub resources: Value,
    pub means_passed: bool,
    pub densities_passed: bool,
}
const METHOD: &str = "proper-Gaussian-repeated-joint-predictive/retained-covariance/v1";
fn values(b: &F64Buffer, expected: usize) -> Result<&[f64], String> {
    if b.struct_size != size_of::<F64Buffer>() as u32
        || b.abi_version != ABI_VERSION
        || b.element_type != 2
        || b.reserved != 0
        || b.length != expected as u64
        || b.byte_length != b.length.checked_mul(8).ok_or("INVALID_RESULT_LENGTH")?
        || expected > 1000000
    {
        return Err("INVALID_PREDICTIVE_BUFFER".into());
    }
    if expected == 0 {
        if !b.data.is_null() {
            return Err("FAILED_PREDICTIVE_PAYLOAD".into());
        }
        return Ok(&[]);
    }
    if b.data.is_null() || b.data as usize % align_of::<f64>() != 0 {
        return Err("INVALID_PREDICTIVE_BUFFER".into());
    }
    let x = unsafe { std::slice::from_raw_parts(b.data, expected) };
    if x.iter().any(|x| !x.is_finite()) {
        return Err("INVALID_FINITE_PAYLOAD".into());
    }
    Ok(x)
}
fn groups(r: &Request, m: bool, d: bool) -> Value {
    json!({"predictive_means":if !r.outputs.means{"not_requested"}else if m{"checks_passed"}else{"failed"},"joint_predictive_log_densities":if !r.outputs.joint_log_densities{"not_requested"}else if d{"checks_passed"}else{"failed"}})
}
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
    let columns: Vec<_> = r.design.column_units.iter().map(|x| bytes(x)).collect();
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
    let future_rows: Vec<_> = r
        .future_noise
        .ordered_row_ids
        .iter()
        .map(|x| bytes(x))
        .collect();
    let future_events: Vec<_> = r.future_noise.event_ids.iter().map(|x| bytes(x)).collect();
    let response_rows: Vec<_> = r
        .future_response
        .ordered_row_ids
        .iter()
        .map(|x| bytes(x))
        .collect();
    let response_ids: Vec<_> = r
        .future_response
        .ordered_parameter_ids
        .iter()
        .map(|x| bytes(x))
        .collect();
    let response_units: Vec<_> = r
        .future_response
        .parameter_units
        .iter()
        .map(|x| bytes(x))
        .collect();
    let future_columns: Vec<_> = r
        .future_response
        .column_units
        .iter()
        .map(|x| bytes(x))
        .collect();
    let vector_rows: Vec<_> = r.future_vectors.as_ref().map_or_else(Vec::new, |f| {
        f.ordered_row_ids.iter().map(|x| bytes(x)).collect()
    });
    let vector_events: Vec<_> = r
        .future_vectors
        .as_ref()
        .map_or_else(Vec::new, |f| f.event_ids.iter().map(|x| bytes(x)).collect());
    let training = GaussianPosteriorBatch {
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
        column_units: strings(&columns),
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
    let mask = r.outputs.means as u32 + 2 * r.outputs.joint_log_densities as u32;
    let batch = GaussianPredictiveBatch {
        struct_size: size_of::<GaussianPredictiveBatch>() as u32,
        abi_version: ABI_VERSION,
        requested_outputs: mask,
        reserved: 0,
        future_noise_independence_declared: r.prediction.future_noise_independence_declared as u32,
        noise_conditional_on_parameters_declared: r
            .prediction
            .noise_conditional_on_parameters_declared
            as u32,
        training,
        future_noise_covariance: doubles(&r.future_noise.covariance_row_major),
        future_response: doubles(&r.future_response.values_row_major),
        future_vectors: doubles(
            r.future_vectors
                .as_ref()
                .map_or(&[], |f| f.vectors.values.as_slice()),
        ),
        future_noise_row_ids: strings(&future_rows),
        future_noise_event_ids: strings(&future_events),
        future_response_row_ids: strings(&response_rows),
        future_parameter_ids: strings(&response_ids),
        future_parameter_units: strings(&response_units),
        future_column_units: strings(&future_columns),
        future_vector_row_ids: strings(&vector_rows),
        future_vector_event_ids: strings(&vector_events),
        future_unit: bytes(&r.future_noise.residual_unit),
        future_noise_identity: bytes(&r.future_noise.noise_identity),
        future_calibration_identity: bytes(&r.future_noise.calibration_identity),
        future_noise_dependence_identity: bytes(&r.future_noise.dependence_identity),
        future_ordering_provenance: bytes(&r.future_noise.ordering_provenance),
        future_response_identity: bytes(&r.future_response.design_identity),
        future_covariance_unit: bytes(&r.prediction.future_covariance_unit),
        future_measure: bytes(&r.prediction.future_measure),
        conditioning_identity: bytes(&r.prediction.conditioning_identity),
        prediction_dependence_identity: bytes(&r.prediction.dependence_identity),
    };
    let count = r.conditioning.vectors.lengths.len();
    let mut raw = ptr::null_mut();
    let transport = unsafe { irred_gaussian_predictive_evaluate(&batch, &policy, &mut raw) };
    let owner = Owned(raw);
    let requested_method = json!({"requested":METHOD,"actual":null,"executed":false});
    let requested_arithmetic = json!({"requested_source":q.arithmetic,"actual_training_source":null,"actual_future_source":null,"posterior_products":null,"predictive_products":null,"output_storage":null});
    let mut resources = serde_json::to_value(q).map_err(|_| "RECORD_ENCODING")?;
    resources["compute_threads"] = json!(1);
    resources["native_byte_scope"] = json!(
        "result/native owners/conversion/phase scratch/pools/views; excludes borrowed Rust request/wire buffers, JSON serialization, stack, allocator bookkeeping and RSS"
    );
    if transport == GAUSSIAN_PREDICTIVE_QUOTA_REFUSED {
        if !raw.is_null() {
            return Err("INVALID_QUOTA_OWNER".into());
        }
        return Ok(Calculation {
            output: json!({"kind":"failure","phase":"admission","status":"invalid_input","numerical_status":"work_limit","native_payload_present":false,"requested_case_count":count,"groups":groups(r,false,false),"rows":[]}),
            method: requested_method,
            arithmetic: requested_arithmetic,
            resources,
            means_passed: false,
            densities_passed: false,
        });
    }
    if transport != OK {
        return Err(format!("CORE_STATUS_{transport}"));
    }
    if raw.is_null() {
        return Err("NULL_PREDICTIVE_RESULT".into());
    }
    let mut v = std::mem::MaybeUninit::<GaussianPredictiveView>::uninit();
    if unsafe { irred_gaussian_predictive_result_view(owner.0, v.as_mut_ptr()) } != OK {
        return Err("INVALID_PREDICTIVE_VIEW".into());
    }
    let v = unsafe { v.assume_init() };
    let flags = [
        v.training_noise_attempted,
        v.training_noise_prepared,
        v.posterior_attempted,
        v.posterior_prepared,
        v.future_noise_attempted,
        v.future_noise_prepared,
        v.predictive_attempted,
        v.predictive_prepared,
        v.batch_called,
        v.batch_admission_completed,
    ];
    if v.struct_size != size_of::<GaussianPredictiveView>() as u32
        || v.abi_version != ABI_VERSION
        || v.requested_outputs != mask
        || v.phase > 5
        || flags.iter().any(|&x| x > 1)
        || flags.windows(2).any(|w| w[1] > w[0])
        || v.peak_payload_bytes > q.maximum_native_bytes
        || v.retained_payload_bytes > v.peak_payload_bytes
        || v.minimum_result_bytes > v.retained_payload_bytes
    {
        return Err("INVALID_PREDICTIVE_VIEW".into());
    }
    let phase = [
        "admission",
        "training_noise_preparation",
        "posterior_preparation",
        "future_noise_preparation",
        "predictive_preparation",
        "coarse_conditioning",
    ][v.phase as usize];
    let label = gaussian_status_name(v.status).ok_or("INVALID_SCIENTIFIC_STATUS")?;
    let num = numerical_status_name(v.numerical_status).ok_or("INVALID_NUMERICAL_STATUS")?;
    let work = [
        v.declared_training_noise_work,
        v.declared_posterior_work,
        v.declared_future_noise_work,
        v.declared_predictive_work,
        v.declared_batch_work,
    ]
    .into_iter()
    .try_fold(0u64, |a, b| a.checked_add(b))
    .ok_or("INVALID_WORK_RECEIPT")?;
    if work != v.declared_total_work {
        return Err("INVALID_WORK_RECEIPT".into());
    }
    resources["actual"] = json!({"phase":phase,"declared_training_noise_work":v.declared_training_noise_work,"declared_posterior_work":v.declared_posterior_work,"declared_future_noise_work":v.declared_future_noise_work,"declared_predictive_work":v.declared_predictive_work,"declared_batch_work":v.declared_batch_work,"declared_total_work":v.declared_total_work,"reservation_available":v.declared_total_work>0,"pooled_numeric_elements":v.pooled_numeric_elements,"minimum_result_bytes":v.minimum_result_bytes,"peak_payload_bytes":v.peak_payload_bytes,"retained_payload_bytes":v.retained_payload_bytes});
    let decode = |b: &Bytes| -> Result<String, String> {
        String::from_utf8(copied(b.data, b.length, 256)?)
            .map_err(|_| "INVALID_METADATA_ENCODING".into())
    };
    let finite = v.status == GAUSSIAN_STATUS_FINITE;
    let method = if v.training_noise_attempted == 1 {
        let actual = decode(&v.method)?;
        let expected = if v.predictive_attempted == 1 {
            METHOD
        } else if v.posterior_attempted == 1 && v.future_noise_attempted == 0 {
            "proper-Gaussian-parameter-posterior/whitened-precision/v1"
        } else {
            "supplied-Gaussian-covariance-preparation/Cholesky/v1"
        };
        if actual != expected {
            return Err("INVALID_PREDICTIVE_METHOD".into());
        }
        json!({"requested":METHOD,"actual":actual,"executed":true,"meaning":"attempted stage algorithm; flags state completion separately","training_noise_attempted":v.training_noise_attempted==1,"training_noise_prepared":v.training_noise_prepared==1,"posterior_attempted":v.posterior_attempted==1,"posterior_prepared":v.posterior_prepared==1,"future_noise_attempted":v.future_noise_attempted==1,"future_noise_prepared":v.future_noise_prepared==1,"predictive_attempted":v.predictive_attempted==1,"predictive_prepared":v.predictive_prepared==1,"batch_called":v.batch_called==1,"batch_admission_completed":v.batch_admission_completed==1})
    } else {
        requested_method
    };
    let arithmetic = if v.training_noise_attempted == 1 {
        let source = decode(&v.arithmetic)?;
        let expected = if q.arithmetic == "wide" {
            "F02/longdouble-cpu/v1"
        } else {
            "F02/binary64-legacy/v1"
        };
        if source != expected {
            return Err("INVALID_PREDICTIVE_ARITHMETIC".into());
        }
        json!({"requested_source":q.arithmetic,"attempted_source":source,"actual_training_source":if v.training_noise_prepared==1{Some(source.as_str())}else{None},"actual_future_source":if v.future_noise_prepared==1{Some(source.as_str())}else{None},"posterior_products":if v.posterior_prepared==1{Some("longdouble-cpu/v1")}else{None},"predictive_products":if v.predictive_prepared==1{Some("longdouble-cpu/v1")}else{None},"output_storage":if v.batch_admission_completed==1{Some("binary64")}else{None}})
    } else {
        requested_arithmetic
    };
    if !finite {
        if v.case_count != 0 || !v.rows.is_null() {
            return Err("FAILED_PREDICTIVE_PAYLOAD".into());
        }
        return Ok(Calculation {
            output: json!({"kind":"failure","phase":phase,"status":label,"numerical_status":num,"native_payload_present":true,"requested_case_count":count,"groups":groups(r,false,false),"rows":[]}),
            method,
            arithmetic,
            resources,
            means_passed: false,
            densities_passed: false,
        });
    }
    if v.phase != 5
        || flags.iter().any(|&x| x != 1)
        || v.numerical_status != 0
        || v.case_count != count as u64
        || count == 0
        || count > 65536
        || v.rows.is_null()
        || v.rows as usize % align_of::<GaussianPredictiveRow>() != 0
        || v.declared_total_work > q.maximum_work_units
        || v.pooled_numeric_elements > q.maximum_elements
    {
        return Err("INVALID_FINITE_EXECUTION".into());
    }
    let k = r.future_noise.ordered_row_ids.len();
    let mut rows = Vec::with_capacity(count);
    let mut means_passed = true;
    let mut densities_passed = true;
    for (i, row) in unsafe { std::slice::from_raw_parts(v.rows, count) }
        .iter()
        .enumerate()
    {
        if row.struct_size != size_of::<GaussianPredictiveRow>() as u32
            || row.abi_version != ABI_VERSION
            || row.case_index != i as u64
            || row.reserved != 0
            || row.joint_density_available > 1
        {
            return Err("INVALID_PREDICTIVE_ROW".into());
        }
        let status = gaussian_status_name(row.status).ok_or("INVALID_SCIENTIFIC_STATUS")?;
        let num = numerical_status_name(row.numerical_status).ok_or("INVALID_NUMERICAL_STATUS")?;
        if row.status != GAUSSIAN_STATUS_FINITE {
            means_passed = false;
            densities_passed = false;
            values(&row.mean, 0)?;
            values(&row.absolute_error_estimates, 0)?;
            if row.joint_density_available != 0 {
                return Err("FAILED_PREDICTIVE_PAYLOAD".into());
            }
            rows.push(json!({"kind":"failure","case_index":i,"case_id":r.conditioning.case_ids[i],"status":status,"numerical_status":num}));
            continue;
        }
        if row.numerical_status != 0 {
            return Err("INVALID_FINITE_PAYLOAD".into());
        }
        let mut out = json!({"kind":"finite","case_index":i,"case_id":r.conditioning.case_ids[i],"status":status,"numerical_status":num});
        if r.outputs.means {
            let mean = values(&row.mean, k)?;
            let error = values(&row.absolute_error_estimates, k)?;
            if error
                .iter()
                .zip(mean)
                .any(|(&e, &m)| e < 0.0 || e / (1.0 + m.abs()) > q.maximum_forward_sensitivity)
            {
                return Err("INVALID_NUMERICAL_ERROR".into());
            }
            out["future_mean"] = json!(mean);
            out["absolute_error_estimates"] = json!(error);
        } else {
            values(&row.mean, 0)?;
            values(&row.absolute_error_estimates, 0)?;
        }
        if r.outputs.joint_log_densities {
            if row.joint_density_available != 1
                || [
                    row.log_density,
                    row.quadratic,
                    row.log_determinant,
                    row.normalization,
                    row.backward_residual,
                    row.estimated_forward_sensitivity,
                ]
                .iter()
                .any(|x| !x.is_finite())
                || row.quadratic < 0.0
                || row.normalization < 0.0
                || row.backward_residual < 0.0
                || row.estimated_forward_sensitivity < 0.0
                || row.estimated_forward_sensitivity > q.maximum_forward_sensitivity
            {
                return Err("INVALID_FINITE_PAYLOAD".into());
            }
            out["joint_density"] = json!({"log_density":row.log_density,"quadratic":row.quadratic,"log_determinant":row.log_determinant,"normalization":row.normalization,"backward_residual":row.backward_residual,"estimated_forward_sensitivity":row.estimated_forward_sensitivity});
        } else if row.joint_density_available != 0 {
            return Err("UNREQUESTED_PREDICTIVE_PAYLOAD".into());
        }
        rows.push(out);
    }
    Ok(Calculation {
        output: json!({"kind":if means_passed&&densities_passed{"finite"}else{"failure"},"phase":phase,"status":label,"numerical_status":num,"native_payload_present":true,"future_unit":r.future_noise.residual_unit,"future_measure":r.prediction.future_measure,"ordered_future_row_ids":r.future_noise.ordered_row_ids,"future_event_ids":r.future_noise.event_ids,"ordered_parameter_ids":r.parameter_prior.ordered_parameter_ids,"parameter_units":r.parameter_prior.parameter_units,"groups":groups(r,means_passed,densities_passed),"rows":rows,"numerical_error_meaning":"native empirical arithmetic diagnostics; distinct from physical future covariance","batch_meaning":"family of full-vector conditionals; no IID or joint-across-cases law"}),
        method,
        arithmetic,
        resources,
        means_passed,
        densities_passed,
    })
}
