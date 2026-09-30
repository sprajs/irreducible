//! One retained normalized-density owner; no projection or Gaussian Rust math.
use super::{
    supernova::{Projection, arithmetic},
    generated::*,
    observations::{bytes, doubles, strings},
};
use crate::{
    bao_ingestion::Decoded,
    model_spec::{ExpansionSpec, Geometry},
};
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};
use std::{ffi::c_void, mem::size_of, ptr};
#[derive(Clone, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct PreparationPolicy {
    pub arithmetic: String,
    pub maximum_queries: u64,
    pub maximum_matrix_elements: u64,
    pub maximum_string_bytes: u64,
    pub maximum_native_bytes: u64,
    pub maximum_forward_sensitivity: f64,
}
#[derive(Clone, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct Model {
    pub expansion: ExpansionSpec,
    pub geometry: Geometry,
    pub h0_rd_km_s: f64,
}
#[derive(Clone, Copy, Deserialize, Serialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum Output {
    NormalizedDensity,
    Predictions,
    Residuals,
}
impl Output {
    fn bit(self) -> u32 {
        match self {
            Self::NormalizedDensity => 1,
            Self::Predictions => 2,
            Self::Residuals => 4,
        }
    }
    pub(crate) fn id(self) -> &'static str {
        match self {
            Self::NormalizedDensity => "normalized_density",
            Self::Predictions => "predictions",
            Self::Residuals => "residuals",
        }
    }
}
#[derive(Clone, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct EvaluationPolicy {
    pub arithmetic: String,
    pub projection: Projection,
    pub maximum_models: u64,
    pub maximum_array_elements: u64,
    pub maximum_native_bytes: u64,
    pub maximum_forward_sensitivity: f64,
}
pub(crate) struct Prepared(pub(crate) *mut c_void);
impl Drop for Prepared {
    fn drop(&mut self) {
        unsafe {
            irred_bao_destroy(self.0);
        }
    }
}
pub(crate) struct Preparation {
    pub owner: Prepared,
    pub retained_bytes: usize,
    pub status: u32,
    pub numerical_status: u32,
    pub metadata: Value,
}
fn text(value: &Bytes) -> Result<String, String> {
    let n = usize::try_from(value.length).map_err(|_| "INVALID_NATIVE_STRING")?;
    if n > 16 * 1024 * 1024 || (n > 0 && value.data.is_null()) {
        return Err("INVALID_NATIVE_STRING".into());
    }
    let v = if n == 0 {
        &[][..]
    } else {
        unsafe { std::slice::from_raw_parts(value.data, n) }
    };
    String::from_utf8(v.to_vec()).map_err(|_| "INVALID_NATIVE_STRING".into())
}
pub(crate) fn prepare(source: &Decoded, policy: &PreparationPolicy) -> Result<Preparation, String> {
    let queries: Vec<_> = source
        .rows
        .iter()
        .map(|r| {
            Ok(BaoQuery {
                z: r.z,
                observable: match r.observable.as_str() {
                    "DM_over_rs" => 0,
                    "DH_over_rs" => 1,
                    "DV_over_rs" => 2,
                    _ => return Err("UNSUPPORTED_BAO_OBSERVABLE"),
                },
                reserved: 0,
            })
        })
        .collect::<Result<_, _>>()?;
    let observed: Vec<_> = source.rows.iter().map(|r| r.value).collect();
    let ids: Vec<_> = source.rows.iter().map(|r| bytes(&r.id)).collect();
    let axes: Vec<_> = source
        .covariance_axis_ids
        .iter()
        .map(|s| bytes(s))
        .collect();
    let descriptor = BaoSource {
        struct_size: size_of::<BaoSource>() as u32,
        abi_version: ABI_VERSION,
        role: if source.role == "released_fitted_distance_summary" {
            0
        } else {
            1
        },
        covariance_unit: 0,
        queries: queries.as_ptr(),
        query_count: queries.len() as u64,
        query_byte_length: std::mem::size_of_val(&*queries) as u64,
        observed: doubles(&observed),
        covariance: doubles(&source.covariance),
        ordered_ids: strings(&ids),
        covariance_axis_ids: strings(&axes),
        table_identity: bytes(source.table_asset.sha256()),
        covariance_identity: bytes(source.covariance_asset.sha256()),
        ordering_provenance: bytes(&source.ordering_provenance),
        calibration_provenance: bytes(&source.calibration_provenance),
        dependence_provenance: bytes(&source.dependence_provenance),
        redshift_convention: bytes(&source.redshift_convention),
        ruler_convention: bytes(&source.ruler_convention),
    };
    let p = BaoPreparationPolicy {
        struct_size: size_of::<BaoPreparationPolicy>() as u32,
        abi_version: ABI_VERSION,
        arithmetic: arithmetic(&policy.arithmetic)?,
        reserved: 0,
        maximum_queries: policy.maximum_queries,
        maximum_matrix_elements: policy.maximum_matrix_elements,
        maximum_string_bytes: policy.maximum_string_bytes,
        maximum_native_bytes: policy.maximum_native_bytes,
        maximum_forward_sensitivity: policy.maximum_forward_sensitivity,
    };
    let mut raw = ptr::null_mut();
    let code = unsafe { irred_bao_prepare(&descriptor, &p, &mut raw) };
    if code != OK {
        if !raw.is_null() {
            unsafe {
                irred_bao_destroy(raw);
            }
        }
        return Err(format!("CORE_STATUS_{code}"));
    }
    if raw.is_null() {
        return Err("NULL_NATIVE_OWNER".into());
    }
    let owner = Prepared(raw);
    let mut view: BaoView = unsafe { std::mem::zeroed() };
    let code = unsafe { irred_bao_source_view(owner.0, &mut view) };
    if code != OK
        || view.struct_size as usize != size_of::<BaoView>()
        || view.abi_version != ABI_VERSION
    {
        return Err("INVALID_NATIVE_VIEW".into());
    }
    Ok(Preparation {
        owner,
        retained_bytes: usize::try_from(view.retained_bytes).map_err(|_| "RESOURCE_LIMIT")?,
        status: view.status,
        numerical_status: view.numerical_status,
        metadata: json!({"arithmetic_id":text(&view.arithmetic_id)?,"density_id":text(&view.density_id)?,"equation_id":text(&view.equation_id)?,"ruler_convention_id":text(&view.ruler_convention_id)?}),
    })
}
struct ResultOwner(*mut c_void);
impl Drop for ResultOwner {
    fn drop(&mut self) {
        unsafe {
            irred_bao_result_destroy(self.0);
        }
    }
}
fn state(s: &OutputState) -> Value {
    json!({"availability":s.availability,"status":s.status,"numerical_status":s.numerical_status})
}
fn passed(s: &OutputState) -> bool {
    s.availability == 1 && s.status == 0 && s.numerical_status == 0
}
fn vector(v: &F64Buffer) -> Result<Vec<f64>, String> {
    if v.struct_size as usize != size_of::<F64Buffer>()
        || v.abi_version != ABI_VERSION
        || v.element_type != 2
        || v.reserved != 0
        || v.byte_length != v.length.checked_mul(8).ok_or("INVALID_NATIVE_ARRAY")?
    {
        return Err("INVALID_NATIVE_ARRAY".into());
    }
    let values = super::observations::copied(v.data, v.length, 1048576)?;
    if values.iter().any(|x| !x.is_finite()) {
        return Err("NONFINITE_NATIVE_PAYLOAD".into());
    }
    Ok(values)
}
pub(crate) fn evaluate(
    owner: &Prepared,
    models: &[Model],
    requested: &[Output],
    policy: &EvaluationPolicy,
) -> Result<(Value, Vec<(&'static str, bool)>), String> {
    let mut mask = 0;
    for want in requested {
        if mask & want.bit() != 0 {
            return Err("DUPLICATE_OUTPUT".into());
        }
        mask |= want.bit()
    }
    if mask == 0 {
        return Err("EMPTY_REQUESTED_OUTPUTS".into());
    }
    let params: Vec<_> = models
        .iter()
        .map(|m| m.expansion.active_parameters())
        .collect();
    let wire: Vec<_> = models
        .iter()
        .enumerate()
        .map(|(i, m)| BaoModel {
            struct_size: size_of::<BaoModel>() as u32,
            abi_version: ABI_VERSION,
            geometry: 0,
            reserved: 0,
            expansion: m.expansion.descriptor(&params[i]),
            h0_rd_km_s: m.h0_rd_km_s,
        })
        .collect();
    let batch = BaoBatch {
        struct_size: size_of::<BaoBatch>() as u32,
        abi_version: ABI_VERSION,
        models: wire.as_ptr(),
        model_count: wire.len() as u64,
        model_byte_length: std::mem::size_of_val(&*wire) as u64,
    };
    let p = BaoEvaluationPolicy {
        struct_size: size_of::<BaoEvaluationPolicy>() as u32,
        abi_version: ABI_VERSION,
        arithmetic: arithmetic(&policy.arithmetic)?,
        requested: mask,
        projection: policy.projection.descriptor(),
        maximum_models: policy.maximum_models,
        maximum_array_elements: policy.maximum_array_elements,
        maximum_native_bytes: policy.maximum_native_bytes,
        maximum_forward_sensitivity: policy.maximum_forward_sensitivity,
    };
    let mut raw = ptr::null_mut();
    let code = unsafe { irred_bao_evaluate(owner.0, &batch, &p, &mut raw) };
    if code != OK {
        if !raw.is_null() {
            unsafe {
                irred_bao_result_destroy(raw);
            }
        }
        return Err(format!("CORE_STATUS_{code}"));
    }
    if raw.is_null() {
        return Err("NULL_NATIVE_RESULT".into());
    }
    let result = ResultOwner(raw);
    let (mut data, mut count, mut status, mut numerical, mut callbacks, mut segments) =
        (ptr::null(), 0, 0, 0, 0, 0);
    let code = unsafe {
        irred_bao_result_view(
            result.0,
            &mut data,
            &mut count,
            &mut status,
            &mut numerical,
            &mut callbacks,
            &mut segments,
        )
    };
    if code != OK || count > models.len() as u64 || (count > 0 && data.is_null()) {
        return Err("INVALID_NATIVE_RESULT".into());
    }
    let rows = if count == 0 {
        &[][..]
    } else {
        unsafe { std::slice::from_raw_parts(data, count as usize) }
    };
    let complete = status == 0 && numerical == 0 && rows.len() == models.len();
    let mut checks: Vec<_> = requested.iter().map(|x| (x.id(), complete)).collect();
    let mut evaluations = vec![];
    for (i, row) in rows.iter().enumerate() {
        if row.struct_size as usize != size_of::<BaoRow>()
            || row.abi_version != ABI_VERSION
            || row.model_index != i as u64
            || row.source.expansion.model != models[i].expansion.tag()
            || row.source.geometry != 0
        {
            return Err("INVALID_NATIVE_RESULT".into());
        }
        let active = vector(&row.source.expansion.parameters)?;
        if active.len() != params[i].len()
            || active
                .iter()
                .zip(&params[i])
                .any(|(a, b)| a.to_bits() != b.to_bits())
            || row.source.h0_rd_km_s.to_bits() != models[i].h0_rd_km_s.to_bits()
        {
            return Err("NATIVE_MODEL_IDENTITY_MISMATCH".into());
        }
        let mut output = json!({"source":models[i],"background_status":row.background_status,"numerical_status":row.numerical_status,"model_id":text(&row.model_id)?,"arithmetic_id":text(&row.arithmetic_id)?,"work":{"callbacks":row.callbacks,"segment_visits":row.segment_visits}});
        for (j, want) in requested.iter().enumerate() {
            let s = match want {
                Output::NormalizedDensity => &row.density.state,
                Output::Predictions => &row.predictions_state,
                Output::Residuals => &row.residuals_state,
            };
            let ok = passed(s);
            checks[j].1 &= ok;
            let mut payload = json!({"state":state(s)});
            if ok {
                if matches!(want, Output::NormalizedDensity)
                    && [
                        row.density.quadratic,
                        row.density.log_determinant,
                        row.density.log_normalization,
                        row.density.log_density,
                        row.density.backward_residual,
                        row.density.estimated_forward_sensitivity,
                    ]
                    .iter()
                    .any(|v| !v.is_finite())
                {
                    return Err("NONFINITE_NATIVE_PAYLOAD".into());
                }
                payload["value"] = match want {
                    Output::NormalizedDensity => {
                        json!({"quadratic":row.density.quadratic,"log_determinant":row.density.log_determinant,"log_normalization":row.density.log_normalization,"log_density":row.density.log_density,"backward_residual":row.density.backward_residual,"estimated_forward_sensitivity":row.density.estimated_forward_sensitivity})
                    }
                    Output::Predictions => json!(vector(&row.predictions)?),
                    Output::Residuals => json!(vector(&row.residuals)?),
                };
            }
            output[want.id()] = payload;
        }
        evaluations.push(output);
    }
    Ok((
        json!({"kind":if checks.iter().all(|(_,p)|*p){"finite"}else{"failure"},"status":status,"numerical_status":numerical,"evaluations":evaluations,"work":{"callbacks":callbacks,"segment_visits":segments},"interpretation":"unqualified"}),
        checks,
    ))
}
