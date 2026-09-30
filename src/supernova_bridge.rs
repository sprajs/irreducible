//! Current selected-source SN owner. Acquisition is never repeated here.
use super::{
    generated::*,
    observations::{Prepared as Observations, copied},
};
use crate::model_spec::{Observer, ObserverConvention};
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};
use std::{ffi::c_void, mem::size_of, ptr};
#[derive(Clone, Deserialize, Serialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub(crate) enum Selection {
    Explicit {
        source_indices: Vec<u64>,
        coordinates: Vec<Coordinate>,
    },
    PantheonZhdGt001,
}
#[derive(Clone, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct Coordinate {
    pub z_expansion: f64,
    pub observer: Observer,
}
#[derive(Clone, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct PreparationPolicy {
    pub arithmetic: String,
    pub maximum_selected_rows: u64,
    pub maximum_matrix_elements: u64,
    pub maximum_string_bytes: u64,
    pub maximum_native_bytes: u64,
    pub maximum_forward_sensitivity: f64,
}
pub(crate) fn arithmetic(value: &str) -> Result<u32, String> {
    match value {
        "binary64" => Ok(0),
        "wide" => Ok(1),
        _ => Err("INVALID_ARITHMETIC".into()),
    }
}
pub(crate) struct Prepared(pub(crate) *mut c_void);
impl Drop for Prepared {
    fn drop(&mut self) {
        unsafe {
            irred_supernova_destroy(self.0);
        }
    }
}
pub(crate) struct Preparation {
    pub owner: Prepared,
    pub selected: Value,
    pub retained_bytes: usize,
    pub status: u32,
    pub numerical_status: u32,
}
fn text(value: &Bytes) -> Result<String, String> {
    let n = usize::try_from(value.length).map_err(|_| "INVALID_NATIVE_STRING")?;
    if n > 16 * 1024 * 1024 || (n > 0 && value.data.is_null()) {
        return Err("INVALID_NATIVE_STRING".into());
    }
    let bytes = if n == 0 {
        &[][..]
    } else {
        unsafe { std::slice::from_raw_parts(value.data, n) }
    };
    String::from_utf8(bytes.to_vec()).map_err(|_| "INVALID_NATIVE_STRING".into())
}
impl Prepared {
    fn view(&self) -> Result<SupernovaView, String> {
        let mut view: SupernovaView = unsafe { std::mem::zeroed() };
        let code = unsafe { irred_supernova_source_view(self.0, &mut view) };
        if code != OK {
            return Err(format!("CORE_STATUS_{code}"));
        }
        if view.abi_version != ABI_VERSION
            || view.struct_size as usize != size_of::<SupernovaView>()
        {
            return Err("INVALID_NATIVE_VIEW".into());
        }
        Ok(view)
    }
}
pub(crate) fn prepare(
    source: &Observations,
    selection: &Selection,
    policy: &PreparationPolicy,
) -> Result<Preparation, String> {
    let coordinates: Vec<_> = match selection {
        Selection::Explicit { coordinates, .. } => coordinates
            .iter()
            .map(|c| MagnitudeCoordinate {
                struct_size: size_of::<MagnitudeCoordinate>() as u32,
                abi_version: ABI_VERSION,
                z_expansion: c.z_expansion,
                observer_redshift: c.observer.redshift,
                observer_convention: match c.observer.convention {
                    ObserverConvention::GeometricSameRedshift => 0,
                    ObserverConvention::ReleasedZhdZhel => 1,
                },
                reserved: 0,
            })
            .collect(),
        _ => vec![],
    };
    let indices = match selection {
        Selection::Explicit { source_indices, .. } => source_indices.as_slice(),
        _ => &[],
    };
    let descriptor = MagnitudeSelection {
        struct_size: size_of::<MagnitudeSelection>() as u32,
        abi_version: ABI_VERSION,
        kind: if matches!(selection, Selection::Explicit { .. }) {
            0
        } else {
            1
        },
        reserved: 0,
        source_indices: U64Buffer {
            data: indices.as_ptr(),
            length: indices.len() as u64,
            byte_length: std::mem::size_of_val(indices) as u64,
        },
        coordinates: coordinates.as_ptr(),
        coordinate_count: coordinates.len() as u64,
        coordinate_byte_length: std::mem::size_of_val(&*coordinates) as u64,
    };
    let p = SupernovaPreparationPolicy {
        struct_size: size_of::<SupernovaPreparationPolicy>() as u32,
        abi_version: ABI_VERSION,
        arithmetic: arithmetic(&policy.arithmetic)?,
        reserved: 0,
        maximum_selected_rows: policy.maximum_selected_rows,
        maximum_matrix_elements: policy.maximum_matrix_elements,
        maximum_string_bytes: policy.maximum_string_bytes,
        maximum_native_bytes: policy.maximum_native_bytes,
        maximum_forward_sensitivity: policy.maximum_forward_sensitivity,
    };
    let mut raw = ptr::null_mut();
    let code = unsafe { irred_supernova_prepare(source.0, &descriptor, &p, &mut raw) };
    if code != OK {
        if !raw.is_null() {
            unsafe {
                irred_supernova_destroy(raw);
            }
        }
        return Err(format!("CORE_STATUS_{code}"));
    }
    if raw.is_null() {
        return Err("NULL_NATIVE_OWNER".into());
    }
    let owner = Prepared(raw);
    let view = owner.view()?;
    let indices = copied(
        view.selected_source_indices.data,
        view.selected_source_indices.length,
        4096,
    )?;
    if view.selected_source_indices.byte_length
        != view
            .selected_source_indices
            .length
            .checked_mul(8)
            .ok_or("INVALID_NATIVE_VIEW")?
    {
        return Err("INVALID_NATIVE_VIEW".into());
    }
    if view.coordinate_count > 4096
        || view.ordered_ids.length > 4096
        || view.coordinate_count != view.selected_source_indices.length
        || view.ordered_ids.length != view.selected_source_indices.length
        || view.ordered_ids.byte_length
            != view
                .ordered_ids
                .length
                .checked_mul(size_of::<Bytes>() as u64)
                .ok_or("INVALID_NATIVE_VIEW")?
        || (view.coordinates.is_null() && view.coordinate_count != 0)
    {
        return Err("INVALID_NATIVE_VIEW".into());
    }
    if view.coordinate_byte_length
        != view
            .coordinate_count
            .checked_mul(size_of::<MagnitudeCoordinate>() as u64)
            .ok_or("INVALID_NATIVE_VIEW")?
    {
        return Err("INVALID_NATIVE_VIEW".into());
    }
    let coordinates = if view.coordinate_count == 0 {
        &[][..]
    } else {
        unsafe {
            std::slice::from_raw_parts(
                view.coordinates,
                usize::try_from(view.coordinate_count).map_err(|_| "INVALID_NATIVE_VIEW")?,
            )
        }
    };
    if coordinates.len() > 4096 || coordinates.len() != indices.len() {
        return Err("INVALID_NATIVE_VIEW".into());
    }
    let ids = unsafe {
        if view.ordered_ids.length == 0 {
            &[][..]
        } else {
            if view.ordered_ids.data.is_null() {
                return Err("INVALID_NATIVE_VIEW".into());
            }
            std::slice::from_raw_parts(
                view.ordered_ids.data,
                usize::try_from(view.ordered_ids.length).map_err(|_| "INVALID_NATIVE_VIEW")?,
            )
        }
    };
    if ids.len() != indices.len() {
        return Err("INVALID_NATIVE_VIEW".into());
    }
    let ids: Vec<_> = ids.iter().map(text).collect::<Result<_, _>>()?;
    let coordinates:Vec<_>=coordinates.iter().map(|c|json!({"z_expansion":c.z_expansion,"observer":{"redshift":c.observer_redshift,"convention":c.observer_convention}})).collect();
    let selected = json!({"source_indices":indices,"ordered_ids":ids,"coordinates":coordinates,
        "arithmetic_id":text(&view.arithmetic_id)?,"score_id":text(&view.score_id)?,
        "shape_convention":text(&view.shape_convention)?,"offset_convention":text(&view.offset_convention)?,
        "source_role":view.source.role});
    Ok(Preparation {
        owner,
        selected,
        retained_bytes: usize::try_from(view.retained_bytes).map_err(|_| "RESOURCE_LIMIT")?,
        status: view.status,
        numerical_status: view.preparation_numerical_status,
    })
}
#[derive(Clone, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct Projection {
    pub integration: Option<super::expansion::Integration>,
    pub maximum_queries: u64,
    pub maximum_callbacks: u64,
    pub maximum_segment_visits: u64,
}
impl Projection {
    pub(super) fn descriptor(&self) -> ProjectionPolicy {
        ProjectionPolicy {
            struct_size: size_of::<ProjectionPolicy>() as u32,
            abi_version: ABI_VERSION,
            has_integration: self.integration.is_some() as u32,
            maximum_depth: self.integration.as_ref().map_or(0, |p| p.maximum_depth),
            absolute_tolerance: self
                .integration
                .as_ref()
                .map_or(0., |p| p.absolute_tolerance),
            relative_tolerance: self
                .integration
                .as_ref()
                .map_or(0., |p| p.relative_tolerance),
            maximum_evaluations_per_integral: self
                .integration
                .as_ref()
                .map_or(0, |p| p.maximum_evaluations),
            maximum_queries: self.maximum_queries,
            maximum_callbacks: self.maximum_callbacks,
            maximum_segment_visits: self.maximum_segment_visits,
        }
    }
}
#[derive(Clone, Copy, Deserialize, Serialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum Output {
    Score,
    GeometricShape,
    MagnitudeEffect,
    CorrectedResiduals,
    ProfiledResiduals,
    Diagnostics,
}
impl Output {
    fn bit(self) -> u32 {
        match self {
            Self::Score => 1,
            Self::GeometricShape => 2,
            Self::MagnitudeEffect => 4,
            Self::CorrectedResiduals => 8,
            Self::ProfiledResiduals => 16,
            Self::Diagnostics => 32,
        }
    }
    pub(crate) fn id(self) -> &'static str {
        match self {
            Self::Score => "score",
            Self::GeometricShape => "geometric_shape",
            Self::MagnitudeEffect => "magnitude_effect",
            Self::CorrectedResiduals => "corrected_residuals",
            Self::ProfiledResiduals => "profiled_residuals",
            Self::Diagnostics => "diagnostics",
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
struct ResultOwner(*mut c_void);
impl Drop for ResultOwner {
    fn drop(&mut self) {
        unsafe {
            irred_supernova_result_destroy(self.0);
        }
    }
}
fn state(value: &OutputState) -> Value {
    json!({"availability":value.availability,"status":value.status,"numerical_status":value.numerical_status})
}
fn passed(value: &OutputState) -> bool {
    value.availability == 1 && value.status == 0 && value.numerical_status == 0
}
fn vector(value: &F64Buffer) -> Result<Vec<f64>, String> {
    if value.struct_size as usize != size_of::<F64Buffer>()
        || value.abi_version != ABI_VERSION
        || value.element_type != 2
        || value.reserved != 0
        || value.byte_length != value.length.checked_mul(8).ok_or("INVALID_NATIVE_ARRAY")?
    {
        return Err("INVALID_NATIVE_ARRAY".into());
    }
    let v = copied(value.data, value.length, 1048576)?;
    if v.iter().any(|x| !x.is_finite()) {
        return Err("NONFINITE_NATIVE_PAYLOAD".into());
    }
    Ok(v)
}
pub(crate) fn evaluate(
    owner: &Prepared,
    models: &[crate::model_spec::SupernovaModel],
    requested: &[Output],
    policy: &EvaluationPolicy,
) -> Result<(Value, Vec<(&'static str, bool)>), String> {
    let mut mask = 0;
    for output in requested {
        if mask & output.bit() != 0 {
            return Err("DUPLICATE_OUTPUT".into());
        }
        mask |= output.bit();
    }
    if mask == 0 {
        return Err("EMPTY_REQUESTED_OUTPUTS".into());
    }
    let parameters: Vec<_> = models
        .iter()
        .map(|m| m.expansion.active_parameters())
        .collect();
    let effects: Vec<Vec<f64>> = models
        .iter()
        .map(|m| match m.source_effect {
            crate::model_spec::SourceEffect::None => vec![],
            crate::model_spec::SourceEffect::GreyLog1pMagnitude { epsilon_mag } => {
                vec![epsilon_mag]
            }
        })
        .collect();
    let wire: Vec<_> = models
        .iter()
        .enumerate()
        .map(|(i, m)| SupernovaModel {
            struct_size: size_of::<SupernovaModel>() as u32,
            abi_version: ABI_VERSION,
            geometry: 0,
            reserved: 0,
            expansion: m.expansion.descriptor(&parameters[i]),
            source_effect: SourceEffectSpec {
                struct_size: size_of::<SourceEffectSpec>() as u32,
                abi_version: ABI_VERSION,
                effect: if effects[i].is_empty() { 0 } else { 1 },
                reserved: 0,
                parameters: super::observations::doubles(&effects[i]),
            },
        })
        .collect();
    let batch = SupernovaBatch {
        struct_size: size_of::<SupernovaBatch>() as u32,
        abi_version: ABI_VERSION,
        models: wire.as_ptr(),
        model_count: wire.len() as u64,
        model_byte_length: std::mem::size_of_val(&*wire) as u64,
    };
    let policy = SupernovaEvaluationPolicy {
        struct_size: size_of::<SupernovaEvaluationPolicy>() as u32,
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
    let code = unsafe { irred_supernova_evaluate(owner.0, &batch, &policy, &mut raw) };
    if code != OK {
        if !raw.is_null() {
            unsafe {
                irred_supernova_result_destroy(raw);
            }
        }
        return Err(format!("CORE_STATUS_{code}"));
    }
    if raw.is_null() {
        return Err("NULL_NATIVE_RESULT".into());
    }
    let result = ResultOwner(raw);
    let (mut rows, mut count, mut status, mut numerical, mut callbacks, mut segments) =
        (ptr::null(), 0, 0, 0, 0, 0);
    let code = unsafe {
        irred_supernova_result_view(
            result.0,
            &mut rows,
            &mut count,
            &mut status,
            &mut numerical,
            &mut callbacks,
            &mut segments,
        )
    };
    if code != OK || count > models.len() as u64 || (count > 0 && rows.is_null()) {
        return Err("INVALID_NATIVE_RESULT".into());
    }
    let rows = if count == 0 {
        &[][..]
    } else {
        unsafe { std::slice::from_raw_parts(rows, count as usize) }
    };
    let complete = status == 0 && numerical == 0 && rows.len() == models.len();
    let mut checks: Vec<_> = requested.iter().map(|x| (x.id(), complete)).collect();
    let mut outputs = vec![];
    for (i, row) in rows.iter().enumerate() {
        if row.struct_size as usize != size_of::<SupernovaRow>()
            || row.abi_version != ABI_VERSION
            || row.model_index != i as u64
            || row.source.expansion.model != models[i].expansion.tag()
            || row.source.geometry != 0
            || row.source.source_effect.effect != if effects[i].is_empty() { 0 } else { 1 }
        {
            return Err("INVALID_NATIVE_RESULT".into());
        }
        let actual = vector(&row.source.expansion.parameters)?;
        if actual.len() != parameters[i].len()
            || actual
                .iter()
                .zip(&parameters[i])
                .any(|(a, b)| a.to_bits() != b.to_bits())
        {
            return Err("NATIVE_MODEL_IDENTITY_MISMATCH".into());
        }
        let effect = vector(&row.source.source_effect.parameters)?;
        if effect.len() != effects[i].len()
            || effect
                .iter()
                .zip(&effects[i])
                .any(|(a, b)| a.to_bits() != b.to_bits())
        {
            return Err("NATIVE_MODEL_IDENTITY_MISMATCH".into());
        }
        let mut output = json!({"source":models[i],"status":row.status,"background_status":row.background_status,"numerical_status":row.numerical_status,"profile_status":row.profile_status,"model_id":text(&row.model_id)?,"hypothesis_id":text(&row.hypothesis_id)?,"arithmetic_id":text(&row.arithmetic_id)?,"work":{"callbacks":row.callbacks,"segment_visits":row.segment_visits}});
        for (j, want) in requested.iter().enumerate() {
            let s = match want {
                Output::Score => &row.score.state,
                Output::GeometricShape => &row.geometry_state,
                Output::MagnitudeEffect => &row.effect_state,
                Output::CorrectedResiduals => &row.corrected_state,
                Output::ProfiledResiduals => &row.profiled_state,
                Output::Diagnostics => &row.diagnostics.state,
            };
            let ok = passed(s);
            checks[j].1 &= ok;
            let mut payload = json!({"state":state(s)});
            if ok {
                let scalars: &[f64] = match want {
                    Output::Score => &[
                        row.score.relative_profile_score,
                        row.score.quadratic,
                        row.score.offset_coefficient,
                    ],
                    Output::Diagnostics => &[
                        row.diagnostics.backward_residual,
                        row.diagnostics.estimated_forward_sensitivity,
                        row.diagnostics.coefficient_backward_residual,
                        row.diagnostics.coefficient_forward_sensitivity,
                        row.diagnostics.residual_l1,
                        row.diagnostics.solution_norm_inf,
                        row.diagnostics.adjusted_residual_l1,
                        row.diagnostics.adjusted_solution_norm_inf,
                    ],
                    _ => &[],
                };
                if scalars.iter().any(|v| !v.is_finite()) {
                    return Err("NONFINITE_NATIVE_PAYLOAD".into());
                }
                let value = match want {
                    Output::Score => {
                        json!({"relative_profile_score":row.score.relative_profile_score,"quadratic":row.score.quadratic,"offset_coefficient":row.score.offset_coefficient})
                    }
                    Output::GeometricShape => json!(vector(&row.geometric_shape)?),
                    Output::MagnitudeEffect => json!(vector(&row.magnitude_effect)?),
                    Output::CorrectedResiduals => json!(vector(&row.corrected_residuals)?),
                    Output::ProfiledResiduals => json!(vector(&row.profiled_residuals)?),
                    Output::Diagnostics => {
                        json!({"backward_residual":row.diagnostics.backward_residual,"estimated_forward_sensitivity":row.diagnostics.estimated_forward_sensitivity,"coefficient_backward_residual":row.diagnostics.coefficient_backward_residual,"coefficient_forward_sensitivity":row.diagnostics.coefficient_forward_sensitivity,"residual_l1":row.diagnostics.residual_l1,"solution_norm_inf":row.diagnostics.solution_norm_inf,"adjusted_residual_l1":row.diagnostics.adjusted_residual_l1,"adjusted_solution_norm_inf":row.diagnostics.adjusted_solution_norm_inf})
                    }
                };
                payload["value"] = value;
            }
            output[want.id()] = payload;
        }
        outputs.push(output);
    }
    Ok((
        json!({"kind":if checks.iter().all(|(_,p)|*p){"finite"}else{"failure"},"status":status,"numerical_status":numerical,"evaluations":outputs,"work":{"callbacks":callbacks,"segment_visits":segments},"requested_outputs":requested,"interpretation":"unqualified"}),
        checks,
    ))
}
