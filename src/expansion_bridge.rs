//! One coarse current expansion batch; native typed groups own all equations.
use super::{
    generated::*,
    observations::{copied, doubles},
};
use crate::model_spec::{
    BackgroundOutput, BackgroundRequest, ExpansionSpec, Geometry, ObserverConvention,
};
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};
use std::{ffi::c_void, mem::size_of, ptr};

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct Integration {
    pub absolute_tolerance: f64,
    pub relative_tolerance: f64,
    pub maximum_evaluations: u64,
    pub maximum_depth: u32,
}
#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct Policy {
    pub integration: Option<Integration>,
    pub maximum_models: u64,
    pub maximum_queries: u64,
    pub maximum_slots: u64,
    pub maximum_callbacks: u64,
    pub maximum_segment_visits: u64,
    pub maximum_native_bytes: u64,
}
impl Policy {
    pub(super) fn descriptor(&self) -> ExpansionPolicy {
        ExpansionPolicy {
            struct_size: size_of::<ExpansionPolicy>() as u32,
            abi_version: ABI_VERSION,
            has_integration: self.integration.is_some() as u32,
            max_depth: self.integration.as_ref().map_or(0, |p| p.maximum_depth),
            absolute_tolerance: self
                .integration
                .as_ref()
                .map_or(0., |p| p.absolute_tolerance),
            relative_tolerance: self
                .integration
                .as_ref()
                .map_or(0., |p| p.relative_tolerance),
            integration_max_evaluations: self
                .integration
                .as_ref()
                .map_or(0, |p| p.maximum_evaluations),
            maximum_models: self.maximum_models,
            maximum_queries: self.maximum_queries,
            maximum_slots: self.maximum_slots,
            maximum_callbacks: self.maximum_callbacks,
            maximum_segment_visits: self.maximum_segment_visits,
            maximum_native_bytes: self.maximum_native_bytes,
        }
    }
}
impl ExpansionSpec {
    pub(crate) fn active_parameters(&self) -> Vec<f64> {
        match self {
            Self::Lcdm { omega_m } => vec![*omega_m],
            Self::ConstantQ { q } => vec![*q],
            Self::Cpl { omega_m, w0, wa } => vec![*omega_m, *w0, *wa],
            Self::FixedQ5 { q } => q.to_vec(),
        }
    }
    pub(crate) fn tag(&self) -> u32 {
        match self {
            Self::Lcdm { .. } => 0,
            Self::ConstantQ { .. } => 1,
            Self::Cpl { .. } => 2,
            Self::FixedQ5 { .. } => 3,
        }
    }
    pub(super) fn descriptor(&self, parameters: &[f64]) -> super::generated::ExpansionSpec {
        super::generated::ExpansionSpec {
            struct_size: size_of::<super::generated::ExpansionSpec>() as u32,
            abi_version: ABI_VERSION,
            model: self.tag(),
            reserved: 0,
            parameters: doubles(parameters),
        }
    }
}
impl BackgroundOutput {
    pub(crate) fn bit(&self) -> u32 {
        match self {
            Self::Expansion => 32,
            Self::Radial => 1,
            Self::LuminosityShape => 2,
            Self::Clock => 4,
            Self::Physical => 8,
            Self::Kinematics => 16,
        }
    }
    pub(crate) fn id(&self) -> &'static str {
        match self {
            Self::Expansion => "expansion",
            Self::Radial => "radial",
            Self::LuminosityShape => "luminosity_shape",
            Self::Clock => "clock",
            Self::Physical => "physical",
            Self::Kinematics => "kinematics",
        }
    }
}
impl BackgroundRequest {
    fn descriptor(&self) -> Result<ExpansionRequest, String> {
        let mut requested = 0;
        for output in &self.requested_outputs {
            if requested & output.bit() != 0 {
                return Err("DUPLICATE_REQUESTED_OUTPUT".into());
            }
            requested |= output.bit();
        }
        if requested == 0 {
            return Err("EMPTY_REQUESTED_OUTPUTS".into());
        }
        Ok(ExpansionRequest {
            struct_size: size_of::<ExpansionRequest>() as u32,
            abi_version: ABI_VERSION,
            requested,
            presence_flags: self.observer.is_some() as u32
                | ((self.physical_scale.is_some() as u32) << 1),
            z_expansion: self.z_expansion,
            observer_redshift: self.observer.map_or(0., |o| o.redshift),
            observer_convention: self.observer.map_or(0, |o| match o.convention {
                ObserverConvention::GeometricSameRedshift => 0,
                ObserverConvention::ReleasedZhdZhel => 1,
            }),
            reserved: 0,
            h0_km_s_mpc: self.physical_scale.map_or(0., |p| p.h0_km_s_mpc),
        })
    }
}
struct Owner(*mut c_void);
impl Drop for Owner {
    fn drop(&mut self) {
        unsafe {
            cosmo_expansion_result_destroy(self.0);
        }
    }
}
fn state(value: &OutputState) -> Result<Value, String> {
    let availability =
        expansion_tag_name("availability", value.availability).ok_or("INVALID_AVAILABILITY")?;
    let status = expansion_tag_name("status", value.status).ok_or("INVALID_EXPANSION_STATUS")?;
    let numerical =
        numerical_status_name(value.numerical_status).ok_or("INVALID_NUMERICAL_STATUS")?;
    if value.reserved != 0 {
        return Err("INVALID_RESULT_RESERVED".into());
    }
    Ok(json!({"availability":availability,"status":status,"numerical_status":numerical}))
}
fn passed(value: &OutputState) -> bool {
    value.availability == 1 && value.status == 0 && value.numerical_status == 0
}
fn group(value: &OutputState, payload: Value) -> Result<Value, String> {
    let mut output = state(value)?;
    output["value"] = if passed(value) { payload } else { Value::Null };
    Ok(output)
}
fn scalar(value: &ScalarOutcome) -> Result<Value, String> {
    group(&value.state, json!(value.value))
}
pub(crate) struct Calculation {
    pub output: Value,
    pub checks: Vec<(&'static str, bool)>,
}
pub(crate) fn evaluate(
    models: &[ExpansionSpec],
    geometry: Geometry,
    queries: &[BackgroundRequest],
    policy: &Policy,
) -> Result<Calculation, String> {
    if models.len() > 65536 || queries.len() > 65536 {
        return Err("HARD_BATCH_LIMIT".into());
    }
    let active: Vec<_> = models
        .iter()
        .map(ExpansionSpec::active_parameters)
        .collect();
    let descriptors: Vec<_> = models
        .iter()
        .zip(&active)
        .map(|(m, p)| m.descriptor(p))
        .collect();
    let requests: Vec<_> = queries
        .iter()
        .map(BackgroundRequest::descriptor)
        .collect::<Result<_, _>>()?;
    let batch = ExpansionBatch {
        struct_size: size_of::<ExpansionBatch>() as u32,
        abi_version: ABI_VERSION,
        geometry: match geometry {
            Geometry::FlatFlrw => 0,
        },
        reserved: 0,
        models: descriptors.as_ptr(),
        model_count: descriptors.len() as u64,
        model_byte_length: std::mem::size_of_val(descriptors.as_slice()) as u64,
        queries: requests.as_ptr(),
        query_count: requests.len() as u64,
        query_byte_length: std::mem::size_of_val(requests.as_slice()) as u64,
    };
    let mut raw = ptr::null_mut();
    let transport = unsafe { cosmo_expansion_evaluate(&batch, &policy.descriptor(), &mut raw) };
    let owner = Owner(raw);
    if transport != OK {
        return Err(format!("CORE_STATUS_{transport}"));
    }
    if owner.0.is_null() {
        return Err("NULL_RESULT".into());
    }
    let (mut rows, mut count, mut status, mut numerical, mut callbacks, mut segments) =
        (ptr::null(), 0, 0, 0, 0, 0);
    if unsafe {
        cosmo_expansion_result_view(
            owner.0,
            &mut rows,
            &mut count,
            &mut status,
            &mut numerical,
            &mut callbacks,
            &mut segments,
        )
    } != OK
    {
        return Err("INVALID_RESULT_VIEW".into());
    }
    let max = models
        .len()
        .checked_mul(queries.len())
        .ok_or("RESULT_LENGTH_OVERFLOW")?;
    let rows = unsafe_view(rows, count, max)?;
    let (mut model_views, mut model_count, mut available) = (ptr::null(), 0, 0);
    if unsafe {
        cosmo_expansion_result_models(owner.0, &mut model_views, &mut model_count, &mut available)
    } != OK
    {
        return Err("INVALID_MODEL_VIEW".into());
    }
    let model_views = unsafe_view(model_views, model_count, models.len())?;
    if available > 1
        || (available == 1 && model_count != models.len() as u64)
        || (available == 0 && model_count != 0)
    {
        return Err("INVALID_SOURCE_AVAILABILITY".into());
    }
    let (mut source_queries, mut source_query_count, mut query_available) = (ptr::null(), 0, 0);
    if unsafe {
        cosmo_expansion_result_queries(
            owner.0,
            &mut source_queries,
            &mut source_query_count,
            &mut query_available,
        )
    } != OK
    {
        return Err("INVALID_QUERY_SOURCE_VIEW".into());
    }
    let source_queries = unsafe_view(source_queries, source_query_count, requests.len())?;
    if query_available != available
        || (available == 1 && source_queries.len() != requests.len())
        || (available == 0 && !source_queries.is_empty())
    {
        return Err("INVALID_QUERY_SOURCE_AVAILABILITY".into());
    }
    for (actual, expected) in source_queries.iter().zip(&requests) {
        if actual.struct_size != expected.struct_size
            || actual.abi_version != expected.abi_version
            || actual.requested != expected.requested
            || actual.presence_flags != expected.presence_flags
            || actual.reserved != 0
            || actual.observer_convention != expected.observer_convention
            || actual.z_expansion.to_bits() != expected.z_expansion.to_bits()
            || actual.observer_redshift.to_bits() != expected.observer_redshift.to_bits()
            || actual.h0_km_s_mpc.to_bits() != expected.h0_km_s_mpc.to_bits()
        {
            return Err("QUERY_SOURCE_CHANGED".into());
        }
    }
    let mut evaluations = Vec::with_capacity(rows.len());
    let mut checks: Vec<(&'static str, bool)> = Vec::new();
    if !models.is_empty() {
        checks.push((
            "model_admission",
            status == 0
                && available == 1
                && model_views
                    .iter()
                    .all(|model| model.preparation_status == 0),
        ));
    }
    for request in queries {
        for output in &request.requested_outputs {
            if !checks.iter().any(|(id, _)| *id == output.id()) {
                checks.push((output.id(), status == 0 && available == 1));
            }
        }
    }
    let mut model_json = Vec::with_capacity(model_views.len());
    for (index, model) in model_views.iter().enumerate() {
        if model.model_index != index as u64
            || model.struct_size != size_of::<ExpansionModelView>() as u32
            || model.abi_version != ABI_VERSION
            || model.reserved != 0
        {
            return Err("INVALID_MODEL_VIEW".into());
        }
        let start = usize::try_from(model.row_offset).map_err(|_| "INVALID_ROW_RANGE")?;
        let length = usize::try_from(model.row_count).map_err(|_| "INVALID_ROW_RANGE")?;
        if start.checked_add(length).is_none_or(|n| n > rows.len()) {
            return Err("INVALID_ROW_RANGE".into());
        }
        if model.preparation_status != 0 || model.evaluation_status != 0 || length != queries.len()
        {
            for (_, check) in &mut checks {
                *check = false;
            }
        }
        let source = copied(
            model.source.parameters.data,
            model.source.parameters.length,
            5,
        )?;
        if source.len() != active[index].len()
            || source
                .iter()
                .zip(&active[index])
                .any(|(a, b)| a.to_bits() != b.to_bits())
        {
            return Err("SOURCE_PARAMETER_CHANGED".into());
        }
        model_json.push(json!({"model_index":index,"source":models[index],"active_parameters":source,
            "preparation_status":expansion_tag_name("status",model.preparation_status).ok_or("INVALID_EXPANSION_STATUS")?,
            "evaluation_status":expansion_tag_name("status",model.evaluation_status).ok_or("INVALID_EXPANSION_STATUS")?,
            "numerical_status":numerical_status_name(model.numerical_status).ok_or("INVALID_NUMERICAL_STATUS")?,
            "row_offset":model.row_offset,"row_count":model.row_count,"model_id":text(&model.model_id)?,
            "constants_id":text(&model.constants_id)?,"radial_equation_id":text(&model.radial_equation_id)?}));
    }
    for row in rows {
        let q = usize::try_from(row.query_index).map_err(|_| "INVALID_QUERY_INDEX")?;
        if row.model_index >= models.len() as u64
            || q >= queries.len()
            || row.struct_size != size_of::<ExpansionRow>() as u32
            || row.abi_version != ABI_VERSION
            || row.has_node > 1
        {
            return Err("INVALID_RESULT_ROW".into());
        }
        let request = &queries[q];
        let mut groups = serde_json::Map::new();
        for output in &request.requested_outputs {
            let (payload, ok) = match output {
                BackgroundOutput::Expansion => (
                    group(
                        &row.expansion.state,
                        json!({"E":row.expansion.E,"H_km_s_mpc":scalar(&row.expansion.H_km_s_mpc)?}),
                    )?,
                    passed(&row.expansion.state)
                        && (request.physical_scale.is_none()
                            || passed(&row.expansion.H_km_s_mpc.state)),
                ),
                BackgroundOutput::Radial => (
                    group(
                        &row.radial.state,
                        json!({"E":row.radial.E,"integral":row.radial.integral,"error_estimate":row.radial.error_estimate}),
                    )?,
                    passed(&row.radial.state),
                ),
                BackgroundOutput::LuminosityShape => (
                    scalar(&row.luminosity_shape)?,
                    passed(&row.luminosity_shape.state),
                ),
                BackgroundOutput::Clock => (
                    group(
                        &row.clock.state,
                        json!({"integral":row.clock.integral,"error_estimate":row.clock.error_estimate,"lookback_seconds":scalar(&row.clock.lookback_seconds)?}),
                    )?,
                    passed(&row.clock.state)
                        && (request.physical_scale.is_none()
                            || passed(&row.clock.lookback_seconds.state)),
                ),
                BackgroundOutput::Physical => (
                    group(
                        &row.physical.state,
                        json!({"radial_mpc":row.physical.radial_mpc,"transverse_mpc":row.physical.transverse_mpc,"angular_diameter_mpc":row.physical.angular_diameter_mpc,"luminosity_mpc":row.physical.luminosity_mpc,"volume_mpc3_per_sr_per_redshift":row.physical.volume_mpc3_per_sr_per_redshift}),
                    )?,
                    passed(&row.physical.state),
                ),
                BackgroundOutput::Kinematics => (
                    group(
                        &row.kinematics.state,
                        json!({"q":row.kinematics.q,"jerk":if row.kinematics.has_jerk==1{Some(row.kinematics.jerk)}else{None},"q0":if row.kinematics.has_q0==1{Some(row.kinematics.q0_within_piecewise_model)}else{None},"q_convention":expansion_tag_name("q_convention",row.kinematics.q_convention).ok_or("INVALID_Q_CONVENTION")?,"jerk_availability":expansion_tag_name("jerk_availability",row.kinematics.jerk_availability).ok_or("INVALID_JERK_AVAILABILITY")?,"bin":row.kinematics.bin}),
                    )?,
                    passed(&row.kinematics.state),
                ),
            };
            groups.insert(output.id().into(), payload);
            if !ok {
                if let Some((_, check)) = checks.iter_mut().find(|(id, _)| *id == output.id()) {
                    *check = false;
                }
            }
        }
        evaluations.push(json!({"model_index":row.model_index,"query_index":row.query_index,"source":request,"admission_status":expansion_tag_name("status",row.admission_status).ok_or("INVALID_ADMISSION_STATUS")?,"node_index":if row.has_node==1{Some(row.node_index)}else{None},"groups":groups}));
    }
    let (mut nodes, mut node_count) = (ptr::null(), 0);
    if unsafe { cosmo_expansion_result_nodes(owner.0, &mut nodes, &mut node_count) } != OK {
        return Err("INVALID_NODE_VIEW".into());
    }
    let nodes=unsafe_view(nodes,node_count,max)?.iter().map(|node|json!({"model_index":node.model_index,"z_expansion":node.z_expansion,"callbacks":node.callbacks,"segment_visits":node.segment_visits})).collect::<Vec<_>>();
    let all_passed = checks.iter().all(|(_, ok)| *ok) && status == 0;
    Ok(Calculation {
        output: json!({"kind":if all_passed{"finite"}else{"failure"},"status":expansion_tag_name("status",status).ok_or("INVALID_EXPANSION_STATUS")?,"numerical_status":numerical_status_name(numerical).ok_or("INVALID_NUMERICAL_STATUS")?,"source_available":available==1,"models":model_json,"queries":if available==1{Some(queries)}else{None},"evaluations":evaluations,"nodes":nodes,"work":{"callbacks":callbacks,"segment_visits":segments},"interpretation":"unqualified"}),
        checks,
    })
}
fn unsafe_view<'a, T>(pointer: *const T, length: u64, cap: usize) -> Result<&'a [T], String> {
    let n = usize::try_from(length).map_err(|_| "INVALID_VIEW_LENGTH")?;
    if n > cap
        || n > isize::MAX as usize / size_of::<T>()
        || (n > 0 && (pointer.is_null() || (pointer as usize) % std::mem::align_of::<T>() != 0))
    {
        return Err("INVALID_RESULT_VIEW".into());
    }
    Ok(if n == 0 {
        &[]
    } else {
        unsafe { std::slice::from_raw_parts(pointer, n) }
    })
}
fn text(value: &Bytes) -> Result<String, String> {
    let raw = unsafe_view(value.data, value.length, 16384)?;
    String::from_utf8(raw.to_vec()).map_err(|_| "INVALID_RESULT_UTF8".into())
}
