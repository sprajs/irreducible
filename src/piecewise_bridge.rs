//! Fixed compiled piecewise model transport; no Rust physical equations.
use super::generated::*;
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};
use std::{ffi::c_void, ptr};
#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct Parameters {
    pub h0_km_s_mpc: f64,
    pub q: [f64; 5],
}
#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct Query {
    pub z_expansion: f64,
    pub z_observer: f64,
    pub convention: String,
}
#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct Policy {
    pub maximum_parameters: u64,
    pub maximum_queries: u64,
    pub maximum_slots: u64,
    pub maximum_native_output_bytes: u64,
    pub maximum_total_segment_visits: u64,
}
struct Owned(*mut c_void);
impl Drop for Owned {
    fn drop(&mut self) {
        unsafe {
            cosmo_piecewise_result_destroy(self.0);
        }
    }
}
unsafe fn text(v: &Bytes) -> Result<String, String> {
    if v.length > 4096 || (v.length > 0 && v.data.is_null()) {
        return Err("INVALID_PIECEWISE_TEXT".into());
    }
    let bytes = if v.length == 0 {
        &[][..]
    } else {
        unsafe { std::slice::from_raw_parts(v.data, v.length as usize) }
    };
    String::from_utf8(bytes.to_vec()).map_err(|_| "INVALID_PIECEWISE_TEXT".into())
}
pub(crate) fn evaluate(
    parameters: &[Parameters],
    queries: &[Query],
    policy: &Policy,
) -> Result<Value, String> {
    let total = parameters
        .len()
        .checked_mul(queries.len())
        .ok_or("INVALID_BATCH_SHAPE")?;
    if parameters.len() > 4096 || queries.len() > 4096 || total > 65536 {
        return Err("RESOURCE_LIMIT".into());
    }
    let params: Vec<_> = parameters
        .iter()
        .map(|p| PiecewiseParameters {
            h0_km_s_mpc: p.h0_km_s_mpc,
            q: p.q,
        })
        .collect();
    let wire_queries: Vec<_> = queries
        .iter()
        .map(|q| {
            Ok(BackgroundQuery {
                z_expansion: q.z_expansion,
                z_observer: q.z_observer,
                convention: background_tag_id("convention", &q.convention)
                    .ok_or("UNKNOWN_BACKGROUND_CONVENTION")?,
                reserved: 0,
            })
        })
        .collect::<Result<_, String>>()?;
    let batch = PiecewiseBatch {
        struct_size: std::mem::size_of::<PiecewiseBatch>() as u32,
        abi_version: ABI_VERSION,
        parameters: params.as_ptr(),
        parameter_count: params.len() as u64,
        parameter_bytes: std::mem::size_of_val(params.as_slice()) as u64,
        queries: wire_queries.as_ptr(),
        query_count: wire_queries.len() as u64,
        query_bytes: std::mem::size_of_val(wire_queries.as_slice()) as u64,
    };
    let p = PiecewisePolicy {
        struct_size: std::mem::size_of::<PiecewisePolicy>() as u32,
        abi_version: ABI_VERSION,
        reserved: 0,
        maximum_parameters: policy.maximum_parameters,
        maximum_queries: policy.maximum_queries,
        maximum_slots: policy.maximum_slots,
        maximum_native_output_bytes: policy.maximum_native_output_bytes,
        maximum_total_segment_visits: policy.maximum_total_segment_visits,
    };
    let mut raw = ptr::null_mut();
    let code = unsafe { cosmo_piecewise_evaluate(&batch, &p, &mut raw) };
    if code != OK {
        return Err(format!("CORE_STATUS_{code}"));
    }
    if raw.is_null() {
        return Err("NULL_PIECEWISE_RESULT".into());
    }
    let owner = Owned(raw);
    let mut slots = ptr::null();
    let mut count = 0;
    let mut status = INVALID_INPUT;
    let mut numerical = 0;
    let mut segments = 0;
    if unsafe {
        cosmo_piecewise_result_view(
            owner.0,
            &mut slots,
            &mut count,
            &mut status,
            &mut numerical,
            &mut segments,
        )
    } != OK
        || count as usize > total
        || (count > 0
            && (slots.is_null() || (slots as usize) % std::mem::align_of::<PiecewiseSlot>() != 0))
    {
        return Err("INVALID_PIECEWISE_RESULT_VIEW".into());
    }
    if count as usize != total && !(count == 0 && status == BACKGROUND_STATUS_WORK_LIMIT) {
        return Err("INVALID_PIECEWISE_SLOT_COUNT".into());
    }
    if segments > policy.maximum_total_segment_visits {
        return Err("INVALID_PIECEWISE_SEGMENT_ACCOUNTING".into());
    }
    let rows = if count == 0 {
        &[][..]
    } else {
        unsafe { std::slice::from_raw_parts(slots, count as usize) }
    };
    let mut output = Vec::with_capacity(rows.len());
    for (i, s) in rows.iter().enumerate() {
        let pi = i / queries.len();
        let qi = i % queries.len();
        if s.struct_size as usize != std::mem::size_of::<PiecewiseSlot>()
            || s.abi_version != ABI_VERSION
            || s.parameter_index as usize != pi
            || s.query_index as usize != qi
            || s.parameters.h0_km_s_mpc.to_bits() != parameters[pi].h0_km_s_mpc.to_bits()
            || s.parameters
                .q
                .iter()
                .zip(parameters[pi].q)
                .any(|(a, b)| a.to_bits() != b.to_bits())
            || s.query.z_expansion.to_bits() != queries[qi].z_expansion.to_bits()
            || s.query.z_observer.to_bits() != queries[qi].z_observer.to_bits()
            || s.query.convention != wire_queries[qi].convention
            || [s.has_geometry, s.has_q, s.has_jerk, s.has_q0]
                .iter()
                .any(|f| *f > 1)
            || s.q_convention > 4
            || s.jerk_availability > 3
        {
            return Err("INVALID_PIECEWISE_SLOT".into());
        }
        let mut row = json!({"identity":{"parameter_index":pi,"query_index":qi,"parameters":parameters[pi],"query":queries[qi],"model_id":unsafe{text(&s.model_id)?},"radial_equation_id":unsafe{text(&s.radial_equation_id)?},"luminosity_equation_id":unsafe{text(&s.luminosity_equation_id)?},"shape_equation_id":unsafe{text(&s.shape_equation_id)?},"constant_set_id":unsafe{text(&s.constant_set_id)?}},"status":background_status_name(s.status).ok_or("INVALID_PIECEWISE_STATUS")?,"numerical_status":numerical_status_name(s.numerical_status).ok_or("INVALID_PIECEWISE_NUMERICAL_STATUS")?,"segments_processed":s.segments_processed,"derivatives":{"q_convention":piecewise_tag_name("q_convention",s.q_convention).ok_or("INVALID_PIECEWISE_Q_CONVENTION")?,"jerk_availability":piecewise_tag_name("jerk_availability",s.jerk_availability).ok_or("INVALID_PIECEWISE_JERK_AVAILABILITY")?,"has_q":s.has_q==1,"has_jerk":s.has_jerk==1,"has_q0":s.has_q0==1}});
        if s.status == BACKGROUND_STATUS_OK {
            let values = [
                s.expansion_E,
                s.h_km_s_mpc,
                s.radial_integral,
                s.radial_mpc,
                s.transverse_mpc,
                s.angular_diameter_mpc,
                s.luminosity_mpc,
                s.dimensionless_luminosity_shape,
                s.lookback_seconds,
                s.volume_mpc3_per_sr_per_redshift,
            ];
            if s.has_geometry != 1
                || s.numerical_status != 0
                || values.iter().any(|x| !x.is_finite())
            {
                return Err("INVALID_PIECEWISE_FINITE_PAYLOAD".into());
            }
            row["geometry"] = json!({"expansion_E":s.expansion_E,"h_km_s_mpc":s.h_km_s_mpc,"radial_integral":s.radial_integral,"radial_mpc":s.radial_mpc,"transverse_mpc":s.transverse_mpc,"angular_diameter_mpc":s.angular_diameter_mpc,"luminosity_mpc":s.luminosity_mpc,"dimensionless_luminosity_shape":s.dimensionless_luminosity_shape,"lookback_seconds":s.lookback_seconds,"volume_mpc3_per_sr_per_redshift":s.volume_mpc3_per_sr_per_redshift});
            for (name, flag, value) in [
                ("assigned_q", s.has_q, s.assigned_q),
                ("jerk", s.has_jerk, s.jerk),
                (
                    "q0_within_piecewise_model",
                    s.has_q0,
                    s.q0_within_piecewise_model,
                ),
            ] {
                if flag == 1 {
                    if !value.is_finite() {
                        return Err("INVALID_PIECEWISE_DERIVATIVE".into());
                    }
                    row["derivatives"][name] = json!(value);
                }
            }
        } else if s.has_geometry + s.has_q + s.has_jerk + s.has_q0 != 0
            || s.q_convention != 0
            || s.jerk_availability != 0
        {
            return Err("INVALID_PIECEWISE_FAILURE_PAYLOAD".into());
        }
        output.push(row);
    }
    let row_segments = rows
        .iter()
        .try_fold(0u64, |sum, row| sum.checked_add(row.segments_processed))
        .ok_or("INVALID_PIECEWISE_SEGMENT_ACCOUNTING")?;
    if row_segments != segments {
        return Err("INVALID_PIECEWISE_SEGMENT_ACCOUNTING".into());
    }
    Ok(
        json!({"kind":if status==BACKGROUND_STATUS_OK&&output.iter().all(|s|s["status"]=="ok"){"finite"}else{"failure"},"batch_status":background_status_name(status).ok_or("INVALID_PIECEWISE_BATCH_STATUS")?,"numerical_status":numerical_status_name(numerical).ok_or("INVALID_PIECEWISE_BATCH_NUMERICAL_STATUS")?,"segments_processed":segments,"slots":output}),
    )
}
