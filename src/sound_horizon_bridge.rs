//! A single synchronous native batch; Rust owns descriptors and records only.
use super::{generated::*, observations::copied};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::{ffi::c_void, mem::size_of, ptr};
#[derive(Clone, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct Point {
    pub h0_km_s_mpc: f64,
    pub omega_m: f64,
    pub omega_r: f64,
    pub omega_b: f64,
    pub omega_gamma: f64,
    pub z_drag: f64,
    pub drag_origin: String,
}
#[derive(Clone, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct Policy {
    pub absolute_tolerance_mpc: f64,
    pub relative_tolerance: f64,
    pub maximum_callbacks_per_point: u64,
    pub maximum_depth: u64,
    pub maximum_points: u64,
    pub maximum_total_callbacks: u64,
    pub maximum_native_bytes: u64,
}
struct Owned(*mut c_void);
impl Drop for Owned {
    fn drop(&mut self) {
        unsafe {
            irred_sound_horizon_result_destroy(self.0);
        }
    }
}
pub(crate) fn evaluate(points: &[Point], policy: &Policy) -> Result<(Value, bool), String> {
    if points.len() > 65536 || policy.maximum_native_bytes > 1024 * 1024 * 1024 {
        return Err("HARD_RESOURCE_LIMIT".into());
    }
    let inputs: Vec<_> = points
        .iter()
        .map(|p| SoundHorizonInput {
            struct_size: size_of::<SoundHorizonInput>() as u32,
            abi_version: ABI_VERSION,
            h0_km_s_mpc: p.h0_km_s_mpc,
            omega_m: p.omega_m,
            omega_r: p.omega_r,
            omega_b: p.omega_b,
            omega_gamma: p.omega_gamma,
            z_drag: p.z_drag,
            drag_origin: Bytes {
                data: p.drag_origin.as_ptr(),
                length: p.drag_origin.len() as u64,
            },
        })
        .collect();
    let batch = SoundHorizonBatch {
        struct_size: size_of::<SoundHorizonBatch>() as u32,
        abi_version: ABI_VERSION,
        points: inputs.as_ptr(),
        count: inputs.len() as u64,
        byte_length: std::mem::size_of_val(inputs.as_slice()) as u64,
    };
    let p = SoundHorizonPolicy {
        struct_size: size_of::<SoundHorizonPolicy>() as u32,
        abi_version: ABI_VERSION,
        absolute_tolerance_mpc: policy.absolute_tolerance_mpc,
        relative_tolerance: policy.relative_tolerance,
        maximum_callbacks_per_point: policy.maximum_callbacks_per_point,
        maximum_depth: policy.maximum_depth,
        maximum_points: policy.maximum_points,
        maximum_total_callbacks: policy.maximum_total_callbacks,
        maximum_native_bytes: policy.maximum_native_bytes,
    };
    let mut raw = ptr::null_mut();
    let status = unsafe { irred_sound_horizon_evaluate(&batch, &p, &mut raw) };
    if status != OK {
        return Err(format!("CORE_STATUS_{status}"));
    }
    if raw.is_null() {
        return Err("NULL_RESULT".into());
    }
    let owner = Owned(raw);
    let mut view: SoundHorizonResultView = unsafe { std::mem::zeroed() };
    let status = unsafe { irred_sound_horizon_result_view(owner.0, &mut view) };
    if status != OK
        || view.struct_size != size_of::<SoundHorizonResultView>() as u32
        || view.abi_version != ABI_VERSION
        || view.reserved != 0
        || (view.numerical_status == 0 && view.count != points.len() as u64)
        || (view.numerical_status != 0 && view.count != 0)
        || view.count > points.len() as u64
        || (view.count > 0
            && (view.rows.is_null()
                || (view.rows as usize) % std::mem::align_of::<SoundHorizonRow>() != 0))
    {
        return Err("INVALID_RESULT_VIEW".into());
    }
    let text = |v: Bytes| -> Result<String, String> {
        String::from_utf8(copied(v.data, v.length, 256)?)
            .map_err(|_| "INVALID_NATIVE_METADATA".into())
    };
    let metadata = json!({"model_id":text(view.model_id)?,"equation_id":text(view.equation_id)?,
        "constants_id":text(view.constant_set_id)?,"coordinate":text(view.coordinate_id)?,"arithmetic_id":text(view.arithmetic_id)?});
    let mut passed = view.numerical_status == 0;
    let mut rows = Vec::new();
    if view.count > 0 {
        for (index, r) in unsafe { std::slice::from_raw_parts(view.rows, view.count as usize) }
            .iter()
            .enumerate()
        {
            if r.struct_size != size_of::<SoundHorizonRow>() as u32
                || r.abi_version != ABI_VERSION
                || r.source.struct_size != size_of::<SoundHorizonInput>() as u32
                || r.source.abi_version != ABI_VERSION
                || r.has_value > 1
                || r.reserved != 0
                || r.reserved2 != 0
                || r.has_value == 1
                    && (r.numerical_status != 0
                        || !r.sound_horizon_mpc.is_finite()
                        || !r.error_estimate_mpc.is_finite())
            {
                return Err("INVALID_RESULT_VIEW".into());
            }
            let original = &points[index];
            let source_numbers = [
                r.source.h0_km_s_mpc,
                r.source.omega_m,
                r.source.omega_r,
                r.source.omega_b,
                r.source.omega_gamma,
                r.source.z_drag,
            ];
            let input_numbers = [
                original.h0_km_s_mpc,
                original.omega_m,
                original.omega_r,
                original.omega_b,
                original.omega_gamma,
                original.z_drag,
            ];
            let origin = String::from_utf8(copied(
                r.source.drag_origin.data,
                r.source.drag_origin.length,
                original.drag_origin.len(),
            )?)
            .map_err(|_| "INVALID_NATIVE_SOURCE")?;
            if source_numbers
                .iter()
                .zip(input_numbers)
                .any(|(a, b)| a.to_bits() != b.to_bits())
                || origin != original.drag_origin
            {
                return Err("INVALID_NATIVE_SOURCE".into());
            }
            let source = json!({"h0_km_s_mpc":r.source.h0_km_s_mpc,"omega_m":r.source.omega_m,"omega_r":r.source.omega_r,"omega_b":r.source.omega_b,"omega_gamma":r.source.omega_gamma,"z_drag":r.source.z_drag,"drag_origin":origin});
            let good = r.numerical_status == 0 && r.has_value == 1;
            passed &= good;
            rows.push(json!({"source":source,"numerical_status":numerical_status_name(r.numerical_status).ok_or("UNKNOWN_NUMERICAL_STATUS")?,
   "sound_horizon_mpc":if good{Some(r.sound_horizon_mpc)}else{None},
   "error_estimate_mpc":if good{Some(r.error_estimate_mpc)}else{None},"callbacks":r.callbacks}));
        }
    }
    Ok((
        json!({"metadata":metadata,"numerical_status":numerical_status_name(view.numerical_status).ok_or("UNKNOWN_NUMERICAL_STATUS")?,"rows":rows,"callbacks":view.callbacks}),
        passed,
    ))
}
