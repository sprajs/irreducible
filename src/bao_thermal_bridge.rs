//! One coarse borrowed C call and owned result; no Rust BAO equations.
use super::generated::*;
use crate::{
    bao_thermal_run::{Envelope, Request},
    outcome::{NumericalCheck, Outcome, OutputCheck},
};
use serde_json::{Value, json};
use std::{
    mem::{size_of, zeroed},
    ptr, slice,
};
struct Owned(*mut std::ffi::c_void);
impl Drop for Owned {
    fn drop(&mut self) {
        unsafe {
            irred_bao_thermal_result_destroy(self.0);
        }
    }
}
struct Wire {
    queries: [BaoQuery; 64],
    ids: [Bytes; 64],
    axes: [Bytes; 64],
    models: [BaoThermalModel; 16],
    species: [BaoThermalSpecies; 256],
}
pub(crate) fn wire_size() -> usize {
    size_of::<Wire>()
}
fn bytes(x: &str) -> Bytes {
    Bytes {
        data: x.as_ptr(),
        length: x.len() as u64,
    }
}
fn doubles(x: &[f64]) -> F64Buffer {
    F64Buffer {
        struct_size: size_of::<F64Buffer>() as u32,
        abi_version: ABI_VERSION,
        element_type: 2,
        reserved: 0,
        data: x.as_ptr(),
        length: x.len() as u64,
        byte_length: std::mem::size_of_val(x) as u64,
    }
}
fn strings(x: &[Bytes]) -> Strings {
    Strings {
        data: x.as_ptr(),
        length: x.len() as u64,
        byte_length: std::mem::size_of_val(x) as u64,
    }
}
fn view_text<'a>(x: Bytes, _owner: &'a Owned) -> Result<&'a str, String> {
    if x.length == 0 {
        return Ok("");
    }
    if x.length > 4096 || x.data.is_null() {
        return Err("THERMAL_BAO_TEXT_VIEW".into());
    }
    std::str::from_utf8(unsafe { slice::from_raw_parts(x.data, x.length as usize) })
        .map_err(|_| "THERMAL_BAO_TEXT_VIEW".into())
}
fn eq(x: f64, y: f64) -> bool {
    x.to_bits() == y.to_bits()
}
fn state(x: OutputState) -> Result<Value, String> {
    if x.reserved != 0 || x.availability > 3 || x.status > 5 || x.numerical_status > 8 {
        return Err("THERMAL_BAO_STATE_VIEW".into());
    }
    Ok(
        json!({"availability":x.availability,"status":x.status,"numerical_status":x.numerical_status}),
    )
}
fn numeric(x: F64Buffer, n: usize, available: bool) -> Result<Value, String> {
    if x.struct_size != size_of::<F64Buffer>() as u32
        || x.abi_version != ABI_VERSION
        || x.element_type != 2
        || x.reserved != 0
        || x.length != if available { n as u64 } else { 0 }
        || x.byte_length != x.length * 8
        || (x.length > 0
            && (x.data.is_null() || x.data as usize % std::mem::align_of::<f64>() != 0))
    {
        return Err("THERMAL_BAO_ARRAY_VIEW".into());
    }
    if !available {
        return Ok(Value::Null);
    }
    let a = unsafe { slice::from_raw_parts(x.data, n) };
    if a.iter().any(|x| !x.is_finite()) {
        return Err("THERMAL_BAO_NONFINITE_VIEW".into());
    }
    // Exactly one copy into final Value storage; no intermediate Vec<f64>.
    Ok(Value::Array(a.iter().map(|x| Value::from(*x)).collect()))
}
fn source_check(v: &BaoThermalView, r: &Request, owner: &Owned) -> Result<(), String> {
    let text = |x| view_text(x, owner);
    if v.source_available == 0 {
        return Ok(());
    }
    let a = &v.source;
    let o = &r.observation;
    let n = o.queries.len();
    if a.struct_size != size_of::<BaoThermalSource>() as u32
        || a.abi_version != ABI_VERSION
        || a.role != 1
        || a.covariance_unit != 0
        || a.query_count != n as u64
        || a.query_byte_length != (n * size_of::<BaoQuery>()) as u64
        || a.queries.is_null()
        || a.queries as usize % std::mem::align_of::<BaoQuery>() != 0
    {
        return Err("THERMAL_BAO_SOURCE_VIEW".into());
    }
    for (q, t) in unsafe { slice::from_raw_parts(a.queries, n) }
        .iter()
        .zip(o.queries.iter())
    {
        let tag = match &*t.observable {
            "DM_over_rs" => 0,
            "DH_over_rs" => 1,
            "DV_over_rs" => 2,
            _ => return Err("THERMAL_BAO_OBSERVABLE".into()),
        };
        if q.reserved != 0 || q.observable != tag || !eq(q.z, t.redshift) {
            return Err("THERMAL_BAO_SOURCE_ORDER".into());
        }
    }
    for (actual, expected) in [
        (a.ordered_ids, &*o.ordered_row_ids),
        (a.covariance_axis_ids, &*o.covariance_axis_ids),
    ] {
        if actual.length != n as u64
            || actual.byte_length != (n * size_of::<Bytes>()) as u64
            || actual.data.is_null()
            || actual.data as usize % std::mem::align_of::<Bytes>() != 0
        {
            return Err("THERMAL_BAO_ID_VIEW".into());
        }
        for (x, y) in unsafe { slice::from_raw_parts(actual.data, n) }
            .iter()
            .zip(expected)
        {
            if text(*x)? != &**y {
                return Err("THERMAL_BAO_SOURCE_ORDER".into());
            }
        }
    }
    // Verify raw native ownership against original input directly. No readback
    // array, clone or substitution into the scientific specification.
    for (actual, expected) in [
        (a.observed, &*o.observed_ratios),
        (a.covariance, &*o.covariance_row_major),
    ] {
        if actual.struct_size != size_of::<F64Buffer>() as u32
            || actual.abi_version != ABI_VERSION
            || actual.element_type != 2
            || actual.reserved != 0
            || actual.length != expected.len() as u64
            || actual.byte_length != (expected.len() * 8) as u64
            || actual.data.is_null()
            || actual.data as usize % 8 != 0
        {
            return Err("THERMAL_BAO_SOURCE_ARRAY".into());
        }
        if unsafe { slice::from_raw_parts(actual.data, expected.len()) }
            .iter()
            .zip(expected)
            .any(|(x, y)| !eq(*x, *y))
        {
            return Err("THERMAL_BAO_SOURCE_BITS".into());
        }
    }
    for (x, y) in [
        (a.table_identity, &*o.table_identity),
        (a.covariance_identity, &*o.covariance_identity),
        (a.ordering_provenance, &*o.ordering_provenance),
        (a.calibration_provenance, &*o.calibration_provenance),
        (a.dependence_provenance, &*o.dependence_provenance),
        (a.redshift_convention, &*o.redshift_convention),
    ] {
        if text(x)? != y {
            return Err("THERMAL_BAO_PROVENANCE_VIEW".into());
        }
    }
    Ok(())
}
pub(crate) fn evaluate(r: Box<Request>, env: Envelope) -> Result<Outcome, String> {
    let mask = r.mask()?;
    let a = &r.observation;
    let q = &r.resource_policy;
    let pp = &q.predictions;
    let tp = &pp.thermal;
    // Generated structs contain only trivial scalar/pointer members. Fixed
    // zeroed descriptor arrays allocate one Box; request buffers remain borrowed.
    let mut w: Box<Wire> = Box::new(unsafe { zeroed() });
    let mut sp = 0;
    for (i, x) in a.queries.iter().enumerate() {
        w.queries[i] = BaoQuery {
            z: x.redshift,
            observable: match &*x.observable {
                "DM_over_rs" => 0,
                "DH_over_rs" => 1,
                "DV_over_rs" => 2,
                _ => return Err("THERMAL_BAO_OBSERVABLE".into()),
            },
            reserved: 0,
        };
        w.ids[i] = bytes(&a.ordered_row_ids[i]);
        w.axes[i] = bytes(&a.covariance_axis_ids[i]);
    }
    for (i, x) in r.models.iter().enumerate() {
        let m = &x.physical_model;
        let start = sp;
        for x in m.species.iter() {
            w.species[sp] = BaoThermalSpecies {
                struct_size: size_of::<BaoThermalSpecies>() as u32,
                abi_version: ABI_VERSION,
                reserved0: 0,
                reserved1: 0,
                mass_ev: x.mass_ev,
                temperature_today_kelvin: x.temperature_today_kelvin,
                statistical_weight: x.statistical_weight,
            };
            sp += 1;
        }
        w.models[i] = BaoThermalModel {
            struct_size: size_of::<BaoThermalModel>() as u32,
            abi_version: ABI_VERSION,
            reserved0: 0,
            reserved1: 0,
            id: bytes(&x.id),
            h0_km_s_mpc: m.h0_km_s_mpc,
            physical_baryon_density: m.physical_baryon_density,
            physical_cdm_density: m.physical_cdm_density,
            tcmb_kelvin: m.tcmb_kelvin,
            physical_massless_nonphoton_density: m.physical_massless_nonphoton_density,
            species: if start == sp {
                ptr::null()
            } else {
                w.species[start..].as_ptr()
            },
            species_count: (sp - start) as u64,
            species_byte_length: ((sp - start) * size_of::<BaoThermalSpecies>()) as u64,
            z_drag: x.z_drag,
            drag_origin: bytes(&x.drag_origin),
            source_origin: bytes(&x.source_origin),
        };
    }
    let source = BaoThermalSource {
        struct_size: size_of::<BaoThermalSource>() as u32,
        abi_version: ABI_VERSION,
        role: 1,
        covariance_unit: 0,
        queries: w.queries.as_ptr(),
        query_count: a.queries.len() as u64,
        query_byte_length: (a.queries.len() * size_of::<BaoQuery>()) as u64,
        observed: doubles(&a.observed_ratios),
        covariance: doubles(&a.covariance_row_major),
        ordered_ids: strings(&w.ids[..a.queries.len()]),
        covariance_axis_ids: strings(&w.axes[..a.queries.len()]),
        table_identity: bytes(&a.table_identity),
        covariance_identity: bytes(&a.covariance_identity),
        ordering_provenance: bytes(&a.ordering_provenance),
        calibration_provenance: bytes(&a.calibration_provenance),
        dependence_provenance: bytes(&a.dependence_provenance),
        redshift_convention: bytes(&a.redshift_convention),
    };
    let batch = BaoThermalBatch {
        struct_size: size_of::<BaoThermalBatch>() as u32,
        abi_version: ABI_VERSION,
        models: w.models.as_ptr(),
        model_count: r.models.len() as u64,
        model_byte_length: (r.models.len() * size_of::<BaoThermalModel>()) as u64,
    };
    let policy = BaoThermalPolicy {
        struct_size: size_of::<BaoThermalPolicy>() as u32,
        abi_version: ABI_VERSION,
        arithmetic: if &*q.arithmetic == "wide" { 1 } else { 0 },
        requested: mask,
        maximum_rows: q.maximum_rows,
        maximum_matrix_elements: q.maximum_matrix_elements,
        maximum_models: q.maximum_models,
        maximum_species_per_model: q.maximum_species_per_model,
        maximum_string_bytes: q.maximum_string_bytes,
        maximum_output_array_elements: q.maximum_output_array_elements,
        maximum_native_bytes: q.maximum_native_bytes,
        maximum_preparation_native_bytes: q.maximum_preparation_native_bytes,
        maximum_evaluation_native_bytes: q.maximum_evaluation_native_bytes,
        maximum_total_callbacks: q.maximum_total_callbacks,
        maximum_forward_sensitivity: q.maximum_forward_sensitivity,
        maximum_projection_log_density_error: q.maximum_projection_log_density_error,
        predictions: BaoThermalPredictionPolicy {
            struct_size: size_of::<BaoThermalPredictionPolicy>() as u32,
            abi_version: ABI_VERSION,
            momentum_method: 0,
            maximum_depth: pp.maximum_depth,
            thermal_maximum_depth: tp.maximum_depth,
            reserved: 0,
            absolute_tolerance_mpc: pp.absolute_tolerance_mpc,
            relative_tolerance: pp.relative_tolerance,
            absolute_tolerance_ratio: pp.absolute_tolerance_ratio,
            relative_tolerance_ratio: pp.relative_tolerance_ratio,
            thermal_absolute_tolerance: tp.absolute_tolerance,
            thermal_relative_tolerance: tp.relative_tolerance,
            maximum_callbacks_per_point: pp.maximum_callbacks_per_point,
            maximum_total_callbacks: pp.maximum_total_callbacks,
            maximum_native_bytes: pp.maximum_native_bytes,
            thermal_maximum_callbacks_per_evaluation: tp.maximum_callbacks_per_evaluation,
            thermal_maximum_total_callbacks: tp.maximum_total_callbacks,
            thermal_maximum_native_bytes: tp.maximum_native_bytes,
        },
    };
    let mut raw = ptr::null_mut();
    let transport = unsafe { irred_bao_thermal_evaluate(&source, &batch, &policy, &mut raw) };
    if transport == BAO_THERMAL_QUOTA_REFUSED {
        if !raw.is_null() {
            let _owned = Owned(raw);
            return Err("THERMAL_BAO_QUOTA_OWNER".into());
        }
        drop(w);
        let spec = serde_json::to_value(&*r).map_err(|_| "THERMAL_BAO_SPEC")?;
        drop(r);
        return Ok(Outcome::scientific(
            spec,
            json!({"kind":"failure","error_id":"NUMERICAL_QUALIFICATION_REQUIRED",
            "transport_status":transport,"native_payload_absent":true,"source_prepare_call_attempted":false,"thermal_batch_call_attempted":false}),
            "thermal_bao_requested_outputs",
            json!({"requested":"retained thermal BAO density","actual":null}),
            json!({"actual":null}),
            json!({"rust_conditional_peak":env.peak}),
            "global cap below minimal owner; no native execution",
        ));
    }
    if transport != OK {
        if !raw.is_null() {
            let _owned = Owned(raw);
            return Err("THERMAL_BAO_TRANSPORT_OWNER".into());
        }
        return Err("THERMAL_BAO_NATIVE_TRANSPORT".into());
    }
    if raw.is_null() {
        return Err("THERMAL_BAO_NULL_RESULT".into());
    }
    let owned = Owned(raw);
    let text = |x| view_text(x, &owned);
    let mut v: BaoThermalView = unsafe { zeroed() };
    if unsafe { irred_bao_thermal_result_view(owned.0, &mut v) } != OK
        || v.struct_size != size_of::<BaoThermalView>() as u32
        || v.abi_version != ABI_VERSION
        || v.reserved != 0
        || v.phase > 3
        || v.source_available > 1
        || v.source_prepare_call_attempted > 1
        || v.source_factor_completed > 1
        || v.thermal_batch_call_attempted > 1
        || v.requested != mask
        || v.arithmetic_requested != policy.arithmetic
        || v.numerical_status > 8
        || v.density_status > 5
        || v.model_count > r.models.len() as u64
        || v.row_byte_length != v.model_count * size_of::<BaoThermalRow>() as u64
        || v.retained_payload_bytes > q.maximum_native_bytes
        || Some(v.callbacks) != v.outer_callbacks.checked_add(v.momentum_callbacks)
        || v.preparation_callbacks > v.momentum_callbacks
        || v.callbacks > q.maximum_total_callbacks
    {
        return Err("THERMAL_BAO_RESULT_VIEW".into());
    }
    if v.source_factor_completed != 0
        && (v.source_prepare_call_attempted == 0 || v.source_available == 0)
    {
        return Err("THERMAL_BAO_PHASE_VIEW".into());
    }
    if v.thermal_batch_call_attempted != 0 && v.source_factor_completed == 0 {
        return Err("THERMAL_BAO_PHASE_VIEW".into());
    }
    if v.source_factor_completed == 0 && v.actual_source_arithmetic_id.length != 0 {
        return Err("THERMAL_BAO_ARITHMETIC_VIEW".into());
    }
    if v.source_factor_completed != 0
        && text(v.actual_source_arithmetic_id)?
            != if policy.arithmetic == 1 {
                "F02/longdouble-cpu/v1"
            } else {
                "F02/binary64-legacy/v1"
            }
    {
        return Err("THERMAL_BAO_ARITHMETIC_VIEW".into());
    }
    source_check(&v, &r, &owned)?;
    let rows = if v.model_count == 0 {
        &[][..]
    } else {
        if v.rows.is_null() || v.rows as usize % std::mem::align_of::<BaoThermalRow>() != 0 {
            return Err("THERMAL_BAO_ROW_VIEW".into());
        }
        unsafe { slice::from_raw_parts(v.rows, v.model_count as usize) }
    };
    let mut outputs = Vec::with_capacity(rows.len());
    let mut passed = [mask & 1 == 0, mask & 2 == 0, mask & 4 == 0];
    if rows.len() == r.models.len() {
        for (i, bit) in [1, 2, 4].iter().enumerate() {
            if mask & bit != 0 {
                passed[i] = true;
            }
        }
    }
    for (i, x) in rows.iter().enumerate() {
        let original = &r.models[i];
        let m = &original.physical_model;
        let src = x.source;
        if x.struct_size != size_of::<BaoThermalRow>() as u32
            || x.abi_version != ABI_VERSION
            || x.model_index != i as u64
            || src.struct_size != size_of::<BaoThermalModel>() as u32
            || src.abi_version != ABI_VERSION
            || src.reserved0 != 0
            || src.reserved1 != 0
            || text(src.id)? != &*original.id
            || text(src.drag_origin)? != &*original.drag_origin
            || text(src.source_origin)? != &*original.source_origin
            || !eq(src.h0_km_s_mpc, m.h0_km_s_mpc)
            || !eq(src.physical_baryon_density, m.physical_baryon_density)
            || !eq(src.physical_cdm_density, m.physical_cdm_density)
            || !eq(src.tcmb_kelvin, m.tcmb_kelvin)
            || !eq(
                src.physical_massless_nonphoton_density,
                m.physical_massless_nonphoton_density,
            )
            || !eq(src.z_drag, original.z_drag)
            || src.species_count != m.species.len() as u64
            || src.species_byte_length != src.species_count * size_of::<BaoThermalSpecies>() as u64
            || Some(x.callbacks) != x.outer_callbacks.checked_add(x.momentum_callbacks)
            || x.preparation_callbacks > x.momentum_callbacks
        {
            return Err("THERMAL_BAO_MODEL_VIEW".into());
        }
        if src.species_count > 0 {
            if src.species.is_null()
                || src.species as usize % std::mem::align_of::<BaoThermalSpecies>() != 0
            {
                return Err("THERMAL_BAO_SPECIES_VIEW".into());
            }
            for (a, b) in unsafe { slice::from_raw_parts(src.species, m.species.len()) }
                .iter()
                .zip(m.species.iter())
            {
                if a.struct_size != size_of::<BaoThermalSpecies>() as u32
                    || a.abi_version != ABI_VERSION
                    || a.reserved0 != 0
                    || a.reserved1 != 0
                    || !eq(a.mass_ev, b.mass_ev)
                    || !eq(a.temperature_today_kelvin, b.temperature_today_kelvin)
                    || !eq(a.statistical_weight, b.statistical_weight)
                {
                    return Err("THERMAL_BAO_SPECIES_VIEW".into());
                }
            }
        }
        if x.density.struct_size != size_of::<BaoThermalDensity>() as u32
            || x.density.abi_version != ABI_VERSION
            || x.density.reserved != 0
            || x.density.has_projection_estimate > 1
            || x.preparation_status > 8
            || x.numerical_status > 8
        {
            return Err("THERMAL_BAO_DIAGNOSTIC_VIEW".into());
        }
        let states = [x.density.state, x.prediction_state, x.residual_state];
        for (j, bit) in [1, 2, 4].iter().enumerate() {
            state(states[j])?;
            if (mask & bit == 0) != (states[j].availability == 0) {
                return Err("THERMAL_BAO_OMISSION_VIEW".into());
            }
            if mask & bit != 0 {
                passed[j] &= states[j].availability == 1;
            }
        }
        let mut row = json!({"model_index":i,"model_id":text(src.id)?,"preparation_status":x.preparation_status,"numerical_status":x.numerical_status,
            "density_state":state(x.density.state)?,"prediction_state":state(x.prediction_state)?,"residual_state":state(x.residual_state)?,
            "callbacks":x.callbacks,"preparation_callbacks":x.preparation_callbacks,"outer_callbacks":x.outer_callbacks,"momentum_callbacks":x.momentum_callbacks});
        if mask & 2 != 0 {
            row["predictions"] = numeric(
                x.predictions,
                a.queries.len(),
                x.prediction_state.availability == 1,
            )?;
        } else {
            numeric(x.predictions, a.queries.len(), false)?;
        }
        if mask & 4 != 0 {
            row["residuals"] = numeric(
                x.residuals,
                a.queries.len(),
                x.residual_state.availability == 1,
            )?;
        } else {
            numeric(x.residuals, a.queries.len(), false)?;
        }
        let d = x.density;
        if d.has_projection_estimate != 0 {
            if !d.projection_log_density_error_estimate.is_finite()
                || d.projection_log_density_error_estimate < 0.
            {
                return Err("THERMAL_BAO_PROJECTION_VIEW".into());
            }
            row["projection_log_density_error_estimate"] =
                json!(d.projection_log_density_error_estimate);
        }
        if d.state.availability == 1 {
            if [
                d.quadratic,
                d.log_determinant,
                d.log_normalization,
                d.log_density,
                d.backward_residual,
                d.estimated_forward_sensitivity,
            ]
            .iter()
            .any(|x| !x.is_finite())
            {
                return Err("THERMAL_BAO_DENSITY_VIEW".into());
            }
            row["density"] = json!({"quadratic":d.quadratic,"log_determinant":d.log_determinant,"log_normalization":d.log_normalization,
                "log_density":d.log_density,"backward_residual":d.backward_residual,"estimated_forward_sensitivity":d.estimated_forward_sensitivity});
        }
        outputs.push(row);
    }
    let complete = passed.iter().all(|x| *x) && v.thermal_batch_call_attempted != 0;
    let output = json!({"kind":if complete{"finite"}else{"failure"},"phase":v.phase,"numerical_status":v.numerical_status,"density_status":v.density_status,
        "source_available":v.source_available!=0,"source_prepare_call_attempted":v.source_prepare_call_attempted!=0,
        "source_factor_completed":v.source_factor_completed!=0,"thermal_batch_call_attempted":v.thermal_batch_call_attempted!=0,"native_payload_absent":false,"models":outputs});
    let method = json!({"requested":text(v.thermal_density_id)?,"source_preparation_attempted":v.source_prepare_call_attempted!=0,
        "source_factor_completed":v.source_factor_completed!=0,"thermal_evaluation_attempted":v.thermal_batch_call_attempted!=0,
        "actual_thermal_equation":if v.thermal_batch_call_attempted!=0{Some(text(v.thermal_equation_id)?)}else{None},
        "physical_mapping":if v.thermal_batch_call_attempted!=0{Some(text(v.physical_mapping_id)?)}else{None}});
    let arithmetic = json!({"requested":&*q.arithmetic,"actual_source_factor":if v.source_factor_completed!=0{Some(text(v.actual_source_arithmetic_id)?)}else{None}});
    let resources = json!({"requested":q,"rust_operation_payload_cap":crate::bao_thermal_run::RUST_CAP,"conditional_decoder_peak":env.decoder,
        "conditional_bridge_peak":env.bridge,"conditional_completion_peak":env.completion,"conditional_recorder_peak":env.recorder,
        "native_retained_payload":v.retained_payload_bytes,"native_preparation_peak_bound":v.preparation_peak_bound_bytes,
        "native_evaluation_peak_bound":v.evaluation_peak_bound_bytes,"native_required_global_peak":v.required_global_peak_bytes,
        "callbacks":v.callbacks,"preparation_callbacks":v.preparation_callbacks,"outer_callbacks":v.outer_callbacks,"momentum_callbacks":v.momentum_callbacks,
        "common_framework_allocator_stack_rss_excluded":true});
    drop(w);
    let specification = serde_json::to_value(&*r).map_err(|_| "THERMAL_BAO_SPEC")?;
    drop(r);
    drop(owned);
    let mut outcome = Outcome::scientific(
        specification,
        output,
        "thermal_bao_requested_outputs",
        method,
        arithmetic,
        resources,
        "conditional supplied-drag native thermal BAO; no observational/inference qualification",
    );
    if !complete {
        outcome.output["error_id"] = json!("NUMERICAL_QUALIFICATION_REQUIRED");
    }
    outcome.outputs.clear();
    for (i, (bit, id)) in [
        (1, "normalized_density"),
        (2, "predictions"),
        (4, "residuals"),
    ]
    .iter()
    .enumerate()
    {
        if mask & bit != 0 {
            outcome.outputs.push(OutputCheck{check_kind:"numerical_contract",id,required:true,numerical:if passed[i]{NumericalCheck::ChecksPassed}else{NumericalCheck::Failed},
            inference:"not_applicable",interpretation:"unqualified",evidence:vec![],validation_coverage:"ordered native thermal group status; separate density projection"});
        }
    }
    Ok(outcome)
}
