// Shared structural transport implementation; no physical equations.
macro_rules! define_background_bridge {
    ($($extra:ident),*) => {
        use serde::{Deserialize, Serialize};
        use serde_json::{Value, json};
        use std::{ffi::c_void, ptr};
        #[derive(Deserialize, Serialize)]
        #[serde(deny_unknown_fields)]
        pub(crate) struct Parameters {
            pub model: String,
            pub h0_km_s_mpc: f64,
            pub omega_m: f64,
            pub constant_q: f64,
            $(pub $extra: f64,)*
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
            pub maximum_total_evaluations: u64,
            pub maximum_evaluations_per_integral: u64,
            pub maximum_depth: u32,
            pub absolute_tolerance: f64,
            pub relative_tolerance: f64,
        }
        struct ResultOwner(*mut c_void);
        impl Drop for ResultOwner {
            fn drop(&mut self) {
                unsafe {
                    destroy_ffi(self.0);
                }
            }
        }
        pub(crate) fn evaluate(
            parameters: &[Parameters],
            queries: &[Query],
            policy: &Policy,
        ) -> Result<Value, String> {
            let attempted_parameters = parameters;
            let total = parameters
                .len()
                .checked_mul(queries.len())
                .ok_or("INVALID_BATCH_SHAPE")?;
            // Bound JSON row construction as well as the separately declared native slot storage.
            if parameters.len() as u64 > policy.maximum_parameters
                || queries.len() as u64 > policy.maximum_queries
                || parameters.len() > 4096
                || queries.len() > 4096
                || total > 65536
                || total as u64 > policy.maximum_slots
                || total
                    .checked_mul(std::mem::size_of::<WireSlot>())
                    .ok_or("INVALID_BATCH_SHAPE")? as u64
                    > policy.maximum_native_output_bytes
            {
                return Err("RESOURCE_LIMIT".into());
            }
            let parameters: Result<Vec<_>, String> = parameters
                .iter()
                .map(|p| {
                    Ok(WireParameters {
                        model: model_tag(&p.model).ok_or("UNKNOWN_BACKGROUND_MODEL")?,
                        reserved: 0,
                        h0_km_s_mpc: p.h0_km_s_mpc,
                        omega_m: p.omega_m,
                        constant_q: p.constant_q,
                        $($extra: p.$extra,)*
                    })
                })
                .collect();
            let parameters = parameters?;
            let queries: Result<Vec<_>, String> = queries
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
                .collect();
            let queries = queries?;
            let descriptor = WireBatch {
                struct_size: std::mem::size_of::<WireBatch>() as u32,
                abi_version: ABI_VERSION,
                parameters: parameters.as_ptr(),
                parameter_count: parameters.len() as u64,
                parameter_byte_length: std::mem::size_of_val(parameters.as_slice()) as u64,
                queries: queries.as_ptr(),
                query_count: queries.len() as u64,
                query_byte_length: std::mem::size_of_val(queries.as_slice()) as u64,
            };
            let policy = BackgroundPolicy {
                struct_size: std::mem::size_of::<BackgroundPolicy>() as u32,
                abi_version: ABI_VERSION,
                maximum_depth: policy.maximum_depth,
                reserved: 0,
                maximum_parameters: policy.maximum_parameters,
                maximum_queries: policy.maximum_queries,
                maximum_slots: policy.maximum_slots,
                maximum_native_output_bytes: policy.maximum_native_output_bytes,
                maximum_total_evaluations: policy.maximum_total_evaluations,
                maximum_evaluations_per_integral: policy.maximum_evaluations_per_integral,
                absolute_tolerance: policy.absolute_tolerance,
                relative_tolerance: policy.relative_tolerance,
            };
            let mut raw = ptr::null_mut();
            let status = unsafe { evaluate_ffi(&descriptor, &policy, &mut raw) };
            let owner = ResultOwner(raw);
            if status != OK {
                return Err(format!("CORE_STATUS_{status}"));
            }
            if raw.is_null() {
                return Err("NULL_BACKGROUND_RESULT".into());
            }
            let mut slots = ptr::null();
            let mut length = 0;
            if unsafe { view_ffi(owner.0, &mut slots, &mut length) } != OK
                || length != total as u64
                || (length > 0
                    && (slots.is_null() || (slots as usize) % std::mem::align_of::<WireSlot>() != 0))
            {
                return Err("INVALID_RESULT_VIEW".into());
            }
            let slots = if length == 0 {
                &[][..]
            } else {
                unsafe { std::slice::from_raw_parts(slots, length as usize) }
            };
            let text = |value: &Bytes| -> Result<String, String> {
                if value.length > 256 || (value.length > 0 && value.data.is_null()) {
                    return Err("INVALID_EQUATION_VIEW".into());
                }
                if value.length == 0 {
                    return Ok(String::new());
                }
                String::from_utf8(
                    unsafe { std::slice::from_raw_parts(value.data, value.length as usize) }.to_vec(),
                )
                .map_err(|_| "INVALID_EQUATION_ENCODING".into())
            };
            let mut evaluations = 0u64;
            let mut failed = false;
            let mut rows = Vec::with_capacity(total);
            for (index, slot) in slots.iter().enumerate() {
                if slot.parameter_index as usize != index / queries.len()
                    || slot.query_index as usize != index % queries.len()
                {
                    return Err("INVALID_RESULT_ORDER".into());
                }
                let p = &parameters[slot.parameter_index as usize];
                let q = &queries[slot.query_index as usize];
                if slot.parameters.model != p.model
                    || slot.parameters.reserved != 0
                    || slot.parameters.h0_km_s_mpc.to_bits() != p.h0_km_s_mpc.to_bits()
                    || slot.parameters.omega_m.to_bits() != p.omega_m.to_bits()
                    || slot.parameters.constant_q.to_bits() != p.constant_q.to_bits()
                    $(|| slot.parameters.$extra.to_bits() != p.$extra.to_bits())*
                    || slot.query.z_expansion.to_bits() != q.z_expansion.to_bits()
                    || slot.query.z_observer.to_bits() != q.z_observer.to_bits()
                    || slot.query.convention != q.convention
                    || slot.query.reserved != 0
                {
                    return Err("INVALID_SOURCE_VIEW".into());
                }
                evaluations = evaluations
                    .checked_add(slot.evaluations)
                    .ok_or("INVALID_EVALUATION_COUNT")?;
                let label = background_status_name(slot.status).ok_or("INVALID_SCIENTIFIC_STATUS")?;
                let numerical =
                    numerical_status_name(slot.numerical_status).ok_or("INVALID_NUMERICAL_STATUS")?;
                let diagnostics = json!({"evaluations":slot.evaluations,"numerical_status":numerical,"error_estimate_kind":"empirical quadrature diagnostic; not certified bound"});
                let mut identity = json!({"parameter_index":slot.parameter_index,"query_index":slot.query_index,"model_id":text(&slot.model_id)?,"constants_id":text(&slot.constants_id)?,"radial_equation_id":text(&slot.radial_equation_id)?,"luminosity_equation_id":text(&slot.luminosity_equation_id)?,"shape_equation_id":text(&slot.shape_equation_id)?});
                identity["source_parameters"] = serde_json::to_value(&attempted_parameters[slot.parameter_index as usize]).map_err(|e|e.to_string())?;
                if slot.status != BACKGROUND_STATUS_OK {
                    failed = true;
                    rows.push(json!({"identity":identity,"result":{"kind":"failure","status":label},"diagnostics":diagnostics}));
                    continue;
                }
                let physical = [
                    slot.expansion_e,
                    slot.h_km_s_mpc,
                    slot.radial_integral,
                    slot.radial_mpc,
                    slot.transverse_mpc,
                    slot.angular_diameter_mpc,
                    slot.luminosity_mpc,
                    slot.dimensionless_luminosity_shape,
                    slot.lookback_seconds,
                    slot.volume_mpc3_per_sr_per_redshift,
                    slot.deceleration_q,
                    slot.jerk,
                    slot.radial_integral_error,
                    slot.lookback_integral_error,
                ];
                if slot.numerical_status != 0 || physical.iter().any(|x| !x.is_finite()) {
                    return Err("INVALID_FINITE_PAYLOAD".into());
                }
                rows.push(json!({"identity":identity,"result":{"kind":"finite","expansion_E":slot.expansion_e,"h_km_s_mpc":slot.h_km_s_mpc,"radial_integral":slot.radial_integral,"radial_mpc":slot.radial_mpc,"transverse_mpc":slot.transverse_mpc,"angular_diameter_mpc":slot.angular_diameter_mpc,"luminosity_mpc":slot.luminosity_mpc,"dimensionless_luminosity_shape":slot.dimensionless_luminosity_shape,"lookback_seconds":slot.lookback_seconds,"volume_mpc3_per_sr_per_redshift":slot.volume_mpc3_per_sr_per_redshift,"deceleration_q":slot.deceleration_q,"jerk":slot.jerk},"diagnostics":{"evaluations":slot.evaluations,"numerical_status":numerical,"radial_integral_error":slot.radial_integral_error,"lookback_integral_error":slot.lookback_integral_error,"error_estimate_kind":"empirical quadrature diagnostic; not certified bound"}}));
            }
            if evaluations > policy.maximum_total_evaluations {
                return Err("INVALID_GLOBAL_WORK_ACCOUNTING".into());
            }
            Ok(
                json!({"kind":if failed {"failure"}else{"finite"},"error_id":if failed {Some("BACKGROUND_EVALUATION_FAILURE")}else{None},"ordering":"parameter major, query minor","evaluations":evaluations,"native_output_bytes":total*std::mem::size_of::<WireSlot>(),"rows":rows}),
            )
        }
    };
}
