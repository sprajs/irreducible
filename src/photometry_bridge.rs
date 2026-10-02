//! One borrowed input batch; all radiometric equations remain native.
use super::{generated::*, observations::copied};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::{mem::size_of, ptr};
#[derive(Clone, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct Input {
    pub luminosity_watt_per_metre: f64,
    pub rest_lower_metre: f64,
    pub rest_upper_metre: f64,
    pub observed_lower_metre: f64,
    pub observed_upper_metre: f64,
    pub luminosity_distance_metre: f64,
    pub redshift: f64,
    pub collecting_area_square_metre: f64,
    pub optical_transmission: f64,
    pub observer_exposure_second: f64,
}
impl Input {
    fn descriptor(&self) -> PhotometryInput {
        PhotometryInput {
            struct_size: size_of::<PhotometryInput>() as u32,
            abi_version: ABI_VERSION,
            luminosity_watt_per_metre: self.luminosity_watt_per_metre,
            rest_lower_metre: self.rest_lower_metre,
            rest_upper_metre: self.rest_upper_metre,
            observed_lower_metre: self.observed_lower_metre,
            observed_upper_metre: self.observed_upper_metre,
            luminosity_distance_metre: self.luminosity_distance_metre,
            redshift: self.redshift,
            collecting_area_square_metre: self.collecting_area_square_metre,
            optical_transmission: self.optical_transmission,
            observer_exposure_second: self.observer_exposure_second,
        }
    }
}
#[derive(Clone, Copy, Deserialize, Serialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum Output {
    IncidentBandFlux,
    CollectedEnergy,
    ExpectedTransmittedPhotons,
}
impl Output {
    pub fn id(self) -> &'static str {
        match self {
            Self::IncidentBandFlux => "incident_band_flux",
            Self::CollectedEnergy => "collected_energy",
            Self::ExpectedTransmittedPhotons => "expected_transmitted_photons",
        }
    }
    pub(super) fn bit(self) -> u32 {
        match self {
            Self::IncidentBandFlux => PHOTOMETRY_OUTPUT_INCIDENT_BAND_FLUX,
            Self::CollectedEnergy => PHOTOMETRY_OUTPUT_COLLECTED_ENERGY,
            Self::ExpectedTransmittedPhotons => PHOTOMETRY_OUTPUT_EXPECTED_TRANSMITTED_PHOTONS,
        }
    }
    pub(super) fn unit(self) -> &'static str {
        match self {
            Self::IncidentBandFlux => "watt_per_square_metre",
            Self::CollectedEnergy => "joule",
            Self::ExpectedTransmittedPhotons => "expected_photon_count",
        }
    }
}
#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct Policy {
    pub maximum_rows: u64,
    pub maximum_native_bytes: u64,
}
struct Owned(*mut std::ffi::c_void);
impl Drop for Owned {
    fn drop(&mut self) {
        unsafe {
            irred_photometry_result_destroy(self.0);
        }
    }
}
pub(super) fn scalar(x: &PhotometryScalar, unit: &str) -> Result<(Value, bool), String> {
    let status = numerical_status_name(x.numerical_status).ok_or("UNKNOWN_NUMERICAL_STATUS")?;
    let (state, value, passed) = match x.availability {
        PHOTOMETRY_AVAILABILITY_AVAILABLE if x.numerical_status == 0 && x.value.is_finite() => {
            ("available", Some(x.value), true)
        }
        PHOTOMETRY_AVAILABILITY_FAILED if x.numerical_status != 0 => ("failed", None, false),
        _ => return Err("INVALID_PHOTOMETRY_VIEW".into()),
    };
    Ok((
        json!({"availability":state,"numerical_status":status,"value":value,"unit":unit}),
        passed,
    ))
}
pub(super) fn text(value: &Bytes) -> Result<String, String> {
    String::from_utf8(copied(value.data, value.length, 256)?)
        .map_err(|_| "INVALID_NATIVE_METADATA".into())
}
pub(crate) struct Calculation {
    pub output: Value,
    pub checks: Vec<(&'static str, bool)>,
    pub model: String,
    pub constants: String,
    pub arithmetic: String,
    pub propagation: String,
}
pub(crate) fn evaluate(
    inputs: &[Input],
    requested: &[Output],
    policy: &Policy,
) -> Result<Calculation, String> {
    if inputs.len() > 65536
        || policy.maximum_rows > 65536
        || policy.maximum_native_bytes > (1 << 30)
    {
        return Err("PHOTOMETRY_HARD_LIMIT".into());
    }
    let mut mask = 0;
    for group in requested {
        if mask & group.bit() != 0 {
            return Err("DUPLICATE_REQUESTED_OUTPUT".into());
        }
        mask |= group.bit();
    }
    if mask == 0 {
        return Err("EMPTY_REQUESTED_OUTPUTS".into());
    }
    let wire: Vec<_> = inputs.iter().map(Input::descriptor).collect();
    let batch = PhotometryBatch {
        struct_size: size_of::<PhotometryBatch>() as u32,
        abi_version: ABI_VERSION,
        data: wire.as_ptr(),
        length: wire.len() as u64,
        byte_length: std::mem::size_of_val(wire.as_slice()) as u64,
    };
    let native_policy = PhotometryPolicy {
        struct_size: size_of::<PhotometryPolicy>() as u32,
        abi_version: ABI_VERSION,
        requested_outputs: mask,
        reserved: 0,
        maximum_rows: policy.maximum_rows,
        maximum_native_bytes: policy.maximum_native_bytes,
    };
    let mut raw = ptr::null_mut();
    let transport = unsafe { irred_photometry_evaluate(&batch, &native_policy, &mut raw) };
    if transport != OK {
        return Err(format!("CORE_STATUS_{transport}"));
    }
    if raw.is_null() {
        return Err("NULL_RESULT".into());
    }
    let owned = Owned(raw);
    let mut view: PhotometryView = unsafe { std::mem::zeroed() };
    let transport = unsafe { irred_photometry_result_view(owned.0, &mut view) };
    if transport != OK {
        return Err(format!("CORE_STATUS_{transport}"));
    }
    if view.struct_size != size_of::<PhotometryView>() as u32
        || view.abi_version != ABI_VERSION
        || view.reserved != 0
        || view.row_count > inputs.len() as u64
        || (view.row_count != 0 && view.rows.is_null())
    {
        return Err("INVALID_PHOTOMETRY_VIEW".into());
    }
    let status = numerical_status_name(view.numerical_status).ok_or("UNKNOWN_NUMERICAL_STATUS")?;
    if (view.numerical_status == 0 && view.row_count != inputs.len() as u64)
        || (view.numerical_status != 0 && view.row_count != 0)
    {
        return Err("INVALID_PHOTOMETRY_VIEW".into());
    }
    let rows = if view.row_count == 0 {
        &[][..]
    } else {
        unsafe { std::slice::from_raw_parts(view.rows, view.row_count as usize) }
    };
    let mut checks: Vec<_> = requested
        .iter()
        .map(|g| (g.id(), view.numerical_status == 0))
        .collect();
    let mut values = Vec::with_capacity(rows.len());
    for (index, row) in rows.iter().enumerate() {
        if row.struct_size != size_of::<PhotometryRow>() as u32
            || row.abi_version != ABI_VERSION
            || row.reserved != 0
        {
            return Err("INVALID_PHOTOMETRY_VIEW".into());
        }
        let admission =
            numerical_status_name(row.admission_status).ok_or("UNKNOWN_NUMERICAL_STATUS")?;
        let mut outputs = serde_json::Map::new();
        for (j, group) in requested.iter().enumerate() {
            let item = match group {
                Output::IncidentBandFlux => &row.incident_band_flux,
                Output::CollectedEnergy => &row.collected_energy,
                Output::ExpectedTransmittedPhotons => &row.expected_transmitted_photons,
            };
            let (value, passed) = scalar(item, group.unit())?;
            checks[j].1 &= passed && row.admission_status == 0;
            outputs.insert(group.id().into(), value);
        }
        let x = &row.source;
        let source = Input {
            luminosity_watt_per_metre: x.luminosity_watt_per_metre,
            rest_lower_metre: x.rest_lower_metre,
            rest_upper_metre: x.rest_upper_metre,
            observed_lower_metre: x.observed_lower_metre,
            observed_upper_metre: x.observed_upper_metre,
            luminosity_distance_metre: x.luminosity_distance_metre,
            redshift: x.redshift,
            collecting_area_square_metre: x.collecting_area_square_metre,
            optical_transmission: x.optical_transmission,
            observer_exposure_second: x.observer_exposure_second,
        };
        values.push(
            json!({"index":index,"source":source,"admission_status":admission,"outputs":outputs}),
        );
    }
    let passed = checks.iter().all(|(_, passed)| *passed);
    Ok(Calculation {
        output: json!({"kind":if passed{"finite"}else{"failure"},"batch_numerical_status":status,"requested_rows":inputs.len(),"evaluations":values}),
        checks,
        model: text(&view.model_id)?,
        constants: text(&view.constants_id)?,
        arithmetic: text(&view.arithmetic_id)?,
        propagation: text(&view.propagation_id)?,
    })
}
