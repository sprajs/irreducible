//! Bounded pooled transport; compiled native code owns scientific admission.
use super::{
    generated::*,
    photometry::{self, Calculation, Output},
};
use serde::{
    de::{self, SeqAccess, Visitor},
    Deserialize, Deserializer, Serialize,
};
use serde_json::json;
use std::{collections::HashSet, fmt, mem::size_of, ops::Deref, ptr};
#[derive(Serialize)]
#[serde(transparent)]
pub(crate) struct Bounded<T, const N: usize>(pub Vec<T>);
impl<T, const N: usize> Deref for Bounded<T, N> {
    type Target = [T];
    fn deref(&self) -> &[T] {
        &self.0
    }
}
impl<'de, T: Deserialize<'de>, const N: usize> Deserialize<'de> for Bounded<T, N> {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        struct V<T, const N: usize>(std::marker::PhantomData<T>);
        impl<'de, T: Deserialize<'de>, const N: usize> Visitor<'de> for V<T, N> {
            type Value = Bounded<T, N>;
            fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
                write!(f, "array with at most {N} entries")
            }
            fn visit_seq<A: SeqAccess<'de>>(self, mut a: A) -> Result<Self::Value, A::Error> {
                if a.size_hint().is_some_and(|n| n > N) {
                    return Err(de::Error::custom("SAMPLED_PHOTOMETRY_HARD_LIMIT"));
                }
                let mut v = Vec::new();
                while v.len() < N {
                    match a.next_element()? {
                        Some(x) => v.push(x),
                        None => return Ok(Bounded(v)),
                    }
                }
                if a.next_element::<de::IgnoredAny>()?.is_some() {
                    return Err(de::Error::custom("SAMPLED_PHOTOMETRY_HARD_LIMIT"));
                }
                Ok(Bounded(v))
            }
        }
        d.deserialize_seq(V::<T, N>(std::marker::PhantomData))
    }
}
#[derive(Serialize, PartialEq, Eq, Hash)]
#[serde(transparent)]
pub(crate) struct Text<const N: usize>(String);
impl<const N: usize> Deref for Text<N> {
    type Target = str;
    fn deref(&self) -> &str {
        &self.0
    }
}
impl<'de, const N: usize> Deserialize<'de> for Text<N> {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        struct V<const N: usize>;
        impl<'de, const N: usize> Visitor<'de> for V<N> {
            type Value = Text<N>;
            fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
                write!(f, "nonempty string with at most {N} UTF-8 bytes")
            }
            fn visit_str<E: de::Error>(self, s: &str) -> Result<Self::Value, E> {
                if s.is_empty() || s.len() > N {
                    return Err(de::Error::custom("SAMPLED_METADATA_LIMIT"));
                }
                Ok(Text(s.to_owned()))
            }
        }
        d.deserialize_str(V::<N>)
    }
}
#[derive(Deserialize, Serialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum Role {
    Measured,
    FittedSummary,
    CalibrationAsset,
    SyntheticControl,
}
#[derive(Deserialize, Serialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum Calibration {
    Fixed,
    DeclaredUncertaintyExcluded,
}
#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct Spectrum {
    pub id: Text<256>,
    pub source_role: Role,
    pub provenance: Text<1024>,
    pub rest_wavelength_metre: Bounded<f64, 65536>,
    pub luminosity_watt_per_metre: Bounded<f64, 65536>,
}
#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct Passband {
    pub id: Text<256>,
    pub source_role: Role,
    pub provenance: Text<1024>,
    pub calibration: Calibration,
    pub calibration_provenance: Text<1024>,
    pub observed_wavelength_metre: Bounded<f64, 65536>,
    pub optical_transmission: Bounded<f64, 65536>,
}
#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct Exposure {
    pub spectrum_index: u64,
    pub passband_index: u64,
    pub luminosity_distance_metre: f64,
    pub redshift: f64,
    pub collecting_area_square_metre: f64,
    pub observer_exposure_second: f64,
}
impl Exposure {
    fn descriptor(&self) -> SampledPhotometryExposure {
        SampledPhotometryExposure {
            struct_size: size_of::<SampledPhotometryExposure>() as u32,
            abi_version: ABI_VERSION,
            spectrum_index: self.spectrum_index,
            passband_index: self.passband_index,
            luminosity_distance_metre: self.luminosity_distance_metre,
            redshift: self.redshift,
            collecting_area_square_metre: self.collecting_area_square_metre,
            observer_exposure_second: self.observer_exposure_second,
        }
    }
}
#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct Policy {
    pub maximum_rows: u64,
    pub maximum_native_bytes: u64,
    pub maximum_total_samples: u64,
    pub maximum_samples: u64,
    pub maximum_segments: u64,
}
struct Owned(*mut std::ffi::c_void);
impl Drop for Owned {
    fn drop(&mut self) {
        unsafe {
            irred_sampled_photometry_result_destroy(self.0);
        }
    }
}
fn valid_text(s: &str, max: usize) -> bool {
    !s.is_empty() && s.len() <= max
}
fn curve(id: &str, x: &[f64], y: &[f64]) -> SampledPhotometryCurve {
    SampledPhotometryCurve {
        struct_size: size_of::<SampledPhotometryCurve>() as u32,
        abi_version: ABI_VERSION,
        id: Bytes {
            data: id.as_ptr(),
            length: id.len() as u64,
        },
        wavelength_metre: x.as_ptr(),
        values: y.as_ptr(),
        length: x.len() as u64,
        byte_length: std::mem::size_of_val(x) as u64,
    }
}
unsafe fn slice<'a, T>(p: *const T, n: u64, max: usize) -> Result<&'a [T], String> {
    if n > max as u64 || (n != 0 && (p.is_null() || (p as usize) % std::mem::align_of::<T>() != 0))
    {
        return Err("INVALID_SAMPLED_PHOTOMETRY_VIEW".into());
    }
    if n == 0 {
        Ok(&[])
    } else {
        Ok(unsafe { std::slice::from_raw_parts(p, n as usize) })
    }
}
fn pool_retained(
    view: &[SampledPhotometryCurve],
    sources: &[SampledPhotometryCurve],
) -> Result<(), String> {
    if view.len() != sources.len() {
        return Err("INVALID_SAMPLED_PHOTOMETRY_VIEW".into());
    }
    for (v, s) in view.iter().zip(sources) {
        if v.struct_size != size_of::<SampledPhotometryCurve>() as u32
            || v.abi_version != ABI_VERSION
            || v.length != s.length
            || v.byte_length != s.byte_length
            || photometry::text(&v.id)? != photometry::text(&s.id)?
        {
            return Err("INVALID_SAMPLED_PHOTOMETRY_VIEW".into());
        }
        // No scientific-array readback merely to cross the boundary. Shape,
        // identity and pointer admission suffice here; permanent native ABI
        // controls prove exact retained source bits and owner lifetime.
        if v.length != 0
            && (v.wavelength_metre.is_null()
                || v.values.is_null()
                || (v.wavelength_metre as usize) % std::mem::align_of::<f64>() != 0
                || (v.values as usize) % std::mem::align_of::<f64>() != 0)
        {
            return Err("INVALID_SAMPLED_PHOTOMETRY_VIEW".into());
        }
    }
    Ok(())
}
pub(crate) struct SampledCalculation {
    pub calculation: Calculation,
    pub method: String,
}
pub(crate) fn evaluate(
    spectra: &[Spectrum],
    passbands: &[Passband],
    exposures: &[Exposure],
    requested: &[Output],
    q: &Policy,
) -> Result<SampledCalculation, String> {
    if spectra.len() > 4096
        || passbands.len() > 4096
        || exposures.len() > 65536
        || q.maximum_rows > 65536
        || q.maximum_native_bytes > (1 << 30)
        || q.maximum_total_samples > 65536
        || q.maximum_samples > 65536
        || q.maximum_segments > 65536
    {
        return Err("SAMPLED_PHOTOMETRY_HARD_LIMIT".into());
    }
    let mut total = 0usize;
    let mut ids = HashSet::new();
    for s in spectra {
        if !valid_text(&s.id, 256)
            || !valid_text(&s.provenance, 1024)
            || !ids.insert(&s.id)
            || s.rest_wavelength_metre.len() != s.luminosity_watt_per_metre.len()
        {
            return Err("INVALID_SPECTRUM_POOL".into());
        }
        total = total
            .checked_add(s.rest_wavelength_metre.len())
            .ok_or("SAMPLED_PHOTOMETRY_HARD_LIMIT")?;
    }
    ids.clear();
    for t in passbands {
        if !valid_text(&t.id, 256)
            || !valid_text(&t.provenance, 1024)
            || !valid_text(&t.calibration_provenance, 1024)
            || !ids.insert(&t.id)
            || t.observed_wavelength_metre.len() != t.optical_transmission.len()
        {
            return Err("INVALID_PASSBAND_POOL".into());
        }
        total = total
            .checked_add(t.observed_wavelength_metre.len())
            .ok_or("SAMPLED_PHOTOMETRY_HARD_LIMIT")?;
    }
    if total > 65536 {
        return Err("SAMPLED_PHOTOMETRY_HARD_LIMIT".into());
    }
    let mut mask = 0;
    for g in requested {
        if mask & g.bit() != 0 {
            return Err("DUPLICATE_REQUESTED_OUTPUT".into());
        }
        mask |= g.bit();
    }
    if mask == 0 {
        return Err("EMPTY_REQUESTED_OUTPUTS".into());
    }
    let spectra_wire: Vec<_> = spectra
        .iter()
        .map(|s| {
            curve(
                &s.id,
                &s.rest_wavelength_metre,
                &s.luminosity_watt_per_metre,
            )
        })
        .collect();
    let passbands_wire: Vec<_> = passbands
        .iter()
        .map(|t| curve(&t.id, &t.observed_wavelength_metre, &t.optical_transmission))
        .collect();
    let exposure_wire: Vec<_> = exposures.iter().map(Exposure::descriptor).collect();
    let batch = SampledPhotometryBatch {
        struct_size: size_of::<SampledPhotometryBatch>() as u32,
        abi_version: ABI_VERSION,
        spectra: spectra_wire.as_ptr(),
        spectrum_count: spectra_wire.len() as u64,
        spectrum_byte_length: std::mem::size_of_val(spectra_wire.as_slice()) as u64,
        passbands: passbands_wire.as_ptr(),
        passband_count: passbands_wire.len() as u64,
        passband_byte_length: std::mem::size_of_val(passbands_wire.as_slice()) as u64,
        exposures: exposure_wire.as_ptr(),
        exposure_count: exposure_wire.len() as u64,
        exposure_byte_length: std::mem::size_of_val(exposure_wire.as_slice()) as u64,
    };
    let policy = SampledPhotometryPolicy {
        struct_size: size_of::<SampledPhotometryPolicy>() as u32,
        abi_version: ABI_VERSION,
        requested_outputs: mask,
        reserved: 0,
        maximum_rows: q.maximum_rows,
        maximum_native_bytes: q.maximum_native_bytes,
        maximum_total_samples: q.maximum_total_samples,
        maximum_samples: q.maximum_samples,
        maximum_segments: q.maximum_segments,
    };
    let mut raw = ptr::null_mut();
    let status = unsafe { irred_sampled_photometry_evaluate(&batch, &policy, &mut raw) };
    if status != OK {
        return Err(format!("CORE_STATUS_{status}"));
    }
    if raw.is_null() {
        return Err("NULL_RESULT".into());
    }
    let owner = Owned(raw);
    let mut v: SampledPhotometryView = unsafe { std::mem::zeroed() };
    let status = unsafe { irred_sampled_photometry_result_view(owner.0, &mut v) };
    if status != OK {
        return Err(format!("CORE_STATUS_{status}"));
    }
    if v.struct_size != size_of::<SampledPhotometryView>() as u32
        || v.abi_version != ABI_VERSION
        || v.reserved != 0
    {
        return Err("INVALID_SAMPLED_PHOTOMETRY_VIEW".into());
    }
    let batch_status =
        numerical_status_name(v.numerical_status).ok_or("UNKNOWN_NUMERICAL_STATUS")?;
    if v.numerical_status == 0 {
        pool_retained(
            unsafe { slice(v.spectra, v.spectrum_count, 4096)? },
            &spectra_wire,
        )?;
        pool_retained(
            unsafe { slice(v.passbands, v.passband_count, 4096)? },
            &passbands_wire,
        )?;
        if v.row_count != exposures.len() as u64 {
            return Err("INVALID_SAMPLED_PHOTOMETRY_VIEW".into());
        }
    } else if v.row_count != 0 || v.spectrum_count != 0 || v.passband_count != 0 {
        return Err("INVALID_SAMPLED_PHOTOMETRY_VIEW".into());
    }
    let rows = unsafe { slice(v.rows, v.row_count, exposures.len())? };
    let mut checks: Vec<_> = requested
        .iter()
        .map(|g| (g.id(), v.numerical_status == 0))
        .collect();
    let mut values = Vec::with_capacity(rows.len());
    for (index, row) in rows.iter().enumerate() {
        if row.struct_size != size_of::<SampledPhotometryRow>() as u32
            || row.abi_version != ABI_VERSION
            || row.reserved != 0
        {
            return Err("INVALID_SAMPLED_PHOTOMETRY_VIEW".into());
        }
        let x = &row.source;
        let source = Exposure {
            spectrum_index: x.spectrum_index,
            passband_index: x.passband_index,
            luminosity_distance_metre: x.luminosity_distance_metre,
            redshift: x.redshift,
            collecting_area_square_metre: x.collecting_area_square_metre,
            observer_exposure_second: x.observer_exposure_second,
        };
        let expected = &exposures[index];
        if x.struct_size != size_of::<SampledPhotometryExposure>() as u32
            || x.abi_version != ABI_VERSION
            || source.spectrum_index != expected.spectrum_index
            || source.passband_index != expected.passband_index
            || ![
                source.luminosity_distance_metre,
                source.redshift,
                source.collecting_area_square_metre,
                source.observer_exposure_second,
            ]
            .iter()
            .zip([
                expected.luminosity_distance_metre,
                expected.redshift,
                expected.collecting_area_square_metre,
                expected.observer_exposure_second,
            ])
            .all(|(a, b)| a.to_bits() == b.to_bits())
        {
            return Err("SOURCE_RETENTION_MISMATCH".into());
        }
        let admission =
            numerical_status_name(row.admission_status).ok_or("UNKNOWN_NUMERICAL_STATUS")?;
        let mut output = serde_json::Map::new();
        for (j, g) in requested.iter().enumerate() {
            let x = match g {
                Output::IncidentBandFlux => &row.incident_band_flux,
                Output::CollectedEnergy => &row.collected_energy,
                Output::ExpectedTransmittedPhotons => &row.expected_transmitted_photons,
            };
            let (value, passed) = photometry::scalar(x, g.unit())?;
            checks[j].1 &= passed && row.admission_status == 0;
            output.insert(g.id().into(), value);
        }
        values.push(json!({"index":index,"source":source,"spectrum_id":spectra.get(x.spectrum_index as usize).map(|s|&s.id),"passband_id":passbands.get(x.passband_index as usize).map(|t|&t.id),"admission_status":admission,"outputs":output}));
    }
    let passed = checks.iter().all(|(_, p)| *p);
    Ok(SampledCalculation {
        method: photometry::text(&v.method_id)?,
        calculation: Calculation {
            output: json!({"kind":if passed{"finite"}else{"failure"},"batch_numerical_status":batch_status,"requested_rows":exposures.len(),"evaluations":values}),
            checks,
            model: photometry::text(&v.model_id)?,
            constants: photometry::text(&v.constants_id)?,
            arithmetic: photometry::text(&v.arithmetic_id)?,
            propagation: photometry::text(&v.propagation_id)?,
        },
    })
}

#[cfg(test)]
mod tests {
    use super::{Bounded, Text};
    #[test]
    fn structural_bounds_before_owned_field_growth() {
        assert!(serde_json::from_slice::<Bounded<u64, 2>>(b"[1,2]").is_ok());
        assert!(serde_json::from_slice::<Bounded<u64, 2>>(b"[1,2,3]").is_err());
        assert!(serde_json::from_slice::<Text<3>>(br#""abc""#).is_ok());
        assert!(serde_json::from_slice::<Text<3>>(br#""abcd""#).is_err());
        assert!(serde_json::from_slice::<Text<3>>(br#""""#).is_err());
        assert!(serde_json::from_str::<Text<3>>("\"éé\"").is_err());
    }
}
