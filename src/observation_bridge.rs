//! Synchronous coarse prepared-object bridge. Native scientific semantics only.
pub(crate) use super::generated::ObservationMetadataFields as ObservationMetadata;
use super::generated::*;
use std::{ffi::c_void, ptr};
pub(crate) struct ObservationInput {
    pub metadata: ObservationMetadata,
    pub table_sha256: String,
    pub uncertainty_sha256: String,
    pub calibration_provenance: String,
    pub dependence_provenance: String,
    pub quality_dictionary: String,
    pub ordering_provenance: String,
    pub measurement_ids: Vec<String>,
    pub event_ids: Vec<String>,
    pub uncertainty_axis_ids: Vec<String>,
    pub values: Vec<f64>,
    pub zhd: Vec<f64>,
    pub zcmb: Vec<f64>,
    pub zhel: Vec<f64>,
    pub uncertainty_matrix: Vec<f64>,
    pub missing: Vec<u8>,
    pub zhd_missing: Vec<u8>,
    pub zcmb_missing: Vec<u8>,
    pub zhel_missing: Vec<u8>,
    pub source_selection: Vec<u8>,
    pub quality: Vec<u64>,
}
pub(crate) struct Prepared(pub(crate) *mut c_void);
impl Drop for Prepared {
    fn drop(&mut self) {
        unsafe {
            cosmo_observation_destroy(self.0);
        }
    }
}
impl Prepared {
    pub(crate) fn retained_bytes(&self) -> Result<usize, String> {
        let mut bytes = 0;
        let transport = unsafe { cosmo_observation_retained_bytes(self.0, &mut bytes) };
        if transport != OK {
            return Err(format!("CORE_STATUS_{transport}"));
        }
        usize::try_from(bytes).map_err(|_| "RETAINED_BYTES_OVERFLOW".into())
    }
    pub(crate) fn source_view(&self) -> Result<SourceView<'_>, String> {
        let mut raw = std::mem::MaybeUninit::<ObservationDescriptor>::uninit();
        if unsafe { cosmo_observation_source_view(self.0, raw.as_mut_ptr()) } != OK {
            return Err("INVALID_SOURCE_VIEW".into());
        }
        let descriptor = unsafe { raw.assume_init() };
        if descriptor.struct_size != std::mem::size_of::<ObservationDescriptor>() as u32
            || descriptor.abi_version != ABI_VERSION
        {
            return Err("INVALID_SOURCE_VIEW".into());
        }
        Ok(SourceView {
            descriptor,
            _owner: std::marker::PhantomData,
        })
    }
}
pub(crate) struct SourceView<'a> {
    pub(crate) descriptor: ObservationDescriptor,
    _owner: std::marker::PhantomData<&'a Prepared>,
}
impl SourceView<'_> {
    // Borrow for hashing/export only. No full-matrix readback Vec or repeated
    // preparation is needed by a retained consumer evaluation.
    pub(crate) fn values(&self) -> Result<&[f64], String> {
        self.doubles(&self.descriptor.values)
    }
    pub(crate) fn covariance(&self) -> Result<&[f64], String> {
        self.doubles(&self.descriptor.uncertainty_matrix)
    }
    fn doubles<'a>(&'a self, buffer: &F64Buffer) -> Result<&'a [f64], String> {
        let n = usize::try_from(buffer.length).map_err(|_| "INVALID_SOURCE_LENGTH")?;
        if buffer.struct_size != std::mem::size_of::<F64Buffer>() as u32
            || buffer.abi_version != ABI_VERSION
            || buffer.element_type != 2
            || buffer.reserved != 0
            || n > isize::MAX as usize / std::mem::size_of::<f64>()
            || buffer.byte_length != n as u64 * 8
            || (n > 0
                && (buffer.data.is_null()
                    || (buffer.data as usize) % std::mem::align_of::<f64>() != 0))
        {
            return Err("INVALID_SOURCE_VIEW".into());
        }
        Ok(if n == 0 {
            &[]
        } else {
            unsafe { std::slice::from_raw_parts(buffer.data, n) }
        })
    }
}
pub(super) fn bytes(s: &str) -> Bytes {
    Bytes {
        data: s.as_ptr(),
        length: s.len() as u64,
    }
}
pub(super) fn strings(v: &[Bytes]) -> Strings {
    Strings {
        data: v.as_ptr(),
        length: v.len() as u64,
        byte_length: std::mem::size_of_val(v) as u64,
    }
}
pub(super) fn doubles(v: &[f64]) -> F64Buffer {
    F64Buffer {
        struct_size: std::mem::size_of::<F64Buffer>() as u32,
        abi_version: ABI_VERSION,
        element_type: 2,
        reserved: 0,
        data: v.as_ptr(),
        length: v.len() as u64,
        byte_length: std::mem::size_of_val(v) as u64,
    }
}
fn mask(v: &[u8]) -> U8Buffer {
    U8Buffer {
        data: v.as_ptr(),
        length: v.len() as u64,
        byte_length: v.len() as u64,
    }
}
fn semantic(status: u32) -> Result<(), String> {
    if status == OBSERVATION_STATUS_OK {
        Ok(())
    } else {
        Err(format!(
            "OBSERVATION_{}",
            observation_status_name(status).ok_or("INVALID_SCIENTIFIC_STATUS")?
        ))
    }
}
pub(super) fn copied<T: Copy>(p: *const T, n: u64, cap: usize) -> Result<Vec<T>, String> {
    let n = usize::try_from(n).map_err(|_| "INVALID_RESULT_LENGTH")?;
    if n > cap
        || n > isize::MAX as usize / std::mem::size_of::<T>()
        || (n > 0 && (p.is_null() || (p as usize) % std::mem::align_of::<T>() != 0))
    {
        return Err("INVALID_RESULT_VIEW".into());
    }
    Ok(if n == 0 {
        vec![]
    } else {
        unsafe { std::slice::from_raw_parts(p, n) }.to_vec()
    })
}
pub(crate) fn prepare_shared_bounded(
    s: &ObservationInput,
    policy: [u64; 3],
    maximum_native_bytes: usize,
) -> Result<Prepared, String> {
    let tag = |g: &str, v: &str| {
        observation_tag_id(g, v).ok_or_else(|| format!("UNKNOWN_OBSERVATION_TAG_{g}:{v}"))
    };
    let ids: Vec<_> = s.measurement_ids.iter().map(|x| bytes(x)).collect();
    let events: Vec<_> = s.event_ids.iter().map(|x| bytes(x)).collect();
    let axes: Vec<_> = s.uncertainty_axis_ids.iter().map(|x| bytes(x)).collect();
    let m = &s.metadata;
    let d = ObservationDescriptor {
        struct_size: std::mem::size_of::<ObservationDescriptor>() as u32,
        abi_version: ABI_VERSION,
        reserved: 0,
        profile: tag("profile", &m.profile)?,
        role: tag("role", &m.role)?,
        unit: tag("unit", &m.unit)?,
        calibration: tag("calibration", &m.calibration)?,
        uncertainty: tag("uncertainty", &m.uncertainty)?,
        uncertainty_unit: tag("uncertainty_unit", &m.uncertainty_unit)?,
        component: tag("component", &m.component)?,
        table_sha256: bytes(&s.table_sha256),
        uncertainty_sha256: bytes(&s.uncertainty_sha256),
        calibration_provenance: bytes(&s.calibration_provenance),
        dependence_provenance: bytes(&s.dependence_provenance),
        quality_dictionary: bytes(&s.quality_dictionary),
        ordering_provenance: bytes(&s.ordering_provenance),
        measurement_ids: strings(&ids),
        event_ids: strings(&events),
        uncertainty_axis_ids: strings(&axes),
        values: doubles(&s.values),
        zhd: doubles(&s.zhd),
        zcmb: doubles(&s.zcmb),
        zhel: doubles(&s.zhel),
        uncertainty_matrix: doubles(&s.uncertainty_matrix),
        missing: mask(&s.missing),
        zhd_missing: mask(&s.zhd_missing),
        zcmb_missing: mask(&s.zcmb_missing),
        zhel_missing: mask(&s.zhel_missing),
        source_selection: mask(&s.source_selection),
        quality: U64Buffer {
            data: s.quality.as_ptr(),
            length: s.quality.len() as u64,
            byte_length: std::mem::size_of_val(&*s.quality) as u64,
        },
    };
    let p = ObservationPolicy {
        maximum_rows: policy[0],
        maximum_matrix_elements: policy[1],
        maximum_string_bytes: policy[2],
    };
    let mut preparation_bytes = 0;
    let preflight = unsafe { cosmo_observation_preparation_bytes(&d, &mut preparation_bytes) };
    if preflight != OK {
        return Err(format!("CORE_STATUS_{preflight}"));
    }
    if preparation_bytes > maximum_native_bytes as u64 {
        return Err("PREPARATION_BYTE_LIMIT".into());
    }
    let mut raw = ptr::null_mut();
    let mut status = u32::MAX;
    let transport = unsafe { cosmo_prepare_observations(&d, &p, &mut raw, &mut status) };
    // Defensive adoption even if a faulty native implementation returns an owner on failure.
    let owned = Prepared(raw);
    if transport != OK {
        return Err(format!("CORE_STATUS_{transport}"));
    }
    semantic(status)?;
    if owned.0.is_null() {
        return Err("NULL_PREPARED".into());
    }
    Ok(owned)
}
