//! Synchronous coarse prepared-object bridge. Native scientific semantics only.
pub(crate) use super::generated::ObservationMetadataFields as ObservationMetadata;
use super::{generated::*, Owned};
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
struct Prepared(*mut c_void);
impl Drop for Prepared {
    fn drop(&mut self) {
        unsafe {
            cosmo_observation_destroy(self.0);
        }
    }
}
fn bytes(s: &str) -> Bytes {
    Bytes {
        data: s.as_ptr(),
        length: s.len() as u64,
    }
}
fn strings(v: &[Bytes]) -> Strings {
    Strings {
        data: v.as_ptr(),
        length: v.len() as u64,
        byte_length: std::mem::size_of_val(v) as u64,
    }
}
fn doubles(v: &[f64]) -> F64Buffer {
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
fn copied<T: Copy>(p: *const T, n: u64, cap: usize) -> Result<Vec<T>, String> {
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
pub(crate) struct ObservationSelection {
    pub mask: Vec<u8>,
    pub source_indices: Vec<u64>,
    pub source_values: Vec<f64>,
    pub retained_matrix: Vec<f64>,
}
pub(crate) fn prepare_observations(
    s: &ObservationInput,
    selection: &str,
    policy: [u64; 3],
) -> Result<ObservationSelection, String> {
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
    let mut view = std::mem::MaybeUninit::<ObservationDescriptor>::uninit();
    if unsafe { cosmo_observation_source_view(owned.0, view.as_mut_ptr()) } != OK {
        return Err("INVALID_SOURCE_VIEW".into());
    }
    let view = unsafe { view.assume_init() };
    if view.struct_size != std::mem::size_of::<ObservationDescriptor>() as u32
        || view.abi_version != ABI_VERSION
    {
        return Err("INVALID_SOURCE_VIEW".into());
    }
    if view.values.length != s.values.len() as u64
        || view.uncertainty_matrix.length != s.uncertainty_matrix.len() as u64
    {
        return Err("INVALID_SOURCE_VIEW_LENGTH".into());
    }
    let source_values = copied(view.values.data, view.values.length, s.values.len())?;
    let retained_matrix = copied(
        view.uncertainty_matrix.data,
        view.uncertainty_matrix.length,
        s.uncertainty_matrix.len(),
    )?;
    let mut raw = ptr::null_mut();
    let transport = unsafe {
        cosmo_observation_select(owned.0, tag("selection", selection)?, &mut raw, &mut status)
    };
    let selection = Owned(raw);
    if transport != OK {
        return Err(format!("CORE_STATUS_{transport}"));
    }
    semantic(status)?;
    if selection.0.is_null() {
        return Err("NULL_SELECTION".into());
    }
    let mut mask_ptr = ptr::null();
    let mut indices_ptr = ptr::null();
    let mut n = 0;
    let mut k = 0;
    if unsafe {
        cosmo_result_selection_view(selection.0, &mut mask_ptr, &mut n, &mut indices_ptr, &mut k)
    } != OK
        || n != s.values.len() as u64
    {
        return Err("INVALID_SELECTION_VIEW".into());
    }
    let mask = copied(mask_ptr, n, s.values.len())?;
    let source_indices = copied(indices_ptr, k, s.values.len())?;
    if mask.iter().any(|x| *x > 1)
        || source_indices.windows(2).any(|w| w[0] >= w[1])
        || source_indices.iter().any(|i| *i >= n)
        || source_indices.len() != mask.iter().filter(|x| **x == 1).count()
        || source_indices.iter().any(|i| mask[*i as usize] != 1)
    {
        return Err("INVALID_SELECTION_VIEW".into());
    }
    Ok(ObservationSelection {
        mask,
        source_indices,
        source_values,
        retained_matrix,
    })
}
