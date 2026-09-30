#[allow(dead_code)]
#[path = "abi_generated.rs"]
mod generated;
use generated::*;
pub use generated::{ABI_VERSION, MAX_BATCH_ELEMENTS};
use std::ptr;
struct Owned(*mut std::ffi::c_void);
impl Drop for Owned {
    fn drop(&mut self) {
        unsafe {
            cosmo_result_destroy(self.0);
        }
    }
}
pub(crate) fn add(a: &[i64], b: &[i64], fault: u32) -> Result<Vec<i64>, String> {
    let desc = |v: &[i64]| Buffer {
        struct_size: std::mem::size_of::<Buffer>() as u32,
        abi_version: ABI_VERSION,
        element_type: 1,
        reserved: 0,
        data: v.as_ptr(),
        length: v.len() as u64,
        byte_length: std::mem::size_of_val(v) as u64,
    };
    let mut raw = ptr::null_mut();
    let status = unsafe { cosmo_add(&desc(a), &desc(b), fault, &mut raw) };
    if status != OK {
        return Err(format!("CORE_STATUS_{status}"));
    }
    if raw.is_null() {
        return Err("NULL_RESULT".into());
    }
    let owned = Owned(raw);
    let mut p = ptr::null();
    let mut n = 0;
    let status = unsafe { cosmo_result_view(owned.0, &mut p, &mut n) };
    if status != OK || n != a.len() as u64 || (n > 0 && (p.is_null() || (p as usize) % 8 != 0)) {
        return Err("INVALID_RESULT_VIEW".into());
    }
    if n == 0 {
        Ok(vec![])
    } else {
        Ok(unsafe { std::slice::from_raw_parts(p, n as usize) }.to_vec())
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn batch() {
        assert_eq!(add(&[1, -3], &[2, 5], 0).unwrap(), [3, 2]);
        assert!(add(&[], &[], 0).unwrap().is_empty());
        assert!(add(&[i64::MAX], &[1], 0).is_err());
        assert!(add(&[1], &[], 0).is_err());
        assert!(add(&[], &[], 1).is_err());
        assert!(add(&[], &[], 2).is_err());
    }
}

pub(crate) use generated::MetadataFields as Metadata;
impl Metadata {
    fn descriptor(&self) -> Result<QuantityMetadata, String> {
        let tag = |group: &str, label: &str| {
            tag_id(group, label).ok_or_else(|| format!("UNKNOWN_METADATA_TAG_{group}:{label}"))
        };
        Ok(QuantityMetadata {
            struct_size: std::mem::size_of::<QuantityMetadata>() as u32,
            abi_version: ABI_VERSION,
            unit: tag("unit", &self.unit)?,
            role: tag("role", &self.role)?,
            frame: tag("frame", &self.frame)?,
            convention: tag("convention", &self.convention)?,
            constant_set: tag("constant_set", &self.constant_set)?,
            reserved: 0,
        })
    }
}
pub(crate) struct QuantitySlot {
    pub status: &'static str,
    pub value: Option<f64>,
}
pub(crate) fn convert_quantities(
    values: &[f64],
    source: &Metadata,
    target: &Metadata,
) -> Result<Vec<QuantitySlot>, String> {
    let source = source.descriptor()?;
    let target = target.descriptor()?;
    let buffer = F64Buffer {
        struct_size: std::mem::size_of::<F64Buffer>() as u32,
        abi_version: ABI_VERSION,
        element_type: 2,
        reserved: 0,
        data: values.as_ptr(),
        length: values.len() as u64,
        byte_length: std::mem::size_of_val(values) as u64,
    };
    let mut raw = ptr::null_mut();
    // All descriptors/contiguous input remain borrowed and immobile until synchronous return.
    let status = unsafe { cosmo_convert_quantities(&buffer, &source, &target, &mut raw) };
    if status != OK {
        return Err(format!("CORE_STATUS_{status}"));
    }
    if raw.is_null() {
        return Err("NULL_RESULT".into());
    }
    let owned = Owned(raw);
    let mut data = ptr::null();
    let mut statuses = ptr::null();
    let mut length = 0;
    // The view is core owned, immutable, and only used while the matching owner lives.
    let status = unsafe { cosmo_result_f64_view(owned.0, &mut data, &mut statuses, &mut length) };
    if status != OK
        || length != values.len() as u64
        || (length > 0
            && (data.is_null()
                || statuses.is_null()
                || (data as usize) % std::mem::align_of::<f64>() != 0
                || (statuses as usize) % std::mem::align_of::<u32>() != 0))
    {
        return Err("INVALID_RESULT_VIEW".into());
    }
    if length == 0 {
        return Ok(vec![]);
    }
    // Length equals the original Rust allocation's length; the ABI supplies same-size arrays.
    let converted = unsafe { std::slice::from_raw_parts(data, length as usize) };
    let statuses = unsafe { std::slice::from_raw_parts(statuses, length as usize) };
    converted
        .iter()
        .zip(statuses)
        .map(|(&value, &status)| {
            let name = quantity_status_name(status).ok_or("UNKNOWN_QUANTITY_STATUS")?;
            if status == 0 && !value.is_finite() {
                return Err("INVALID_FINITE_RESULT".into());
            }
            Ok(QuantitySlot {
                status: name,
                value: (status == 0).then_some(value),
            })
        })
        .collect()
}

pub(crate) struct NumericalSlot {
    pub status: &'static str,
    pub value: Option<f64>,
    pub error_estimate: Option<f64>,
    pub evaluations: u64,
}
pub(crate) fn numerics_evaluate(
    method: &str,
    values: &[f64],
) -> Result<Vec<NumericalSlot>, String> {
    let operation = numerical_operation_id(method).ok_or("UNKNOWN_NUMERICAL_METHOD")?;
    let buffer = F64Buffer {
        struct_size: std::mem::size_of::<F64Buffer>() as u32,
        abi_version: ABI_VERSION,
        element_type: 2,
        reserved: 0,
        data: values.as_ptr(),
        length: values.len() as u64,
        byte_length: std::mem::size_of_val(values) as u64,
    };
    let mut raw = ptr::null_mut();
    // One synchronous compiled batch; no Rust numerical inner loop.
    let status = unsafe { cosmo_numerics_evaluate(operation, &buffer, &mut raw) };
    if status != OK {
        return Err(format!("CORE_STATUS_{status}"));
    }
    if raw.is_null() {
        return Err("NULL_RESULT".into());
    }
    let owner = Owned(raw);
    let mut data = ptr::null();
    let mut statuses = ptr::null();
    let mut errors = ptr::null();
    let mut evaluations = ptr::null();
    let mut length = 0;
    let status = unsafe {
        cosmo_result_numerics_view(
            owner.0,
            &mut data,
            &mut statuses,
            &mut errors,
            &mut evaluations,
            &mut length,
        )
    };
    let expected =
        numerical_output_length(method, values.len()).ok_or("UNKNOWN_NUMERICAL_METHOD")?;
    if status != OK
        || length != expected as u64
        || (length > 0
            && (data.is_null()
                || statuses.is_null()
                || errors.is_null()
                || evaluations.is_null()
                || (data as usize) % 8 != 0
                || (errors as usize) % 8 != 0
                || (statuses as usize) % 4 != 0
                || (evaluations as usize) % 8 != 0))
    {
        return Err("INVALID_NUMERICAL_RESULT_VIEW".into());
    }
    if length == 0 {
        return Ok(vec![]);
    }
    // All four immutable views have matched lengths and remain live until owner drops.
    let data = unsafe { std::slice::from_raw_parts(data, length as usize) };
    let statuses = unsafe { std::slice::from_raw_parts(statuses, length as usize) };
    let errors = unsafe { std::slice::from_raw_parts(errors, length as usize) };
    let evaluations = unsafe { std::slice::from_raw_parts(evaluations, length as usize) };
    (0..expected)
        .map(|i| {
            let name = numerical_status_name(statuses[i]).ok_or("UNKNOWN_NUMERICAL_STATUS")?;
            if statuses[i] == 0 && !data[i].is_finite() {
                return Err("INVALID_FINITE_RESULT".into());
            }
            // Diagnostics are not certified error bounds. Missing finite diagnostics stay absent.
            Ok(NumericalSlot {
                status: name,
                value: (statuses[i] == 0).then_some(data[i]),
                error_estimate: errors[i].is_finite().then_some(errors[i]),
                evaluations: evaluations[i],
            })
        })
        .collect()
}

#[path = "observation_bridge.rs"]
mod observations;
pub(crate) use observations::{prepare_observations, ObservationInput, ObservationMetadata};
