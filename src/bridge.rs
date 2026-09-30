#[allow(dead_code)]
#[path = "abi_generated.rs"]
mod generated;
pub use generated::ABI_VERSION;
use generated::*;
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
