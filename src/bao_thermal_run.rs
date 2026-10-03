//! Strict bounded structural decoder. All scientific equations remain native.
use crate::{bridge::bao_thermal, outcome::Outcome};
use serde::de::{self, SeqAccess, Visitor};
use serde::ser::SerializeSeq;
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use serde_json::{Value, json};
use sha2::{Digest, Sha256};
use std::{
    fmt,
    mem::{MaybeUninit, size_of},
    ops::Deref,
};
pub(crate) const INPUT_MAX: usize = 1 << 20;
pub(crate) const RUST_CAP: usize = 64 << 20;
const Q: usize = 2 << 20;
// Source/default library identity, private node/growth, decoder/JSON/recorder
// allocation and fallback proofs must be earned together before this changes.
// This is an execution-admission gate, not a second implementation route.
const ALLOCATION_PROFILE_ACCEPTED: bool = false;
pub(crate) struct Text256 {
    bytes: [u8; 256],
    len: usize,
}
impl Deref for Text256 {
    type Target = str;
    fn deref(&self) -> &str {
        unsafe { std::str::from_utf8_unchecked(&self.bytes[..self.len]) }
    }
}
impl Serialize for Text256 {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        s.serialize_str(self)
    }
}
impl<'de> Deserialize<'de> for Text256 {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        struct V;
        impl<'de> Visitor<'de> for V {
            type Value = Text256;
            fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
                f.write_str("1..256 UTF-8 bytes")
            }
            fn visit_str<E: de::Error>(self, v: &str) -> Result<Text256, E> {
                if v.is_empty() || v.len() > 256 {
                    return Err(E::custom("THERMAL_BAO_TEXT_LIMIT"));
                }
                let mut out = Text256 {
                    bytes: [0; 256],
                    len: v.len(),
                };
                out.bytes[..v.len()].copy_from_slice(v.as_bytes());
                Ok(out)
            }
        }
        d.deserialize_str(V)
    }
}
pub(crate) struct FixedPool<T, const N: usize> {
    data: [MaybeUninit<T>; N],
    len: usize,
}
impl<T, const N: usize> FixedPool<T, N> {
    fn new() -> Self {
        Self {
            data: [const { MaybeUninit::uninit() }; N],
            len: 0,
        }
    }
}
impl<T, const N: usize> Deref for FixedPool<T, N> {
    type Target = [T];
    fn deref(&self) -> &[T] {
        unsafe { std::slice::from_raw_parts(self.data.as_ptr().cast(), self.len) }
    }
}
impl<T, const N: usize> Drop for FixedPool<T, N> {
    fn drop(&mut self) {
        for x in &mut self.data[..self.len] {
            unsafe {
                x.assume_init_drop();
            }
        }
    }
}
impl<T: Serialize, const N: usize> Serialize for FixedPool<T, N> {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        let mut seq = s.serialize_seq(Some(self.len))?;
        for x in self.iter() {
            seq.serialize_element(x)?;
        }
        seq.end()
    }
}
impl<'de, T: Deserialize<'de>, const N: usize> Deserialize<'de> for FixedPool<T, N> {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        struct V<T, const N: usize>(std::marker::PhantomData<T>);
        impl<'de, T: Deserialize<'de>, const N: usize> Visitor<'de> for V<T, N> {
            type Value = FixedPool<T, N>;
            fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
                f.write_str("bounded inline array")
            }
            fn visit_seq<A: SeqAccess<'de>>(self, mut a: A) -> Result<Self::Value, A::Error> {
                let mut out = FixedPool::new();
                while out.len < N {
                    match a.next_element::<T>()? {
                        Some(x) => {
                            out.data[out.len].write(x);
                            out.len += 1;
                        }
                        None => return Ok(out),
                    }
                }
                if a.next_element::<de::IgnoredAny>()?.is_some() {
                    return Err(de::Error::custom("THERMAL_BAO_ARRAY_LIMIT"));
                }
                Ok(out)
            }
        }
        d.deserialize_seq(V::<T, N>(std::marker::PhantomData))
    }
}
#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct Query {
    pub redshift: f64,
    pub observable: Text256,
}
#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct Observation {
    pub ordered_row_ids: FixedPool<Text256, 64>,
    pub queries: FixedPool<Query, 64>,
    pub observed_ratios: FixedPool<f64, 64>,
    pub covariance_axis_ids: FixedPool<Text256, 64>,
    pub covariance_row_major: FixedPool<f64, 4096>,
    pub ratio_unit: Text256,
    pub covariance_unit: Text256,
    pub table_identity: Text256,
    pub covariance_identity: Text256,
    pub ordering_provenance: Text256,
    pub calibration_provenance: Text256,
    pub dependence_provenance: Text256,
    pub redshift_convention: Text256,
}
#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct Species {
    pub mass_ev: f64,
    pub temperature_today_kelvin: f64,
    pub statistical_weight: f64,
}
#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct PhysicalModel {
    pub h0_km_s_mpc: f64,
    pub physical_baryon_density: f64,
    pub physical_cdm_density: f64,
    pub tcmb_kelvin: f64,
    pub physical_massless_nonphoton_density: f64,
    pub species: FixedPool<Species, 16>,
}
#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct Model {
    pub id: Text256,
    pub physical_model: PhysicalModel,
    pub z_drag: f64,
    pub drag_origin: Text256,
    pub source_origin: Text256,
}
#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct ThermalPolicy {
    pub absolute_tolerance: f64,
    pub relative_tolerance: f64,
    pub maximum_callbacks_per_evaluation: u64,
    pub maximum_total_callbacks: u64,
    pub maximum_depth: u32,
    pub maximum_native_bytes: u64,
    pub momentum_method: Text256,
}
#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct PredictionPolicy {
    pub absolute_tolerance_mpc: f64,
    pub relative_tolerance: f64,
    pub absolute_tolerance_ratio: f64,
    pub relative_tolerance_ratio: f64,
    pub maximum_callbacks_per_point: u64,
    pub maximum_total_callbacks: u64,
    pub maximum_depth: u32,
    pub maximum_native_bytes: u64,
    pub thermal: ThermalPolicy,
}
#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct Policy {
    pub arithmetic: Text256,
    pub maximum_rows: u64,
    pub maximum_matrix_elements: u64,
    pub maximum_models: u64,
    pub maximum_species_per_model: u64,
    pub maximum_string_bytes: u64,
    pub maximum_output_array_elements: u64,
    pub maximum_native_bytes: u64,
    pub maximum_preparation_native_bytes: u64,
    pub maximum_evaluation_native_bytes: u64,
    pub maximum_total_callbacks: u64,
    pub maximum_forward_sensitivity: f64,
    pub maximum_projection_log_density_error: f64,
    pub predictions: PredictionPolicy,
}
#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct Request {
    pub schema_version: u32,
    pub operation: Text256,
    pub source_semantics: Text256,
    pub observation: Observation,
    pub models: FixedPool<Model, 16>,
    pub requested_outputs: FixedPool<Text256, 3>,
    pub resource_policy: Policy,
}
impl Request {
    pub(crate) fn mask(&self) -> Result<u32, String> {
        let mut mask = 0;
        for x in self.requested_outputs.iter() {
            let bit = match &**x {
                "normalized_density" => 1,
                "predictions" => 2,
                "residuals" => 4,
                _ => return Err("THERMAL_BAO_OUTPUT".into()),
            };
            if mask & bit != 0 {
                return Err("THERMAL_BAO_DUPLICATE_OUTPUT".into());
            }
            mask |= bit;
        }
        if mask == 0 {
            return Err("THERMAL_BAO_OUTPUT".into());
        }
        Ok(mask)
    }
    fn validate(&self) -> Result<usize, String> {
        let a = &self.observation;
        let n = a.queries.len();
        if self.schema_version != 2
            || &*self.operation != "bao.thermal_density"
            || &*self.source_semantics != "synthetic_controls"
            || n == 0
            || self.models.is_empty()
            || a.ordered_row_ids.len() != n
            || a.observed_ratios.len() != n
            || a.covariance_axis_ids.len() != n
            || a.covariance_row_major.len() != n * n
            || &*a.ratio_unit != "one"
            || &*a.covariance_unit != "ratio_squared"
            || !matches!(&*self.resource_policy.arithmetic, "wide" | "binary64")
            || &*self.resource_policy.predictions.thermal.momentum_method != "direct_adaptive"
        {
            return Err("THERMAL_BAO_STRUCTURE".into());
        }
        let p = &self.resource_policy;
        let pr = &p.predictions;
        let tp = &pr.thermal;
        if p.maximum_rows > 64
            || p.maximum_matrix_elements > 4096
            || p.maximum_models > 16
            || p.maximum_species_per_model > 16
            || p.maximum_string_bytes > 65536
            || p.maximum_output_array_elements > 2048
            || p.maximum_native_bytes > 64 << 20
            || p.maximum_preparation_native_bytes > p.maximum_native_bytes
            || p.maximum_evaluation_native_bytes > p.maximum_native_bytes
            || p.maximum_total_callbacks > 200000000
            || pr.maximum_callbacks_per_point > 20000000
            || pr.maximum_total_callbacks > 100000000
            || pr.maximum_depth > 64
            || pr.maximum_native_bytes > 64 << 20
            || tp.maximum_callbacks_per_evaluation > 200000
            || tp.maximum_total_callbacks > 100000000
            || tp.maximum_depth > 64
            || tp.maximum_native_bytes > 64 << 20
        {
            return Err("THERMAL_BAO_POLICY".into());
        }
        for x in [
            p.maximum_forward_sensitivity,
            p.maximum_projection_log_density_error,
            pr.absolute_tolerance_mpc,
            pr.relative_tolerance,
            pr.absolute_tolerance_ratio,
            pr.relative_tolerance_ratio,
            tp.absolute_tolerance,
            tp.relative_tolerance,
        ] {
            if !x.is_finite() || x <= 0. {
                return Err("THERMAL_BAO_POLICY".into());
            }
        }
        self.mask()?;
        for i in 0..n {
            if &*a.ordered_row_ids[i] != &*a.covariance_axis_ids[i]
                || a.ordered_row_ids[..i]
                    .iter()
                    .any(|x| &**x == &*a.ordered_row_ids[i])
                || !matches!(
                    &*a.queries[i].observable,
                    "DM_over_rs" | "DH_over_rs" | "DV_over_rs"
                )
            {
                return Err("THERMAL_BAO_ORDER".into());
            }
        }
        for i in 0..self.models.len() {
            if self.models[..i]
                .iter()
                .any(|x| &*x.id == &*self.models[i].id)
            {
                return Err("THERMAL_BAO_MODEL_ID".into());
            }
        }
        let mut chars = 0usize;
        let mut add = |x: &str| {
            chars = chars.checked_add(x.len()).unwrap_or(usize::MAX);
        };
        for x in [
            &*self.operation,
            &*self.source_semantics,
            &*a.ratio_unit,
            &*a.covariance_unit,
            &*a.table_identity,
            &*a.covariance_identity,
            &*a.ordering_provenance,
            &*a.calibration_provenance,
            &*a.dependence_provenance,
            &*a.redshift_convention,
            &*self.resource_policy.arithmetic,
            &*self.resource_policy.predictions.thermal.momentum_method,
        ] {
            add(x);
        }
        for x in a
            .ordered_row_ids
            .iter()
            .chain(a.covariance_axis_ids.iter())
            .chain(self.requested_outputs.iter())
        {
            add(x);
        }
        for q in a.queries.iter() {
            add(&q.observable);
        }
        for m in self.models.iter() {
            add(&m.id);
            add(&m.drag_origin);
            add(&m.source_origin);
        }
        if chars > 65536 {
            return Err("THERMAL_BAO_TEXT_LIMIT".into());
        }
        Ok(chars)
    }
}
pub(crate) struct Envelope {
    pub decoder: usize,
    pub bridge: usize,
    pub completion: usize,
    pub recorder: usize,
    pub peak: usize,
}
fn sum(terms: &[(usize, usize)]) -> Option<usize> {
    terms
        .iter()
        .try_fold(0usize, |a, (n, w)| a.checked_add(n.checked_mul(*w)?))
}
pub(crate) fn envelope(
    bytes: usize,
    n: usize,
    m: usize,
    species: usize,
    chars: usize,
) -> Option<Envelope> {
    let k = sum(&[
        (n, n),
        (n, 12),
        (m, 96),
        (species, 12),
        (m.checked_mul(n)?, 2),
        (1, 512),
    ])?;
    let e = sum(&[(n, 16), (m, 96), (species, 16), (1, 256)])?;
    let maps = sum(&[(n, 4), (m, 16), (species, 8), (1, 128)])?;
    let l = sum(&[(n, 4), (m, 16), (species, 4), (1, 256)])?;
    let t = sum(&[(chars, 2), (e, 64), (1, 32768)])?;
    let a = sum(&[
        (k, 256),
        (maps, 512),
        (e, 1024),
        (t, 4),
        (e.checked_add(l)?, 16),
    ])?;
    let j = sum(&[(k, 64), (t, 8), (1, 65536)])?;
    let i = size_of::<Request>();
    let w = bao_thermal::wire_size();
    let decoder = sum(&[(i, 2), (bytes, 6), (1, Q)])?;
    let bridge = sum(&[(i, 1), (w, 1), (a, 1), (j, 1), (1, Q)])?;
    let completion = sum(&[(i, 1), (a, 2), (j, 1), (1, Q)])?;
    let recorder = sum(&[(a, 5), (j, 3), (1, Q)])?;
    Some(Envelope {
        decoder,
        bridge,
        completion,
        recorder,
        peak: decoder.max(bridge).max(completion).max(recorder),
    })
}
fn refusal(
    raw: &[u8],
    cause: &'static str,
    env: Option<&Envelope>,
    requested: Option<(&'static str, u32)>,
) -> Outcome {
    let digest = format!("{:x}", Sha256::digest(raw));
    let resources = json!({"rust_operation_payload_cap":RUST_CAP,"short_refusal_conditional_bound":Q,
        "allocation_profile":"unqualified","bounds_measured":false,"outer_cli_logical_limit":16<<20,
        "common_framework_allocator_stack_rss_excluded":true,
        "conditional_peak":env.map(|x|x.peak)});
    Outcome::scientific(
        json!({"schema_version":2,"operation":"bao.thermal_density","raw_input_sha256":digest}),
        json!({"kind":"failure","error_id":"NUMERICAL_QUALIFICATION_REQUIRED","cause":cause,
            "native_payload_absent":true,"source_prepare_call_attempted":false,"thermal_batch_call_attempted":false}),
        "thermal_bao_requested_outputs",
        json!({"requested":"retained thermal BAO normalized ratio density","requested_output_mask":requested.map(|x|x.1),"actual":null}),
        json!({"requested":requested.map(|x|x.0),"requested_unresolved_before_decode":requested.is_none(),"actual":null}),
        resources,
        "source-only resource/profile refusal; no native equation executed; fallback bound remains conditional",
    )
}
pub(crate) fn execute(raw: &[u8]) -> Result<Outcome, String> {
    if raw.len() > INPUT_MAX {
        return Err("THERMAL_BAO_INPUT_LIMIT".into());
    }
    if size_of::<Request>() > 256 << 10
        || bao_thermal::wire_size() > 128 << 10
        || size_of::<Value>() > 64
        || size_of::<String>() > 32
        || size_of::<serde_json::Map<String, Value>>() > 64
    {
        return Ok(refusal(raw, "unsupported_public_layout", None, None));
    }
    let decoder = sum(&[(size_of::<Request>(), 2), (raw.len(), 6), (1, Q)])
        .ok_or("THERMAL_BAO_RESOURCE_OVERFLOW")?;
    if decoder > RUST_CAP {
        return Ok(refusal(raw, "decoder_resource_limit", None, None));
    }
    // Typed visitors own one final inline object; no Value request tree or
    // scientific array flattening. Derive's object field masks reject duplicates.
    let r: Box<Request> =
        Box::new(serde_json::from_slice(raw).map_err(|_| "THERMAL_BAO_REQUEST".to_string())?);
    let chars = r.validate()?;
    let requested = Some((
        if &*r.resource_policy.arithmetic == "wide" {
            "wide"
        } else {
            "binary64"
        },
        r.mask()?,
    ));
    let species = r
        .models
        .iter()
        .map(|x| x.physical_model.species.len())
        .sum();
    let env = envelope(
        raw.len(),
        r.observation.queries.len(),
        r.models.len(),
        species,
        chars,
    );
    if env.as_ref().is_none_or(|x| x.peak > RUST_CAP) {
        drop(r);
        return Ok(refusal(raw, "rust_resource_limit", env.as_ref(), requested));
    }
    if !ALLOCATION_PROFILE_ACCEPTED {
        drop(r);
        return Ok(refusal(
            raw,
            "unsupported_allocation_profile",
            env.as_ref(),
            requested,
        ));
    }
    bao_thermal::evaluate(r, env.unwrap())
}
