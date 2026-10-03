//! Strict inline structural resolution; all posterior equations remain native.
use crate::{
    bridge::gaussian_posterior,
    outcome::{NumericalCheck, Outcome, OutputCheck},
};
use serde::de::{self, DeserializeSeed, SeqAccess, Visitor};
use serde::ser::SerializeSeq;
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use serde_json::json;
use std::{fmt, ops::Deref};
#[derive(Serialize)]
#[serde(transparent)]
pub(crate) struct Bounded<T, const N: usize>(Vec<T>);
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
                write!(f, "at most {N} entries")
            }
            fn visit_seq<A: SeqAccess<'de>>(self, mut a: A) -> Result<Self::Value, A::Error> {
                let mut v = Vec::new();
                while v.len() < N {
                    match a.next_element()? {
                        Some(x) => v.push(x),
                        None => return Ok(Bounded(v)),
                    }
                }
                if a.next_element::<de::IgnoredAny>()?.is_some() {
                    return Err(de::Error::custom("POSTERIOR_ARRAY_LIMIT"));
                }
                Ok(Bounded(v))
            }
        }
        d.deserialize_seq(V::<T, N>(std::marker::PhantomData))
    }
}
#[derive(Serialize)]
#[serde(transparent)]
pub(crate) struct Text(String);
impl Deref for Text {
    type Target = str;
    fn deref(&self) -> &str {
        &self.0
    }
}
impl<'de> Deserialize<'de> for Text {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        struct V;
        impl<'de> Visitor<'de> for V {
            type Value = Text;
            fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
                f.write_str("1..256 UTF-8 bytes")
            }
            fn visit_str<E: de::Error>(self, s: &str) -> Result<Text, E> {
                if s.is_empty() || s.len() > 256 {
                    Err(E::custom("POSTERIOR_TEXT_LIMIT"))
                } else {
                    Ok(Text(s.into()))
                }
            }
        }
        d.deserialize_str(V)
    }
}
type IDs = Bounded<Text, 65536>;
type Numbers = Bounded<f64, 1000000>;
// Decode nested original vectors directly into one contiguous borrowed ABI pool.
// No row-by-row FFI or second flattening copy is needed.
pub(crate) struct VectorPool {
    pub values: Vec<f64>,
    pub lengths: Vec<usize>,
}
impl Serialize for VectorPool {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        let mut seq = s.serialize_seq(Some(self.lengths.len()))?;
        let mut off = 0;
        for &n in &self.lengths {
            seq.serialize_element(&self.values[off..off + n])?;
            off += n;
        }
        seq.end()
    }
}
impl<'de> Deserialize<'de> for VectorPool {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        struct Pool;
        struct Row<'a>(&'a mut Vec<f64>);
        impl<'de> DeserializeSeed<'de> for Row<'_> {
            type Value = usize;
            fn deserialize<D: Deserializer<'de>>(self, d: D) -> Result<usize, D::Error> {
                impl<'de> Visitor<'de> for Row<'_> {
                    type Value = usize;
                    fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
                        f.write_str("bounded numerical vector")
                    }
                    fn visit_seq<A: SeqAccess<'de>>(self, mut a: A) -> Result<usize, A::Error> {
                        let start = self.0.len();
                        while let Some(x) = a.next_element::<f64>()? {
                            if self.0.len() == 1000000 {
                                return Err(de::Error::custom("POSTERIOR_POOL_LIMIT"));
                            }
                            self.0.push(x);
                        }
                        Ok(self.0.len() - start)
                    }
                }
                d.deserialize_seq(self)
            }
        }
        impl<'de> Visitor<'de> for Pool {
            type Value = VectorPool;
            fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
                f.write_str("bounded case vectors")
            }
            fn visit_seq<A: SeqAccess<'de>>(self, mut a: A) -> Result<VectorPool, A::Error> {
                let mut p = VectorPool {
                    values: vec![],
                    lengths: vec![],
                };
                while p.lengths.len() < 65536 {
                    match a.next_element_seed(Row(&mut p.values))? {
                        Some(n) => p.lengths.push(n),
                        None => return Ok(p),
                    }
                }
                if a.next_element::<de::IgnoredAny>()?.is_some() {
                    return Err(de::Error::custom("POSTERIOR_CASE_LIMIT"));
                }
                Ok(p)
            }
        }
        d.deserialize_seq(Pool)
    }
}
#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct Noise {
    pub ordered_row_ids: IDs,
    pub event_ids: IDs,
    pub residual_unit: Text,
    pub covariance_row_major: Numbers,
    pub noise_identity: Text,
    pub calibration_identity: Text,
    pub dependence_identity: Text,
    pub ordering_provenance: Text,
}
#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct Design {
    pub ordered_row_ids: IDs,
    pub ordered_parameter_ids: IDs,
    pub parameter_units: IDs,
    pub column_units: IDs,
    pub values_row_major: Numbers,
    pub design_identity: Text,
}
#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct Prior {
    pub ordered_parameter_ids: IDs,
    pub parameter_units: IDs,
    pub shared_nuisance_ids: IDs,
    pub mean: Numbers,
    pub covariance_row_major: Numbers,
    pub prior_identity: Text,
    pub parameter_measure: Text,
    pub dependence_identity: Text,
    pub noise_independence_declared: bool,
}
#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct Conditioning {
    pub ordered_row_ids: IDs,
    pub event_ids: IDs,
    pub case_ids: IDs,
    pub vectors: VectorPool,
}
#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct Policy {
    pub arithmetic: String,
    pub maximum_elements: u64,
    pub maximum_cases: u64,
    pub maximum_string_bytes: u64,
    pub maximum_native_bytes: u64,
    pub maximum_work_units: u64,
    pub maximum_forward_sensitivity: f64,
}
#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct Request {
    pub schema_version: u32,
    pub operation: String,
    pub source_semantics: String,
    pub noise: Noise,
    pub design: Design,
    pub parameter_prior: Prior,
    pub conditioning: Conditioning,
    pub resource_policy: Policy,
}
pub(crate) fn execute(bytes: &[u8]) -> Result<Outcome, String> {
    if bytes.len() > 16 * 1024 * 1024 {
        return Err("INPUT_LIMIT".into());
    }
    crate::strict_json::validate(bytes).map_err(|_| "INVALID_JSON")?;
    let r: Request = serde_json::from_slice(bytes).map_err(|e| e.to_string())?;
    if r.schema_version != 2
        || r.operation != "statistics.gaussian_posterior"
        || r.source_semantics != "synthetic_controls"
    {
        return Err("UNSUPPORTED_SPECIFICATION".into());
    }
    if r.conditioning
        .vectors
        .lengths
        .iter()
        .any(|&n| n != r.noise.ordered_row_ids.len())
    {
        return Err("INVALID_BATCH_SHAPE".into());
    }
    let c = gaussian_posterior::evaluate(&r)?;
    let mut spec = serde_json::to_value(&r).map_err(|_| "RECORD_ENCODING")?;
    spec["requested_outputs"] = json!([{"id":"posterior_covariance","required":true,"numerical_gate":"required","inference_gate":"not_applicable"},{"id":"posterior_means","required":true,"numerical_gate":"required","inference_gate":"not_applicable"}]);
    spec["batch_meaning"] = json!(
        "family of conditionals under one fixed noise/design/independent proper prior; not a joint parameter posterior or IID generated draws"
    );
    let mut outcome = Outcome::scientific(
        spec,
        c.output,
        "posterior_means",
        c.method,
        c.arithmetic,
        c.resources,
        "named native proper-Gaussian posterior controls; supplied synthetic prior/data applicability unqualified",
    );
    outcome.outputs=[("posterior_covariance",c.covariance_passed),("posterior_means",c.means_passed)].into_iter().map(|(id,pass)|OutputCheck{check_kind:"numerical_contract",id,required:true,numerical:if pass{NumericalCheck::ChecksPassed}else{NumericalCheck::Failed},inference:"not_applicable",interpretation:"unqualified",evidence:vec![],validation_coverage:"original two-dimensional rational/cofactor and native/ABI/CLI parity; inherited scientific ancestry"}).collect();
    Ok(outcome)
}
