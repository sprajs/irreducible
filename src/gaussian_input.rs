//! Shared bounded Gaussian wire inputs; no statistical equations.
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
pub(crate) type IDs = Bounded<Text, 65536>;
pub(crate) type Numbers = Bounded<f64, 1000000>;
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
        impl<'de> DeserializeSeed<'de> for Row<'_> {
            type Value = usize;
            fn deserialize<D: Deserializer<'de>>(self, d: D) -> Result<usize, D::Error> {
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
