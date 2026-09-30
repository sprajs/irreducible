//! Validate raw JSON before intermediate Values can erase duplicate members.
//! Input size is bounded by the caller; serde_json's recursion limit is retained.
use serde::de::{self, DeserializeSeed, MapAccess, SeqAccess, Visitor};
use std::{collections::HashSet, fmt};

struct Unique;
impl<'de> DeserializeSeed<'de> for Unique {
    type Value = ();
    fn deserialize<D: de::Deserializer<'de>>(self, d: D) -> Result<(), D::Error> {
        d.deserialize_any(self)
    }
}
impl<'de> Visitor<'de> for Unique {
    type Value = ();
    fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.write_str("JSON with unique object members")
    }
    fn visit_bool<E: de::Error>(self, _: bool) -> Result<(), E> {
        Ok(())
    }
    fn visit_i64<E: de::Error>(self, _: i64) -> Result<(), E> {
        Ok(())
    }
    fn visit_u64<E: de::Error>(self, _: u64) -> Result<(), E> {
        Ok(())
    }
    fn visit_f64<E: de::Error>(self, _: f64) -> Result<(), E> {
        Ok(())
    }
    fn visit_str<E: de::Error>(self, _: &str) -> Result<(), E> {
        Ok(())
    }
    fn visit_unit<E: de::Error>(self) -> Result<(), E> {
        Ok(())
    }
    fn visit_seq<A: SeqAccess<'de>>(self, mut a: A) -> Result<(), A::Error> {
        while a.next_element_seed(Unique)?.is_some() {}
        Ok(())
    }
    fn visit_map<A: MapAccess<'de>>(self, mut a: A) -> Result<(), A::Error> {
        let mut keys = HashSet::new();
        while let Some(key) = a.next_key::<String>()? {
            if !keys.insert(key) {
                return Err(de::Error::custom("duplicate JSON member"));
            }
            a.next_value_seed(Unique)?;
        }
        Ok(())
    }
}
pub(crate) fn validate(input: &[u8]) -> Result<(), serde_json::Error> {
    let mut d = serde_json::Deserializer::from_slice(input);
    Unique.deserialize(&mut d)?;
    d.end()
}
#[cfg(test)]
mod tests {
    use super::validate;
    #[test]
    fn recursive_and_decoded_names() {
        for raw in [
            r#"{"q":3,"q":0}"#,
            r#"{"x":[{"q":1,"\u0071":2}]}"#,
            r#"{"unused":{"a":1,"a":1}}"#,
        ] {
            assert!(validate(raw.as_bytes()).is_err());
        }
        assert!(validate(br#"{"x":[{"q":1},{"q":2}],"s":"q"}"#).is_ok());
        assert!(validate(b"{} {}").is_err());
    }
}
