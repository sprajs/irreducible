//! Strict fixed synthetic predictive request; equations remain native.
use crate::gaussian_input::{Conditioning, Design, IDs, Noise, Policy, Prior, Text, VectorPool};
use serde::de::{MapAccess, Visitor, value::MapAccessDeserializer};
use serde::{Deserialize, Deserializer, Serialize};
use std::{fmt, marker::PhantomData};
#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct RequestedOutputs {
    pub means: bool,
    pub joint_log_densities: bool,
}
#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct Prediction {
    pub conditioning_identity: Text,
    pub dependence_identity: Text,
    pub future_covariance_unit: Text,
    pub future_measure: Text,
    pub future_noise_independence_declared: bool,
    pub noise_conditional_on_parameters_declared: bool,
}
#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct FutureVectors {
    pub ordered_row_ids: IDs,
    pub event_ids: IDs,
    pub vectors: VectorPool,
}
#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct Request {
    pub schema_version: u32,
    pub operation: String,
    pub source_semantics: String,
    #[serde(deserialize_with = "object_only")]
    pub noise: Noise,
    #[serde(deserialize_with = "object_only")]
    pub design: Design,
    #[serde(deserialize_with = "object_only")]
    pub parameter_prior: Prior,
    #[serde(deserialize_with = "object_only")]
    pub conditioning: Conditioning,
    #[serde(deserialize_with = "object_only")]
    pub future_noise: Noise,
    #[serde(deserialize_with = "object_only")]
    pub future_response: Design,
    #[serde(deserialize_with = "object_only")]
    pub prediction: Prediction,
    #[serde(deserialize_with = "object_only")]
    pub outputs: RequestedOutputs,
    #[serde(
        default,
        deserialize_with = "present_future_vectors",
        skip_serializing_if = "Option::is_none"
    )]
    pub future_vectors: Option<FutureVectors>,
    #[serde(deserialize_with = "object_only")]
    pub resource_policy: Policy,
}
// This operation requires object syntax even though serde derive may also
// accept positional sequences for structs. Reused posterior types are unchanged.
fn object_only<'de, T: Deserialize<'de>, D: Deserializer<'de>>(d: D) -> Result<T, D::Error> {
    struct Object<T>(PhantomData<T>);
    impl<'de, T: Deserialize<'de>> Visitor<'de> for Object<T> {
        type Value = T;
        fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
            f.write_str("a JSON object")
        }
        fn visit_map<A: MapAccess<'de>>(self, a: A) -> Result<T, A::Error> {
            T::deserialize(MapAccessDeserializer::new(a))
        }
    }
    d.deserialize_map(Object::<T>(PhantomData))
}
// Missing member defaults to None. Every present value must be an object:
// explicit null or positional sequence fails the map-only decoder.
fn present_future_vectors<'de, D: Deserializer<'de>>(
    d: D,
) -> Result<Option<FutureVectors>, D::Error> {
    object_only::<FutureVectors, D>(d).map(Some)
}

pub(crate) fn execute(bytes: &[u8]) -> Result<crate::outcome::Outcome, String> {
    use crate::outcome::{NumericalCheck, Outcome, OutputCheck};
    use serde_json::json;
    if bytes.len() > 16 * 1024 * 1024 {
        return Err("INPUT_LIMIT".into());
    }
    crate::strict_json::validate(bytes).map_err(|_| "INVALID_JSON")?;
    let mut d = serde_json::Deserializer::from_slice(bytes);
    let r: Request = object_only(&mut d).map_err(|e| e.to_string())?;
    d.end().map_err(|e| e.to_string())?;
    if r.schema_version != 2
        || r.operation != "statistics.gaussian_predictive"
        || r.source_semantics != "synthetic_controls"
    {
        return Err("UNSUPPORTED_SPECIFICATION".into());
    }
    if !r.outputs.means && !r.outputs.joint_log_densities {
        return Err("EMPTY_OUTPUT_MASK".into());
    }
    if r.outputs.joint_log_densities != r.future_vectors.is_some() {
        return Err("FUTURE_VECTOR_MASK_MISMATCH".into());
    }
    let n = r.noise.ordered_row_ids.len();
    let k = r.future_noise.ordered_row_ids.len();
    let b = r.conditioning.vectors.lengths.len();
    if r.conditioning.vectors.lengths.iter().any(|&x| x != n)
        || r.conditioning.case_ids.len() != b
        || r.future_vectors.as_ref().is_some_and(|f| {
            f.vectors.lengths.len() != b || f.vectors.lengths.iter().any(|&x| x != k)
        })
    {
        return Err("INVALID_BATCH_SHAPE".into());
    }
    // Hard request bounds precede allocation of borrowed ABI descriptors.
    // Smaller original quotas are adjudicated by the single native call.
    let mut elements = 0usize;
    for len in [
        r.noise.covariance_row_major.len(),
        r.design.values_row_major.len(),
        r.parameter_prior.mean.len(),
        r.parameter_prior.covariance_row_major.len(),
        r.conditioning.vectors.values.len(),
        r.future_noise.covariance_row_major.len(),
        r.future_response.values_row_major.len(),
        r.future_vectors
            .as_ref()
            .map_or(0, |f| f.vectors.values.len()),
    ] {
        elements = elements.checked_add(len).ok_or("REQUEST_ELEMENT_LIMIT")?;
    }
    if r.outputs.means {
        elements = elements
            .checked_add(
                b.checked_mul(k)
                    .and_then(|x| x.checked_mul(2))
                    .ok_or("REQUEST_ELEMENT_LIMIT")?,
            )
            .ok_or("REQUEST_ELEMENT_LIMIT")?;
    }
    if r.outputs.joint_log_densities {
        elements = elements
            .checked_add(b.checked_mul(6).ok_or("REQUEST_ELEMENT_LIMIT")?)
            .ok_or("REQUEST_ELEMENT_LIMIT")?;
    }
    if elements > 1000000 {
        return Err("REQUEST_ELEMENT_LIMIT".into());
    }
    let mut text_bytes = 0usize;
    for ids in [
        &r.noise.ordered_row_ids,
        &r.noise.event_ids,
        &r.design.ordered_row_ids,
        &r.design.ordered_parameter_ids,
        &r.design.parameter_units,
        &r.design.column_units,
        &r.parameter_prior.ordered_parameter_ids,
        &r.parameter_prior.parameter_units,
        &r.parameter_prior.shared_nuisance_ids,
        &r.conditioning.ordered_row_ids,
        &r.conditioning.event_ids,
        &r.conditioning.case_ids,
        &r.future_noise.ordered_row_ids,
        &r.future_noise.event_ids,
        &r.future_response.ordered_row_ids,
        &r.future_response.ordered_parameter_ids,
        &r.future_response.parameter_units,
        &r.future_response.column_units,
    ] {
        for text in ids.iter() {
            text_bytes = text_bytes
                .checked_add(text.len())
                .ok_or("REQUEST_TEXT_LIMIT")?;
        }
    }
    if let Some(f) = &r.future_vectors {
        for ids in [&f.ordered_row_ids, &f.event_ids] {
            for text in ids.iter() {
                text_bytes = text_bytes
                    .checked_add(text.len())
                    .ok_or("REQUEST_TEXT_LIMIT")?;
            }
        }
    }
    for text in [
        &r.noise.residual_unit,
        &r.noise.noise_identity,
        &r.noise.calibration_identity,
        &r.noise.dependence_identity,
        &r.noise.ordering_provenance,
        &r.design.design_identity,
        &r.parameter_prior.prior_identity,
        &r.parameter_prior.parameter_measure,
        &r.parameter_prior.dependence_identity,
        &r.future_noise.residual_unit,
        &r.future_noise.noise_identity,
        &r.future_noise.calibration_identity,
        &r.future_noise.dependence_identity,
        &r.future_noise.ordering_provenance,
        &r.future_response.design_identity,
        &r.prediction.future_covariance_unit,
        &r.prediction.future_measure,
        &r.prediction.conditioning_identity,
        &r.prediction.dependence_identity,
    ] {
        text_bytes = text_bytes
            .checked_add(text.len())
            .ok_or("REQUEST_TEXT_LIMIT")?;
    }
    if text_bytes > 1048576 {
        return Err("REQUEST_TEXT_LIMIT".into());
    }
    let c = crate::bridge::gaussian_predictive::evaluate(&r)?;
    let mut spec = serde_json::to_value(&r).map_err(|_| "RECORD_ENCODING")?;
    let requested: Vec<_> = [
        ("predictive_means", r.outputs.means, c.means_passed),
        (
            "joint_predictive_log_densities",
            r.outputs.joint_log_densities,
            c.densities_passed,
        ),
    ]
    .into_iter()
    .filter(|(_, requested, _)| *requested)
    .collect();
    spec["requested_outputs"]=json!(requested.iter().map(|(id,_,_)|json!({"id":id,"required":true,"numerical_gate":"required","inference_gate":"not_applicable"})).collect::<Vec<_>>());
    spec["batch_meaning"] = json!(
        "ordered family of conditionals under one fixed synthetic model; no IID or joint-across-cases law"
    );
    let mut outcome = Outcome::scientific(
        spec,
        c.output,
        requested[0].0,
        c.method,
        c.arithmetic,
        c.resources,
        "inherited native joint Gaussian controls; generic fixed synthetic responses only",
    );
    outcome.outputs=requested.into_iter().map(|(id,_,pass)|OutputCheck{check_kind:"numerical_contract",id,required:true,numerical:if pass{NumericalCheck::ChecksPassed}else{NumericalCheck::Failed},inference:"not_applicable",interpretation:"unqualified",evidence:vec![],validation_coverage:"original rational/cofactor and direct-product quadrature ancestry; native/ABI/CLI parity is interface evidence"}).collect();
    Ok(outcome)
}
