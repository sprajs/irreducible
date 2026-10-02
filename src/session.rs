//! Current process-local session. Complete replies precede handle publication.
use crate::{
    observation_run::{self, Acquired},
    retained_context::{Context, Error, Limits},
    stream_io::{self, Line},
};
use serde::Deserialize;
use serde_json::{Value, json};
use std::{
    io::{BufRead, Write},
    path::Path,
};
#[derive(Deserialize)]
#[serde(tag = "action", rename_all = "snake_case", deny_unknown_fields)]
enum Command {
    PrepareObservations {
        request: Value,
    },
    PrepareSupernova {
        source_handle: u64,
        selection: crate::bridge::supernova::Selection,
        preparation_policy: crate::bridge::supernova::PreparationPolicy,
    },
    PrepareBao {
        source: crate::bao_ingestion::Source,
        preparation_policy: crate::bridge::bao::PreparationPolicy,
    },
    BaoDensity {
        handle: u64,
        models: Vec<crate::bridge::bao::Model>,
        requested_outputs: Vec<crate::bridge::bao::Output>,
        numerical_policy: crate::bridge::bao::EvaluationPolicy,
    },
    PhotometryPredict {
        request: Value,
    },
    SoundHorizonEvaluate {
        request: Value,
    },
    BackgroundEvaluate {
        request: Value,
    },
    SupernovaProfile {
        handle: u64,
        models: Vec<crate::model_spec::SupernovaModel>,
        requested_outputs: Vec<crate::bridge::supernova::Output>,
        numerical_policy: crate::bridge::supernova::EvaluationPolicy,
    },
    Release {
        handle: u64,
    },
}
// Consumer variants are added when their native retained bridges freeze.
pub(crate) enum Consumer {
    Supernova {
        owner: crate::bridge::supernova::Prepared,
        specification: Value,
    },
    Bao {
        owner: crate::bridge::bao::Prepared,
        specification: Value,
    },
}
fn error_frame(error: &str, limit: usize) -> Result<Vec<u8>, String> {
    // Parser diagnostics can contain a user-controlled field name. Keep the
    // rejection frame bounded independently of that text; raw admitted input
    // is retained by the caller's record path rather than copied into errors.
    let error = if error.len() > 128 {
        "REQUEST_REJECTED"
    } else {
        error
    };
    stream_io::encode(
        &json!({"execution":"failed","error_id":error,"accepted":false}),
        limit,
    )
}
fn outcome_reply(outcome: crate::outcome::Outcome) -> Value {
    let passed = outcome.numerical_passed();
    json!({"execution":outcome.execution,"specification":outcome.specification,
        "output":outcome.output,"method":outcome.method,"arithmetic":outcome.arithmetic,
        "resources":outcome.resources,"outputs":outcome.outputs,
        "accepted":passed,"accepted_scope":"numerical_contract","numerical":if passed{"checks_passed"}else{"failed"},
        "interpretation":"unqualified"})
}
pub(crate) fn run<R: BufRead, W: Write>(
    reader: &mut R,
    writer: &mut W,
    store: &Path,
    limits: Limits,
) -> Result<(), String> {
    if limits.output_line_bytes < 1024 || limits.input_line_bytes == 0 || limits.commands == 0 {
        return Err("INVALID_SESSION_LIMITS".into());
    }
    std::fs::create_dir_all(store.join("objects")).map_err(|_| "STORE_IO")?;
    let mut context: Context<Acquired, Consumer> = Context::new(limits);
    let build: Value = serde_json::from_str(include_str!(env!("IRRED_BUILD_MANIFEST")))
        .map_err(|_| "INVALID_BUILD_MANIFEST")?;
    let executable_digest = crate::records::hash(
        &std::fs::read(std::env::current_exe().map_err(|_| "EXECUTABLE_IO")?)
            .map_err(|_| "EXECUTABLE_IO")?,
    );
    let configuration = json!({"limits":limits,"assurance":"numerical_contract","backend":"portable_cpu","compute_threads":1,"io_threads":1,"rng":"not_applicable"});
    let configuration_bytes = serde_json::to_vec(&configuration).map_err(|_| "RECORD_ENCODING")?;
    let configuration_digest = crate::records::hash(&configuration_bytes);
    crate::records::publish(
        &store.join("objects").join(&configuration_digest),
        &configuration_bytes,
    )?;
    let mut admitted = 0u64;
    loop {
        // Exactly one reply per admitted request. The session closes normally
        // after its Nth reply, without reading or replying to an N+1 frame.
        if admitted == limits.commands {
            return Ok(());
        }
        let line = stream_io::read_line(reader, limits.input_line_bytes).map_err(|_| "INPUT_IO")?;
        if matches!(line, Line::End) {
            return Ok(());
        }
        admitted += 1;
        let input_digest = match &line {
            Line::Complete(bytes) => Some(crate::records::hash(bytes)),
            _ => None,
        };
        let mut receipt = json!({"schema_version":2,"build_id":build["build_id"],"executable_digest":executable_digest,
            "execution_configuration_digest":configuration_digest,"input_digest":input_digest,"output_digest":"0".repeat(64),"command_index":admitted,
            "requested_assurance":"numerical_contract","record_scope":"execution_binding; disposition and output checks in result"});
        let receipt_bytes = serde_json::to_vec(&receipt).map_err(|_| "RECORD_ENCODING")?;
        // Reserve the exact fixed framing/receipt length BEFORE transactional
        // dispatch. All substituted digests have the same encoded length.
        let overhead = receipt_bytes
            .len()
            .checked_add(24)
            .ok_or("OUTPUT_LINE_LIMIT")?;
        let payload_limit = limits
            .output_line_bytes
            .checked_sub(overhead)
            .filter(|n| *n >= 128)
            .ok_or("INVALID_SESSION_LIMITS")?;
        let reply = match line {
            Line::TooLarge => error_frame("INPUT_LINE_LIMIT", payload_limit),
            Line::Complete(bytes) => {
                let input_digest = crate::records::hash(&bytes);
                crate::records::publish(&store.join("objects").join(input_digest), &bytes)?;
                dispatch(&mut context, &bytes, store, payload_limit)
            }
            Line::End => unreachable!(),
        };
        let frame = match reply {
            Ok(frame) => frame,
            Err(error) => error_frame(&error, payload_limit)?,
        };
        let output_digest = crate::records::hash(&frame);
        crate::records::publish(&store.join("objects").join(&output_digest), &frame)?;
        receipt["output_digest"] = json!(output_digest);
        let receipt_bytes = serde_json::to_vec(&receipt).map_err(|_| "RECORD_ENCODING")?;
        let receipt_digest = crate::records::hash(&receipt_bytes);
        crate::records::publish(&store.join("objects").join(receipt_digest), &receipt_bytes)?;
        let mut encoded = Vec::with_capacity(overhead + frame.len());
        encoded.extend_from_slice(b"{\"receipt\":");
        encoded.extend_from_slice(&receipt_bytes);
        encoded.extend_from_slice(b",\"result\":");
        encoded.extend_from_slice(frame.strip_suffix(b"\n").unwrap_or(&frame));
        encoded.extend_from_slice(b"}\n");
        debug_assert!(encoded.len() <= limits.output_line_bytes);
        writer.write_all(&encoded).map_err(|_| "OUTPUT_IO")?;
        writer.flush().map_err(|_| "OUTPUT_IO")?;
    }
}
fn dispatch(
    context: &mut Context<Acquired, Consumer>,
    bytes: &[u8],
    store: &Path,
    limit: usize,
) -> Result<Vec<u8>, String> {
    crate::strict_json::validate(bytes).map_err(|_| "INVALID_SESSION_COMMAND")?;
    let command: Command = serde_json::from_slice(bytes).map_err(|_| "INVALID_SESSION_COMMAND")?;
    match command {
        Command::PrepareBao {
            source,
            preparation_policy,
        } => {
            let mut effective = preparation_policy.clone();
            effective.maximum_native_bytes = effective
                .maximum_native_bytes
                .min(u64::try_from(context.remaining_bytes()).map_err(|_| "RESOURCE_LIMIT")?);
            let peak = context
                .reserve_peak(
                    usize::try_from(effective.maximum_native_bytes)
                        .map_err(|_| "RESOURCE_LIMIT")?,
                )
                .map_err(|_| "RETAINED_BYTE_LIMIT")?;
            let acquired = crate::bao_run::acquire(
                source,
                &preparation_policy,
                effective.clone(),
                store,
                peak.bytes(),
            )?;
            let reply = json!({"execution":"completed","preparation_digest":acquired.preparation_digest,"specification":acquired.specification,"effective_runtime_policy":acquired.effective_policy,
                "status":acquired.status,"numerical_status":acquired.numerical_status,"accepted":acquired.status==0&&acquired.numerical_status==0,"accepted_scope":"numerical_contract"});
            let id = context
                .insert_independent_consumer_reserved(
                    Consumer::Bao {
                        owner: acquired.owner,
                        specification: acquired.specification,
                    },
                    acquired.preparation_digest,
                    acquired.retained_bytes,
                    peak,
                )
                .map_err(|_| "RETAINED_BYTE_LIMIT")?;
            context
                .publish_reply(id, |id| {
                    let mut reply = reply;
                    reply["handle"] = json!(id);
                    stream_io::encode(&reply, limit).map_err(|_| Error::OutputLine)
                })
                .map_err(|_| "OUTPUT_LINE_LIMIT".into())
        }
        Command::BaoDensity {
            handle,
            models,
            requested_outputs,
            mut numerical_policy,
        } => {
            let requested_policy = numerical_policy.clone();
            numerical_policy.maximum_native_bytes = numerical_policy
                .maximum_native_bytes
                .min(u64::try_from(context.remaining_bytes()).map_err(|_| "RESOURCE_LIMIT")?);
            let consumer = context.consumer(handle).map_err(|error| match error {
                Error::WrongKind => "WRONG_HANDLE_KIND",
                _ => "UNKNOWN_HANDLE",
            })?;
            let Consumer::Bao {
                owner,
                specification,
            } = &consumer.value
            else {
                return Err("WRONG_HANDLE_KIND".into());
            };
            let (output, checks) = crate::bridge::bao::evaluate(
                owner,
                &models,
                &requested_outputs,
                &numerical_policy,
            )?;
            let spec = json!({"operation":"bao.density","preparation_digest":consumer.preparation_digest,"preparation":specification,"models":models,"requested_outputs":requested_outputs,"numerical_policy":requested_policy});
            let encoded = serde_json::to_vec(&spec).map_err(|_| "RECORD_ENCODING")?;
            let digest = crate::records::hash(&encoded);
            crate::records::publish(&store.join("objects").join(&digest), &encoded)?;
            let encoded = serde_json::to_vec(&output).map_err(|_| "RECORD_ENCODING")?;
            let output_digest = crate::records::hash(&encoded);
            crate::records::publish(&store.join("objects").join(&output_digest), &encoded)?;
            let passed = checks.iter().all(|(_, p)| *p);
            stream_io::encode(
                &json!({"execution":"completed","scientific_specification_digest":digest,"output_digest":output_digest,"output":output,"method":"normalized_conditional_free_ruler_gaussian","arithmetic":requested_policy.arithmetic,"effective_runtime_policy":numerical_policy,
                "outputs":checks.iter().map(|(id,p)|json!({"id":id,"required":true,"numerical":if *p{"checks_passed"}else{"failed"}})).collect::<Vec<_>>(),"accepted":passed,"accepted_scope":"numerical_contract","interpretation":"unqualified"}),
                limit,
            )
        }
        Command::SupernovaProfile {
            handle,
            models,
            requested_outputs,
            mut numerical_policy,
        } => {
            let requested_policy = numerical_policy.clone();
            numerical_policy.maximum_native_bytes = numerical_policy
                .maximum_native_bytes
                .min(u64::try_from(context.remaining_bytes()).map_err(|_| "RESOURCE_LIMIT")?);
            let consumer = context.consumer(handle).map_err(|error| match error {
                Error::WrongKind => "WRONG_HANDLE_KIND",
                _ => "UNKNOWN_HANDLE",
            })?;
            let Consumer::Supernova {
                owner,
                specification,
            } = &consumer.value
            else {
                return Err("WRONG_HANDLE_KIND".into());
            };
            let (output, checks) = crate::bridge::supernova::evaluate(
                owner,
                &models,
                &requested_outputs,
                &numerical_policy,
            )?;
            let spec = json!({"operation":"supernova.profile","preparation_digest":consumer.preparation_digest,
                "preparation":specification,"models":models,"requested_outputs":requested_outputs,"numerical_policy":requested_policy});
            let encoded = serde_json::to_vec(&spec).map_err(|_| "RECORD_ENCODING")?;
            let spec_digest = crate::records::hash(&encoded);
            crate::records::publish(&store.join("objects").join(&spec_digest), &encoded)?;
            let encoded = serde_json::to_vec(&output).map_err(|_| "RECORD_ENCODING")?;
            let output_digest = crate::records::hash(&encoded);
            crate::records::publish(&store.join("objects").join(&output_digest), &encoded)?;
            let passed = checks.iter().all(|(_, p)| *p);
            stream_io::encode(
                &json!({"execution":"completed","scientific_specification_digest":spec_digest,
                "output_digest":output_digest,"output":output,"method":"conditional_single_offset_profile","arithmetic":requested_policy.arithmetic,"effective_runtime_policy":numerical_policy,"outputs":checks.iter().map(|(id,p)|json!({"id":id,"required":true,"numerical":if *p{"checks_passed"}else{"failed"}})).collect::<Vec<_>>(),
                "accepted":passed,"accepted_scope":"numerical_contract","interpretation":"unqualified"}),
                limit,
            )
        }
        Command::PrepareSupernova {
            source_handle,
            selection,
            mut preparation_policy,
        } => {
            let requested_policy = preparation_policy.clone();
            let source = context.source(source_handle).map_err(|error| match error {
                Error::WrongKind => "WRONG_HANDLE_KIND",
                _ => "UNKNOWN_HANDLE",
            })?;
            let remaining =
                u64::try_from(context.remaining_bytes()).map_err(|_| "RESOURCE_LIMIT")?;
            preparation_policy.maximum_native_bytes =
                preparation_policy.maximum_native_bytes.min(remaining);
            let peak = context
                .reserve_peak(
                    usize::try_from(preparation_policy.maximum_native_bytes)
                        .map_err(|_| "RESOURCE_LIMIT")?,
                )
                .map_err(|_| "RETAINED_BYTE_LIMIT")?;
            let source_spec_bytes = serde_json::to_vec(&source.value.specification)
                .map_err(|_| "RECORD_ENCODING")?
                .len();
            let selected_request_bytes = serde_json::to_vec(&selection)
                .map_err(|_| "RECORD_ENCODING")?
                .len();
            let metadata_peak = source_spec_bytes
                .checked_add(selected_request_bytes)
                .and_then(|n| n.checked_mul(128))
                .and_then(|n| {
                    n.checked_add(
                        usize::try_from(preparation_policy.maximum_string_bytes)
                            .ok()?
                            .checked_mul(64)?,
                    )
                })
                .and_then(|n| n.checked_add(4096))
                .ok_or("RESOURCE_LIMIT")?;
            let native_bytes = peak
                .bytes()
                .checked_sub(metadata_peak)
                .ok_or("PREPARATION_BYTE_LIMIT")?;
            preparation_policy.maximum_native_bytes =
                u64::try_from(native_bytes).map_err(|_| "RESOURCE_LIMIT")?;
            let prepared = crate::bridge::supernova::prepare(
                &source.value.owner,
                &selection,
                &preparation_policy,
            )?;
            let spec = json!({"operation":"supernova.profile","preparation":{"source_preparation_digest":source.preparation_digest,
                "source_specification":source.value.specification,"selection":prepared.selected,"policy":requested_policy}});
            let encoded = serde_json::to_vec(&spec).map_err(|_| "RECORD_ENCODING")?;
            let digest = crate::records::hash(&encoded);
            crate::records::publish(&store.join("objects").join(&digest), &encoded)?;
            // Additional context metadata/container/bookkeeping envelope, not
            // an allocator/RSS promise. Source matrix remains charged once.
            let bytes = prepared
                .retained_bytes
                .checked_add(encoded.len().checked_mul(64).ok_or("RESOURCE_LIMIT")?)
                .and_then(|n| n.checked_add(4096))
                .ok_or("RESOURCE_LIMIT")?;
            let reply = json!({"execution":"completed","preparation_digest":digest,"specification":spec,"effective_runtime_policy":preparation_policy,
                "preparation_status":prepared.status,"numerical_status":prepared.numerical_status,
                "accepted":prepared.status==0,"accepted_scope":"numerical_contract"});
            let id = context
                .insert_consumer_reserved(
                    source,
                    Consumer::Supernova {
                        owner: prepared.owner,
                        specification: spec,
                    },
                    digest,
                    bytes,
                    peak,
                )
                .map_err(|_| "RETAINED_BYTE_LIMIT")?;
            context
                .publish_reply(id, |id| {
                    let mut reply = reply;
                    reply["handle"] = json!(id);
                    stream_io::encode(&reply, limit).map_err(|_| Error::OutputLine)
                })
                .map_err(|_| "OUTPUT_LINE_LIMIT".into())
        }
        Command::PrepareObservations { request } => {
            let input = serde_json::to_vec(&request).map_err(|_| "INVALID_REQUEST")?;
            let peak_bytes = observation_run::peak_bound(&input)?;
            let peak = context
                .reserve_peak(peak_bytes)
                .map_err(|_| "RETAINED_BYTE_LIMIT")?;
            let acquired = observation_run::acquire(&input, store, peak.bytes())?;
            let spec_bytes =
                serde_json::to_vec(&acquired.specification).map_err(|_| "RECORD_ENCODING")?;
            crate::records::publish(
                &store.join("objects").join(&acquired.preparation_digest),
                &spec_bytes,
            )?;
            let output_bytes =
                serde_json::to_vec(&acquired.output).map_err(|_| "RECORD_ENCODING")?;
            let output_digest = crate::records::hash(&output_bytes);
            crate::records::publish(&store.join("objects").join(&output_digest), &output_bytes)?;
            let reply = json!({"execution":"completed","preparation_digest":acquired.preparation_digest,
                "specification":acquired.specification,"output":acquired.output,"output_digest":output_digest,
                "accepted":true,"accepted_scope":"structural_source_contract"});
            let digest = acquired.preparation_digest.clone();
            let charge = acquired.retained_bytes;
            let id = context
                .insert_source_reserved(acquired, digest, charge, peak)
                .map_err(|_| "HANDLE_LIMIT")?;
            context
                .publish_reply(id, |id| {
                    let mut reply = reply;
                    reply["handle"] = json!(id);
                    stream_io::encode(&reply, limit).map_err(|_| Error::OutputLine)
                })
                .map_err(|_| "OUTPUT_LINE_LIMIT".into())
        }
        Command::PhotometryPredict { mut request } => {
            let requested_policy = request["resource_policy"].clone();
            let requested = request["resource_policy"]["maximum_native_bytes"]
                .as_u64()
                .ok_or("INVALID_REQUEST")?;
            let remaining =
                u64::try_from(context.remaining_bytes()).map_err(|_| "RESOURCE_LIMIT")?;
            request["resource_policy"]["maximum_native_bytes"] = json!(requested.min(remaining));
            let effective_policy = request["resource_policy"].clone();
            let input = serde_json::to_vec(&request).map_err(|_| "INVALID_REQUEST")?;
            let mut outcome = crate::photometry_run::execute(&input)?;
            outcome.specification["resource_policy"] = requested_policy;
            let mut reply = outcome_reply(outcome);
            reply["effective_runtime_policy"] = effective_policy;
            stream_io::encode(&reply, limit)
        }
        Command::SoundHorizonEvaluate { mut request } => {
            let requested_policy = request["numerical_policy"].clone();
            // Bound transient evaluation against currently retained owners.
            // Effective context allowance is execution evidence, not scientific identity.
            let requested = request["numerical_policy"]["maximum_native_bytes"]
                .as_u64()
                .ok_or("INVALID_REQUEST")?;
            let remaining =
                u64::try_from(context.remaining_bytes()).map_err(|_| "RESOURCE_LIMIT")?;
            request["numerical_policy"]["maximum_native_bytes"] = json!(requested.min(remaining));
            let input = serde_json::to_vec(&request).map_err(|_| "INVALID_REQUEST")?;
            let effective_policy = request["numerical_policy"].clone();
            let mut outcome = crate::sound_horizon_run::execute(&input)?;
            outcome.specification["numerical_policy"] = requested_policy;
            let mut reply = outcome_reply(outcome);
            reply["effective_runtime_policy"] = effective_policy;
            stream_io::encode(&reply, limit)
        }
        Command::BackgroundEvaluate { mut request } => {
            let requested_policy = request["numerical_policy"].clone();
            // Bound transient evaluation against currently retained owners.
            // Effective context allowance is execution evidence, not scientific identity.
            let requested = request["numerical_policy"]["maximum_native_bytes"]
                .as_u64()
                .ok_or("INVALID_REQUEST")?;
            let remaining =
                u64::try_from(context.remaining_bytes()).map_err(|_| "RESOURCE_LIMIT")?;
            request["numerical_policy"]["maximum_native_bytes"] = json!(requested.min(remaining));
            let input = serde_json::to_vec(&request).map_err(|_| "INVALID_REQUEST")?;
            let effective_policy = request["numerical_policy"].clone();
            let mut outcome = crate::background_run::execute(&input)?;
            outcome.specification["numerical_policy"] = requested_policy;
            let mut reply = outcome_reply(outcome);
            reply["effective_runtime_policy"] = effective_policy;
            stream_io::encode(&reply, limit)
        }
        Command::Release { handle } => context
            .release_reply(handle, |id| {
                stream_io::encode(
                    &json!({"execution":"completed","released_handle":id,"disposition":"released"}),
                    limit,
                )
                .map_err(|_| Error::OutputLine)
            })
            .map_err(|error| {
                match error {
                    Error::UnknownHandle => "UNKNOWN_HANDLE",
                    Error::OutputLine => "OUTPUT_LINE_LIMIT",
                    _ => "HANDLE_ERROR",
                }
                .into()
            }),
    }
}

// One-shot callers use exactly the retained dispatch path and then drop the
// context. No second reader, observation preparation, or likelihood factory.
pub(crate) fn execute_once(input: &[u8], store: &Path) -> Result<crate::outcome::Outcome, String> {
    #[derive(Deserialize)]
    #[serde(deny_unknown_fields)]
    struct Request {
        schema_version: u32,
        operation: String,
        observations: Option<Value>,
        source: Option<Value>,
        selection: Option<Value>,
        preparation_policy: Value,
        models: Value,
        requested_outputs: Value,
        numerical_policy: Value,
    }
    let request: Request = serde_json::from_slice(input).map_err(|_| "INVALID_REQUEST")?;
    if request.schema_version != 2 {
        return Err("UNSUPPORTED_SPECIFICATION".into());
    }
    let arithmetic = request.preparation_policy["arithmetic"].clone();
    let limits = Limits::default();
    let mut context = Context::new(limits);
    let mut call = |command: Value| -> Result<Value, String> {
        let bytes = serde_json::to_vec(&command).map_err(|_| "REQUEST_ENCODING")?;
        let reply = dispatch(&mut context, &bytes, store, limits.output_line_bytes)?;
        serde_json::from_slice(&reply).map_err(|_| "REPLY_ENCODING".into())
    };
    let (prepared, action, method) = match request.operation.as_str() {
        "supernova.profile" => {
            if request.source.is_some() {
                return Err("INCOMPATIBLE_SOURCE_FIELDS".into());
            }
            let source = call(
                json!({"action":"prepare_observations","request":request.observations.ok_or("MISSING_OBSERVATIONS")?}),
            )?;
            let prepared = call(
                json!({"action":"prepare_supernova","source_handle":source["handle"],
                "selection":request.selection.ok_or("MISSING_SELECTION")?,"preparation_policy":request.preparation_policy}),
            )?;
            (
                prepared,
                "supernova_profile",
                "conditional_single_offset_profile",
            )
        }
        "bao.density" => {
            if request.observations.is_some() || request.selection.is_some() {
                return Err("INCOMPATIBLE_SOURCE_FIELDS".into());
            }
            let prepared = call(
                json!({"action":"prepare_bao","source":request.source.ok_or("MISSING_SOURCE")?,
                "preparation_policy":request.preparation_policy}),
            )?;
            (
                prepared,
                "bao_density",
                "normalized_conditional_free_ruler_gaussian",
            )
        }
        _ => return Err("UNSUPPORTED_SPECIFICATION".into()),
    };
    let reply = call(
        json!({"action":action,"handle":prepared["handle"],"models":request.models,
        "requested_outputs":request.requested_outputs,"numerical_policy":request.numerical_policy}),
    )?;
    let digest = reply["scientific_specification_digest"]
        .as_str()
        .ok_or("MISSING_SPECIFICATION")?;
    let specification: Value = serde_json::from_slice(
        &std::fs::read(store.join("objects").join(digest)).map_err(|_| "STORE_IO")?,
    )
    .map_err(|_| "INVALID_STORED_SPECIFICATION")?;
    let mut outcome = crate::outcome::Outcome::scientific(
        specification,
        reply["output"].clone(),
        "batch",
        json!(method),
        arithmetic,
        json!({"effective_runtime_policy":reply["effective_runtime_policy"]}),
        "named native reference cases; arbitrary request interpretation remains unqualified",
    );
    outcome.outputs = reply["outputs"].as_array().ok_or("INVALID_OUTPUT_CHECKS")?.iter().map(|check| {
        let id = match check["id"].as_str() {
            Some("score")=>"score",Some("geometric_shape")=>"geometric_shape",Some("magnitude_effect")=>"magnitude_effect",
            Some("corrected_residuals")=>"corrected_residuals",Some("profiled_residuals")=>"profiled_residuals",Some("diagnostics")=>"diagnostics",
            Some("normalized_density")=>"normalized_density",Some("predictions")=>"predictions",Some("residuals")=>"residuals",
            _=>return Err("INVALID_OUTPUT_CHECK_ID"),
        };
        Ok(crate::outcome::OutputCheck{check_kind:"numerical_contract",id,required:true,
            numerical:if check["numerical"]=="checks_passed"{crate::outcome::NumericalCheck::ChecksPassed}else{crate::outcome::NumericalCheck::Failed},
            inference:"not_applicable",interpretation:"unqualified",evidence:vec![],validation_coverage:"named native reference cases; request applicability not established"})
    }).collect::<Result<Vec<_>,_>>()?;
    Ok(outcome)
}
