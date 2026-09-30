use crate::{
    bridge::{
        ABI_VERSION, MAX_BATCH_ELEMENTS, Metadata, add, convert_quantities, numerics_evaluate,
    },
    records::{hash, publish, runtime_libraries},
};
use serde::Deserialize;
use serde_json::{Value, json};
use std::{fs, io::Read, path::PathBuf};
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Request {
    schema_version: u32,
    operation: String,
    a: Vec<i64>,
    b: Vec<i64>,
    #[serde(default)]
    fault: u32,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct QuantityRequest {
    schema_version: u32,
    operation: String,
    values: Vec<f64>,
    source: Metadata,
    target: Metadata,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct NumericalRequest {
    schema_version: u32,
    operation: String,
    method: String,
    values: Vec<f64>,
}
#[derive(Deserialize)]
struct OperationHeader {
    operation: String,
}
pub(crate) fn execute() -> Result<(), String> {
    let args: Vec<_> = std::env::args().collect();
    let manifest: Value = serde_json::from_str(include_str!(env!("IRRED_BUILD_MANIFEST")))
        .map_err(|e| e.to_string())?;
    match args.get(1).map(String::as_str) {
        Some("stream") if args.len() == 3 || args.len() == 4 => {
            let store = PathBuf::from(&args[2]);
            let limits = if let Some(path) = args.get(3) {
                serde_json::from_slice(&fs::read(path).map_err(|_| "LIMITS_IO")?).map_err(|_| "INVALID_SESSION_LIMITS")?
            } else { crate::retained_context::Limits::default() };
            crate::current_session::run(&mut std::io::stdin().lock(), &mut std::io::stdout().lock(), &store, limits)
        }
        Some("describe") | Some("version") if args.len() == 3 && args[2] == "--json" => {
            println!(
                "{}",
                json!({"schema_version":2,"product":"Irreducible","executable":"irred","version":env!("CARGO_PKG_VERSION"),"abi_version":ABI_VERSION,"build":manifest,
                    "capabilities":[
                        {"id":"fixture.checked_i64_add.v1","implementation":"implemented","scientific":false,"qualification":"unqualified"},
                        {"id":"quantity.convert.v1","implementation":"implemented","scientific":true,"qualification":"unqualified"},
                        {"id":"numerics.scalar_batch.v1","implementation":"implemented","scientific":true,"qualification":"unqualified"},
                        {"id":"background.evaluate","implementation":"implemented","scientific":true,"qualification":"unqualified","requested_groups":["radial","luminosity_shape","clock","physical","kinematics","expansion"],"models":["lcdm","constant_q","cpl","fixed_q5"],"node_reuse":"exact z bits within each model batch; no cross-model cache"},
                        {"id":"observations.prepare","implementation":"implemented","scientific":true,"qualification":"unqualified","profiles":["pantheon_plus_released_v1","gaussian_fixture_v1","typed_magnitude_covariance"],"ownership":"immutable shared native source"},
                        {"id":"statistics.gaussian","implementation":"implemented","scientific":true,"qualification":"unqualified","modes":["normalized_density","profile_offset_score"]},
                        {"id":"supernova.profile","implementation":"implemented","scientific":true,"qualification":"unqualified","models":["lcdm","constant_q","cpl","fixed_q5"],"source_effects":["none","grey_log1p_magnitude"],"target":"conditional single-offset relative profile"},
                        {"id":"bao.density","implementation":"implemented","scientific":true,"qualification":"unqualified","models":["lcdm","constant_q","cpl","fixed_q5"],"target":"conditional free-H0rd normalized Gaussian density"}],
                    "abi_schema":serde_json::from_str::<Value>(include_str!("../schema/abi.json")).map_err(|e|e.to_string())?,
                    "commands":["describe --json","version --json","run REQUEST STORE [--assurance numerical_contract|qualified]","stream STORE [LIMITS_JSON]"],
                    "interface_policy":"one current ABI revision; no compatibility aliases",
                    "scientific_qualifications":[]})
            );
            Ok(())
        }
        Some("run") if args.len()==4 || (args.len()==6 && args[4]=="--assurance" && matches!(args[5].as_str(),"numerical_contract"|"qualified")) => {
            let assurance=args.get(5).map(String::as_str).unwrap_or("numerical_contract");
            let mut input = Vec::new();
            fs::File::open(args.get(2).ok_or("missing request path")?).map_err(|_| "INPUT_IO")?
                .take((16 << 20) + 1).read_to_end(&mut input).map_err(|_| "INPUT_IO")?;
            if input.len() > 16 << 20 { return Err("INPUT_BYTE_LIMIT".into()); }
            let store = PathBuf::from(args.get(3).ok_or("missing store path")?);
            fs::create_dir_all(store.join("objects")).map_err(|e| e.to_string())?;
            fs::create_dir_all(store.join("attempts")).map_err(|e| e.to_string())?;
            let attempt = format!(
                "{}-{}",
                std::process::id(),
                std::time::SystemTime::now()
                    .duration_since(std::time::UNIX_EPOCH)
                    .unwrap()
                    .as_nanos()
            );
            let executable = fs::read(std::env::current_exe().map_err(|e| e.to_string())?)
                .map_err(|e| e.to_string())?;
            let input_hash = hash(&input);
            publish(&store.join("objects").join(&input_hash), &input)?;
            let initial = json!({"schema_version":1,"attempt_id":attempt,"execution":"incomplete","input_digest":input_hash,"build_id":manifest["build_id"],"source_revision":manifest["git_head"],"source_status":manifest["git_status"],"executable_digest":hash(&executable),"numerical":"not_assessed","inference":"not_applicable","interpretation":"not_assessed","runtime_libraries":runtime_libraries()?,"resource_budget":{"compute_threads":1,"io_threads":1},"rng":"not_applicable","backend":"portable_cpu","precision":"unresolved","fixture_fault_controls":{"abort":std::env::var_os("COSMOLOGY_TEST_ABORT").is_some(),"panic":std::env::var_os("COSMOLOGY_TEST_PANIC").is_some(),"abort_after_output":std::env::var_os("COSMOLOGY_TEST_ABORT_AFTER_OUTPUT").is_some()}});
            publish(
                &store
                    .join("attempts")
                    .join(format!("{attempt}.incomplete.json")),
                &serde_json::to_vec(&initial).unwrap(),
            )?;
            if std::env::var_os("COSMOLOGY_TEST_ABORT").is_some() {
                std::process::abort()
            }
            let mut resolved = None;
            let mut fault = 0;
            let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(
                || -> Result<crate::outcome::Outcome, String> {
                    crate::strict_json::validate(&input).map_err(|e| e.to_string())?;
                    let header: OperationHeader = serde_json::from_slice(&input).map_err(|e|e.to_string())?;
                    match header.operation.as_str() {
                        "fixture.checked_i64_add.v1" => {
                            let request: Request = serde_json::from_slice(&input).map_err(|e|e.to_string())?;
                            if request.schema_version != 2 { return Err("UNSUPPORTED_SPECIFICATION".into()); }
                            fault = request.fault;
                            resolved = Some(json!({"schema_version":2,"operation":request.operation,
                                "equation_id":"EQ-fixture.checked_i64_add.v1","a":request.a,"b":request.b,
                                "requested_outputs":[{"id":"sum","required":true,"numerical_gate":"required","inference_gate":"not_applicable"}]}));
                            if std::env::var_os("COSMOLOGY_TEST_PANIC").is_some() {panic!("injected Rust panic");}
                            add(&request.a,&request.b,request.fault).map(|v|crate::outcome::Outcome::exact(resolved.clone().unwrap(),json!({"kind":"finite","values":v})))
                        },
                        "quantity.convert.v1" => {
                            let request: QuantityRequest = serde_json::from_slice(&input).map_err(|e|e.to_string())?;
                            if request.schema_version != 2 { return Err("UNSUPPORTED_SPECIFICATION".into()); }
                            resolved = Some(json!({"schema_version":2,"operation":request.operation,
                                "equation_id":"EQ-quantity-conversion-v1","values":request.values,
                                "source":request.source,"target":request.target,
                                "requested_outputs":[{"id":"converted","required":true,"numerical_gate":"required","inference_gate":"not_applicable"}]}));
                            if std::env::var_os("COSMOLOGY_TEST_PANIC").is_some() {panic!("injected Rust panic");}
                            let slots = convert_quantities(&request.values,&request.source,&request.target)?;
                            let failed = slots.iter().any(|slot|slot.value.is_none());
                            let evaluations: Vec<_> = slots.into_iter().map(|slot|match slot.value {
                                Some(value)=>json!({"kind":"finite","value":value}),
                                None=>json!({"kind":"failure","error_id":slot.status})
                            }).collect();
                            let output=json!({"kind":if failed {"failure"} else {"finite"},
                                "error_id":if failed {Some("QUANTITY_DOMAIN_FAILURE")} else {None},
                                "source_values":request.values,"source":request.source,
                                "target":request.target,"evaluations":evaluations});
                            Ok(crate::outcome::Outcome::scientific(resolved.clone().unwrap(),output,"converted",json!("typed_physical_conversion"),json!("binary64_storage_longdouble_intermediate"),json!({"max_batch_elements":MAX_BATCH_ELEMENTS}),"bounded native/interface checks; request applicability not established"))
                        },
                        "numerics.scalar_batch.v1" => {
                            let request: NumericalRequest = serde_json::from_slice(&input).map_err(|e|e.to_string())?;
                            if request.schema_version!=2 {return Err("UNSUPPORTED_SPECIFICATION".into());}
                            resolved=Some(json!({"schema_version":2,"operation":request.operation,"method":request.method,
                                "equation_id":format!("F02/scalar/{}/v1",request.method),"values":request.values,
                                "requested_outputs":[{"id":"evaluations","required":true,"numerical_gate":"required","inference_gate":"not_applicable"}]}));
                            if std::env::var_os("COSMOLOGY_TEST_PANIC").is_some(){panic!("injected Rust panic");}
                            let slots=numerics_evaluate(&request.method,&request.values)?;
                            let failed=slots.iter().any(|slot|slot.value.is_none());
                            let evaluations:Vec<_>=slots.into_iter().map(|slot| {
                                let result=match slot.value {
                                    Some(value)=>json!({"kind":"finite","value":value}),
                                    None=>json!({"kind":"failure","error_id":slot.status})
                                };
                                json!({"result":result,"diagnostics":{"error_estimate":slot.error_estimate,
                                    "error_estimate_kind":"empirical_or_arithmetic_diagnostic_not_certified_bound",
                                    "evaluations":slot.evaluations}})
                            }).collect();
                            let output=json!({"kind":if failed {"failure"}else{"finite"},
                                "error_id":if failed {Some("NUMERICAL_EVALUATION_FAILURE")}else{None},
                                "method":request.method,"source_values":request.values,"evaluations":evaluations});
                            Ok(crate::outcome::Outcome::scientific(resolved.clone().unwrap(),output,"evaluations",json!(request.method),json!("binary64_storage_method_declared_intermediate"),json!({"max_batch_elements":MAX_BATCH_ELEMENTS}),"bounded native/interface checks; request applicability not established"))
                        },
                        "background.evaluate" => crate::current_background_run::execute(&input),
                        "supernova.profile" | "bao.density" => crate::current_session::execute_once(&input,&store),
                        "statistics.gaussian" => crate::statistics_run::execute(&input,&store),
                        "observations.prepare" => crate::current_observation_run::execute(&input,&store),
                        _ => Err("UNSUPPORTED_SPECIFICATION".into()),
                    }
                },
            )).unwrap_or_else(|_| Err("RUST_PANIC".into()));
            let success=result.as_ref().is_ok_and(|outcome|outcome.numerical_passed());
            let qualified=result.as_ref().is_ok_and(|outcome|outcome.qualification_passed());
            let assurance_satisfied=success && (assurance=="numerical_contract" || qualified);
            let output=match &result {Ok(o)=>o.output.clone(),Err(error)=>json!({"kind":"failure","error_id":error,"failure_class":if error=="RUST_PANIC" {"internal"}else{"request_or_transport"}})};
            let bytes=serde_json::to_vec(&output).map_err(|e|e.to_string())?;
            let output_hash=hash(&bytes);publish(&store.join("objects").join(&output_hash),&bytes)?;
            if std::env::var_os("COSMOLOGY_TEST_ABORT_AFTER_OUTPUT").is_some(){std::process::abort()}
            let mut final_record=initial;
            final_record["assurance"]=json!({"requested":assurance,"satisfied":assurance_satisfied,"numerical_contract":success,"qualification":if qualified {"not_applicable"}else{"unqualified"}});
            final_record["accepted"]=json!(assurance_satisfied);
            final_record["accepted_scope"]=json!(assurance);
            final_record["numerical"]=json!(if result.is_err(){"not_assessed"}else if success{"checks_passed"}else{"failed"});
            final_record["execution"]=match &result { Ok(outcome)=>serde_json::to_value(&outcome.execution).map_err(|e|e.to_string())?, Err(_)=>json!("failed") };
            if let Ok(outcome)=&result {
                resolved=Some(outcome.specification.clone());
                final_record["method"]=outcome.method.clone();
                final_record["precision"]=outcome.arithmetic.clone();
                final_record["resource_budget"]["operation"]=outcome.resources.clone();
                final_record["outputs"]=serde_json::to_value(&outcome.outputs).map_err(|e|e.to_string())?;
                final_record["interpretation"]=json!(if outcome.scientific {"unqualified"}else{"not_applicable"});
            } else {
                final_record["outputs"]=json!([]);
                final_record["failure_class"]=output["failure_class"].clone();
            }
            if let Some(spec) = resolved {
                let spec_bytes = serde_json::to_vec(&spec).unwrap();
                let spec_digest = hash(&spec_bytes);
                publish(&store.join("objects").join(&spec_digest), &spec_bytes)?;
                let execution_spec = json!({"scientific_specification_digest":spec_digest,"input_digest":input_hash,"build_id":manifest["build_id"],"executable_digest":hash(&executable),"resource_budget":final_record["resource_budget"],"rng":final_record["rng"],"runtime_libraries":final_record["runtime_libraries"],"fixture_fault_controls":final_record["fixture_fault_controls"],"fault":fault,"backend":"portable_cpu","precision":final_record["precision"],"method":final_record["method"],"assurance":assurance});
                final_record["scientific_specification_digest"] = json!(spec_digest);
                final_record["execution_identity"] =
                    json!(hash(&serde_json::to_vec(&execution_spec).unwrap()));
            }
            final_record["output_digest"] = json!(output_hash);
            publish(
                &store
                    .join("attempts")
                    .join(format!("{attempt}.complete.json")),
                &serde_json::to_vec(&final_record).unwrap(),
            )?;
            println!("{}", json!({"receipt":final_record,"result":output}));
            if !success {Err(output["error_id"].as_str().unwrap_or("EVALUATION_FAILURE").to_owned())}
            else if !assurance_satisfied {Err("NUMERICAL_QUALIFICATION_REQUIRED".into())}
            else {Ok(())}
        }
        _ => Err("usage: irred describe --json | version --json | run REQUEST STORE [--assurance numerical_contract|qualified]".into()),
    }
}
