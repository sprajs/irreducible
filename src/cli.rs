use crate::{
    bridge::{
        ABI_VERSION, MAX_BATCH_ELEMENTS, Metadata, add, convert_quantities, numerics_evaluate,
    },
    records::{hash, publish, runtime_libraries},
};
use serde::Deserialize;
use serde_json::{Value, json};
use std::{fs, path::PathBuf};
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
        Some("describe") | Some("version") if args.len() == 3 && args[2] == "--json" => {
            println!(
                "{}",
                json!({"schema_version":1,"product":"Irreducible","executable":"irred","version":env!("CARGO_PKG_VERSION"),"abi_version":ABI_VERSION,"build":manifest,"capabilities":[{"id":"fixture.checked_i64_add.v1","implementation":"implemented","qualification":"unqualified","scientific":false},{"id":"quantity.convert.v1","implementation":"implemented","qualification":"unqualified","scientific":true},{"id":"numerics.scalar_batch.v1","implementation":"implemented","qualification":"unqualified","scientific":true},{"id":"background.parameter_query_batch.v1","implementation":"implemented","qualification":"unqualified","scientific":true,"interface_gate":"bounded_native_and_cli_passed","models":["flat_lcdm_late_v1","constant_q_flat_v1"]},{"id":"supernova.profile_batch.v1","implementation":"implemented","qualification":"unqualified","scientific":true,"interface_gate":"bounded_native_and_cli_passed","models":["flat_lcdm_late_v1","constant_q_flat_v1"],"target":"relative single-offset profile score; density and evidence not applicable","validation_coverage":"named original-input seven-point native regression; request applicability not established"},{"id":"background.parameter_query_batch.v2","implementation":"implemented","qualification":"unqualified","scientific":true,"interface_gate":"bounded_native_and_cli_passed","models":["flat_lcdm_late_v1","constant_q_flat_v1","flat_cpl_late_v1"],"parameters":"all fields explicit; canonical inactive values required","validation_coverage":"bounded native model and synthetic interface fixtures; request applicability not established"},{"id":"supernova.profile_batch.v2","implementation":"implemented","qualification":"unqualified","scientific":true,"interface_gate":"bounded_native_and_cli_passed","models":["flat_lcdm_late_v1","constant_q_flat_v1","flat_cpl_late_v1"],"target":"relative single-offset profile score; density and evidence not applicable","validation_coverage":"named four-point original-input native comparison and synthetic interface cases; request applicability not established"},{"id":"background.piecewise_query_batch.v1","implementation":"implemented","qualification":"unqualified","scientific":true,"interface_gate":"bounded_native_and_cli_passed","parameters":"required H0 and five ordered q coefficients","work":"one global analytic segment budget; no quadrature callbacks","validation_coverage":"bounded native fixed-five-bin and hostile transport cases; request applicability not established"},{"id":"bao.gaussian_batch.v1","implementation":"implemented","qualification":"unqualified","scientific":true,"interface_gate":"bounded_native_and_cli_passed","target":"normalized conditional free-ruler Gaussian density","computational_h0_km_s_mpc":70,"validation_coverage":"named native eleven-point original-input regression; request applicability not established"},{"id":"statistics.gaussian_batch.v1","implementation":"implemented","qualification":"unqualified","scientific":true,"interface_gate":"bounded_native_and_cli_passed","modes":["normalized_density","profile_offset_score"]},{"id":"observations.prepare.v1","implementation":"implemented","qualification":"unqualified","scientific":true,"profiles":["pantheon_plus_released_v1","gaussian_fixture_v1"],"fits_codec":"unavailable_in_product"}],"abi_schema":serde_json::from_str::<Value>(include_str!("../schema/abi.json")).map_err(|e|e.to_string())?,"commands":["describe --json","version --json","run REQUEST STORE"],"scientific_qualifications":[]})
            );
            Ok(())
        }
        Some("run") if args.len() == 4 => {
            let input =
                fs::read(args.get(2).ok_or("missing request path")?).map_err(|e| e.to_string())?;
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
                || -> Result<Value, String> {
                    let header: OperationHeader = serde_json::from_slice(&input).map_err(|e|e.to_string())?;
                    match header.operation.as_str() {
                        "fixture.checked_i64_add.v1" => {
                            let request: Request = serde_json::from_slice(&input).map_err(|e|e.to_string())?;
                            if request.schema_version != 1 { return Err("UNSUPPORTED_SPECIFICATION".into()); }
                            fault = request.fault;
                            resolved = Some(json!({"schema_version":1,"operation":request.operation,
                                "equation_id":"EQ-fixture.checked_i64_add.v1","a":request.a,"b":request.b,
                                "requested_outputs":[{"id":"sum","required":true,"numerical_gate":"not_applicable","inference_gate":"not_applicable"}]}));
                            if std::env::var_os("COSMOLOGY_TEST_PANIC").is_some() {panic!("injected Rust panic");}
                            add(&request.a,&request.b,request.fault).map(|v|json!({"kind":"finite","values":v}))
                        },
                        "quantity.convert.v1" => {
                            let request: QuantityRequest = serde_json::from_slice(&input).map_err(|e|e.to_string())?;
                            if request.schema_version != 1 { return Err("UNSUPPORTED_SPECIFICATION".into()); }
                            resolved = Some(json!({"schema_version":1,"operation":request.operation,
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
                            Ok(json!({"kind":if failed {"failure"} else {"finite"},
                                "error_id":if failed {Some("QUANTITY_DOMAIN_FAILURE")} else {None},
                                "source_values":request.values,"source":request.source,
                                "target":request.target,"evaluations":evaluations}))
                        },
                        "numerics.scalar_batch.v1" => {
                            let request: NumericalRequest = serde_json::from_slice(&input).map_err(|e|e.to_string())?;
                            if request.schema_version!=1 {return Err("UNSUPPORTED_SPECIFICATION".into());}
                            resolved=Some(json!({"schema_version":1,"operation":request.operation,"method":request.method,
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
                            Ok(json!({"kind":if failed {"failure"}else{"finite"},
                                "error_id":if failed {Some("NUMERICAL_EVALUATION_FAILURE")}else{None},
                                "method":request.method,"source_values":request.values,"evaluations":evaluations}))
                        },
                        "background.piecewise_query_batch.v1" => {
                            let (spec,output)=crate::piecewise_run::execute(&input,&store)?;resolved=Some(spec);Ok(output)
                        }
                        "bao.gaussian_batch.v1" => {
                            let (spec,output)=crate::bao_run::execute(&input,&store)?;resolved=Some(spec);Ok(output)
                        }
                        "background.parameter_query_batch.v2" => {
                            let (spec,output)=crate::background_v2_run::execute(&input,&store)?;resolved=Some(spec);Ok(output)
                        }
                        "supernova.profile_batch.v2" => {
                            let (spec,output)=crate::supernova_v2_run::execute(&input,&store)?;resolved=Some(spec);Ok(output)
                        }
                        "background.parameter_query_batch.v1" => {
                            let (spec,output)=crate::background_run::execute(&input,&store)?;resolved=Some(spec);Ok(output)
                        }
                        "supernova.profile_batch.v1" => {
                            let (spec,output)=crate::supernova_run::execute(&input,&store)?;resolved=Some(spec);Ok(output)
                        }
                        "statistics.gaussian_batch.v1" => {
                            let (spec,output)=crate::statistics_run::execute(&input,&store)?;resolved=Some(spec);Ok(output)
                        }
                        "observations.prepare.v1" => {
                            let (spec,output)=crate::observation_run::execute(&input,&store)?;resolved=Some(spec);Ok(output)
                        },
                        _ => Err("UNSUPPORTED_SPECIFICATION".into()),
                    }
                },
            )).unwrap_or_else(|_| Err("RUST_PANIC".into()));
            let output = match &result {
                Ok(value) => value.clone(),
                Err(error) => json!({"kind":"failure","error_id":error}),
            };
            let success = result.is_ok() && output["kind"] != "failure";
            let bytes = serde_json::to_vec(&output).unwrap();
            let output_hash = hash(&bytes);
            publish(&store.join("objects").join(&output_hash), &bytes)?;
            if std::env::var_os("COSMOLOGY_TEST_ABORT_AFTER_OUTPUT").is_some() {
                std::process::abort()
            }
            let is_numerical = resolved
                .as_ref()
                .is_some_and(|spec| spec["operation"] == "numerics.scalar_batch.v1");
            let is_quantity = resolved
                .as_ref()
                .is_some_and(|spec| spec["operation"] == "quantity.convert.v1");
            let is_observation = resolved.as_ref().is_some_and(|spec| {
                spec["operation"] == "background.piecewise_query_batch.v1"
                    || spec["operation"] == "bao.gaussian_batch.v1"
                    || spec["operation"] == "observations.prepare.v1"
                    || spec["operation"] == "statistics.gaussian_batch.v1"
                    || spec["operation"] == "background.parameter_query_batch.v1"
                    || spec["operation"] == "background.parameter_query_batch.v2"
                    || spec["operation"] == "supernova.profile_batch.v1"
                    || spec["operation"] == "supernova.profile_batch.v2"
            });
            let is_supernova = resolved.as_ref().is_some_and(|spec| {
                spec["operation"] == "supernova.profile_batch.v1"
                    || spec["operation"] == "supernova.profile_batch.v2"
            });
            let is_piecewise = resolved.as_ref().is_some_and(|spec|spec["operation"] == "background.piecewise_query_batch.v1");
            let is_bao = resolved
                .as_ref()
                .is_some_and(|spec| spec["operation"] == "bao.gaussian_batch.v1");
            let mut final_record = initial;
            final_record["precision"] = json!(if is_observation {
                "binary64_source_parse_no_scientific_transformation"
            } else if is_numerical {
                "binary64_storage_method_declared_intermediate"
            } else if is_quantity {
                "binary64_storage_host_long_double_intermediate"
            } else {
                "exact_i64"
            });
            if is_supernova {
                final_record["precision"] =
                    resolved.as_ref().unwrap()["policy"]["arithmetic"].clone();
                final_record["resource_budget"]["supernova_policy"] =
                    resolved.as_ref().unwrap()["policy"].clone();
            }
            if is_piecewise {
                final_record["precision"]=json!("binary64_storage_native_longdouble_analytic_intermediate");
                final_record["resource_budget"]["piecewise_policy"]=resolved.as_ref().unwrap()["policy"].clone();
            }
            if is_bao {
                final_record["precision"] =
                    resolved.as_ref().unwrap()["evaluation_policy"]["arithmetic"].clone();
                final_record["resource_budget"]["bao_prepare_policy"] =
                    resolved.as_ref().unwrap()["prepare_policy"].clone();
                final_record["resource_budget"]["bao_evaluation_policy"] =
                    resolved.as_ref().unwrap()["evaluation_policy"].clone();
            }
            final_record["accepted"] =
                json!(success && !is_quantity && !is_numerical && !is_observation);
            if is_numerical {
                final_record["resource_budget"]["max_batch_elements"] = json!(MAX_BATCH_ELEMENTS);
            }
            if is_quantity || is_numerical || is_observation {
                final_record["outputs"] = json!([{"id":if is_observation {"prepared_observations"} else if is_numerical {"evaluations"} else {"converted"},"required":true,"numerical":"not_assessed","inference":"not_applicable","evidence":[]}]);
            }
            if is_supernova {
                let validation_coverage = if resolved
                    .as_ref()
                    .is_some_and(|spec| spec["operation"] == "supernova.profile_batch.v2")
                {
                    "named original-input four-point CPL regression; request applicability not established"
                } else {
                    "named original-input seven-point regression; request applicability not established"
                };
                final_record["outputs"] = json!([{ "id":"profile_scores", "required":true, "numerical":output["calculation"]["numerical_status"].as_str().unwrap_or("not_evaluated"), "inference":"not_applicable", "interpretation":"unqualified", "evidence":[], "validation_coverage":validation_coverage }]);
            }
            if is_piecewise {
                final_record["outputs"]=json!([{ "id":"piecewise_background_slots","required":true,"numerical":if success {"checks_passed"}else{"failed"},"inference":"not_applicable","interpretation":"unqualified","evidence":[],"validation_coverage":"named native fixed-five-bin comparisons; request applicability not established" }]);
            }
            if is_bao {
                final_record["outputs"] = json!([{ "id":"bao_densities", "required":true, "numerical":output["calculation"]["numerical_status"].as_str().unwrap_or("not_evaluated"), "inference":"not_applicable", "interpretation":"unqualified", "evidence":[], "validation_coverage":"named native eleven-point original-input regression; request applicability not established" }]);
            }
            if let Some(spec) = resolved {
                let spec_bytes = serde_json::to_vec(&spec).unwrap();
                let spec_digest = hash(&spec_bytes);
                publish(&store.join("objects").join(&spec_digest), &spec_bytes)?;
                let execution_spec = json!({"scientific_specification_digest":spec_digest,"input_digest":input_hash,"build_id":manifest["build_id"],"executable_digest":hash(&executable),"resource_budget":final_record["resource_budget"],"rng":final_record["rng"],"runtime_libraries":final_record["runtime_libraries"],"fixture_fault_controls":final_record["fixture_fault_controls"],"fault":fault,"backend":"portable_cpu","precision":final_record["precision"]});
                final_record["scientific_specification_digest"] = json!(spec_digest);
                final_record["execution_identity"] =
                    json!(hash(&serde_json::to_vec(&execution_spec).unwrap()));
            }
            final_record["execution"] =
                json!(if success || ((is_supernova || is_bao || is_piecewise) && result.is_ok()) {
                    "completed"
                } else {
                    "failed"
                });
            final_record["output_digest"] = json!(output_hash);
            publish(
                &store
                    .join("attempts")
                    .join(format!("{attempt}.complete.json")),
                &serde_json::to_vec(&final_record).unwrap(),
            )?;
            println!("{}", json!({"receipt":final_record,"result":output}));
            if success && (is_quantity || is_numerical || is_observation) {
                Err("NUMERICAL_QUALIFICATION_REQUIRED".into())
            } else if success {
                Ok(())
            } else {
                Err(output["error_id"]
                    .as_str()
                    .unwrap_or("EVALUATION_FAILURE")
                    .to_string())
            }
        }
        _ => Err("usage: irred describe --json | version --json | run REQUEST STORE".into()),
    }
}
