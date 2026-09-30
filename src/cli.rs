use crate::{
    bridge::{ABI_VERSION, add},
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
pub(crate) fn execute() -> Result<(), String> {
    let args: Vec<_> = std::env::args().collect();
    let manifest: Value = serde_json::from_str(include_str!("../build/build-manifest.json"))
        .map_err(|e| e.to_string())?;
    match args.get(1).map(String::as_str) {
        Some("describe") | Some("version") if args.len() == 3 && args[2] == "--json" => {
            println!(
                "{}",
                json!({"schema_version":1,"version":env!("CARGO_PKG_VERSION"),"abi_version":ABI_VERSION,"build":manifest,"capabilities":[{"id":"fixture.checked_i64_add.v1","implementation":"implemented","qualification":"unqualified","scientific":false}],"commands":["describe --json","version --json","run REQUEST STORE"],"scientific_qualifications":[]})
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
            let initial = json!({"schema_version":1,"attempt_id":attempt,"execution":"incomplete","input_digest":input_hash,"build_id":manifest["build_id"],"executable_digest":hash(&executable),"numerical":"not_assessed","inference":"not_applicable","interpretation":"not_assessed","runtime_libraries":runtime_libraries()?,"resource_budget":{"compute_threads":1,"io_threads":1},"rng":"not_applicable","backend":"portable_cpu","precision":"exact_i64","fixture_fault_controls":{"abort":std::env::var_os("COSMOLOGY_TEST_ABORT").is_some(),"panic":std::env::var_os("COSMOLOGY_TEST_PANIC").is_some(),"abort_after_output":std::env::var_os("COSMOLOGY_TEST_ABORT_AFTER_OUTPUT").is_some()}});
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
            let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(
                || -> Result<Vec<i64>, String> {
                    let request: Request =
                        serde_json::from_slice(&input).map_err(|e| e.to_string())?;
                    if request.schema_version != 1
                        || request.operation != "fixture.checked_i64_add.v1"
                    {
                        return Err("UNSUPPORTED_SPECIFICATION".into());
                    }
                    resolved = Some(json!({
                        "schema_version": 1,
                        "operation": request.operation,
                        "equation_id": "EQ-fixture.checked_i64_add.v1",
                        "a": request.a,
                        "b": request.b,
                        "requested_outputs": [{
                            "id": "sum", "required": true,
                            "numerical_gate": "not_applicable",
                            "inference_gate": "not_applicable"
                        }]
                    }));
                    if std::env::var_os("COSMOLOGY_TEST_PANIC").is_some() {
                        panic!("injected Rust panic");
                    }
                    add(&request.a, &request.b, request.fault)
                },
            ))
            .unwrap_or_else(|_| Err("RUST_PANIC".into()));
            let output = match &result {
                Ok(v) => json!({"kind":"finite","values":v}),
                Err(e) => json!({"kind":"failure","error_id":e}),
            };
            let bytes = serde_json::to_vec(&output).unwrap();
            let output_hash = hash(&bytes);
            publish(&store.join("objects").join(&output_hash), &bytes)?;
            if std::env::var_os("COSMOLOGY_TEST_ABORT_AFTER_OUTPUT").is_some() {
                std::process::abort()
            }
            let mut final_record = initial;
            if let Some(spec) = resolved {
                let spec_bytes = serde_json::to_vec(&spec).unwrap();
                let spec_digest = hash(&spec_bytes);
                publish(&store.join("objects").join(&spec_digest), &spec_bytes)?;
                let execution_spec = json!({"scientific_specification_digest":spec_digest,"input_digest":input_hash,"build_id":manifest["build_id"],"executable_digest":hash(&executable),"resource_budget":final_record["resource_budget"],"rng":final_record["rng"],"runtime_libraries":final_record["runtime_libraries"],"fixture_fault_controls":final_record["fixture_fault_controls"],"fault":serde_json::from_slice::<Request>(&input).unwrap().fault,"backend":"portable_cpu","precision":"exact_i64"});
                final_record["scientific_specification_digest"] = json!(spec_digest);
                final_record["execution_identity"] =
                    json!(hash(&serde_json::to_vec(&execution_spec).unwrap()));
            }
            final_record["execution"] = json!(if result.is_ok() {
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
            result.map(|_| ())
        }
        _ => Err("usage: cosmology describe --json | version --json | run REQUEST STORE".into()),
    }
}
