//! Original inline fixture and operation-owned records; no Rust posterior equations.
use serde_json::{json, Value};
use std::{
    fs,
    path::PathBuf,
    process::Command,
    time::{SystemTime, UNIX_EPOCH},
};
struct Scratch(PathBuf);
impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}
fn fixture() -> Value {
    serde_json::from_str(include_str!("fixtures/gaussian-posterior.json")).unwrap()
}
fn run_raw(bytes: &[u8]) -> (Value, i32, Option<Value>) {
    let d = Scratch(std::env::temp_dir().join(format!(
            "irred-proper-posterior-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        )));
    fs::create_dir(&d.0).unwrap();
    let p = d.0.join("request.json");
    let store = d.0.join("store");
    fs::write(&p, bytes).unwrap();
    let out = Command::new(env!("CARGO_BIN_EXE_irred"))
        .args(["run", p.to_str().unwrap(), store.to_str().unwrap()])
        .output()
        .unwrap();
    let v: Value = serde_json::from_slice(&out.stdout).unwrap();
    let digest = v["receipt"]["input_digest"].as_str().unwrap();
    assert_eq!(fs::read(store.join("objects").join(digest)).unwrap(), bytes);
    let spec = v["receipt"]["scientific_specification_digest"]
        .as_str()
        .map(|x| {
            serde_json::from_slice(&fs::read(store.join("objects").join(x)).unwrap()).unwrap()
        });
    (v, out.status.code().unwrap(), spec)
}
fn run(v: &Value) -> (Value, i32, Option<Value>) {
    run_raw(&serde_json::to_vec(v).unwrap())
}
fn near(v: &Value, x: f64) {
    assert!((v.as_f64().unwrap() - x).abs() <= 2e-12 * (1. + x.abs()));
}
#[test]
fn original_two_parameter_family_and_record_parity() {
    let input = fixture();
    let (v, code, spec) = run(&input);
    assert_eq!(code, 0);
    assert_eq!(v["receipt"]["execution"], "completed");
    assert_eq!(v["receipt"]["accepted_scope"], "numerical_contract");
    let o = &v["result"];
    near(&o["rows"][0]["mean"][0], 20604. / 25511.);
    near(&o["rows"][0]["mean"][1], -7125. / 51022.);
    near(&o["rows"][1]["mean"][0], 0.5);
    near(&o["rows"][1]["mean"][1], -0.25);
    for (i, x) in [15160., 4466., 4466., 8592.].iter().enumerate() {
        near(&o["covariance"][i], x / 25511.);
    }
    assert_eq!(o["event_ids"], input["noise"]["event_ids"]);
    assert_eq!(
        o["ordered_parameter_ids"],
        input["parameter_prior"]["ordered_parameter_ids"]
    );
    assert_eq!(v["receipt"]["outputs"][0]["id"], "posterior_covariance");
    assert_eq!(v["receipt"]["outputs"][1]["id"], "posterior_means");
    assert_eq!(
        v["receipt"]["method"]["actual"],
        "proper-Gaussian-parameter-posterior/whitened-precision/v1"
    );
    assert_eq!(
        v["receipt"]["precision"]["actual_source"],
        "F02/longdouble-cpu/v1"
    );
    let spec = spec.unwrap();
    assert_eq!(
        spec["conditioning"]["vectors"],
        input["conditioning"]["vectors"]
    );
    assert_eq!(spec["parameter_prior"], input["parameter_prior"]);
    assert_eq!(spec["noise"], input["noise"]);
    assert!(o["covariance_relative_error_estimate"].as_f64().unwrap() >= 0.);
    assert!(o.get("log_density").is_none());
    assert!(o.get("evidence").is_none());
}
#[test]
fn case_refusal_keeps_covariance_and_original_neighbors() {
    let mut input = fixture();
    input["conditioning"]["vectors"][1] = json!([1e-320, 0.]);
    let (v, code, spec) = run(&input);
    assert_ne!(code, 0);
    assert_eq!(v["receipt"]["execution"], "completed");
    assert_eq!(v["receipt"]["outputs"][0]["numerical"], "checks_passed");
    assert_eq!(v["receipt"]["outputs"][1]["numerical"], "failed");
    assert_eq!(v["result"]["rows"][0]["kind"], "finite");
    assert_eq!(v["result"]["rows"][1]["kind"], "failure");
    assert!(v["result"]["rows"][1].get("mean").is_none());
    assert!(v["result"]["rows"][1]
        .get("absolute_error_estimates")
        .is_none());
    assert_eq!(spec.unwrap()["conditioning"], input["conditioning"]);
}
#[test]
fn quota_zero_is_completed_without_claiming_native_execution() {
    let mut input = fixture();
    input["resource_policy"]["maximum_native_bytes"] = json!(0);
    let (v, code, spec) = run(&input);
    assert_ne!(code, 0);
    assert_eq!(v["receipt"]["execution"], "completed");
    assert_eq!(v["result"]["native_payload_present"], false);
    assert_eq!(v["result"]["numerical_status"], "work_limit");
    assert_eq!(v["receipt"]["method"]["actual"], Value::Null);
    assert_eq!(v["receipt"]["method"]["executed"], false);
    assert_eq!(v["result"]["covariance"], Value::Null);
    assert_eq!(v["result"]["rows"], json!([]));
    assert_eq!(spec.unwrap()["resource_policy"]["maximum_native_bytes"], 0);
}
#[test]
fn unchanged_quota_and_semantic_refusals() {
    for (path, value) in [
        ("maximum_work_units", json!(3583)),
        ("maximum_elements", json!(29)),
        ("maximum_cases", json!(1)),
        ("maximum_string_bytes", json!(0)),
    ] {
        let mut input = fixture();
        input["resource_policy"][path] = value;
        let (v, code, _) = run(&input);
        assert_ne!(code, 0);
        assert_eq!(v["receipt"]["execution"], "completed");
        assert_eq!(v["result"]["numerical_status"], "work_limit");
        assert_eq!(v["result"]["covariance"], Value::Null);
    }
    let mut input = fixture();
    input["parameter_prior"]["noise_independence_declared"] = json!(false);
    let (v, code, _) = run(&input);
    assert_ne!(code, 0);
    assert_eq!(v["receipt"]["execution"], "completed");
    assert_eq!(v["result"]["status"], "incompatible_metadata");
    input = fixture();
    input["conditioning"]["event_ids"] = json!(["synthetic-event-1", "synthetic-event-0"]);
    let (v, code, _) = run(&input);
    assert_ne!(code, 0);
    assert_eq!(v["result"]["status"], "incompatible_metadata");
    input = fixture();
    input["parameter_prior"]["covariance_row_major"] = json!([1., 2., 2., 1.]);
    let (v, code, _) = run(&input);
    assert_ne!(code, 0);
    assert_eq!(v["receipt"]["execution"], "completed");
    assert_eq!(v["result"]["numerical_status"], "not_positive_definite");
    assert_eq!(v["receipt"]["method"]["executed"], true);
    assert_eq!(v["receipt"]["method"]["noise_prepared"], true);
    assert_eq!(
        v["receipt"]["method"]["posterior_preparation_attempted"],
        true
    );
    assert_eq!(v["receipt"]["method"]["posterior_prepared"], false);
    assert_eq!(v["receipt"]["method"]["conditioning_cases_attempted"], 0);
    assert_eq!(
        v["receipt"]["precision"]["actual_source"],
        "F02/longdouble-cpu/v1"
    );
    assert_eq!(v["receipt"]["precision"]["posterior_products"], Value::Null);
    input = fixture();
    input["noise"]["covariance_row_major"] = json!([1., 2., 2., 1.]);
    let (v, code, _) = run(&input);
    assert_ne!(code, 0);
    assert_eq!(v["result"]["phase"], "noise_preparation");
    assert_eq!(v["receipt"]["method"]["executed"], true);
    assert_eq!(
        v["receipt"]["method"]["actual"],
        "supplied-Gaussian-covariance-preparation/Cholesky/v1"
    );
    assert_eq!(v["receipt"]["method"]["noise_prepared"], false);
    assert_eq!(v["receipt"]["precision"]["actual_source"], Value::Null);
    assert_eq!(
        v["receipt"]["precision"]["attempted_source"],
        "F02/longdouble-cpu/v1"
    );
}
#[test]
fn proper_prior_retains_unmeasured_parameter_for_rank_deficient_and_underdetermined_designs() {
    let mut input = fixture();
    input["parameter_prior"]["covariance_row_major"] = json!([1., 0., 0., 0.5]);
    input["design"]["values_row_major"] = json!([1., 0., 2., 0.]);
    let (v, code, _) = run(&input);
    assert_eq!(code, 0);
    near(&v["result"]["covariance"][3], 0.5);
    near(&v["result"]["covariance"][1], 0.);
    near(&v["result"]["rows"][0]["mean"][1], -0.25);
    input["noise"]["ordered_row_ids"] = json!(["r0"]);
    input["noise"]["event_ids"] = json!(["synthetic-event-0"]);
    input["noise"]["covariance_row_major"] = json!([2.]);
    input["design"]["ordered_row_ids"] = json!(["r0"]);
    input["design"]["values_row_major"] = json!([1., 0.]);
    input["conditioning"]["ordered_row_ids"] = json!(["r0"]);
    input["conditioning"]["event_ids"] = json!(["synthetic-event-0"]);
    input["conditioning"]["vectors"] = json!([[1.25], [0.5]]);
    let (v, code, _) = run(&input);
    assert_eq!(code, 0);
    near(&v["result"]["covariance"][0], 2. / 3.);
    near(&v["result"]["covariance"][3], 0.5);
    near(&v["result"]["rows"][0]["mean"][0], 0.75);
    near(&v["result"]["rows"][0]["mean"][1], -0.25);
    near(&v["result"]["rows"][1]["mean"][0], 0.5);
}
#[test]
fn strict_nested_fields_duplicates_and_discovery() {
    let mut input = fixture();
    input["noise"]["ignored"] = json!(true);
    let (v, code, _) = run(&input);
    assert_ne!(code, 0);
    assert_eq!(v["receipt"]["execution"], "failed");
    let raw = serde_json::to_string(&fixture()).unwrap().replace(
        "\"noise_independence_declared\":true",
        "\"noise_independence_declared\":false,\"noise_independence_declared\":true",
    );
    let (v, code, spec) = run_raw(raw.as_bytes());
    assert_ne!(code, 0);
    assert!(spec.is_none());
    assert_eq!(v["receipt"]["execution"], "failed");
    let out = Command::new(env!("CARGO_BIN_EXE_irred"))
        .args(["describe", "--json"])
        .output()
        .unwrap();
    let v: Value = serde_json::from_slice(&out.stdout).unwrap();
    assert!(v["capabilities"]
        .as_array()
        .unwrap()
        .iter()
        .any(|x| x["id"] == "statistics.gaussian_posterior"));
}
