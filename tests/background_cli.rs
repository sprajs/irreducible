//! Bridge-only background transport/status checks. Mathematical qualification
//! comes from independent native analytic/fixed-panel P01 fixtures.
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
fn request() -> Value {
    json!({
    "schema_version":1,"operation":"background.parameter_query_batch.v1",
    "parameters":[{"model":"flat_lcdm_late_v1","h0_km_s_mpc":70.,"omega_m":0.3,"constant_q":0.},{"model":"constant_q_flat_v1","h0_km_s_mpc":70.,"omega_m":0.,"constant_q":-1.}],
    "queries":[{"z_expansion":0.,"z_observer":0.,"convention":"geometric_same_redshift"},{"z_expansion":1.,"z_observer":1.,"convention":"geometric_same_redshift"},{"z_expansion":1.,"z_observer":0.9,"convention":"released_zhd_zhel"}],
    "policy":{"maximum_parameters":2,"maximum_queries":3,"maximum_slots":6,"maximum_native_output_bytes":1048576,"maximum_total_evaluations":2000000,"maximum_evaluations_per_integral":100000,"maximum_depth":30,"absolute_tolerance":1e-13,"relative_tolerance":1e-12}
    })
}
fn run_raw(bytes: &[u8]) -> (Value, i32, Option<Value>) {
    let scratch = Scratch(std::env::temp_dir().join(format!(
            "irred-background-cli-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        )));
    fs::create_dir(&scratch.0).unwrap();
    let input = scratch.0.join("request.json");
    let store = scratch.0.join("store");
    fs::write(&input, bytes).unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_irred"))
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .args(["run", input.to_str().unwrap(), store.to_str().unwrap()])
        .output()
        .unwrap();
    let response: Value = serde_json::from_slice(&output.stdout).unwrap();
    let resolved = response["receipt"]["scientific_specification_digest"]
        .as_str()
        .map(|d| {
            serde_json::from_slice(&fs::read(store.join("objects").join(d)).unwrap()).unwrap()
        });
    (response, output.status.code().unwrap(), resolved)
}
fn run(spec: &Value) -> (Value, i32, Option<Value>) {
    run_raw(&serde_json::to_vec(spec).unwrap())
}
fn calculation(v: &Value) -> &Value {
    &v["result"]
}
#[test]
fn batch_order_units_identity_and_zero() {
    let spec = request();
    let (response, exit, resolved) = run(&spec);
    assert_eq!(exit, 6);
    assert_eq!(response["receipt"]["accepted"], false);
    assert_eq!(response["receipt"]["execution"], "completed");
    let c = calculation(&response);
    assert_eq!(c["kind"], "finite");
    assert_eq!(c["rows"].as_array().unwrap().len(), 6);
    for (index, row) in c["rows"].as_array().unwrap().iter().enumerate() {
        assert_eq!(row["identity"]["parameter_index"], index / 3);
        assert_eq!(row["identity"]["query_index"], index % 3);
        assert_eq!(row["result"]["kind"], "finite");
        assert_eq!(row["diagnostics"]["numerical_status"], "ok");
        assert!(row["identity"]["constants_id"].as_str().unwrap().len() > 0);
    }
    let zero = &c["rows"][0]["result"];
    assert_eq!(zero["radial_mpc"], 0.);
    assert_eq!(zero["luminosity_mpc"], 0.);
    assert_eq!(zero["jerk"], 1.);
    assert!((zero["deceleration_q"].as_f64().unwrap() + 0.55).abs() < 1e-14);
    let de_sitter = &c["rows"][4]["result"];
    assert_eq!(de_sitter["expansion_E"], 1.);
    assert_eq!(de_sitter["h_km_s_mpc"], 70.);
    assert_eq!(de_sitter["deceleration_q"], -1.);
    assert_eq!(de_sitter["jerk"], 1.);
    assert!((de_sitter["luminosity_mpc"].as_f64().unwrap() - 2. * 299792.458 / 70.).abs() < 1e-7);
    assert!(
        (c["rows"][5]["result"]["dimensionless_luminosity_shape"]
            .as_f64()
            .unwrap()
            - 1.9)
            .abs()
            < 1e-12
    );
    assert_ne!(
        c["rows"][4]["identity"]["luminosity_equation_id"],
        c["rows"][5]["identity"]["luminosity_equation_id"]
    );
    let resolved = resolved.unwrap();
    assert_eq!(resolved["parameters"], spec["parameters"]);
    assert_eq!(resolved["queries"], spec["queries"]);
    assert_eq!(resolved["policy"], spec["policy"]);
    assert_eq!(resolved["units"]["distance"], "Mpc");
    assert_eq!(resolved["units"]["lookback"], "s");
    assert_eq!(resolved["units"]["volume"], "Mpc^3/sr/redshift");
    assert!(c["evaluations"].as_u64().unwrap() <= 2000000);
    assert_eq!(
        response["receipt"]["outputs"][0]["numerical"],
        "not_assessed"
    );
}
#[test]
fn scientific_failures_omit_payload_and_preserve_work() {
    let mut spec = request();
    spec["policy"]["maximum_total_evaluations"] = json!(1);
    let (response, exit, _) = run(&spec);
    assert_ne!(exit, 0);
    assert_ne!(exit, 6);
    let c = calculation(&response);
    assert_eq!(c["kind"], "failure");
    assert!(c["evaluations"].as_u64().unwrap() <= 1);
    for (i, row) in c["rows"].as_array().unwrap().iter().enumerate() {
        if i % 3 == 0 {
            assert_eq!(row["result"]["kind"], "finite");
        } else {
            assert_eq!(row["result"]["status"], "work_limit");
            assert_eq!(row["diagnostics"]["numerical_status"], "work_limit");
            assert!(row["result"].get("luminosity_mpc").is_none());
        }
    }
    let mut spec = request();
    spec["queries"][1]["z_observer"] = json!(0.9);
    let (response, exit, _) = run(&spec);
    assert_ne!(exit, 6);
    assert_eq!(
        calculation(&response)["rows"][1]["result"]["status"],
        "incompatible_convention"
    );
    assert!(calculation(&response)["rows"][1]["result"]
        .get("expansion_E")
        .is_none());
    let mut spec = request();
    spec["queries"][1]["z_expansion"] = json!(6.);
    spec["queries"][1]["z_observer"] = json!(6.);
    let (response, _, _) = run(&spec);
    assert_eq!(
        calculation(&response)["rows"][1]["result"]["status"],
        "unsupported_domain"
    );
    let mut spec = request();
    spec["parameters"][0]["h0_km_s_mpc"] = json!(f64::MAX);
    let (response, _, _) = run(&spec);
    assert_eq!(
        calculation(&response)["rows"][1]["result"]["status"],
        "numerical_failure"
    );
}
#[test]
fn strict_metadata_and_resource_caps() {
    for (field, value) in [
        ("maximum_slots", json!(5)),
        ("maximum_native_output_bytes", json!(1)),
        ("maximum_parameters", json!(1)),
        ("maximum_queries", json!(2)),
    ] {
        let mut spec = request();
        spec["policy"][field] = value;
        let (response, exit, _) = run(&spec);
        assert_ne!(exit, 0);
        assert_ne!(exit, 6);
        assert_eq!(calculation(&response)["error_id"], "RESOURCE_LIMIT");
    }
    let mut spec = request();
    spec["parameters"][0]["model"] = json!("unknown_model");
    let (response, exit, _) = run(&spec);
    assert_ne!(exit, 6);
    assert_eq!(
        calculation(&response)["error_id"],
        "UNKNOWN_BACKGROUND_MODEL"
    );
    let mut spec = request();
    spec["queries"][0]["convention"] = json!("implicit_prefactor");
    let (response, exit, _) = run(&spec);
    assert_ne!(exit, 6);
    assert_eq!(
        calculation(&response)["error_id"],
        "UNKNOWN_BACKGROUND_CONVENTION"
    );
    let mut spec = request();
    spec["parameters"][0]["secret_default"] = json!(1);
    let (_, exit, _) = run(&spec);
    assert_ne!(exit, 0);
    assert_ne!(exit, 6);
    let text = serde_json::to_string(&request()).unwrap();
    let duplicate = text.replacen(
        "\"schema_version\":1",
        "\"schema_version\":1,\"schema_version\":1",
        1,
    );
    let (_, exit, _) = run_raw(duplicate.as_bytes());
    assert_ne!(exit, 0);
    assert_ne!(exit, 6);
}
