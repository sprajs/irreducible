//! Gaussian interface/status/identity tests, core-author role; no Rust equations.
use serde_json::{Value, json};
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
    serde_json::from_slice(&fs::read("tests/fixtures/gaussian-density.json").unwrap()).unwrap()
}
fn run(spec: &Value) -> (Value, i32, Option<Value>) {
    let scratch = Scratch(std::env::temp_dir().join(format!(
            "irred-gaussian-cli-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        )));
    fs::create_dir(&scratch.0).unwrap();
    let input = scratch.0.join("request.json");
    let store = scratch.0.join("store");
    fs::write(&input, serde_json::to_vec(spec).unwrap()).unwrap();
    let out = Command::new(env!("CARGO_BIN_EXE_irred"))
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .args(["run", input.to_str().unwrap(), store.to_str().unwrap()])
        .output()
        .unwrap();
    let response: Value = serde_json::from_slice(&out.stdout).unwrap();
    let resolved = response["receipt"]["scientific_specification_digest"]
        .as_str()
        .map(|d| {
            serde_json::from_slice(&fs::read(store.join("objects").join(d)).unwrap()).unwrap()
        });
    (response, out.status.code().unwrap(), resolved)
}
fn calc(r: &Value) -> &Value {
    &r["result"]["calculation"]
}
#[test]
fn normalized_profile_and_proper_prior_modes() {
    let mut spec = request();
    spec["residuals"] = json!([[2., -3.], [0., 0.]]);
    let (r, exit, resolved) = run(&spec);
    assert_eq!(exit, 0);
    assert_eq!(r["receipt"]["accepted"], true);
    assert_eq!(r["receipt"]["execution"], "completed");
    let c = calc(&r);
    assert_eq!(c["kind"], "finite");
    assert_eq!(c["matrix_validation_scope"], "selected_covariance_only");
    assert_eq!(c["ordered_ids"], json!(["r0", "r1"]));
    assert_eq!(c["selection_history"][0]["kept_row_ids"], c["ordered_ids"]);
    assert_eq!(c["calibration_provenance"], "");
    assert_eq!(c["dependence_provenance"], "");
    assert!((c["rows"][0]["quadratic"].as_f64().unwrap() - 2.4).abs() < 4e-11);
    assert_eq!(c["rows"][1]["quadratic"], 0.);
    assert!(c["rows"][0]["log_density"].as_f64().unwrap().is_finite());
    let resolved = resolved.unwrap();
    assert_eq!(resolved["ordered_ids"], spec["ordered_ids"]);
    assert_eq!(resolved["residuals"], spec["residuals"]);
    assert_eq!(resolved["equation_id"], "F03/prepared-Gaussian/v1");
    assert_eq!(r["receipt"]["outputs"][0]["numerical"], "checks_passed");
    assert_eq!(r["receipt"]["outputs"][0]["id"], "gaussian_batch");
    assert_eq!(r["receipt"]["method"], "normalized_density");
    assert_eq!(r["receipt"]["precision"], "binary64_legacy_v1");
    assert_eq!(r["receipt"]["accepted_scope"], "numerical_contract");
    spec["mode"] = json!("profile_offset_score");
    spec["response"] = json!([1., 1.]);
    let (r, exit, _) = run(&spec);
    assert_eq!(exit, 0);
    let row = &calc(&r)["rows"][0];
    assert_eq!(row["density"], Value::Null);
    assert_eq!(row["normalization"], "not_applicable");
    assert!(row.get("log_density").is_none());
    assert!((row["profile_coefficient"].as_f64().unwrap() - 7. / 11.).abs() < 1e-10);
    assert!((row["quadratic"].as_f64().unwrap() - 25. / 11.).abs() < 1e-10);
    assert_eq!(row["numerical_status"], "ok");
    assert!(row["estimated_forward_sensitivity"].as_f64().unwrap() > 0.);
    let mut tighter = spec.clone();
    tighter["maximum_forward_sensitivity"] = json!(1e-30);
    let (failed, code, _) = run(&tighter);
    assert_ne!(code, 0);
    assert_eq!(
        calc(&failed)["rows"][0]["numerical_status"],
        "conditioning_budget_exceeded"
    );
    assert!(calc(&failed)["rows"][0].get("quadratic").is_none());
    spec["mode"] = json!("normalized_density");
    spec["response"] = Value::Null;
    spec["proper_prior"] = json!({"response":[1.,1.],"ordered_ids":["r0","r1"],"mean":1./3.,"variance":2.,"latent_identity":"proper-offset","independence_declared":true});
    let (r, exit, _) = run(&spec);
    assert_eq!(exit, 0);
    let c = calc(&r);
    assert_eq!(c["priors"][0]["mean"], spec["proper_prior"]["mean"]);
    assert_eq!(c["priors"][0]["variance"], 2.);
    assert_eq!(c["priors"][0]["applied_row_ids"], spec["ordered_ids"]);
    assert_eq!(c["priors"][0]["independence_declared"], true);
    assert_eq!(c["dependence_provenance"], "");
    assert!((c["rows"][0]["quadratic"].as_f64().unwrap() - 1175. / 513.).abs() < 4e-11);
}
#[test]
fn row_failures_no_sentinel_and_order_checks() {
    let mut spec = request();
    spec["residuals"] = json!([[1., 2.], [f64::MAX, -f64::MAX]]);
    let (r, exit, _) = run(&spec);
    assert_ne!(exit, 0);
    assert_ne!(exit, 6);
    assert_eq!(r["receipt"]["execution"], "completed");
    assert_eq!(calc(&r)["kind"], "failure");
    assert_eq!(calc(&r)["rows"][0]["kind"], "finite");
    let row = &calc(&r)["rows"][1];
    assert_eq!(row["kind"], "failure");
    for field in [
        "log_density",
        "quadratic",
        "log_determinant",
        "normalization",
        "coefficient",
    ] {
        assert!(
            row.get(field).is_none(),
            "failure payload must not expose {field}"
        );
    }
    spec = request();
    spec["ordered_ids"] = json!(["r1", "r0"]);
    let (r, exit, _) = run(&spec);
    assert_ne!(exit, 0);
    assert_eq!(calc(&r)["rows"][0]["kind"], "failure");
    assert_eq!(calc(&r)["rows"][0]["status"], "incompatible_metadata");
    spec = request();
    spec["residuals"] = json!([]);
    let (r, exit, _) = run(&spec);
    assert_eq!(exit, 0);
    assert_eq!(calc(&r)["rows"], json!([]));
}
#[test]
fn malformed_units_shapes_modes_and_priors() {
    let baseline = request();
    for (field, value) in [
        ("schema_version", json!(2)),
        ("mode", json!("invented")),
        ("residual_unit", json!("metre")),
        ("response_unit", json!("magnitude")),
        ("response", json!([1., 1.])),
        ("residuals", json!([[1.]])),
        ("maximum_batch_elements", json!(1)),
        ("maximum_forward_sensitivity", json!(0)),
        ("ignored", json!(true)),
    ] {
        let mut s = baseline.clone();
        s[field] = value;
        let (r, exit, _) = run(&s);
        assert_ne!(exit, 0, "{field}");
        assert_eq!(r["result"]["kind"], "failure", "{field}");
        assert_eq!(r["receipt"]["accepted"], false);
    }
    let mut s = baseline.clone();
    s["mode"] = json!("profile_offset_score");
    s["response"] = json!([1., 1.]);
    s["proper_prior"] = json!({"response":[1.,1.],"ordered_ids":["r0","r1"],"mean":0.,"variance":1.,"latent_identity":"x","independence_declared":true});
    let (r, exit, _) = run(&s);
    assert_ne!(exit, 0);
    assert_eq!(r["result"]["kind"], "failure");
    s = baseline;
    s["proper_prior"] = json!({"response":[1.,1.],"ordered_ids":["r0","r1"],"mean":0.,"variance":1.,"latent_identity":"x","independence_declared":false});
    let (r, exit, _) = run(&s);
    assert_ne!(exit, 0);
    assert_eq!(r["result"]["kind"], "failure");
}
