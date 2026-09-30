//! Independent v2 boundary tests, synthetic input only. Native physics peer
//! fixtures remain separate; successful execution never autoqualifies a request.
//! Independent interface tests: generated release-profile-shaped fixture only.
//! Analytic deSitter shape=(1+zHEL)*zHD, C=[[4,1],[1,9]] gives
//! q=(r0-r1)^2/11 and offset=(8*r0+3*r1)/11. No inference claim.
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
fn scratch() -> Scratch {
    let p = std::env::temp_dir().join(format!(
        "irred-supernova-{}-{}",
        std::process::id(),
        SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    fs::create_dir(&p).unwrap();
    fs::write(p.join("table.dat"),"CID IDSURVEY zHD zCMB zHEL m_b_corr EXTRA\nA 51 0.009 0.009 0.009 17 excluded\nA 52 0.1 0.1 0.1 2 generated\nB 51 0.2 0.2 0.2 -3 generated\n").unwrap();
    // Full raw matrix intentionally invalid; selected principal block valid.
    fs::write(p.join("cov.dat"), "3\n-5 8 7 6 4 1 0 1 9\n").unwrap();
    Scratch(p)
}
fn request(s: &Scratch) -> Value {
    json!({"schema_version":1,"operation":"supernova.profile_batch.v1","observations":{"schema_version":1,"operation":"observations.prepare.v1","table":s.0.join("table.dat"),"uncertainty":s.0.join("cov.dat"),"metadata":{"profile":"pantheon_plus_released_v1","role":"released_fitted_summary","unit":"magnitude","calibration":"unknown","uncertainty":"covariance","uncertainty_unit":"magnitude_squared","component":"total"},"selection":"pantheon_zhd_gt_001","ordering_provenance":"generated interface fixture supplied axis order; not actual release or independently verified matrix linkage","resources":{"maximum_asset_bytes":10000,"maximum_rows":3,"maximum_matrix_elements":9,"maximum_string_bytes":4096}},"models":[{"model":"flat_lcdm_late_v1","omega_m":0,"constant_q":0},{"model":"constant_q_flat_v1","omega_m":0,"constant_q":-1}],"policy":{"arithmetic":"longdouble_cpu_v1","include_residual_arrays":true,"maximum_models":4,"maximum_source_rows":3,"maximum_matrix_elements":9,"maximum_array_elements":100,"maximum_native_output_bytes":100000,"maximum_total_evaluations":100000,"maximum_evaluations_per_integral":10000,"maximum_depth":24,"absolute_tolerance":1e-12,"relative_tolerance":1e-12,"maximum_forward_sensitivity":1e-10}})
}
fn run(s: &Scratch, v: &Value) -> (Value, i32) {
    let p = s.0.join("request.json");
    fs::write(&p, serde_json::to_vec(v).unwrap()).unwrap();
    let o = Command::new(env!("CARGO_BIN_EXE_irred"))
        .arg("run")
        .arg(p)
        .arg(s.0.join("store")).args(["--assurance","qualified"])
        .output()
        .unwrap();
    (
        serde_json::from_slice(&o.stdout).unwrap_or_else(|_| {
            panic!(
                "stdout={} stderr={}",
                String::from_utf8_lossy(&o.stdout),
                String::from_utf8_lossy(&o.stderr)
            )
        }),
        o.status.code().unwrap(),
    )
}
fn sn_v2(s: &Scratch) -> Value {
    let mut v = request(s);
    v["operation"] = json!("supernova.profile_batch.v2");
    v["models"] = json!([
      {"model":"flat_cpl_late_v1","omega_m":0.,"constant_q":0.,"w0":-1.,"wa":0.},
      {"model":"flat_lcdm_late_v1","omega_m":0.,"constant_q":0.,"w0":-1.,"wa":0.},
      {"model":"constant_q_flat_v1","omega_m":0.,"constant_q":-1.,"w0":-1.,"wa":0.}]);
    v
}
fn bg_v2() -> Value {
    json!({"schema_version":1,"operation":"background.parameter_query_batch.v2",
 "parameters":[{"model":"flat_cpl_late_v1","h0_km_s_mpc":70.,"omega_m":0.3,"constant_q":0.,"w0":-1.,"wa":0.},{"model":"flat_lcdm_late_v1","h0_km_s_mpc":70.,"omega_m":0.3,"constant_q":0.,"w0":-1.,"wa":0.}],
 "queries":[{"z_expansion":0.,"z_observer":0.,"convention":"geometric_same_redshift"},{"z_expansion":1.,"z_observer":1.,"convention":"geometric_same_redshift"},{"z_expansion":1.,"z_observer":0.9,"convention":"released_zhd_zhel"}],
 "policy":{"maximum_parameters":4,"maximum_queries":3,"maximum_slots":12,"maximum_native_output_bytes":1048576,"maximum_total_evaluations":2000000,"maximum_evaluations_per_integral":100000,"maximum_depth":30,"absolute_tolerance":1e-13,"relative_tolerance":1e-12}})
}
#[test]
fn background_v2_explicit_sources_zero_lambda_and_version_strictness() {
    let s = scratch();
    let req = bg_v2();
    let (r, exit) = run(&s, &req);
    assert_eq!(exit, 6);
    assert_eq!(r["receipt"]["accepted"], false);
    assert_eq!(r["receipt"]["execution"], "completed");
    let rows = r["result"]["rows"].as_array().unwrap();
    assert_eq!(rows.len(), 6);
    for i in 0..3 {
        assert_eq!(rows[i]["result"], rows[i + 3]["result"]);
        assert_eq!(
            rows[i]["identity"]["source_parameters"],
            req["parameters"][0]
        );
        assert_eq!(
            rows[i + 3]["identity"]["source_parameters"],
            req["parameters"][1]
        );
    }
    assert_eq!(rows[0]["result"]["radial_mpc"], 0.);
    assert_eq!(rows[0]["result"]["jerk"], 1.);
    let mut legacy = req.clone();
    legacy["operation"] = json!("background.parameter_query_batch.v1");
    legacy["parameters"] =
        json!([{"model":"flat_lcdm_late_v1","h0_km_s_mpc":70.,"omega_m":0.3,"constant_q":0.}]);
    let (old, oldexit) = run(&s, &legacy);
    assert_eq!(oldexit, 6);
    for i in 0..3 {
        assert_eq!(rows[i + 3]["result"], old["result"]["rows"][i]["result"]);
    }
    let mut varied = req.clone();
    varied["parameters"][0]["w0"] = json!(-0.9);
    varied["parameters"][0]["wa"] = json!(0.4);
    let (a, exit) = run(&s, &varied);
    assert_eq!(exit, 6);
    assert!(
        (a["result"]["rows"][0]["result"]["deceleration_q"]
            .as_f64()
            .unwrap()
            + 0.445)
            .abs()
            < 1e-12
    );
    assert!((a["result"]["rows"][0]["result"]["jerk"].as_f64().unwrap() - 1.1365).abs() < 1e-12);
    assert_ne!(
        r["receipt"]["scientific_specification_digest"],
        a["receipt"]["scientific_specification_digest"]
    );
    for field in ["w0", "wa"] {
        let mut missing = req.clone();
        missing["parameters"][0]
            .as_object_mut()
            .unwrap()
            .remove(field);
        assert_eq!(run(&s, &missing).1, 2);
    }
    let mut forbidden = legacy.clone();
    forbidden["parameters"][0]["model"] = json!("flat_cpl_late_v1");
    assert_eq!(run(&s, &forbidden).1, 2);
    let mut unknown = req.clone();
    unknown["parameters"][0]["w_runtime_expression"] = json!("-1");
    assert_eq!(run(&s, &unknown).1, 2);
}
#[test]
fn supernova_v2_model_fields_and_legacy_profile_parity() {
    let s = scratch();
    let req = sn_v2(&s);
    let (r, exit) = run(&s, &req);
    assert_eq!(exit, 6);
    assert_eq!(r["receipt"]["accepted"], false);
    assert_eq!(r["receipt"]["outputs"][0]["numerical"], "checks_passed");
    assert_eq!(
        r["receipt"]["outputs"][0]["validation_coverage"],
        "named original-input four-point CPL regression; request applicability not established"
    );
    let rows = r["result"]["calculation"]["rows"].as_array().unwrap();
    assert_eq!(rows.len(), 3);
    for i in 0..3 {
        assert_eq!(rows[i]["identity"]["source_parameters"], req["models"][i]);
        assert_eq!(rows[i]["result"], rows[0]["result"]);
        assert_eq!(rows[i]["result"]["density_status"], "not_applicable");
        assert_eq!(rows[i]["result"]["normalization_status"], "not_applicable");
    }
    let (old, oldexit) = run(&s, &request(&s));
    assert_eq!(oldexit, 6);
    assert_eq!(
        rows[0]["result"],
        old["result"]["calculation"]["rows"][0]["result"]
    );
    let mut varied = req.clone();
    varied["models"][0]["omega_m"] = json!(0.3);
    varied["models"][0]["w0"] = json!(-0.9);
    varied["models"][0]["wa"] = json!(0.4);
    let (a, exit) = run(&s, &varied);
    assert_eq!(exit, 6);
    assert_eq!(
        a["result"]["calculation"]["rows"][0]["identity"]["source_parameters"],
        varied["models"][0]
    );
    assert_ne!(
        r["receipt"]["scientific_specification_digest"],
        a["receipt"]["scientific_specification_digest"]
    );
    for field in ["w0", "wa"] {
        let mut missing = req.clone();
        missing["models"][0].as_object_mut().unwrap().remove(field);
        assert_eq!(run(&s, &missing).1, 2);
    }
    let mut legacy = request(&s);
    legacy["models"][0]["model"] = json!("flat_cpl_late_v1");
    assert_eq!(run(&s, &legacy).1, 2);
}
#[test]
fn v2_failures_keep_attempts_omit_payload_and_global_work_is_bounded() {
    let s = scratch();
    let mut bg = bg_v2();
    bg["parameters"][1]["wa"] = json!(0.1);
    let (r, exit) = run(&s, &bg);
    assert_ne!(exit, 0);
    let row = &r["result"]["rows"][3];
    assert_eq!(row["identity"]["source_parameters"]["wa"], 0.1);
    assert_eq!(row["result"]["kind"], "failure");
    assert!(row["result"]["radial_mpc"].is_null());
    let mut bg = bg_v2();
    bg["policy"]["maximum_total_evaluations"] = json!(3);
    let (r, _) = run(&s, &bg);
    assert!(r["result"]["evaluations"].as_u64().unwrap() <= 3);
    assert_eq!(r["result"]["rows"][1]["result"]["kind"], "failure");
    assert_eq!(r["result"]["rows"][4]["result"]["kind"], "failure");
    let mut sn = sn_v2(&s);
    sn["models"][1]["wa"] = json!(0.1);
    let (r, _) = run(&s, &sn);
    let row = &r["result"]["calculation"]["rows"][1];
    assert_eq!(row["identity"]["source_parameters"]["wa"], 0.1);
    assert_eq!(row["result"]["kind"], "failure");
    assert!(row["result"]["quadratic"].is_null());
    assert!(row["result"]["arrays"].is_null());
    assert_eq!(r["receipt"]["execution"], "completed");
    assert_eq!(r["receipt"]["outputs"][0]["numerical"], "failed");
    let mut sn = sn_v2(&s);
    sn["policy"]["maximum_total_evaluations"] = json!(3);
    let (r, _) = run(&s, &sn);
    assert!(r["result"]["calculation"]["evaluations"].as_u64().unwrap() <= 3);
    let mut sn = sn_v2(&s);
    sn["policy"]["maximum_array_elements"] = json!(1);
    assert_eq!(run(&s, &sn).1, 2);
    let mut sn = sn_v2(&s);
    sn["policy"]["maximum_forward_sensitivity"] = json!(1e-30);
    let (r, _) = run(&s, &sn);
    assert_eq!(
        r["result"]["calculation"]["source"]["preparation_numerical_status"],
        "conditioning_budget_exceeded"
    );
    assert_eq!(r["receipt"]["accepted"], false);
}
