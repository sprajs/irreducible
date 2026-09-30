//! Independent interface tests: generated release-profile-shaped fixture only.
//! Analytic deSitter shape=(1+zHEL)*zHD, C=[[4,1],[1,9]] gives
//! q=(r0-r1)^2/11 and offset=(8*r0+3*r1)/11. No inference claim.
use serde_json::{Value, json};
use sha2::{Digest, Sha256};
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
    fs::write(p.join("table.dat"),"CID IDSURVEY zHD zCMB zHEL m_b_corr EXTRA\nA 51 .009 .009 .009 17 excluded\nA 52 .1 .1 .1 2 generated\nB 51 .2 .2 .2 -3 generated\n").unwrap();
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
        .arg(s.0.join("store"))
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
#[test]
fn analytic_scores_order_assets_and_unqualified_status() {
    let s = scratch();
    let req = request(&s);
    let (r, exit) = run(&s, &req);
    assert_eq!(exit, 6);
    assert_eq!(r["receipt"]["accepted"], false);
    assert_eq!(r["receipt"]["execution"], "completed");
    assert_eq!(r["receipt"]["numerical"], "not_assessed");
    assert_eq!(r["receipt"]["outputs"][0]["numerical"], "checks_passed");
    let c = &r["result"]["calculation"];
    let o = &r["result"]["observations"];
    assert_eq!(c["kind"], "finite");
    assert_eq!(c["interpretation_status"], "unqualified");
    assert_eq!(c["execution_status"], "completed");
    assert_eq!(c["source"]["selected_source_indices"], json!([1, 2]));
    assert_eq!(c["source"]["z_expansion"], json!([0.1, 0.2]));
    assert_eq!(c["source"]["z_observer"], json!([0.1, 0.2]));
    assert_eq!(
        c["source"]["ordered_ids"],
        json!([o["measurement_ids"][1], o["measurement_ids"][2]])
    );
    assert_eq!(c["source"]["matrix_validation_assessed"], true);
    assert_eq!(c["source"]["table_identity"], o["table_digest"]);
    for (i, row) in c["rows"].as_array().unwrap().iter().enumerate() {
        assert_eq!(row["identity"]["model_index"], i);
        let x = &row["result"];
        assert_eq!(x["kind"], "finite");
        assert!((x["quadratic"].as_f64().unwrap() - 4.073716198571997).abs() < 1e-10);
        assert!((x["offset_coefficient"].as_f64().unwrap() - 4.967374906181537).abs() < 1e-10);
        assert!((x["relative_profile_score"].as_f64().unwrap() + 2.0368580992859986).abs() < 1e-10);
        assert_eq!(x["density_status"], "not_applicable");
        assert_eq!(x["normalization_status"], "not_applicable");
        assert_eq!(x["arrays"]["shape_magnitudes"].as_array().unwrap().len(), 2);
        assert_eq!(row["diagnostics"]["numerical_status"], "ok");
    }
    for name in ["table", "uncertainty"] {
        let bytes = fs::read(req["observations"][name].as_str().unwrap()).unwrap();
        let digest = format!("{:x}", Sha256::digest(&bytes));
        assert_eq!(
            fs::read(s.0.join("store/objects").join(digest)).unwrap(),
            bytes
        );
    }
    let matrix_digest = o["full_uncertainty_matrix"]["object_digest"]
        .as_str()
        .unwrap();
    let bytes = fs::read(s.0.join("store/objects").join(matrix_digest)).unwrap();
    assert_eq!(format!("{:x}", Sha256::digest(&bytes)), matrix_digest);
    let matrix: Value = serde_json::from_slice(&bytes).unwrap();
    assert_eq!(
        matrix["values"],
        json!([-5., 8., 7., 6., 4., 1., 0., 1., 9.])
    );
}
#[test]
fn failures_caps_causes_and_mixed_slots() {
    let s = scratch();
    let mut v = request(&s);
    v["models"]
        .as_array_mut()
        .unwrap()
        .push(json!({"model":"flat_lcdm_late_v1","omega_m":2,"constant_q":0}));
    let (r, exit) = run(&s, &v);
    assert_eq!(exit, 2);
    assert_eq!(r["receipt"]["execution"], "completed");
    assert_eq!(r["receipt"]["outputs"][0]["numerical"], "failed");
    let rows = &r["result"]["calculation"]["rows"];
    assert_eq!(rows[0]["result"]["kind"], "finite");
    assert_eq!(rows[2]["result"]["kind"], "failure");
    assert!(rows[2]["result"].get("quadratic").is_none());
    assert!(rows[2]["result"].get("arrays").is_none());
    v = request(&s);
    v["policy"]["maximum_total_evaluations"] = json!(3);
    let (r, exit) = run(&s, &v);
    assert_eq!(exit, 2);
    let c = &r["result"]["calculation"];
    assert!(c["evaluations"].as_u64().unwrap() <= 3);
    for row in c["rows"].as_array().unwrap() {
        assert_eq!(row["result"]["kind"], "failure");
        assert_eq!(row["diagnostics"]["numerical_status"], "work_limit");
    }
    v = request(&s);
    v["policy"]["maximum_forward_sensitivity"] = json!(1e-30);
    let (r, exit) = run(&s, &v);
    assert_eq!(exit, 2);
    assert_eq!(
        r["result"]["calculation"]["source"]["preparation_numerical_status"],
        "conditioning_budget_exceeded"
    );
    for key in ["maximum_array_elements", "maximum_native_output_bytes"] {
        v = request(&s);
        v["policy"][key] = json!(1);
        let (_, exit) = run(&s, &v);
        assert_eq!(exit, 2);
    }
    v = request(&s);
    v["models"] = json!([]);
    let (empty, exit) = run(&s, &v);
    assert_eq!(exit, 6);
    assert_eq!(empty["result"]["calculation"]["rows"], json!([]));
    assert_eq!(empty["result"]["calculation"]["evaluations"], 0);
    v = request(&s);
    v["policy"]["include_residual_arrays"] = json!(false);
    let (r, exit) = run(&s, &v);
    assert_eq!(exit, 6);
    assert_eq!(
        r["result"]["calculation"]["rows"][0]["result"]["arrays"]["shape_magnitudes"],
        json!([])
    );
}
#[test]
fn strict_source_identity_and_changed_same_path_inputs() {
    let s = scratch();
    for key in ["w0", "unknown"] {
        let mut v = request(&s);
        v["models"][0][key] = json!(0);
        let (_, exit) = run(&s, &v);
        assert_eq!(exit, 2);
    }
    let mut v = request(&s);
    v["observations"]["metadata"]["unit"] = json!("physical_length");
    let (_, exit) = run(&s, &v);
    assert_eq!(exit, 2);
    v = request(&s);
    v["observations"]["ordering_provenance"] = json!("");
    let (_, exit) = run(&s, &v);
    assert_eq!(exit, 2);
    v = request(&s);
    v["policy"]["arithmetic"] = json!("unknown");
    let (_, exit) = run(&s, &v);
    assert_eq!(exit, 2);
    v = request(&s);
    let (a, exit) = run(&s, &v);
    assert_eq!(exit, 6);
    let digest = a["result"]["observations"]["table_digest"]
        .as_str()
        .unwrap()
        .to_owned();
    let original = fs::read(s.0.join("store/objects").join(&digest)).unwrap();
    let changed = String::from_utf8(original.clone())
        .unwrap()
        .replace(".1 .1 .1 2 generated", ".1 .1 .1 3 generated");
    fs::write(s.0.join("table.dat"), changed).unwrap();
    let (b, exit) = run(&s, &v);
    assert_eq!(exit, 6);
    assert_ne!(
        a["result"]["observations"]["table_digest"],
        b["result"]["observations"]["table_digest"]
    );
    assert_ne!(
        a["receipt"]["scientific_specification_digest"],
        b["receipt"]["scientific_specification_digest"]
    );
    assert_eq!(
        fs::read(s.0.join("store/objects").join(digest)).unwrap(),
        original
    );
    let matrix_digest = b["result"]["observations"]["uncertainty_digest"]
        .as_str()
        .unwrap()
        .to_owned();
    let matrix_original = fs::read(s.0.join("store/objects").join(&matrix_digest)).unwrap();
    fs::write(s.0.join("cov.dat"), "3\n-5 8 7 6 5 1 0 1 9\n").unwrap();
    let (c, exit) = run(&s, &v);
    assert_eq!(exit, 6);
    assert_ne!(
        b["result"]["observations"]["uncertainty_digest"],
        c["result"]["observations"]["uncertainty_digest"]
    );
    assert_ne!(
        b["receipt"]["scientific_specification_digest"],
        c["receipt"]["scientific_specification_digest"]
    );
    assert_ne!(
        b["result"]["calculation"]["rows"][0]["result"]["quadratic"],
        c["result"]["calculation"]["rows"][0]["result"]["quadratic"]
    );
    assert_eq!(
        fs::read(s.0.join("store/objects").join(matrix_digest)).unwrap(),
        matrix_original
    );
    v["policy"]["arithmetic"] = json!("binary64_legacy_v1");
    let (d, exit) = run(&s, &v);
    assert_eq!(exit, 6);
    assert_eq!(
        c["result"]["observations"]["table_digest"],
        d["result"]["observations"]["table_digest"]
    );
    assert_ne!(
        c["receipt"]["scientific_specification_digest"],
        d["receipt"]["scientific_specification_digest"]
    );
    assert_ne!(
        c["result"]["calculation"]["source"]["arithmetic_id"],
        d["result"]["calculation"]["source"]["arithmetic_id"]
    );
}
