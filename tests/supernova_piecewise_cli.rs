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
    fs::write(p.join("table.dat"),"CID IDSURVEY zHD zCMB zHEL m_b_corr EXTRA\nA 51 0.009 0.009 0.009 17 excluded\nA 52 0.1 0.1 0.1 2 generated\nB 51 0.2 0.2 0.2 -3 generated\n").unwrap();
    // Full raw matrix intentionally invalid; selected principal block valid.
    fs::write(p.join("cov.dat"), "3\n-5 8 7 6 4 1 0 1 9\n").unwrap();
    Scratch(p)
}
fn request(s: &Scratch) -> Value {
    json!({"schema_version":1,"operation":"supernova.piecewise_profile_batch.v1","observations":{"schema_version":1,"operation":"observations.prepare.v1","table":s.0.join("table.dat"),"uncertainty":s.0.join("cov.dat"),"metadata":{"profile":"pantheon_plus_released_v1","role":"released_fitted_summary","unit":"magnitude","calibration":"unknown","uncertainty":"covariance","uncertainty_unit":"magnitude_squared","component":"total"},"selection":"pantheon_zhd_gt_001","ordering_provenance":"generated interface fixture supplied axis order; not actual release or independently verified matrix linkage","resources":{"maximum_asset_bytes":10000,"maximum_rows":3,"maximum_matrix_elements":9,"maximum_string_bytes":4096}},"models":[{"q":[-1.,-1.,-1.,-1.,-1.]},{"q":[0.,0.,0.,0.,0.]}],"policy":{"arithmetic":"longdouble_cpu_v1","include_residual_arrays":true,"maximum_models":4,"maximum_source_rows":3,"maximum_matrix_elements":9,"maximum_array_elements":100,"maximum_native_output_bytes":100000,"maximum_queries":3,"maximum_total_segment_visits":100,"maximum_forward_sensitivity":1e-10}})
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
fn analytic_source_and_immutable_status() {
    let s = scratch();
    let req = request(&s);
    let (v, e) = run(&s, &req);
    assert_eq!(e, 6);
    let c = &v["result"]["calculation"];
    let o = &v["result"]["observations"];
    assert_eq!(c["kind"], "finite");
    assert_eq!(c["batch_status"], "ok");
    assert_eq!(c["source"]["selected_source_indices"], json!([1, 2]));
    assert_eq!(
        c["source"]["ordered_ids"],
        json!([o["measurement_ids"][1], o["measurement_ids"][2]])
    );
    assert_eq!(c["source"]["z_expansion"], json!([0.1, 0.2]));
    let x = &c["rows"][0]["result"];
    let r0 = 2. - 5. * (0.11_f64).log10();
    let r1 = -3. - 5. * (0.24_f64).log10();
    assert!((x["quadratic"].as_f64().unwrap() - (r0 - r1).powi(2) / 11.).abs() < 1e-10);
    assert_eq!(x["density_status"], "not_applicable");
    assert_eq!(
        c["rows"][0]["identity"]["source_parameters"],
        req["models"][0]
    );
    assert_eq!(
        c["validation_coverage"]["profile"],
        "original_inputs_piecewise_six_points_v1"
    );
    assert_eq!(v["receipt"]["accepted"], false);
    assert_eq!(v["receipt"]["execution"], "completed");
    assert_eq!(v["receipt"]["outputs"][0]["numerical"], "checks_passed");
    assert_eq!(c["interpretation_status"], "unqualified");
    for key in ["scientific_specification_digest", "output_digest"] {
        let h = v["receipt"][key].as_str().unwrap();
        let b = fs::read(s.0.join("store/objects").join(h)).unwrap();
        assert_eq!(format!("{:x}", Sha256::digest(&b)), h);
    }
    for key in ["table", "uncertainty"] {
        let b = fs::read(req["observations"][key].as_str().unwrap()).unwrap();
        let h = format!("{:x}", Sha256::digest(&b));
        assert_eq!(fs::read(s.0.join("store/objects").join(h)).unwrap(), b);
    }
    let mut changed = req.clone();
    changed["models"][0]["q"][0] = json!(-0.5);
    let (w, _) = run(&s, &changed);
    assert_ne!(
        v["receipt"]["scientific_specification_digest"],
        w["receipt"]["scientific_specification_digest"]
    );
}
#[test]
fn analytic_global_work_and_owned_caps() {
    let s = scratch();
    let mut q = request(&s);
    q["policy"]["maximum_total_segment_visits"] = json!(3);
    let (v, e) = run(&s, &q);
    assert_eq!(e, 2);
    let c = &v["result"]["calculation"];
    assert_eq!(c["batch_status"], "ok");
    assert_eq!(c["segment_visits"], 3);
    assert_eq!(c["rows"][0]["result"]["kind"], "finite");
    assert_eq!(c["rows"][1]["result"]["kind"], "failure");
    assert!(c["rows"][1]["result"].get("quadratic").is_none());
    assert_eq!(v["receipt"]["outputs"][0]["numerical"], "failed");
    assert_eq!(v["receipt"]["execution"], "completed");
    assert_eq!(v["receipt"]["accepted"], false);
    assert_eq!(c["interpretation_status"], "unqualified");
    for field in [
        "maximum_models",
        "maximum_array_elements",
        "maximum_native_output_bytes",
    ] {
        let mut q = request(&s);
        q["policy"][field] = json!(0);
        let (v, e) = run(&s, &q);
        assert_eq!(e, 2);
        let c = &v["result"]["calculation"];
        assert_eq!(c["batch_status"], "work_limit");
        assert_eq!(c["rows"], json!([]));
        assert_eq!(c["segment_visits"], 0);
    }
    let mut q = request(&s);
    q["models"] = json!([]);
    q["policy"]["maximum_models"] = json!(0);
    let (v, e) = run(&s, &q);
    assert_eq!(e, 6);
    assert_eq!(v["result"]["calculation"]["batch_status"], "ok");
    let mut q = request(&s);
    q["models"][0]["q"][0] = json!(3.);
    let (v, e) = run(&s, &q);
    assert_eq!(e, 2);
    assert_eq!(
        v["result"]["calculation"]["rows"][0]["identity"]["source_parameters"],
        q["models"][0]
    );
    assert_eq!(
        v["result"]["calculation"]["rows"][0]["result"]["kind"],
        "failure"
    );
}
#[test]
fn strict_schema_and_export_suppression() {
    let s = scratch();
    for mode in 0..4 {
        let mut q = request(&s);
        match mode {
            0 => q["models"][0]["q"] = json!([0., 0.]),
            1 => q["policy"]["maximum_depth"] = json!(24),
            2 => {
                q["policy"].as_object_mut().unwrap().remove("arithmetic");
            }
            _ => q["policy"]["arithmetic"] = json!("unknown"),
        };
        let (v, e) = run(&s, &q);
        assert_eq!(e, 2);
        assert!(v["result"]["calculation"]["rows"].is_null());
    }
    let mut q = request(&s);
    q["policy"]["include_residual_arrays"] = json!(false);
    let (v, e) = run(&s, &q);
    assert_eq!(e, 6);
    assert_eq!(
        v["result"]["calculation"]["rows"][0]["result"]["arrays"]["shape_magnitudes"],
        json!([])
    );
    q["policy"]["maximum_array_elements"] = json!(15);
    let (v, e) = run(&s, &q);
    assert_eq!(e, 2);
    assert_eq!(v["result"]["calculation"]["batch_status"], "work_limit");
}
