//! Independent retained BAO CLI/reader review. Synthetic controls and generated
//! release-format bytes are not published observations. Analytic q=.6/det35
//! validates normalized terms; wrapper equivalence remains boundary evidence.
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
        "irred-bao-piecewise-peer-{}-{}",
        std::process::id(),
        SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    fs::create_dir(&p).unwrap();
    Scratch(p)
}
fn hash(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}
fn policy() -> Value {
    json!({"arithmetic":"longdouble_cpu_v1","include_predictions":true,"include_residuals":true,"maximum_models":4,"maximum_rows":16,"maximum_matrix_elements":256,"maximum_string_bytes":4096,"maximum_native_bytes":1000000,"maximum_array_elements":256,"maximum_native_output_bytes":1000000,"maximum_queries":16,"maximum_total_segment_visits":2000,"maximum_forward_sensitivity":1e-10})
}
fn provenance(source: &mut Value) {
    source["ordering_provenance"] =
        json!("supplied control row order, covariance linkage not independently verified");
    source["calibration_provenance"] = json!("synthetic fixture; released calibration unknown");
    source["dependence_provenance"] =
        json!("declared correlated covariance, cross-probe dependence unknown");
    source["redshift_convention"] = json!("P01/released-effective-redshift/v1");
    source["ruler_convention"] = json!("P01/free-H0rd-km-s-no-early-physics/v1");
    source["computational_h0_convention"] = json!("P01/computational-H0-fixed-70-km-s-Mpc/v1");
}
fn request() -> Value {
    let mut source = json!({"profile":"synthetic_inline_v1","rows":[{"id":"a","z":0.1,"observable":"DM_over_rs","value":299792.458/10000.0*0.1+1.0},{"id":"b","z":0.2,"observable":"DH_over_rs","value":299792.458/10000.0+2.0}],"covariance":[4,1,1,9],"covariance_axis_ids":["a","b"]});
    provenance(&mut source);
    json!({"schema_version":1,"operation":"bao.piecewise_gaussian_batch.v1","source":source,"prepare_policy":policy(),"evaluation_policy":policy(),"models":[{"q":[-1.0,-1.0,-1.0,-1.0,-1.0],"h0_rd_km_s":10000.0},{"q":[-1.0,-1.0,-1.0,-1.0,-1.0],"h0_rd_km_s":10000.0}]})
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
    let r = serde_json::from_slice(&o.stdout).unwrap_or_else(|_| {
        panic!(
            "{} / {}",
            String::from_utf8_lossy(&o.stdout),
            String::from_utf8_lossy(&o.stderr)
        )
    });
    (r, o.status.code().unwrap())
}
fn c(v: &Value) -> &Value {
    &v["result"]["calculation"]
}
fn object(s: &Scratch, d: &str) -> Value {
    let b = fs::read(s.0.join("store/objects").join(d)).unwrap();
    assert_eq!(hash(&b), d);
    serde_json::from_slice(&b).unwrap()
}
#[test]
fn normalized_analytic_source_and_immutable_unqualified_records() {
    let s = scratch();
    let req = request();
    let (r, code) = run(&s, &req);
    assert_eq!(code, 6, "{r}");
    assert_eq!(r["receipt"]["accepted"], false);
    assert_eq!(r["receipt"]["execution"], "completed");
    assert_eq!(r["receipt"]["outputs"][0]["numerical"], "checks_passed");
    assert_eq!(c(&r)["interpretation_status"], "unqualified");
    assert_eq!(r["result"]["source"]["role"], "synthetic_control");
    assert_eq!(r["result"]["source"]["computational_h0_km_s_mpc"], 70);
    for (i, row) in c(&r)["rows"].as_array().unwrap().iter().enumerate() {
        assert_eq!(row["identity"]["source_parameters"], req["models"][i]);
        assert_eq!(row["identity"]["arithmetic_id"], "F02/longdouble-cpu/v1");
        assert_eq!(row["result"]["kind"], "finite");
        let q = row["result"]["quadratic"].as_f64().unwrap();
        let det = row["result"]["log_determinant"].as_f64().unwrap();
        let norm = row["result"]["normalization"].as_f64().unwrap();
        let density = row["result"]["log_density"].as_f64().unwrap();
        assert!((q - 0.6).abs() < 1e-12);
        assert!((det - 35f64.ln()).abs() < 1e-12);
        assert!((norm - 2.0 * (2.0 * std::f64::consts::PI).ln()).abs() < 1e-12);
        assert!((density + 0.5 * (q + det + norm)).abs() < 1e-12);
        assert!(row["result"].get("relative_profile_score").is_none());
        assert_eq!(row["result"]["predictions"].as_array().unwrap().len(), 2);
    }
    let digest = r["result"]["source"]["full_covariance"]["object_digest"]
        .as_str()
        .unwrap();
    let matrix = object(&s, digest);
    assert_eq!(matrix["axis_ids"], json!(["a", "b"]));
    assert_eq!(matrix["values"], json!([4.0, 1.0, 1.0, 9.0]));
    let spec = object(
        &s,
        r["receipt"]["scientific_specification_digest"]
            .as_str()
            .unwrap(),
    );
    assert_eq!(spec["models"], req["models"]);
    assert_eq!(spec["prepare_policy"], req["prepare_policy"]);
    assert_eq!(spec["evaluation_policy"], req["evaluation_policy"]);
    object(&s, r["receipt"]["output_digest"].as_str().unwrap());
    let (again, again_code) = run(&s, &req);
    assert_eq!(again_code, 6);
    assert_eq!(
        again["receipt"]["scientific_specification_digest"],
        r["receipt"]["scientific_specification_digest"]
    );
    assert_ne!(again["receipt"]["attempt_id"], r["receipt"]["attempt_id"]);
}
#[test]
fn phase_quotas_mixed_attempts_and_strict_analytic_packet() {
    let s = scratch();
    let mut req = request();
    req["evaluation_policy"]["maximum_rows"] = json!(0);
    req["evaluation_policy"]["maximum_matrix_elements"] = json!(0);
    req["evaluation_policy"]["maximum_string_bytes"] = json!(0);
    req["evaluation_policy"]["maximum_native_bytes"] = json!(0);
    let (r, code) = run(&s, &req);
    assert_eq!(code, 6, "{r}");
    req["models"][1]["q"][0] = json!(3.0);
    let (r, code) = run(&s, &req);
    assert_eq!(code, 2, "{r}");
    assert_eq!(r["receipt"]["execution"], "completed");
    assert_eq!(r["receipt"]["accepted"], false);
    let rows = c(&r)["rows"].as_array().unwrap();
    assert_eq!(rows.len(), 2);
    assert_eq!(rows[0]["result"]["kind"], "finite");
    assert_ne!(rows[1]["result"]["kind"], "finite");
    assert_eq!(rows[1]["identity"]["source_parameters"], req["models"][1]);
    for key in ["predictions", "residuals", "log_density", "quadratic"] {
        assert!(rows[1]["result"].get(key).is_none());
    }
    let mut req = request();
    req["evaluation_policy"]["maximum_total_segment_visits"] = json!(3);
    let (r, code) = run(&s, &req);
    assert_eq!(code, 2, "{r}");
    assert_eq!(c(&r)["rows"][0]["result"]["kind"], "finite");
    assert_ne!(c(&r)["rows"][1]["result"]["kind"], "finite");
    for field in [
        "maximum_models",
        "maximum_array_elements",
        "maximum_native_output_bytes",
    ] {
        let mut req = request();
        req["evaluation_policy"][field] = json!(0);
        req["evaluation_policy"]["include_predictions"] = json!(false);
        req["evaluation_policy"]["include_residuals"] = json!(false);
        let (r, code) = run(&s, &req);
        assert_eq!(code, 2, "{r}");
        assert!(c(&r)["rows"].as_array().unwrap().is_empty());
        assert_eq!(r["receipt"]["execution"], "completed");
    }
    let mut req = request();
    req["models"] = json!([]);
    let (r, code) = run(&s, &req);
    assert_eq!(code, 6, "{r}");
    assert!(c(&r)["rows"].as_array().unwrap().is_empty());
    for kind in 0..4 {
        let mut req = request();
        match kind {
            0 => {
                req["models"][0]["q"] = json!([0.0, 0.0]);
            }
            1 => {
                req["evaluation_policy"]["background"] = json!({});
            }
            2 => {
                req["models"][0]["epsilon"] = json!(0.0);
            }
            _ => {
                req["evaluation_policy"]
                    .as_object_mut()
                    .unwrap()
                    .remove("arithmetic");
            }
        }
        let (r, code) = run(&s, &req);
        assert_eq!(code, 2, "{r}");
    }
}
