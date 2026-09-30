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
        "irred-bao-peer-{}-{}",
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
    json!({"background":{"maximum_parameters":4,"maximum_queries":16,"maximum_slots":64,"maximum_native_output_bytes":1000000,"maximum_total_evaluations":2000000,"maximum_evaluations_per_integral":100000,"maximum_depth":30,"absolute_tolerance":1e-13,"relative_tolerance":1e-12},"arithmetic":"longdouble_cpu_v1","include_predictions":true,"include_residuals":true,"maximum_models":4,"maximum_rows":16,"maximum_matrix_elements":256,"maximum_string_bytes":4096,"maximum_native_bytes":1000000,"maximum_array_elements":256,"maximum_native_output_bytes":1000000,"maximum_forward_sensitivity":1e-10})
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
    json!({"schema_version":1,"operation":"bao.gaussian_batch.v1","source":source,"prepare_policy":policy(),"evaluation_policy":policy(),"models":[{"model":"constant_q_flat_v1","omega_m":0.0,"constant_q":-1.0,"w0":-1.0,"wa":0.0,"h0_rd_km_s":10000.0},{"model":"flat_lcdm_late_v1","omega_m":0.0,"constant_q":0.0,"w0":-1.0,"wa":0.0,"h0_rd_km_s":10000.0},{"model":"flat_cpl_late_v1","omega_m":0.0,"constant_q":0.0,"w0":-1.0,"wa":0.0,"h0_rd_km_s":10000.0}]})
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
fn normalized_terms_owned_source_and_unqualified_execution() {
    let s = scratch();
    let req = request();
    let (r, code) = run(&s, &req);
    assert_eq!(code, 6, "{r}");
    assert_eq!(r["result"]["source"]["role"], "synthetic_control");
    assert_eq!(r["result"]["source"]["computational_h0_km_s_mpc"], 70);
    assert_eq!(r["receipt"]["accepted"], false);
    assert_eq!(r["receipt"]["execution"], "completed");
    assert_eq!(r["receipt"]["outputs"][0]["numerical"], "checks_passed");
    assert_eq!(c(&r)["interpretation_status"], "unqualified");
    for (i, row) in c(&r)["rows"].as_array().unwrap().iter().enumerate() {
        assert_eq!(row["identity"]["source_parameters"], req["models"][i]);
        assert_eq!(row["identity"]["arithmetic_id"], "F02/longdouble-cpu/v1");
        assert_eq!(row["result"]["kind"], "finite");
        let q = row["result"]["quadratic"].as_f64().unwrap();
        let ld = row["result"]["log_determinant"].as_f64().unwrap();
        let norm = row["result"]["normalization"].as_f64().unwrap();
        let log = row["result"]["log_density"].as_f64().unwrap();
        assert!((q - 0.6).abs() < 1e-12);
        assert!((ld - 35f64.ln()).abs() < 1e-12);
        assert!((norm - 2.0 * (2.0 * std::f64::consts::PI).ln()).abs() < 1e-12);
        assert!((log + 0.5 * (q + ld + norm)).abs() < 1e-12);
        assert_eq!(row["result"]["predictions"].as_array().unwrap().len(), 2);
        assert_eq!(row["result"]["residuals"].as_array().unwrap().len(), 2);
        assert!(row["result"].get("relative_profile_score").is_none());
    }
    assert_eq!(c(&r)["preparation"]["policy"], req["prepare_policy"]);
    assert_eq!(c(&r)["evaluation_policy"], req["evaluation_policy"]);
    let md = r["result"]["source"]["full_covariance"]["object_digest"]
        .as_str()
        .unwrap();
    let matrix = object(&s, md);
    assert_eq!(matrix["values"], json!([4.0, 1.0, 1.0, 9.0]));
    assert_eq!(matrix["axis_ids"], json!(["a", "b"]));
    assert_eq!(
        matrix["ordering_assessment"],
        "supplied declaration; not independently verified from unlabeled bytes"
    );
    for key in ["scientific_specification_digest", "output_digest"] {
        object(&s, r["receipt"][key].as_str().unwrap());
    }
}
#[test]
fn mixed_failures_explicit_policy_global_work_and_payload_omission() {
    let s = scratch();
    let mut req = request();
    req["models"][1]["w0"] = json!(-0.9);
    let (r, code) = run(&s, &req);
    assert_eq!(code, 2, "{r}");
    assert_eq!(r["receipt"]["execution"], "completed");
    assert_eq!(r["receipt"]["outputs"][0]["numerical"], "failed");
    assert_eq!(c(&r)["rows"][0]["result"]["kind"], "finite");
    assert_eq!(
        c(&r)["rows"][1]["identity"]["source_parameters"],
        req["models"][1]
    );
    assert_eq!(c(&r)["rows"][1]["result"], json!({"kind":"failure"}));
    let mut req = request();
    req["evaluation_policy"]["maximum_forward_sensitivity"] = json!(1e-30);
    let (r, code) = run(&s, &req);
    assert_eq!(code, 2, "{r}");
    assert_eq!(
        c(&r)["rows"][0]["numerical_status"],
        "conditioning_budget_exceeded"
    );
    assert_eq!(c(&r)["rows"][0]["result"], json!({"kind":"failure"}));
    let mut req = request();
    req["evaluation_policy"]["background"]["maximum_total_evaluations"] = json!(20);
    let (r, code) = run(&s, &req);
    assert_eq!(code, 2, "{r}");
    assert!(c(&r)["evaluations"].as_u64().unwrap() <= 20);
    assert!(
        c(&r)["rows"]
            .as_array()
            .unwrap()
            .iter()
            .any(|x| x["numerical_status"] == "work_limit")
    );
    let mut req = request();
    req["evaluation_policy"]["include_predictions"] = json!(false);
    req["evaluation_policy"]["include_residuals"] = json!(false);
    let (r, code) = run(&s, &req);
    assert_eq!(code, 6, "{r}");
    assert_eq!(c(&r)["rows"][0]["result"]["predictions"], json!([]));
    req["evaluation_policy"]["maximum_array_elements"] = json!(11);
    assert_eq!(run(&s, &req).1, 2);
    let mut req = request();
    req["models"] = json!([]);
    let (r, code) = run(&s, &req);
    assert_eq!(code, 6, "{r}");
    assert_eq!(c(&r)["rows"], json!([]));
}
fn release_files(s: &Scratch) -> Value {
    let mut mean =
        String::from("# Generated format fixture, not original published observations\n");
    for i in 0..13 {
        let tag = if i == 11 {
            "DH_over_rs"
        } else if i == 12 {
            "DM_over_rs"
        } else if i % 2 == 0 {
            "DM_over_rs"
        } else {
            "DH_over_rs"
        };
        let z = if i >= 11 { 2.33 } else { 0.1 + i as f64 * 0.02 };
        mean += &format!("{z} {} {tag}\n", 10 + i);
    }
    fs::write(s.0.join("mean.txt"), mean).unwrap();
    let rows: Vec<String> = (0..13)
        .map(|i| {
            (0..13)
                .map(|j| if i == j { "1" } else { "0" })
                .collect::<Vec<_>>()
                .join(" ")
        })
        .collect();
    fs::write(s.0.join("cov.txt"), rows.join("\n") + "\n").unwrap();
    let mut req = request();
    req["source"] = json!({"profile":"desi_dr2_all_gccomb_13_v1","mean":s.0.join("mean.txt"),"covariance":s.0.join("cov.txt"),"maximum_asset_bytes":10000});
    provenance(&mut req["source"]);
    req
}
#[test]
fn strict_reader_ragged_equal_count_unknown_tags_and_source_shape() {
    let s = scratch();
    let req = release_files(&s);
    let (r, code) = run(&s, &req);
    assert_eq!(code, 6, "{r}");
    assert_eq!(
        r["result"]["source"]["rows"][11]["observable"],
        "DH_over_rs"
    );
    assert_eq!(
        r["result"]["source"]["rows"][12]["observable"],
        "DM_over_rs"
    );
    let correct = fs::read(s.0.join("cov.txt")).unwrap();
    let text = String::from_utf8(correct.clone()).unwrap();
    let mut lines: Vec<Vec<&str>> = text
        .lines()
        .map(|l| l.split_whitespace().collect())
        .collect();
    let item = lines[0].pop().unwrap();
    lines[1].push(item);
    assert_eq!(lines.iter().map(Vec::len).sum::<usize>(), 169);
    fs::write(
        s.0.join("cov.txt"),
        lines
            .iter()
            .map(|l| l.join(" "))
            .collect::<Vec<_>>()
            .join("\n"),
    )
    .unwrap();
    let (r, code) = run(&s, &req);
    assert_eq!(code, 2, "{r}");
    assert_eq!(r["result"]["error_id"], "INVALID_BAO_COVARIANCE_SHAPE");
    fs::write(s.0.join("cov.txt"), correct).unwrap();
    let bytes = fs::read(s.0.join("mean.txt")).unwrap();
    let changed = String::from_utf8(bytes.clone())
        .unwrap()
        .replacen("DM_over_rs", "UNKNOWN", 1);
    fs::write(s.0.join("mean.txt"), changed).unwrap();
    let (r, code) = run(&s, &req);
    assert_eq!(code, 2);
    assert_eq!(r["result"]["error_id"], "UNSUPPORTED_BAO_OBSERVABLE");
    fs::write(s.0.join("mean.txt"), bytes).unwrap();
    let mut bad = request();
    bad["source"]["covariance_axis_ids"] = json!(["b", "a"]);
    assert_eq!(run(&s, &bad).1, 2);
    let mut bad = request();
    bad["source"]["rows"][0]["unit"] = json!("magnitude");
    assert_eq!(run(&s, &bad).1, 2);
    let mut bad = request();
    bad["models"][0].as_object_mut().unwrap().remove("w0");
    assert_eq!(run(&s, &bad).1, 2);
    let mut bad = request();
    bad["prepare_policy"]["arithmetic"] = json!("unknown");
    assert_eq!(run(&s, &bad).1, 2);
}
#[test]
fn changed_same_path_assets_and_policy_change_identity_preserve_objects() {
    let s = scratch();
    let req = release_files(&s);
    let (old, code) = run(&s, &req);
    assert_eq!(code, 6, "{old}");
    let old_table = old["result"]["source"]["table_asset_sha256"]
        .as_str()
        .unwrap()
        .to_owned();
    let old_cov = old["result"]["source"]["covariance_asset_sha256"]
        .as_str()
        .unwrap()
        .to_owned();
    let old_table_bytes = fs::read(s.0.join("mean.txt")).unwrap();
    let old_cov_bytes = fs::read(s.0.join("cov.txt")).unwrap();
    let mut text = String::from_utf8(old_table_bytes.clone()).unwrap();
    text = text.replacen("0.1 10 DM_over_rs", "0.1 11 DM_over_rs", 1);
    fs::write(s.0.join("mean.txt"), text).unwrap();
    let (table_changed, code) = run(&s, &req);
    assert_eq!(code, 6, "{table_changed}");
    assert_ne!(
        old["receipt"]["scientific_specification_digest"],
        table_changed["receipt"]["scientific_specification_digest"]
    );
    assert_ne!(
        old_table,
        table_changed["result"]["source"]["table_asset_sha256"]
            .as_str()
            .unwrap()
    );
    let matrix = String::from_utf8(old_cov_bytes.clone())
        .unwrap()
        .replacen("1 0 0", "2 0 0", 1);
    fs::write(s.0.join("cov.txt"), matrix).unwrap();
    let (cov_changed, code) = run(&s, &req);
    assert_eq!(code, 6, "{cov_changed}");
    assert_ne!(
        table_changed["receipt"]["scientific_specification_digest"],
        cov_changed["receipt"]["scientific_specification_digest"]
    );
    assert_ne!(
        old_cov,
        cov_changed["result"]["source"]["covariance_asset_sha256"]
            .as_str()
            .unwrap()
    );
    assert_eq!(
        fs::read(s.0.join("store/objects").join(&old_table)).unwrap(),
        old_table_bytes
    );
    assert_eq!(
        fs::read(s.0.join("store/objects").join(&old_cov)).unwrap(),
        old_cov_bytes
    );
    let mut p = request();
    let (wide, code) = run(&s, &p);
    assert_eq!(code, 6);
    p["prepare_policy"]["arithmetic"] = json!("binary64_legacy_v1");
    p["evaluation_policy"]["arithmetic"] = json!("binary64_legacy_v1");
    let (legacy, code) = run(&s, &p);
    assert_eq!(code, 6, "{legacy}");
    assert_ne!(
        wide["receipt"]["scientific_specification_digest"],
        legacy["receipt"]["scientific_specification_digest"]
    );
    assert_eq!(
        c(&legacy)["rows"][0]["identity"]["arithmetic_id"],
        "F02/binary64-legacy/v1"
    );
}
