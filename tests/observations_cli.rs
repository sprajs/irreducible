//! Exact structural/semantic cases, independently expected selections. No Rust science arithmetic.
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
    let s = Scratch(std::env::temp_dir().join(format!(
            "irred-observation-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        )));
    fs::create_dir(&s.0).unwrap();
    s
}
fn request(s: &Scratch) -> Value {
    json!({"schema_version":1,"operation":"observations.prepare.v1","table":s.0.join("table.dat"),"uncertainty":s.0.join("cov.dat"),"metadata":{"profile":"pantheon_plus_released_v1","role":"released_fitted_summary","unit":"magnitude","calibration":"unknown","uncertainty":"covariance","uncertainty_unit":"magnitude_squared","component":"total"},"selection":"pantheon_zhd_gt_001","ordering_provenance":"supplied release table row order; not independently verified","resources":{"maximum_asset_bytes":10000,"maximum_rows":3,"maximum_matrix_elements":9,"maximum_string_bytes":4096}})
}
fn fixture(s: &Scratch) {
    fs::write(s.0.join("table.dat"),"CID IDSURVEY zHD zCMB zHEL m_b_corr EXTRA\nA 51 .009 .02 .03 17 keep\nA 52 .020 .03 .04 18 keep2\nB 51 .010 .04 .05 19 keep3\n").unwrap();
    fs::write(s.0.join("cov.dat"), "3\n4 1.00000003 0 1 9 0 0 0 0\n").unwrap();
}
fn run(s: &Scratch, v: &Value) -> (Value, i32) {
    let path = s.0.join("request.json");
    fs::write(&path, serde_json::to_vec(v).unwrap()).unwrap();
    let out = Command::new(env!("CARGO_BIN_EXE_irred"))
        .arg("run")
        .arg(path)
        .arg(s.0.join("store"))
        .output()
        .unwrap();
    (
        serde_json::from_slice(&out.stdout)
            .unwrap_or_else(|_| panic!("{}", String::from_utf8_lossy(&out.stderr))),
        out.status.code().unwrap(),
    )
}
#[test]
fn released_semantics_order_masks_and_raw_retention() {
    let s = scratch();
    fixture(&s);
    let req = request(&s);
    let (r, exit) = run(&s, &req);
    assert_eq!(exit, 6);
    let o = &r["result"];
    assert_eq!(o["values"], json!([17., 18., 19.]));
    assert_eq!(o["event_ids"], json!(["A", "A", "B"]));
    assert_ne!(o["measurement_ids"][0], o["measurement_ids"][1]);
    assert_eq!(o["selected_mask"], json!([0, 1, 0]));
    assert_eq!(o["selected_source_indices"], json!([1]));
    let matrix_digest = o["full_uncertainty_matrix"]["object_digest"]
        .as_str()
        .unwrap();
    let matrix_bytes = fs::read(s.0.join("store/objects").join(matrix_digest)).unwrap();
    assert_eq!(
        format!("{:x}", Sha256::digest(&matrix_bytes)),
        matrix_digest
    );
    let matrix: Value = serde_json::from_slice(&matrix_bytes).unwrap();
    assert_eq!(
        matrix["values"],
        json!([4., 1.00000003, 0., 1., 9., 0., 0., 0., 0.])
    );
    assert_eq!(matrix["axis_ids"], o["measurement_ids"]);
    assert_eq!(matrix["source_asset_digest"], o["uncertainty_digest"]);
    assert_eq!(o["full_uncertainty_matrix"]["elements"], 9);
    assert_eq!(o["uncertainty_axis_ids"], o["measurement_ids"]);
    assert_eq!(o["metadata"], req["metadata"]);
    assert_eq!(o["original_fields"][1][6], "keep2");
    assert_eq!(o["quality_flags_available"], false);
    assert_eq!(o["ordering_verified_from_unlabeled_matrix"], false);
    assert_eq!(o["inference_independence"], "not_asserted");
    assert_eq!(r["receipt"]["accepted"], false);
    for name in ["table", "uncertainty"] {
        let bytes = fs::read(req[name].as_str().unwrap()).unwrap();
        let digest = format!("{:x}", Sha256::digest(&bytes));
        assert_eq!(
            fs::read(s.0.join("store/objects").join(&digest)).unwrap(),
            bytes
        );
        assert_eq!(o[format!("{name}_digest")], digest);
    }
}
#[test]
fn same_path_mutation_changes_scientific_identity_without_overwrite() {
    let s = scratch();
    fixture(&s);
    let req = request(&s);
    let (a, _) = run(&s, &req);
    let table = s.0.join("table.dat");
    let old = fs::read(&table).unwrap();
    let changed = String::from_utf8(old.clone())
        .unwrap()
        .replace("18 keep2", "28 keep2");
    assert_eq!(old.len(), changed.len());
    fs::write(&table, changed).unwrap();
    let (b, _) = run(&s, &req);
    assert_ne!(
        a["receipt"]["scientific_specification_digest"],
        b["receipt"]["scientific_specification_digest"]
    );
    assert_ne!(a["result"]["table_digest"], b["result"]["table_digest"]);
    assert_eq!(a["result"]["values"][1], 18.);
    assert_eq!(b["result"]["values"][1], 28.);
    let old_digest = a["result"]["table_digest"].as_str().unwrap();
    assert_eq!(
        fs::read(s.0.join("store/objects").join(old_digest)).unwrap(),
        old
    );
}
#[test]
fn invalid_semantics_and_resource_contract_fail_closed() {
    let s = scratch();
    fixture(&s);
    let req = request(&s);
    for (field, value) in [
        ("calibration", "not_applicable"),
        ("role", "posterior_summary"),
        ("uncertainty_unit", "inverse_magnitude_squared"),
    ] {
        let mut q = req.clone();
        q["metadata"][field] = json!(value);
        let (r, exit) = run(&s, &q);
        assert_ne!(exit, 0);
        assert_eq!(r["result"]["kind"], "failure");
        assert!(r["result"].get("selected_source_indices").is_none());
    }
    let mut q = req.clone();
    q["ordering_provenance"] = json!("");
    assert_eq!(run(&s, &q).0["result"]["kind"], "failure");
    q = req.clone();
    q["resources"]["maximum_rows"] = json!(2);
    assert_eq!(run(&s, &q).0["result"]["kind"], "failure");
    q = req.clone();
    q["source_selection"] = json!([1, 2, 1]);
    assert_eq!(run(&s, &q).0["result"]["kind"], "failure");
}

#[test]
fn same_path_matrix_mutation_changes_identity_and_preserves_original() {
    let s = scratch();
    fixture(&s);
    let req = request(&s);
    let (a, _) = run(&s, &req);
    let path = s.0.join("cov.dat");
    let old = fs::read(&path).unwrap();
    let changed = String::from_utf8(old.clone())
        .unwrap()
        .replace("4 1.00000003", "5 1.00000003");
    assert_eq!(old.len(), changed.len());
    fs::write(path, changed).unwrap();
    let (b, _) = run(&s, &req);
    assert_ne!(
        a["receipt"]["scientific_specification_digest"],
        b["receipt"]["scientific_specification_digest"]
    );
    assert_ne!(
        a["result"]["uncertainty_digest"],
        b["result"]["uncertainty_digest"]
    );
    assert_ne!(
        a["result"]["full_uncertainty_matrix"]["object_digest"],
        b["result"]["full_uncertainty_matrix"]["object_digest"]
    );
    assert_eq!(
        fs::read(
            s.0.join("store/objects")
                .join(a["result"]["uncertainty_digest"].as_str().unwrap())
        )
        .unwrap(),
        old
    );
}
#[test]
fn unused_nonfinite_source_is_distinct_from_missing_and_selected_failure() {
    let s = scratch();
    fixture(&s);
    let req = request(&s);
    let path = s.0.join("table.dat");
    let text = fs::read_to_string(&path)
        .unwrap()
        .replace(".009 .02 .03 17 keep", ".009 NaN .03 NaN keep");
    fs::write(&path, &text).unwrap();
    let (r, exit) = run(&s, &req);
    assert_eq!(exit, 6);
    assert_eq!(r["result"]["selected_source_indices"], json!([1]));
    assert_eq!(
        r["result"]["source_nonfinite"]["values"],
        json!([true, false, false])
    );
    assert_eq!(
        r["result"]["source_nonfinite"]["zcmb"],
        json!([true, false, false])
    );
    assert_eq!(r["result"]["source_missing"]["values"], json!([0, 0, 0]));
    assert_eq!(r["result"]["original_fields"][0][5], "NaN");
    let mut q = req;
    q["selection"] = json!("all");
    let (r, exit) = run(&s, &q);
    assert_eq!(exit, 2);
    assert_eq!(r["result"]["kind"], "failure");
    assert!(r["result"].get("values").is_none());
}
