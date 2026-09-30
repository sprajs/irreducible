//! Independent fixed-bin transport checks; analytic de Sitter is a small
//! control. Native split-reference gates provide separate physics evidence.
use serde_json::{json, Value};
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
        "irred-piecewise-cli-{}-{}",
        std::process::id(),
        SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    fs::create_dir(&p).unwrap();
    Scratch(p)
}
fn query(z: f64) -> Value {
    json!({"z_expansion":z,"z_observer":z,"convention":"geometric_same_redshift"})
}
fn request() -> Value {
    json!({"schema_version":1,"operation":"background.piecewise_query_batch.v1","parameters":[{"h0_km_s_mpc":70.0,"q":[-1.0,-1.0,-1.0,-1.0,-1.0]},{"h0_km_s_mpc":100.0,"q":[-0.4,-0.4,-0.2,0.1,0.3]}],"queries":[query(0.0),query(0.1),query(0.3),query(2.5),{"z_expansion":0.7,"z_observer":0.71,"convention":"released_zhd_zhel"}],"policy":{"maximum_parameters":8,"maximum_queries":16,"maximum_slots":128,"maximum_native_output_bytes":1000000,"maximum_total_segment_visits":1000}})
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
    let r = serde_json::from_slice(&o.stdout)
        .unwrap_or_else(|_| panic!("{}", String::from_utf8_lossy(&o.stderr)));
    (r, o.status.code().unwrap())
}
fn object(s: &Scratch, d: &str) -> Value {
    let b = fs::read(s.0.join("store/objects").join(d)).unwrap();
    assert_eq!(format!("{:x}", Sha256::digest(&b)), d);
    serde_json::from_slice(&b).unwrap()
}
fn no_geometry(row: &Value) {
    assert!(row.get("geometry").is_none());
    let d = &row["derivatives"];
    assert_eq!(d["has_q"], false);
    assert_eq!(d["has_jerk"], false);
    assert_eq!(d["has_q0"], false);
    assert_eq!(d["q_convention"], "not_assessed");
    assert_eq!(d["jerk_availability"], "not_assessed");
    for k in ["assigned_q", "jerk", "q0_within_piecewise_model"] {
        assert!(d.get(k).is_none());
    }
}
#[test]
fn analytic_geometry_derivative_absence_source_and_records() {
    let s = scratch();
    let req = request();
    let (r, exit) = run(&s, &req);
    assert_eq!(exit, 6, "{r}");
    assert_eq!(r["receipt"]["execution"], "completed");
    assert_eq!(r["receipt"]["accepted"], false);
    assert_eq!(r["receipt"]["outputs"][0]["numerical"], "checks_passed");
    assert_eq!(r["receipt"]["outputs"][0]["interpretation"], "unqualified");
    assert_eq!(
        r["receipt"]["precision"],
        "binary64_storage_longdouble_analytic_intermediate"
    );
    let rows = r["result"]["slots"].as_array().unwrap();
    assert_eq!(rows.len(), 10);
    for (i, row) in rows.iter().enumerate() {
        assert_eq!(row["identity"]["parameter_index"], i / 5);
        assert_eq!(row["identity"]["query_index"], i % 5);
        assert_eq!(row["identity"]["parameters"], req["parameters"][i / 5]);
        assert_eq!(row["identity"]["query"], req["queries"][i % 5]);
        assert_eq!(row["status"], "ok");
        assert_eq!(row["numerical_status"], "ok");
        if i < 5 {
            let z = req["queries"][i]["z_expansion"].as_f64().unwrap();
            let g = &row["geometry"];
            assert_eq!(g["expansion_E"], 1.0);
            assert!((g["radial_integral"].as_f64().unwrap() - z).abs() < 2e-14);
            assert!((g["radial_mpc"].as_f64().unwrap() - 299792.458 / 70.0 * z).abs() < 1e-10);
            assert_eq!(row["derivatives"]["assigned_q"], -1.0);
            assert_eq!(row["derivatives"]["jerk"], 1.0);
        }
    }
    assert_eq!(
        rows[0]["derivatives"]["q_convention"],
        "right_limit_at_zero"
    );
    assert_eq!(rows[0]["derivatives"]["has_q0"], true);
    assert_eq!(rows[0]["geometry"]["luminosity_mpc"], 0.0);
    assert_eq!(
        rows[3]["derivatives"]["q_convention"],
        "left_limit_at_final_endpoint"
    );
    assert_eq!(rows[6]["derivatives"]["has_jerk"], true);
    assert_eq!(
        rows[6]["derivatives"]["jerk_availability"],
        "ordinary_within_bin"
    );
    assert_eq!(
        rows[7]["derivatives"]["q_convention"],
        "right_limit_at_internal_jump"
    );
    assert_eq!(rows[7]["derivatives"]["has_jerk"], false);
    assert!(rows[7]["derivatives"].get("jerk").is_none());
    assert_eq!(rows[7]["derivatives"]["assigned_q"], -0.2);
    let spec = object(
        &s,
        r["receipt"]["scientific_specification_digest"]
            .as_str()
            .unwrap(),
    );
    assert_eq!(spec["parameters"], req["parameters"]);
    assert_eq!(spec["policy"], req["policy"]);
    assert_eq!(spec["units"]["volume"], "Mpc^3/sr/redshift");
    assert_eq!(
        r["receipt"]["resource_budget"]["operation"],
        req["policy"]
    );
    object(&s, r["receipt"]["output_digest"].as_str().unwrap());
}
#[test]
fn global_atomic_admission_failed_work_and_attempted_sources() {
    let s = scratch();
    let mut req = request();
    req["queries"] = json!([query(1.0), query(0.1)]);
    req["policy"]["maximum_total_segment_visits"] = json!(6);
    let (r, exit) = run(&s, &req);
    assert_eq!(exit, 2, "{r}");
    assert_eq!(r["receipt"]["execution"], "completed");
    assert_eq!(r["receipt"]["outputs"][0]["numerical"], "failed");
    let rows = r["result"]["slots"].as_array().unwrap();
    assert_eq!(r["result"]["segments_processed"], 6);
    assert_eq!(r["result"]["batch_status"], "ok");
    assert_eq!(r["result"]["numerical_status"], "ok");
    assert_eq!(rows[2]["status"], "work_limit");
    assert_eq!(rows[2]["numerical_status"], "work_limit");
    assert_eq!(rows[2]["segments_processed"], 0);
    no_geometry(&rows[2]);
    assert_eq!(rows[3]["status"], "ok");
    assert_eq!(rows[3]["segments_processed"], 1);
    req["parameters"] = json!([{"h0_km_s_mpc":70.0,"q":[-1.0,-1.0,-1.0,-1.0,-1.0]}]);
    req["queries"] = json!([query(1e-200), query(0.1)]);
    req["policy"]["maximum_total_segment_visits"] = json!(1);
    let (r, exit) = run(&s, &req);
    assert_eq!(exit, 2);
    let rows = r["result"]["slots"].as_array().unwrap();
    assert_ne!(rows[0]["status"], "ok");
    assert_eq!(rows[0]["segments_processed"], 1);
    assert_eq!(r["result"]["segments_processed"], 1);
    no_geometry(&rows[0]);
    assert_eq!(rows[1]["status"], "work_limit");
    no_geometry(&rows[1]);
    req = request();
    req["parameters"][1]["q"][0] = json!(3.0);
    let (r, exit) = run(&s, &req);
    assert_eq!(exit, 2);
    assert_eq!(
        r["result"]["slots"][5]["identity"]["parameters"],
        req["parameters"][1]
    );
    no_geometry(&r["result"]["slots"][5]);
    let original = r["receipt"]["scientific_specification_digest"].clone();
    req["parameters"][1]["q"][0] = json!(-0.4);
    let (r, exit) = run(&s, &req);
    assert_eq!(exit, 6);
    assert_ne!(r["receipt"]["scientific_specification_digest"], original);
    object(&s, original.as_str().unwrap());
}
#[test]
fn strict_schema_caps_empty_batches_and_json_product_limit() {
    let s = scratch();
    for field in [
        "maximum_parameters",
        "maximum_queries",
        "maximum_slots",
        "maximum_native_output_bytes",
    ] {
        let mut req = request();
        req["policy"][field] = json!(0);
        let (r, exit) = run(&s, &req);
        assert_eq!(exit, 2, "{field}: {r}");
        assert_eq!(r["result"]["batch_status"], "work_limit");
        assert_eq!(r["result"]["numerical_status"], "work_limit");
        assert_eq!(r["result"]["slots"], json!([]));
    }
    let mut req = request();
    req["parameters"] = json!([]);
    let (r, exit) = run(&s, &req);
    assert_eq!(exit, 6, "{r}");
    assert_eq!(r["result"]["slots"], json!([]));
    for mode in 0..5 {
        let mut req = request();
        match mode {
            0 => req["parameters"][0]["q"] = json!([0.0, 0.0, 0.0, 0.0]),
            1 => req["parameters"][0]["arithmetic"] = json!("implicit"),
            2 => req["queries"][0]["convention"] = json!("unknown"),
            3 => {
                req["policy"]
                    .as_object_mut()
                    .unwrap()
                    .remove("maximum_total_segment_visits");
            }
            _ => req["policy"]["maximum_slots"] = json!(1048577),
        }
        assert_eq!(run(&s, &req).1, 2);
    }
    let mut invalid_convention = request();
    invalid_convention["queries"][0]["convention"] = json!("unknown");
    let (r, exit) = run(&s, &invalid_convention);
    assert_eq!(exit, 2);
    assert_eq!(r["result"]["error_id"], "UNKNOWN_BACKGROUND_CONVENTION");
    assert!(r["result"].get("slots").is_none());
    let mut req = request();
    req["parameters"] = json!(vec![req["parameters"][0].clone(); 257]);
    req["queries"] = json!(vec![query(0.1); 256]);
    let (r, exit) = run(&s, &req);
    assert_eq!(exit, 2);
    assert_eq!(r["result"]["error_id"], "RESOURCE_LIMIT");
}
