//! Transport/immutable records. Existing independent native peers own scientific evidence.
use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use std::{
    fs,
    io::Write,
    process::{Command, Stdio},
};
struct Store(std::path::PathBuf);
impl Store {
    fn new() -> Self {
        let p = std::env::temp_dir().join(format!(
            "irred-sampled-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir(&p).unwrap();
        Self(p)
    }
    fn object(&self, digest: &Value) -> Value {
        let h = digest.as_str().unwrap();
        let b = fs::read(self.0.join("store/objects").join(h)).unwrap();
        assert_eq!(format!("{:x}", Sha256::digest(&b)), h);
        serde_json::from_slice(&b).unwrap()
    }
}
impl Drop for Store {
    fn drop(&mut self) {
        fs::remove_dir_all(&self.0).unwrap();
    }
}
fn request() -> Value {
    serde_json::from_str(include_str!("fixtures/sampled_photometry.json")).unwrap()
}
fn run_raw(s: &Store, b: &[u8], qualified: bool) -> (Value, i32) {
    let p = s.0.join("request.json");
    fs::write(&p, b).unwrap();
    let mut c = Command::new(env!("CARGO_BIN_EXE_irred"));
    c.arg("run").arg(p).arg(s.0.join("store"));
    if qualified {
        c.args(["--assurance", "qualified"]);
    }
    let o = c.output().unwrap();
    let v: Value = serde_json::from_slice(&o.stdout)
        .unwrap_or_else(|_| panic!("{}", String::from_utf8_lossy(&o.stderr)));
    let hash = format!("{:x}", Sha256::digest(b));
    assert_eq!(v["receipt"]["input_digest"], hash);
    assert_eq!(fs::read(s.0.join("store/objects").join(hash)).unwrap(), b);
    assert_eq!(s.object(&v["receipt"]["output_digest"]), v["result"]);
    (v, o.status.code().unwrap())
}
fn run(s: &Store, r: &Value) -> (Value, i32) {
    run_raw(s, &serde_json::to_vec(r).unwrap(), false)
}
#[test]
fn constant_rectangular_parity_and_source_records() {
    let s = Store::new();
    let r = request();
    let (v, c) = run(&s, &r);
    assert_eq!(c, 0);
    assert_eq!(
        v["receipt"]["method"],
        "analytic_piecewise_linear_merged_bernstein"
    );
    assert_eq!(
        v["receipt"]["precision"],
        "binary64_storage_longdouble_intermediate"
    );
    let spec = s.object(&v["receipt"]["scientific_specification_digest"]);
    for field in ["spectra", "passbands", "exposures", "resource_policy"] {
        assert_eq!(spec[field], r[field]);
    }
    let rect: Value = serde_json::from_str(include_str!("fixtures/photometry.json")).unwrap();
    let (a, c) = run(&s, &rect);
    assert_eq!(c, 0);
    for id in [
        "incident_band_flux",
        "collected_energy",
        "expected_transmitted_photons",
    ] {
        let x = v["result"]["evaluations"][0]["outputs"][id]["value"]
            .as_f64()
            .unwrap();
        let y = a["result"]["evaluations"][0]["outputs"][id]["value"]
            .as_f64()
            .unwrap();
        assert!((x - y).abs() <= 2e-12 * x.abs() + 1e-300);
    }
    let (q, c) = run_raw(&s, &serde_json::to_vec(&r).unwrap(), true);
    assert_eq!(c, 6);
    assert_eq!(q["result"], v["result"]);
    assert_eq!(q["receipt"]["accepted"], false);
    let mut measured = r.clone();
    measured["spectra"][0]["source_role"] = json!("measured");
    measured["spectra"][0]["provenance"] =
        json!("caller-declared measured source, no release verification");
    measured["passbands"][0]["source_role"] = json!("calibration_asset");
    measured["passbands"][0]["calibration"] = json!("declared_uncertainty_excluded");
    let (m, c) = run(&s, &measured);
    assert_eq!(c, 0);
    let spec = s.object(&m["receipt"]["scientific_specification_digest"]);
    assert_eq!(spec["spectra"], measured["spectra"]);
    assert_eq!(spec["passbands"], measured["passbands"]);
    assert_ne!(
        m["receipt"]["scientific_specification_digest"],
        v["receipt"]["scientific_specification_digest"]
    );
    assert_eq!(m["result"], v["result"]);
}
#[test]
fn mixed_rows_order_omission_and_quotas() {
    let s = Store::new();
    let mut r = request();
    let mut e = r["exposures"][0].clone();
    e["redshift"] = json!(-1);
    r["exposures"].as_array_mut().unwrap().push(e);
    let mut e = r["exposures"][0].clone();
    e["spectrum_index"] = json!(99);
    r["exposures"].as_array_mut().unwrap().push(e);
    r["resource_policy"]["maximum_rows"] = json!(3);
    let (v, c) = run(&s, &r);
    assert_eq!(c, 2);
    let rows = v["result"]["evaluations"].as_array().unwrap();
    assert_eq!(rows.len(), 3);
    for (i, row) in rows.iter().enumerate() {
        for field in ["spectrum_index", "passband_index"] {
            assert_eq!(row["source"][field], r["exposures"][i][field]);
        }
        for field in [
            "luminosity_distance_metre",
            "redshift",
            "collecting_area_square_metre",
            "observer_exposure_second",
        ] {
            assert_eq!(
                row["source"][field].as_f64().unwrap().to_bits(),
                r["exposures"][i][field].as_f64().unwrap().to_bits()
            );
        }
        assert_eq!(row["index"], i);
    }
    assert_eq!(rows[0]["admission_status"], "ok");
    assert_eq!(rows[1]["admission_status"], "outside_domain");
    assert_eq!(rows[2]["admission_status"], "invalid_input");
    assert!(rows[2]["spectrum_id"].is_null());
    for row in &rows[1..] {
        assert_eq!(row["outputs"]["collected_energy"]["availability"], "failed");
        assert!(row["outputs"]["collected_energy"]["value"].is_null());
    }
    r = request();
    r["requested_outputs"] = json!(["collected_energy"]);
    r["exposures"][0]["collecting_area_square_metre"] = json!(0);
    r["exposures"][0]["luminosity_distance_metre"] = json!(1e308);
    let (v, c) = run(&s, &r);
    assert_eq!(c, 0);
    assert_eq!(
        v["result"]["evaluations"][0]["outputs"]
            .as_object()
            .unwrap()
            .len(),
        1
    );
    assert_eq!(
        v["result"]["evaluations"][0]["outputs"]["collected_energy"]["value"],
        0.0
    );
    for (key, n) in [
        ("maximum_rows", 0),
        ("maximum_total_samples", 3),
        ("maximum_native_bytes", 0),
    ] {
        let mut r = request();
        r["resource_policy"][key] = json!(n);
        let (v, c) = run(&s, &r);
        assert_eq!(c, 2);
        assert_eq!(v["result"]["batch_numerical_status"], "work_limit");
        assert_eq!(v["result"]["evaluations"], json!([]));
    }
    let mut r = request();
    r["resource_policy"]["maximum_segments"] = json!(0);
    let (v, c) = run(&s, &r);
    assert_eq!(c, 2);
    assert_eq!(v["result"]["batch_numerical_status"], "ok");
    assert_eq!(
        v["result"]["evaluations"][0]["admission_status"],
        "work_limit"
    );
    assert_eq!(
        v["result"]["evaluations"][0]["outputs"]["collected_energy"]["availability"],
        "failed"
    );
}
#[test]
fn raw_duplicate_and_structural_rejections() {
    let s = Store::new();
    for case in 0..5 {
        let mut r = request();
        match case {
            0 => {
                let c = r["spectra"][0].clone();
                r["spectra"].as_array_mut().unwrap().push(c);
            }
            1 => r["spectra"][0]["luminosity_watt_per_metre"] = json!([1]),
            2 => r["passbands"][0]["calibration"] = json!("inferred"),
            3 => r["resource_policy"]["maximum_samples"] = json!(65537),
            _ => r["spectra"][0]["unexpected"] = json!(1),
        }
        let (v, c) = run(&s, &r);
        assert_eq!(c, 2);
        assert_eq!(v["receipt"]["execution"], "failed");
        assert!(v["receipt"]["scientific_specification_digest"].is_null());
    }
    let text = serde_json::to_string(&request()).unwrap().replace(
        "\"spectrum_index\":0",
        "\"spectrum_index\":0,\"spectrum_index\":99",
    );
    let (v, c) = run_raw(&s, text.as_bytes(), false);
    assert_eq!(c, 2);
    assert_eq!(v["receipt"]["execution"], "failed");
    assert!(v["receipt"]["scientific_specification_digest"].is_null());
}
#[test]
fn linear_polynomial_fact_and_refinement() {
    let s = Store::new();
    let mut r = request();
    r["spectra"][0]["rest_wavelength_metre"] = json!([1e-6, 2e-6]);
    r["spectra"][0]["luminosity_watt_per_metre"] = json!([2, 4]);
    r["passbands"][0]["optical_transmission"] = json!([0, 0.5]);
    let (v, c) = run(&s, &r);
    assert_eq!(c, 0);
    // Independent polynomial antiderivatives with high-precision decimal pi/h/c:
    // integral L =3e-6; integral L*T =5e-7*5/3;
    // integral lambda*L*T =5e-13*17/6; multiply A*t=6 then 1/(4pi), photon 1/hc.
    for (id, expected) in [
        ("incident_band_flux", 2.3873241463784300365e-7),
        ("collected_energy", 3.9788735772973833942e-7),
        ("expected_transmitted_photons", 3405119277257.2164),
    ] {
        let x = v["result"]["evaluations"][0]["outputs"][id]["value"]
            .as_f64()
            .unwrap();
        assert!(
            x > 0. && (x - expected).abs() <= 2e-12 * expected + 1e-300,
            "{id} actual {x} expected {expected}"
        );
    }
    r["spectra"][0]["rest_wavelength_metre"] = json!([1e-6, 1.5e-6, 2e-6]);
    r["spectra"][0]["luminosity_watt_per_metre"] = json!([2, 3, 4]);
    r["resource_policy"]["maximum_total_samples"] = json!(5);
    r["resource_policy"]["maximum_samples"] = json!(5);
    r["resource_policy"]["maximum_segments"] = json!(2);
    let (w, c) = run(&s, &r);
    assert_eq!(c, 0);
    for id in [
        "incident_band_flux",
        "collected_energy",
        "expected_transmitted_photons",
    ] {
        let a = v["result"]["evaluations"][0]["outputs"][id]["value"]
            .as_f64()
            .unwrap();
        let b = w["result"]["evaluations"][0]["outputs"][id]["value"]
            .as_f64()
            .unwrap();
        assert!((a - b).abs() <= 2e-12 * a.abs() + 1e-300);
    }
}
#[test]
fn stateless_stream_parity_and_raw_frame_receipts() {
    let s = Store::new();
    let r = request();
    let (one, c) = run(&s, &r);
    assert_eq!(c, 0);
    let command = json!({"action":"photometry_predict","request":r});
    let mut raw = serde_json::to_vec(&command).unwrap();
    raw.push(b'\n');
    let mut child = Command::new(env!("CARGO_BIN_EXE_irred"))
        .arg("stream")
        .arg(s.0.join("stream-store"))
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .unwrap();
    child.stdin.take().unwrap().write_all(&raw).unwrap();
    let o = child.wait_with_output().unwrap();
    assert!(o.status.success());
    let v: Value = serde_json::from_slice(&o.stdout).unwrap();
    assert_eq!(v["result"]["output"], one["result"]);
    assert_eq!(v["result"]["method"], one["receipt"]["method"]);
    assert_eq!(v["result"]["accepted"], true);
    assert_eq!(
        v["receipt"]["input_digest"],
        format!("{:x}", Sha256::digest(&raw))
    );
    let digest = v["receipt"]["input_digest"].as_str().unwrap();
    assert_eq!(
        fs::read(s.0.join("stream-store/objects").join(digest)).unwrap(),
        raw
    );
    let spec = s.object(&one["receipt"]["scientific_specification_digest"]);
    assert_eq!(v["result"]["specification"], spec);
}

#[test]
fn partial_group_failure_and_signed_zero_bits() {
    let s = Store::new();
    let mut r = request();
    r["exposures"][0]["collecting_area_square_metre"] = json!(-0.0);
    r["exposures"][0]["luminosity_distance_metre"] = json!(1e308);
    let (v, c) = run(&s, &r);
    assert_eq!(c, 2);
    let row = &v["result"]["evaluations"][0];
    assert_eq!(
        row["source"]["collecting_area_square_metre"]
            .as_f64()
            .unwrap()
            .to_bits(),
        (-0.0f64).to_bits()
    );
    assert_eq!(row["admission_status"], "ok");
    assert_eq!(
        row["outputs"]["incident_band_flux"]["availability"],
        "failed"
    );
    assert_eq!(
        row["outputs"]["collected_energy"]["availability"],
        "available"
    );
    assert_eq!(
        row["outputs"]["expected_transmitted_photons"]["availability"],
        "available"
    );
    let checks = v["receipt"]["outputs"].as_array().unwrap();
    assert_eq!(checks[0]["numerical"], "failed");
    assert_eq!(checks[1]["numerical"], "checks_passed");
    assert_eq!(checks[2]["numerical"], "checks_passed");
}
