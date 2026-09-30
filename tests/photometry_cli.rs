//! Boundary/record controls; native independent frequency oracle owns science evidence.
use serde_json::Value;
use sha2::{Digest, Sha256};
use std::{fs, process::Command};
fn run(v: &Value, qualified: bool) -> (Value, i32) {
    let dir = std::env::temp_dir().join(format!(
        "irred-photo-cli-{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    fs::create_dir(&dir).unwrap();
    let bytes = serde_json::to_vec(v).unwrap();
    let request = dir.join("request.json");
    fs::write(&request, &bytes).unwrap();
    let store = dir.join("store");
    let mut c = Command::new(env!("CARGO_BIN_EXE_irred"));
    c.arg("run").arg(request).arg(&store);
    if qualified {
        c.args(["--assurance", "qualified"]);
    }
    let o = c.output().unwrap();
    let out: Value = serde_json::from_slice(&o.stdout)
        .unwrap_or_else(|_| panic!("{}", String::from_utf8_lossy(&o.stderr)));
    let digest = format!("{:x}", Sha256::digest(&bytes));
    assert_eq!(out["receipt"]["input_digest"], digest);
    assert_eq!(fs::read(store.join("objects").join(digest)).unwrap(), bytes);
    let hash = out["receipt"]["output_digest"].as_str().unwrap();
    let recorded = fs::read(store.join("objects").join(hash)).unwrap();
    assert_eq!(format!("{:x}", Sha256::digest(&recorded)), hash);
    assert_eq!(
        serde_json::from_slice::<Value>(&recorded).unwrap(),
        out["result"]
    );
    fs::remove_dir_all(dir).unwrap();
    (out, o.status.code().unwrap())
}
fn request() -> Value {
    serde_json::from_str(include_str!("fixtures/photometry.json")).unwrap()
}
#[test]
fn analytic_si_outputs_and_assurance() {
    let r = request();
    let (v, code) = run(&r, false);
    assert_eq!(code, 0);
    assert_eq!(v["receipt"]["accepted"], true);
    let row = &v["result"]["evaluations"][0];
    assert_eq!(row["source"], r["inputs"][0]);
    let outputs = &row["outputs"];
    let flux = r["inputs"][0]["luminosity_watt_per_metre"]
        .as_f64()
        .unwrap()
        / (4.0 * std::f64::consts::PI)
        * 1e-6;
    let expected = [
        ("incident_band_flux", flux),
        ("collected_energy", 3.0 * flux),
        (
            "expected_transmitted_photons",
            3.0 * flux * 1.5e-6 / (6.62607015e-34 * 299792458.0),
        ),
    ];
    for (id, x) in expected {
        assert_eq!(outputs[id]["availability"], "available");
        let y = outputs[id]["value"].as_f64().unwrap();
        assert!((y - x).abs() <= 2e-12 * x.abs() + 1e-300);
    }
    let (q, code) = run(&r, true);
    assert_eq!(code, 6);
    assert_eq!(q["receipt"]["accepted"], false);
    assert_eq!(q["result"], v["result"]);
}
#[test]
fn independent_groups_mixed_failure_and_quotas() {
    let mut r = request();
    let mut bad = r["inputs"][0].clone();
    bad["redshift"] = (-1.0).into();
    r["inputs"].as_array_mut().unwrap().push(bad);
    r["resource_policy"]["maximum_rows"] = 2.into();
    let (v, c) = run(&r, false);
    assert_eq!(c, 2);
    assert_eq!(v["result"]["evaluations"].as_array().unwrap().len(), 2);
    assert_eq!(
        v["result"]["evaluations"][0]["outputs"]["collected_energy"]["availability"],
        "available"
    );
    assert_eq!(
        v["result"]["evaluations"][1]["outputs"]["collected_energy"]["availability"],
        "failed"
    );
    assert!(v["result"]["evaluations"][1]["outputs"]["collected_energy"]["value"].is_null());
    let mut r = request();
    r["requested_outputs"] = serde_json::json!(["collected_energy"]);
    r["inputs"][0]["luminosity_distance_metre"] = 1e308.into();
    r["inputs"][0]["collecting_area_square_metre"] = 0.into();
    let (v, c) = run(&r, false);
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
    r["resource_policy"]["maximum_native_bytes"] = 0.into();
    let (v, c) = run(&r, false);
    assert_eq!(c, 2);
    assert_eq!(v["result"]["batch_numerical_status"], "work_limit");
    assert!(v["result"]["evaluations"].as_array().unwrap().is_empty());
}
