//! Current typed projection/assurance boundary, not independent science qualification.
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
        "irred-current-bg-{}-{}",
        std::process::id(),
        SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    fs::create_dir(&p).unwrap();
    Scratch(p)
}
fn digest(b: &[u8]) -> String {
    format!("{:x}", Sha256::digest(b))
}
fn request() -> Value {
    json!({"schema_version":2,"operation":"background.evaluate","geometry":{"kind":"flat_flrw"},"models":[{"kind":"constant_q","q":0.0}],"queries":[{"z_expansion":1.0,"requested_outputs":["expansion"]},{"z_expansion":1.0,"requested_outputs":["expansion"]}],"numerical_policy":{"maximum_models":4,"maximum_queries":16,"maximum_slots":64,"maximum_callbacks":0,"maximum_segment_visits":0,"maximum_native_bytes":1000000}})
}
fn run(s: &Scratch, v: &Value, qualified: bool) -> (Value, i32) {
    let bytes = serde_json::to_vec(v).unwrap();
    let p = s.0.join("request.json");
    fs::write(&p, &bytes).unwrap();
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_irred"));
    cmd.args([
        "run",
        p.to_str().unwrap(),
        s.0.join("store").to_str().unwrap(),
    ]);
    if qualified {
        cmd.args(["--assurance", "qualified"]);
    }
    let o = cmd.output().unwrap();
    let value: Value = serde_json::from_slice(&o.stdout)
        .unwrap_or_else(|_| panic!("{}", String::from_utf8_lossy(&o.stderr)));
    assert_eq!(value["receipt"]["input_digest"], digest(&bytes));
    assert_eq!(
        fs::read(s.0.join("store/objects").join(digest(&bytes))).unwrap(),
        bytes
    );
    let d = value["receipt"]["output_digest"].as_str().unwrap();
    let out = fs::read(s.0.join("store/objects").join(d)).unwrap();
    assert_eq!(digest(&out), d);
    assert_eq!(
        serde_json::from_slice::<Value>(&out).unwrap(),
        value["result"]
    );
    (value, o.status.code().unwrap())
}
#[test]
fn expansion_only_exact_identity_and_assurance() {
    let s = scratch();
    let (v, code) = run(&s, &request(), false);
    assert_eq!(code, 0);
    assert_eq!(v["receipt"]["accepted"], true);
    let r = &v["result"];
    assert_eq!(r["work"]["callbacks"], 0);
    assert_eq!(r["work"]["segment_visits"], 0);
    assert_eq!(r["nodes"].as_array().unwrap().len(), 1);
    for row in r["evaluations"].as_array().unwrap() {
        assert_eq!(row["groups"].as_object().unwrap().len(), 1);
        assert_eq!(
            row["groups"]["expansion"]["value"]["E"]
                .as_f64()
                .unwrap()
                .to_bits(),
            2.0f64.to_bits()
        );
    }
    let (q, code) = run(&s, &request(), true);
    assert_eq!(code, 6);
    assert_eq!(q["result"], v["result"]);
    assert_eq!(q["receipt"]["accepted"], false);
}
#[test]
fn nested_dimensional_failures_preserve_dimensionless_outputs() {
    let s = scratch();
    let mut v = request();
    v["queries"] = json!([{"z_expansion":0.5,"requested_outputs":["expansion","clock"],"physical_scale":{"h0_km_s_mpc":0.0}}]);
    v["numerical_policy"]["integration"] = json!({"absolute_tolerance":1e-12,"relative_tolerance":1e-12,"maximum_evaluations":100000,"maximum_depth":24});
    v["numerical_policy"]["maximum_callbacks"] = json!(100000);
    let (r, code) = run(&s, &v, false);
    assert_eq!(code, 2);
    let g = &r["result"]["evaluations"][0]["groups"];
    assert_eq!(g["expansion"]["availability"], "available");
    assert_eq!(
        g["expansion"]["value"]["H_km_s_mpc"]["availability"],
        "failed"
    );
    assert_eq!(g["clock"]["availability"], "available");
    assert_eq!(
        g["clock"]["value"]["lookback_seconds"]["availability"],
        "failed"
    );
    assert!(g["clock"]["value"]["lookback_seconds"]["value"].is_null());
    assert_eq!(r["receipt"]["execution"], "completed");
    assert_eq!(r["receipt"]["accepted"], false);
}
#[test]
fn invalid_model_is_owned_empty_not_false_success_and_schema_is_strict() {
    let s = scratch();
    let mut v = request();
    v["models"] = json!([{"kind":"constant_q","q":0.0},{"kind":"constant_q","q":3.0}]);
    let (r, code) = run(&s, &v, false);
    assert_eq!(code, 2);
    assert_eq!(r["result"]["models"][1]["row_count"], 0);
    assert_eq!(
        r["result"]["models"][1]["source"]["q"]
            .as_f64()
            .unwrap()
            .to_bits(),
        3.0f64.to_bits()
    );
    assert_eq!(r["result"]["evaluations"].as_array().unwrap().len(), 2);
    assert_eq!(r["receipt"]["execution"], "completed");
    v["models"][0]["omega_m"] = json!(0.3);
    let (r, code) = run(&s, &v, false);
    assert_eq!(code, 2);
    assert_eq!(r["receipt"]["execution"], "failed");
}

#[test]
fn final_units_survive_unrequested_metre_overflow() {
    let s = scratch();
    let mut v = request();
    v["models"] = json!([{"kind":"constant_q","q":-1.0}]);
    v["numerical_policy"]["integration"] = json!({"absolute_tolerance":1e-14,"relative_tolerance":1e-12,"maximum_evaluations":10000,"maximum_depth":30});
    v["numerical_policy"]["maximum_callbacks"] = json!(10000);
    for (z, outputs) in [(0.5, json!(["clock"])), (1e-300, json!(["physical"]))] {
        v["queries"] = json!([{"z_expansion":z,"requested_outputs":outputs,"observer":{"redshift":z,"convention":"geometric_same_redshift"},"physical_scale":{"h0_km_s_mpc":1e-285}}]);
        let (r, code) = run(&s, &v, false);
        assert_eq!(code, 0);
        assert_eq!(r["receipt"]["accepted"], true);
        let groups = &r["result"]["evaluations"][0]["groups"];
        if z == 0.5 {
            assert_eq!(
                groups["clock"]["value"]["lookback_seconds"]["availability"],
                "available"
            );
            let seconds = groups["clock"]["value"]["lookback_seconds"]["value"]
                .as_f64()
                .unwrap();
            // Independent numeric fixture from the separately reviewed unit
            // definition and exact de Sitter clock log(1+z), budget unchanged.
            assert!((seconds / 1.2511345941663363e304 - 1.0).abs() < 1e-13);
        } else {
            assert_eq!(groups["physical"]["availability"], "available");
            for key in [
                "radial_mpc",
                "transverse_mpc",
                "angular_diameter_mpc",
                "luminosity_mpc",
            ] {
                let value = groups["physical"]["value"][key].as_f64().unwrap();
                assert!((value / 2.99792458e-10 - 1.0).abs() < 1e-13);
            }
            let volume = groups["physical"]["value"]["volume_mpc3_per_sr_per_redshift"]
                .as_f64()
                .unwrap();
            assert!((volume / 2.69440024173739849e271 - 1.0).abs() < 1e-13);
        }
    }
    v["queries"] = json!([{"z_expansion":0.5,"requested_outputs":["physical","expansion"],"observer":{"redshift":0.5,"convention":"geometric_same_redshift"},"physical_scale":{"h0_km_s_mpc":1e-285}}]);
    let (r, code) = run(&s, &v, false);
    assert_eq!(code, 2);
    let groups = &r["result"]["evaluations"][0]["groups"];
    assert_eq!(groups["physical"]["availability"], "failed");
    assert!(groups["physical"]["value"].is_null());
    assert_eq!(groups["expansion"]["availability"], "available");
}
