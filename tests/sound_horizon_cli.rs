//! Boundary/records controls; independent scientific equations live in native tests.
use serde_json::{json, Value};
use sha2::Digest;
use std::{
    fs,
    process::Command,
    time::{SystemTime, UNIX_EPOCH},
};
fn request() -> Value {
    json!({"schema_version":2,"operation":"cosmology.sound_horizon","points":[{"h0_km_s_mpc":70.0,"omega_m":0.0,"omega_r":1.0,"omega_b":0.0,"omega_gamma":1.0,"z_drag":999.0,"drag_origin":"analytic radiation-only control"}],"numerical_policy":{"absolute_tolerance_mpc":1e-9,"relative_tolerance":2e-11,"maximum_callbacks_per_point":1000,"maximum_depth":40,"maximum_points":8,"maximum_total_callbacks":1000,"maximum_native_bytes":1048576}})
}
fn run(r: &Value, qualified: bool) -> (Value, i32) {
    let dir = std::env::temp_dir().join(format!(
        "irred-sound-{}-{}",
        std::process::id(),
        SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    fs::create_dir(&dir).unwrap();
    let path = dir.join("input.json");
    fs::write(&path, serde_json::to_vec(r).unwrap()).unwrap();
    let mut c = Command::new(env!("CARGO_BIN_EXE_irred"));
    c.args([
        "run",
        path.to_str().unwrap(),
        dir.join("store").to_str().unwrap(),
    ]);
    if qualified {
        c.args(["--assurance", "qualified"]);
    }
    let o = c.output().unwrap();
    let v: Value = serde_json::from_slice(&o.stdout)
        .unwrap_or_else(|_| panic!("{}", String::from_utf8_lossy(&o.stderr)));
    let input_bytes = fs::read(&path).unwrap();
    let input_digest = format!("{:x}", sha2::Sha256::digest(&input_bytes));
    assert_eq!(v["receipt"]["input_digest"], input_digest);
    assert_eq!(
        fs::read(dir.join("store/objects").join(&input_digest)).unwrap(),
        input_bytes
    );
    let output_digest = v["receipt"]["output_digest"].as_str().unwrap();
    let output_bytes = fs::read(dir.join("store/objects").join(output_digest)).unwrap();
    assert_eq!(
        format!("{:x}", sha2::Sha256::digest(&output_bytes)),
        output_digest
    );
    assert_eq!(
        serde_json::from_slice::<Value>(&output_bytes).unwrap(),
        v["result"]
    );
    for (i, row) in v["result"]["rows"].as_array().unwrap().iter().enumerate() {
        assert_eq!(row["source"], r["points"][i]);
    }
    assert_eq!(
        v["result"]["metadata"]["model_id"],
        "flat-pressureless-matter-massless-radiation-lambda"
    );
    assert_eq!(
        v["result"]["metadata"]["equation_id"],
        "conditional-tight-coupling-sound-horizon-scale-factor"
    );
    assert_eq!(
        v["receipt"]["precision"],
        v["result"]["metadata"]["arithmetic_id"]
    );
    fs::remove_dir_all(dir).unwrap();
    (v, o.status.code().unwrap())
}
#[test]
fn radiation_control_and_qualification() {
    let (v, c) = run(&request(), false);
    assert_eq!(c, 0);
    let row = &v["result"]["rows"][0];
    let expected = 299792.458 / 70.0 / 3f64.sqrt() / 1000.0;
    assert!((row["sound_horizon_mpc"].as_f64().unwrap() - expected).abs() < 1e-9);
    assert_eq!(row["numerical_status"], "ok");
    assert_eq!(v["receipt"]["accepted"], true);
    let (q, c) = run(&request(), true);
    assert_eq!(c, 6);
    assert_eq!(q["result"], v["result"]);
}
#[test]
fn mixed_invalid_model_and_work_cap_preserve_causes() {
    let mut r = request();
    let mut invalid = r["points"][0].clone();
    invalid["omega_m"] = json!(1.0);
    r["points"].as_array_mut().unwrap().push(invalid);
    let (v, c) = run(&r, false);
    assert_eq!(c, 2);
    assert_eq!(v["result"]["rows"][0]["numerical_status"], "ok");
    assert_eq!(v["result"]["rows"][1]["numerical_status"], "outside_domain");
    assert!(v["result"]["rows"][1]["sound_horizon_mpc"].is_null());
    r["numerical_policy"]["maximum_total_callbacks"] = json!(0);
    let (v, c) = run(&r, false);
    assert_eq!(c, 2);
    assert_eq!(v["result"]["rows"][0]["numerical_status"], "work_limit");
    assert_eq!(v["result"]["callbacks"], 0);
}

#[test]
fn stateless_stream_matches_one_shot_and_rejects_duplicates() {
    use std::{io::Write, process::Stdio};
    let r = request();
    let (one, code) = run(&r, false);
    assert_eq!(code, 0);
    let dir = std::env::temp_dir().join(format!(
        "irred-sound-stream-{}-{}",
        std::process::id(),
        SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    fs::create_dir(&dir).unwrap();
    let mut child = Command::new(env!("CARGO_BIN_EXE_irred"))
        .arg("stream")
        .arg(&dir)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .unwrap();
    let mut input = child.stdin.take().unwrap();
    input.write_all(b"{\"action\":\"sound_horizon_evaluate\",\"action\":\"sound_horizon_evaluate\",\"request\":{}}\n").unwrap();
    writeln!(
        input,
        "{}",
        json!({"action":"sound_horizon_evaluate","request":r})
    )
    .unwrap();
    drop(input);
    let o = child.wait_with_output().unwrap();
    let lines: Vec<Value> = String::from_utf8(o.stdout)
        .unwrap()
        .lines()
        .map(|l| serde_json::from_str(l).unwrap())
        .collect();
    assert_eq!(lines.len(), 2);
    assert_eq!(lines[0]["result"]["execution"], "failed");
    assert_eq!(lines[1]["result"]["output"], one["result"]);
    assert_eq!(
        format!(
            "{:x}",
            sha2::Sha256::digest(serde_json::to_vec(&lines[1]["result"]["specification"]).unwrap())
        ),
        one["receipt"]["scientific_specification_digest"]
    );
    fs::remove_dir_all(dir).unwrap();
}
