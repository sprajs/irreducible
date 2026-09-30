//! Handwritten raw bytes exercise parsing before JSON Value loses duplicates.
use serde_json::Value;
use sha2::{Digest, Sha256};
use std::{
    fs,
    io::Write,
    process::{Command, Stdio},
    time::{SystemTime, UNIX_EPOCH},
};
const VALID: &str = r#"{"schema_version":2,"operation":"background.evaluate","geometry":{"kind":"flat_flrw"},"models":[{"kind":"constant_q","q":0}],"queries":[{"z_expansion":1,"requested_outputs":["expansion"]}],"numerical_policy":{"maximum_models":1,"maximum_queries":1,"maximum_slots":1,"maximum_callbacks":0,"maximum_segment_visits":0,"maximum_native_bytes":1000000}}"#;
#[test]
fn raw_duplicates_reject_both_routes_and_stream_recovers() {
    let dir = std::env::temp_dir().join(format!(
        "irred-strict-json-{}-{}",
        std::process::id(),
        SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    fs::create_dir_all(&dir).unwrap();
    let duplicates = [
        VALID.replace("\"q\":0", "\"q\":3,\"q\":0"),
        VALID.replace(
            "\"maximum_callbacks\":0",
            "\"maximum_callbacks\":1,\"maximum_callbacks\":0",
        ),
        VALID.replace(
            "\"geometry\":",
            "\"irrelevant\":{\"x\":1,\"x\":2},\"geometry\":",
        ),
        VALID.replace(
            "\"geometry\":",
            "\"source\":{\"metadata\":{\"role\":1,\"role\":2}},\"geometry\":",
        ),
    ];
    let mut frames = String::new();
    for (i, raw) in duplicates.iter().enumerate() {
        let request = dir.join(format!("request-{i}.json"));
        fs::write(&request, raw).unwrap();
        let output = Command::new(env!("CARGO_BIN_EXE_irred"))
            .args([
                "run",
                request.to_str().unwrap(),
                dir.join(format!("one-{i}")).to_str().unwrap(),
            ])
            .output()
            .unwrap();
        assert!(!output.status.success());
        let v: Value = serde_json::from_slice(&output.stdout).unwrap();
        assert_eq!(v["receipt"]["execution"], "failed");
        let digest = format!("{:x}", Sha256::digest(raw.as_bytes()));
        assert_eq!(v["receipt"]["input_digest"], digest);
        assert!(v["receipt"]["scientific_specification_digest"].is_null());
        assert_eq!(
            fs::read(dir.join(format!("one-{i}/objects/{digest}"))).unwrap(),
            raw.as_bytes()
        );
        frames.push_str(&format!(
            "{{\"action\":\"background_evaluate\",\"request\":{raw}}}\n"
        ));
    }
    frames
        .push_str("{\"action\":\"release\",\"action\":\"background_evaluate\",\"request\":{} }\n");
    frames.push_str(&format!(
        "{{\"action\":\"background_evaluate\",\"request\":{VALID}}}\n"
    ));
    let mut child = Command::new(env!("CARGO_BIN_EXE_irred"))
        .args(["stream", dir.join("stream").to_str().unwrap()])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .unwrap();
    child
        .stdin
        .take()
        .unwrap()
        .write_all(frames.as_bytes())
        .unwrap();
    let output = child.wait_with_output().unwrap();
    assert!(output.status.success());
    let rows: Vec<Value> = String::from_utf8(output.stdout)
        .unwrap()
        .lines()
        .map(|s| serde_json::from_str(s).unwrap())
        .collect();
    assert_eq!(rows.len(), 6);
    for row in &rows[..5] {
        assert_eq!(row["result"]["execution"], "failed");
        assert_eq!(row["result"]["error_id"], "INVALID_SESSION_COMMAND");
    }
    assert_eq!(rows[5]["result"]["execution"], "completed");
    fs::remove_dir_all(dir).unwrap();
}
