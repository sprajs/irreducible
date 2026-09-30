//! Exact IEEE binary64 request identity, independent of tolerant physics checks.
//! Regression: default serde_json decoder rounded -15.803224250000001 one ULP
//! toward zero. Source ASCII std::parse and native equations are separate.
use serde_json::Value;
use std::{
    fs,
    process::Command,
    time::{SystemTime, UNIX_EPOCH},
};
fn values() -> Vec<f64> {
    let mut v = vec!["-15.803224250000001".parse::<f64>().unwrap()];
    // Selection threshold, piecewise endpoint and q admission endpoints.
    for base in [0.01_f64, 2.5, -3.0, 2.0] {
        v.extend([
            f64::from_bits(base.to_bits() - 1),
            base,
            f64::from_bits(base.to_bits() + 1),
        ]);
    }
    v
}
#[test]
fn decimal_roundtrip_preserves_bits_at_scientific_boundaries() {
    for expected in values() {
        let text = expected.to_string();
        assert_eq!(text.parse::<f64>().unwrap().to_bits(), expected.to_bits());
        let parsed: f64 = serde_json::from_str(&text).unwrap();
        assert_eq!(parsed.to_bits(), expected.to_bits(), "lexeme {text}");
        // A parser must not change which side of a named domain/selection edge
        // this source lies on, even when downstream numerical tolerance passes.
        for edge in [0.01, 2.5, -3.0, 2.0] {
            assert_eq!(parsed.partial_cmp(&edge), expected.partial_cmp(&edge));
        }
    }
}
#[test]
fn actual_scalar_request_spec_source_and_native_sum_preserve_decimal() {
    let path = std::env::temp_dir().join(format!(
        "irred-json-precision-{}-{}",
        std::process::id(),
        SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    fs::create_dir(&path).unwrap();
    struct Scratch(std::path::PathBuf);
    impl Drop for Scratch {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.0);
        }
    }
    let path = Scratch(path);
    for (i, expected) in values().into_iter().enumerate() {
        let text = expected.to_string();
        let raw = format!("{{\"schema_version\":1,\"operation\":\"numerics.scalar_batch.v1\",\"method\":\"compensated_sum\",\"values\":[{text}]}}");
        let request = path.0.join("request.json");
        fs::write(&request, &raw).unwrap();
        let store = path.0.join(format!("store-{i}"));
        // Use the real noninteractive run entry point; no alternate test parser.
        let output = Command::new(env!("CARGO_BIN_EXE_irred"))
            .arg("run")
            .arg(&request)
            .arg(&store)
            .output()
            .unwrap();
        assert_eq!(
            output.status.code(),
            Some(6),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        let result: Value = serde_json::from_slice(&output.stdout).unwrap();
        let hash = result["receipt"]["scientific_specification_digest"]
            .as_str()
            .unwrap();
        let spec: Value =
            serde_json::from_slice(&fs::read(store.join("objects").join(hash)).unwrap()).unwrap();
        assert_eq!(
            spec["values"][0].as_f64().unwrap().to_bits(),
            expected.to_bits(),
            "spec {text}"
        );
        assert_eq!(
            result["result"]["source_values"][0]
                .as_f64()
                .unwrap()
                .to_bits(),
            expected.to_bits(),
            "source {text}"
        );
        assert_eq!(
            result["result"]["evaluations"][0]["result"]["value"]
                .as_f64()
                .unwrap()
                .to_bits(),
            expected.to_bits(),
            "native sum {text}"
        );
    }
}
