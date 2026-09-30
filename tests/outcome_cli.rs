//! Recorder/assurance regression controls. These are interface assertions,
//! not independent physics qualification or new numerical reference values.
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
        "irred-outcome-{}-{}",
        std::process::id(),
        SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    fs::create_dir(&p).unwrap();
    fs::write(p.join("table.dat"),"CID IDSURVEY zHD zCMB zHEL m_b_corr EXTRA\nA 51 .009 .009 .009 17 excluded\nA 52 .1 .1 .1 2 generated\nB 51 .2 .2 .2 -3 generated\n").unwrap();
    fs::write(p.join("cov.dat"), "3\n-5 8 7 6 4 1 0 1 9\n").unwrap();
    Scratch(p)
}
fn fixture(name: &str) -> Value {
    let mut value: Value =
        serde_json::from_slice(&fs::read(format!("tests/fixtures/{name}.json")).unwrap()).unwrap();
    value["schema_version"] = json!(2);
    value
}
fn digest(b: &[u8]) -> String {
    format!("{:x}", Sha256::digest(b))
}
fn object(s: &Scratch, d: &Value) -> Value {
    let d = d.as_str().unwrap();
    let bytes = fs::read(s.0.join("store/objects").join(d)).unwrap();
    assert_eq!(digest(&bytes), d);
    serde_json::from_slice(&bytes).unwrap()
}
fn run(s: &Scratch, v: &Value, assurance: Option<&str>, panic: bool) -> (Value, i32) {
    let bytes = serde_json::to_vec(v).unwrap();
    let path = s.0.join("request.json");
    fs::write(&path, &bytes).unwrap();
    let mut c = Command::new(env!("CARGO_BIN_EXE_irred"));
    c.args([
        "run",
        path.to_str().unwrap(),
        s.0.join("store").to_str().unwrap(),
    ]);
    if let Some(a) = assurance {
        c.args(["--assurance", a]);
    }
    if panic {
        c.env("IRRED_TEST_PANIC", "1");
    }
    let o = c.output().unwrap();
    let r: Value = serde_json::from_slice(&o.stdout)
        .unwrap_or_else(|_| panic!("{}", String::from_utf8_lossy(&o.stderr)));
    assert_eq!(r["receipt"]["input_digest"], digest(&bytes));
    assert_eq!(
        fs::read(s.0.join("store/objects").join(digest(&bytes))).unwrap(),
        bytes
    );
    assert_eq!(object(s, &r["receipt"]["output_digest"]), r["result"]);
    assert_eq!(
        r["receipt"]["executable_digest"],
        digest(&fs::read(env!("CARGO_BIN_EXE_irred")).unwrap())
    );
    (r, o.status.code().unwrap())
}
fn gaussian() -> Value {
    let mut r = fixture("gaussian-density");
    r["schema_version"] = json!(2);
    r["operation"] = json!("statistics.gaussian");
    r["selection"] = json!("all");
    r["maximum_preparation_bytes"] = json!(16777216);
    r["maximum_evaluation_bytes"] = json!(16777216);
    let o = &mut r["observations"];
    o["schema_version"] = json!(2);
    o["operation"] = json!("observations.prepare");
    o.as_object_mut().unwrap().remove("selection");
    o["exports"] = json!([]);
    o["calibration_provenance"] = json!("");
    o["dependence_provenance"] = json!("");
    o["quality_dictionary"] = json!("");
    o["source_selection"] = json!([1, 1]);
    o["resources"]["maximum_preparation_bytes"] = json!(16777216);
    r
}
fn background() -> Value {
    json!({"schema_version":2,"operation":"background.evaluate","geometry":{"kind":"flat_flrw"},"models":[{"kind":"constant_q","q":0.0}],"queries":[{"z_expansion":1.0,"requested_outputs":["expansion"]}],"numerical_policy":{"maximum_models":4,"maximum_queries":16,"maximum_slots":64,"maximum_callbacks":0,"maximum_segment_visits":0,"maximum_native_bytes":1000000}})
}
fn cases(_s: &Scratch) -> Vec<(Value, &'static str, &'static str, &'static str)> {
    let g = gaussian();
    vec![
        (
            fixture("quantity-length"),
            "converted",
            "typed_physical_conversion",
            "binary64_storage_longdouble_intermediate",
        ),
        (
            fixture("numerical-sum"),
            "evaluations",
            "compensated_sum",
            "binary64_storage_method_declared_intermediate",
        ),
        (
            g["observations"].clone(),
            "prepared_observations",
            "immutable_typed_source_preparation",
            "binary64_source_storage",
        ),
        (
            g,
            "gaussian_batch",
            "normalized_density",
            "binary64_legacy_v1",
        ),
        (
            background(),
            "model_admission",
            "compiled_flat_flrw_requested_projections",
            "binary64_storage_longdouble_intermediate",
        ),
    ]
}
#[test]
fn operation_owned_metadata_and_assurance_identity() {
    let s = scratch();
    for (request, id, method, precision) in cases(&s) {
        let (r, code) = run(&s, &request, None, false);
        assert_eq!(code, 0, "{}: {r}", request["operation"]);
        let receipt = &r["receipt"];
        assert_eq!(receipt["execution"], "completed");
        assert_eq!(receipt["accepted"], true);
        assert_eq!(receipt["accepted_scope"], "numerical_contract");
        assert_eq!(receipt["assurance"]["satisfied"], true);
        assert_eq!(receipt["outputs"][0]["id"], id);
        assert_eq!(
            receipt["outputs"][0]["check_kind"],
            if id == "prepared_observations" {
                "structural_source_contract"
            } else {
                "numerical_contract"
            }
        );
        assert_eq!(receipt["method"], method);
        assert_eq!(receipt["precision"], precision);
        assert_eq!(receipt["interpretation"], "unqualified");
        let spec = object(&s, &receipt["scientific_specification_digest"]);
        assert_eq!(spec["operation"], request["operation"]);
        assert_eq!(
            spec["requested_outputs"][0]["id"],
            receipt["outputs"][0]["id"]
        );
        assert_eq!(
            spec["requested_outputs"][0]["required"],
            receipt["outputs"][0]["required"]
        );
        assert_eq!(spec["requested_outputs"][0]["numerical_gate"], "required");
        let (strict, code) = run(&s, &request, Some("qualified"), false);
        assert_eq!(code, 6, "{strict}");
        assert_eq!(strict["receipt"]["accepted"], false);
        assert_eq!(strict["receipt"]["assurance"]["numerical_contract"], true);
        assert_eq!(
            strict["receipt"]["assurance"]["qualification"],
            "unqualified"
        );
        assert_eq!(
            strict["receipt"]["scientific_specification_digest"],
            receipt["scientific_specification_digest"]
        );
        assert_eq!(strict["receipt"]["output_digest"], receipt["output_digest"]);
        assert_ne!(
            strict["receipt"]["execution_identity"],
            receipt["execution_identity"]
        );
    }
    let exact = fixture("exact-add");
    for assurance in [None, Some("qualified")] {
        let (r, code) = run(&s, &exact, assurance, false);
        assert_eq!(code, 0);
        assert_eq!(r["receipt"]["outputs"][0]["numerical"], "checks_passed");
        assert_eq!(
            r["receipt"]["outputs"][0]["check_kind"],
            "exact_integer_contract"
        );
        assert_eq!(r["receipt"]["interpretation"], "not_applicable");
        let spec = object(&s, &r["receipt"]["scientific_specification_digest"]);
        assert_eq!(
            spec["requested_outputs"][0]["id"],
            r["receipt"]["outputs"][0]["id"]
        );
        assert_eq!(
            spec["requested_outputs"][0]["required"],
            r["receipt"]["outputs"][0]["required"]
        );
        assert_eq!(spec["requested_outputs"][0]["numerical_gate"], "required");
    }
}
#[test]
fn completed_mixed_failures_are_not_transport_or_internal_failures() {
    let s = scratch();
    let mut q = gaussian();
    q["residuals"] = json!([[1., 2.], [f64::MAX, -f64::MAX]]);
    for assurance in [None, Some("qualified")] {
        let (r, code) = run(&s, &q, assurance, false);
        assert_eq!(code, 2);
        assert_eq!(r["receipt"]["execution"], "completed");
        assert_eq!(r["receipt"]["numerical"], "failed");
        assert_eq!(r["receipt"]["assurance"]["satisfied"], false);
        assert_eq!(r["result"]["rows"][0]["kind"], "finite");
        assert_eq!(r["result"]["rows"][1]["kind"], "failure");
        assert!(r["result"]["rows"][1].get("quadratic").is_none());
    }
    let mut malformed = background();
    malformed["ignored"] = json!(true);
    let (r, code) = run(&s, &malformed, None, false);
    assert_eq!(code, 2);
    assert_eq!(r["receipt"]["execution"], "failed");
    assert_eq!(r["receipt"]["failure_class"], "request_or_transport");
    assert_eq!(r["receipt"]["numerical"], "not_assessed");
    assert_eq!(r["receipt"]["outputs"], json!([]));
    let (r, code) = run(&s, &fixture("exact-add"), None, true);
    assert_eq!(code, 2);
    assert_eq!(r["receipt"]["execution"], "failed");
    assert_eq!(r["receipt"]["failure_class"], "internal");
    assert_eq!(r["result"]["error_id"], "RUST_PANIC");
}

#[test]
fn well_formed_json_with_native_descriptor_rejection_is_transport_failure() {
    let s = scratch();
    let mut q = background();
    // Structurally valid JSON, but above the native policy hard cap. This is
    // descriptor admission, unlike an owned per-row numerical failure.
    q["numerical_policy"]["maximum_models"] = json!(65537);
    let (r, code) = run(&s, &q, None, false);
    assert_eq!(code, 2);
    assert_eq!(r["receipt"]["execution"], "failed");
    assert_eq!(r["receipt"]["failure_class"], "request_or_transport");
    assert_eq!(r["receipt"]["numerical"], "not_assessed");
    assert_eq!(r["result"]["error_id"], "CORE_STATUS_2");
    assert_eq!(r["receipt"]["outputs"], json!([]));
}

#[test]
fn request_schema_two_accepts_primitives_and_rejects_retired_one() {
    let s = scratch();
    for name in ["exact-add", "quantity-length", "numerical-sum"] {
        let current = fixture(name);
        let (r, code) = run(&s, &current, None, false);
        assert_eq!(code, 0, "{r}");
        let mut retired = current;
        retired["schema_version"] = json!(1);
        let (r, code) = run(&s, &retired, None, false);
        assert_eq!(code, 2, "{r}");
        assert_eq!(r["receipt"]["execution"], "failed");
        assert_eq!(r["receipt"]["accepted"], false);
    }
}
