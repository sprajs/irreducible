//! Original fixed future law and operation-owned records; no Rust predictive equations.
use serde_json::{Value, json};
use sha2::{Digest, Sha256};
use std::{
    collections::BTreeMap,
    fs,
    path::PathBuf,
    process::Command,
    time::{SystemTime, UNIX_EPOCH},
};
struct Scratch(PathBuf);
fn fixture() -> Value {
    serde_json::from_str(include_str!("fixtures/gaussian-predictive.json")).unwrap()
}
fn run_raw(bytes: &[u8]) -> (Value, i32, Option<Value>) {
    // Controller may choose an ignored evidence root. The CI default is also
    // ignored and retains this bounded suite's per-call stores even on panic.
    let parent = std::env::var_os("IRRED_TEST_ARTIFACT_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                .join("build/test-artifacts/predictive-gaussian-cli")
        });
    fs::create_dir_all(&parent).unwrap();
    let d = Scratch(parent.join(format!(
            "irred-predictive-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        )));
    fs::create_dir(&d.0).unwrap();
    let p = d.0.join("request.json");
    let store = d.0.join("store");
    fs::write(&p, bytes).unwrap();
    let out = Command::new(env!("CARGO_BIN_EXE_irred"))
        .args(["run", p.to_str().unwrap(), store.to_str().unwrap()])
        .output()
        .unwrap();
    fs::write(d.0.join("stdout.json"), &out.stdout).unwrap();
    fs::write(d.0.join("stderr.log"), &out.stderr).unwrap();
    let exit_status = json!({
        "exit_status_display": out.status.to_string(),
        "success": out.status.success(),
        "exit_code": out.status.code(),
    });
    #[cfg(unix)]
    let exit_status = {
        use std::os::unix::process::ExitStatusExt;
        let mut exit_status = exit_status;
        exit_status["unix_signal"] = json!(out.status.signal());
        exit_status["unix_raw_status"] = json!(out.status.into_raw());
        exit_status
    };
    fs::write(
        d.0.join("exit-status.json"),
        serde_json::to_vec(&exit_status).unwrap(),
    )
    .unwrap();
    fn inventory(
        root: &std::path::Path,
        dir: &std::path::Path,
        files: &mut BTreeMap<String, String>,
    ) {
        for entry in fs::read_dir(dir).unwrap() {
            let path = entry.unwrap().path();
            if path.is_dir() {
                inventory(root, &path, files);
            } else if path.is_file() {
                files.insert(
                    path.strip_prefix(root)
                        .unwrap()
                        .to_string_lossy()
                        .into_owned(),
                    format!("{:x}", Sha256::digest(fs::read(path).unwrap())),
                );
            }
        }
    }
    let mut hashes = BTreeMap::new();
    inventory(&d.0, &d.0, &mut hashes);
    fs::write(d.0.join("capture.json"), serde_json::to_vec_pretty(&json!({"request_sha256":format!("{:x}", Sha256::digest(bytes)),"files":hashes,"exit_status":exit_status})).unwrap()).unwrap();
    let v: Value = serde_json::from_slice(&out.stdout).unwrap();
    let digest = v["receipt"]["input_digest"].as_str().unwrap();
    assert_eq!(fs::read(store.join("objects").join(digest)).unwrap(), bytes);
    let spec = v["receipt"]["scientific_specification_digest"]
        .as_str()
        .map(|x| {
            serde_json::from_slice(&fs::read(store.join("objects").join(x)).unwrap()).unwrap()
        });
    (v, out.status.code().unwrap(), spec)
}
fn run(v: &Value) -> (Value, i32, Option<Value>) {
    run_raw(&serde_json::to_vec(v).unwrap())
}
fn near(v: &Value, x: f64) {
    assert!((v.as_f64().unwrap() - x).abs() <= 2e-12 * (1. + x.abs()));
}

fn failed_request(input: &Value) {
    let (v, code, spec) = run(input);
    assert_ne!(code, 0);
    assert_eq!(v["receipt"]["execution"], "failed");
    assert!(spec.is_none());
}
#[test]
fn named_joint_law_preserves_mean_density_order_and_receipts() {
    let input = fixture();
    let (v, code, spec) = run(&input);
    assert_eq!(code, 0);
    assert_eq!(v["receipt"]["execution"], "completed");
    assert_eq!(v["receipt"]["accepted_scope"], "numerical_contract");
    assert_eq!(v["receipt"]["rng"], "not_applicable");
    let o = &v["result"];
    near(&o["rows"][0]["future_mean"][0], -258. / 668.);
    near(&o["rows"][0]["future_mean"][1], 855. / 668.);
    near(&o["rows"][1]["future_mean"][0], 0.);
    near(&o["rows"][1]["future_mean"][1], 0.75);
    near(&o["rows"][1]["joint_density"]["quadratic"], 0.);
    let scalar_product =
        -0.5 * ((4444_f64 * 2618. / (1336. * 1336.)).ln() + 2. * (2. * std::f64::consts::PI).ln());
    assert!(
        o["rows"][1]["joint_density"]["log_density"]
            .as_f64()
            .unwrap()
            - scalar_product
            > 1e-3
    );
    assert_eq!(
        o["ordered_future_row_ids"],
        input["future_noise"]["ordered_row_ids"]
    );
    assert_eq!(o["future_event_ids"], input["future_noise"]["event_ids"]);
    assert_eq!(
        o["ordered_parameter_ids"],
        input["parameter_prior"]["ordered_parameter_ids"]
    );
    assert!(o.get("covariance").is_none());
    assert_eq!(v["receipt"]["outputs"].as_array().unwrap().len(), 2);
    assert_eq!(v["receipt"]["outputs"][0]["id"], "predictive_means");
    assert_eq!(
        v["receipt"]["outputs"][1]["id"],
        "joint_predictive_log_densities"
    );
    for field in [
        "training_noise_attempted",
        "training_noise_prepared",
        "posterior_attempted",
        "posterior_prepared",
        "future_noise_attempted",
        "future_noise_prepared",
        "predictive_attempted",
        "predictive_prepared",
        "batch_called",
        "batch_admission_completed",
    ] {
        assert_eq!(v["receipt"]["method"][field], true);
    }
    assert_eq!(
        v["receipt"]["resource_budget"]["operation"]["actual"]["declared_total_work"],
        7296
    );
    assert_eq!(
        v["receipt"]["resource_budget"]["operation"]["actual"]["pooled_numeric_elements"],
        50
    );
    let spec = spec.unwrap();
    assert_eq!(
        spec["conditioning"]["vectors"],
        input["conditioning"]["vectors"]
    );
    assert_eq!(
        spec["future_vectors"]["vectors"],
        input["future_vectors"]["vectors"]
    );
    assert_eq!(spec["prediction"], input["prediction"]);
}
#[test]
fn masks_omit_unused_payload_and_keep_core_request_replayable() {
    for (means, density) in [(true, false), (false, true), (true, true)] {
        let mut input = fixture();
        input["outputs"] = json!({"means":means,"joint_log_densities":density});
        if !density {
            input.as_object_mut().unwrap().remove("future_vectors");
        }
        let (v, code, spec) = run(&input);
        assert_eq!(code, 0);
        assert_eq!(
            v["receipt"]["outputs"].as_array().unwrap().len(),
            usize::from(means) + usize::from(density)
        );
        for row in v["result"]["rows"].as_array().unwrap() {
            assert_eq!(row.get("future_mean").is_some(), means);
            assert_eq!(row.get("absolute_error_estimates").is_some(), means);
            assert_eq!(row.get("joint_density").is_some(), density);
        }
        assert_eq!(
            v["result"]["groups"]["predictive_means"],
            if means {
                "checks_passed"
            } else {
                "not_requested"
            }
        );
        assert_eq!(
            v["result"]["groups"]["joint_predictive_log_densities"],
            if density {
                "checks_passed"
            } else {
                "not_requested"
            }
        );
        let mut replay = spec.unwrap();
        assert_eq!(replay.get("future_vectors").is_some(), density);
        replay.as_object_mut().unwrap().remove("requested_outputs");
        replay.as_object_mut().unwrap().remove("batch_meaning");
        let (_, code, _) = run(&replay);
        assert_eq!(code, 0);
    }
}
#[test]
fn per_case_failure_is_retained_without_zero_payload_or_renormalization() {
    let mut input = fixture();
    input["conditioning"]["vectors"] = json!([[1.25, -0.5], [1e-310, 0.]]);
    let (v, code, spec) = run(&input);
    assert_eq!(code, 2);
    assert_eq!(v["receipt"]["execution"], "completed");
    assert_eq!(v["result"]["rows"].as_array().unwrap().len(), 2);
    assert_eq!(v["result"]["rows"][0]["kind"], "finite");
    assert_eq!(v["result"]["rows"][1]["kind"], "failure");
    assert_eq!(v["result"]["rows"][1]["case_index"], 1);
    assert!(v["result"]["rows"][1].get("future_mean").is_none());
    assert!(v["result"]["rows"][1].get("joint_density").is_none());
    for o in v["receipt"]["outputs"].as_array().unwrap() {
        assert_eq!(o["numerical"], "failed");
    }
    assert_eq!(
        spec.unwrap()["conditioning"]["vectors"],
        input["conditioning"]["vectors"]
    );
}
#[test]
fn global_spd_failures_name_actual_stage_without_claiming_completion() {
    for (field, phase, completed) in [
        (
            "noise",
            "training_noise_preparation",
            "training_noise_prepared",
        ),
        (
            "parameter_prior",
            "posterior_preparation",
            "posterior_prepared",
        ),
        (
            "future_noise",
            "future_noise_preparation",
            "future_noise_prepared",
        ),
    ] {
        let mut input = fixture();
        input[field]["covariance_row_major"] = json!([1., 2., 2., 1.]);
        let (v, code, _) = run(&input);
        assert_eq!(code, 2);
        assert_eq!(v["receipt"]["execution"], "completed");
        assert_eq!(v["result"]["phase"], phase);
        assert_eq!(v["result"]["rows"], json!([]));
        assert_eq!(v["receipt"]["method"][completed], false);
        assert_eq!(v["receipt"]["method"]["executed"], true);
        assert_eq!(v["receipt"]["method"]["batch_called"], false);
    }
}
#[test]
fn original_smaller_quotas_are_completed_causal_refusals() {
    for (field, value) in [
        ("maximum_native_bytes", 0),
        ("maximum_work_units", 7295),
        ("maximum_elements", 49),
        ("maximum_cases", 1),
        ("maximum_string_bytes", 0),
    ] {
        let mut input = fixture();
        input["resource_policy"][field] = json!(value);
        let (v, code, spec) = run(&input);
        assert_eq!(code, 2);
        assert_eq!(v["receipt"]["execution"], "completed");
        assert_eq!(v["result"]["numerical_status"], "work_limit");
        assert_eq!(v["result"]["rows"], json!([]));
        assert_eq!(v["receipt"]["method"]["executed"], false);
        assert_eq!(v["receipt"]["method"]["actual"], Value::Null);
        assert_eq!(spec.unwrap()["resource_policy"][field], value);
    }
    let mut input = fixture();
    input["conditioning"]["vectors"] = json!([]);
    input["conditioning"]["case_ids"] = json!([]);
    input["future_vectors"]["vectors"] = json!([]);
    let (v, code, _) = run(&input);
    assert_eq!(code, 2);
    assert_eq!(v["result"]["phase"], "admission");
    assert_eq!(v["receipt"]["method"]["executed"], false);
}
#[test]
fn supplied_independence_units_and_disjoint_lineage_are_mandatory() {
    for variant in 0..6 {
        let mut input = fixture();
        match variant {
            0 => input["prediction"]["future_noise_independence_declared"] = json!(false),
            1 => input["prediction"]["noise_conditional_on_parameters_declared"] = json!(false),
            2 => input["future_response"]["ordered_parameter_ids"] = json!(["slope", "offset"]),
            3 => input["prediction"]["future_covariance_unit"] = json!("other^2"),
            4 => {
                input["future_noise"]["event_ids"] = input["noise"]["event_ids"].clone();
                input["future_vectors"]["event_ids"] = input["noise"]["event_ids"].clone();
            }
            _ => {
                input["future_noise"]["ordered_row_ids"] =
                    input["noise"]["ordered_row_ids"].clone();
                input["future_response"]["ordered_row_ids"] =
                    input["noise"]["ordered_row_ids"].clone();
                input["future_vectors"]["ordered_row_ids"] =
                    input["noise"]["ordered_row_ids"].clone();
            }
        }
        let (v, code, _) = run(&input);
        assert_eq!(code, 2);
        assert_eq!(v["result"]["status"], "incompatible_metadata");
        assert_eq!(v["receipt"]["execution"], "completed");
    }
}
#[test]
fn null_design_or_response_is_lawful_with_proper_noise_and_prior() {
    for field in ["design", "future_response"] {
        let mut input = fixture();
        input[field]["values_row_major"] = json!([0., 0., 0., 0.]);
        let (v, code, _) = run(&input);
        assert_eq!(code, 0);
        near(&v["result"]["rows"][1]["future_mean"][0], 0.);
        near(
            &v["result"]["rows"][1]["future_mean"][1],
            if field == "design" { 0.75 } else { 0. },
        );
    }
}
#[test]
fn strict_objects_presence_duplicates_and_hard_request_bounds() {
    let original = fixture();
    for field in [
        "noise",
        "design",
        "parameter_prior",
        "conditioning",
        "future_noise",
        "future_response",
        "prediction",
        "outputs",
        "future_vectors",
        "resource_policy",
    ] {
        for replacement in [Value::Null, json!([]), json!([1, 2, 3])] {
            let mut input = original.clone();
            input[field] = replacement;
            failed_request(&input);
        }
        let mut input = original.clone();
        input[field]["unused"] = json!(true);
        failed_request(&input);
    }
    failed_request(&json!([original]));
    let mut input = fixture();
    input["outputs"] = json!({"means":false,"joint_log_densities":false});
    failed_request(&input);
    input = fixture();
    input.as_object_mut().unwrap().remove("future_vectors");
    failed_request(&input);
    input = fixture();
    input["outputs"]["joint_log_densities"] = json!(false);
    failed_request(&input);
    input = fixture();
    input["future_noise"]["noise_identity"] = json!("x".repeat(257));
    failed_request(&input);
    input = fixture();
    input["resource_policy"]["maximum_native_bytes"] = json!(268435457);
    failed_request(&input);
    input = fixture();
    input["conditioning"]["vectors"] = json!([[1.25], [0.5, 0.25]]);
    failed_request(&input);
    for field in [
        "noise",
        "design",
        "parameter_prior",
        "conditioning",
        "future_noise",
        "future_response",
        "prediction",
        "outputs",
        "future_vectors",
        "resource_policy",
    ] {
        let input = fixture();
        let object = &input[field];
        let (key, value) = object.as_object().unwrap().iter().next().unwrap();
        let original_object = serde_json::to_string(object).unwrap();
        let duplicate = format!(
            "{{\"{key}\":{},{}",
            serde_json::to_string(value).unwrap(),
            &original_object[1..]
        );
        let raw = serde_json::to_string(&input)
            .unwrap()
            .replacen(&original_object, &duplicate, 1);
        let (v, code, spec) = run_raw(raw.as_bytes());
        assert_ne!(code, 0);
        assert_eq!(v["receipt"]["execution"], "failed");
        assert!(spec.is_none());
    }
    input = fixture();
    input["future_noise"]["covariance_row_major"] = json!(vec![0.; 1000000]);
    failed_request(&input);
    input = fixture();
    input["conditioning"]["case_ids"] = json!(vec!["case"; 65537]);
    failed_request(&input);
    let raw = serde_json::to_string(&fixture()).unwrap().replace(
        "\"future_noise_independence_declared\":true",
        "\"future_noise_independence_declared\":false,\"future_noise_independence_declared\":true",
    );
    let (v, code, spec) = run_raw(raw.as_bytes());
    assert_ne!(code, 0);
    assert_eq!(v["receipt"]["execution"], "failed");
    assert!(spec.is_none());
    let raw = serde_json::to_string(&fixture())
        .unwrap()
        .replace("0.125", "1e400");
    let (v, code, spec) = run_raw(raw.as_bytes());
    assert_ne!(code, 0);
    assert_eq!(v["receipt"]["execution"], "failed");
    assert!(spec.is_none());
}
#[test]
fn current_discovery_exposes_single_predictive_operation() {
    let out = Command::new(env!("CARGO_BIN_EXE_irred"))
        .args(["describe", "--json"])
        .output()
        .unwrap();
    assert!(out.status.success());
    let v: Value = serde_json::from_slice(&out.stdout).unwrap();
    let matches: Vec<_> = v["capabilities"]
        .as_array()
        .unwrap()
        .iter()
        .filter(|x| x["id"] == "statistics.gaussian_predictive")
        .collect();
    assert_eq!(matches.len(), 1);
    assert_eq!(matches[0]["qualification"], "unqualified");
    assert_eq!(
        matches[0]["outputs"],
        json!(["predictive_means", "joint_predictive_log_densities"])
    );
}
