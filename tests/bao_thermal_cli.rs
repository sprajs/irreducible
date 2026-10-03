//! Original strict thermal input and durable per-call evidence; no Rust physics.
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
    serde_json::from_str(include_str!("fixtures/bao-thermal.json")).unwrap()
}
fn run_raw(bytes: &[u8]) -> (Value, i32, Option<Value>) {
    // Controller may choose an ignored evidence root. The CI default is also
    // ignored and retains this bounded suite's per-call stores even on panic.
    let parent = std::env::var_os("IRRED_TEST_ARTIFACT_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("build/test-artifacts/thermal-bao-cli")
        });
    fs::create_dir_all(&parent).unwrap();
    let d = Scratch(parent.join(format!(
            "irred-thermal-bao-{}-{}",
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
    let launched = Command::new(env!("CARGO_BIN_EXE_irred"))
        .args(["run", p.to_str().unwrap(), store.to_str().unwrap()])
        .output();
    let out = match launched {
        Ok(out) => out,
        Err(error) => {
            fs::write(d.0.join("launch-error.json"), serde_json::to_vec(&json!({"error":error.to_string(),"request_sha256":format!("{:x}",Sha256::digest(bytes))})).unwrap()).unwrap();
            panic!(
                "child launch failed; original request and error retained at {}",
                d.0.display()
            );
        }
    };
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
fn typed_resolution(actual: &Value, original: &Value, numeric: bool) {
    match (actual, original) {
        (Value::Number(a), Value::Number(b)) if numeric => {
            assert_eq!(a.as_f64().unwrap().to_bits(), b.as_f64().unwrap().to_bits())
        }
        (Value::Array(a), Value::Array(b)) => {
            assert_eq!(a.len(), b.len());
            for (a, b) in a.iter().zip(b) {
                typed_resolution(a, b, numeric);
            }
        }
        (Value::Object(a), Value::Object(b)) => {
            assert_eq!(a.len(), b.len());
            for (k, b) in b {
                let is_numeric = matches!(
                    k.as_str(),
                    "redshift"
                        | "observed_ratios"
                        | "covariance_row_major"
                        | "h0_km_s_mpc"
                        | "physical_baryon_density"
                        | "physical_cdm_density"
                        | "tcmb_kelvin"
                        | "physical_massless_nonphoton_density"
                        | "mass_ev"
                        | "temperature_today_kelvin"
                        | "statistical_weight"
                        | "z_drag"
                        | "absolute_tolerance_mpc"
                        | "relative_tolerance"
                        | "absolute_tolerance_ratio"
                        | "relative_tolerance_ratio"
                        | "absolute_tolerance"
                        | "maximum_forward_sensitivity"
                        | "maximum_projection_log_density_error"
                );
                typed_resolution(a.get(k).unwrap(), b, is_numeric);
            }
        }
        _ => assert_eq!(actual, original),
    }
}
// These active tests earn structural/refusal lineage only. The profile gate
// cannot be treated as a successful numerical CLI consumer.
#[test]
fn unsupported_profile_preserves_original_bytes_and_completed_refusal() {
    let raw = include_bytes!("fixtures/bao-thermal.json");
    let (v, status, spec) = run_raw(raw);
    assert_eq!(status, 6);
    assert_eq!(v["result"]["kind"], "failure");
    assert_eq!(v["result"]["error_id"], "NUMERICAL_QUALIFICATION_REQUIRED");
    assert_eq!(v["result"]["cause"], "unsupported_allocation_profile");
    assert_eq!(v["result"]["native_payload_absent"], true);
    assert_eq!(v["result"]["source_prepare_call_attempted"], false);
    assert_eq!(v["result"]["thermal_batch_call_attempted"], false);
    assert_eq!(v["receipt"]["precision"]["requested"], "wide");
    assert_eq!(v["receipt"]["method"]["requested_output_mask"], 7);
    assert_eq!(
        spec.unwrap()["raw_input_sha256"],
        format!("{:x}", Sha256::digest(raw))
    );
}
#[test]
fn wrong_order_and_units_are_structural() {
    for key in ["covariance_axis_ids", "ratio_unit"] {
        let mut r = fixture();
        if key == "ratio_unit" {
            r["observation"][key] = json!("Mpc");
        } else {
            r["observation"][key][0] = json!("synthetic-DV-z1");
        }
        let (v, status, _) = run(&r);
        assert_eq!(status, 2);
        assert_eq!(
            v["result"]["error_id"],
            if key == "ratio_unit" {
                "THERMAL_BAO_STRUCTURE"
            } else {
                "THERMAL_BAO_ORDER"
            }
        );
    }
}
#[test]
fn duplicates_unknown_fields_and_oversized_typed_arrays_refuse() {
    let mut r = fixture();
    r["unexpected"] = json!(true);
    let (_, status, _) = run(&r);
    assert_eq!(status, 2);
    let mut r = fixture();
    r["models"] = json!([r["models"][0].clone(), r["models"][0].clone()]);
    let (_, status, _) = run(&r);
    assert_eq!(status, 2);
    let mut r = fixture();
    r["observation"]["observed_ratios"] = json!(vec![0; 65]);
    let (_, status, _) = run(&r);
    assert_eq!(status, 2);
}
#[test]
fn finite_positive_subnormal_weight_is_not_a_wrapper_normal_gate() {
    let mut r = fixture();
    r["models"][0]["physical_model"]["species"][0]["statistical_weight"] = json!(f64::from_bits(1));
    let (v, status, _) = run(&r);
    assert_eq!(status, 6);
    assert_eq!(v["result"]["cause"], "unsupported_allocation_profile");
}
#[test]
fn discovery_preserves_free_ruler_and_names_profile_gate() {
    let out = Command::new(env!("CARGO_BIN_EXE_irred"))
        .args(["describe", "--json"])
        .output()
        .unwrap();
    assert!(out.status.success());
    let v: Value = serde_json::from_slice(&out.stdout).unwrap();
    let all = v["capabilities"].as_array().unwrap();
    assert!(all.iter().any(|x| x["id"] == "bao.density"));
    let t = all
        .iter()
        .find(|x| x["id"] == "bao.thermal_density")
        .unwrap();
    assert_eq!(t["qualification"], "unqualified");
    assert_eq!(t["source_semantics"], "synthetic_controls");
    assert_eq!(
        t["execution_gate"],
        "unsupported allocation profile refuses before native execution"
    );
}
// Explicit remaining numerical gates are retained as unexecuted controls,
// separate from active source/profile refusal tests. Root must earn profile
// admission before enabling them; an ignored test is never a numerical pass.
#[test]
#[ignore = "allocation-profile proof and root runtime admission pending"]
fn original_native_cli_named_cofactor_parity() {
    let r = fixture();
    let (v, status, spec) = run(&r);
    assert_eq!(status, 0);
    let spec = spec.unwrap();
    typed_resolution(&spec, &r, false);
    let models = v["result"]["models"].as_array().unwrap();
    assert_eq!(models.len(), 2);
    let inverse = [29., -6., 1., -6., 24., -4., 1., -4., 19.];
    for (i, m) in models.iter().enumerate() {
        assert_eq!(m["model_id"], r["models"][i]["id"]);
        let p = m["predictions"].as_array().unwrap();
        let e = m["residuals"].as_array().unwrap();
        assert_eq!(p.len(), 3);
        assert_eq!(e.len(), 3);
        let mut q = 0.;
        for j in 0..3 {
            assert_eq!(
                e[j].as_f64().unwrap().to_bits(),
                (r["observation"]["observed_ratios"][j].as_f64().unwrap() - p[j].as_f64().unwrap())
                    .to_bits()
            );
            for k in 0..3 {
                q += e[j].as_f64().unwrap() * inverse[j * 3 + k] * e[k].as_f64().unwrap() / 1760.;
            }
        }
        assert!((m["density"]["quadratic"].as_f64().unwrap() - q).abs() <= 1e-8);
        assert!(
            (m["density"]["log_determinant"].as_f64().unwrap() - 450560_f64.ln()).abs() <= 1e-8
        );
        let norm = 3. * (2. * std::f64::consts::PI).ln();
        assert!((m["density"]["log_normalization"].as_f64().unwrap() - norm).abs() <= 1e-8);
        assert!(
            (m["density"]["log_density"].as_f64().unwrap() + (q + 450560_f64.ln() + norm) / 2.)
                .abs()
                <= 1e-8
        );
    }
}
#[test]
#[ignore = "allocation-profile proof and root runtime admission pending"]
fn completed_native_global_and_projection_refusals() {
    for cap in [
        "maximum_native_bytes",
        "maximum_preparation_native_bytes",
        "maximum_evaluation_native_bytes",
    ] {
        let mut r = fixture();
        r["resource_policy"][cap] = json!(0);
        if cap == "maximum_native_bytes" {
            r["resource_policy"]["maximum_preparation_native_bytes"] = json!(0);
            r["resource_policy"]["maximum_evaluation_native_bytes"] = json!(0);
        }
        let (v, status, _) = run(&r);
        assert_eq!(status, 6);
        assert_eq!(v["result"]["error_id"], "NUMERICAL_QUALIFICATION_REQUIRED");
    }
    let mut r = fixture();
    r["resource_policy"]["maximum_projection_log_density_error"] = json!(1e-30);
    let (v, status, _) = run(&r);
    assert_eq!(status, 6);
    assert!(v["result"]["models"][0]["predictions"].is_array());
    assert!(v["result"]["models"][0]["residuals"].is_array());
    assert!(v["result"]["models"][0]["density"].is_null());
}
