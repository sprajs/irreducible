//! Current retained-consumer boundary controls. Analytic values are unchanged;
//! transport parity is separate from the directed original-data science gates.
use serde_json::{Value, json};
use sha2::{Digest, Sha256};
use std::{
    fs,
    path::PathBuf,
    process::Command,
    time::{SystemTime, UNIX_EPOCH},
};
fn root() -> PathBuf {
    let p = std::env::temp_dir().join(format!(
        "irred-current-consumer-{}-{}",
        std::process::id(),
        SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    fs::create_dir_all(&p).unwrap();
    p
}
fn run(p: &PathBuf, r: &Value) -> (Value, i32) {
    let path = p.join("request.json");
    fs::write(&path, serde_json::to_vec(r).unwrap()).unwrap();
    let o = Command::new(env!("CARGO_BIN_EXE_irred"))
        .arg("run")
        .arg(path)
        .arg(p.join("store"))
        .output()
        .unwrap();
    (
        serde_json::from_slice(&o.stdout)
            .unwrap_or_else(|_| panic!("{}", String::from_utf8_lossy(&o.stderr))),
        o.status.code().unwrap(),
    )
}
fn record(p: &PathBuf, o: &Value, r: &Value, method: &str) {
    let receipt = &o["receipt"];
    assert_eq!(receipt["method"], method);
    assert_eq!(receipt["precision"], "wide");
    assert_eq!(receipt["interpretation"], "unqualified");
    let requested = r["requested_outputs"].as_array().unwrap();
    let outputs = receipt["outputs"].as_array().unwrap();
    assert_eq!(outputs.len(), requested.len());
    for (a, b) in requested.iter().zip(outputs) {
        assert_eq!(a, &b["id"]);
        assert_eq!(b["required"], true);
        assert_eq!(b["check_kind"], "numerical_contract");
    }
    for key in ["scientific_specification_digest", "output_digest"] {
        let d = receipt[key].as_str().unwrap();
        let bytes = fs::read(p.join("store/objects").join(d)).unwrap();
        assert_eq!(format!("{:x}", Sha256::digest(&bytes)), d);
        let value: Value = serde_json::from_slice(&bytes).unwrap();
        if key == "output_digest" {
            assert_eq!(value, o["result"]);
        } else {
            assert_eq!(value["operation"], r["operation"]);
            assert_eq!(value["requested_outputs"], r["requested_outputs"]);
        }
    }
}
fn bao() -> Value {
    let scale = 299792.458 / 10000.;
    json!({"schema_version":2,"operation":"bao.density","source":{"profile":"synthetic_inline_v1","rows":[{"id":"a","z":0.1,"observable":"DM_over_rs","value":scale*0.1+1.0},{"id":"b","z":0.2,"observable":"DH_over_rs","value":scale+2.0}],"covariance":[4.0,1.0,1.0,9.0],"covariance_axis_ids":["a","b"],"ordering_provenance":"supplied synthetic order","calibration_provenance":"synthetic fixture; unknown calibration","dependence_provenance":"declared correlated covariance; cross-probe dependence unknown","redshift_convention":"P01/released-effective-redshift/v1","ruler_convention":"P01/free-H0rd-km-s-no-early-physics/v1"},"preparation_policy":{"arithmetic":"wide","maximum_queries":16,"maximum_matrix_elements":256,"maximum_string_bytes":4096,"maximum_native_bytes":8388608,"maximum_forward_sensitivity":1e-10},"models":[{"expansion":{"kind":"constant_q","q":-1.0},"geometry":{"kind":"flat_flrw"},"h0_rd_km_s":10000.0}],"requested_outputs":["normalized_density","predictions","residuals"],"numerical_policy":{"arithmetic":"wide","projection":{"maximum_queries":16,"maximum_callbacks":100000,"maximum_segment_visits":2000,"integration":{"absolute_tolerance":1e-13,"relative_tolerance":1e-12,"maximum_evaluations":100000,"maximum_depth":30}},"maximum_models":4,"maximum_array_elements":256,"maximum_native_bytes":8388608,"maximum_forward_sensitivity":1e-10}})
}
#[test]
fn current_bao_normalized_analytic_and_mixed_failure() {
    let p = root();
    let r = bao();
    let (o, code) = run(&p, &r);
    assert_eq!(code, 0, "{o}");
    record(&p, &o, &r, "normalized_conditional_free_ruler_gaussian");
    let d = &o["result"]["evaluations"][0]["normalized_density"]["value"];
    assert!((d["quadratic"].as_f64().unwrap() - 0.6).abs() < 4e-11);
    assert!((d["log_determinant"].as_f64().unwrap() - 35f64.ln()).abs() < 4e-11);
    assert!(
        (d["log_density"].as_f64().unwrap()
            + 0.5 * (0.6 + 35f64.ln() + 2.0 * (2.0 * std::f64::consts::PI).ln()))
        .abs()
            < 4e-11
    );
    assert_eq!(o["receipt"]["execution"], "completed");
    assert_eq!(o["receipt"]["accepted"], true);
    let mut mixed = r.clone();
    mixed["models"].as_array_mut().unwrap().push(json!({"expansion":{"kind":"constant_q","q":3.0},"geometry":{"kind":"flat_flrw"},"h0_rd_km_s":10000.0}));
    let (o, code) = run(&p, &mixed);
    assert_eq!(code, 2, "{o}");
    assert_eq!(o["receipt"]["execution"], "completed");
    assert_eq!(o["receipt"]["accepted"], false);
    assert!(o["result"]["evaluations"][0]["normalized_density"]["value"].is_object());
    assert!(o["result"]["evaluations"][1]["normalized_density"]["value"].is_null());
    let mut predictions = r;
    predictions["requested_outputs"] = json!(["predictions"]);
    let (o, code) = run(&p, &predictions);
    assert_eq!(code, 0, "{o}");
    assert!(o["result"]["evaluations"][0].get("residuals").is_none());
    assert!(
        o["result"]["evaluations"][0]
            .get("normalized_density")
            .is_none()
    );
    let _ = fs::remove_dir_all(p);
}
fn sn(p: &PathBuf) -> Value {
    fs::write(
        p.join("table.txt"),
        format!(
            "ROW_ID EVENT_ID MAG\nb eb {:.17}\na ea {:.17}\n",
            5.0 * 2f64.log10() + 1.0,
            5.0 * 6f64.log10() + 2.0
        ),
    )
    .unwrap();
    fs::write(p.join("cov.txt"), "2\n3 .25 .25 2\n").unwrap();
    json!({"schema_version":2,"operation":"supernova.profile","observations":{"schema_version":2,"operation":"observations.prepare","table":p.join("table.txt"),"uncertainty":p.join("cov.txt"),"metadata":{"profile":"typed_magnitude_covariance","role":"observed_measurement","unit":"magnitude","calibration":"unknown","uncertainty":"covariance","uncertainty_unit":"magnitude_squared","component":"unknown"},"resources":{"maximum_asset_bytes":4096,"maximum_rows":8,"maximum_matrix_elements":64,"maximum_string_bytes":4096,"maximum_preparation_bytes":8388608},"exports":[],"calibration_provenance":"explicit uncalibrated synthetic measured-scale control","dependence_provenance":"unknown","quality_dictionary":"not supplied","ordering_provenance":"supplied row order","source_selection":[1,1]},"selection":{"kind":"explicit","source_indices":[0,1],"coordinates":[{"z_expansion":1.0,"observer":{"redshift":1.0,"convention":"geometric_same_redshift"}},{"z_expansion":2.0,"observer":{"redshift":2.0,"convention":"geometric_same_redshift"}}]},"preparation_policy":{"arithmetic":"wide","maximum_selected_rows":2,"maximum_matrix_elements":4,"maximum_string_bytes":4096,"maximum_native_bytes":8388608,"maximum_forward_sensitivity":1e-10},"models":[{"expansion":{"kind":"constant_q","q":-1.0},"geometry":{"kind":"flat_flrw"},"source_effect":{"kind":"none"}}],"requested_outputs":["score","geometric_shape","magnitude_effect","corrected_residuals","profiled_residuals","diagnostics"],"numerical_policy":{"arithmetic":"wide","projection":{"maximum_queries":2,"maximum_callbacks":10000,"maximum_segment_visits":2000,"integration":{"absolute_tolerance":1e-12,"relative_tolerance":1e-12,"maximum_evaluations":10000,"maximum_depth":30}},"maximum_models":4,"maximum_array_elements":32,"maximum_native_bytes":8388608,"maximum_forward_sensitivity":1e-10}})
}
#[test]
fn current_sn_analytic_measured_source_and_effect_only() {
    let p = root();
    let r = sn(&p);
    let (o, code) = run(&p, &r);
    assert_eq!(code, 0, "{o}");
    assert_eq!(o["receipt"]["execution"], "completed");
    assert_eq!(o["receipt"]["accepted"], true);
    record(&p, &o, &r, "conditional_single_offset_profile");
    // For a common offset, exact 2x2 profile q=(r0-r1)^2/(C00+C11-2C01)=2/9.
    assert!(
        (o["result"]["evaluations"][0]["score"]["value"]["quadratic"]
            .as_f64()
            .unwrap()
            - 2.0 / 9.0)
            .abs()
            < 1e-10
    );
    let serialized = o["result"].to_string();
    assert!(
        !serialized.contains("released_fitted_summary"),
        "measured source must not be upgraded: {o}"
    );
    let mut effect = r.clone();
    effect["requested_outputs"] = json!(["magnitude_effect"]);
    effect["models"][0]["source_effect"] = json!({"kind":"grey_log1p_magnitude","epsilon_mag":0.2});
    effect["numerical_policy"]["projection"]["maximum_callbacks"] = json!(0);
    effect["numerical_policy"]["projection"]["integration"] = Value::Null;
    let (o, code) = run(&p, &effect);
    assert_eq!(code, 0, "{o}");
    let mut invalid = effect;
    invalid["models"][0]["expansion"]["q"] = json!(3.0);
    let (o, code) = run(&p, &invalid);
    assert_eq!(code, 2, "{o}");
    assert_eq!(o["receipt"]["execution"], "completed");
    let _ = fs::remove_dir_all(p);
}

#[test]
fn current_bao_released_reader_shape_and_tags_preserve_original_bytes() {
    let p = root();
    let mut r = bao();
    let mut mean = String::from("# Generated format fixture; not published observations\n");
    for i in 0..13 {
        let tag = if i == 11 {
            "DH_over_rs"
        } else if i == 12 || i % 2 == 0 {
            "DM_over_rs"
        } else {
            "DH_over_rs"
        };
        let z = if i >= 11 { 2.33 } else { 0.1 + i as f64 * 0.02 };
        mean += &format!("{z} {} {tag}\n", 10 + i);
    }
    let rows: Vec<String> = (0..13)
        .map(|i| {
            (0..13)
                .map(|j| if i == j { "1" } else { "0" })
                .collect::<Vec<_>>()
                .join(" ")
        })
        .collect();
    let covariance = rows.join("\n") + "\n";
    fs::write(p.join("mean.txt"), &mean).unwrap();
    fs::write(p.join("cov.txt"), &covariance).unwrap();
    let old = r["source"].clone();
    r["source"] = json!({"profile":"desi_dr2_all_gccomb_13_v1","mean":p.join("mean.txt"),"covariance":p.join("cov.txt"),"maximum_asset_bytes":10000});
    for key in [
        "ordering_provenance",
        "calibration_provenance",
        "dependence_provenance",
        "redshift_convention",
        "ruler_convention",
    ] {
        r["source"][key] = old[key].clone();
    }
    let (o, code) = run(&p, &r);
    assert_eq!(code, 0, "{o}");
    assert_eq!(fs::read(p.join("mean.txt")).unwrap(), mean.as_bytes());
    assert_eq!(fs::read(p.join("cov.txt")).unwrap(), covariance.as_bytes());
    let mut tokens: Vec<Vec<&str>> = covariance
        .lines()
        .map(|l| l.split_whitespace().collect())
        .collect();
    let last = tokens[0].pop().unwrap();
    tokens[1].push(last);
    assert_eq!(tokens.iter().map(Vec::len).sum::<usize>(), 169);
    fs::write(
        p.join("cov.txt"),
        tokens
            .iter()
            .map(|l| l.join(" "))
            .collect::<Vec<_>>()
            .join("\n"),
    )
    .unwrap();
    let (o, code) = run(&p, &r);
    assert_eq!(code, 2, "{o}");
    assert_eq!(o["result"]["error_id"], "INVALID_BAO_COVARIANCE_SHAPE");
    assert_eq!(o["receipt"]["execution"], "failed");
    fs::write(p.join("cov.txt"), covariance).unwrap();
    fs::write(
        p.join("mean.txt"),
        mean.replacen("DM_over_rs", "UNKNOWN", 1),
    )
    .unwrap();
    let (o, code) = run(&p, &r);
    assert_eq!(code, 2, "{o}");
    assert_eq!(o["result"]["error_id"], "UNSUPPORTED_BAO_OBSERVABLE");
    let mut axis = bao();
    axis["source"]["covariance_axis_ids"] = json!(["b", "a"]);
    let (o, code) = run(&p, &axis);
    assert_eq!(code, 2, "{o}");
    let _ = fs::remove_dir_all(p);
}
