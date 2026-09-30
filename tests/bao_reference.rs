//! Optional original-input check. No equations or external data in ordinary CI.
//! Frozen hashes/history are owned by the permanent native fixture header.
use sha2::{Digest, Sha256};
use std::{fs::File, io::Read, path::PathBuf, process::Command};
fn pinned_hash(label: &str) -> String {
    let source = include_str!("../cpp/tests/fixtures/bao_reference.hpp");
    let marker = format!(" {label} =");
    let tail = source
        .split_once(&marker)
        .expect("missing frozen hash declaration")
        .1;
    let hash = tail.split('"').nth(1).expect("missing frozen hash literal");
    assert_eq!(hash.len(), 64);
    assert!(hash.bytes().all(|b| b.is_ascii_hexdigit()));
    hash.to_owned()
}
fn digest(path: &PathBuf) -> String {
    let mut file = File::open(path).expect("open explicit original asset");
    let mut hash = Sha256::new();
    let mut buffer = [0u8; 65536];
    loop {
        let n = file.read(&mut buffer).expect("read original asset");
        if n == 0 {
            break;
        }
        hash.update(&buffer[..n]);
    }
    format!("{:x}", hash.finalize())
}
#[test]
#[ignore = "requires explicit original assets and prebuilt optional native harness"]
fn released_assets_and_native_eleven_point_comparison() {
    let required =
        |key| PathBuf::from(std::env::var_os(key).unwrap_or_else(|| panic!("required {key}")));
    let mean = required("IRRED_BAO_MEAN");
    let covariance = required("IRRED_BAO_COVARIANCE");
    let executable = required("IRRED_BAO_NATIVE_HARNESS");
    assert_eq!(
        digest(&mean),
        pinned_hash("bao_mean_sha256"),
        "mean identity mismatch; native comparison not executed"
    );
    assert_eq!(
        digest(&covariance),
        pinned_hash("bao_cov_sha256"),
        "covariance identity mismatch; native comparison not executed"
    );
    println!("native_harness_sha256={}", digest(&executable));
    let output = Command::new(&executable)
        .arg("--verified-original-assets")
        .arg(&mean)
        .arg(&covariance)
        .output()
        .expect("execute prebuilt native comparison");
    assert_eq!(
        digest(&mean),
        pinned_hash("bao_mean_sha256"),
        "mean changed during native comparison"
    );
    assert_eq!(
        digest(&covariance),
        pinned_hash("bao_cov_sha256"),
        "covariance changed during native comparison"
    );
    println!("{}", String::from_utf8_lossy(&output.stdout));
    eprintln!("{}", String::from_utf8_lossy(&output.stderr));
    assert!(
        output.status.success(),
        "native original-input comparison failed: {:?}",
        output.status
    );
}

// The default serde_json float decoder may round a printed f64 one ULP away.
// Exact output-bit checks use Rust's correctly rounded decimal parser on the
// flat numeric arrays; serde remains responsible for JSON structure/status.
fn printed_arrays(text: &str, key: &str) -> Vec<Vec<f64>> {
    let marker = format!("\"{key}\":[");
    text.split(&marker)
        .skip(1)
        .map(|tail| {
            tail.split_once(']')
                .unwrap()
                .0
                .split(',')
                .filter(|x| !x.trim().is_empty())
                .map(|x| x.trim().parse::<f64>().unwrap())
                .collect()
        })
        .collect()
}
#[test]
#[ignore = "requires exact original BAO assets; bounded eleven-point CLI transport comparison"]
fn released_assets_and_cli_eleven_points() {
    use serde_json::{Value, json};
    use std::fs;
    let mean = PathBuf::from(std::env::var_os("IRRED_BAO_MEAN").expect("IRRED_BAO_MEAN"));
    let covariance =
        PathBuf::from(std::env::var_os("IRRED_BAO_COVARIANCE").expect("IRRED_BAO_COVARIANCE"));
    assert_eq!(digest(&mean), pinned_hash("bao_mean_sha256"));
    assert_eq!(digest(&covariance), pinned_hash("bao_cov_sha256"));
    let header = include_str!("../cpp/tests/fixtures/bao_reference.hpp");
    let body = header
        .split_once("bao_points{{")
        .unwrap()
        .1
        .split_once("}};")
        .unwrap()
        .0;
    let numbers = body
        .split(',')
        .map(|s| s.trim_matches(|c: char| c.is_whitespace() || c == '{' || c == '}'))
        .filter(|s| !s.is_empty())
        .map(|s| s.parse::<f64>().unwrap())
        .collect::<Vec<_>>();
    assert_eq!(numbers.len(), 11 * 22);
    let points = numbers.chunks_exact(22).collect::<Vec<_>>();
    let models=points.iter().map(|p|json!({"model":if p[1]==-1.0&&p[2]==0.0{"flat_lcdm_late_v1"}else{"flat_cpl_late_v1"},"omega_m":p[0],"constant_q":0.0,"w0":p[1],"wa":p[2],"h0_rd_km_s":p[3]})).collect::<Vec<_>>();
    let policy = json!({"background":{"maximum_parameters":64,"maximum_queries":13,"maximum_slots":143,"maximum_native_output_bytes":1048576,"maximum_total_evaluations":20000000,"maximum_evaluations_per_integral":100000,"maximum_depth":30,"absolute_tolerance":1e-14,"relative_tolerance":1e-14},"arithmetic":"longdouble_cpu_v1","include_predictions":true,"include_residuals":true,"maximum_models":64,"maximum_rows":13,"maximum_matrix_elements":169,"maximum_string_bytes":65536,"maximum_native_bytes":1048576,"maximum_array_elements":286,"maximum_native_output_bytes":1048576,"maximum_forward_sensitivity":1e-10});
    let request = json!({"schema_version":1,"operation":"bao.gaussian_batch.v1","source":{"profile":"desi_dr2_all_gccomb_13_v1","mean":mean,"covariance":covariance,"maximum_asset_bytes":1048576,"ordering_provenance":"supplied released original row order; not independently verified from unlabeled covariance","calibration_provenance":"released fitted-distance Gaussian; free empirical ruler","dependence_provenance":"internal covariance supplied; external probe dependence not assessed","redshift_convention":"P01/released-effective-redshift/v1","ruler_convention":"P01/free-H0rd-km-s-no-early-physics/v1","computational_h0_convention":"P01/computational-H0-fixed-70-km-s-Mpc/v1"},"prepare_policy":policy,"evaluation_policy":policy,"models":models});
    let scratch = std::env::temp_dir().join(format!(
        "irred-bao-original-cli-{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    fs::create_dir(&scratch).unwrap();
    struct Scratch(PathBuf);
    impl Drop for Scratch {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.0);
        }
    }
    let scratch = Scratch(scratch);
    let executable = PathBuf::from(env!("CARGO_BIN_EXE_irred"));
    let discovery = Command::new(&executable)
        .args(["describe", "--json"])
        .output()
        .unwrap();
    assert!(discovery.status.success());
    let manifest: Value = serde_json::from_slice(&discovery.stdout).unwrap();
    let destination = std::env::var_os("IRRED_BAO_CLI_RECORD_DIRECTORY").map(PathBuf::from);
    if let Some(d) = &destination {
        fs::create_dir(d).expect("capture directory must be new");
        fs::write(d.join("discovery-before.json"), &discovery.stdout).unwrap();
        fs::copy(&executable, d.join("executed-irred")).unwrap();
    }
    let request_path = scratch.0.join("request.json");
    let request_bytes = serde_json::to_vec(&request).unwrap();
    fs::write(&request_path, &request_bytes).unwrap();
    let output = Command::new(&executable)
        .arg("run")
        .arg(&request_path)
        .arg(scratch.0.join("store")).args(["--assurance","qualified"])
        .output()
        .unwrap();
    assert_eq!(digest(&mean), pinned_hash("bao_mean_sha256"));
    assert_eq!(digest(&covariance), pinned_hash("bao_cov_sha256"));
    assert_eq!(
        output.status.code(),
        Some(6),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let r: Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(r["receipt"]["accepted"], false);
    assert_eq!(r["receipt"]["execution"], "completed");
    assert_eq!(r["receipt"]["outputs"][0]["numerical"], "checks_passed");
    assert_eq!(r["receipt"]["interpretation"], "not_assessed");
    assert_eq!(r["receipt"]["outputs"][0]["interpretation"], "unqualified");
    assert_eq!(r["receipt"]["build_id"], manifest["build"]["build_id"]);
    // Original bytes, not merely a matching profile name, own these row identities.
    let source = &r["result"]["source"];
    let source_rows = source["rows"].as_array().unwrap();
    let original_mean = fs::read_to_string(&mean).unwrap();
    let fields = original_mean
        .lines()
        .filter(|line| !line.trim().is_empty() && !line.trim_start().starts_with('#'))
        .map(|line| line.split_whitespace().collect::<Vec<_>>())
        .collect::<Vec<_>>();
    assert_eq!(fields.len(), 13);
    assert_eq!(source_rows.len(), 13);
    for (j, (row, field)) in source_rows.iter().zip(&fields).enumerate() {
        assert_eq!(field.len(), 3);
        assert_eq!(
            row["z"].as_f64().unwrap().to_bits(),
            field[0].parse::<f64>().unwrap().to_bits()
        );
        assert_eq!(
            row["value"].as_f64().unwrap().to_bits(),
            field[1].parse::<f64>().unwrap().to_bits()
        );
        assert_eq!(row["observable"], field[2]);
        assert_eq!(
            row["id"],
            format!("{}:row:{j}:{}", pinned_hash("bao_mean_sha256"), field[2])
        );
    }
    assert_eq!(source_rows[11]["observable"], "DH_over_rs");
    assert_eq!(source_rows[12]["observable"], "DM_over_rs");
    assert_eq!(source_rows[11]["z"], 2.33);
    assert_eq!(source_rows[12]["z"], 2.33);
    let output_text = std::str::from_utf8(&output.stdout).unwrap();
    let printed_predictions = printed_arrays(output_text, "predictions");
    let printed_residuals = printed_arrays(output_text, "residuals");
    assert_eq!(printed_predictions.len(), 11);
    assert_eq!(printed_residuals.len(), 11);
    let rows = r["result"]["calculation"]["rows"].as_array().unwrap();
    assert_eq!(rows.len(), 11);
    let mut maximum_error = 0.0f64;
    for (i, (row, p)) in rows.iter().zip(&points).enumerate() {
        assert_eq!(row["identity"]["model_index"], i);
        assert_eq!(row["identity"]["source_parameters"], request["models"][i]);
        assert_eq!(row["numerical_status"], "ok");
        assert_eq!(row["status"], "finite");
        let values = &row["result"];
        let error = (values["log_density"].as_f64().unwrap() - p[4]).abs();
        maximum_error = maximum_error.max(error);
        assert!(error <= 1e-8, "point{i} density error {error}");
        assert!((values["quadratic"].as_f64().unwrap() - p[5]).abs() <= 2e-8);
        assert!((values["log_determinant"].as_f64().unwrap() - p[6]).abs() <= 1e-10);
        assert!((values["normalization"].as_f64().unwrap() - p[7]).abs() <= 1e-10);
        let predicted = values["predictions"].as_array().unwrap();
        assert_eq!(predicted.len(), 13);
        assert_eq!(values["residuals"].as_array().unwrap().len(), 13);
        for (j, x) in predicted.iter().enumerate() {
            let expected = p[9 + j];
            assert!((x.as_f64().unwrap() - expected).abs() <= 2e-12 + 2e-10 * expected.abs());
            let observed = source_rows[j]["value"].as_f64().unwrap();
            assert_eq!(
                printed_residuals[i][j].to_bits(),
                (observed - printed_predictions[i][j]).to_bits()
            );
        }
    }
    let objects = scratch.0.join("store/objects");
    assert_eq!(
        fs::read(objects.join(pinned_hash("bao_mean_sha256"))).unwrap(),
        fs::read(&mean).unwrap()
    );
    assert_eq!(
        fs::read(objects.join(pinned_hash("bao_cov_sha256"))).unwrap(),
        fs::read(&covariance).unwrap()
    );
    let matrix_text = fs::read_to_string(
        objects.join(source["full_covariance"]["object_digest"].as_str().unwrap()),
    )
    .unwrap();
    let matrix: Value = serde_json::from_str(&matrix_text).unwrap();
    let covariance_values = fs::read_to_string(&covariance)
        .unwrap()
        .split_whitespace()
        .map(|value| value.parse::<f64>().unwrap())
        .collect::<Vec<_>>();
    assert_eq!(covariance_values.len(), 169);
    let printed_values = printed_arrays(&matrix_text, "values");
    assert_eq!(printed_values.len(), 1);
    for (retained, original) in printed_values[0].iter().zip(&covariance_values) {
        assert_eq!(retained.to_bits(), original.to_bits());
    }
    assert_eq!(matrix["values"].as_array().unwrap().len(), 169);
    assert_eq!(matrix["axis_ids"], source["full_covariance"]["axis_ids"]);
    assert_eq!(
        matrix["axis_ids"],
        json!(
            source_rows
                .iter()
                .map(|row| row["id"].clone())
                .collect::<Vec<_>>()
        )
    );
    for entry in fs::read_dir(&objects).unwrap() {
        let p = entry.unwrap().path();
        assert_eq!(digest(&p), p.file_name().unwrap().to_str().unwrap());
    }
    for key in [
        "scientific_specification_digest",
        "output_digest",
        "input_digest",
    ] {
        let h = r["receipt"][key].as_str().unwrap();
        assert_eq!(digest(&objects.join(h)), h);
    }
    if let Some(d) = &destination {
        fn copy_tree(a: &std::path::Path, b: &std::path::Path) {
            fs::create_dir(b).unwrap();
            for e in fs::read_dir(a).unwrap() {
                let e = e.unwrap();
                let dest = b.join(e.file_name());
                if e.file_type().unwrap().is_dir() {
                    copy_tree(&e.path(), &dest);
                } else {
                    fs::copy(e.path(), dest).unwrap();
                }
            }
        }
        fs::write(d.join("request.json"), &request_bytes).unwrap();
        fs::write(d.join("stdout.json"), &output.stdout).unwrap();
        fs::write(d.join("stderr.txt"), &output.stderr).unwrap();
        copy_tree(&scratch.0.join("store"), &d.join("store"));
        assert_eq!(
            digest(&d.join("executed-irred")),
            r["receipt"]["executable_digest"].as_str().unwrap()
        );
    }
    println!(
        "{}",
        json!({"build_id":r["receipt"]["build_id"],"points":11,"maximum_reference_fixture_density_difference":maximum_error,"assembled_budget":1e-8,"accepted":false})
    );
}

#[test]
#[ignore = "requires exact original BAO assets and matching current native direct24 harness"]
fn released_assets_and_cli_piecewise_twentyfour_points() {
    use serde_json::{Value, json};
    use std::fs;
    let mean = PathBuf::from(std::env::var_os("IRRED_BAO_MEAN").expect("IRRED_BAO_MEAN"));
    let covariance =
        PathBuf::from(std::env::var_os("IRRED_BAO_COVARIANCE").expect("IRRED_BAO_COVARIANCE"));
    assert_eq!(digest(&mean), pinned_hash("bao_mean_sha256"));
    assert_eq!(digest(&covariance), pinned_hash("bao_cov_sha256"));
    let harness = PathBuf::from(
        std::env::var_os("IRRED_BAO_PIECEWISE_NATIVE_HARNESS")
            .expect("IRRED_BAO_PIECEWISE_NATIVE_HARNESS"),
    );
    let harness_hash = digest(&harness);
    let direct = Command::new(&harness)
        .arg("--verified-original-assets")
        .arg(&mean)
        .arg(&covariance)
        .arg("--direct-native-transcript")
        .output()
        .unwrap();
    assert!(
        direct.status.success(),
        "{}",
        String::from_utf8_lossy(&direct.stderr)
    );
    assert_eq!(digest(&mean), pinned_hash("bao_mean_sha256"));
    assert_eq!(digest(&covariance), pinned_hash("bao_cov_sha256"));
    let native: Value = serde_json::from_slice(&direct.stdout).unwrap();
    assert_eq!(native["suite"], "piecewise_BAO_direct_native_transport");
    assert_eq!(native["mean_sha256"], pinned_hash("bao_mean_sha256"));
    assert_eq!(native["covariance_sha256"], pinned_hash("bao_cov_sha256"));
    assert_eq!(native["status"], 0);
    assert_eq!(native["numerical_status"], 0);
    let q_grid: [[f64; 5]; 8] = [
        [0.0; 5],
        [-1.0; 5],
        [0.5; 5],
        [-0.4, -0.4, -0.2, 0.1, 0.3],
        [-1.0, 0.0, -1.0, 0.0, -1.0],
        [-3.0, 2.0, -3.0, 2.0, -3.0],
        [-3.0; 5],
        [2.0; 5],
    ];
    let models = q_grid
        .into_iter()
        .flat_map(|q| [5000.0, 10000.0, 15000.0].map(|r| json!({"q":q,"h0_rd_km_s":r})))
        .collect::<Vec<_>>();
    let points = native["rows"].as_array().unwrap();
    assert_eq!(points.len(), 24);
    for (i, p) in points.iter().enumerate() {
        assert_eq!(p["point"], i);
        for k in 0..5 {
            assert_eq!(
                p["q"][k].as_f64().unwrap().to_bits(),
                models[i]["q"][k].as_f64().unwrap().to_bits()
            );
        }
        assert_eq!(
            p["h0rd"].as_f64().unwrap().to_bits(),
            models[i]["h0_rd_km_s"].as_f64().unwrap().to_bits()
        );
    }
    let policy = json!({"arithmetic":"longdouble_cpu_v1","include_predictions":true,"include_residuals":true,"maximum_models":24,"maximum_rows":13,"maximum_matrix_elements":169,"maximum_string_bytes":65536,"maximum_native_bytes":1048576,"maximum_array_elements":624,"maximum_native_output_bytes":1048576,"maximum_queries":13,"maximum_total_segment_visits":2000,"maximum_forward_sensitivity":1e-10});
    let request = json!({"schema_version":1,"operation":"bao.piecewise_gaussian_batch.v1","source":{"profile":"desi_dr2_all_gccomb_13_v1","mean":mean,"covariance":covariance,"maximum_asset_bytes":1048576,"ordering_provenance":"supplied released original row order; not independently verified from unlabeled covariance","calibration_provenance":"released fitted-distance Gaussian; free empirical ruler","dependence_provenance":"internal covariance supplied; external probe dependence not assessed","redshift_convention":"P01/released-effective-redshift/v1","ruler_convention":"P01/free-H0rd-km-s-no-early-physics/v1","computational_h0_convention":"P01/computational-H0-fixed-70-km-s-Mpc/v1"},"prepare_policy":policy,"evaluation_policy":policy,"models":models});
    let scratch = std::env::temp_dir().join(format!(
        "irred-bao-piecewise-original-cli-{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    fs::create_dir(&scratch).unwrap();
    struct Scratch(PathBuf);
    impl Drop for Scratch {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.0);
        }
    }
    let scratch = Scratch(scratch);
    let executable = PathBuf::from(env!("CARGO_BIN_EXE_irred"));
    let discovery = Command::new(&executable)
        .args(["describe", "--json"])
        .output()
        .unwrap();
    assert!(discovery.status.success());
    let manifest: Value = serde_json::from_slice(&discovery.stdout).unwrap();
    let destination =
        std::env::var_os("IRRED_BAO_PIECEWISE_CLI_RECORD_DIRECTORY").map(PathBuf::from);
    if let Some(d) = &destination {
        fs::create_dir(d).expect("capture directory must be new");
        fs::write(d.join("discovery-before.json"), &discovery.stdout).unwrap();
        fs::copy(&executable, d.join("executed-irred")).unwrap();
        fs::copy(&harness, d.join("executed-native-harness")).unwrap();
        fs::write(d.join("direct-native.json"), &direct.stdout).unwrap();
        fs::write(d.join("native-harness-sha256.txt"), &harness_hash).unwrap();
    }
    let request_path = scratch.0.join("request.json");
    let request_bytes = serde_json::to_vec(&request).unwrap();
    fs::write(&request_path, &request_bytes).unwrap();
    let output = Command::new(&executable)
        .arg("run")
        .arg(&request_path)
        .arg(scratch.0.join("store")).args(["--assurance","qualified"])
        .output()
        .unwrap();
    assert_eq!(digest(&mean), pinned_hash("bao_mean_sha256"));
    assert_eq!(digest(&covariance), pinned_hash("bao_cov_sha256"));
    assert_eq!(
        output.status.code(),
        Some(6),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let r: Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(r["receipt"]["accepted"], false);
    assert_eq!(r["receipt"]["execution"], "completed");
    assert_eq!(r["receipt"]["outputs"][0]["numerical"], "checks_passed");
    assert_eq!(r["receipt"]["interpretation"], "not_assessed");
    assert_eq!(r["receipt"]["outputs"][0]["interpretation"], "unqualified");
    assert_eq!(r["receipt"]["build_id"], manifest["build"]["build_id"]);
    // Original bytes, not merely a matching profile name, own these row identities.
    let source = &r["result"]["source"];
    let source_rows = source["rows"].as_array().unwrap();
    let original_mean = fs::read_to_string(&mean).unwrap();
    let fields = original_mean
        .lines()
        .filter(|line| !line.trim().is_empty() && !line.trim_start().starts_with('#'))
        .map(|line| line.split_whitespace().collect::<Vec<_>>())
        .collect::<Vec<_>>();
    assert_eq!(fields.len(), 13);
    assert_eq!(source_rows.len(), 13);
    for (j, (row, field)) in source_rows.iter().zip(&fields).enumerate() {
        assert_eq!(field.len(), 3);
        assert_eq!(
            row["z"].as_f64().unwrap().to_bits(),
            field[0].parse::<f64>().unwrap().to_bits()
        );
        assert_eq!(
            row["value"].as_f64().unwrap().to_bits(),
            field[1].parse::<f64>().unwrap().to_bits()
        );
        assert_eq!(row["observable"], field[2]);
        assert_eq!(
            row["id"],
            format!("{}:row:{j}:{}", pinned_hash("bao_mean_sha256"), field[2])
        );
    }
    assert_eq!(source_rows[11]["observable"], "DH_over_rs");
    assert_eq!(source_rows[12]["observable"], "DM_over_rs");
    assert_eq!(source_rows[11]["z"], 2.33);
    assert_eq!(source_rows[12]["z"], 2.33);
    let native_source = &native["source"];
    assert_eq!(native_source["role"], 0);
    assert_eq!(native_source["covariance_unit"], 0);
    assert_eq!(
        native_source["table_identity"],
        source["table_asset_sha256"]
    );
    assert_eq!(
        native_source["covariance_identity"],
        source["covariance_asset_sha256"]
    );
    for key in [
        "ordering_provenance",
        "calibration_provenance",
        "dependence_provenance",
    ] {
        assert_eq!(native_source[key], source[key]);
    }
    assert_eq!(native_source["ordered_ids"].as_array().unwrap().len(), 13);
    assert_eq!(native_source["queries"].as_array().unwrap().len(), 13);
    assert_eq!(native_source["observed"].as_array().unwrap().len(), 13);
    for j in 0..13 {
        assert_eq!(native_source["ordered_ids"][j], source_rows[j]["id"]);
        assert_eq!(
            native_source["queries"][j]["z"].as_f64().unwrap().to_bits(),
            source_rows[j]["z"].as_f64().unwrap().to_bits()
        );
        assert_eq!(
            native_source["observed"][j].as_f64().unwrap().to_bits(),
            source_rows[j]["value"].as_f64().unwrap().to_bits()
        );
        let tag = match source_rows[j]["observable"].as_str().unwrap() {
            "DM_over_rs" => 0,
            "DH_over_rs" => 1,
            "DV_over_rs" => 2,
            _ => panic!("unknown original observable"),
        };
        assert_eq!(native_source["queries"][j]["observable"], tag);
    }
    let output_text = std::str::from_utf8(&output.stdout).unwrap();
    let printed_predictions = printed_arrays(output_text, "predictions");
    let printed_residuals = printed_arrays(output_text, "residuals");
    assert_eq!(printed_predictions.len(), 24);
    assert_eq!(printed_residuals.len(), 24);
    let rows = r["result"]["calculation"]["rows"].as_array().unwrap();
    assert_eq!(rows.len(), 24);
    assert_eq!(
        r["result"]["calculation"]["segment_visits"],
        native["segments"]
    );
    let mut maximum_error = 0.0f64;
    for (i, (row, p)) in rows.iter().zip(points).enumerate() {
        assert_eq!(row["identity"]["model_index"], i);
        for k in 0..5 {
            assert_eq!(
                row["identity"]["source_parameters"]["q"][k]
                    .as_f64()
                    .unwrap()
                    .to_bits(),
                request["models"][i]["q"][k].as_f64().unwrap().to_bits()
            );
        }
        assert_eq!(
            row["identity"]["source_parameters"]["h0_rd_km_s"]
                .as_f64()
                .unwrap()
                .to_bits(),
            request["models"][i]["h0_rd_km_s"]
                .as_f64()
                .unwrap()
                .to_bits()
        );
        assert_eq!(row["identity"]["arithmetic_id"], "F02/longdouble-cpu/v1");
        assert_eq!(row["numerical_status"], "ok");
        assert_eq!(row["status"], "finite");
        assert_eq!(p["numerical_status"], 0);
        assert_eq!(p["density_status"], 0);
        assert_eq!(row["segment_visits"], p["segments"]);
        let values = &row["result"];
        for (cli, direct) in [
            ("log_density", "density"),
            ("quadratic", "quadratic"),
            ("log_determinant", "logdet"),
            ("normalization", "normalization"),
            ("backward_residual", "backward_residual"),
            ("estimated_forward_sensitivity", "forward_sensitivity"),
        ] {
            let a = values[cli].as_f64().unwrap();
            let b = p[direct].as_f64().unwrap();
            maximum_error = maximum_error.max((a - b).abs());
            assert_eq!(a.to_bits(), b.to_bits(), "point{i} {cli}");
        }
        for j in 0..13 {
            assert_eq!(
                printed_predictions[i][j].to_bits(),
                p["prediction"][j].as_f64().unwrap().to_bits()
            );
            assert_eq!(
                printed_residuals[i][j].to_bits(),
                p["residual"][j].as_f64().unwrap().to_bits()
            );
        }
        assert!(values.get("relative_profile_score").is_none());
    }
    let objects = scratch.0.join("store/objects");
    assert_eq!(
        fs::read(objects.join(pinned_hash("bao_mean_sha256"))).unwrap(),
        fs::read(&mean).unwrap()
    );
    assert_eq!(
        fs::read(objects.join(pinned_hash("bao_cov_sha256"))).unwrap(),
        fs::read(&covariance).unwrap()
    );
    let matrix_text = fs::read_to_string(
        objects.join(source["full_covariance"]["object_digest"].as_str().unwrap()),
    )
    .unwrap();
    let matrix: Value = serde_json::from_str(&matrix_text).unwrap();
    let covariance_values = fs::read_to_string(&covariance)
        .unwrap()
        .split_whitespace()
        .map(|value| value.parse::<f64>().unwrap())
        .collect::<Vec<_>>();
    assert_eq!(covariance_values.len(), 169);
    assert_eq!(native_source["covariance"].as_array().unwrap().len(), 169);
    for (j, expected) in covariance_values.iter().enumerate() {
        assert_eq!(
            native_source["covariance"][j].as_f64().unwrap().to_bits(),
            expected.to_bits()
        );
    }

    let printed_values = printed_arrays(&matrix_text, "values");
    assert_eq!(printed_values.len(), 1);
    for (retained, original) in printed_values[0].iter().zip(&covariance_values) {
        assert_eq!(retained.to_bits(), original.to_bits());
    }
    assert_eq!(matrix["values"].as_array().unwrap().len(), 169);
    assert_eq!(matrix["axis_ids"], source["full_covariance"]["axis_ids"]);
    assert_eq!(
        matrix["axis_ids"],
        json!(
            source_rows
                .iter()
                .map(|row| row["id"].clone())
                .collect::<Vec<_>>()
        )
    );
    for entry in fs::read_dir(&objects).unwrap() {
        let p = entry.unwrap().path();
        assert_eq!(digest(&p), p.file_name().unwrap().to_str().unwrap());
    }
    for key in [
        "scientific_specification_digest",
        "output_digest",
        "input_digest",
    ] {
        let h = r["receipt"][key].as_str().unwrap();
        assert_eq!(digest(&objects.join(h)), h);
    }
    if let Some(d) = &destination {
        fn copy_tree(a: &std::path::Path, b: &std::path::Path) {
            fs::create_dir(b).unwrap();
            for e in fs::read_dir(a).unwrap() {
                let e = e.unwrap();
                let dest = b.join(e.file_name());
                if e.file_type().unwrap().is_dir() {
                    copy_tree(&e.path(), &dest);
                } else {
                    fs::copy(e.path(), dest).unwrap();
                }
            }
        }
        fs::write(d.join("request.json"), &request_bytes).unwrap();
        fs::write(d.join("stdout.json"), &output.stdout).unwrap();
        fs::write(d.join("stderr.txt"), &output.stderr).unwrap();
        copy_tree(&scratch.0.join("store"), &d.join("store"));
        assert_eq!(
            digest(&d.join("executed-irred")),
            r["receipt"]["executable_digest"].as_str().unwrap()
        );
    }
    assert_eq!(digest(&harness), harness_hash);
    println!(
        "{}",
        json!({"build_id":r["receipt"]["build_id"],"points":24,"maximum_direct_native_difference":maximum_error,"ancestry":"exact transport parity; independent scientific reference is directed interval gate","accepted":false})
    );
}
