//! Optional original-input check. No equations or external data in ordinary CI.
//! Frozen hashes/history are owned by the permanent native fixture header.
use sha2::{Digest, Sha256};
use std::{fs::File, io::Read, path::PathBuf, process::Command};
fn pinned_hash(label: &str) -> String {
    let source = include_str!("../cpp/tests/fixtures/w01_historical.hpp");
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
fn released_assets_and_native_historical_comparison() {
    let required =
        |key| PathBuf::from(std::env::var_os(key).unwrap_or_else(|| panic!("required {key}")));
    let table = required("IRRED_W01_TABLE");
    let covariance = required("IRRED_W01_COVARIANCE");
    let executable = required("IRRED_W01_NATIVE_HARNESS");
    assert_eq!(
        digest(&table),
        pinned_hash("w01_table_sha256"),
        "table identity mismatch; native comparison not executed"
    );
    assert_eq!(
        digest(&covariance),
        pinned_hash("w01_cov_sha256"),
        "covariance identity mismatch; native comparison not executed"
    );
    println!("native_harness_sha256={}", digest(&executable));
    let output = Command::new(&executable)
        .arg(&covariance)
        .arg(&table)
        .arg("--verified-original-assets")
        .output()
        .expect("execute prebuilt native comparison");
    assert_eq!(
        digest(&table),
        pinned_hash("w01_table_sha256"),
        "table changed during native comparison"
    );
    assert_eq!(
        digest(&covariance),
        pinned_hash("w01_cov_sha256"),
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

#[test]
#[ignore = "requires explicit original table/covariance; bounded seven-point CLI transport regression"]
fn released_assets_and_cli_seven_points() {
    use serde_json::{Value, json};
    use std::fs;
    let table = PathBuf::from(std::env::var_os("IRRED_W01_TABLE").expect("IRRED_W01_TABLE"));
    let covariance =
        PathBuf::from(std::env::var_os("IRRED_W01_COVARIANCE").expect("IRRED_W01_COVARIANCE"));
    assert_eq!(digest(&table), pinned_hash("w01_table_sha256"));
    assert_eq!(digest(&covariance), pinned_hash("w01_cov_sha256"));
    // Parse frozen reference constants; Rust implements no prediction/profile equation.
    let header = include_str!("../cpp/tests/fixtures/w01_historical.hpp");
    let rows = header
        .split_once("w01_historical{{")
        .unwrap()
        .1
        .split_once("}};")
        .unwrap()
        .0;
    let mut models = Vec::new();
    let mut expected = Vec::new();
    for row in rows.split("},") {
        let fields: Vec<_> = row
            .trim_matches(|c: char| c.is_whitespace() || c == '{' || c == '}')
            .split(',')
            .map(str::trim)
            .filter(|x| !x.is_empty())
            .collect();
        if fields.is_empty() {
            continue;
        }
        assert_eq!(fields.len(), 5);
        let lcdm = fields[0] == "true";
        assert!(lcdm || fields[0] == "false");
        let parameter: f64 = fields[1].parse().unwrap();
        models.push(
            json!({"model":if lcdm {"flat_lcdm_late_v1"}else{"constant_q_flat_v1"},
            "omega_m":if lcdm {parameter}else{0.0},"constant_q":if lcdm {0.0}else{parameter}}),
        );
        expected.push(fields[3].parse::<f64>().unwrap());
    }
    assert_eq!(models.len(), 7);
    let path = std::env::temp_dir().join(format!(
        "irred-original-cli-{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    fs::create_dir(&path).unwrap();
    struct Scratch(PathBuf);
    impl Drop for Scratch {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.0);
        }
    }
    let scratch = Scratch(path);
    let request = json!({"schema_version":1,"operation":"supernova.profile_batch.v1",
        "observations":{"schema_version":1,"operation":"observations.prepare.v1","table":table,"uncertainty":covariance,
        "metadata":{"profile":"pantheon_plus_released_v1","role":"released_fitted_summary","unit":"magnitude","calibration":"unknown",
        "uncertainty":"covariance","uncertainty_unit":"magnitude_squared","component":"total"},"selection":"pantheon_zhd_gt_001",
        "ordering_provenance":"supplied original release row order; unlabeled matrix linkage not independently verified",
        "calibration_provenance":"released fitted summaries; absolute calibration unknown and free offset profiled",
        "dependence_provenance":"released covariance; repeated events retained; no extra independence assumed",
        "resources":{"maximum_asset_bytes":268435456,"maximum_rows":4096,"maximum_matrix_elements":16777216,"maximum_string_bytes":16777216}},
        "models":models,"policy":{"arithmetic":"longdouble_cpu_v1","include_residual_arrays":false,"maximum_models":64,"maximum_source_rows":4096,
        "maximum_matrix_elements":16777216,"maximum_array_elements":1000000,"maximum_native_output_bytes":536870912,"maximum_total_evaluations":20000000,
        "maximum_evaluations_per_integral":100000,"maximum_depth":30,"absolute_tolerance":1e-12,"relative_tolerance":1e-12,"maximum_forward_sensitivity":1e-10}});
    let request_path = scratch.0.join("request.json");
    fs::write(&request_path, serde_json::to_vec(&request).unwrap()).unwrap();
    let executable = PathBuf::from(env!("CARGO_BIN_EXE_irred"));
    println!("cli_sha256={}", digest(&executable));
    let output = Command::new(&executable)
        .arg("run")
        .arg(&request_path)
        .arg(scratch.0.join("store"))
        .output()
        .unwrap();
    assert_eq!(digest(&table), pinned_hash("w01_table_sha256"));
    assert_eq!(digest(&covariance), pinned_hash("w01_cov_sha256"));
    assert_eq!(
        output.status.code(),
        Some(6),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let response: Value = serde_json::from_slice(&output.stdout).unwrap();
    let c = &response["result"]["calculation"];
    assert_eq!(c["kind"], "finite");
    assert_eq!(c["source"]["source_row_count"], 1701);
    assert_eq!(c["source"]["ordered_ids"].as_array().unwrap().len(), 1590);
    assert_eq!(c["source"]["arithmetic_id"], "F02/longdouble-cpu/v1");
    assert_eq!(response["receipt"]["accepted"], false);
    assert_eq!(response["receipt"]["execution"], "completed");
    assert_eq!(
        response["receipt"]["outputs"][0]["numerical"],
        "checks_passed"
    );
    assert_eq!(c["interpretation_status"], "unqualified");
    let rows = c["rows"].as_array().unwrap();
    assert_eq!(rows.len(), 7);
    let mut maximum_error = 0.0_f64;
    for (i, row) in rows.iter().enumerate() {
        assert_eq!(row["identity"]["model_index"], i);
        assert_eq!(row["identity"]["arithmetic_id"], "F02/longdouble-cpu/v1");
        assert_eq!(row["diagnostics"]["numerical_status"], "ok");
        assert_eq!(row["result"]["kind"], "finite");
        let error = (row["result"]["relative_profile_score"].as_f64().unwrap() - expected[i]).abs();
        assert!(
            error <= 2e-7,
            "point{i} fixed historical allocation exceeded: {error}"
        );
        maximum_error = maximum_error.max(error);
        assert_eq!(row["result"]["arrays"]["shape_magnitudes"], json!([]));
        assert_eq!(row["result"]["density_status"], "not_applicable");
    }
    assert!(c["evaluations"].as_u64().unwrap() <= 20_000_000);
    let objects = scratch.0.join("store/objects");
    for key in ["w01_table_sha256", "w01_cov_sha256"] {
        let hash = pinned_hash(key);
        assert_eq!(digest(&objects.join(&hash)), hash);
    }
    for key in ["scientific_specification_digest", "output_digest"] {
        let hash = response["receipt"][key].as_str().unwrap();
        assert_eq!(digest(&objects.join(hash)), hash);
    }
    println!(
        "{}",
        json!({"build_id":response["receipt"]["build_id"],"points":7,"selected_rows":1590,
        "maximum_historical_stable_score_error":maximum_error,"fixed_allocation":2e-7,"accepted":false})
    );
}
