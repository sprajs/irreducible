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
        .arg(scratch.0.join("store")).args(["--assurance","qualified"])
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

#[test]
#[ignore = "requires explicit original table/covariance; bounded four-point explicit CPL v2 CLI transport regression"]
fn released_assets_and_cli_cpl_four_points() {
    use serde_json::{Value, json};
    use std::fs;
    let table = PathBuf::from(std::env::var_os("IRRED_W01_TABLE").expect("IRRED_W01_TABLE"));
    let covariance =
        PathBuf::from(std::env::var_os("IRRED_W01_COVARIANCE").expect("IRRED_W01_COVARIANCE"));
    assert_eq!(digest(&table), pinned_hash("w01_table_sha256"));
    assert_eq!(digest(&covariance), pinned_hash("w01_cov_sha256"));
    // Parse frozen reference constants; Rust implements no prediction/profile equation.
    let header = include_str!("../cpp/tests/fixtures/w01_cpl.hpp");
    for (cpl_key, historical_key) in [
        ("w01_cpl_table_sha256", "w01_table_sha256"),
        ("w01_cpl_covariance_sha256", "w01_cov_sha256"),
    ] {
        let marker = format!(" {cpl_key} =");
        let hash = header
            .split_once(&marker)
            .unwrap()
            .1
            .split('"')
            .nth(1)
            .unwrap();
        assert_eq!(
            hash,
            pinned_hash(historical_key),
            "CPL fixture must name the same exact original assets"
        );
    }
    let rows = header
        .split_once("w01_cpl_points{")
        .unwrap()
        .1
        .split_once("};")
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
        models.push(json!({"model":"flat_cpl_late_v1","omega_m":fields[0].parse::<f64>().unwrap(),"constant_q":fields[1].parse::<f64>().unwrap(),"w0":fields[2].parse::<f64>().unwrap(),"wa":fields[3].parse::<f64>().unwrap()}));
        expected.push(fields[4].parse::<f64>().unwrap());
    }
    assert_eq!(models.len(), 4);
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
    let request = json!({"schema_version":1,"operation":"supernova.profile_batch.v2",
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
        .arg(scratch.0.join("store")).args(["--assurance","qualified"])
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
    assert_eq!(rows.len(), 4);
    let mut maximum_error = 0.0_f64;
    for (i, row) in rows.iter().enumerate() {
        assert_eq!(row["identity"]["model_index"], i);
        assert_eq!(row["identity"]["source_parameters"], request["models"][i]);
        assert_eq!(row["identity"]["arithmetic_id"], "F02/longdouble-cpu/v1");
        assert_eq!(row["diagnostics"]["numerical_status"], "ok");
        assert_eq!(row["result"]["kind"], "finite");
        let error = (row["result"]["relative_profile_score"].as_f64().unwrap() - expected[i]).abs();
        assert!(
            error <= 1e-6,
            "point{i} frozen native CPL comparison allocation exceeded: {error}"
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
    // Optional developer capture: exclusive new destination preserves this one
    // acceptance run. Ordinary optional checks still clean their scratch data.
    if let Some(destination) = std::env::var_os("IRRED_W01_CPL_RECORD_DIRECTORY") {
        fn copy_tree(source: &std::path::Path, target: &std::path::Path) {
            fs::create_dir(target).unwrap();
            for entry in fs::read_dir(source).unwrap() {
                let entry = entry.unwrap();
                let to = target.join(entry.file_name());
                if entry.file_type().unwrap().is_dir() {
                    copy_tree(&entry.path(), &to);
                } else {
                    fs::copy(entry.path(), to).unwrap();
                }
            }
        }
        let destination = PathBuf::from(destination);
        fs::create_dir(&destination).expect("acceptance capture destination must be new");
        fs::copy(&request_path, destination.join("request.json")).unwrap();
        fs::write(destination.join("stdout.json"), &output.stdout).unwrap();
        fs::write(destination.join("stderr.txt"), &output.stderr).unwrap();
        copy_tree(&scratch.0.join("store"), &destination.join("store"));
        for key in ["scientific_specification_digest", "output_digest"] {
            let hash = response["receipt"][key].as_str().unwrap();
            assert_eq!(digest(&destination.join("store/objects").join(hash)), hash);
        }
        println!("acceptance_capture={}", destination.display());
    }
    println!(
        "{}",
        json!({"build_id":response["receipt"]["build_id"],"points":4,"selected_rows":1590,
        "maximum_frozen_native_score_error":maximum_error,"fixed_allocation":1e-6,"accepted":false})
    );
}

#[test]
#[ignore = "requires original assets and matching direct-native piecewise transcript harness"]
fn released_assets_and_cli_piecewise_six_points() {
    use serde_json::{Value, json};
    use std::fs;
    let required = |k| PathBuf::from(std::env::var_os(k).unwrap_or_else(|| panic!("required {k}")));
    let table = required("IRRED_W01_TABLE");
    let covariance = required("IRRED_W01_COVARIANCE");
    let native = required("IRRED_W01_PIECEWISE_NATIVE_HARNESS");
    assert_eq!(digest(&table), pinned_hash("w01_table_sha256"));
    assert_eq!(digest(&covariance), pinned_hash("w01_cov_sha256"));
    let n = Command::new(&native)
        .args([covariance.as_os_str(), table.as_os_str()])
        .arg("--verified-original-assets")
        .arg("--direct-native-transcript")
        .output()
        .unwrap();
    assert!(n.status.success(), "{}", String::from_utf8_lossy(&n.stderr));
    let transcript: Value = serde_json::from_slice(&n.stdout).unwrap();
    assert_eq!(transcript["selected_count"], 1590);
    let models: Vec<Value> = transcript["rows"]
        .as_array()
        .unwrap()
        .iter()
        .map(|r| json!({"q":r["q"].as_array().unwrap().iter().map(|x|x.as_f64().unwrap()).collect::<Vec<_>>()}))
        .collect();
    assert_eq!(models.len(), 6);
    let frozen = [
        [0.; 5],
        [-1.; 5],
        [0.5; 5],
        [-0.4, -0.4, -0.2, 0.1, 0.3],
        [-1., 0., -1., 0., -1.],
        [-3., 2., -3., 2., -3.],
    ];
    for (i, q) in frozen.iter().enumerate() {
        for (j, value) in q.iter().enumerate() {
            assert_eq!(
                models[i]["q"][j].as_f64().unwrap().to_bits(),
                (*value as f64).to_bits()
            );
        }
    }

    let path = std::env::temp_dir().join(format!(
        "irred-piecewise-original-{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    fs::create_dir(&path).unwrap();
    let request = json!({"schema_version":1,"operation":"supernova.piecewise_profile_batch.v1","observations":{"schema_version":1,"operation":"observations.prepare.v1","table":table,"uncertainty":covariance,"metadata":{"profile":"pantheon_plus_released_v1","role":"released_fitted_summary","unit":"magnitude","calibration":"unknown","uncertainty":"covariance","uncertainty_unit":"magnitude_squared","component":"total"},"selection":"pantheon_zhd_gt_001","ordering_provenance":"supplied original release order; unlabeled axes declared","calibration_provenance":"unknown calibration; free offset","dependence_provenance":"released covariance; repeated events retained","resources":{"maximum_asset_bytes":268435456,"maximum_rows":4096,"maximum_matrix_elements":16777216,"maximum_string_bytes":16777216}},"models":models,"policy":{"arithmetic":"longdouble_cpu_v1","include_residual_arrays":false,"maximum_models":6,"maximum_source_rows":4096,"maximum_matrix_elements":16777216,"maximum_queries":4096,"maximum_array_elements":1000000,"maximum_native_output_bytes":536870912,"maximum_total_segment_visits":50000,"maximum_forward_sensitivity":1e-10}});
    let exe = PathBuf::from(env!("CARGO_BIN_EXE_irred"));
    let exe_hash = digest(&exe);
    let discovery = Command::new(&exe)
        .args(["describe", "--json"])
        .output()
        .unwrap();
    assert!(discovery.status.success());
    fs::copy(&exe, path.join("executed-irred")).unwrap();
    fs::write(path.join("discovery-before.json"), &discovery.stdout).unwrap();
    fs::write(path.join("native-transcript.json"), &n.stdout).unwrap();
    fs::write(
        path.join("request.json"),
        serde_json::to_vec(&request).unwrap(),
    )
    .unwrap();
    let o = Command::new(&exe)
        .arg("run")
        .arg(path.join("request.json"))
        .arg(path.join("store")).args(["--assurance","qualified"])
        .output()
        .unwrap();
    fs::write(path.join("stdout.json"), &o.stdout).unwrap();
    fs::write(path.join("stderr.txt"), &o.stderr).unwrap();
    assert_eq!(o.status.code(), Some(6));
    let v: Value = serde_json::from_slice(&o.stdout).unwrap();
    let c = &v["result"]["calculation"];
    let obs = &v["result"]["observations"];
    assert_eq!(c["source"]["source_row_count"], 1701);
    assert_eq!(c["source"]["ordered_ids"].as_array().unwrap().len(), 1590);
    assert_eq!(c["segment_visits"], transcript["segment_visits"]);
    let fields = [
        "backward_residual",
        "estimated_forward_sensitivity",
        "coefficient_solve_backward_residual",
        "coefficient_solve_forward_sensitivity",
        "residual_l1",
        "solution_norm_inf",
        "adjusted_residual_l1",
        "adjusted_solution_norm_inf",
    ];
    for i in 0..6 {
        let r = &c["rows"][i];
        let t = &transcript["rows"][i];
        for j in 0..5 {
            assert_eq!(
                r["identity"]["source_parameters"]["q"][j]
                    .as_f64()
                    .unwrap()
                    .to_bits(),
                request["models"][i]["q"][j].as_f64().unwrap().to_bits()
            );
        }
        assert_eq!(r["result"]["kind"], "finite");
        for (a, b) in [
            ("relative_profile_score", "score"),
            ("quadratic", "quadratic"),
            ("offset_coefficient", "offset"),
        ] {
            assert_eq!(
                r["result"][a].as_f64().unwrap().to_bits(),
                t[b].as_f64().unwrap().to_bits()
            );
        }
        for (j, key) in fields.iter().enumerate() {
            assert_eq!(
                r["diagnostics"][key].as_f64().unwrap().to_bits(),
                t["diagnostics"][j].as_f64().unwrap().to_bits()
            );
        }
        assert_eq!(r["diagnostics"]["segment_visits"], t["segment_visits"]);
        assert_eq!(r["diagnostics"]["numerical_status"], "ok");
        assert_eq!(r["result"]["density_status"], "not_applicable");
    }
    let objects = path.join("store/objects");
    for (p, k) in [
        (&table, "w01_table_sha256"),
        (&covariance, "w01_cov_sha256"),
    ] {
        let h = pinned_hash(k);
        assert_eq!(digest(p), h);
        assert_eq!(fs::read(objects.join(&h)).unwrap(), fs::read(p).unwrap());
    }
    let ids = obs["measurement_ids"].as_array().unwrap();
    let indices = c["source"]["selected_source_indices"].as_array().unwrap();
    for (i, index) in indices.iter().enumerate() {
        assert_eq!(
            c["source"]["ordered_ids"][i],
            ids[index.as_u64().unwrap() as usize]
        );
    }
    let selected_expected: Vec<usize> = obs["zhd"]
        .as_array()
        .unwrap()
        .iter()
        .enumerate()
        .filter_map(|(i, z)| (z.as_f64().unwrap() > 0.01).then_some(i))
        .collect();
    assert_eq!(
        indices
            .iter()
            .map(|x| x.as_u64().unwrap() as usize)
            .collect::<Vec<_>>(),
        selected_expected
    );
    for (i, index) in indices.iter().enumerate() {
        let j = index.as_u64().unwrap() as usize;
        assert_eq!(
            c["source"]["z_expansion"][i].as_f64().unwrap().to_bits(),
            obs["zhd"][j].as_f64().unwrap().to_bits()
        );
        assert_eq!(
            c["source"]["z_observer"][i].as_f64().unwrap().to_bits(),
            obs["zhel"][j].as_f64().unwrap().to_bits()
        );
    }
    let mh = obs["full_uncertainty_matrix"]["object_digest"]
        .as_str()
        .unwrap();
    assert_eq!(digest(&objects.join(mh)), mh);
    let matrix: Value = serde_json::from_slice(&fs::read(objects.join(mh)).unwrap()).unwrap();
    let values = matrix["values"].as_array().unwrap();
    assert_eq!(values.len(), 1701 * 1701);
    let mut selected = Sha256::new();
    for a in indices {
        for b in indices {
            selected.update(
                values[a.as_u64().unwrap() as usize * 1701 + b.as_u64().unwrap() as usize]
                    .as_f64()
                    .unwrap()
                    .to_le_bytes(),
            );
        }
    }
    assert_eq!(
        format!("{:x}", selected.finalize()),
        transcript["expected_selected_f64le_sha256"]
            .as_str()
            .unwrap()
    );
    assert_eq!(v["receipt"]["accepted"], false);
    assert_eq!(v["receipt"]["execution"], "completed");
    assert_eq!(v["receipt"]["outputs"][0]["numerical"], "checks_passed");
    assert_eq!(c["interpretation_status"], "unqualified");
    for key in ["scientific_specification_digest", "output_digest"] {
        let h = v["receipt"][key].as_str().unwrap();
        assert_eq!(digest(&objects.join(h)), h);
    }
    assert_eq!(v["receipt"]["executable_digest"], exe_hash);
    let described: Value = serde_json::from_slice(&discovery.stdout).unwrap();
    assert_eq!(v["receipt"]["build_id"], described["build"]["build_id"]);
    assert_eq!(digest(&exe), exe_hash);
    assert_eq!(digest(&path.join("executed-irred")), exe_hash);
    println!(
        "{}",
        json!({"capture":path,"native_sha256":digest(&native),"cli_sha256":exe_hash,"build_id":v["receipt"]["build_id"],"points":6,"exact_transport":true,"independent_science":"separate native reference qualification"})
    );
}
