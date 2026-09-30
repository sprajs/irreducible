//! Optional minimal current original-asset ingestion/transport guards.
//! Exact production-native transcript agreement is not independent science.
//! The full retired original CLI campaign is not repeated here.
use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use std::{
    fs,
    path::PathBuf,
    process::Command,
    time::{SystemTime, UNIX_EPOCH},
};
fn digest(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}
fn path(key: &str) -> PathBuf {
    std::env::var_os(key)
        .map(PathBuf::from)
        .unwrap_or_else(|| panic!("required {key}"))
}
fn object(store: &PathBuf, d: &Value) -> Value {
    let d = d.as_str().unwrap();
    let b = fs::read(store.join("objects").join(d)).unwrap();
    assert_eq!(digest(&b), d);
    serde_json::from_slice(&b).unwrap()
}
fn run(request: &Value) -> (Value, PathBuf) {
    let p = std::env::temp_dir().join(format!(
        "irred-original-{}-{}",
        std::process::id(),
        SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    fs::create_dir_all(&p).unwrap();
    let bytes = serde_json::to_vec(request).unwrap();
    fs::write(p.join("request.json"), &bytes).unwrap();
    let executable = PathBuf::from(env!("CARGO_BIN_EXE_irred"));
    let executable_hash = digest(&fs::read(&executable).unwrap());
    let discovery = Command::new(&executable)
        .args(["describe", "--json"])
        .output()
        .unwrap();
    assert!(discovery.status.success());
    fs::write(p.join("discovery-before.json"), &discovery.stdout).unwrap();
    let o = Command::new(&executable)
        .arg("run")
        .arg(p.join("request.json"))
        .arg(p.join("store"))
        .args(["--assurance", "qualified"])
        .output()
        .unwrap();
    fs::write(p.join("stdout.json"), &o.stdout).unwrap();
    let r: Value = serde_json::from_slice(&o.stdout)
        .unwrap_or_else(|_| panic!("{}", String::from_utf8_lossy(&o.stderr)));
    assert_eq!(o.status.code(), Some(6), "{r}");
    let q = &r["receipt"];
    assert_eq!(q["execution"], "completed");
    assert_eq!(q["numerical"], "checks_passed");
    assert_eq!(q["accepted"], false);
    assert_eq!(q["assurance"]["numerical_contract"], true);
    assert_eq!(q["interpretation"], "unqualified");
    assert_eq!(q["executable_digest"], executable_hash);
    assert_eq!(q["input_digest"], digest(&bytes));
    assert_eq!(
        fs::read(p.join("store/objects").join(digest(&bytes))).unwrap(),
        bytes
    );
    assert_eq!(object(&p.join("store"), &q["output_digest"]), r["result"]);
    println!("capture={}", p.display());
    (r, p)
}
fn preparation(rows: usize, matrix: usize) -> Value {
    json!({"arithmetic":"wide","maximum_selected_rows":rows,"maximum_matrix_elements":matrix,"maximum_string_bytes":1048576,"maximum_native_bytes":536870912,"maximum_forward_sensitivity":1e-10})
}
fn evaluation(rows: usize) -> Value {
    json!({"arithmetic":"wide","projection":{"maximum_queries":rows,"maximum_callbacks":20000000,"maximum_segment_visits":50000,"integration":{"absolute_tolerance":1e-12,"relative_tolerance":1e-12,"maximum_evaluations":100000,"maximum_depth":30}},"maximum_models":1,"maximum_array_elements":4*rows,"maximum_native_bytes":536870912,"maximum_forward_sensitivity":1e-10})
}
#[test]
#[ignore = "requires exact original assets and matching recorded production-native stdout"]
fn original_pantheon_reader_one_lcdm_point() {
    let table = path("IRRED_W01_TABLE");
    let covariance = path("IRRED_W01_COVARIANCE");
    let native_log = path("IRRED_SN_NATIVE_OUTPUT");
    let native_harness = path("IRRED_W01_NATIVE_HARNESS");
    let tb = fs::read(&table).unwrap();
    let cb = fs::read(&covariance).unwrap();
    let th = digest(&tb);
    let ch = digest(&cb);
    assert_eq!(
        th,
        "1cb0fc379ef066afdc2ffd1857681cc478024570d8a3eba284fb645775198cf8"
    );
    assert_eq!(
        ch,
        "abf806d966485e64afdb359c87bffc0ecc00d05eff0a31ced66f247385df0fdc"
    );
    let log = fs::read_to_string(native_log).unwrap();
    assert!(log.contains(&format!(
        "native_harness_sha256={}",
        digest(&fs::read(native_harness).unwrap())
    )));
    let expected: Value = log
        .lines()
        .filter_map(|l| serde_json::from_str::<Value>(l).ok())
        .find(|v| v["point"] == 1 && v["epsilon_mag"] == 0)
        .unwrap();
    assert_eq!(
        expected["omega_m"].as_f64().unwrap().to_bits(),
        0.3f64.to_bits()
    );
    assert_eq!(expected["passed"], true);
    let lines: Vec<_> = std::str::from_utf8(&tb)
        .unwrap()
        .lines()
        .filter(|l| !l.trim().is_empty() && !l.starts_with('#'))
        .collect();
    let header: Vec<_> = lines[0].split_whitespace().collect();
    let hd = header.iter().position(|s| *s == "zHD").unwrap();
    let hel = header.iter().position(|s| *s == "zHEL").unwrap();
    let rows: Vec<Vec<_>> = lines[1..]
        .iter()
        .map(|l| l.split_whitespace().collect())
        .collect();
    assert_eq!(rows.len(), 1701);
    let indices: Vec<_> = rows
        .iter()
        .enumerate()
        .filter(|(_, r)| r[hd].parse::<f64>().unwrap() > 0.01)
        .map(|(i, _)| i)
        .collect();
    assert_eq!(indices.len(), 1590);
    let request = json!({"schema_version":2,"operation":"supernova.profile","observations":{"schema_version":2,"operation":"observations.prepare","table":table,"uncertainty":covariance,"metadata":{"profile":"pantheon_plus_released_v1","role":"released_fitted_summary","unit":"magnitude","calibration":"released_corrected","uncertainty":"covariance","uncertainty_unit":"magnitude_squared","component":"total"},"resources":{"maximum_asset_bytes":134217728,"maximum_rows":1701,"maximum_matrix_elements":2893401,"maximum_string_bytes":1048576,"maximum_preparation_bytes":536870912},"exports":[],"source_selection":vec![1u8;1701],"calibration_provenance":"named release corrected magnitude declaration; no new calibration","dependence_provenance":"supplied STAT+SYS covariance; cross-probe dependence unknown","quality_dictionary":"not supplied","ordering_provenance":"supplied original file order; unlabeled covariance axes not independently verified"},"selection":{"kind":"pantheon_zhd_gt001"},"preparation_policy":preparation(1590,2528100),"models":[{"expansion":{"kind":"lcdm","omega_m":0.3},"geometry":{"kind":"flat_flrw"},"source_effect":{"kind":"none"}}],"requested_outputs":["score"],"numerical_policy":evaluation(1590)});
    let (r, p) = run(&request);
    let spec = object(
        &p.join("store"),
        &r["receipt"]["scientific_specification_digest"],
    );
    let prep = &spec["preparation"]["preparation"];
    let src = &prep["source_specification"];
    assert_eq!(src["table_identity"], th);
    assert_eq!(src["covariance_identity"], ch);
    assert_eq!(src["ordered_ids"].as_array().unwrap().len(), 1701);
    assert_eq!(src["ordered_ids"], src["covariance_axis_ids"]);
    let selection = &prep["selection"];
    assert_eq!(selection["source_indices"], json!(indices));
    for (j, &i) in indices.iter().enumerate() {
        assert_eq!(selection["ordered_ids"][j], format!("{th}:row:{i}"));
        assert_eq!(
            selection["coordinates"][j]["z_expansion"]
                .as_f64()
                .unwrap()
                .to_bits(),
            rows[i][hd].parse::<f64>().unwrap().to_bits()
        );
        assert_eq!(
            selection["coordinates"][j]["observer"]["redshift"]
                .as_f64()
                .unwrap()
                .to_bits(),
            rows[i][hel].parse::<f64>().unwrap().to_bits()
        );
        assert_eq!(selection["coordinates"][j]["observer"]["convention"], 1);
    }
    let matrix: Vec<f64> = std::str::from_utf8(&cb)
        .unwrap()
        .split_whitespace()
        .skip(1)
        .map(|s| s.parse().unwrap())
        .collect();
    assert_eq!(matrix.len(), 1701 * 1701);
    let mut full = Sha256::new();
    for x in &matrix {
        full.update(x.to_le_bytes());
    }
    assert_eq!(
        src["covariance_f64le_digest"],
        format!("{:x}", full.finalize())
    );
    let mut selected = Sha256::new();
    for &i in &indices {
        for &j in &indices {
            selected.update(matrix[i * 1701 + j].to_le_bytes());
        }
    }
    assert_eq!(
        format!("{:x}", selected.finalize()),
        "64339b81cfd28998ca04f589394c9209b44c93dd63112704eb21bdeb6f300ee4"
    );
    assert_eq!(
        r["result"]["evaluations"][0]["score"]["value"]["relative_profile_score"]
            .as_f64()
            .unwrap()
            .to_bits(),
        expected["relative_score"].as_f64().unwrap().to_bits()
    );
    for (raw, hash) in [(&tb, &th), (&cb, &ch)] {
        assert_eq!(fs::read(p.join("store/objects").join(hash)).unwrap(), *raw);
    }
    assert_eq!(digest(&fs::read(&table).unwrap()), th);
    assert_eq!(digest(&fs::read(&covariance).unwrap()), ch);
    // Native reference driver uses CID/survey/original-row IDs, while the current
    // reader uses raw-SHA row IDs. Exact score comparison does not claim identical
    // metadata objects. Numeric rows/selection/covariance/z roles are audited above.
}
#[test]
#[ignore = "requires exact original assets and matching recorded production-native stdout"]
fn original_desi_reader_one_normalized_point() {
    let mean = path("IRRED_BAO_MEAN");
    let covariance = path("IRRED_BAO_COVARIANCE");
    let mb = fs::read(&mean).unwrap();
    let cb = fs::read(&covariance).unwrap();
    let mh = digest(&mb);
    let ch = digest(&cb);
    assert_eq!(
        mh,
        "9ac154ab583ce759c0f7eef3c978c7c70a6ead2d18774caceadf1a350a640585"
    );
    assert_eq!(
        ch,
        "252a143274c8a07c78694c119617d36594f6d7965d00319ca611c6ffb886e509"
    );
    let log = fs::read_to_string(path("IRRED_BAO_NATIVE_OUTPUT")).unwrap();
    assert!(log.contains(&format!(
        "native_harness_sha256={}",
        digest(&fs::read(path("IRRED_BAO_NATIVE_HARNESS")).unwrap())
    )));
    let native: Vec<_> = log
        .lines()
        .map(|l| l.split_whitespace().collect::<Vec<_>>())
        .find(|v| {
            v.first() == Some(&"BAO_ORIGINAL")
                && v[1].parse::<f64>().unwrap() == 1e-12
                && v[2] == "1"
        })
        .unwrap();
    assert_eq!(native[3], "0");
    assert_eq!(native[4], "0");
    let request = json!({"schema_version":2,"operation":"bao.density","source":{"profile":"desi_dr2_all_gccomb_13_v1","mean":mean,"covariance":covariance,"maximum_asset_bytes":10000,"ordering_provenance":"declared original DESI13 release row order; hash-gated caller","calibration_provenance":"released free ruler distance summaries; no early ruler calibration","dependence_provenance":"supplied full BAO covariance; crossprobe dependence unknown","redshift_convention":"P01/released-effective-redshift/v1","ruler_convention":"P01/free-H0rd-km-s-no-early-physics/v1"},"preparation_policy":{"arithmetic":"wide","maximum_queries":13,"maximum_matrix_elements":169,"maximum_string_bytes":10000,"maximum_native_bytes":1000000,"maximum_forward_sensitivity":1e-10},"models":[{"expansion":{"kind":"cpl","omega_m":0.3,"w0":-1.0,"wa":0.0},"geometry":{"kind":"flat_flrw"},"h0_rd_km_s":10000.0}],"requested_outputs":["normalized_density","predictions","residuals"],"numerical_policy":evaluation(13)});
    let (r, p) = run(&request);
    let spec = object(
        &p.join("store"),
        &r["receipt"]["scientific_specification_digest"],
    );
    let src = &spec["preparation"]["source"];
    assert_eq!(src["table_identity"], mh);
    assert_eq!(src["covariance_identity"], ch);
    let rows: Vec<_> = std::str::from_utf8(&mb)
        .unwrap()
        .lines()
        .filter(|l| !l.trim().is_empty() && !l.starts_with('#'))
        .map(|l| l.split_whitespace().collect::<Vec<_>>())
        .collect();
    assert_eq!(rows.len(), 13);
    for (i, row) in rows.iter().enumerate() {
        let actual = &src["rows"][i];
        assert_eq!(actual["id"], format!("{mh}:row:{i}:{}", row[2]));
        assert_eq!(
            actual["z"].as_f64().unwrap().to_bits(),
            row[0].parse::<f64>().unwrap().to_bits()
        );
        assert_eq!(
            actual["value"].as_f64().unwrap().to_bits(),
            row[1].parse::<f64>().unwrap().to_bits()
        );
        assert_eq!(actual["observable"], row[2]);
        assert_eq!(src["covariance_axis_ids"][i], actual["id"]);
    }
    assert_eq!(src["rows"][11]["observable"], "DH_over_rs");
    assert_eq!(src["rows"][12]["observable"], "DM_over_rs");
    let e = &r["result"]["evaluations"][0];
    let d = &e["normalized_density"]["value"];
    for (key, index) in [
        ("log_density", 5),
        ("quadratic", 6),
        ("log_determinant", 7),
        ("log_normalization", 8),
    ] {
        assert_eq!(
            d[key].as_f64().unwrap().to_bits(),
            native[index].parse::<f64>().unwrap().to_bits(),
            "{key}"
        );
    }
    for i in 0..13 {
        let prediction = e["predictions"]["value"][i].as_f64().unwrap();
        assert_eq!(
            prediction.to_bits(),
            native[11 + i].parse::<f64>().unwrap().to_bits()
        );
        let observed = rows[i][1].parse::<f64>().unwrap();
        assert_eq!(
            e["residuals"]["value"][i].as_f64().unwrap().to_bits(),
            (observed - prediction).to_bits()
        );
    }
    for (raw, hash) in [(&mb, &mh), (&cb, &ch)] {
        assert_eq!(fs::read(p.join("store/objects").join(hash)).unwrap(), *raw);
    }
    assert_eq!(digest(&fs::read(&mean).unwrap()), mh);
    assert_eq!(digest(&fs::read(&covariance).unwrap()), ch);
}
