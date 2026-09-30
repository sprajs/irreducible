//! Generic magnitude covariance admission is structural, not release verification.
use serde_json::{Value, json};
use std::{fs, process::Command};
#[test]
fn explicit_generic_source_preserves_role_axes_and_unknown_provenance() {
    let root = std::env::temp_dir().join(format!("irred-generic-source-{}", std::process::id()));
    fs::create_dir_all(&root).unwrap();
    let table = root.join("table.txt");
    let cov = root.join("cov.txt");
    fs::write(
        &table,
        "ROW_ID EVENT_ID MAG\nrow-b event-b 23\nrow-a event-a 21\n",
    )
    .unwrap();
    fs::write(&cov, "2\n3 .25 .25 2\n").unwrap();
    let mut request = json!({"schema_version":2,"operation":"observations.prepare","table":table,"uncertainty":cov,"metadata":{"profile":"typed_magnitude_covariance","role":"released_fitted_summary","unit":"magnitude","calibration":"unknown","uncertainty":"covariance","uncertainty_unit":"magnitude_squared","component":"unknown"},"exports":["source_values","covariance","original_fields"],"source_selection":[1,1],"calibration_provenance":"unknown supplied calibration","dependence_provenance":"unknown overlap","quality_dictionary":"unavailable","ordering_provenance":"declared table row order; not independently verified","resources":{"maximum_asset_bytes":4096,"maximum_rows":8,"maximum_matrix_elements":64,"maximum_string_bytes":4096,"maximum_preparation_bytes":16777216}});
    let run = |v: &Value| {
        let path = root.join("request.json");
        fs::write(&path, serde_json::to_vec(v).unwrap()).unwrap();
        let out = Command::new(env!("CARGO_BIN_EXE_irred"))
            .arg("run")
            .arg(path)
            .arg(root.join("store"))
            .output()
            .unwrap();
        (
            out.status.code().unwrap(),
            serde_json::from_slice::<Value>(&out.stdout).unwrap(),
        )
    };
    let (exit, response) = run(&request);
    assert_eq!(exit, 0);
    let result = &response["result"];
    assert_eq!(result["source_type"], request["metadata"]);
    assert_eq!(result["values"], json!([23., 21.]));
    assert_eq!(result["measurement_ids"], json!(["row-b", "row-a"]));
    assert!(result.get("selected_source_indices").is_none());
    assert_eq!(
        response["receipt"]["outputs"][0]["check_kind"],
        "structural_source_contract"
    );
    assert_eq!(response["receipt"]["interpretation"], "unqualified");
    let matrix = result["full_uncertainty_matrix"]["object_digest"]
        .as_str()
        .unwrap();
    let stored: Value =
        serde_json::from_slice(&fs::read(root.join("store/objects").join(matrix)).unwrap())
            .unwrap();
    assert_eq!(stored["axis_ids"], result["measurement_ids"]);
    assert_eq!(stored["values"], json!([3., 0.25, 0.25, 2.]));
    request["metadata"]["role"] = json!("posterior_summary");
    assert_ne!(run(&request).0, 0);
    request["metadata"]["role"] = json!("released_fitted_summary");
    request["selection"] = json!("pantheon_zhd_gt_001");
    assert_ne!(run(&request).0, 0);
    fs::remove_dir_all(root).unwrap();
}
