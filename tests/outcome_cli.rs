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
    serde_json::from_slice(&fs::read(format!("tests/fixtures/{name}.json")).unwrap()).unwrap()
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
        c.env("COSMOLOGY_TEST_PANIC", "1");
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
fn sn_request(s: &Scratch) -> Value {
    json!({"schema_version":1,"operation":"supernova.profile_batch.v1","observations":{"schema_version":1,"operation":"observations.prepare.v1","table":s.0.join("table.dat"),"uncertainty":s.0.join("cov.dat"),"metadata":{"profile":"pantheon_plus_released_v1","role":"released_fitted_summary","unit":"magnitude","calibration":"unknown","uncertainty":"covariance","uncertainty_unit":"magnitude_squared","component":"total"},"selection":"pantheon_zhd_gt_001","ordering_provenance":"generated interface fixture supplied axis order; not actual release or independently verified matrix linkage","resources":{"maximum_asset_bytes":10000,"maximum_rows":3,"maximum_matrix_elements":9,"maximum_string_bytes":4096}},"models":[{"model":"flat_lcdm_late_v1","omega_m":0,"constant_q":0},{"model":"constant_q_flat_v1","omega_m":0,"constant_q":-1}],"policy":{"arithmetic":"longdouble_cpu_v1","include_residual_arrays":true,"maximum_models":4,"maximum_source_rows":3,"maximum_matrix_elements":9,"maximum_array_elements":100,"maximum_native_output_bytes":100000,"maximum_total_evaluations":100000,"maximum_evaluations_per_integral":10000,"maximum_depth":24,"absolute_tolerance":1e-12,"relative_tolerance":1e-12,"maximum_forward_sensitivity":1e-10}})
}
fn bao_policy() -> Value {
    json!({"background":{"maximum_parameters":4,"maximum_queries":16,"maximum_slots":64,"maximum_native_output_bytes":1000000,"maximum_total_evaluations":2000000,"maximum_evaluations_per_integral":100000,"maximum_depth":30,"absolute_tolerance":1e-13,"relative_tolerance":1e-12},"arithmetic":"longdouble_cpu_v1","include_predictions":true,"include_residuals":true,"maximum_models":4,"maximum_rows":16,"maximum_matrix_elements":256,"maximum_string_bytes":4096,"maximum_native_bytes":1000000,"maximum_array_elements":256,"maximum_native_output_bytes":1000000,"maximum_forward_sensitivity":1e-10})
}
fn bao_provenance(source: &mut Value) {
    source["ordering_provenance"] =
        json!("supplied control row order, covariance linkage not independently verified");
    source["calibration_provenance"] = json!("synthetic fixture; released calibration unknown");
    source["dependence_provenance"] =
        json!("declared correlated covariance, cross-probe dependence unknown");
    source["redshift_convention"] = json!("P01/released-effective-redshift/v1");
    source["ruler_convention"] = json!("P01/free-H0rd-km-s-no-early-physics/v1");
    source["computational_h0_convention"] = json!("P01/computational-H0-fixed-70-km-s-Mpc/v1");
}
fn bao_request() -> Value {
    let mut source = json!({"profile":"synthetic_inline_v1","rows":[{"id":"a","z":0.1,"observable":"DM_over_rs","value":299792.458/10000.0*0.1+1.0},{"id":"b","z":0.2,"observable":"DH_over_rs","value":299792.458/10000.0+2.0}],"covariance":[4,1,1,9],"covariance_axis_ids":["a","b"]});
    bao_provenance(&mut source);
    json!({"schema_version":1,"operation":"bao.gaussian_batch.v1","source":source,"prepare_policy":bao_policy(),"evaluation_policy":bao_policy(),"models":[{"model":"constant_q_flat_v1","omega_m":0.0,"constant_q":-1.0,"w0":-1.0,"wa":0.0,"h0_rd_km_s":10000.0},{"model":"flat_lcdm_late_v1","omega_m":0.0,"constant_q":0.0,"w0":-1.0,"wa":0.0,"h0_rd_km_s":10000.0},{"model":"flat_cpl_late_v1","omega_m":0.0,"constant_q":0.0,"w0":-1.0,"wa":0.0,"h0_rd_km_s":10000.0}]})
}
fn pw_query(z: f64) -> Value {
    json!({"z_expansion":z,"z_observer":z,"convention":"geometric_same_redshift"})
}
fn pw_request() -> Value {
    json!({"schema_version":1,"operation":"background.piecewise_query_batch.v1","parameters":[{"h0_km_s_mpc":70.0,"q":[-1.0,-1.0,-1.0,-1.0,-1.0]},{"h0_km_s_mpc":100.0,"q":[-0.4,-0.4,-0.2,0.1,0.3]}],"queries":[pw_query(0.0),pw_query(0.1),pw_query(0.3),pw_query(2.5),{"z_expansion":0.7,"z_observer":0.71,"convention":"released_zhd_zhel"}],"policy":{"maximum_parameters":8,"maximum_queries":16,"maximum_slots":128,"maximum_native_output_bytes":1000000,"maximum_total_segment_visits":1000}})
}
fn cases(s: &Scratch) -> Vec<(Value, &'static str, &'static str, &'static str)> {
    let bg = fixture("background-late");
    let mut bg2 = bg.clone();
    bg2["operation"] = json!("background.parameter_query_batch.v2");
    for a in bg2["parameters"].as_array_mut().unwrap() {
        a["w0"] = json!(-1.0);
        a["wa"] = json!(0.0);
    }
    let sn = sn_request(s);
    let obs = sn["observations"].clone();
    let mut sn2 = sn.clone();
    sn2["operation"] = json!("supernova.profile_batch.v2");
    for a in sn2["models"].as_array_mut().unwrap() {
        a["w0"] = json!(-1.0);
        a["wa"] = json!(0.0);
    }
    let mut sp = sn.clone();
    sp["operation"] = json!("supernova.piecewise_profile_batch.v1");
    sp["models"] = json!([{"q":[-1.,-1.,-1.,-1.,-1.]}]);
    sp["policy"] = json!({"arithmetic":"longdouble_cpu_v1","include_residual_arrays":true,"maximum_models":4,"maximum_source_rows":3,"maximum_matrix_elements":9,"maximum_array_elements":100,"maximum_native_output_bytes":100000,"maximum_queries":3,"maximum_total_segment_visits":100,"maximum_forward_sensitivity":1e-10});
    let mut bp = bao_request();
    bp["operation"] = json!("bao.piecewise_gaussian_batch.v1");
    bp["models"] = json!([{"q":[-1.,-1.,-1.,-1.,-1.],"h0_rd_km_s":10000.0}]);
    let ap = json!({"arithmetic":"longdouble_cpu_v1","include_predictions":true,"include_residuals":true,"maximum_models":4,"maximum_rows":16,"maximum_matrix_elements":256,"maximum_string_bytes":4096,"maximum_native_bytes":1000000,"maximum_array_elements":256,"maximum_native_output_bytes":1000000,"maximum_queries":16,"maximum_total_segment_visits":100,"maximum_forward_sensitivity":1e-10});
    bp["prepare_policy"] = ap.clone();
    bp["evaluation_policy"] = ap;
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
            obs,
            "prepared_observations",
            "typed_observation_preparation",
            "binary64_source_storage_no_scientific_transformation",
        ),
        (
            fixture("gaussian-density"),
            "gaussian_batch",
            "normalized_density",
            "binary64_legacy_v1",
        ),
        (
            bg,
            "background_slots",
            "adaptive_simpson",
            "binary64_storage_longdouble_intermediate",
        ),
        (
            bg2,
            "background_slots",
            "adaptive_simpson",
            "binary64_storage_longdouble_intermediate",
        ),
        (
            pw_request(),
            "piecewise_background_slots",
            "analytic_piecewise_segments",
            "binary64_storage_longdouble_analytic_intermediate",
        ),
        (
            sn,
            "profile_scores",
            "retained_single_offset_profile",
            "longdouble_cpu_v1",
        ),
        (
            sn2,
            "profile_scores",
            "retained_single_offset_profile",
            "longdouble_cpu_v1",
        ),
        (
            sp,
            "profile_scores",
            "retained_single_offset_profile_analytic_segments",
            "longdouble_cpu_v1",
        ),
        (
            bao_request(),
            "bao_densities",
            "retained_normalized_gaussian",
            "longdouble_cpu_v1",
        ),
        (
            bp,
            "bao_densities",
            "retained_normalized_gaussian_analytic_segments",
            "longdouble_cpu_v1",
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
    let mut q = fixture("gaussian-density");
    q["residuals"] = json!([[1., 2.], [f64::MAX, -f64::MAX]]);
    for assurance in [None, Some("qualified")] {
        let (r, code) = run(&s, &q, assurance, false);
        assert_eq!(code, 2);
        assert_eq!(r["receipt"]["execution"], "completed");
        assert_eq!(r["receipt"]["numerical"], "failed");
        assert_eq!(r["receipt"]["assurance"]["satisfied"], false);
        assert_eq!(r["result"]["calculation"]["rows"][0]["kind"], "finite");
        assert_eq!(r["result"]["calculation"]["rows"][1]["kind"], "failure");
        assert!(
            r["result"]["calculation"]["rows"][1]
                .get("quadratic")
                .is_none()
        );
    }
    let mut malformed = fixture("background-late");
    malformed["ignored"] = json!(true);
    let (r, code) = run(&s, &malformed, None, false);
    assert_eq!(code, 2);
    assert_eq!(r["receipt"]["execution"], "failed");
    assert_eq!(r["receipt"]["failure_class"], "request_or_transport");
    assert_eq!(r["receipt"]["numerical"], "not_assessed");
    assert_eq!(r["receipt"]["outputs"], json!([]));
    let (r, code) = run(&s, &fixture("background-late"), None, true);
    assert_eq!(code, 2);
    assert_eq!(r["receipt"]["execution"], "failed");
    assert_eq!(r["receipt"]["failure_class"], "internal");
    assert_eq!(r["result"]["error_id"], "RUST_PANIC");
}

#[test]
fn well_formed_json_with_native_descriptor_rejection_is_transport_failure() {
    let s = scratch();
    let mut q = fixture("background-late");
    // Structurally valid JSON, but above the native policy hard cap. This is
    // descriptor admission, unlike an owned per-row numerical failure.
    q["policy"]["maximum_parameters"] = json!(4097);
    let (r, code) = run(&s, &q, None, false);
    assert_eq!(code, 2);
    assert_eq!(r["receipt"]["execution"], "failed");
    assert_eq!(r["receipt"]["failure_class"], "request_or_transport");
    assert_eq!(r["receipt"]["numerical"], "not_assessed");
    assert_eq!(r["result"]["error_id"], "CORE_STATUS_2");
    assert_eq!(r["receipt"]["outputs"], json!([]));
}
