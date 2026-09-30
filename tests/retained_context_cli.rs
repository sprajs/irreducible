//! Permanent public retained-context controls, including exact one-shot parity.
//! Boundary evidence is separate from original-data scientific qualification.
use serde_json::{Value, json};
use sha2::{Digest, Sha256};
use std::{
    fs,
    io::{BufRead, BufReader, Write},
    path::PathBuf,
    process::{Child, ChildStdin, ChildStdout, Command, Stdio},
    time::{Duration, Instant, SystemTime, UNIX_EPOCH},
};
fn root() -> PathBuf {
    let p = std::env::temp_dir().join(format!(
        "irred-retained-peer-{}-{}",
        std::process::id(),
        SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    fs::create_dir_all(&p).unwrap();
    p
}
fn fixtures() -> Vec<Value> {
    let p = root();
    fs::write(
        p.join("table.txt"),
        "ROW_ID EVENT_ID MAG\nrow-b event-b 23\nrow-a event-a 21\n",
    )
    .unwrap();
    fs::write(p.join("covariance.txt"), "2\n3 .25 .25 2\n").unwrap();
    vec![
        json!({"action":"prepare_observations","request":{"schema_version":2,
        "operation":"observations.prepare","table":p.join("table.txt"),"uncertainty":p.join("covariance.txt"),
        "metadata":{"profile":"typed_magnitude_covariance","role":"released_fitted_summary",
            "unit":"magnitude","calibration":"unknown","uncertainty":"covariance",
            "uncertainty_unit":"magnitude_squared","component":"unknown"},
        "resources":{"maximum_asset_bytes":4096,"maximum_rows":8,"maximum_matrix_elements":64,
            "maximum_string_bytes":4096,"maximum_preparation_bytes":8388608},"exports":[],
        "calibration_provenance":"unknown explicit fitted summary","dependence_provenance":"unknown",
        "quality_dictionary":"not supplied","ordering_provenance":"supplied table row order; axes not independently verified",
        "source_selection":[1,1]}}),
        json!({"action":"prepare_supernova","source_handle":1,
            "selection":{"kind":"explicit","source_indices":[0,1],"coordinates":[
                {"z_expansion":1.0,"observer":{"redshift":1.0,"convention":"geometric_same_redshift"}},
                {"z_expansion":2.0,"observer":{"redshift":2.0,"convention":"geometric_same_redshift"}}]},
            "preparation_policy":{"arithmetic":"wide","maximum_selected_rows":2,"maximum_matrix_elements":4,
                "maximum_string_bytes":4096,"maximum_native_bytes":8388608,"maximum_forward_sensitivity":1e-10}}),
        json!({"action":"release","handle":1}),
        json!({"action":"supernova_profile","handle":2,"models":[{"expansion":{"kind":"constant_q","q":-1.0},
            "geometry":{"kind":"flat_flrw"},"source_effect":{"kind":"none"}}],
            "requested_outputs":["score","geometric_shape","magnitude_effect","corrected_residuals","profiled_residuals","diagnostics"],
            "numerical_policy":{"arithmetic":"wide","projection":{"maximum_queries":2,"maximum_callbacks":10000,
                "maximum_segment_visits":0,"integration":{"absolute_tolerance":1e-12,"relative_tolerance":1e-12,
                    "maximum_evaluations":10000,"maximum_depth":30}},"maximum_models":2,"maximum_array_elements":8,
                "maximum_native_bytes":8388608,"maximum_forward_sensitivity":1e-10}}),
        json!({"action":"release","handle":2}),
    ]
}
struct Session {
    child: Child,
    input: ChildStdin,
    output: BufReader<ChildStdout>,
    store: PathBuf,
    commands: u64,
    cap: usize,
}
impl Session {
    fn new(commands: u64, cap: usize) -> Self {
        let store = root();
        let executable = std::env::var("IRRED_PUBLIC_EXECUTABLE")
            .unwrap_or_else(|_| env!("CARGO_BIN_EXE_irred").into());
        let mut command = Command::new(executable);
        let limits = store.join("limits.json");
        fs::write(
            &limits,
            serde_json::to_vec(
                &json!({"sources":16,"consumers":16,"retained_bytes":8388608,
            "commands":commands,"input_line_bytes":16*1024*1024,"output_line_bytes":cap}),
            )
            .unwrap(),
        )
        .unwrap();
        command.arg("stream").arg(&store).arg(limits);
        let mut child = command
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::inherit())
            .spawn()
            .unwrap();
        let input = child.stdin.take().unwrap();
        let output = BufReader::new(child.stdout.take().unwrap());
        Self {
            child,
            input,
            output,
            store,
            commands,
            cap,
        }
    }
    fn send(&mut self, x: &Value) -> Value {
        let mut bytes = serde_json::to_vec(x).unwrap();
        bytes.push(b'\n');
        self.send_raw(&bytes)
    }
    fn send_raw(&mut self, bytes: &[u8]) -> Value {
        self.input.write_all(bytes).unwrap();
        self.input.flush().unwrap();
        let mut line = String::new();
        assert!(self.output.read_line(&mut line).unwrap() > 0);
        let envelope: Value = serde_json::from_str(&line).unwrap();
        {
            let receipt = &envelope["receipt"];
            assert_eq!(receipt["requested_assurance"], "numerical_contract");
            let config_digest = receipt["execution_configuration_digest"].as_str().unwrap();
            let config_bytes = fs::read(self.store.join("objects").join(config_digest)).unwrap();
            assert_eq!(
                format!("{:x}", Sha256::digest(&config_bytes)),
                config_digest
            );
            let config: Value = serde_json::from_slice(&config_bytes).unwrap();
            assert_eq!(config["limits"]["commands"], self.commands);
            assert_eq!(config["limits"]["output_line_bytes"], self.cap);
            assert_eq!(config["limits"]["retained_bytes"], 8388608);
            assert_eq!(config["assurance"], "numerical_contract");
            let digest = receipt["output_digest"].as_str().unwrap();
            let stored = fs::read(self.store.join("objects").join(digest)).unwrap();
            assert_eq!(format!("{:x}", Sha256::digest(&stored)), digest);
            assert_eq!(
                serde_json::from_slice::<Value>(&stored).unwrap(),
                envelope["result"]
            );
            if bytes.len() <= 16 * 1024 * 1024 {
                assert_eq!(
                    receipt["input_digest"],
                    format!("{:x}", Sha256::digest(&bytes))
                );
            }
        }
        let reply = envelope["result"].clone();
        if bytes.len() <= 16 * 1024 * 1024 {
            let digest = format!("{:x}", Sha256::digest(&bytes));
            assert_eq!(
                fs::read(self.store.join("objects").join(digest)).unwrap(),
                bytes
            );
        }
        if let Some(digest) = reply["output_digest"].as_str() {
            let stored = fs::read(self.store.join("objects").join(digest)).unwrap();
            assert_eq!(format!("{:x}", Sha256::digest(&stored)), digest);
            assert_eq!(
                serde_json::from_slice::<Value>(&stored).unwrap(),
                reply["output"]
            );
        }
        reply
    }
    fn close(mut self) {
        drop(self.input);
        assert!(self.child.wait().unwrap().success());
        let mut extra = String::new();
        assert_eq!(self.output.read_line(&mut extra).unwrap(), 0);
    }
}
#[test]
fn actual_retained_identity_lifetime_and_no_reacquisition() {
    let f = fixtures();
    let mut a = Session::new(100, 1 << 20);
    let source = a.send(&f[0]);
    let prep = a.send(&f[1]);
    assert_eq!(a.send(&f[2])["disposition"], "released");
    let output = a.send(&f[3]);
    assert_eq!(output["accepted"], true);
    {
        let request = json!({"schema_version":2,"operation":"supernova.profile",
            "observations":f[0]["request"],"selection":f[1]["selection"],
            "preparation_policy":f[1]["preparation_policy"],"models":f[3]["models"],
            "requested_outputs":f[3]["requested_outputs"],"numerical_policy":f[3]["numerical_policy"]});
        let one_shot = root();
        let path = one_shot.join("request.json");
        fs::write(&path, serde_json::to_vec(&request).unwrap()).unwrap();
        let executed = Command::new(
            std::env::var("IRRED_PUBLIC_EXECUTABLE")
                .unwrap_or_else(|_| env!("CARGO_BIN_EXE_irred").into()),
        )
        .arg("run")
        .arg(path)
        .arg(one_shot.join("store"))
        .output()
        .unwrap();
        assert_eq!(
            executed.status.code(),
            Some(0),
            "{}",
            String::from_utf8_lossy(&executed.stdout)
        );
        let recorded: Value = serde_json::from_slice(&executed.stdout).unwrap();
        assert_eq!(recorded["result"], output["output"]);
        assert_eq!(
            recorded["receipt"]["scientific_specification_digest"],
            output["scientific_specification_digest"]
        );
        assert_eq!(recorded["receipt"]["accepted"], true);
    }

    assert_eq!(a.send(&f[4])["disposition"], "released");
    a.close();
    let mut b = Session::new(80, 2 << 20);
    b.send(&f[0]);
    b.send(&f[0]);
    let occupied = b.send(&f[1]);
    assert_eq!(prep["preparation_digest"], occupied["preparation_digest"]);
    assert_eq!(prep["specification"], occupied["specification"]);
    assert_ne!(
        prep["effective_runtime_policy"],
        occupied["effective_runtime_policy"]
    );
    let mut wrong_kind = f[3].clone();
    wrong_kind["handle"] = json!(1);
    let rejected = b.send(&wrong_kind);
    assert_eq!(rejected["execution"], "failed");
    assert_eq!(rejected["error_id"], "WRONG_HANDLE_KIND");
    let mut eval = f[3].clone();
    eval["handle"] = json!(3);
    assert_eq!(output["output"], b.send(&eval)["output"]);
    b.close();
    let mut c = Session::new(100, 1 << 20);
    c.send(&f[0]);
    c.send(&f[1]);
    c.send(&f[1]);
    c.send(&f[2]);
    let first = c.send(&f[3]);
    assert_eq!(first["output"], c.send(&eval)["output"]);
    c.send(&f[4]);
    assert_eq!(first["output"], c.send(&eval)["output"]);
    assert_eq!(
        c.send(&json!({"action":"release","handle":3}))["disposition"],
        "released"
    );
    let recreated = c.send(&f[0]);
    assert_eq!(recreated["handle"], 4);
    assert_eq!(
        source["preparation_digest"],
        recreated["preparation_digest"]
    );
    assert_eq!(c.send(&f[2])["execution"], "failed");
    c.close();
    // Only disposable copies are made unavailable. Original inputs untouched.
    let p = root();
    let mut copied = f[0].clone();
    for key in ["table", "uncertainty"] {
        let original = copied["request"][key].as_str().unwrap();
        let target = p.join(key);
        fs::copy(original, &target).unwrap();
        copied["request"][key] = json!(target);
    }
    let mut d = Session::new(100, 1 << 20);
    d.send(&copied);
    d.send(&f[1]);
    fs::remove_file(p.join("table")).unwrap();
    fs::remove_file(p.join("uncertainty")).unwrap();
    d.send(&f[2]);
    assert_eq!(output["output"], d.send(&f[3])["output"]);
    assert_eq!(output["output"], d.send(&f[3])["output"]);
    d.send(&f[4]);
    d.close();
}
#[test]
fn actual_reply_rollback_frame_drain_and_command_eof() {
    let f = fixtures();
    let mut s = Session::new(100, 1024);
    for _ in 0..3 {
        let reply = s.send(&f[0]);
        assert_eq!(reply["execution"], "failed");
        assert_eq!(reply["error_id"], "OUTPUT_LINE_LIMIT");
        assert!(reply.get("handle").is_none());
    }
    assert_eq!(s.send(&f[2])["execution"], "failed");
    s.close();
    let mut s = Session::new(100, 1 << 20);
    let reply = s.send(&json!({"action":"X".repeat(17*1024*1024)}));
    assert_eq!(reply["error_id"], "INPUT_LINE_LIMIT");
    assert_eq!(
        s.send(&json!({"action":"release","handle":999}))["execution"],
        "failed"
    );
    s.close();
    let mut s = Session::new(1, 1 << 20);
    writeln!(s.input, "{}", json!({"action":"release","handle":999})).unwrap();
    s.input
        .write_all(b"{\"action\":\"release\",\"handle\":")
        .unwrap();
    s.input.flush().unwrap();
    let mut reply = String::new();
    assert!(s.output.read_line(&mut reply).unwrap() > 0);
    assert_eq!(
        {
            let x: Value = serde_json::from_str(&reply).unwrap();
            x["result"]["execution"].clone()
        },
        "failed"
    );
    let deadline = Instant::now() + Duration::from_secs(3);
    loop {
        if let Some(status) = s.child.try_wait().unwrap() {
            assert!(status.success());
            break;
        }
        assert!(
            Instant::now() < deadline,
            "N+1 unterminated frame must not be read"
        );
        std::thread::sleep(Duration::from_millis(10));
    }
    let mut extra = String::new();
    assert_eq!(s.output.read_line(&mut extra).unwrap(), 0);
}
#[test]
fn actual_nested_required_dimensional_failures_preserve_dimensionless_groups() {
    let request = json!({"schema_version":2,"operation":"background.evaluate",
        "geometry":{"kind":"flat_flrw"},"models":[{"kind":"constant_q","q":0.0}],
        "queries":[{"z_expansion":0.5,"requested_outputs":["expansion","clock"],
                    "physical_scale":{"h0_km_s_mpc":0.0}}],
        "numerical_policy":{"maximum_models":4,"maximum_queries":16,"maximum_slots":64,
            "maximum_callbacks":100000,"maximum_segment_visits":0,"maximum_native_bytes":1000000,
            "integration":{"absolute_tolerance":1e-12,"relative_tolerance":1e-12,
                           "maximum_evaluations":100000,"maximum_depth":24}}});
    let mut s = Session::new(10, 1 << 20);
    let reply = s.send(&json!({"action":"background_evaluate","request":request}));
    assert_eq!(reply["execution"], "completed");
    assert_eq!(reply["accepted"], false);
    assert_eq!(reply["numerical"], "failed");
    let groups = &reply["output"]["evaluations"][0]["groups"];
    assert_eq!(groups["expansion"]["availability"], "available");
    assert_eq!(
        groups["expansion"]["value"]["H_km_s_mpc"]["availability"],
        "failed"
    );
    assert_eq!(groups["clock"]["availability"], "available");
    assert_eq!(
        groups["clock"]["value"]["lookback_seconds"]["availability"],
        "failed"
    );
    assert!(groups["clock"]["value"]["lookback_seconds"]["value"].is_null());
    s.close();
}

#[test]
fn raw_duplicate_rejection_preserves_live_owner_and_quota() {
    let f = fixtures();
    let mut s = Session::new(100, 1 << 20);
    // Duplicate a legitimate nested acquisition field before any handle exists.
    let raw = serde_json::to_string(&f[0]).unwrap().replace(
        "\"maximum_rows\":8",
        "\"maximum_rows\":0,\"maximum_rows\":8",
    );
    assert_ne!(raw, serde_json::to_string(&f[0]).unwrap());
    let rejected = s.send_raw(format!("{raw}\n").as_bytes());
    assert_eq!(rejected["error_id"], "INVALID_SESSION_COMMAND");
    assert_eq!(s.send(&f[0])["handle"], 1);
    let prepared = s.send(&f[1]);
    assert_eq!(prepared["handle"], 2);
    let original = s.send(&f[3]);
    assert_eq!(original["accepted"], true);
    for raw in [
        "{\"action\":\"release\",\"handle\":999,\"handle\":2}\n".to_string(),
        serde_json::to_string(&f[1]).unwrap().replace(
            "\"source_handle\":1",
            "\"source_handle\":999,\"source_handle\":1",
        ) + "\n",
        serde_json::to_string(&f[3])
            .unwrap()
            .replace("\"handle\":2", "\"handle\":999,\"handle\":2")
            + "\n",
    ] {
        let rejected = s.send_raw(raw.as_bytes());
        assert_eq!(rejected["execution"], "failed");
        assert_eq!(rejected["error_id"], "INVALID_SESSION_COMMAND");
        assert_eq!(s.send(&f[3])["output"], original["output"]);
    }
    assert_eq!(s.send(&f[4])["disposition"], "released");
    assert_eq!(s.send(&f[1])["handle"], 3);
    s.send(&json!({"action":"release","handle":3}));
    assert_eq!(s.send(&f[2])["disposition"], "released");
    assert_eq!(s.send(&f[0])["handle"], 4);
    s.close();
}
