//! Independent CLI contracts: exact integer oracle and durable immutable records.
use serde_json::{Value, json};
use sha2::{Digest, Sha256};
use std::{
    collections::BTreeMap,
    fs,
    path::{Path, PathBuf},
    process::{Command, Output},
    time::{SystemTime, UNIX_EPOCH},
};
fn digest(b: &[u8]) -> String {
    format!("{:x}", Sha256::digest(b))
}
struct Scratch(PathBuf);
impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}
fn files(dir: &Path) -> BTreeMap<String, Vec<u8>> {
    let mut out = BTreeMap::new();
    for sub in ["objects", "attempts"] {
        for entry in fs::read_dir(dir.join(sub)).unwrap() {
            let p = entry.unwrap().path();
            if p.is_file() {
                out.insert(
                    format!("{sub}/{}", p.file_name().unwrap().to_str().unwrap()),
                    fs::read(p).unwrap(),
                );
            }
        }
    }
    out
}
fn run(request: &Path, store: &Path, flag: Option<&str>) -> Output {
    // Shell only disables core files for the deliberately aborting child.
    let mut cmd = Command::new("sh");
    cmd.args([
        "-c",
        "ulimit -c 0; exec \"$@\"",
        "verifier",
        env!("CARGO_BIN_EXE_irred"),
        "run",
    ]);
    cmd.arg(request).arg(store);
    if let Some(f) = flag {
        cmd.env(f, "1");
    }
    cmd.output().unwrap()
}
fn decode(out: &Output) -> Value {
    serde_json::from_slice(&out.stdout).unwrap()
}
#[test]
fn independent_records_adversarial() {
    let dir = Scratch(
        std::env::temp_dir().join(format!(
            "cosmology-native-verifier-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        )),
    );
    fs::create_dir(&dir.0).unwrap();
    let request = dir.0.join("request.json");
    let store = dir.0.join("store");
    let mut payload = json!({"schema_version":2,"operation":"fixture.checked_i64_add.v1","a":[0,1,-2,1048576,-1048576],"b":[0,-1,3,-1048576,1048576],"fault":0});
    fs::write(&request, serde_json::to_vec(&payload).unwrap()).unwrap();
    let first = run(&request, &store, None);
    assert!(first.status.success());
    let first = decode(&first);
    assert_eq!(first["result"]["values"], json!([0, 0, 1, 0, 0]));
    let receipt = &first["receipt"];
    assert_eq!(receipt["execution"], "completed");
    assert_eq!(receipt["numerical"], "checks_passed");
    let old = files(&store);
    let second = run(&request, &store, None);
    assert!(second.status.success());
    let second = decode(&second);
    assert_ne!(receipt["attempt_id"], second["receipt"]["attempt_id"]);
    for (name, bytes) in old {
        assert_eq!(
            fs::read(store.join(name)).unwrap(),
            bytes,
            "duplicate run overwrote immutable bytes"
        );
    }
    for key in ["input_digest", "output_digest"] {
        let name = receipt[key].as_str().unwrap();
        assert_eq!(
            digest(&fs::read(store.join("objects").join(name)).unwrap()),
            name
        );
    }
    assert_eq!(
        digest(&fs::read(env!("CARGO_BIN_EXE_irred")).unwrap()),
        receipt["executable_digest"]
    );
    let before = fs::read(&request).unwrap();
    payload["a"][1] = json!(2);
    fs::write(&request, serde_json::to_vec(&payload).unwrap()).unwrap();
    assert_eq!(before.len(), fs::read(&request).unwrap().len());
    let changed = decode(&run(&request, &store, None));
    assert_ne!(receipt["input_digest"], changed["receipt"]["input_digest"]);
    assert_eq!(changed["result"]["values"][1], 1);
    for fault in [1, 2] {
        payload["fault"] = json!(fault);
        fs::write(&request, serde_json::to_vec(&payload).unwrap()).unwrap();
        let out = run(&request, &store, None);
        assert!(!out.status.success());
        let v = decode(&out);
        assert_eq!(v["receipt"]["execution"], "failed");
        assert_eq!(v["result"]["kind"], "failure");
    }
    payload["fault"] = json!(0);
    fs::write(&request, serde_json::to_vec(&payload).unwrap()).unwrap();
    let out = run(&request, &store, Some("COSMOLOGY_TEST_PANIC"));
    assert!(!out.status.success());
    assert_eq!(decode(&out)["receipt"]["execution"], "failed");
    let before = files(&store);
    let out = run(&request, &store, Some("COSMOLOGY_TEST_ABORT"));
    assert!(!out.status.success());
    let after = files(&store);
    let added: Vec<_> = after
        .iter()
        .filter(|(name, _)| !before.contains_key(*name))
        .collect();
    assert_eq!(added.len(), 1);
    assert!(added[0].0.ends_with(".incomplete.json"));
    let incomplete: Value = serde_json::from_slice(added[0].1).unwrap();
    assert_eq!(incomplete["execution"], "incomplete");
    let path = store
        .join("objects")
        .join(digest(&fs::read(&request).unwrap()));
    fs::write(path, b"corrupted").unwrap();
    let out = run(&request, &store, None);
    assert!(!out.status.success());
    assert!(out.stdout.is_empty());
    assert!(String::from_utf8_lossy(&out.stderr).contains("immutable object conflict"));
}

#[allow(dead_code)]
#[path = "../src/abi_generated.rs"]
mod independently_checked_layout;
#[test]
fn host_layout_and_status_values() {
    use independently_checked_layout::*;
    assert_eq!(
        std::mem::size_of::<usize>(),
        8,
        "current 64-bit host evidence only"
    );
    assert_eq!(std::mem::size_of::<Buffer>(), 40);
    assert_eq!(std::mem::align_of::<Buffer>(), 8);
    assert_eq!(
        [
            std::mem::offset_of!(Buffer, struct_size),
            std::mem::offset_of!(Buffer, abi_version),
            std::mem::offset_of!(Buffer, element_type),
            std::mem::offset_of!(Buffer, reserved),
            std::mem::offset_of!(Buffer, data),
            std::mem::offset_of!(Buffer, length),
            std::mem::offset_of!(Buffer, byte_length)
        ],
        [0, 4, 8, 12, 16, 24, 32]
    );
    assert_eq!(
        [
            ABI_VERSION,
            OK,
            ABI_MISMATCH,
            INVALID_INPUT,
            OVERFLOW,
            ALLOCATION_FAILURE,
            EXCEPTION
        ],
        [2, 0, 1, 2, 3, 4, 5]
    );
}

#[test]
fn scientific_identity_strict_parse_and_late_abort() {
    let dir = Scratch(
        std::env::temp_dir().join(format!(
            "cosmology-native-spec-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        )),
    );
    fs::create_dir(&dir.0).unwrap();
    let request = dir.0.join("request.json");
    let store = dir.0.join("store");
    let bytes=b"{\"schema_version\":2,\"operation\":\"fixture.checked_i64_add.v1\",\"a\":[2],\"b\":[3],\"fault\":0}";
    fs::write(&request, bytes).unwrap();
    let first = decode(&run(&request, &store, None));
    let receipt = &first["receipt"];
    let spec = receipt["scientific_specification_digest"].as_str().unwrap();
    let contents = fs::read(store.join("objects").join(spec)).unwrap();
    assert_eq!(digest(&contents), spec);
    let parsed: Value = serde_json::from_slice(&contents).unwrap();
    assert_eq!(parsed["equation_id"], "EQ-fixture.checked_i64_add.v1");
    for (path, hash) in receipt["runtime_libraries"].as_object().unwrap() {
        assert_eq!(digest(&fs::read(path).unwrap()), hash.as_str().unwrap());
    }
    fs::write(&request,b"{ \"b\": [3], \"a\": [2], \"operation\": \"fixture.checked_i64_add.v1\", \"schema_version\": 2 }").unwrap();
    let reordered = decode(&run(&request, &store, None));
    assert_eq!(
        reordered["receipt"]["scientific_specification_digest"],
        spec
    );
    assert_ne!(
        reordered["receipt"]["input_digest"],
        receipt["input_digest"]
    );
    assert_ne!(
        reordered["receipt"]["execution_identity"],
        receipt["execution_identity"]
    );
    fs::write(&request,b"{\"schema_version\":2,\"operation\":\"fixture.checked_i64_add.v1\",\"a\":[2],\"b\":[3],\"fault\":1}").unwrap();
    let fault = decode(&run(&request, &store, None));
    assert_eq!(fault["receipt"]["scientific_specification_digest"], spec);
    assert_ne!(
        fault["receipt"]["execution_identity"],
        receipt["execution_identity"]
    );
    for invalid in [b"{\"schema_version\":2,\"schema_version\":2,\"operation\":\"fixture.checked_i64_add.v1\",\"a\":[2],\"b\":[3]}".as_slice(),b"{\"schema_version\":2,\"operation\":\"fixture.checked_i64_add.v1\",\"a\":[2],\"b\":[3],\"ignored\":1}".as_slice()] {
  fs::write(&request,invalid).unwrap();let out=run(&request,&store,None);assert!(!out.status.success());let v=decode(&out);assert_eq!(v["receipt"]["execution"],"failed");assert_eq!(v["result"]["kind"],"failure");
 }
    for args in [
        vec!["describe"],
        vec!["describe", "--json", "extra"],
        vec!["version", "--unknown"],
        vec!["run", "x", "y", "extra"],
    ] {
        let out = Command::new(env!("CARGO_BIN_EXE_irred"))
            .args(args)
            .output()
            .unwrap();
        assert!(!out.status.success());
        assert!(out.stdout.is_empty());
    }
    fs::write(&request, bytes).unwrap();
    let late_store = dir.0.join("late-store");
    let out = run(
        &request,
        &late_store,
        Some("COSMOLOGY_TEST_ABORT_AFTER_OUTPUT"),
    );
    assert!(!out.status.success());
    let records = files(&late_store);
    let attempts: Vec<_> = records
        .keys()
        .filter(|k| k.starts_with("attempts/"))
        .collect();
    assert_eq!(attempts.len(), 1);
    assert!(attempts[0].ends_with(".incomplete.json"));
    assert_eq!(
        records.keys().filter(|k| k.starts_with("objects/")).count(),
        2,
        "input and output retained without final receipt"
    );
}
