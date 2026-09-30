//! Run this ignored, exclusive build test by executing its compiled test binary.
//! It temporarily edits source comments, restores bytes on unwind, and rebuilds.
use serde_json::Value;
use std::{
    fs,
    path::{Path, PathBuf},
    process::Command,
};
struct Restore {
    path: PathBuf,
    bytes: Vec<u8>,
}
impl Drop for Restore {
    fn drop(&mut self) {
        fs::write(&self.path, &self.bytes).expect("restore original source");
    }
}
struct Remove(PathBuf);
impl Drop for Remove {
    fn drop(&mut self) {
        let _ = fs::remove_file(&self.0);
    }
}
fn build(root: &Path) -> String {
    let out = Command::new("python3")
        .arg("tools/build.py")
        .current_dir(root)
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let m: Value =
        serde_json::from_slice(&fs::read(root.join("build/build-manifest.json")).unwrap()).unwrap();
    m["build_id"].as_str().unwrap().to_owned()
}
#[test]
#[ignore = "exclusive build/mutation test: execute compiled test binary directly after cargo test --no-run"]
fn source_receipt_flags_and_cache_identity() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let base = build(&root);
    let probe = Remove(root.join("build/verifier-native-mutable-probe.txt"));
    fs::write(&probe.0, b"mutable receipt\n").unwrap();
    assert_eq!(
        build(&root),
        base,
        "mutable receipt affected scientific identity"
    );
    drop(probe);
    for path in ["src/main.rs", "cpp/src/abi.cpp"] {
        let path = root.join(path);
        let restore = Restore {
            bytes: fs::read(&path).unwrap(),
            path,
        };
        let mut changed = restore.bytes.clone();
        changed.extend_from_slice(b"\n// independent native build-identity mutation\n");
        fs::write(&restore.path, changed).unwrap();
        assert_ne!(build(&root), base, "source mutation missing from identity");
        drop(restore);
        assert_eq!(build(&root), base, "restored source identity differs");
    }
    for name in [
        "RUSTFLAGS",
        "CARGO_ENCODED_RUSTFLAGS",
        "CXXFLAGS",
        "CXX",
        "LDFLAGS",
    ] {
        let out = Command::new("python3")
            .arg("tools/build.py")
            .current_dir(&root)
            .env(name, "verifier_unsupported_override")
            .output()
            .unwrap();
        assert!(!out.status.success());
        assert!(String::from_utf8_lossy(&out.stderr).contains(name));
    }
    let out = Command::new(root.join(".build-tools/bin/cmake"))
        .args([
            "-S",
            "cpp",
            "-B",
            "build/native",
            "-DCMAKE_CXX_FLAGS=-ffast-math",
        ])
        .current_dir(&root)
        .output()
        .unwrap();
    assert!(out.status.success());
    assert_eq!(build(&root), base, "cache reset failed to restore identity");
    let m: Value =
        serde_json::from_slice(&fs::read(root.join("build/build-manifest.json")).unwrap()).unwrap();
    for row in m["native_compile_commands"].as_array().unwrap() {
        assert!(!row["command"].as_str().unwrap().contains("-ffast-math"));
    }
}
