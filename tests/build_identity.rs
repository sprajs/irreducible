//! Run this ignored, exclusive build test by executing its compiled test binary.
//! It temporarily edits sources and an included CMake fragment, restoring on unwind.
use serde_json::Value;
use sha2::{Digest, Sha256};
use std::{
    fs,
    io::Write,
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
        "stdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&out.stdout),
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
    let fragment_path = "cpp/tests/build_identity_probe.cmake";
    let mut fragment_file = fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(root.join(fragment_path))
        .expect("exclusive fragment path must not replace an existing file");
    let fragment = Remove(root.join(fragment_path));
    let marker_path = root.join("build/native/build-identity-probe.txt");
    assert!(
        !marker_path.exists(),
        "exclusive marker path already exists"
    );
    let marker = Remove(marker_path);
    let original_fragment =
        b"file(WRITE \"${CMAKE_BINARY_DIR}/build-identity-probe.txt\" \"original\")\n";
    fragment_file.write_all(original_fragment).unwrap();
    drop(fragment_file);
    let cmake = Restore {
        path: root.join("cpp/CMakeLists.txt"),
        bytes: fs::read(root.join("cpp/CMakeLists.txt")).unwrap(),
    };
    let mut included = cmake.bytes.clone();
    included.extend_from_slice(
        b"\ninclude(\"${CMAKE_CURRENT_LIST_DIR}/tests/build_identity_probe.cmake\")\n",
    );
    fs::write(&cmake.path, included).unwrap();
    let included_base = build(&root);
    let original_digest = format!("{:x}", Sha256::digest(original_fragment));
    let manifest = || -> Value {
        serde_json::from_slice(&fs::read(root.join("build/build-manifest.json")).unwrap()).unwrap()
    };
    assert_eq!(manifest()["sources"][fragment_path], original_digest);
    assert_eq!(fs::read(&marker.0).unwrap(), b"original");
    let fragment_restore = Restore {
        path: fragment.0.clone(),
        bytes: original_fragment.to_vec(),
    };
    let changed_fragment =
        b"file(WRITE \"${CMAKE_BINARY_DIR}/build-identity-probe.txt\" \"changed\")\n";
    fs::write(&fragment.0, changed_fragment).unwrap();
    assert_ne!(
        build(&root),
        included_base,
        "included CMake fragment mutation missing from identity"
    );
    let changed_digest = format!("{:x}", Sha256::digest(changed_fragment));
    assert_ne!(changed_digest, original_digest);
    assert_eq!(manifest()["sources"][fragment_path], changed_digest);
    assert_eq!(fs::read(&marker.0).unwrap(), b"changed");
    drop(fragment_restore);
    assert_eq!(
        build(&root),
        included_base,
        "restored CMake fragment identity differs"
    );
    assert_eq!(manifest()["sources"][fragment_path], original_digest);
    assert_eq!(fs::read(&marker.0).unwrap(), b"original");
    drop(cmake);
    drop(fragment);
    drop(marker);
    assert_eq!(build(&root), base, "removed CMake probe identity differs");
    assert!(manifest()["sources"].get(fragment_path).is_none());
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
