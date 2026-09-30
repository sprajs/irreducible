//! Execute in an isolated checkout after both profiles have been built.
//! This checks build identity/embedding and that Release native checks execute.
use serde_json::Value;
use sha2::{Digest, Sha256};
use std::{fs, path::PathBuf, process::Command};
#[test]
#[ignore = "requires isolated checkout and prebuilt Debug/Release profiles"]
fn profile_manifests_and_release_negative_witness() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let mut manifests = Vec::new();
    for (profile, file) in [
        ("debug", "build/build-manifest.json"),
        ("release", "build/build-manifest-release.json"),
    ] {
        let manifest: Value = serde_json::from_slice(&fs::read(root.join(file)).unwrap()).unwrap();
        assert_eq!(manifest["profile"], profile);
        for command in ["describe", "version"] {
            let output = Command::new(root.join(format!("target/{profile}/irred")))
                .args([command, "--json"])
                .output()
                .unwrap();
            assert!(output.status.success());
            let discovery: Value = serde_json::from_slice(&output.stdout).unwrap();
            assert_eq!(
                discovery["build"], manifest,
                "wrong embedded active-profile manifest"
            );
        }
        manifests.push(manifest);
    }
    assert_eq!(manifests[0]["sources"], manifests[1]["sources"]);
    assert_ne!(manifests[0]["build_id"], manifests[1]["build_id"]);
    assert_ne!(manifests[0]["flags"], manifests[1]["flags"]);
    assert_ne!(manifests[0]["rust_flags"], manifests[1]["rust_flags"]);
    for (path, digest) in manifests[0]["sources"].as_object().unwrap() {
        assert_eq!(
            format!("{:x}", Sha256::digest(fs::read(root.join(path)).unwrap())),
            digest.as_str().unwrap(),
            "source identity mismatch: {path}"
        );
    }
    assert_ne!(
        Sha256::digest(fs::read(root.join("build/native/libirred_core.a")).unwrap()),
        Sha256::digest(fs::read(root.join("build/native-release/libirred_core.a")).unwrap())
    );
    let scratch = root.join("build/profile-negative-witness");
    fs::create_dir_all(&scratch).unwrap();
    let original = fs::read_to_string(root.join("cpp/tests/test_abi.cpp")).unwrap();
    assert!(original.contains("p[1] == 2"));
    let source = scratch.join("test.cpp");
    fs::write(&source, original.replace("p[1] == 2", "p[1] == 99")).unwrap();
    let executable = scratch.join("test");
    let compiler = manifests[1]["native_compile_commands"][0]["command"]
        .as_str()
        .unwrap()
        .split_whitespace()
        .next()
        .unwrap();
    let compile = Command::new(compiler)
        .args([
            "-std=c++20",
            "-O3",
            "-DNDEBUG",
            "-fno-fast-math",
            "-ffp-contract=off",
            "-I",
        ])
        .arg(root.join("cpp/include"))
        .arg(&source)
        .arg(root.join("build/native-release/libirred_core.a"))
        .arg("-o")
        .arg(&executable)
        .output()
        .unwrap();
    assert!(
        compile.status.success(),
        "{}",
        String::from_utf8_lossy(&compile.stderr)
    );
    let output = Command::new("bash")
        .args(["-c", "ulimit -c 0; exec \"$@\"", "--"])
        .arg(&executable)
        .output()
        .unwrap();
    assert!(
        !output.status.success(),
        "Release compiled away a required check"
    );
    assert!(String::from_utf8_lossy(&output.stderr).contains("ABI fixture check failed"));
    fs::remove_dir_all(scratch).unwrap();
}
