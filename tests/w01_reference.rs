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
