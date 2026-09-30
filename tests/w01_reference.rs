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
fn check_native(all_named: bool) {
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
    let mut command = Command::new(&executable);
    command.env_remove("IRRED_SN_ALL_NAMED");
    if all_named {
        command.env("IRRED_SN_ALL_NAMED", "1");
    }
    let output = command
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
    let marker = String::from_utf8_lossy(&output.stdout)
        .lines()
        .filter_map(|line| serde_json::from_str::<serde_json::Value>(line).ok())
        .find(|v| {
            v["suite"]
                == if all_named {
                    "current1590_named29"
                } else {
                    "grey12_current1590_conditional_not1820"
                }
        })
        .expect("expected named native suite marker");
    assert_eq!(marker["models"], if all_named { 29 } else { 12 });
    assert_eq!(marker["n"], 1590);
    assert_eq!(marker["passed"], true);
    assert!(
        output.status.success(),
        "native original-input comparison failed: {:?}",
        output.status
    );
}

#[test]
#[ignore = "requires explicit original assets and prebuilt optional native harness"]
fn released_assets_and_native_historical_comparison() {
    check_native(false);
}
#[test]
#[ignore = "requires explicit original assets and prebuilt optional native harness"]
fn released_assets_and_native_current_twenty_nine_points() {
    check_native(true);
}
