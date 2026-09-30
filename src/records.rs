use serde_json::{Value, json};
use sha2::{Digest, Sha256};
use std::{
    fs::{self, File, OpenOptions},
    io::Write,
    path::Path,
};
pub(crate) fn hash(b: &[u8]) -> String {
    format!("{:x}", Sha256::digest(b))
}
pub(crate) fn publish(path: &Path, bytes: &[u8]) -> Result<(), String> {
    if path.exists() {
        if fs::read(path).map_err(|e| e.to_string())? == bytes {
            return Ok(());
        }
        return Err("immutable object conflict".into());
    }
    let tmp = path.with_extension(format!("stage-{}", std::process::id()));
    let mut f = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&tmp)
        .map_err(|e| e.to_string())?;
    f.write_all(bytes)
        .and_then(|_| f.sync_all())
        .map_err(|e| e.to_string())?;
    if let Err(e) = fs::hard_link(&tmp, path) {
        if e.kind() != std::io::ErrorKind::AlreadyExists
            || fs::read(path).map_err(|e| e.to_string())? != bytes
        {
            return Err(e.to_string());
        }
    }
    fs::remove_file(tmp).map_err(|e| e.to_string())?;
    File::open(path.parent().unwrap())
        .and_then(|f| f.sync_all())
        .map_err(|e| e.to_string())?;
    Ok(())
}
pub(crate) fn runtime_libraries() -> Result<Value, String> {
    let maps = fs::read_to_string("/proc/self/maps").map_err(|e| e.to_string())?;
    let paths: std::collections::BTreeSet<_> = maps
        .lines()
        .filter_map(|line| line.split_whitespace().last())
        .filter(|p| p.starts_with('/') && p.contains(".so"))
        .collect();
    let mut libs = serde_json::Map::new();
    for p in paths {
        let bytes = fs::read(p).map_err(|e| e.to_string())?;
        libs.insert(p.to_string(), json!(hash(&bytes)));
    }
    Ok(Value::Object(libs))
}
