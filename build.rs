fn main() {
    let root = std::env::var("CARGO_MANIFEST_DIR").unwrap();
    let profile = std::env::var("PROFILE").unwrap();
    let (native, manifest) = match profile.as_str() {
        "debug" => ("build/native", "build/build-manifest.json"),
        "release" => ("build/native-release", "build/build-manifest-release.json"),
        _ => panic!("unsupported build profile; use tools/build.py"),
    };
    println!("cargo:rustc-link-search=native={root}/{native}");
    println!("cargo:rustc-link-lib=static=irred_core");
    println!("cargo:rustc-link-lib=dylib=stdc++");
    println!("cargo:rustc-env=IRRED_BUILD_MANIFEST={root}/{manifest}");
    println!("cargo:rerun-if-changed={manifest}");
    println!("cargo:rerun-if-changed={native}/libirred_core.a");
}
