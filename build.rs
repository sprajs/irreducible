fn main() {
    let root = std::env::var("CARGO_MANIFEST_DIR").unwrap();
    println!("cargo:rustc-link-search=native={root}/build/native");
    println!("cargo:rustc-link-lib=static=cosmology_core");
    println!("cargo:rustc-link-lib=dylib=stdc++");
    println!("cargo:rerun-if-changed=build/build-manifest.json");
    println!("cargo:rerun-if-changed=build/native/libcosmology_core.a");
}
