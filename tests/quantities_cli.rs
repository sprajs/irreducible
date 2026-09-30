//! Independent CLI contract for C++-owned quantity conversions; no Rust unit algebra.
use serde_json::{json,Value};
use std::{fs,path::PathBuf,process::Command,time::{SystemTime,UNIX_EPOCH}};
struct Scratch(PathBuf);impl Drop for Scratch{fn drop(&mut self){let _=fs::remove_dir_all(&self.0);}}
fn metadata(unit:&str,role:&str,frame:&str,convention:&str)->Value {json!({"unit":unit,"role":role,"frame":frame,"convention":convention,"constant_set":"SI-IAU-definitions-v1"})}
fn run(values:Value,source:Value,target:Value)->(Value,i32){
 let scratch=Scratch(std::env::temp_dir().join(format!("cosmology-native-quantity-{}-{}",std::process::id(),SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_nanos())));fs::create_dir(&scratch.0).unwrap();let request=scratch.0.join("request.json");let store=scratch.0.join("store");
 let spec=json!({"schema_version":1,"operation":"quantity.convert.v1","values":values,"source":source,"target":target});fs::write(&request,serde_json::to_vec(&spec).unwrap()).unwrap();let out=Command::new(env!("CARGO_BIN_EXE_irred")).arg("run").arg(request).arg(&store).output().unwrap();let response:Value=serde_json::from_slice(&out.stdout).unwrap();
 let digest=response["receipt"]["scientific_specification_digest"].as_str().unwrap();let resolved:Value=serde_json::from_slice(&fs::read(store.join("objects").join(digest)).unwrap()).unwrap();assert_eq!(resolved["source"],spec["source"]);assert_eq!(resolved["target"],spec["target"]);assert_eq!(resolved["values"].as_array().unwrap().iter().map(|v|v.as_f64().unwrap()).collect::<Vec<_>>(),spec["values"].as_array().unwrap().iter().map(|v|v.as_f64().unwrap()).collect::<Vec<_>>());assert_eq!(resolved["equation_id"],"EQ-quantity-conversion-v1");(response,out.status.code().unwrap())
}
#[test]
fn quantity_values_statuses_and_metadata() {
 let km=metadata("kilometre","physical_length","none","physical");let m=metadata("metre","physical_length","none","physical");
 let (r,exit)=run(json!([1,-1,0]),km.clone(),m.clone());assert_eq!(r["result"]["source"],km);assert_eq!(r["result"]["target"],m);assert_eq!(r["result"]["source_values"],json!([1.,-1.,0.]));assert_eq!(r["result"]["evaluations"],json!([{"kind":"finite","value":1000.},{"kind":"finite","value":-1000.},{"kind":"finite","value":0.}]));assert_eq!(r["receipt"]["execution"],"completed");
 // The unqualified bootstrap CLI withholds accepted science even if arithmetic completes.
 if r["receipt"]["outputs"][0]["numerical"]=="not_assessed" {assert_eq!(r["receipt"]["accepted"],false);assert_eq!(exit,6);}
 let lum=metadata("metre","luminosity_distance","none","physical");let (r,exit)=run(json!([1,-1,0]),lum.clone(),lum);assert_ne!(exit,0);assert_eq!(r["receipt"]["execution"],"failed");assert_eq!(r["result"]["evaluations"][0]["value"],1.);assert_eq!(r["result"]["evaluations"][1]["kind"],"failure");assert_eq!(r["result"]["evaluations"][1]["error_id"],"invalid_domain");assert!(r["result"]["evaluations"][1].get("value").is_none());assert_eq!(r["result"]["evaluations"][2]["value"],0.);
 let z=metadata("one","redshift","cmb","none");let zh=metadata("one","redshift","heliocentric","none");let(r,exit)=run(json!([0.1]),z,zh);assert_ne!(exit,0);assert_eq!(r["result"]["evaluations"][0]["error_id"],"unsupported_transform");assert!(r["result"]["evaluations"][0].get("value").is_none());
 let(r,_)=run(json!([]),metadata("kilometre","physical_length","none","physical"),m);assert_eq!(r["result"]["evaluations"],json!([]));
}
#[test]
fn frozen_independent_constant_values() {
 let(r,_)=run(json!([1]),metadata("parsec","physical_length","none","physical"),metadata("metre","physical_length","none","physical"));let pc=r["result"]["evaluations"][0]["value"].as_f64().unwrap();let expected=30856775814913672.789139379577964716107319211604092_f64;assert!(pc.to_bits().abs_diff(expected.to_bits())<=1,"independent Machin/Decimal fixture");
 let(r,_)=run(json!([70]),metadata("km_per_s_per_mpc","expansion_rate","none","none"),metadata("inverse_second","expansion_rate","none","none"));let h=r["result"]["evaluations"][0]["value"].as_f64().unwrap();let expected=2.2685455026110555162663814014470700829931485893862e-18_f64;assert!(h.to_bits().abs_diff(expected.to_bits())<=2,"independent H0 Decimal fixture");
}
