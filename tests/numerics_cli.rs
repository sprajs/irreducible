//! Boundary/status/record evidence only; core-author test, no Rust numerical engine.
use serde_json::{json,Value};
use std::{fs,path::PathBuf,process::Command,time::{SystemTime,UNIX_EPOCH}};
struct Scratch(PathBuf);
impl Drop for Scratch { fn drop(&mut self){let _=fs::remove_dir_all(&self.0);} }
fn run_bytes(bytes:&[u8])->(Value,i32,Option<Value>){
 let scratch=Scratch(std::env::temp_dir().join(format!("irred-numerical-cli-{}-{}",std::process::id(),SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_nanos())));
 fs::create_dir(&scratch.0).unwrap();let request=scratch.0.join("request.json");let store=scratch.0.join("store");fs::write(&request,bytes).unwrap();
 let out=Command::new(env!("CARGO_BIN_EXE_irred")).arg("run").arg(&request).arg(&store).output().unwrap();
 let response:Value=serde_json::from_slice(&out.stdout).unwrap();
 let resolved=response["receipt"]["scientific_specification_digest"].as_str().map(|d|serde_json::from_slice(&fs::read(store.join("objects").join(d)).unwrap()).unwrap());
 (response,out.status.code().unwrap(),resolved)
}
fn run(method:&str,values:Value)->(Value,i32,Value){
 let spec=json!({"schema_version":1,"operation":"numerics.scalar_batch.v1","method":method,"values":values});
 let (response,exit,resolved)=run_bytes(&serde_json::to_vec(&spec).unwrap());let resolved=resolved.unwrap();
 assert_eq!(resolved["method"],method);assert_eq!(resolved["values"],spec["values"]);assert_eq!(resolved["equation_id"],format!("F02/scalar/{method}/v1"));
 assert_eq!(resolved["requested_outputs"],json!([{"id":"evaluations","required":true,"numerical_gate":"required","inference_gate":"not_applicable"}]));
 assert_eq!(response["result"]["source_values"],spec["values"]);
 assert_eq!(response["receipt"]["backend"],"portable_cpu");assert_eq!(response["receipt"]["resource_budget"]["max_batch_elements"],1000000);
 assert_eq!(response["receipt"]["precision"],"binary64_storage_method_declared_intermediate");
 (response,exit,resolved)
}
fn finite(response:&Value,i:usize)->f64 {let s=&response["result"]["evaluations"][i];assert_eq!(s["result"]["kind"],"finite");assert_eq!(s["diagnostics"]["error_estimate_kind"],"empirical_or_arithmetic_diagnostic_not_certified_bound");s["result"]["value"].as_f64().unwrap()}
fn failed(response:&Value,i:usize,status:&str){let s=&response["result"]["evaluations"][i]["result"];assert_eq!(s["kind"],"failure");assert_eq!(s["error_id"],status);assert!(s.get("value").is_none(),"failed row must not expose finite sentinel");}
#[test]
fn scalar_status_values_and_unqualified_gate(){
 let (r,exit,_)=run("compensated_sum",json!([1e16,1.,-1e16]));assert_eq!(finite(&r,0),1.);assert_eq!(exit,6);assert_eq!(r["receipt"]["execution"],"completed");assert_eq!(r["receipt"]["accepted"],false);assert_eq!(r["receipt"]["outputs"][0]["numerical"],"not_assessed");assert_eq!(r["result"]["evaluations"][0]["diagnostics"]["evaluations"],3);
 let(r,exit,_)=run("log1p",json!([0.,0.5,-1.,-2.]));assert_ne!(exit,0);assert_eq!(finite(&r,0),0.);let expected=0.40546510810816438197801311546434913657199042346249_f64;assert!((finite(&r,1)-expected).abs()<1e-15);failed(&r,2,"outside_domain");failed(&r,3,"outside_domain");assert_eq!(r["receipt"]["execution"],"failed");assert_eq!(r["receipt"]["accepted"],false);
 let(r,exit,_)=run("expm1",json!([0.,1e-12,1000.]));assert_ne!(exit,0);assert_eq!(finite(&r,0),0.);failed(&r,2,"overflow");
 let(r,_,_)=run("log_gamma_positive",json!([0.5,5.,0.]));let expected=0.57236494292470008707171367567652935582364740645766_f64;assert!((finite(&r,0)-expected).abs()<1e-14);failed(&r,2,"outside_domain");
 let(r,exit,_)=run("log_sum_exp",json!([1000.,1000.]));assert_eq!(exit,6);assert!((finite(&r,0)-1000.6931471805599).abs()<1e-12);
}
#[test]
fn empty_shapes_and_failed_reduction(){
 let(r,exit,_)=run("compensated_sum",json!([]));assert_eq!(exit,6);assert_eq!(finite(&r,0),0.);
 let(r,exit,_)=run("log_sum_exp",json!([]));assert_ne!(exit,0);failed(&r,0,"invalid_input");
 for method in ["log1p","expm1","log_gamma_positive"] {let(r,exit,_)=run(method,json!([]));assert_eq!(exit,6);assert_eq!(r["result"]["evaluations"],json!([]));}
 let(r,exit,_)=run("compensated_sum",json!([f64::MAX,f64::MAX]));assert_ne!(exit,0);failed(&r,0,"overflow");
}
#[test]
fn malformed_scalar_specifications_fail_closed(){
 for spec in [json!({"schema_version":2,"operation":"numerics.scalar_batch.v1","method":"log1p","values":[0.]}),json!({"schema_version":1,"operation":"numerics.scalar_batch.v1","method":"invented","values":[0.]}),json!({"schema_version":1,"operation":"numerics.scalar_batch.v1","method":"log1p","values":[0.],"ignored":true}),json!({"schema_version":1,"operation":"numerics.scalar_batch.v1","method":"log1p","values":["0"]})] {
  let(r,exit,_)=run_bytes(&serde_json::to_vec(&spec).unwrap());assert_ne!(exit,0);assert_eq!(r["result"]["kind"],"failure");assert_eq!(r["receipt"]["execution"],"failed");
 }
 for bytes in [br#"{"schema_version":1,"operation":"numerics.scalar_batch.v1","method":"log1p","method":"expm1","values":[0]}"#.as_slice(),br#"{"schema_version":1,"operation":"numerics.scalar_batch.v1","method":"log1p","values":[NaN]}"#.as_slice(),br#"{"schema_version":1,"operation":"numerics.scalar_batch.v1","method":"log1p","values":[1e999]}"#.as_slice()] {
  let(r,exit,_)=run_bytes(bytes);assert_ne!(exit,0);assert_eq!(r["result"]["kind"],"failure");
 }
}
