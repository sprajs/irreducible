mod background_run;
mod background_v2_run;
mod bao_ingestion;
mod bao_run;
mod bridge;
mod cli;
mod ingestion;
mod observation_run;
mod records;
mod piecewise_run;
mod statistics_run;
mod supernova_run;
mod supernova_v2_run;
mod supernova_piecewise_run;
fn main() {
    if let Err(e) = cli::execute() {
        eprintln!("{}", serde_json::json!({"kind":"failure","error_id":e}));
        std::process::exit(if e == "NUMERICAL_QUALIFICATION_REQUIRED" {
            6
        } else {
            2
        })
    }
}
