mod bao_ingestion;
mod bridge;
mod cli;
mod current_background_run;
mod current_bao_run;
mod current_observation_run;
mod current_session;
mod ingestion;
mod model_spec;
mod outcome;
mod records;
mod retained_context;
mod statistics_run;
mod stream_io;
mod strict_json;
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
