mod gaussian_predictive_run;
mod gaussian_input;
mod gaussian_posterior_run;
mod sound_horizon_run;
mod bao_ingestion;
mod bridge;
mod cli;
mod background_run;
mod photometry_run;
mod sampled_photometry_run;
mod bao_run;
mod observation_run;
mod session;
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
