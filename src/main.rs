mod bridge;
mod cli;
mod ingestion;
mod observation_run;
mod records;
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
