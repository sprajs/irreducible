mod bridge;
mod cli;
mod records;
fn main() {
    if let Err(e) = cli::execute() {
        eprintln!("{}", serde_json::json!({"kind":"failure","error_id":e}));
        std::process::exit(2)
    }
}
