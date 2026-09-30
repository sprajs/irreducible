//! Operation-owned execution and assurance metadata; no scientific equations.
use serde::Serialize;
use serde_json::{Value, json};
#[derive(Clone, Copy, Serialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum NumericalCheck {
    ChecksPassed,
    Failed,
}
#[derive(Serialize)]
pub(crate) struct OutputCheck {
    pub check_kind: &'static str,
    pub id: &'static str,
    pub required: bool,
    pub numerical: NumericalCheck,
    pub inference: &'static str,
    pub interpretation: &'static str,
    pub evidence: Vec<Value>,
    pub validation_coverage: &'static str,
}
#[derive(Serialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum ExecutionDisposition {
    Completed,
}
pub(crate) struct Outcome {
    pub execution: ExecutionDisposition,
    pub specification: Value,
    pub output: Value,
    pub method: Value,
    pub arithmetic: Value,
    pub resources: Value,
    pub outputs: Vec<OutputCheck>,
    pub scientific: bool,
}
impl Outcome {
    pub fn scientific(
        specification: Value,
        output: Value,
        id: &'static str,
        method: Value,
        arithmetic: Value,
        resources: Value,
        coverage: &'static str,
    ) -> Self {
        let numerical = if output["kind"] == "finite" {
            NumericalCheck::ChecksPassed
        } else {
            NumericalCheck::Failed
        };
        Self {
            execution: ExecutionDisposition::Completed,
            specification,
            output,
            method,
            arithmetic,
            resources,
            scientific: true,
            outputs: vec![OutputCheck {
                check_kind: "numerical_contract",
                id,
                required: true,
                numerical,
                inference: "not_applicable",
                interpretation: "unqualified",
                evidence: vec![],
                validation_coverage: coverage,
            }],
        }
    }
    pub fn exact(specification: Value, output: Value) -> Self {
        Self {
            execution: ExecutionDisposition::Completed,
            specification,
            output,
            method: json!("checked_i64_add"),
            arithmetic: json!("exact_i64"),
            resources: json!({}),
            scientific: false,
            outputs: vec![OutputCheck {
                check_kind: "exact_integer_contract",
                id: "sum",
                required: true,
                numerical: NumericalCheck::ChecksPassed,
                inference: "not_applicable",
                interpretation: "not_applicable",
                evidence: vec![],
                validation_coverage: "checked integer contract",
            }],
        }
    }
    pub fn structural(
        specification: Value,
        output: Value,
        id: &'static str,
        method: Value,
        arithmetic: Value,
        resources: Value,
        coverage: &'static str,
    ) -> Self {
        let mut out = Self::scientific(
            specification,
            output,
            id,
            method,
            arithmetic,
            resources,
            coverage,
        );
        for check in &mut out.outputs {
            check.check_kind = "structural_source_contract";
        }
        out
    }
    pub fn numerical_passed(&self) -> bool {
        self.outputs
            .iter()
            .filter(|o| o.required)
            .all(|o| matches!(o.numerical, NumericalCheck::ChecksPassed))
    }
    pub fn qualification_passed(&self) -> bool {
        !self.scientific
    }
}
