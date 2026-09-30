//! Narrow structural BAO adapters. No distance, ruler or probability equations.
use crate::ingestion::{Asset, read_asset};
use serde::{Deserialize, Serialize};
use std::path::Path;
#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct Row {
    pub id: String,
    pub z: f64,
    pub observable: String,
    pub value: f64,
}
#[derive(Deserialize, Serialize)]
#[serde(tag = "profile", deny_unknown_fields)]
pub(crate) enum Source {
    #[serde(rename = "desi_dr2_all_gccomb_13_v1")]
    Released {
        mean: String,
        covariance: String,
        maximum_asset_bytes: usize,
        ordering_provenance: String,
        calibration_provenance: String,
        dependence_provenance: String,
        redshift_convention: String,
        ruler_convention: String,
        computational_h0_convention: String,
    },
    #[serde(rename = "synthetic_inline_v1")]
    Synthetic {
        rows: Vec<Row>,
        covariance: Vec<f64>,
        covariance_axis_ids: Vec<String>,
        ordering_provenance: String,
        calibration_provenance: String,
        dependence_provenance: String,
        redshift_convention: String,
        ruler_convention: String,
        computational_h0_convention: String,
    },
}
pub(crate) struct Decoded {
    pub rows: Vec<Row>,
    pub covariance: Vec<f64>,
    pub covariance_axis_ids: Vec<String>,
    pub table_asset: Asset,
    pub covariance_asset: Asset,
    pub original_fields: Vec<Vec<String>>,
    pub role: &'static str,
    pub ordering_provenance: String,
    pub calibration_provenance: String,
    pub dependence_provenance: String,
    pub redshift_convention: String,
    pub ruler_convention: String,
    pub computational_h0_convention: String,
}
fn observable(s: &str) -> Result<(), String> {
    match s {
        "DM_over_rs" | "DH_over_rs" | "DV_over_rs" => Ok(()),
        _ => Err("UNSUPPORTED_BAO_OBSERVABLE".into()),
    }
}
fn lines(a: &Asset, rows: usize, columns: usize) -> Result<Vec<Vec<String>>, String> {
    let text = std::str::from_utf8(a.bytes()).map_err(|_| "INVALID_UTF8")?;
    let mut result = Vec::with_capacity(rows);
    for line in text
        .lines()
        .map(str::trim)
        .filter(|s| !s.is_empty() && !s.starts_with('#'))
    {
        if result.len() == rows {
            return Err("INVALID_BAO_SOURCE_SHAPE".into());
        }
        let mut fields = Vec::with_capacity(columns);
        for token in line.split_whitespace() {
            if fields.len() == columns {
                return Err("INVALID_BAO_SOURCE_SHAPE".into());
            }
            fields.push(token.to_owned());
        }
        if fields.len() != columns {
            return Err("INVALID_BAO_SOURCE_SHAPE".into());
        }
        result.push(fields);
    }
    if result.len() != rows {
        return Err("INVALID_BAO_SOURCE_SHAPE".into());
    }
    Ok(result)
}
impl Source {
    pub(crate) fn decode(
        &self,
        maximum_rows: usize,
        maximum_matrix_elements: usize,
    ) -> Result<Decoded, String> {
        let (
            rows,
            covariance,
            axis,
            table,
            matrix,
            original,
            role,
            ordering,
            calibration,
            dependence,
            redshift,
            ruler,
            h0,
        ) = match self {
            Self::Released {
                mean,
                covariance,
                maximum_asset_bytes,
                ordering_provenance,
                calibration_provenance,
                dependence_provenance,
                redshift_convention,
                ruler_convention,
                computational_h0_convention,
            } => {
                if maximum_rows < 13
                    || maximum_matrix_elements < 169
                    || *maximum_asset_bytes > 16 * 1024 * 1024
                {
                    return Err("RESOURCE_LIMIT".into());
                }
                let table = read_asset(Path::new(mean), *maximum_asset_bytes)?;
                let matrix = read_asset(Path::new(covariance), *maximum_asset_bytes)?;
                let original = lines(&table, 13, 3).map_err(|_| "INVALID_BAO_RELEASE_SHAPE")?;
                if original.len() != 13 || original.iter().any(|r| r.len() != 3) {
                    return Err("INVALID_BAO_RELEASE_SHAPE".into());
                }
                let rows = original
                    .iter()
                    .enumerate()
                    .map(|(i, r)| {
                        observable(&r[2])?;
                        Ok(Row {
                            id: format!("{}:row:{i}:{}", table.sha256(), r[2]),
                            z: r[0].parse().map_err(|_| "INVALID_NUMBER")?,
                            value: r[1].parse().map_err(|_| "INVALID_NUMBER")?,
                            observable: r[2].clone(),
                        })
                    })
                    .collect::<Result<Vec<_>, String>>()?;
                let values = lines(&matrix, 13, 13)
                    .map_err(|_| "INVALID_BAO_COVARIANCE_SHAPE")?
                    .into_iter()
                    .flatten()
                    .map(|s| s.parse::<f64>().map_err(|_| "INVALID_NUMBER".to_string()))
                    .collect::<Result<Vec<_>, _>>()?;
                if values.len() != 169 {
                    return Err("INVALID_BAO_COVARIANCE_SHAPE".into());
                }
                let axis = rows.iter().map(|r| r.id.clone()).collect();
                (
                    rows,
                    values,
                    axis,
                    table,
                    matrix,
                    original,
                    "released_fitted_distance_summary",
                    ordering_provenance,
                    calibration_provenance,
                    dependence_provenance,
                    redshift_convention,
                    ruler_convention,
                    computational_h0_convention,
                )
            }
            Self::Synthetic {
                rows,
                covariance,
                covariance_axis_ids,
                ordering_provenance,
                calibration_provenance,
                dependence_provenance,
                redshift_convention,
                ruler_convention,
                computational_h0_convention,
            } => {
                let count = rows.len();
                let square = count.checked_mul(count).ok_or("RESOURCE_LIMIT")?;
                if count == 0
                    || count > maximum_rows
                    || square > maximum_matrix_elements
                    || covariance.len() != square
                    || covariance_axis_ids.len() != count
                {
                    return Err("INVALID_BAO_SOURCE_SHAPE".into());
                }
                for row in rows {
                    observable(&row.observable)?;
                }
                let table_bytes = serde_json::to_vec(rows).map_err(|e| e.to_string())?;
                let covariance_bytes = serde_json::to_vec(covariance).map_err(|e| e.to_string())?;
                let table = Asset::from_bytes(table_bytes, 16 * 1024 * 1024)?;
                let matrix = Asset::from_bytes(covariance_bytes, 128 * 1024 * 1024)?;
                let copied =
                    serde_json::from_slice::<Vec<Row>>(table.bytes()).map_err(|e| e.to_string())?;
                (
                    copied,
                    covariance.clone(),
                    covariance_axis_ids.clone(),
                    table,
                    matrix,
                    Vec::new(),
                    "synthetic_control",
                    ordering_provenance,
                    calibration_provenance,
                    dependence_provenance,
                    redshift_convention,
                    ruler_convention,
                    computational_h0_convention,
                )
            }
        };
        Ok(Decoded {
            rows,
            covariance,
            covariance_axis_ids: axis,
            table_asset: table,
            covariance_asset: matrix,
            original_fields: original,
            role,
            ordering_provenance: ordering.clone(),
            calibration_provenance: calibration.clone(),
            dependence_provenance: dependence.clone(),
            redshift_convention: redshift.clone(),
            ruler_convention: ruler.clone(),
            computational_h0_convention: h0.clone(),
        })
    }
}
