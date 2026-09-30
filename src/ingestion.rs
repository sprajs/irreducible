//! Structural, compile-selected local readers. Scientific selection belongs to C++.
use sha2::{Digest, Sha256};
use std::{collections::HashMap, fs, path::Path};
#[derive(Debug, Clone)]
pub struct Asset {
    sha256: String,
    bytes: Vec<u8>,
}
impl Asset {
    pub fn sha256(&self) -> &str {
        &self.sha256
    }
    pub fn bytes(&self) -> &[u8] {
        &self.bytes
    }
    pub fn from_bytes(bytes: Vec<u8>, maximum_bytes: usize) -> Result<Self, String> {
        if bytes.len() > maximum_bytes {
            return Err("RESOURCE_LIMIT".into());
        }
        Ok(Self {
            sha256: format!("{:x}", Sha256::digest(&bytes)),
            bytes,
        })
    }
}
#[derive(Debug, Clone)]
pub struct Table {
    pub asset: Asset,
    pub measurement_ids: Vec<String>,
    pub event_ids: Vec<String>,
    pub values: Vec<f64>,
    pub zhd: Vec<f64>,
    pub zcmb: Vec<f64>,
    pub zhel: Vec<f64>,
    pub quality: Vec<u64>,
    pub original_columns: Vec<String>,
    pub original_fields: Vec<Vec<String>>,
}
#[derive(Debug, Clone)]
pub struct Matrix {
    pub ordering_provenance: String,
    pub asset: Asset,
    pub values: Vec<f64>,
    pub dimension: usize,
    pub axis_ids: Vec<String>,
}
pub fn read_asset(path: &Path, maximum_bytes: usize) -> Result<Asset, String> {
    use std::io::Read;
    let f = fs::File::open(path).map_err(|_| "INPUT_OPEN")?;
    let mut bytes = Vec::new();
    f.take(
        (maximum_bytes as u64)
            .checked_add(1)
            .ok_or("RESOURCE_LIMIT")?,
    )
    .read_to_end(&mut bytes)
    .map_err(|_| "INPUT_READ")?;
    if bytes.len() > maximum_bytes {
        return Err("RESOURCE_LIMIT".into());
    }
    Asset::from_bytes(bytes, maximum_bytes)
}
fn number(s: &str) -> Result<f64, String> {
    s.parse::<f64>().map_err(|_| "INVALID_NUMBER".into())
}
pub fn pantheon_plus(asset: Asset, maximum_rows: usize) -> Result<Table, String> {
    let text = std::str::from_utf8(&asset.bytes).map_err(|_| "INVALID_UTF8")?;
    let mut lines = text.lines().filter(|x| !x.trim().is_empty());
    let header = lines
        .next()
        .ok_or("MISSING_HEADER")?
        .trim()
        .trim_start_matches('#');
    let columns: Vec<String> = header.split_whitespace().map(str::to_owned).collect();
    let mut map = HashMap::new();
    for (i, k) in columns.iter().enumerate() {
        if map.insert(k.as_str(), i).is_some() {
            return Err("DUPLICATE_COLUMN".into());
        }
    }
    let required = ["CID", "IDSURVEY", "zHD", "zCMB", "zHEL", "m_b_corr"];
    let positions: Vec<usize> = required
        .iter()
        .map(|k| map.get(k).copied().ok_or("MISSING_COLUMN"))
        .collect::<Result<_, _>>()?;
    let mut rows = Vec::new();
    let mut events = Vec::new();
    let mut vals = Vec::new();
    let mut hd = Vec::new();
    let mut cmb = Vec::new();
    let mut hel = Vec::new();
    let mut fields = Vec::new();
    for line in lines {
        if line.trim_start().starts_with('#') {
            continue;
        }
        if rows.len() >= maximum_rows {
            return Err("RESOURCE_LIMIT".into());
        }
        let f: Vec<String> = line.split_whitespace().map(str::to_owned).collect();
        if f.len() != columns.len() {
            return Err("ROW_WIDTH".into());
        }
        let i = rows.len();
        rows.push(format!("{}:row:{}", asset.sha256, i));
        events.push(f[positions[0]].clone());
        // IDSURVEY remains in original_fields; CID is event identity, not measurement identity.
        hd.push(number(&f[positions[2]])?);
        cmb.push(number(&f[positions[3]])?);
        hel.push(number(&f[positions[4]])?);
        vals.push(number(&f[positions[5]])?);
        fields.push(f);
    }
    if rows.is_empty() {
        return Err("EMPTY_TABLE".into());
    }
    let n = rows.len();
    Ok(Table {
        asset,
        measurement_ids: rows,
        event_ids: events,
        values: vals,
        zhd: hd,
        zcmb: cmb,
        zhel: hel,
        quality: vec![0; n],
        original_columns: columns,
        original_fields: fields,
    })
}
pub fn pantheon_covariance(
    asset: Asset,
    axis_ids: Vec<String>,
    maximum_elements: usize,
) -> Result<Matrix, String> {
    let text = std::str::from_utf8(&asset.bytes).map_err(|_| "INVALID_UTF8")?;
    let mut tokens = text.split_whitespace();
    let n = tokens
        .next()
        .ok_or("MISSING_DIMENSION")?
        .parse::<usize>()
        .map_err(|_| "INVALID_DIMENSION")?;
    let count = n.checked_mul(n).ok_or("RESOURCE_LIMIT")?;
    if n == 0 || count > maximum_elements {
        return Err("RESOURCE_LIMIT".into());
    }
    if axis_ids.len() != n {
        return Err("AXIS_SHAPE".into());
    }
    let mut values = Vec::with_capacity(count);
    for _ in 0..count {
        values.push(number(tokens.next().ok_or("TRUNCATED_MATRIX")?)?);
    }
    if tokens.next().is_some() {
        return Err("TRAILING_MATRIX".into());
    }
    // Preserve source asymmetry/precision; interpretation and normalization are native policies.
    Ok(Matrix {
        ordering_provenance: "supplied axis IDs in original unlabeled row-major source order; declaration not independently verified".into(),
        asset,
        values,
        dimension: n,
        axis_ids,
    })
}
/// Synthetic magnitude fixture: exact ROW_ID EVENT_ID MAG header. This is a
/// separately identified control, never a substituted survey release.
pub fn magnitude_table(asset: Asset, maximum_rows: usize) -> Result<Table, String> {
    let text = std::str::from_utf8(&asset.bytes).map_err(|_| "INVALID_UTF8")?;
    let mut lines = text.lines().filter(|x| !x.trim().is_empty());
    if lines.next().ok_or("MISSING_HEADER")?.trim() != "ROW_ID EVENT_ID MAG" {
        return Err("WRONG_PROFILE".into());
    }
    let mut ids = Vec::new();
    let mut events = Vec::new();
    let mut values = Vec::new();
    let mut fields = Vec::new();
    for line in lines {
        if ids.len() >= maximum_rows {
            return Err("RESOURCE_LIMIT".into());
        }
        let f: Vec<String> = line.split_whitespace().map(str::to_owned).collect();
        if f.len() != 3 {
            return Err("ROW_WIDTH".into());
        }
        ids.push(f[0].clone());
        events.push(f[1].clone());
        values.push(number(&f[2])?);
        fields.push(f);
    }
    if ids.is_empty() {
        return Err("EMPTY_TABLE".into());
    }
    let n = ids.len();
    Ok(Table {
        asset,
        measurement_ids: ids,
        event_ids: events,
        values,
        zhd: Vec::new(),
        zcmb: Vec::new(),
        zhel: Vec::new(),
        quality: vec![0; n],
        original_columns: vec!["ROW_ID".into(), "EVENT_ID".into(), "MAG".into()],
        original_fields: fields,
    })
}
