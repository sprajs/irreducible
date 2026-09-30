//! Current compiled hypotheses and projection inputs. No inactive parameters,
//! runtime expressions or source-acquisition logic live here.
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(tag = "kind", deny_unknown_fields)]
pub(crate) enum ExpansionSpec {
    #[serde(rename = "lcdm")]
    Lcdm { omega_m: f64 },
    #[serde(rename = "constant_q")]
    ConstantQ { q: f64 },
    #[serde(rename = "cpl")]
    Cpl { omega_m: f64, w0: f64, wa: f64 },
    #[serde(rename = "fixed_q5")]
    FixedQ5 { q: [f64; 5] },
}

#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
#[serde(tag = "kind", deny_unknown_fields)]
pub(crate) enum Geometry {
    #[serde(rename = "flat_flrw")]
    FlatFlrw,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(tag = "kind", deny_unknown_fields)]
pub(crate) enum SourceEffect {
    #[serde(rename = "none")]
    None,
    #[serde(rename = "grey_log1p_magnitude")]
    GreyLog1pMagnitude { epsilon_mag: f64 },
}

#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum ObserverConvention {
    GeometricSameRedshift,
    ReleasedZhdZhel,
}
#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct Observer {
    pub redshift: f64,
    pub convention: ObserverConvention,
}
#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct PhysicalScale {
    pub h0_km_s_mpc: f64,
}

#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum BackgroundOutput {
    Expansion,
    Radial,
    LuminosityShape,
    Clock,
    Physical,
    Kinematics,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct BackgroundRequest {
    pub z_expansion: f64,
    pub requested_outputs: Vec<BackgroundOutput>,
    pub observer: Option<Observer>,
    pub physical_scale: Option<PhysicalScale>,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct SupernovaModel {
    pub expansion: ExpansionSpec,
    pub geometry: Geometry,
    pub source_effect: SourceEffect,
}
