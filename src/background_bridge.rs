use super::generated::*;
type WireParameters = BackgroundParameters;
type WireBatch = BackgroundBatch;
type WireSlot = BackgroundSlot;
use super::generated::{
    cosmo_background_evaluate as evaluate_ffi, cosmo_background_result_destroy as destroy_ffi,
    cosmo_background_result_view as view_ffi,
};
fn model_tag(label: &str) -> Option<u32> {
    super::generated::background_tag_id("model", label)
}
define_background_bridge!();
