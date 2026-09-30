use super::generated::*;
type WireParameters = BackgroundParametersV2;
type WireBatch = BackgroundBatchV2;
type WireSlot = BackgroundSlotV2;
use super::generated::{
    cosmo_background_evaluate_v2 as evaluate_ffi,
    cosmo_background_result_v2_destroy as destroy_ffi, cosmo_background_result_v2_view as view_ffi,
};
fn model_tag(label: &str) -> Option<u32> {
    super::generated::background_v2_tag_id("model", label)
}

define_background_bridge!(w0, wa);
