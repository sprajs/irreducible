use super::generated::*;
type WireModel = SupernovaModelV2;
type WireBatch = SupernovaBatchV2;
type WireSlot = SupernovaSlotV2;
use super::generated::{
    cosmo_supernova_evaluate_v2 as evaluate_ffi, cosmo_supernova_result_v2_destroy as destroy_ffi,
    cosmo_supernova_result_v2_view as view_ffi,
};
fn model_tag(label: &str) -> Option<u32> {
    super::generated::background_v2_tag_id("model", label)
}

define_supernova_bridge!(w0, wa);
