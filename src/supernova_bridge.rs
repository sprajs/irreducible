use super::generated::*;
type WireModel = SupernovaModel;
type WireBatch = SupernovaBatch;
type WireSlot = SupernovaSlot;
use super::generated::{
    cosmo_supernova_evaluate as evaluate_ffi, cosmo_supernova_result_destroy as destroy_ffi,
    cosmo_supernova_result_view as view_ffi,
};
fn model_tag(label: &str) -> Option<u32> {
    super::generated::background_tag_id("model", label)
}
define_supernova_bridge!();
