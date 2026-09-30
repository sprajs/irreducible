#include "irred/supernova.hpp"
namespace irred::supernova {
SelectedMagnitudeSource
pantheon_zhd_gt_001(std::shared_ptr<const observations::Prepared> p) {
  SelectedMagnitudeSource out;
  out.source = std::move(p);
  if (!out.source || out.source->status() != observations::Status::ok)
    return out;
  const auto &s = out.source->source();
  if (s.profile != observations::Profile::pantheon_plus_released_v1 ||
      s.role != observations::Role::released_fitted_summary)
    return out;
  auto selection =
      out.source->select(observations::Selection::pantheon_zhd_gt_001);
  if (selection.status != observations::Status::ok)
    return out;
  out.source_indices = std::move(selection.source_indices);
  out.ordered_ids.reserve(out.source_indices.size());
  out.coordinates.reserve(out.source_indices.size());
  for (auto i : out.source_indices) {
    out.ordered_ids.push_back(s.measurement_ids[i]);
    out.coordinates.push_back(
        {s.zhd[i], {s.zhel[i], cosmology::Convention::released_zhd_zhel}});
  }
  return out;
}
} // namespace irred::supernova
