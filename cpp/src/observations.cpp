#include "irred/observations.hpp"
#include <cmath>
#include <limits>
#include <unordered_set>
#include <utility>
namespace irred::observations {
namespace {
bool digest(const std::string& s) { if(s.size()!=64)return false; for(char c:s)if(!((c>='0'&&c<='9')||(c>='a'&&c<='f')))return false; return true; }
bool mask(const std::vector<std::uint8_t>& v,std::size_t n) { if(v.size()!=n)return false; for(auto x:v)if(x>1)return false;return true; }
}
Prepared prepare(Input in,Policy p) {
 Prepared out; out.source_=std::move(in); const auto& s=out.source_; const auto n=s.values.size();
 if((s.profile!=Profile::pantheon_plus_released_v1&&s.profile!=Profile::gaussian_fixture_v1&&s.profile!=Profile::fits_length_fixture_v1&&s.profile!=Profile::typed_magnitude_covariance)||(s.calibration!=Calibration::unknown&&s.calibration!=Calibration::released_corrected&&s.calibration!=Calibration::not_applicable)||(s.uncertainty!=Uncertainty::none&&s.uncertainty!=Uncertainty::covariance&&s.uncertainty!=Uncertainty::precision)){out.status_=Status::incompatible_semantics;return out;}
 if(s.component!=Component::unknown&&s.component!=Component::statistical&&s.component!=Component::systematic&&s.component!=Component::total){out.status_=Status::incompatible_semantics;return out;}
 const auto expected=s.uncertainty==Uncertainty::none?UncertaintyUnit::none:s.unit==Unit::magnitude?(s.uncertainty==Uncertainty::covariance?UncertaintyUnit::magnitude_squared:UncertaintyUnit::inverse_magnitude_squared):(s.uncertainty==Uncertainty::covariance?UncertaintyUnit::metre_squared:UncertaintyUnit::inverse_metre_squared);
 if(s.uncertainty_unit!=expected){out.status_=Status::incompatible_semantics;return out;}
 if(n==0||n>p.maximum_rows){out.status_=Status::resource_limit;return out;}
 if(s.measurement_ids.size()!=n||s.event_ids.size()!=n||s.quality.size()!=n||!mask(s.missing,n)||(!s.source_selection.empty()&&!mask(s.source_selection,n)))return out;
 if(!digest(s.table_sha256)){out.status_=Status::invalid_identity;return out;}
 std::size_t string_bytes=0;
 auto admit=[&](const std::string& v){if(v.size()>p.maximum_string_bytes-string_bytes)return false;string_bytes+=v.size();return true;};
 if(!admit(s.table_sha256)||!admit(s.uncertainty_sha256)||!admit(s.ordering_provenance)||!admit(s.calibration_provenance)||!admit(s.dependence_provenance)||!admit(s.quality_dictionary)){out.status_=Status::resource_limit;return out;}
 for(const auto* v:{&s.measurement_ids,&s.event_ids,&s.uncertainty_axis_ids})for(const auto& id:*v)if(!admit(id)){out.status_=Status::resource_limit;return out;}
 std::unordered_set<std::string> ids;
 for(std::size_t i=0;i<n;++i)if(s.measurement_ids[i].empty()||s.event_ids[i].empty()||!ids.insert(s.measurement_ids[i]).second){out.status_=Status::invalid_identity;return out;}
 if(s.profile==Profile::pantheon_plus_released_v1){
  if(s.calibration==Calibration::not_applicable||s.role!=Role::released_fitted_summary||s.unit!=Unit::magnitude||s.zhd.size()!=n||s.zcmb.size()!=n||s.zhel.size()!=n||!mask(s.zhd_missing,n)||!mask(s.zcmb_missing,n)||!mask(s.zhel_missing,n)){out.status_=Status::incompatible_semantics;return out;}
 }else if(s.profile==Profile::typed_magnitude_covariance){
  if((s.role!=Role::observed_measurement&&s.role!=Role::released_fitted_summary&&s.role!=Role::synthetic_control)||s.unit!=Unit::magnitude||s.uncertainty!=Uncertainty::covariance){out.status_=Status::incompatible_semantics;return out;}
  for(const auto& column:{std::pair{&s.zhd,&s.zhd_missing},std::pair{&s.zcmb,&s.zcmb_missing},std::pair{&s.zhel,&s.zhel_missing}}){
   if((column.first->empty()&&!column.second->empty())||(!column.first->empty()&&(column.first->size()!=n||!mask(*column.second,n)))){out.status_=Status::incompatible_semantics;return out;}
  }
 }else if(s.role!=Role::synthetic_control||(s.profile==Profile::fits_length_fixture_v1?s.unit!=Unit::metre:s.unit!=Unit::magnitude)||!s.zhd.empty()||!s.zcmb.empty()||!s.zhel.empty()||!s.zhd_missing.empty()||!s.zcmb_missing.empty()||!s.zhel_missing.empty()) {out.status_=Status::incompatible_semantics;return out;}
 if(s.uncertainty==Uncertainty::none){if(s.component!=Component::unknown||!s.uncertainty_matrix.empty()||!s.uncertainty_axis_ids.empty()||!s.uncertainty_sha256.empty())return out;}
 else {
  if(n>std::numeric_limits<std::size_t>::max()/n||n*n>p.maximum_matrix_elements){out.status_=Status::resource_limit;return out;}
  if(s.ordering_provenance.empty()){out.status_=Status::incompatible_semantics;return out;}
  if(s.uncertainty_matrix.size()!=n*n||s.uncertainty_axis_ids!=s.measurement_ids)return out;
  if(!digest(s.uncertainty_sha256)){out.status_=Status::invalid_identity;return out;}
  for(double x:s.uncertainty_matrix)if(!std::isfinite(x)){out.status_=Status::nonfinite_required_value;return out;}
 }
 out.status_=Status::ok;return out;
}
SelectionResult Prepared::select(Selection choice) const {
 SelectionResult r; if(choice!=Selection::all&&choice!=Selection::pantheon_zhd_gt_001){r.status=Status::incompatible_semantics;return r;} if(status_!=Status::ok){r.status=status_;return r;}const auto& s=source_;const auto n=s.values.size();r.mask.assign(n,0);
 auto fail=[&](Status status){r.status=status;r.mask.clear();r.source_indices.clear();return r;};
 if(choice==Selection::pantheon_zhd_gt_001&&s.profile!=Profile::pantheon_plus_released_v1){return fail(Status::incompatible_semantics);}
 for(std::size_t i=0;i<n;++i){
  if(!s.source_selection.empty()&&!s.source_selection[i])continue;
  if(choice==Selection::pantheon_zhd_gt_001){if(s.zhd_missing[i]){return fail(Status::missing_required_value);}if(!std::isfinite(s.zhd[i])){return fail(Status::nonfinite_required_value);}if(!(s.zhd[i]>0.01))continue;}
  if(s.missing[i]||(s.profile==Profile::pantheon_plus_released_v1&&(s.zhd_missing[i]||s.zhel_missing[i]))){return fail(Status::missing_required_value);}
  if(!std::isfinite(s.values[i])||(s.profile==Profile::pantheon_plus_released_v1&&(!std::isfinite(s.zhd[i])||!std::isfinite(s.zhel[i])))){return fail(Status::nonfinite_required_value);}
  r.mask[i]=1;r.source_indices.push_back(i);
 }
 r.status=Status::ok;return r;
}
std::optional<std::size_t> retained_source_payload_bound(const Prepared& prepared) noexcept {
 const auto& s=prepared.source(); std::size_t bytes=0;
 auto add=[&](std::size_t count,std::size_t width){
  if(count>std::numeric_limits<std::size_t>::max()/width)return false;
  const auto n=count*width;
  if(n>std::numeric_limits<std::size_t>::max()-bytes)return false;
  bytes+=n;return true;
 };
 auto string=[&](const std::string& value){
  // Conservatively charge SSO capacity too; implementation-independent upper bound.
  return value.capacity()!=std::numeric_limits<std::size_t>::max()&&add(value.capacity()+1,1);
 };
 auto strings=[&](const std::vector<std::string>& values){
  if(!add(values.capacity(),sizeof(std::string)))return false;
  for(const auto& value:values)if(!string(value))return false;
  return true;
 };
 for(const auto* value:{&s.ordering_provenance,&s.calibration_provenance,&s.dependence_provenance,&s.quality_dictionary,&s.table_sha256,&s.uncertainty_sha256})if(!string(*value))return std::nullopt;
 for(const auto* values:{&s.measurement_ids,&s.event_ids,&s.uncertainty_axis_ids})if(!strings(*values))return std::nullopt;
 for(const auto* values:{&s.values,&s.zhd,&s.zcmb,&s.zhel,&s.uncertainty_matrix})if(!add(values->capacity(),sizeof(double)))return std::nullopt;
 for(const auto* values:{&s.missing,&s.zhd_missing,&s.zcmb_missing,&s.zhel_missing,&s.source_selection})if(!add(values->capacity(),sizeof(std::uint8_t)))return std::nullopt;
 if(!add(s.quality.capacity(),sizeof(std::uint64_t)))return std::nullopt;
 return bytes;
}

}
