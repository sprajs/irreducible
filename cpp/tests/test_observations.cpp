#include "cosmology/observations.hpp"
#include <stdexcept>
#define CHECK(...) do { if(!(__VA_ARGS__)) throw std::runtime_error("contract failed: " #__VA_ARGS__); } while(false)
#include <cmath>
#include <limits>
using namespace cosmology::observations;
Input fixture(){Input x{};x.profile=Profile::pantheon_plus_released_v1;x.role=Role::released_fitted_summary;x.unit=Unit::magnitude;x.calibration=Calibration::unknown;x.uncertainty=Uncertainty::covariance;x.uncertainty_unit=UncertaintyUnit::magnitude_squared;x.table_sha256=std::string(64,'a');x.uncertainty_sha256=std::string(64,'b');x.measurement_ids={"a:0","a:1","a:2"};x.event_ids={"A","A","B"};x.uncertainty_axis_ids=x.measurement_ids;x.ordering_provenance="supplied original release row order; not independently verified";x.values={17,18,19};x.zhd={.009,.02,.01};x.zcmb={.02,.03,.04};x.zhel={.03,.04,.05};x.missing=x.zhd_missing=x.zcmb_missing=x.zhel_missing={0,0,0};x.quality={0,2,128};x.uncertainty_matrix={4,1.00000003,0,1,9,0,0,0,0};return x;}
int main(){
 auto x=fixture();auto p=prepare(x,{3,9});CHECK(p.status()==Status::ok);x.values[1]=999;CHECK(p.source().values[1]==18);CHECK(!p.calibration_state_declared());CHECK(p.source().event_ids[0]==p.source().event_ids[1]);CHECK(p.source().uncertainty_matrix[1]!=p.source().uncertainty_matrix[3]);
 CHECK(!p.calibration_provenance_known());
 auto s=p.select(Selection::pantheon_zhd_gt_001);CHECK(s.status==Status::ok);CHECK(s.source_indices==std::vector<std::size_t>{1});CHECK(s.mask==std::vector<std::uint8_t>({0,1,0}));
 x=fixture();std::swap(x.uncertainty_axis_ids[0],x.uncertainty_axis_ids[1]);CHECK(prepare(x,{3,9}).status()==Status::invalid_shape);
 x=fixture();x.measurement_ids[1]=x.measurement_ids[0];CHECK(prepare(x,{3,9}).status()==Status::invalid_identity);
 x=fixture();x.role=Role::posterior_summary;CHECK(prepare(x,{3,9}).status()==Status::incompatible_semantics);
 x=fixture();x.missing[1]=1;p=prepare(x,{3,9});CHECK(p.status()==Status::ok);CHECK(p.select(Selection::pantheon_zhd_gt_001).status==Status::missing_required_value);CHECK(p.source().missing[1]==1);
 x.source_selection={1,0,1};p=prepare(x,{3,9});CHECK(p.select(Selection::pantheon_zhd_gt_001).status==Status::ok);
 x=fixture();x.values[1]=std::numeric_limits<double>::quiet_NaN();p=prepare(x,{3,9});CHECK(p.select(Selection::all).status==Status::nonfinite_required_value);
 x=fixture();CHECK(prepare(x,{2,9}).status()==Status::resource_limit);CHECK(prepare(x,{3,8}).status()==Status::resource_limit);
 x=fixture();x.uncertainty=Uncertainty::precision;x.uncertainty_unit=UncertaintyUnit::inverse_magnitude_squared;p=prepare(x,{3,9});CHECK(p.status()==Status::ok);CHECK(p.source().uncertainty_matrix==x.uncertainty_matrix);
 x=fixture();x.ordering_provenance.clear();CHECK(prepare(x,{3,9}).status()==Status::incompatible_semantics);
 x=fixture();x.calibration=Calibration::not_applicable;CHECK(prepare(x,{3,9}).status()==Status::incompatible_semantics);
 x=fixture();x.profile=static_cast<Profile>(99);CHECK(prepare(x,{3,9}).status()==Status::incompatible_semantics);
 x=fixture();x.calibration=static_cast<Calibration>(99);CHECK(prepare(x,{3,9}).status()==Status::incompatible_semantics);
 x=fixture();x.uncertainty=static_cast<Uncertainty>(99);CHECK(prepare(x,{3,9}).status()==Status::incompatible_semantics);
 x=fixture();x.unit=Unit::metre;CHECK(prepare(x,{3,9}).status()==Status::incompatible_semantics);
 x=fixture();x.uncertainty_unit=UncertaintyUnit::inverse_magnitude_squared;CHECK(prepare(x,{3,9}).status()==Status::incompatible_semantics);
 x=fixture();x.component=static_cast<Component>(99);CHECK(prepare(x,{3,9}).status()==Status::incompatible_semantics);
 x=fixture();CHECK(prepare(x,{3,9,1}).status()==Status::resource_limit);
 x=fixture();p=prepare(x,{3,9});CHECK(p.select(static_cast<Selection>(99)).status==Status::incompatible_semantics);
 x=fixture();x.missing[2]=1;p=prepare(x,{3,9});auto failed=p.select(Selection::all);CHECK(failed.status==Status::missing_required_value);CHECK(failed.mask.empty()&&failed.source_indices.empty());CHECK(p.source().values.size()==3);
}

