// Independent semantic adversaries, no scientific arithmetic/legacy package oracle.
#include "cosmology/observations.hpp"
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
using namespace cosmology::observations;
namespace {
int checks=0;void check(bool ok,const char* name){++checks;if(!ok)throw std::runtime_error(name);}
Input fixture(){Input x{};x.profile=Profile::gaussian_fixture_v1;x.role=Role::synthetic_control;x.unit=Unit::magnitude;x.calibration=Calibration::unknown;x.uncertainty=Uncertainty::covariance;x.uncertainty_unit=UncertaintyUnit::magnitude_squared;x.component=Component::unknown;x.ordering_provenance="synthetic explicit row axis";x.table_sha256=std::string(64,'a');x.uncertainty_sha256=std::string(64,'b');x.measurement_ids={"row1","row2","row3"};x.event_ids={"eventA","eventA","eventB"};x.uncertainty_axis_ids=x.measurement_ids;x.values={-2,0,std::numeric_limits<double>::quiet_NaN()};x.missing={0,0,1};x.quality={0,2,128};x.source_selection={1,1,0};x.uncertainty_matrix={4,1,2,1.00000003,9,3,2,3,-1};return x;}
}
int main(){try{
 auto x=fixture();const auto original=x.uncertainty_matrix;auto p=prepare(x,{3,9});check(p.status()==Status::ok,"preserve asymmetric indefinite source without repair");
 check(!p.calibration_state_declared()&&!p.calibration_provenance_known(),"unknown calibration stays unknown");
 check(p.source().calibration_provenance.empty()&&p.source().dependence_provenance.empty(),"unknown lineage not invented");
 check(p.source().uncertainty_matrix==original,"full covariance retained");check(p.source().uncertainty_axis_ids==x.measurement_ids,"explicit covariance order retained");check(p.source().event_ids[0]==p.source().event_ids[1]&&p.source().measurement_ids[0]!=p.source().measurement_ids[1],"repeated event is not deduplicated");
 auto selected=p.select(Selection::all);check(selected.status==Status::ok&&selected.source_indices==std::vector<std::size_t>({0,1})&&selected.mask==std::vector<std::uint8_t>({1,1,0}),"exact source map excluded missing retained");check(std::isnan(p.source().values[2])&&p.source().missing[2]==1,"source missing payload retained not selected");check(p.source().quality==std::vector<std::uint64_t>({0,2,128})&&p.source().quality_dictionary.empty(),"raw flags unknown dictionary retained");
 auto wrong=p.select(Selection::pantheon_zhd_gt_001);check(wrong.status==Status::incompatible_semantics&&wrong.mask.empty()&&wrong.source_indices.empty(),"all failed derived selection maps empty");
 x.source_selection={1,1,1};p=prepare(x,{3,9});selected=p.select(Selection::all);check(selected.status==Status::missing_required_value&&selected.mask.empty()&&selected.source_indices.empty(),"late missing row no partial map");check(p.source().source_selection==x.source_selection,"original source mask preserved after failure");
 x=fixture();x.uncertainty=Uncertainty::precision;x.uncertainty_unit=UncertaintyUnit::inverse_magnitude_squared;p=prepare(x,{3,9});check(p.status()==Status::ok&&p.source().uncertainty==Uncertainty::precision&&p.source().uncertainty_matrix==original,"precision tag not covariance reinterpretation");
 x=fixture();std::swap(x.uncertainty_axis_ids[0],x.uncertainty_axis_ids[1]);check(prepare(x,{3,9}).status()!=Status::ok,"permuted uncertainty axis rejects");
 x=fixture();x.ordering_provenance.clear();check(prepare(x,{3,9}).status()==Status::incompatible_semantics,"missing ordering declaration rejects");
 x=fixture();x.role=Role::posterior_summary;check(prepare(x,{3,9}).status()==Status::incompatible_semantics,"posterior cannot become observation");
 x=fixture();x.calibration=Calibration::released_corrected;p=prepare(x,{3,9});check(p.status()==Status::ok&&p.calibration_state_declared()&&!p.calibration_provenance_known(),"declared correction state not known recipe");
 x=fixture();x.measurement_ids[1]=x.measurement_ids[0];check(prepare(x,{3,9}).status()==Status::invalid_identity,"duplicate measurement identity rejects");
 x=fixture();x.uncertainty_matrix[0]=std::numeric_limits<double>::infinity();check(prepare(x,{3,9}).status()==Status::nonfinite_required_value,"nonfinite uncertainty always rejects");
 x=fixture();check(prepare(x,{2,9}).status()==Status::resource_limit,"row budget");check(prepare(x,{3,8}).status()==Status::resource_limit,"matrix budget");check(prepare(x,{3,9,64}).status()==Status::resource_limit,"metadata byte budget");
 x=fixture();x.profile=Profile::pantheon_plus_released_v1;x.role=Role::released_fitted_summary;x.zhd={std::nextafter(.01,0.),.01,std::nextafter(.01,1.)};x.zcmb={.2,.3,.4};x.zhel={.21,.31,.41};x.zhd_missing={0,0,0};x.zcmb_missing={0,0,0};x.zhel_missing={0,0,0};x.values={1,2,3};x.missing={0,0,0};x.source_selection={1,1,1};
 p=prepare(x,{3,9});check(p.status()==Status::ok,"explicit released summary profile");selected=p.select(Selection::pantheon_zhd_gt_001);check(selected.status==Status::ok&&selected.source_indices==std::vector<std::size_t>({2}),"strict observed zHD threshold neighbors");check(p.source().zcmb==x.zcmb&&p.source().zhel==x.zhel,"redshift frames remain distinct");
 x.calibration=Calibration::not_applicable;check(prepare(x,{3,9}).status()==Status::incompatible_semantics,"released calibration cannot be notapplicable");
 std::printf("{\"suite\":\"independent_observation_semantics\",\"checks\":%d,\"passed\":true,\"scope\":\"metadata, immutable full uncertainty and selection only; no evaluated likelihood\"}\n",checks);return 0;
 }catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 1;}}
