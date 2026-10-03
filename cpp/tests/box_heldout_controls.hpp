#pragma once
#include "irred/gaussian_box_heldout.hpp"
#include <algorithm>
#include <cmath>
#include <utility>
#include <stdexcept>
namespace box_heldout_controls {
using namespace irred::statistics;
struct Input {
  Metadata rows;
  DesignMetadata design;
  BoxSupport box;
  BoxHeldoutMetadata heldout;
  std::size_t heldout_index=2;
  std::vector<double> covariance={1,0,.5,0,1,.25,.5,.25,1},x={1,0,0,1,.5,.25},offsets={1,2,3};
  Input() {
    rows.ordered_ids={"r0","r1","r2"};rows.measure="product d(u)";
    rows.table_identity="original dyadic correlated control/v1";
    rows.uncertainty_identity="full supplied dyadic3 covariance/v1";
    rows.ordering_provenance="literal r0,r1,r2";rows.calibration_provenance="synthetic fixed offsets";
    rows.dependence_provenance="full C including k=(1/2,1/4)";rows.source_semantics="synthetic controls";
    design={{"original0","original2"},{"u0","u2"},{},"u","dyadic2 active source coordinates/v1",rows.dependence_provenance};
    box={{"original0","original2"},{-8,-9},{8,9},"product d(original0)d(original2)","normalized original finite box/v1",
        "synthetic original1 literal0 point mass outside active2"};
    heldout.source_contract_id="synthetic-original-dyadic/v1";
    heldout.candidate_identity="withhold literal original r2/v1";
    heldout.conditioning_identity="fixed C/X/offsets/box; synthetic conditional law";
    heldout.original_row_lineage="r0,r1,r2 exact original order";
    heldout.event_lineage="synthetic controls; no observational event claim";
    heldout.calibration_dependence_identity=rows.dependence_provenance;
    heldout.heldout_unit="u";heldout.heldout_covariance_unit="u^2";heldout.heldout_measure="d(u)";
    heldout.ordered_original_parameter_ids={"original0","original1","original2"};
    heldout.active_original_parameter_indices={0,2};heldout.fixed_original_parameter_index=1;
    heldout.fixed_original_parameter_value=0;heldout.fixed_coordinate_provenance=box.fixed_coordinate_provenance;
  }
  Gaussian gaussian(MatrixKind kind=MatrixKind::covariance)const {
    return prepare_gaussian(covariance,kind,rows,covariance.size(),1e-10,irred::numerics::Arithmetic::longdouble_cpu_v1);
  }
  BoxHeldoutPolicy policy(const Gaussian&g,std::size_t trainings=1,std::size_t candidates=1,std::size_t requests=1)const {
    BoxHeldoutPolicy p;p.maximum_training_vectors=trainings;p.maximum_candidate_values=candidates;p.maximum_requests=requests;
    const auto setup=GaussianBoxHeldout::preparation_work_bound(g,design,box,heldout);
    std::size_t identity_bytes=0;for(std::size_t j=0;j<rows.ordered_ids.size();++j)if(j!=heldout_index)identity_bytes+=rows.ordered_ids[j].size();
    const auto evaluation=GaussianBoxHeldout::evaluation_work_bound(rows.ordered_ids.size(),design.ordered_parameter_ids.size(),trainings,candidates,requests,identity_bytes);
    if(!setup||!evaluation)throw std::runtime_error("source work bound unavailable");
    p.maximum_preparation_work_units=*setup;p.maximum_evaluation_work_units=*evaluation;return p;
  }
  GaussianBoxHeldout prepare(Gaussian&&g,BoxHeldoutPolicy p)const {
    return GaussianBoxHeldout::prepare(std::move(g),x,offsets,rows.ordered_ids,design,box,heldout_index,heldout,p);
  }
};
} // namespace box_heldout_controls
