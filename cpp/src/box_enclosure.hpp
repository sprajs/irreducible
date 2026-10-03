#pragma once
#include "irred/gaussian_box.hpp"
namespace irred::detail {
struct BoxRefusal { numerics::Status status; std::optional<long double> witness={}; };
struct BoxCdfCounters { std::size_t nodes=0,evaluations=0,bisections=0; };
bool box_supported() noexcept;
bool box_valid_policy(const statistics::BoxPolicy &) noexcept;
double box_scalar(long double);
statistics::BoxInterval box_sqrt(double);
long double box_normalize(statistics::GaussianBoxResult &,const statistics::BoxSupport &,
                   const statistics::BoxPolicy &,std::size_t observation_dimension);
void box_quantile(statistics::GaussianBoxResult &,const statistics::BoxSupport &,
                  statistics::BoxMarginalRequest,const statistics::BoxPolicy &,
                  BoxCdfCounters &,long double excluded_mass);
statistics::BoxInterval box_observation_constant(double full_logdet,double training_logdet);
statistics::BoxInterval box_density(statistics::BoxInterval joint,
    statistics::BoxInterval training,statistics::BoxInterval log_constant,double maximum_width);
} // namespace irred::detail
