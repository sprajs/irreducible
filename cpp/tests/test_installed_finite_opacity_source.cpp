// Fresh public-header/static-archive consumer. No private equation includes.
#include "irred/finite_opacity_source.hpp"
#include <iostream>
#include <type_traits>
#include <utility>
int main() {
  using namespace irred::cosmology;
  using S=irred::numerics::Status;
  static_assert(!std::is_copy_constructible_v<FiniteOpacitySourceProducer>);
  static_assert(!std::is_copy_constructible_v<FiniteOpacitySourceResult>);
  FiniteOpacityRequest invalid;
  auto refused=prepare_finite_opacity_source(std::move(invalid));
  if(refused.status()==S::ok || refused.identity())return 1;
  auto result=refused.produce();
  if(result.status==S::ok || result.source || result.source_numerically_admitted || result.boundary())return 2;
  auto moved=std::move(refused);
  if(refused.status()!=S::invalid_input || refused.identity() || moved.status()==S::ok)return 3;
  if(finite_opacity_model_id.empty() || finite_opacity_method_id.empty() || finite_opacity_arithmetic_id.empty())return 4;
  std::cout<<"PASS installed finite opacity owned source public contract\n";
}
