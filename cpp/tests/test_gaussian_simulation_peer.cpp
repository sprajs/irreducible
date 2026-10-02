// Independent dyadic triangular facts and 90-digit Decimal transcendental facts.
#include "irred/gaussian_simulation.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>
namespace {
namespace n=irred::numerics;using W=long double;unsigned checks=0;
void need(bool b,const char*w){++checks;if(!b)throw std::runtime_error(w);}
void near(W a,W b,W budget){need(std::abs(a-b)<=budget*(1+std::abs(b)),"independent frozen reference");}
void run(){// L=[2,0,0;.5,1,0;-.25,.125,.5], exact supplied C=L L^T.
 const std::array<double,9>C{4,1,-.5,1,1.25,0,-.5,0,.328125};
 for(auto arithmetic:{n::Arithmetic::binary64_legacy_v1,n::Arithmetic::longdouble_cpu_v1}){
 auto f=n::cholesky(C,3,100,arithmetic);need(f.status()==n::Status::ok,"dyadic covariance admitted");
 for(std::array<W,3>z: {std::array<W,3>{1,-2,3},std::array<W,3>{0,0,0},std::array<W,3>{-4,.5,-1}}){
 auto v=n::colour(f,z,100,1024,1e-10);need(v.status==n::Status::ok,"retained colour");
 const std::array<W,3> expected{2*z[0],.5L*z[0]+z[1],-.25L*z[0]+.125L*z[1]+.5L*z[2]};
 for(size_t i=0;i<3;++i){near(v.value[i],expected[i],2e-18L);need(v.absolute_error_estimates[i]>=std::abs(v.value[i]-expected[i]),"arithmetic estimate covers exact fixture");}}
 const std::array<W,3>z{1,2,3};auto bound=n::colouring_payload_bound(3);need(bound.has_value(),"colour bound");
 need(n::colour(f,z,100,*bound-1,1e-10).status==n::Status::work_limit,"colour payload boundary");
 }
 const std::array<std::array<std::uint32_t,4>,4>words{{{0,0,0,0},{0,UINT32_MAX,UINT32_MAX,0},{0,2147483647,1073741823,0},{0,123456789,987654321,0}}};
 const std::array<W,4> expected{6.76370563500189520164273164794584827896L,.0000152587890629440892057960934192990747549L,8.61227204468765399962667142034704743609e-10L,.334655772280110387830408546412615361708L};
 for(size_t i=0;i<4;++i)near(irred::random::normal_cosine(words[i]),expected[i],2e-18L);
}
}
int main(){try{run();std::cout<<"PASS "<<checks<<" independent Gaussian simulation controls\n";}catch(const std::exception&e){std::cerr<<"FAIL "<<checks<<": "<<e.what()<<'\n';return 1;}}
