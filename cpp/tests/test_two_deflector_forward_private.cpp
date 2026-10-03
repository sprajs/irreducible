// Shared-source regression for F02/F03; NOT independent numerical evidence.
// Include the actual private owner body without adding unused public outputs.
#include "../src/two_deflector_forward.cpp"
#include <iostream>
#include <stdexcept>
using namespace irred::lensing;
int main() {
  try {
    using W=long double;
    auto need=[](bool ok,const char *why) { if(!ok) throw std::runtime_error(why); };
    const double denormal=std::numeric_limits<double>::denorm_min();
    need(returned_integral_rounding(denormal)==W(denormal)/2&&
         returned_integral_rounding(denormal)>0,"subnormal returned-integral cast remains a positive wide diagnostic");
    need(returned_integral_rounding(1)==std::ldexp(1.L,-53),"binary64 binade cast uses larger adjacent half spacing");
    TwoDeflectorScene scene;
    scene.deflectors={SoftenedPotentialComponent{{0,0},.6,.3,.8,.125},
                      SoftenedPotentialComponent{{.4,-.15},.4,.2,.9,-.25}};
    scene.source={{0,0},2,2,0,1e13}; scene.theta_scale_radians=1e-5;
    scene.exposure_seconds=600; scene.origin="private source-regression only";
    const AffineDetectorCutout cutout{{0,0},{1,1,1,1.08}};
    const std::array<ShiftPSFComponent,1> psf{{{{0,0},1}}};
    const ForwardPixelRectangle pixel{49.99999,50.00001,-50.00001,-49.99999,0};
    const ForwardPixelPolicy policy;
    std::array<Matrix,3> derived{};
    for(std::size_t j=0;j<2;++j) {
      const auto &d=scene.deflectors[j];
      derived[j]=rotate_diagonal(d.angle_radians,1,1/(W(d.axis_ratio)*d.axis_ratio));
    }
    derived[2]=rotate_diagonal(0,.25,.25);
    ForwardWork batch,row;
    Context context{scene,cutout,derived,psf,pixel,policy,batch,row};
    const auto coordinate=context.detector(50,-50);
    need(coordinate.coordinate_error>64*eps*(1+norm(coordinate.value)),
         "affine cancellation retains original product magnitude");
    W map_error=0,unperturbed_error=0;
    const auto mapped=context.map(coordinate,map_error);
    const auto unperturbed=context.map({coordinate.value,0},unperturbed_error);
    need(mapped.x==unperturbed.x&&mapped.y==unperturbed.y&&map_error>unperturbed_error,
         "affine sensitivity charges diagnostics without changing the physical map");
    W sensitivity=1;
    for(const auto &d:scene.deflectors) sensitivity+=W(d.strength)/(W(d.axis_ratio)*d.axis_ratio*d.core);
    need(map_error>=unperturbed_error+sensitivity*coordinate.coordinate_error,
         "global map derivative bound propagates per-sample coordinate error");
    std::cout<<"PASS private source regressions F02/F03 (shared ancestry)\n";
    return 0;
  } catch(const std::exception &e) { std::cerr<<"FAIL private "<<e.what()<<'\n'; return 1; }
}
