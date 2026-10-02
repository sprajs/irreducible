#pragma once
// Original direct-SI reference for the frozen pure-H F=1/B2 closure.
// No engine headers, production helper/readback, external solver or rate table.
// Physical constants and long-double/system-libm ancestry are shared, declared
// inputs; Radau IIA, finite-difference Newton and GL integration are independent.
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <vector>
namespace hydrogen_thermal_reference {
using W = long double;
using State = std::array<W, 2>;
constexpr W kb=1.380649e-23L, h=6.62607015e-34L, c=299792458,
 ev=1.602176634e-19L, me=9.1093837139e-31L, mp=1.67262192595e-27L,
 G=6.67430e-11L, sigma=6.6524587051e-29L, lya=121.5682e-9L,
 binding=13.598434599702L*ev, excitation=h*c/lya;
struct Species { W temperature, weight; };
struct Model {
 W H0=67.4, baryon=.02237, cdm=.12, T0=2.7255, other=1.7e-5;
 std::vector<Species> species;
};
struct Stats {
 size_t steps=0, rhs=0, newton=0;
 W residual=0, min_x=1, min_neutral=1, min_theta=1,
   initial_stiffness=0, maximum_initial_step_stiffness=0;
};
W mpc() { return 1e6L*648000/std::numbers::pi_v<W>*149597870700.L; }
W radiation_constant() {
 return 8*std::pow(std::numbers::pi_v<W>,5)*std::pow(kb,4)/(15*h*h*h*c*c*c);
}
struct Physics {
 Model model;
 W initial, late, n0, photon, radiation, matter, lambda, H100;
 Physics(Model m, W zi, W zl):model(std::move(m)),initial(zi),late(zl) {
  H100=100000/mpc();
  W critical=3*H100*H100/(8*std::numbers::pi_v<W>*G);
  n0=critical*model.baryon/(mp+me-binding/(c*c));
  photon=radiation_constant()*std::pow(model.T0,4)/(c*c*critical);
  radiation=photon+model.other;
  for(auto a:model.species)
   radiation+=7.L/16*a.weight*radiation_constant()*std::pow(a.temperature,4)/(c*c*critical);
  matter=model.baryon+model.cdm;
  lambda=std::pow(model.H0/100,2)-radiation-matter;
  if(!(lambda>=0)) throw std::runtime_error("reference flat closure");
 }
 W hubble(W z) const {
  W u=1+z; return H100*std::sqrt(radiation*u*u*u*u+matter*u*u*u+lambda);
 }
 W density(W z) const { return n0*std::pow(1+z,3); }
 W gamma_ratio(W z,W x) const {
  return 8*sigma*radiation_constant()*std::pow(model.T0*(1+z),4)/(3*me*c*hubble(z))*x/(1+x);
 }
 W initial_x() const {
  W T=model.T0*(1+initial);
  W logS=1.5L*std::log(2*std::numbers::pi_v<W>*me*kb*T/(h*h))
    -binding/(kb*T)-std::log(density(initial));
  // Original monotone logit equation; no production quadratic branch.
  W a=-1000,b=1000;
  for(unsigned i=0;i<256;++i) {
   W w=(a+b)/2, soft=w>0 ? w+std::log1p(std::exp(-w)) : std::log1p(std::exp(w));
   if(2*w-soft>logS) b=w; else a=w;
  }
  W w=(a+b)/2; return w>=0 ? 1/(1+std::exp(-w)) : std::exp(w)/(1+std::exp(w));
 }
 State operator()(W s,const State& y) const {
  W u=(1+initial)*std::exp(-s),z=u-1,T=y[1]*model.T0*u,n=density(z),H=hubble(z);
  W t=T/10000, alpha=1e-19L*4.309L*std::pow(t,-.6166L)/(1+.6703L*std::pow(t,.53L));
  W beta=alpha*std::pow(2*std::numbers::pi_v<W>*me*kb*T/(h*h),1.5L)*std::exp(-(binding-excitation)/(kb*T));
  W K=lya*lya*lya/(8*std::numbers::pi_v<W>*H),neutral=1-y[0];
  W inhibition=(1+K*8.22458L*n*neutral)/(1+K*(8.22458L+beta)*n*neutral);
  return {-inhibition*(n*alpha*y[0]*y[0]-beta*neutral*std::exp(-excitation/(kb*T)))/H,
          -y[1]-gamma_ratio(z,y[0])*(y[1]-1)};
 }
 std::array<W,2> opacity_s(W s,W x) const {
  W u=(1+initial)*std::exp(-s),z=u-1,q=c*sigma*density(z)*x/hubble(z);
  return {q,q/(3*model.baryon/(4*photon*u))};
 }
};
bool positive(const State& y) { return y[0]>0 && y[0]<1 && y[1]>0 && std::isfinite(y[1]); }
W norm(const std::array<W,4>& a) {
 W n=0;for(W v:a) {
  if(!std::isfinite(v))return std::numeric_limits<W>::infinity();
  n=std::max(n,std::abs(v));
 }
 return n;
}
std::array<W,4> solve4(std::array<std::array<W,4>,4> a,std::array<W,4> b) {
 for(unsigned k=0;k<4;++k) {
  unsigned p=k;for(unsigned i=k+1;i<4;++i)if(std::abs(a[i][k])>std::abs(a[p][k]))p=i;
  if(!(std::abs(a[p][k])>1e-28L))throw std::runtime_error("reference Newton pivot");
  std::swap(a[p],a[k]);std::swap(b[p],b[k]);
  for(unsigned i=k+1;i<4;++i) {
   W r=a[i][k]/a[k][k];for(unsigned j=k;j<4;++j)a[i][j]-=r*a[k][j];b[i]-=r*b[k];
  }
 }
 std::array<W,4> x{};
 for(int i=3;i>=0;--i){W v=b[i];for(unsigned j=i+1;j<4;++j)v-=a[i][j]*x[j];x[i]=v/a[i][i];}
 return x;
}
struct Step { State end, first, last; };
template<class F> Step radau(F&& f,W s,W ds,State y,Stats& stats) {
 constexpr W A[2][2]{{5.L/12,-1.L/12},{3.L/4,1.L/4}};
 std::array<W,4> v{y[0],y[1],y[0],y[1]};
 auto residual=[&](const auto& t){
  State k[]{f(s+ds/3,{t[0],t[1]}),f(s+ds,{t[2],t[3]})};stats.rhs+=2;
  std::array<W,4> r{};
  for(unsigned i=0;i<2;++i)for(unsigned j=0;j<2;++j)
   r[2*i+j]=t[2*i+j]-y[j]-ds*(A[i][0]*k[0][j]+A[i][1]*k[1][j]);
  return r;
 };
 constexpr W tolerance=8e-17L;
 for(unsigned iteration=0;iteration<48;++iteration) {
  auto r=residual(v); W n=norm(r);++stats.newton;
  if(n<=tolerance) {
   stats.residual=std::max(stats.residual,n);++stats.steps;
   for(unsigned i=0;i<2;++i){stats.min_x=std::min(stats.min_x,v[2*i]);stats.min_neutral=std::min(stats.min_neutral,1-v[2*i]);stats.min_theta=std::min(stats.min_theta,v[2*i+1]);}
   return {{v[2],v[3]},{v[0],v[1]},{v[2],v[3]}};
  }
  if(!std::isfinite(n))throw std::runtime_error("reference nonfinite residual");
  std::array<std::array<W,4>,4> jac{};
  for(unsigned j=0;j<4;++j) {
   W delta=std::sqrt(std::numeric_limits<W>::epsilon())*(1+std::abs(v[j]));
   auto a=v,b=v;a[j]+=delta;b[j]-=delta;auto ra=residual(a),rb=residual(b);
   for(unsigned i=0;i<4;++i)jac[i][j]=(ra[i]-rb[i])/(2*delta);
  }
  for(auto& q:r)q=-q;
  auto direction=solve4(jac,r);W damping=1;bool accepted=false;
  for(unsigned i=0;i<64;++i) {
   auto trial=v;for(unsigned j=0;j<4;++j)trial[j]+=damping*direction[j];
   if(positive({trial[0],trial[1]}) && positive({trial[2],trial[3]}) && norm(residual(trial))<n){v=trial;accepted=true;break;}
   damping/=2;
  }
  if(!accepted)throw std::runtime_error("reference Newton line search");
 }
 throw std::runtime_error("reference Newton iteration ceiling");
}
struct Node { W s,z;State state,derivative;W thomson=0,drag=0; };
struct Value { W z,x,temperature,thomson,drag,opacity,visibility,survival; };
struct Result {
 Physics physics;Stats stats;std::vector<Node> nodes;std::vector<size_t> queries;
 W root=0;
 Value query(size_t i) const {
  auto n=nodes.at(queries.at(i));W tau=nodes.back().thomson-n.thomson,drag=nodes.back().drag-n.drag;
  auto q=physics.opacity_s(n.s,n.state[0]);W survival=std::exp(-tau),opacity=q[0]/(1+n.z);
  return {n.z,n.state[0],n.state[1]*physics.model.T0*(1+n.z),tau,drag,opacity,opacity*survival,survival};
 }
};
Result integrate(Model m,W initial,W late,const std::vector<W>& requested,unsigned refinement,W maximum_s_step=4e-4L) {
 Result out{Physics(m,initial,late),{}, {},{},0};
 std::vector<W> stops=requested;stops.push_back(initial);stops.push_back(late);stops.push_back(initial-.1L);
 std::sort(stops.begin(),stops.end(),std::greater<W>());stops.erase(std::unique(stops.begin(),stops.end()),stops.end());
 if(stops.front()!=initial||stops.back()!=late||!refinement)throw std::runtime_error("reference query domain");
 State state{out.physics.initial_x(),1};W s=0,at=0,ad=0;
 auto append=[&](W z){auto d=out.physics(s,state);++out.stats.rhs;out.nodes.push_back({s,z,state,d,at,ad});};append(initial);
 out.stats.initial_stiffness=1+out.physics.gamma_ratio(initial,1.L); // conservative x<=1 bound
 for(size_t j=1;j<stops.size();++j) {
  W end=std::log((1+initial)/(1+stops[j])),length=end-s;
  size_t n=std::max<size_t>(1,std::ceil(length/maximum_s_step));
  // The whole cell must be in the first .1 redshift. Testing its upper
  // endpoint also selected the following large cell and over-resolved it.
  if(stops[j]>=initial-.1L)n=std::max<size_t>(n,std::ceil(length*out.stats.initial_stiffness*8));
  n*=refinement;W begin=s,ds=length/n;
  for(size_t k=0;k<n;++k) {
   if(stops[j]>=initial-.1L)out.stats.maximum_initial_step_stiffness=std::max(out.stats.maximum_initial_step_stiffness,ds*out.stats.initial_stiffness);
   auto step=radau(out.physics,s,ds,state,out.stats);
   auto a=out.physics.opacity_s(s+ds/3,step.first[0]),b=out.physics.opacity_s(s+ds,step.last[0]);
   at+=ds*(3*a[0]+b[0])/4;ad+=ds*(3*a[1]+b[1])/4;
   state=step.end;s=begin+(k+1)*ds;append(k+1==n ? stops[j] : (1+initial)*std::exp(-s)-1);
  }
 }
 for(W z:requested){auto i=std::find_if(out.nodes.begin(),out.nodes.end(),[z](const Node& n){return n.z==z;});if(i==out.nodes.end())throw std::runtime_error("reference exact stop absent");out.queries.push_back(i-out.nodes.begin());}
 for(size_t i=1;i<out.nodes.size();++i) {
  W a=out.nodes.back().drag-out.nodes[i-1].drag,b=out.nodes.back().drag-out.nodes[i].drag;
  if(a>=1 && b<=1){out.root=out.nodes[i].z+(1-b)*(out.nodes[i-1].z-out.nodes[i].z)/(a-b);break;}
 }
 return out;
}
struct Gauss8 {
 std::array<W,8> x{},w{};
 Gauss8(){for(unsigned i=0;i<8;++i){W z=std::cos(std::numbers::pi_v<W>*(i+.75L)/8.5L),derivative=0;
  for(unsigned j=0;j<32;++j){W a=1,b=z;for(unsigned k=2;k<=8;++k){W next=((2*k-1)*z*b-(k-1)*a)/k;a=b;b=next;}derivative=8*(z*b-a)/(z*z-1);W d=b/derivative;z-=d;if(std::abs(d)<4*std::numeric_limits<W>::epsilon())break;}
  x[i]=z;w[i]=2/((1-z*z)*derivative*derivative);}}
};
template<class F> W gl8(F&& f,W a,W b){static const Gauss8 q;W sum=0;for(unsigned i=0;i<8;++i)sum+=q.w[i]*f((a+b)/2+(b-a)/2*q.x[i]);return (b-a)/2*sum;}
W dense_x(const Node& a,const Node& b,W s){W h=b.s-a.s,t=(s-a.s)/h;
 W x=(2*t*t*t-3*t*t+1)*a.state[0]+(t*t*t-2*t*t+t)*h*a.derivative[0]+(-2*t*t*t+3*t*t)*b.state[0]+(t*t*t-t*t)*h*b.derivative[0];
 if(!(x>0&&x<1))throw std::runtime_error("reference dense positivity");
 return x;}
struct Quadrature { W thomson=0,drag=0,mass=0,survival=0; };
Quadrature opacity_gl(const Result& r){
 std::vector<W> cumulative(r.nodes.size());W drag=0;
 for(size_t i=1;i<r.nodes.size();++i){auto a=r.nodes[i-1],b=r.nodes[i];
  auto q=[&](W s){return r.physics.opacity_s(s,dense_x(a,b,s));};
  cumulative[i]=cumulative[i-1]+gl8([&](W s){return q(s)[0];},a.s,b.s);
  drag+=gl8([&](W s){return q(s)[1];},a.s,b.s);}
 W total=cumulative.back(),mass=0;
 for(size_t i=1;i<r.nodes.size();++i){auto a=r.nodes[i-1],b=r.nodes[i];auto q=[&](W s){return r.physics.opacity_s(s,dense_x(a,b,s))[0];};
  mass+=gl8([&](W s){W part=gl8(q,a.s,s);return q(s)*std::exp(-(total-cumulative[i-1]-part));},a.s,b.s);}
 return {total,drag,mass,std::exp(-total)};
}
} // namespace hydrogen_thermal_reference
