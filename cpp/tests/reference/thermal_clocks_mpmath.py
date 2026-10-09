# Independent high-precision validation only; never production shared physics.
# Distinct fixed Gauss-Legendre full-support strategy, no production tail brackets.
import mpmath as mp,json,sys
from pathlib import Path
out=Path(sys.argv[1] if len(sys.argv)>1 else ".");out.mkdir(parents=True,exist_ok=True)
results={}
for digits,nq,na in ((80,160,256),(100,240,384)):
 mp.mp.dps=digits
 qnodes,qweights=mp.gauss_quadrature(nq,'legendre')
 anodes,aweights=mp.gauss_quadrature(na,'legendre')
 q=[40*(v+1) for v in qnodes]
 w=[40*qweights[i]*q[i]**2/(mp.exp(q[i])+1) for i in range(nq)]
 pi=mp.pi;c=mp.mpf(299792458);ev=mp.mpf('1.602176634e-19');kb=mp.mpf('1.380649e-23')/ev
 hbar=mp.mpf('6.62607015e-34')/(2*pi);mpc=mp.mpf(648000000000)*mp.mpf(149597870700)/pi
 H=mp.mpf('67.4');h2=(H/100)**2;T=mp.mpf('2.7255');ratio=mp.mpf('.71611');tnu=T*ratio*kb
 critical=3*(H*1000/mpc)**2*c*c/(8*pi*mp.mpf('6.67430e-11'))/ev*(hbar*c/ev)**3
 coefficient=2*tnu**4/(2*pi*pi*critical);M=mp.mpf('.06')/tnu
 radiation=pi*pi/15*(kb*T)**4/critical*(1+mp.mpf(7)/8*(mp.mpf(4)/11)**(mp.mpf(4)/3)*(mp.mpf('3.046')-(ratio/(mp.mpf(4)/11)**(mp.mpf(1)/3))**4))
 matter=(mp.mpf('.02237')+mp.mpf('.12'))/h2
 def rho(a):
  y=M*a
  return mp.fsum(w[i]*mp.sqrt(q[i]*q[i]+y*y) for i in range(nq))
 lam=1-radiation-matter-coefficient*rho(1)
 def D(a):return radiation+matter*a+lam*a**4+coefficient*rho(a)
 time_scale=mpc/(H*1000)/(mp.mpf('365.25')*86400*10**9);distance_scale=c/(H*1000)
 def integrate(a,age):
  # Full support a'=a*x^2: resolves dust/radiation endpoint independently.
  vals=[]
  for i,v in enumerate(anodes):
   x=(v+1)/2;ap=a*x*x
   vals.append(aweights[i]/2*2*a*x*(ap if age else 1)/mp.sqrt(D(ap)))
  return mp.fsum(vals)
 def lookback(a):
  # Full finite interval; direct support preserves tiny near-one clock.
  return mp.fsum(aweights[i]/2*(1-a)*(a+(1-a)*(anodes[i]+1)/2)/mp.sqrt(D(a+(1-a)*(anodes[i]+1)/2)) for i in range(na))
 values=[]
 # Decimal string of actual binary64 supplied coordinates for near-one reference.
 for astr in ('0.001','0.01','0.1','0.5','1.0','0.99999999999999988897769753748434595763683319091796875'):
  a=mp.mpf(astr);eta=distance_scale*integrate(a,False)
  values.append({'a':astr,'age_gyr':mp.nstr(time_scale*integrate(a,True),digits),
      'lookback_gyr':mp.nstr(time_scale*lookback(a),digits),
      'comoving_particle_horizon_mpc':mp.nstr(eta,digits),
      'proper_particle_horizon_mpc':mp.nstr(a*eta,digits)})
 # Positive omitted momentum tail, uniform y<=M; Lambda normalization retained.
 exp_tail=lambda n:mp.exp(-80)*mp.factorial(n)*mp.fsum(mp.mpf(80)**k/mp.factorial(k) for k in range(n+1))
 deltaD=2*coefficient*(exp_tail(3)+M*exp_tail(2))
 D0=D(0);radius=deltaD/(mp.sqrt(D0-deltaD)*(mp.sqrt(D0)+mp.sqrt(D0-deltaD)))
 results[str(digits)]={'momentum_order':nq,'outer_order':na,'values':values,
    'uniform_inverse_sqrt_momentum_tail_relative_bound':mp.nstr(radius,digits),
    'lambda':mp.nstr(lam,digits)}
 print('completed precision',digits,flush=True)
open(out/'clock-reference.json','w').write(json.dumps(results,indent=2)+'\n')
# Allocation refinement gate separate from comparison against native output.
worst=mp.mpf(0)
for low,high in zip(results['80']['values'],results['100']['values']):
 for key in ('age_gyr','lookback_gyr','comoving_particle_horizon_mpc','proper_particle_horizon_mpc'):
  v=mp.mpf(high[key]);e=abs(v-mp.mpf(low[key]))+abs(v)*mp.mpf(results['100']['uniform_inverse_sqrt_momentum_tail_relative_bound'])
  fraction=e/(mp.mpf('1e-8')+mp.mpf('2e-10')*abs(v));worst=max(worst,fraction)
open(out/'clock-reference-refinement.json','w').write(json.dumps({'max_fraction':str(worst),'required_max_fraction':.05,'passed':bool(worst<=.05)},indent=2)+'\n')
print('maximum reference refinement fraction',mp.nstr(worst,10),flush=True)
assert worst<=.05
