# Independent validation only; never a production shared-physics route.
# Fixed high-precision GL in log(1+z), independently fixed FD momentum GL.
import mpmath as mp, math, json, sys
from pathlib import Path
out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=True)
models=[('eds',70.,1.,0.,0.),('milne',70.,0.,0.,0.),
        ('open',67.4,.3,.0001,.5999),('closed',67.4,.3,.0001,.7999),
        ('closed-dust',70.,10.,0.,0.),('fd',)]
pairs=[(.1,.5),(.5,2.),(2.,10.),(1.,math.nextafter(1.,2.))]
results={}
for digits,nq,nz in ((60,160,64),(80,240,96)):
 mp.mp.dps=digits;pi=mp.pi;c=mp.mpf(299792458);G=mp.mpf('6.67430e-11')
 mpc=mp.mpf(648000000000)*mp.mpf(149597870700)/pi
 zn,zw=mp.gauss_quadrature(nz,'legendre');qn,qw=mp.gauss_quadrature(nq,'legendre')
 q=[40*(v+1) for v in qn];weights=[40*qw[i]*q[i]**2/(mp.exp(q[i])+1) for i in range(nq)]
 rows=[]
 for model in models:
  name=model[0];tail=mp.mpf(0)
  if name=='fd':
   H=mp.mpf(67.4);h2=(H/100)**2;T=mp.mpf(2.7255);ratio=mp.mpf(.71611)
   ev=mp.mpf('1.602176634e-19');kb=mp.mpf('1.380649e-23')/ev
   hbar=mp.mpf('6.62607015e-34')/(2*pi);tnu=T*ratio*kb
   critical=3*(H*1000/mpc)**2*c*c/(8*pi*G)/ev*(hbar*c/ev)**3
   coef=2*tnu**4/(2*pi*pi*critical);M=mp.mpf(.06)/tnu
   og=pi*pi/15*(kb*T)**4/critical;r0=(mp.mpf(4)/11)**(mp.mpf(1)/3)
   radiation=og*(1+mp.mpf(7)/8*r0**4*(mp.mpf(3.046)-(ratio/r0)**4))
   matter=(mp.mpf(.02237)+mp.mpf(.12))/h2;k=mp.mpf(0)
   def rho(a):return mp.fsum(weights[i]*mp.sqrt(q[i]**2+(M*a)**2) for i in range(nq))
   lam=1-radiation-matter-coef*rho(1)
   def inverse(z):
    a=1/(1+z);return a*a/mp.sqrt(radiation+matter*a+lam*a**4+coef*rho(a))
   exp_tail=lambda n:mp.exp(-80)*mp.factorial(n)*mp.fsum(mp.mpf(80)**j/mp.factorial(j) for j in range(n+1))
   delta=2*coef*(exp_tail(3)+M*exp_tail(2));D0=radiation+coef*rho(0)
   tail=delta/(mp.sqrt(D0-delta)*(mp.sqrt(D0)+mp.sqrt(D0-delta)))
  else:
   _,H0,m0,r0,l0=model;H,m,r,lam=map(mp.mpf,(H0,m0,r0,l0))
   # Actual retained binary64 closure, not a silently re-rounded alternative.
   k=mp.mpf(1.-m0-r0-l0)
   def inverse(z):
    u=1+z;return 1/mp.sqrt(r*u**4+m*u**3+k*u*u+lam)
  unit=c/(H*1000)
  def radial(lo,hi):
   lower=mp.log1p(lo);width=mp.log1p((hi-lo)/(1+lo))
   return unit*width/2*mp.fsum(zw[i]*mp.exp(lower+width*(zn[i]+1)/2)*
       inverse(mp.expm1(lower+width*(zn[i]+1)/2)) for i in range(nz))
  def angular(radial,z):
   x=radial/unit
   if k>0:value=mp.sinh(mp.sqrt(k)*x)/mp.sqrt(k);response=mp.cosh(mp.sqrt(k)*(abs(x)+abs(x)*tail))
   elif k<0:value=mp.sin(mp.sqrt(-k)*x)/mp.sqrt(-k);response=1
   else:value=x;response=1
   return unit*value/(1+z),response*abs(radial)*tail/(1+z)
  def quotient(a,b,d,scale):
   v=scale*a[0]*b[0]/d[0]
   lo=scale*(a[0]-a[1])*(b[0]-b[1])/(d[0]+d[1])
   hi=scale*(a[0]+a[1])*(b[0]+b[1])/(d[0]-d[1])
   return v,max(v-lo,hi-v)
  selected=pairs+([(.5,100.)] if name=='closed-dust' else [])
  for zl,zs in selected:
   l,s=mp.mpf(zl),mp.mpf(zs)
   dl=angular(radial(0,l),l);ds=angular(radial(0,s),s);dls=angular(radial(l,s),s)
   one=(mp.mpf(1),mp.mpf(0));beta=quotient(dls,one,ds,1)
   prod=quotient(dl,dls,one,1);sigma=quotient(ds,one,prod,c*c/(4*pi*G*mpc))
   delay=quotient(dl,ds,dls,1+l);values=[dl,ds,dls,beta,sigma,delay]
   rows.append({'model':name,'lens_redshift':zl,'source_redshift':zs,
       'values':[mp.nstr(v[0],digits) for v in values],
       'momentum_tail_errors':[mp.nstr(v[1],digits) for v in values]})
  print('completed',digits,name,flush=True)
 results[str(digits)]={'momentum_order':nq,'radial_order':nz,'rows':rows}
(out/'geometry-reference.json').write_text(json.dumps(results,indent=2)+'\n')
worst=mp.mpf(0);witness=None
for lo,hi in zip(results['60']['rows'],results['80']['rows']):
 for i,(low,high,tail) in enumerate(zip(lo['values'],hi['values'],hi['momentum_tail_errors'])):
  v=mp.mpf(high);e=abs(v-mp.mpf(low))+mp.mpf(tail)
  absolute=mp.mpf('1e-12' if i==3 else '1e-10' if i==4 else '1e-8')
  fraction=e/(absolute+mp.mpf('1e-8')*abs(v))
  if fraction>worst:worst=fraction;witness={'model':hi['model'],'zl':hi['lens_redshift'],'zs':hi['source_redshift'],'output':i}
receipt={'max_fraction':str(worst),'required_max_fraction':.05,'passed':bool(worst<=.05),'witness':witness}
(out/'geometry-reference-refinement.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(receipt,flush=True);assert worst<=.05
