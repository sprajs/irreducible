#pragma once
#include <cmath>
namespace irred::cosmology::detail {
// Sole retained-curvature transverse law. Its callers own source/path/error
// admission; chi is dimensionless H0*radial/c and k is the actual retained Ok.
inline long double curved_transverse(long double chi,long double k) {
  const long double t=k*chi*chi;
  if(std::abs(t)<1e-4L)
    return chi*(1+t*(1.L/6+t*(1.L/120+t*(1.L/5040+t/362880))));
  const long double q=std::sqrt(std::abs(k));
  return k>0?std::sinh(q*chi)/q:std::sin(q*chi)/q;
}
inline long double curved_transverse_response(long double chi,
                                             long double chi_error,
                                             long double k) {
  return k>0?std::cosh(std::sqrt(k)*(std::abs(chi)+chi_error)):1;
}
} // namespace irred::cosmology::detail
