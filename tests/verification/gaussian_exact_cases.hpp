#pragma once
// Independently derived 2x2 adjugate/scalar Gaussian cases; no production ancestry.
namespace verifier_gaussian {
struct Rational {long long numerator,denominator; constexpr long double value()const{return static_cast<long double>(numerator)/denominator;}};
inline constexpr Rational covariance[4]={{4,1},{1,1},{1,1},{9,1}};
inline constexpr Rational precision[4]={{9,35},{-1,35},{-1,35},{4,35}};
inline constexpr Rational residual[2]={{2,1},{-3,1}};
inline constexpr Rational quadratic{12,5};
inline constexpr Rational offset_precision{11,35};
inline constexpr Rational flat_offset_fit{7,11};
inline constexpr Rational projected_chi2{25,11};
inline constexpr Rational marginal_A_variance{4,1};
inline constexpr Rational conditional_A_variance{35,9};
// Proper nuisance t~N(1/3,2), y-model residual r=(2,-3), design=(1,1).
// y marginal covariance [[6,3],[3,11]], residual around prior mean (5/3,-10/3).
inline constexpr Rational proper_marginal_determinant{57,1};
inline constexpr Rational proper_marginal_quadratic{1175,513};
inline constexpr Rational nuisance_posterior_mean{77,171};
inline constexpr Rational nuisance_posterior_variance{70,57};
// Normalized log marginal = -log(2pi)-log(57)/2-(1175/513)/2.
// With zero residual and C -> 4*C, determinant ->16*det(C), loglike shift=-log(4).
}
