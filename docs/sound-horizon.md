# Conditional sound horizon

`cosmology.sound_horizon` computes the comoving sound horizon in Mpc at an explicitly supplied drag redshift. It does not predict that redshift, recombination, an ionization history or a CMB observable. The calculation is conditional on pressureless matter, massless radiation and a cosmological constant in a flat background, with today's scale factor equal to one.

Supply H0 in km/s/Mpc and today's critical-density fractions Omega_m, Omega_r, Omega_b and Omega_gamma. Baryons are a subset of matter; photons are a subset of radiation. H0 must be finite and positive, radiation and photons strictly positive, and all other fractions nonnegative. The flat closure is Omega_Lambda = 1 − Omega_m − Omega_r >= 0. No temperature or effective-neutrino-number conversion is implicit. A finite nonnegative `z_drag` and its `drag_origin` are supplied with each point. Zero redshift is a mathematical control of this approximation, not a claim that actual late-time baryons remain tightly coupled to radiation.

The equations are

\[
H^2(a)/H_0^2=\Omega_r/a^4+\Omega_m/a^3+\Omega_\Lambda,
\quad B=3\Omega_b/(4\Omega_\gamma),\quad a_d=(1+z_d)^{-1},
\]

\[
r_s=\frac{c\,a_d}{\sqrt3 H_0}\int_0^1
\frac{dt}{\sqrt{\Omega_r+\Omega_m a_dt+\Omega_\Lambda(a_dt)^4}
\sqrt{1+B a_dt}}.
\]

The exact changes of coordinate a = 1/(1+z) and t = a/a_d cover the entire infinite-redshift interval. There is no artificial redshift cutoff or omitted-tail estimate. Positive radiation makes the a=0 endpoint finite; numerical quadrature error still remains.

Each batch owns its attempted points and origin strings. Rows retain individual numerical causes and actual callback counts, including failed work. Successful rows contain a positive normal binary64 `sound_horizon_mpc` and an `error_estimate_mpc`; failed rows have neither usable scalar. A configured batch quota can return an empty work-limit result. Malformed transport descriptors are rejected before copying their payloads.

The numerical policy supplies absolute tolerance in Mpc, relative tolerance, per-point and total callback limits, depth, point count and combined native/wire payload bytes. Empty native batches allocate no dynamic payload; an allocated successful boundary owner still counts towards the combined envelope. A fixed bounded empty work-limit diagnostic owner is permitted outside this admission budget, even when the requested quota is zero. Allocator bookkeeping, stack frames and RSS are outside this payload estimate. The supported arithmetic requires round-to-nearest and long double with at least 64 mantissa bits and exponent range 16384. An explicit arithmetic allowance rejects unattainable tolerances; the quadrature estimate and allowance are not rigorous bounds on libm or quadrature error.

Permanent native controls use radiation-only and Lambda=0 analytic limits, independent fixed Gauss-Legendre refinement, exact closure adversaries, subsets, scaling, monotonicity and failure/resource cases. The named comparison budget is 1e-9 Mpc + 2e-11 times the reference magnitude; reference refinement must occupy at most 5% of that budget. These controls do not automatically qualify arbitrary requests or a downstream BAO likelihood. Runtime interpretation remains unqualified.
