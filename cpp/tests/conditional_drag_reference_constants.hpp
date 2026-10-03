#pragma once
// Auxiliary enclosures of exact declared source numbers. The unchanged nominal
// reference still owns its stored constants and expression order. This adapter
// does not turn central atomic values into a newly measured input law.
#include "conditional_drag_reference_certificate.hpp"

namespace conditional_drag_reference_certificate {
struct AtomicNumbers {
  Interval kb, h, c, ev, me, G, sigma, chiH, chiHe;
  // Automatic/runtime construction only; prevent () from zeroing Wide slots.
  AtomicNumbers() noexcept {}
};
struct ColdImage {
  Interval value;
  bool absolute_cold = false;
  ColdImage() noexcept {} // No floating initialization under () or {}.
};
class SourceNumbers {
  Owner &owner_;
  Arithmetic &arithmetic_;
public:
  explicit SourceNumbers(Arithmetic &a) noexcept : owner_(a.owner()), arithmetic_(a) {}
  SourceNumbers(const SourceNumbers &) = delete;
  SourceNumbers &operator=(const SourceNumbers &) = delete;
  SourceNumbers(SourceNumbers &&) = delete;
  SourceNumbers &operator=(SourceNumbers &&) = delete;

  bool integer(Interval &out, std::int64_t value) noexcept {
    if (!owner_.guard()) return false;
    out.valid = false;
    if (!owner_.arithmetic_ready() || !owner_.action(Action::floating(1, 1)))
      return false;
    W stored; stored = static_cast<W>(value);
    // Every signed64 integer is exactly representable under the admitted
    // digits64 binary80 profile; zero is an exact integer witness.
    return arithmetic_.point(out, stored);
  }
  bool rational(Interval &out, std::int64_t numerator,
                std::int64_t denominator) noexcept {
    if (!owner_.guard()) return false;
    out.valid = false;
    if (denominator == 0) return owner_.fail(Cause::invalid_input);
    Interval n, d;
    return integer(n, numerator) && integer(d, denominator) &&
           arithmetic_.divide(out, n, d);
  }
  bool power(Interval &out, const Interval &base, unsigned exponent) noexcept {
    if (!owner_.guard()) return false;
    if (&out == &base) return owner_.fail(Cause::invalid_ownership);
    out.valid = false;
    if (exponent > 128) return owner_.fail(Cause::invalid_input);
    Interval checked_base, buffers[2];
    if (!arithmetic_.copy(checked_base, base) || !integer(buffers[0], 1)) return false;
    unsigned current = 0, step = 0;
    for (;;) {
      if (!owner_.guard()) return false;
      if (step == exponent) break;
      if (!owner_.iteration() ||
          !arithmetic_.multiply(buffers[1 - current], buffers[current], checked_base))
        return false;
      current = 1 - current; ++step;
    }
    return arithmetic_.copy(out, buffers[current]);
  }
  // Exact coefficient * 10^exponent; this is an auxiliary constant enclosure,
  // never a replacement of the original decimal-literal central arithmetic.
  bool decimal(Interval &out, std::int64_t coefficient, int exponent) noexcept {
    if (!owner_.guard()) return false;
    out.valid = false;
    if (exponent < -64 || exponent > 64) return owner_.fail(Cause::invalid_input);
    Interval c, ten, scale;
    if (!integer(c, coefficient)) return false;
    if (exponent == 0 || coefficient == 0) return arithmetic_.copy(out, c);
    if (!integer(ten, 10) || !power(scale, ten, unsigned(exponent < 0 ? -exponent : exponent)))
      return false;
    return exponent < 0 ? arithmetic_.divide(out, c, scale) :
                          arithmetic_.multiply(out, c, scale);
  }
  bool mpc(Interval &out) noexcept {
    // Preserve ancestor::mpc()'s exact AU/pi identity AND expression grouping:
    // 1e6 * 648000 / pi * 149597870700. The original nominal is untouched.
    if (!owner_.guard()) return false;
    out.valid = false;
    Interval million, arcseconds, au, pi_value, product, quotient;
    return integer(million, 1000000) && integer(arcseconds, 648000) &&
           integer(au, INT64_C(149597870700)) && pi(pi_value) &&
           arithmetic_.multiply(product, million, arcseconds) &&
           arithmetic_.divide(quotient, product, pi_value) &&
           arithmetic_.multiply(out, quotient, au);
  }

  // 32-term alternating arctan(q), q=1/d, and the positive next-term bound.
  // For d>=5: 0<q<1, terms decrease; after32 terms the next term is positive.
  // This exact rational theorem is independent of a libm atan candidate.
  bool arctan_inverse(Interval &out, std::int64_t denominator) noexcept {
    if (!owner_.guard()) return false;
    out.valid = false;
    if (denominator < 5) return owner_.fail(Cause::invalid_input);
    Interval q, q2, powers[2], sums[2], divisor, term, tail, upper_sum;
    if (!rational(q, 1, denominator) || !arithmetic_.multiply(q2, q, q) ||
        !arithmetic_.copy(powers[0], q) || !integer(sums[0], 0)) return false;
    unsigned current = 0, j = 0;
    for (;;) {
      if (!owner_.guard()) return false;
      if (j == 32) break;
      if (!owner_.iteration() || !integer(divisor, 2 * j + 1) ||
          !arithmetic_.divide(term, powers[current], divisor)) return false;
      if (!(j & 1U)) {
        if (!arithmetic_.add(sums[1 - current], sums[current], term)) return false;
      } else if (!arithmetic_.subtract(sums[1 - current], sums[current], term))
        return false;
      if (!arithmetic_.multiply(powers[1 - current], powers[current], q2)) return false;
      current = 1 - current; ++j;
    }
    // powers[current] is q^65. Account its exact remaining alternating tail.
    if (!integer(divisor, 65) || !arithmetic_.divide(tail, powers[current], divisor) ||
        !arithmetic_.add(upper_sum, sums[current], tail)) return false;
    return arithmetic_.enclose(out, sums[current].lower, upper_sum.upper);
  }
  bool pi(Interval &out) noexcept {
    if (!owner_.guard()) return false;
    out.valid = false;
    Interval a5, a239, sixteen, four, left, right;
    return arctan_inverse(a5, 5) && arctan_inverse(a239, 239) &&
           integer(sixteen, 16) && integer(four, 4) &&
           arithmetic_.multiply(left, sixteen, a5) &&
           arithmetic_.multiply(right, four, a239) &&
           arithmetic_.subtract(out, left, right);
  }
  // EXACT original declared central constants, not a newly qualified atomic
  // measurement/uncertainty law. Nominal ancestor symbols remain unchanged.
  bool atomic(AtomicNumbers &out) noexcept {
    Interval chiH_ev, chiHe_ev;
    return decimal(out.kb, INT64_C(1380649), -29) &&
           decimal(out.h, INT64_C(662607015), -42) &&
           integer(out.c, INT64_C(299792458)) &&
           decimal(out.ev, INT64_C(1602176634), -28) &&
           decimal(out.me, INT64_C(91093837139), -41) &&
           decimal(out.G, INT64_C(667430), -16) &&
           decimal(out.sigma, INT64_C(66524587051), -39) &&
           decimal(chiH_ev, INT64_C(13598434599702), -12) &&
           decimal(chiHe_ev, INT64_C(24587389011), -9) &&
           arithmetic_.multiply(out.chiH, chiH_ev, out.ev) &&
           arithmetic_.multiply(out.chiHe, chiHe_ev, out.ev);
  }
  // exp(-128*7/10)<2^-128 because exp(7/10)>Taylor4(7/10)
  // =482921/240000>2. A certified upper argument<=EXACT -448/5
  // therefore has this positive ABSOLUTE enclosure. Lower0 is never a
  // structural-zero function witness or an ordinary exp postcheck result.
  bool cold_exp(ColdImage &out, const Interval &argument) noexcept {
    static_assert(INT64_C(482921) > 2 * INT64_C(240000));
    if (!owner_.guard()) return false;
    if (&out.value == &argument) return owner_.fail(Cause::invalid_ownership);
    out.value.valid = false; out.absolute_cold = false;
    Interval checked, threshold;
    if (!arithmetic_.copy(checked, argument) || !rational(threshold, -448, 5) ||
        !owner_.guard()) return false;
    if (checked.upper > threshold.lower) return owner_.fail(Cause::invalid_input);
    if (!owner_.action(Action::floating(2, 0, 2))) return false;
    W zero, ceiling; zero = 0; ceiling = 0x1p-128L;
    if (!arithmetic_.enclose(out.value, zero, ceiling)) return false;
    out.absolute_cold = true;
    return true;
  }
};
} // namespace conditional_drag_reference_certificate
