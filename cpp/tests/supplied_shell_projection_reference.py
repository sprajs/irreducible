"""Bounded independent angular reference for the supplied finite radial law.

Source implementation; execution needs the repository owner's separate lease.
No native Bessel code, installed table, binary launcher or physical CMB source.
"""
from __future__ import annotations

from dataclasses import dataclass
from decimal import (Decimal, Context, ROUND_FLOOR, ROUND_CEILING,
                     InvalidOperation, DivisionByZero, Overflow, Underflow,
                     Subnormal, Clamped, Inexact, Rounded)
import json
import struct
import sys
import decimal


SCHEMA = "irred/supplied-shell-projection-reference/v1"
MAX_BITS = 16384
MAX_WORK = 40000000
MAX_PAYLOAD = 8388608
MAX_SOURCE_BYTES = 65536
MAX_REPORT_BYTES = 262144
NUMERIC_SLOT = 4096
WRAPPER_SLOT = 256
SCRATCH_BYTES = 256 * NUMERIC_SLOT + 64 * WRAPPER_SLOT
ZERO, ONE = Decimal(0), Decimal(1)
ABS_ALLOC = Decimal("1e-10")
REL_ALLOC = Decimal("1e-10")
TAIL_GOAL = Decimal("1e-50")


class Refusal(Exception):
    """A retained failure; no budget or source changes are implicit."""


def fixed_keys(value, required, optional=()):
    if not isinstance(value, dict) or not set(required) <= value.keys() or any(
            key not in required and key not in optional for key in value):
        raise Refusal("unknown/missing fixed schema member")


class Work:
    limits = {"root_bracket_checks": 8192, "panel_node_visits": 262144,
              "trig_terms": 33554432}

    def __init__(self):
        self.counts = {}
        self.total = 0
        self.rejected = {}
        self.rounding_events = {"Inexact": 0, "Rounded": 0}
        self.trapped_events = {}
        self.source_bytes = 0
        self.rejected_bytes = None
        self.stage = "source-admission"
        self.last_attempt = None

    def charge(self, name, amount=1):
        if not isinstance(amount, int) or amount < 0:
            raise Refusal("invalid work charge")
        previous = self.counts.get(name, 0)
        if amount > MAX_WORK - self.total or amount > self.limits.get(name, MAX_WORK) - previous:
            self.rejected[name] = self.rejected.get(name, 0) + 1
            raise Refusal("work cap: " + name)
        self.total += amount
        self.counts[name] = previous + amount
        self.last_attempt = {"unit": name, "amount": amount, "total_after_charge": self.total}

    def bytes(self, amount):
        if type(amount) is not int or amount < 0 or amount > MAX_SOURCE_BYTES - self.source_bytes:
            self.rejected_bytes = {"requested": amount, "remaining": MAX_SOURCE_BYTES - self.source_bytes}
            raise Refusal("inspected source byte cap")
        self.source_bytes += amount


class Payload:
    """Prospective reservations include old/new objects and serial temporaries."""
    def __init__(self, raw_bytes):
        self.input_bytes = 256 * (raw_bytes + 1) + 65536
        self.persistent = 0
        self.report_bytes = 0
        self.peak = 0
        self.check()

    def check(self, extra=0):
        live = (self.input_bytes + SCRATCH_BYTES + self.persistent +
                8 * self.report_bytes + 65536 + extra)
        if live > MAX_PAYLOAD:
            raise Refusal("whole-live payload cap")
        self.peak = max(self.peak, live)

    def hold(self, numeric=0, wrappers=0):
        amount = numeric * NUMERIC_SLOT + wrappers * WRAPPER_SLOT
        self.check(amount)
        self.persistent += amount
        return amount

    def release(self, amount):
        if amount < 0 or amount > self.persistent:
            raise Refusal("payload ownership underflow")
        self.persistent -= amount

    def prepare_report(self, text_bytes):
        if type(text_bytes) is not int or text_bytes < 0:
            raise Refusal("invalid prospective report text bound")
        predicted = 8192 + 6 * text_bytes
        if predicted > 16384 or predicted > MAX_REPORT_BYTES - self.report_bytes:
            raise Refusal("prospective report record byte cap")
        self.check(8 * predicted)
        return predicted

    def append_report(self, rows, record, text_bytes=None):
        # Bound JSON escaping of the retained input identifier before creating
        # encoder buffers. Other fields have bounded100-digit scalar strings.
        if text_bytes is None:
            text_bytes = len(record.get("case_id", "").encode("utf-8"))
        predicted = self.prepare_report(text_bytes)
        encoded = json.dumps(record, separators=(",", ":"), allow_nan=False)
        size = len(encoded.encode("utf-8")) + 1
        if size > predicted or size > MAX_REPORT_BYTES - self.report_bytes:
            raise Refusal("retained report byte cap")
        self.report_bytes += size
        self.check()
        rows.append(record)


@dataclass(frozen=True, slots=True)
class Rational:
    n: int
    d: int = 1


@dataclass(frozen=True, slots=True)
class Interval:
    lo: Decimal
    hi: Decimal


@dataclass(frozen=True, slots=True)
class Node:
    mu: Interval
    weight: Interval


@dataclass(frozen=True, slots=True)
class Ball:
    center: Decimal
    root: Decimal
    trig: Decimal
    arithmetic: Decimal
    quadrature: Decimal = ZERO
    refinement: Decimal = ZERO


class Arithmetic:
    def __init__(self, work):
        self.work = work
        self.down = Context(prec=100, rounding=ROUND_FLOOR,
                            Emin=-999999, Emax=999999, clamp=0)
        self.up = Context(prec=100, rounding=ROUND_CEILING,
                          Emin=-999999, Emax=999999, clamp=0)
        for context in (self.down, self.up):
            for signal in (InvalidOperation, DivisionByZero, Overflow,
                           Underflow, Subnormal, Clamped):
                context.traps[signal] = True
            context.traps[Inexact] = context.traps[Rounded] = False

    def integer(self, x):
        if not isinstance(x, int) or x.bit_length() > MAX_BITS or sys.getsizeof(x) > NUMERIC_SLOT:
            raise Refusal("integer width/object cap")
        return x

    def int_op(self, operation, a, b):
        self.integer(a)
        self.integer(b)
        if operation == "mul":
            width = 0 if a == 0 or b == 0 else a.bit_length() + b.bit_length()
        elif operation in ("add", "sub"):
            width = max(a.bit_length(), b.bit_length()) + 1
        else:
            width = max(a.bit_length(), b.bit_length())
        if width > MAX_BITS:
            raise Refusal("prospective integer temporary width cap")
        self.work.charge("exact_integer_operations")
        if operation == "mul":
            answer = a * b
        elif operation == "add":
            answer = a + b
        elif operation == "sub":
            answer = a - b
        elif operation == "div":
            if b == 0:
                raise Refusal("integer division by zero")
            answer = a // b
        elif operation == "mod":
            if b == 0:
                raise Refusal("integer remainder by zero")
            answer = a % b
        else:
            raise Refusal("unknown integer operation")
        return self.integer(answer)

    def gcd(self, a, b):
        a, b = abs(self.integer(a)), abs(self.integer(b))
        while b:
            a, b = b, self.int_op("mod", a, b)
        return a

    def rat(self, n, d=1):
        self.integer(n)
        self.integer(d)
        if d == 0:
            raise Refusal("rational zero denominator")
        if d < 0:
            n, d = -n, -d
        if n == 0:
            return Rational(0)
        g = self.gcd(n, d)
        return Rational(self.int_op("div", n, g), self.int_op("div", d, g))

    def rat_op(self, op, a, b):
        if op in ("add", "sub"):
            left = self.int_op("mul", a.n, b.d)
            right = self.int_op("mul", b.n, a.d)
            return self.rat(self.int_op(op, left, right), self.int_op("mul", a.d, b.d))
        if op == "mul":
            return self.rat(self.int_op("mul", a.n, b.n), self.int_op("mul", a.d, b.d))
        if op == "div":
            return self.rat(self.int_op("mul", a.n, b.d), self.int_op("mul", a.d, b.n))
        raise Refusal("unknown rational operation")

    def decimal(self, x):
        if not x.is_finite() or sys.getsizeof(x) > NUMERIC_SLOT:
            raise Refusal("Decimal finite/object cap")
        return x

    def operation(self, name, x, y, upward):
        self.work.charge("directed_interval_operations")
        context = self.up if upward else self.down
        context.clear_flags()
        try:
            answer = getattr(context, name)(x, y)
        except ArithmeticError as error:
            kind = type(error).__name__
            self.work.trapped_events[kind] = self.work.trapped_events.get(kind, 0) + 1
            raise Refusal("directed Decimal trap: " + type(error).__name__) from error
        finally:
            for signal in (Inexact, Rounded):
                if context.flags[signal]:
                    self.work.rounding_events[signal.__name__] += 1
        return self.decimal(answer)

    def point(self, x):
        return Interval(self.decimal(x), x)

    def add(self, a, b):
        return Interval(self.operation("add", a.lo, b.lo, False),
                        self.operation("add", a.hi, b.hi, True))

    def sub(self, a, b):
        return Interval(self.operation("subtract", a.lo, b.hi, False),
                        self.operation("subtract", a.hi, b.lo, True))

    def neg(self, a):
        return Interval(a.hi.copy_negate(), a.lo.copy_negate())

    def mul(self, a, b):
        pairs = ((a.lo, b.lo), (a.lo, b.hi), (a.hi, b.lo), (a.hi, b.hi))
        if a.lo == a.hi == ZERO or b.lo == b.hi == ZERO:
            return self.point(ZERO)
        lower = [self.operation("multiply", x, y, False) for x, y in pairs]
        upper = [self.operation("multiply", x, y, True) for x, y in pairs]
        return Interval(min(lower), max(upper))

    def div(self, a, b):
        if b.lo <= ZERO <= b.hi:
            raise Refusal("interval denominator contains zero")
        pairs = ((a.lo, b.lo), (a.lo, b.hi), (a.hi, b.lo), (a.hi, b.hi))
        lower = [self.operation("divide", x, y, False) for x, y in pairs]
        upper = [self.operation("divide", x, y, True) for x, y in pairs]
        return Interval(min(lower), max(upper))

    def abs_bound(self, a):
        return max(a.lo.copy_abs(), a.hi.copy_abs())

    def center_radius(self, a):
        center = self.operation("divide", self.operation("add", a.lo, a.hi, False), Decimal(2), False)
        # An exact input Decimal can have more than100 digits. Its rounded
        # midpoint can lie just outside even a point interval; signed directed
        # subtraction, rather than abs of an upward-rounded negative value,
        # is therefore mandatory here.
        return center, self.abs_bound(self.sub(self.point(center), a))

    def upper_add(self, *values):
        answer = ZERO
        for value in values:
            if value < ZERO:
                raise Refusal("negative upper radius")
            answer = self.operation("add", answer, value, True)
        return answer

    def upper_mul(self, x, y):
        if x < ZERO or y < ZERO:
            raise Refusal("negative upper-radius factor")
        return self.operation("multiply", x, y, True)

    def rational_interval(self, value):
        return self.div(self.point(Decimal(value.n)), self.point(Decimal(value.d)))

    def profile_controls(self):
        self.work.stage = "directed-Decimal-profile-controls"
        self.work.charge("comparison_cast_operations")
        tiny = Decimal("1e-100")
        next_one = Decimal("1." + "0" * 98 + "1")
        positive_down = self.operation("add", ONE, tiny, False)
        positive_up = self.operation("add", ONE, tiny, True)
        negative_down = self.operation("add", ONE.copy_negate(), tiny.copy_negate(), False)
        negative_up = self.operation("add", ONE.copy_negate(), tiny.copy_negate(), True)
        if ((positive_down, positive_up) != (ONE, next_one) or
                (negative_down, negative_up) != (next_one.copy_negate(), ONE.copy_negate())):
            raise Refusal("signed directed Decimal profile control")
        cancellation = self.operation("subtract", ONE, ONE, False)
        if cancellation != ZERO:
            raise Refusal("exact Decimal cancellation control")
        subnormal_trapped = False
        try:
            self.operation("multiply", Decimal("1e-999999"), Decimal("0.1"), False)
        except Refusal as error:
            if not isinstance(error.__cause__, Subnormal):
                raise
            subnormal_trapped = True
        if not subnormal_trapped:
            raise Refusal("Decimal subnormal trap control")
        return {"positive_down": str(positive_down), "positive_up": str(positive_up),
                "negative_down": str(negative_down), "negative_up": str(negative_up),
                "cancellation_exact_zero": True, "exact_subnormal_trapped": True}


def binary64(bits, arithmetic):
    arithmetic.work.charge("source_scalar_inspections")
    if not isinstance(bits, str) or len(bits) != 16:
        raise Refusal("invalid IEEE64 bit string")
    arithmetic.work.bytes(len(bits))
    if any(c not in "0123456789abcdef" for c in bits):
        raise Refusal("invalid IEEE64 bit string")
    arithmetic.work.charge("comparison_cast_operations")
    value = struct.unpack(">d", bytes.fromhex(bits))[0]
    exact = Decimal.from_float(value)
    arithmetic.decimal(exact)
    numerator, denominator = value.as_integer_ratio()
    return exact, arithmetic.rat(numerator, denominator)


def _trim(poly):
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return poly


def integer_legendre(a):
    a.work.charge("source_tag_inspections")
    basis = [[1], [0, 2]]
    for n in range(1, 16):
        previous, current = basis[-2], basis[-1]
        coefficients = [0] * (n + 2)
        first = a.int_op("mul", 2, a.int_op("add", a.int_op("mul", 2, n), 1))
        second = a.int_op("mul", 4, n)
        divisor = a.int_op("add", n, 1)
        for i in range(n + 2):
            left = a.int_op("mul", first, current[i - 1]) if 0 < i <= len(current) else 0
            right = a.int_op("mul", second, previous[i]) if i < len(previous) else 0
            numerator = a.int_op("sub", left, right)
            if a.int_op("mod", numerator, divisor) != 0:
                raise Refusal("Legendre integer divisibility")
            coefficients[i] = a.int_op("div", numerator, divisor)
        basis.append(coefficients)
    return tuple(tuple(poly) for poly in basis)


def primitive(poly, a):
    common = 1
    for item in poly:
        g = a.gcd(common, item.d)
        common = a.int_op("mul", a.int_op("div", common, g), item.d)
    values = [a.int_op("mul", item.n, a.int_op("div", common, item.d)) for item in poly]
    g = 0
    for value in values:
        g = a.gcd(g, value)
    if g:
        values = [a.int_op("div", value, g) for value in values]
    return _trim(values)


def sturm_sequence(poly, a):
    derivative = [a.int_op("mul", i, poly[i]) for i in range(1, len(poly))]
    sequence = [poly, derivative]
    while len(sequence[-1]) > 1:
        dividend = [Rational(x) for x in sequence[-2]]
        divisor = [Rational(x) for x in sequence[-1]]
        remainder = dividend[:]
        while len(remainder) >= len(divisor) and any(x.n for x in remainder):
            offset = len(remainder) - len(divisor)
            factor = a.rat_op("div", remainder[-1], divisor[-1])
            for i, coefficient in enumerate(divisor):
                product = a.rat_op("mul", factor, coefficient)
                remainder[i + offset] = a.rat_op("sub", remainder[i + offset], product)
            while len(remainder) > 1 and remainder[-1].n == 0:
                remainder.pop()
        negative = [Rational(-x.n, x.d) for x in remainder]
        following = primitive(negative, a)
        if not any(following) or len(following) >= len(sequence[-1]):
            raise Refusal("Sturm nondecreasing/zero remainder")
        sequence.append(following)
    if sequence[-1][0] == 0:
        raise Refusal("Sturm terminal polynomial")
    return sequence


def dyadic_sign(poly, numerator, denominator, a):
    a.work.charge("polynomial_evaluations")
    result, power = poly[-1], denominator
    for coefficient in reversed(poly[:-1]):
        result = a.int_op("add", a.int_op("mul", result, numerator),
                          a.int_op("mul", coefficient, power))
        power = a.int_op("mul", power, denominator)
    return (result > 0) - (result < 0)


def variations(sequence, n, d, a):
    signs = [dyadic_sign(poly, n, d, a) for poly in sequence]
    signs = [value for value in signs if value]
    return sum(signs[i] != signs[i - 1] for i in range(1, len(signs)))


def polynomial(poly, scale, argument, a):
    a.work.charge("polynomial_evaluations")
    result = a.point(ZERO)
    denominator = a.point(Decimal(scale))
    for coefficient in reversed(poly):
        c = a.div(a.point(Decimal(coefficient)), denominator)
        result = a.add(a.mul(result, argument), c)
    return result


def prepare_nodes(basis, a):
    a.work.stage = "Sturm-construction"
    sequence = sturm_sequence(basis[16], a)
    a.work.stage = "Sturm-root-isolation"
    for endpoint in (-1, 1):
        a.work.charge("root_bracket_checks")
        if dyadic_sign(basis[16], endpoint, 1, a) == 0:
            raise Refusal("GL16 full-bracket endpoint root")
    a.work.charge("root_bracket_checks")
    previous = variations(sequence, -1, 1, a)
    terminal = variations(sequence, 1, 1, a)
    if previous - terminal != 16:
        raise Refusal("Sturm full root count is not16")
    brackets = []
    left, denominator = -64, 64
    for right in range(-63, 65):
        a.work.charge("root_bracket_checks")
        current = variations(sequence, right, denominator, a)
        count = previous - current
        if count not in (0, 1):
            raise Refusal("unisolated Sturm panel")
        if count:
            a.work.charge("root_bracket_checks")
            sign_left = dyadic_sign(basis[16], left, denominator, a)
            a.work.charge("root_bracket_checks")
            sign_right = dyadic_sign(basis[16], right, denominator, a)
            a.work.charge("root_bracket_checks")
            if sign_left * sign_right != -1:
                raise Refusal("root panel endpoint/sign failure")
            brackets.append((left, right, denominator, sign_left))
        left, previous = right, current
    if len(brackets) != 16:
        raise Refusal("incomplete disjoint root brackets")
    derivative = [a.int_op("mul", i, basis[16][i]) for i in range(1, 17)]
    nodes = []
    for index, (lo, hi, d, sign_lo) in enumerate(brackets):
        a.work.stage = "root-refinement:" + str(index)
        attempts = 0
        while d.bit_length() - 1 < 320:
            if attempts >= 384:
                raise Refusal("root refinement per-bracket cap")
            a.work.charge("root_bracket_checks")
            attempts += 1
            middle = a.int_op("add", lo, hi)
            next_d = a.int_op("mul", d, 2)
            sign = dyadic_sign(basis[16], middle, next_d, a)
            if sign == 0:
                raise Refusal("unexpected rational GL16 root")
            if sign == sign_lo:
                lo, hi = middle, a.int_op("mul", hi, 2)
            else:
                lo, hi = a.int_op("mul", lo, 2), middle
            d = next_d
        lower = a.rational_interval(Rational(lo, d))
        upper = a.rational_interval(Rational(hi, d))
        z = Interval(lower.lo, upper.hi)
        dp = polynomial(derivative, 65536, z, a)
        if dp.lo <= ZERO <= dp.hi:
            raise Refusal("GL derivative interval contains zero")
        a.work.charge("weight_calculations")
        denominator_i = a.mul(a.sub(a.point(ONE), a.mul(z, z)), a.mul(dp, dp))
        weight = a.div(a.point(Decimal(2)), denominator_i)
        if weight.lo <= ZERO:
            raise Refusal("GL weight nonpositive")
        mu = a.div(a.add(z, a.point(ONE)), a.point(Decimal(2)))
        half_weight = a.div(weight, a.point(Decimal(2)))
        if not ZERO < mu.lo <= mu.hi < ONE:
            raise Refusal("mapped GL node domain")
        nodes.append(Node(mu, half_weight))
    return tuple(nodes)


def trig(argument, odd, a):
    if argument == ZERO:
        return (ZERO if odd else ONE), ZERO, ZERO
    x = a.point(argument)
    square = a.mul(x, x)
    term = x if odd else a.point(ONE)
    total = a.point(ZERO)
    for n in range(128):
        a.work.charge("trig_terms")
        twice = a.int_op("mul", 2, n)
        denominator = a.int_op("mul", a.int_op("add", twice, 2),
                               a.int_op("add", twice, 3 if odd else 1))
        ratio = a.div(square, a.point(Decimal(denominator)))
        # Current term is the charged first omitted term once the decreasing
        # tail is earned. Otherwise it remains an included term.
        if n >= 1 and ratio.hi < ONE and a.abs_bound(term) <= TAIL_GOAL:
            center, arithmetic = a.center_radius(total)
            return center, a.abs_bound(term), arithmetic
        total = a.add(total, term)
        if n == 127:
            raise Refusal("trig term/tail cap")
        term = a.neg(a.mul(term, ratio))
    raise Refusal("trig term/tail cap")


def factorial(n, a):
    result = 1
    for value in range(2, n + 1):
        result = a.int_op("mul", result, value)
    return result


def integer_power(base, exponent, a):
    result = 1
    for _ in range(exponent):
        result = a.int_op("mul", result, base)
    return result


def quadrature_envelope(ell, panels, basis, a):
    numerator = 0
    for r in range(ell + 1):
        derivative_sum = 0
        for m in range(r, ell + 1):
            falling = 1
            for i in range(r):
                falling = a.int_op("mul", falling, m - i)
            derivative_sum = a.int_op("add", derivative_sum,
                                    a.int_op("mul", abs(basis[ell][m]), falling))
        choose = a.int_op("div", factorial(32, a),
                          a.int_op("mul", factorial(r, a), factorial(32 - r, a)))
        component = a.int_op("mul", a.int_op("mul", choose, derivative_sum),
                             integer_power(8, 32 - r, a))
        numerator = a.int_op("add", numerator, component)
    factor = integer_power(factorial(16, a), 4, a)
    numerator = a.int_op("mul", numerator, factor)
    denominator = a.int_op("mul", a.int_op("mul", 33,
                            integer_power(factorial(32, a), 3, a)),
                            a.int_op("mul", integer_power(2, ell, a),
                                     integer_power(panels, 32, a)))
    return a.rational_interval(Rational(numerator, denominator)).hi


def angular(x, ell, panels, nodes, basis, a):
    if x == ZERO:
        return Ball(ONE if ell == 0 else ZERO, ZERO, ZERO, ZERO)
    derivative = 0
    for m in range(1, ell + 1):
        derivative = a.int_op("add", derivative, a.int_op("mul", abs(basis[ell][m]), m))
    scale = integer_power(2, ell, a)
    lipschitz = a.upper_add(a.rational_interval(Rational(derivative, scale)).hi, Decimal(8))
    summed = a.point(ZERO)
    root_error = trig_error = arithmetic_error = ZERO
    for panel in range(panels):
        for node in nodes:
            a.work.charge("panel_node_visits")
            mu = a.div(a.add(a.point(Decimal(panel)), node.mu), a.point(Decimal(panels)))
            weight = a.div(node.weight, a.point(Decimal(panels)))
            muc, mur = a.center_radius(mu)
            wc, wr = a.center_radius(weight)
            if not ZERO <= muc <= ONE:
                raise Refusal("reference point outside angular panel")
            angle = a.mul(a.point(x), a.point(muc))
            tc_arg, arg_error = a.center_radius(angle)
            tc, tail, tround = trig(tc_arg, bool(ell % 2), a)
            p = polynomial(basis[ell], scale, a.point(muc), a)
            pc, pr = a.center_radius(p)
            a.work.charge("weighted_products")
            product = a.mul(a.mul(a.point(wc), a.point(pc)), a.point(tc))
            product_center, product_error = a.center_radius(product)
            a.work.charge("signed_additions")
            summed = a.add(summed, a.point(product_center))
            root_error = a.upper_add(root_error,
                a.upper_mul(a.abs_bound(weight), a.upper_mul(lipschitz, mur)), wr)
            poly_bound = a.upper_add(pc.copy_abs(), pr)
            trig_error = a.upper_add(trig_error,
                a.upper_mul(wc.copy_abs(), a.upper_mul(poly_bound, tail)))
            arithmetic_error = a.upper_add(arithmetic_error, product_error,
                a.upper_mul(wc.copy_abs(), a.upper_add(
                    a.upper_mul(poly_bound, a.upper_add(tround, arg_error)),
                    a.upper_mul(tc.copy_abs(), pr))))
    center, sum_error = a.center_radius(summed)
    if (ell // 2) % 2:
        center = center.copy_negate()
    return Ball(center, root_error, trig_error,
                a.upper_add(arithmetic_error, sum_error),
                quadrature_envelope(ell, panels, basis, a))


def assemble(shells, ell, panels, nodes, basis, a):
    summed = a.point(ZERO)
    root = tail = rounding = quadrature = ZERO
    for amplitude, phase in shells:
        a.work.charge("shell_visits")
        xc, xr = a.center_radius(phase)
        if not ZERO <= xc <= Decimal(8):
            raise Refusal("reference phase center domain")
        if amplitude == ZERO:
            continue
        value = angular(xc, ell, panels, nodes, basis, a)
        a.work.charge("weighted_products")
        product = a.mul(a.point(amplitude), a.point(value.center))
        product_center, product_error = a.center_radius(product)
        a.work.charge("signed_additions")
        summed = a.add(summed, a.point(product_center))
        absolute = amplitude.copy_abs()
        root = a.upper_add(root, a.upper_mul(absolute, value.root))
        tail = a.upper_add(tail, a.upper_mul(absolute, value.trig))
        rounding = a.upper_add(rounding, product_error,
            a.upper_mul(absolute, a.upper_add(value.arithmetic,
                a.operation("divide", xr, Decimal(2), True))))
        quadrature = a.upper_add(quadrature, a.upper_mul(absolute, value.quadrature))
    center, error = a.center_radius(summed)
    return Ball(center, root, tail, a.upper_add(rounding, error), quadrature)


def radius(ball, a):
    return a.upper_add(ball.root, ball.trig, ball.arithmetic,
                       ball.quadrature, ball.refinement)


def scalar_identity(shells, ell, a):
    """Original j0/j1 closed forms, independent of Legendre/GL16."""
    summed = a.point(ZERO)
    argument_error = ZERO
    for amplitude, phase in shells:
        a.work.charge("shell_visits")
        if amplitude == ZERO:
            continue
        x, xr = a.center_radius(phase)
        if x == ZERO:
            value = a.point(ONE if ell == 0 else ZERO)
        else:
            sc, st, sr = trig(x, True, a)
            se = a.upper_add(st, sr)
            sin_i = Interval(a.operation("subtract", sc, se, False),
                             a.operation("add", sc, se, True))
            if ell == 0:
                value = a.div(sin_i, a.point(x))
            else:
                cc, ct, cr = trig(x, False, a)
                ce = a.upper_add(ct, cr)
                cos_i = Interval(a.operation("subtract", cc, ce, False),
                                 a.operation("add", cc, ce, True))
                value = a.div(a.sub(sin_i, a.mul(a.point(x), cos_i)),
                              a.mul(a.point(x), a.point(x)))
        a.work.charge("weighted_products")
        product = a.mul(a.point(amplitude), value)
        a.work.charge("signed_additions")
        summed = a.add(summed, product)
        argument_error = a.upper_add(argument_error,
            a.upper_mul(amplitude.copy_abs(), a.operation("divide", xr, Decimal(2), True)))
    center, error = a.center_radius(summed)
    return center, a.upper_add(error, argument_error)


def compare(native_value, native_radius, refs, a):
    a.work.charge("comparison_cast_operations")
    fine, medium, coarse = refs[2], refs[1], refs[0]
    change1 = a.abs_bound(a.sub(a.point(fine.center), a.point(medium.center)))
    change2 = a.abs_bound(a.sub(a.point(medium.center), a.point(coarse.center)))
    refinement = max(change1, change2)
    reference = Ball(fine.center, fine.root, fine.trig, fine.arithmetic,
                     fine.quadrature, refinement)
    er = radius(reference, a)
    allocation_i = a.add(a.point(ABS_ALLOC),
                         a.mul(a.point(REL_ALLOC), a.point(reference.center.copy_abs())))
    allocation = allocation_i.lo
    reference_limit = a.operation("multiply", Decimal("0.05"), allocation, False)
    native_limit = a.operation("multiply", Decimal("0.90"), allocation, False)
    ref_fractions = ((reference.root, Decimal("0.15")),
                     (reference.trig, Decimal("0.15")),
                     (reference.arithmetic, Decimal("0.20")),
                     (a.upper_add(reference.quadrature, reference.refinement), Decimal("0.50")))
    ingredient_pass = all(value <= a.operation("multiply", share, reference_limit, False)
                          for value, share in ref_fractions)
    discrepancy = a.abs_bound(a.sub(a.point(native_value), a.point(reference.center)))
    complete = a.upper_add(discrepancy, native_radius, er)
    # Upper discrepancy versus a LOWER radius sum proves consistency; an outer
    # endpoint intersection that merely might overlap cannot pass this gate.
    radius_sum_lower = a.operation("add", native_radius, er, False)
    gates = {"reference_allocation": er <= reference_limit,
             "native_allocation": native_radius <= native_limit,
             "combined_allocation": complete <= allocation,
             "interval_consistency": discrepancy <= radius_sum_lower,
             "ingredient_allocations": ingredient_pass}
    return {"status": "ok" if all(gates.values()) else "comparison_failed",
            "gates": gates, "center": str(reference.center),
            "radius": str(er), "allocation_lower": str(allocation),
            "discrepancy_upper": str(discrepancy),
            "root_node_weight": str(reference.root), "trig": str(reference.trig),
            "arithmetic": str(reference.arithmetic),
            "quadrature": str(reference.quadrature), "refinement": str(refinement),
            "q_centers": [str(item.center) for item in refs],
            "q_radii": [str(radius(item, a)) for item in refs]}


def inspect_label(label, a):
    a.work.charge("source_tag_inspections")
    if not isinstance(label, str) or not label or "\x00" in label:
        raise Refusal("empty/NUL/invalid identifier or provenance label")
    try:
        size = len(label.encode("utf-8", errors="strict"))
    except UnicodeError as error:
        raise Refusal("invalid UTF8 identifier or provenance label") from error
    a.work.bytes(size)
    return size


def validate_source(source, ells, a):
    fixed_keys(source, ("identity", "k_ids", "shell_ids", "k_bits", "chi_bits",
                       "amplitude_bits", "mode", "units", "convention", "roles",
                       "mode_origin", "ordering_provenance"),
               ("endpoint",))
    if not isinstance(source, dict) or not isinstance(ells, list) or not 1 <= len(ells) <= 9:
        raise Refusal("source/request shape")
    tags = {"mode": "scalar-unit-zeta", "units": "comoving-Mpc-no-h",
            "convention": "outward-exp-plus-i"}
    for key, expected in tags.items():
        a.work.charge("source_tag_inspections")
        if source.get(key) != expected:
            raise Refusal("source tag: " + key)
    for ell in ells:
        a.work.charge("source_scalar_inspections")
        if type(ell) is not int or not 0 <= ell <= 8:
            raise Refusal("requested multipole domain")
    ks, chis, amplitudes = source.get("k_bits"), source.get("chi_bits"), source.get("amplitude_bits")
    if not all(isinstance(value, list) for value in (ks, chis, amplitudes)):
        raise Refusal("source table types")
    if not 1 <= len(ks) <= 16 or not 1 <= len(chis) <= 16 or len(amplitudes) != len(ks) * len(chis):
        raise Refusal("bounded reference shape")
    for key, count in (("k_ids", len(ks)), ("shell_ids", len(chis)), ("roles", len(chis))):
        if not isinstance(source.get(key), list) or len(source[key]) != count:
            raise Refusal("source axis association")
    labels = ([source["identity"], source["mode_origin"], source["ordering_provenance"]] +
              source["k_ids"] + source["shell_ids"])
    endpoint = source.get("endpoint")
    if endpoint is not None:
        fixed_keys(endpoint, ("producer_identity", "coordinate_identity", "unit_identity",
                              "support_lower_bits", "support_upper_bits", "survival_status",
                              "survival_value_bits", "survival_radius_bits"))
        labels += [endpoint[key] for key in ("producer_identity", "coordinate_identity", "unit_identity")]
    provenance_bytes = sum(inspect_label(label, a) for label in labels)
    if endpoint is not None:
        lower, _ = binary64(endpoint["support_lower_bits"], a)
        upper, _ = binary64(endpoint["support_upper_bits"], a)
        a.work.charge("comparison_cast_operations")
        if lower > upper:
            raise Refusal("endpoint support ordering")
        provenance_bytes += 32
        for key in ("survival_value_bits", "survival_radius_bits"):
            if endpoint[key] is not None:
                value, _ = binary64(endpoint[key], a)
                provenance_bytes += 16
                a.work.charge("comparison_cast_operations")
                if value < ZERO or (key == "survival_value_bits" and value > ONE):
                    raise Refusal("endpoint survival provenance domain")
        a.work.charge("source_scalar_inspections")
        if type(endpoint["survival_status"]) is not int or not 0 <= endpoint["survival_status"] <= 8:
            raise Refusal("endpoint survival status")
    for role in source["roles"]:
        a.work.charge("source_tag_inspections")
        if role not in ("radial-atom", "boundary-response"):
            raise Refusal("source role")
    return provenance_bytes


def decode_source(source, a):
    # Complete structure/provenance was admitted and retained first.
    ks, chis, amplitudes = source["k_bits"], source["chi_bits"], source["amplitude_bits"]
    kd = [binary64(item, a) for item in ks]
    cd = [binary64(item, a) for item in chis]
    ad = [binary64(item, a) for item in amplitudes]
    if any(value < ZERO for value, _ in kd + cd) or any(value.copy_abs() > ONE for value, _ in ad):
        raise Refusal("source scalar domain")
    result = []
    for ki, (k, kr) in enumerate(kd):
        shells = []
        for si, (chi, cr) in enumerate(cd):
            a.work.charge("phase_products")
            exact = a.rat_op("mul", kr, cr)
            if exact.n != 0:
                lower_left = a.int_op("mul", exact.n, 1 << 32)
                upper_right = a.int_op("mul", exact.d, 8)
                if lower_left < exact.d or exact.n > upper_right:
                    raise Refusal("exact phase domain")
            phase = a.mul(a.point(k), a.point(chi))
            shells.append((ad[ki * len(chis) + si][0], phase))
        result.append(shells)
    return tuple(tuple(row) for row in result)


class Reference:
    def __init__(self, work, payload):
        self.work, self.payload = work, payload
        self.arithmetic = Arithmetic(work)
        self.root_reservation = payload.hold(numeric=512, wrappers=128)
        self.basis = self.nodes = None
        self.context = None

    def case(self, case, records, provenance):
        self.work.stage = "source-admission"
        case_id = case.get("id") if isinstance(case, dict) else None
        self.context = {"case_id": case_id if isinstance(case_id, str) else None,
                        "stage": "source-admission"}
        source = case.get("source") if isinstance(case, dict) else None
        if not isinstance(source, dict):
            raise Refusal("source object")
        ks, cs = source.get("k_bits"), source.get("chi_bits")
        if not isinstance(ks, list) or not isinstance(cs, list) or not 1 <= len(ks) <= 16 or not 1 <= len(cs) <= 16:
            raise Refusal("source structural shape")
        # Original graph remains reserved while the once-decoded immutable
        # typed source is acquired. Old/new arrays, all Rational components,
        # phase intervals and variable pointer tables coexist in this bound.
        k, s = len(ks), len(cs)
        held = self.payload.hold(numeric=3 * (k + s + k * s) + 2 * k * s,
                                 wrappers=8 * (k + s + k * s) + 64)
        try:
            self._case(case, records, provenance)
        finally:
            self.payload.release(held)

    def _case(self, case, records, provenance):
        a = self.arithmetic
        fixed_keys(case, ("id", "source", "requested_ell", "native"))
        if not isinstance(case, dict) or not isinstance(case.get("id"), str) or not case["id"]:
            raise Refusal("case identity")
        self.context = {"case_id": case["id"], "stage": "source-admission"}
        case_bytes = inspect_label(case["id"], a)
        ells = case.get("requested_ell")
        source = case["source"]
        provenance_bytes = validate_source(source, ells, a)
        self.payload.prepare_report(case_bytes + provenance_bytes)
        # Shared references to the original object preserve literal failed
        # endpoint metadata without copying the source model into each row.
        self.payload.append_report(provenance, {
            "case_id": case["id"], "source_identity": source["identity"],
            "mode_origin": source["mode_origin"],
            "ordering_provenance": source["ordering_provenance"],
            "endpoint": source.get("endpoint")}, case_bytes + provenance_bytes)
        shells = decode_source(source, a)
        native = case.get("native")
        fixed_keys(native, ("prepare_status", "batch_status", "rows", "work", "payload_bytes"),
                   ("requested_accuracy_reduction_power",))
        if (type(native["prepare_status"]) is not int or native["prepare_status"] != 0 or
                type(native["batch_status"]) is not int or native["batch_status"] != 0):
            raise Refusal("native preparation/batch refused")
        native_work = native.get("work")
        nonbyte = ("series", "radius", "source_scalar_inspections", "source_tag_inspections",
                   "shell_visits", "phase_products", "Bessel_evaluations", "shell_products",
                   "signed_additions", "output_casts")
        fixed_keys(native_work, nonbyte + ("total", "source_bytes_inspected", "rejected_work_requests"))
        total = 0
        for name in nonbyte + ("total", "source_bytes_inspected", "rejected_work_requests"):
            a.work.charge("source_scalar_inspections")
            count = native_work[name]
            if not isinstance(name, str) or type(count) is not int or not 0 <= count <= 8000000:
                raise Refusal("native named work receipt")
            if name in nonbyte:
                a.work.charge("comparison_cast_operations")
                if count > 8000000 - total:
                    raise Refusal("native complete work sum cap")
                total += count
        a.work.charge("comparison_cast_operations")
        if (native_work["series"] > 6000000 or total != native_work["total"] or
                native_work["source_bytes_inspected"] > 65536):
            raise Refusal("native named work receipt domain")
        # The rejected0..8M bound above is stricter saved-corpus admission;
        # native arbitrary refused calls instead use guarded size_t saturation.
        a.work.charge("source_scalar_inspections")
        if type(native.get("payload_bytes")) is not int or not 0 < native["payload_bytes"] <= 4194304:
            raise Refusal("native whole-live payload receipt")
        reduction = native.get("requested_accuracy_reduction_power", 0)
        a.work.charge("source_scalar_inspections")
        if type(reduction) is not int or not 0 <= reduction <= 32:
            raise Refusal("native requested accuracy reduction")
        if measured_size(shells, work=self.work) > self.payload.input_bytes:
            raise Refusal("decoded typed source exceeds prospective input reserve")
        rows = native.get("rows")
        if not isinstance(rows, list) or len(rows) != len(shells) * len(ells):
            raise Refusal("native row shape")
        if self.nodes is None:
            self.basis = integer_legendre(a)
            self.nodes = prepare_nodes(self.basis, a)
        for ki, shell_row in enumerate(shells):
            for li, ell in enumerate(ells):
                self.context = {"case_id": case["id"], "k_index": ki,
                                "ell_index": li, "ell": ell}
                self.work.stage = "angular-output"
                row = rows[ki * len(ells) + li]
                fixed_keys(row, ("k_index", "ell_index", "ell", "status", "value_bits", "radius_bits"))
                a.work.charge("comparison_cast_operations")
                if (any(type(row[key]) is not int for key in ("k_index", "ell_index", "ell", "status")) or
                    row["k_index"] != ki or row["ell_index"] != li or row["ell"] != ell or row["status"] != 0):
                    raise Refusal("native ordered output/refusal")
                value, _ = binary64(row.get("value_bits"), a)
                error, _ = binary64(row.get("radius_bits"), a)
                if error < ZERO:
                    raise Refusal("negative native radius")
                for bits in (row["value_bits"], row["radius_bits"]):
                    a.work.charge("comparison_cast_operations")
                    pattern = int(bits, 16)
                    if pattern & 0x7fffffffffffffff and pattern & 0x7ff0000000000000 == 0:
                        raise Refusal("native nonzero subnormal storage")
                refs = [assemble(shell_row, ell, q, self.nodes, self.basis, a) for q in (1, 2, 4)]
                self.payload.prepare_report(case_bytes)
                record = compare(value, error, refs, a)
                if ell in (0, 1):
                    scalar_center, scalar_error = scalar_identity(shell_row, ell, a)
                    difference = a.abs_bound(a.sub(a.point(refs[2].center), a.point(scalar_center)))
                    agreement = a.operation("add", radius(refs[2], a), scalar_error, False)
                    limit = a.operation("multiply", Decimal("0.05"), Decimal(record["allocation_lower"]), False)
                    scalar_ok = difference <= agreement and scalar_error <= limit
                    record["scalar_identity"] = {"center": str(scalar_center),
                        "radius": str(scalar_error), "consistency": scalar_ok}
                    if not scalar_ok:
                        record["status"] = "comparison_failed"
                record.update(case_id=case["id"], k_index=ki, ell_index=li, ell=ell)
                self.payload.append_report(records, record, case_bytes)


def measured_size(value, seen=None, work=None):
    if work is not None:
        work.charge("source_graph_inspections")
    if seen is None:
        seen = set()
    identity = id(value)
    if identity in seen:
        return 0
    seen.add(identity)
    size = sys.getsizeof(value)
    if isinstance(value, dict):
        size += sum(measured_size(k, seen, work) + measured_size(v, seen, work) for k, v in value.items())
    elif isinstance(value, (tuple, list)):
        size += sum(measured_size(item, seen, work) for item in value)
    elif isinstance(value, Rational):
        size += measured_size(value.n, seen, work) + measured_size(value.d, seen, work)
    elif isinstance(value, Interval):
        size += measured_size(value.lo, seen, work) + measured_size(value.hi, seen, work)
    return size


def check_runtime(work):
    work.charge("comparison_cast_operations")
    if sys.implementation.name != "cpython" or struct.calcsize("P") != 8:
        raise Refusal("unadmitted Python object profile")
    if sys.int_info.bits_per_digit != 30 or sys.int_info.sizeof_digit != 4:
        raise Refusal("unadmitted integer object profile")
    if (sys.float_info.radix, sys.float_info.mant_dig, sys.float_info.min_exp,
        sys.float_info.max_exp) != (2, 53, -1021, 1024):
        raise Refusal("unadmitted binary64 profile")
    if sys.getsizeof(0) + 4 * ((MAX_BITS + 29) // 30) + 128 > NUMERIC_SLOT:
        raise Refusal("integer slot bound")
    # A16384-bit integer has at most4933 decimal digits. The prospective
    # scratch reservation already includes this numeric slot and its4933-byte
    # constructor string before the runtime layout probe is allocated.
    work.charge("comparison_cast_operations")
    if sys.getsizeof(Decimal("9" * 4933)) > NUMERIC_SLOT:
        raise Refusal("maximum exact Decimal slot bound")
    if any(sys.getsizeof(x) > WRAPPER_SLOT for x in ([], {}, Rational(0), Interval(ZERO, ZERO))):
        raise Refusal("wrapper/header bound")


def check_corpus(corpus, raw_size):
    payload = Payload(raw_size)
    work = Work()
    work.bytes(raw_size)
    records = []
    provenance = []
    report = {"schema": SCHEMA, "status": "refused", "outputs": records,
        "source_provenance": provenance,
        "model": "synthetic/supplied-real-unit-zeta-isotropic-finite-radial-atoms/v1",
        "method": "original-rational-Sturm-GL16-Legendre-angular-q1q2q4/v1",
        "arithmetic": "Decimal100-directed-signed-intervals-owned-Taylor/v1",
        "scope": "conditioned exact emitted finite source; no physical source closure or C_l",
        "python": sys.version, "libmpdec": decimal.__libmpdec_version__,
        "decimal_profile": {"precision": 100, "Emin": -999999,
                            "Emax": 999999, "clamp": 0,
                            "rounding": ["ROUND_FLOOR", "ROUND_CEILING"]}}
    reference = None
    try:
        check_runtime(work)
        if measured_size(corpus, work=work) > payload.input_bytes:
            raise Refusal("parsed graph exceeds prospective input reserve")
        if not isinstance(corpus, dict) or corpus.get("schema") != SCHEMA:
            raise Refusal("corpus schema")
        fixed_keys(corpus, ("schema", "cases"))
        cases = corpus.get("cases")
        if not isinstance(cases, list) or not cases:
            raise Refusal("empty/invalid corpus")
        reference = Reference(work, payload)
        report["arithmetic_profile_controls"] = reference.arithmetic.profile_controls()
        for case in cases:
            reference.case(case, records, provenance)
        report["status"] = "ok" if all(row["status"] == "ok" for row in records) else "comparison_failed"
    except Refusal as error:
        report["reason"] = str(error)
    except (ArithmeticError, MemoryError, ValueError, TypeError) as error:
        report["reason"] = "retained implementation refusal: " + type(error).__name__
    report.update(work=dict(work.counts, total=work.total), rejected_work=work.rejected,
                  source_bytes_inspected=work.source_bytes,
                  rejected_source_bytes=work.rejected_bytes,
                  last_attempt=work.last_attempt, stage=work.stage,
                  failed_reference=(reference.context if reference is not None and report["status"] == "refused" else None),
                  rounding_events=work.rounding_events, peak_reserved_payload_bytes=payload.peak,
                  trapped_events=work.trapped_events,
                  payload_scope="owned input/reference/report plus reserved simultaneous temporaries; excludes interpreter/allocator/RSS")
    return report
