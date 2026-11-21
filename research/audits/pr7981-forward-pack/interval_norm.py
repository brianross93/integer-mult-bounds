"""Certified complex sup-norm bounds for this diagnostic's four-lane gate word.

Matrices act on column vectors. The row triangle inequality bounds their
complex sup operator norm by the largest row sum of moduli. Choosing unit
input phases conjugate to that row's coefficient phases attains the bound.
A ball stores a Gaussian dyadic centre and an error in this operator norm.
Spectator fibres are identical blocks. No production code or numerical library.
"""
from fractions import Fraction
from functools import lru_cache
from math import isqrt
from typing import NamedTuple


BITS = 160
UNIT = 1 << BITS
LANES = 4


def ceildiv(x, y):
    return -((-x) // y)


class Ball(NamedTuple):
    rows: tuple
    error: int
    lo: int
    hi: int


def enclose(rows, error):
    lo = hi = 0
    for row in rows:
        lower = upper = 0
        for x, y in row:
            q = x*x + y*y
            n = isqrt(q)
            lower += n
            upper += n + (n*n != q)
        lo, hi = max(lo, lower), max(hi, upper)
    return Ball(rows, error, max(0, lo-error), hi+error)


def diagonal(zs):
    # Rounding one diagonal entry to the centre costs < 2 last-place units
    # in complex modulus. Only one entry is supported per row.
    nums = [(r * UNIT, i * UNIT) for r, i in zs]
    residual = max((x.denominator != 1) + (y.denominator != 1) for x, y in nums)
    zs = [(x.numerator // x.denominator, y.numerator // y.denominator) for x, y in nums]
    return enclose(tuple(tuple(z if i == j else (0, 0) for j in range(LANES))
                         for i, z in enumerate(zs)), residual)


def hadamard(bit):
    if bit not in (0, 1):
        raise ValueError("invalid local significance")
    mask = 1 << bit
    def coefficient(i, j):
        if i == j:
            return -UNIT if i & mask else UNIT
        return UNIT if j == i ^ mask else 0
    return enclose(tuple(tuple((coefficient(i, j), 0) for j in range(LANES))
                         for i in range(LANES)), 0)


IDENTITY = diagonal([(Fraction(1), Fraction(0))]*LANES)


@lru_cache(maxsize=4096)
def multiply(a, b):
    """Enclose exact AB; a.hi already includes a.error."""
    if a is IDENTITY:
        return b
    if b is IDENTITY:
        return a
    rows = []
    lost = 0
    for ar in a.rows:
        row = []
        cost = 0
        for j in range(LANES):
            x = y = 0
            for k, (u, v) in enumerate(ar):
                r, s = b.rows[k][j]
                x += u*r - v*s
                y += u*s + v*r
            xr, xi = divmod(x, UNIT)
            yr, yi = divmod(y, UNIT)
            row.append((xr, yr))
            cost += (xi != 0) + (yi != 0)
        rows.append(tuple(row))
        lost = max(lost, cost)
    # AB - mid(A)mid(B) = A*(B-mid(B)) + (A-mid(A))*mid(B).
    inherited = ceildiv(a.hi*b.error + a.error*(b.hi-b.error), UNIT)
    return enclose(tuple(rows), inherited+lost)


def compile_gates(tables, gates, p, raw):
    one = Fraction(1, 1 << p)
    @lru_cache(maxsize=None)
    def pair(kind, index):
        if kind == "H":
            h = hadamard(index)
            inv = enclose(tuple(tuple((x//2, 0) for x, _ in row) for row in h.rows), 0)
            return h, inv
        if kind != "D":
            raise ValueError("unknown gate kind")
        table = tables[index]
        zs = [(Fraction(z["raw"][0]), Fraction(z["raw"][1])) if raw else
              (one*z["word"][0], one*z["word"][1]) for z in table]
        if len(zs) != LANES:
            raise ValueError("expected four active lanes")
        inverse = [(x/(x*x+y*y), -y/(x*x+y*y)) for x, y in zs]
        return diagonal(zs), diagonal(inverse)
    return [pair(k, s) for k, s in gates]


def sixfold_witness(tables, gates):
    """Two actual stages, one spectator fibre, one output coefficient.

    The unit sup-norm input (1,1,0,0) on this fibre goes to (2,0,0,0).
    The next diagonal has first coefficient -3, so its first output is -6.
    Taking the second input of the interval premise to be zero forces M >= 6.
    """
    kind, tag = gates[22]
    if tuple(gates[21]) != ("H", 0) or kind != "D":
        raise ValueError("two-gate witness does not match the gate sequence")
    first = tuple(Fraction(x) for x in tables[tag][0]["raw"])
    if first != (-3, 0):
        raise ValueError("two-gate witness coefficient differs from -3")
    return 6


def build_certificate(tables, gates, p, raw=False):
    """Return exact rational interval bounds and a certified full-product ball.

    A balanced tree avoids sequential interval dependency growth. All prefixes
    P[j] and exact inverse prefixes are bounded. Every subword is P[j] P[i]^-1.
    For a lower bound, any particular subword supplies its induced row norm.
    """
    pairs = compile_gates(tables, gates, p, raw)
    tree = {}
    best = (UNIT, 0, 0)
    def record(ball, i, j):
        nonlocal best
        if ball.lo > best[0]:
            best = (ball.lo, i, j)
    def build(i, j):
        if j == i+1:
            fwd, inv = pairs[i]
        else:
            k = (i+j)//2
            lf, li = build(i,k)
            rf, ri = build(k,j)
            fwd, inv = multiply(rf, lf), multiply(li, ri)
        tree[i,j] = fwd, inv
        record(fwd, i, j)
        return fwd, inv
    if not gates:
        raise ValueError("empty gate word")
    full, _ = build(0, len(gates))
    max_pre = max_inv = UNIT
    # Inclusion of every leaf start and the end covers all N+1 prefixes.
    def visit(i, j, forward, inverse):
        nonlocal max_pre, max_inv
        if j == i+1:
            max_pre = max(max_pre, forward.hi)
            max_inv = max(max_inv, inverse.hi)
            record(forward, 0, i)
            return
        k = (i+j)//2
        visit(i, k, forward, inverse)
        lf, li = tree[i,k]
        visit(k, j, multiply(lf,forward), multiply(inverse,li))
    visit(0,len(gates),IDENTITY,IDENTITY)
    max_pre = max(max_pre, full.hi)
    max_inv = max(max_inv, tree[0,len(gates)][1].hi)
    return {"full": full.rows, "error": Fraction(full.error, UNIT),
            "M_upper": Fraction(max_pre*max_inv, UNIT**2),
            "M_lower": Fraction(best[0], UNIT), "lower_span": best[1:]}
