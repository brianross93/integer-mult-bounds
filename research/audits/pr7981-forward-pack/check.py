"""Offline checker for the left, forward, ternary completed pack of PR 7981.

Independent witness replay and a re-signed semantic control in the spirit of
Brian Ross's Krenn verifier and Theta Twist. Python standard library only.
"""
import argparse
import gzip
import hashlib
import json
import sys
from collections import Counter
from fractions import Fraction as F
from math import isqrt, lcm
from pathlib import Path


BITS = 160  # outward square-root bounds, independent of the working precision


class Mismatch(ValueError):
    pass


def require(test, message):
    if not test:
        raise Mismatch(message)


def seal(trace):
    payload = {k: v for k, v in trace.items() if k != "payload_sha256"}
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def load(path):
    with gzip.open(path, "rt", encoding="utf-8") as f:
        trace = json.load(f)
    require(seal(trace) == trace["payload_sha256"], "payload SHA-256 mismatch")
    return trace


def root(x):
    """Rational lower and upper bounds for the nonnegative square root."""
    x = F(x)
    scale = 1 << BITS
    n = isqrt((x.numerator << (2*BITS)) // x.denominator)
    lo = F(n, scale)
    return lo, lo if lo*lo == x else F(n+1, scale)


ROOT3 = root(3)


def affine(a, b, d):
    """Enclose (a+b*sqrt(3))/d. All arguments are exact."""
    x, y = (F(a+b*v, d) for v in ROOT3)
    return min(x, y), max(x, y)


def ternary_mean(entry):
    """Lean's false character on Fin 3, identity on the fourth active wire.

    Work in integer-word units. The 8 spectator fibres are independent,
    including fibres with a padding coordinate on the spectator axis.
    """
    out = [None]*32
    for s in range(8):
        (a,b), (c,d), (e,f), (g,h) = entry[s::8]
        out[s] = (affine(a+c+e, 0, 3), affine(b+d+f, 0, 3))
        # exp(-2*pi*i/3) = (-1-i*sqrt(3))/2. Its square is its conjugate.
        for k, sign in ((1,1), (2,-1)):
            out[s+8*k] = (affine(2*a-c-e, sign*(d-f), 6),
                           affine(2*b-d-f, sign*(e-c), 6))
        out[s+24] = ((g,g), (h,h))
    return out


def norm_error(values, ideal, divisor):
    """Field dist: max over the 32 rows of Euclidean complex modulus."""
    worst_lo, worst_hi, index = F(), F(), 0
    for k, (z, box) in enumerate(zip(values, ideal, strict=True)):
        lo2, hi2 = F(), F()
        for v, (a,b) in zip(z, box, strict=True):
            lo2 += max(a-v, v-b, 0)**2
            hi2 += max(abs(v-a), abs(v-b))**2
        lo = root(lo2)[0]/divisor
        hi = root(hi2)[1]/divisor
        worst_lo = max(lo, worst_lo)
        if hi > worst_hi:
            worst_hi, index = hi, k
    return worst_lo, worst_hi, index


def mul_words(x, a):
    u,v = x
    c,d = a
    return (u*c-v*d, u*d+v*c)


def rounded(n, q):
    """Nearest integer, with an exact tie assigned to the greater integer."""
    return (2*n+q)//(2*q)


def replay(trace, q):
    """Check each measured state and bound the least relative local epsilon."""
    compiled = []
    for t, table in enumerate(trace["tables"]):
        require(len(table) == 4, f"table width {t}")
        words = []
        for cell in table:
            raw = [F(s) for s in cell["raw"]]
            word = cell["word"]
            require([round(q*z) for z in raw] == word, f"coefficient quantization at table {t}")
            den = lcm(*(z.denominator for z in raw))
            words.append((word, [int(z*den) for z in raw], den))
        compiled.append(words)
    gates, states = trace["gates"], trace["states"]
    require(len(gates) == len(states) > 2, "gate/state count")
    require(gates[0][0] == gates[-1][0] == "D", "missing entrance or exit diagonal")
    half, rem = divmod(len(gates)-2, 2)
    require(rem == 0 and gates[1:1+half] == gates[1+half:-1], "two root passes differ")
    before = trace["entry"]
    worst_lo, worst_hi, worst_gate, worst_lane = F(), F(), -1, -1
    for j, ((op, tag), after) in enumerate(zip(gates, states, strict=True)):
        require(len(after) == 32, f"gate {j} field extent")
        if op == "H":
            require(tag in (0,1), f"gate {j} H significance")
            for i, val in enumerate(after):
                bit = 8*(1 << tag)
                u,v = before[i & ~bit]
                c,d = before[i | bit]
                expected = [u-c, v-d] if i & bit else [u+c, v+d]
                require(val == expected, f"gate {j}, lane {i}: H value")
        elif op == "D":
            num, den2, lane = 0, 1, 0
            radius2 = max(u*u+v*v for u,v in before)
            for i, (x,y) in enumerate(zip(before, after, strict=True)):
                stored, raw_num, den = compiled[tag][i//8]
                re, im = mul_words(x, stored)
                require(y == [rounded(re,q), rounded(im,q)], f"gate {j}, lane {i}: rounded diagonal")
                # Unrounded *rational* gate: coefficient conversion and the
                # final data-word rounding both contribute to Lean's defect.
                re, im = mul_words(x, raw_num)
                residual = (den*y[0]-re)**2 + (den*y[1]-im)**2
                if residual*den2 > num*den**2:
                    num, den2, lane = residual, den**2, i
            if num:
                rlo, rhi = root(F(num,den2))
                alo, ahi = root(radius2)
                lo, hi = rlo/(q+ahi), rhi/(q+alo)
                worst_lo = max(worst_lo, lo)
                if hi > worst_hi:
                    worst_hi, worst_gate, worst_lane = hi, j, lane
        else:
            raise Mismatch(f"unknown gate {op}")
        before = after
    require(before == trace["exit"], "last actual state differs from exit")
    return worst_lo, worst_hi, worst_gate, worst_lane


def show(name, value):
    print(f"{name} = {value} (~{float(value):.12g})")


def check_origin(trace, q):
    require(trace["commit"] == "9585c1fa2ed5f2316ac09e04b357cb720badfcb8", "revision")
    require((trace["pack_id"],trace["diagnostic"],trace["operand"],trace["pack_index"]) ==
            ("3x5-P72/left-forward/pack-1/physical-axis-1-side-3",
             "measure_product.measure(Product,72)","left",1), "diagnostic or pack identity")
    require((trace["sides"], trace["widths"], trace["axis"], trace["direction"],
             trace["opposite"], trace["P"], trace["p"], trace["b"]) ==
            ([5,3], [3,2], 1, -1, False, 72, 48, 8), "pack geometry or direction")
    require(trace["left_negative"] and trace["left_magnitude"] == 741503 and
            trace["right_magnitude"] == 466013 and trace["transform_returned"], "diagnostic invocation")
    require(trace["pack_spans"] == [[0,1],[1,2]], "axis partition")
    require([(b["true"], b["padded"], b["direction"], b["requested_bits"], b["generated_layers"])
             for b in trace["compiler"]] == [(3,4,-1,80,189),(5,8,-1,80,581)],
            "private axis compilation")
    rootgen = trace["root_generation"]
    require([F(x) for x in rootgen["turns"]] == [-F(j*(j+3),6) for j in range(3)]
            and list(map(F,rootgen["radii"])) == [1,1,1] and F(rootgen["rotation"]) == -F(1,12)
            and F(rootgen["accretivity"]) == F(1,12), "ternary root phases")
    require(tuple(rootgen[k] for k in ("requested_bits","column_bits","phase_bits",
                                       "phase_grid","phase_internal_bits")) == (80,112,116,6,296),
            "private coefficient precisions")
    column = [[F(v) for v in cell] for cell in rootgen["column"]]
    require(len(column) == 3 and column[1] == column[2] and column[0] != column[1] and
            all(len(cell)==2 and all(((1<<112)*v).denominator==1 for v in cell)
                for cell in column), "private circulant column")
    # Packed source: five true x positions per valid y row. Lean's second
    # axis input is 5^-1 k mod 3. The first coordinate has empty prefix.
    mapping = trace["source_slot_at_crt_slot"]
    require(sorted(mapping) == list(range(32)), "CRT tag probe is not a permutation")
    before = [[0,0] for _ in range(32)]
    for k in range(15):
        src = (k//5)*8 + k%5
        dst = k%5 + 8*((2*k)%3)
        require(mapping[dst] == src, f"CRT input address at k={k}")
        digit = (741503 >> (8*k)) % 256
        before[src] = [q*digit//1024,0]
    require(before == trace["before_crt"], "pre-CRT left operand")
    require([before[j] for j in mapping] == trace["operand_field"], "post-CRT left operand")
    require(len(trace["entry"]) == len(trace["exit"]) == 32, "endpoint extent")
    h_count = sum(op == "H" for op, _ in trace["gates"])
    require(0 < h_count == trace["h_after_exit"]-trace["h_before_entry"], "physical H call count")
    require((trace["h_before_entry"],trace["h_after_exit"],trace["full_one_pass_layer_count"],
             trace["storage_bytes_per_component"]) == (41832,52416,6615,23)
            and trace["crt_stats"] == {"single_targets":24}, "execution path")
    require([(k["physical_positions"],k["selected_bits"],k["L"],k["bytes_per_record"],k["calls"])
             for k in trace["h_kernel_shapes"]] == [([3],[0],2,46,72576),([4],[1],2,46,12096)],
            "stopped H kernel path")


def audit(trace, endpoint_only):
    from interval_norm import UNIT, build_certificate, sixfold_witness
    q = 1 << trace["P"]
    check_origin(trace, q)
    delta = F(1,8*(1 << 48)*3)  # GatePacks.at_selected_precision uses K=d+1.
    entry = trace["entry"]
    radius = root(max(a*a+b*b for a,b in entry))[1]/q
    require(radius <= 2, "gate_defects entry radius exceeds two")
    print(trace["pack_id"], "stages", dict(Counter(op for op,_ in trace["gates"])))
    e_lo, e_hi, j, k = replay(trace,q)
    print("all local transitions PASS; largest relative rounding bound at stage", j, "slot", k)
    show("entry radius upper",radius)
    show("relative epsilon lower",e_lo)
    show("relative epsilon upper",e_hi)
    ideal = ternary_mean(entry)
    ylo, yhi, where = norm_error(trace["exit"],ideal,q)
    show("delta",delta)
    print("endpoint maximum at physical slot", where)
    show("endpoint residual lower",ylo)
    show("endpoint residual upper",yhi)
    show("endpoint margin lower",delta-yhi)
    require(ylo <= delta, f"padded ternary mean: physical slot {where}, "
            f"residual >= {float(ylo):.12g} > delta {float(delta):.12g}")
    if yhi > delta:
        print("INCONCLUSIVE: the endpoint enclosure straddles the allowance")
        return 3
    print("actual-to-ideal completed-pack Close PASS")
    if where == 31 and entry[31] == [0,0] and trace["exit"][31] == [60,-30]:
        other_hi = norm_error(trace["exit"][:31],ideal[:31],1)[1]
        if other_hi < 56:
            print("exact worst endpoint residual = 30*sqrt(5)/2^72 at slot 31;"
                  " all other residuals < 56 working ulps")
    if endpoint_only:
        return 0
    n = len(trace["gates"])
    # Minimal obstruction uses only stages 21, 22, 26 and the total stage count.
    kind, t = trace["gates"][26]
    require(kind == "D" and list(map(F,trace["tables"][t][2]["raw"])) == [F(1,8),0],
            "stage 26, slot 17 is not an exact eighth")
    before = trace["states"][25]
    x, y = before[17][0], trace["states"][26][17][0]
    require(abs(F(x,8)-y) == F(1,2), "stage 26 lacks half-unit defect")
    radius2 = max(a*a+b*b for a,b in before)
    require(4*radius2 <= q*q, "rounding witness exceeded radius 1/2")
    print("stage 26, slot 17 real: 1/8 times",x,"rounds to",y)
    print("input sup norm squared =", F(radius2,q*q), "<= 1/4")
    m = sixfold_witness(trace["tables"],trace["gates"])
    show("two-stage M lower",m)
    simple = F(n*m*(2+4*m),3*q)
    show("elementary spend lower",simple)
    show("elementary excess over two-pack budget",simple-F(1,16*(1 << 48)))
    cert = build_certificate(trace["tables"],trace["gates"],trace["P"],raw=True)
    show("all-subword M upper",cert["M_upper"])
    show(f"necessary M lower on [{cert['lower_span'][0]},{cert['lower_span'][1]})",
         cert["M_lower"])
    # Full matrix is quantized only in our enclosure. Apply its centre to
    # the exact tape. The operator error bound multiplies the entry radius.
    centre = [[0,0] for _ in range(32)]
    for s in range(8):
        for k, row in enumerate(cert["full"]):
            p = [mul_words(x,a) for x,a in zip(entry[s::8],row,strict=True)]
            centre[s+8*k] = [F(sum(v[i] for v in p),UNIT) for i in range(2)]
    bias_lo, bias_hi, _ = norm_error(centre,ideal,q)
    matrix_error = cert["error"]*radius
    beta_hi = bias_hi + matrix_error
    beta_lo = max(0,bias_lo-matrix_error)
    show("synthesis beta upper",beta_hi)
    show("matrix enclosure contribution",matrix_error)
    # Substitute this pack's observed geometry in the submitted selector.
    # These caps are about GateData for this one trace, not universal inspect.
    axes = len(trace["sides"])
    public_bits = trace["p"]
    extent = sum(trace["widths"])
    max_width = (extent + axes - 1) // axes
    active = max(1, axes // (1024 * (public_bits + 1).bit_length()))
    E = 4096 * active * (max_width + 6)
    J = 20000 * (axes + 1) * (max_width + 2)**6 + 2 * axes
    K = axes + 1
    required_bits = public_bits + 3 + 8*E + 4*(J*K).bit_length() + K.bit_length() + 64
    A = 1 << E
    require(n <= J and cert["M_upper"] <= A and e_hi*q <= K*K*A*A
            and beta_hi*q <= 2*K*A, "selected-precision numeric cap on this trace")
    print("individual inspect numeric caps PASS; E, J, K =", E, J, K)
    print("selector requires", required_bits, "working bits; observed", trace["P"])
    low = n*cert["M_lower"]*e_lo
    high = n*cert["M_upper"]*e_hi
    spend_lo = low*(2+4*cert["M_lower"])+beta_lo
    spend_hi = high*(2+4*cert["M_upper"])+beta_hi
    show("bootstrap LHS upper",high)
    show("spend LHS lower",spend_lo)
    show("spend LHS upper",spend_hi)
    if spend_lo > delta:
        print("spend FAIL for the captured D/H decomposition")
        show("certified spend excess", spend_lo-delta)
        return 2
    elif high <= F(1,2) and spend_hi <= delta:
        print("GateWitness PASS for the captured D/H decomposition")
        return 0
    else:
        print("GateWitness INCONCLUSIVE: certified upper and lower bounds do not decide spend")
        return 3


def negative(trace, out):
    """Double the last true-wire gain, recompute its rounded output, and reseal."""
    q = 1 << trace["P"]
    op, t = trace["gates"][-1]
    require(op == "D", "final operation is not a diagonal")
    table = []
    for i, cell in enumerate(trace["tables"][t]):
        raw = [F(v)*(2 if i < 3 else 1) for v in cell["raw"]]
        table.append(dict(raw=list(map(str, raw)), word=[round(q*v) for v in raw]))
    trace["gates"][-1] = ["D",len(trace["tables"])]
    trace["tables"].append(table)
    before = trace["states"][-2]
    new = [[rounded(v,q) for v in mul_words(x, table[i//8]["word"])] for i,x in enumerate(before)]
    require(new != trace["exit"], "negative control is inert")
    trace["exit"] = trace["states"][-1] = new
    trace["payload_sha256"] = seal(trace)
    with open(out, "wb") as f:
        f.write(gzip.compress(json.dumps(trace,separators=(",",":")).encode(),mtime=0))
    print("generated semantic negative", trace["payload_sha256"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("trace", type=Path)
    parser.add_argument("--negative", type=Path, metavar="OUTPUT")
    parser.add_argument("--endpoint-only", action="store_true",
                        help="verify all transitions and the local Close inequality; skip witness funding")
    args = parser.parse_args()
    trace = load(args.trace)
    if args.negative:
        negative(trace, args.negative)
        return
    print("integrity PASS", trace["payload_sha256"])
    return audit(trace, args.endpoint_only)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (Mismatch, KeyError, IndexError, ZeroDivisionError) as exc:
        print("REJECT", str(exc), file=sys.stderr)
        sys.exit(1)
