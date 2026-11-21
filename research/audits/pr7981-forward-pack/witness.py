"""Short independent spending obstruction. Run the full checker for all stages."""
import gzip
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path
import sys


def require(test, message):
    if not test:
        raise ValueError(message)


def witness(path):
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        trace = json.load(stream)
    seal = trace.pop("payload_sha256")
    data = json.dumps(trace, sort_keys=True, separators=(",", ":")).encode()
    require(hashlib.sha256(data).hexdigest() == seal, "trace digest")
    gates, states = trace["gates"], trace["states"]
    n, q, d = len(gates), 2**72, 2**48
    require((n, len(states), trace["P"], trace["p"]) == (23816, n, 72, 48), "geometry")
    require(all(x*x+y*y <= 4*q*q for x,y in trace["entry"]), "entry radius for nonnegative bias")
    require([kind for kind, arg in gates[20:23]] == ["D","H","D"], "three-gate interval")
    def coeffs(j):
        return [tuple(map(F, a["raw"])) for a in trace["tables"][gates[j][1]]]
    # Stage 20, physical slot 20: spectator 4, active coordinate 2.
    require(coeffs(20)[2] == (F(-5,24),0), "stage 20 coefficient")
    require((states[19][20][1], states[20][20][1]) == (36,-7), "imaginary component")
    require(max(x*x+y*y for x,y in states[19]) <= (q//16)**2, "stage 20 radius")
    eps = abs(-7-F(-5,24)*36)/(q*(1+F(1,16)))
    require(gates[21] == ["H",0] and coeffs(22) == [(-3,0),(1,0),(-3,0),(1,0)],
            "H and diagonal coefficients")
    # H0 sends (1,1,0,0) to (2,0,0,0); D sends that to (-6,0,0,0).
    m = 6
    bill = n*m*eps*(2+4*m)
    require(bill > F(1,16*d) > F(1,24*d), "spend obstruction")
    print("PROVED: stages 20-22 force spending above both pack budgets")
    print("epsilon >=", eps, " M >=", m, " spending >=", bill)
    print("ratio to K=3 budget >=", bill*24*d, " ; to two-pack budget >=", bill*16*d)


if __name__ == "__main__":
    witness(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).with_name("pack.json.gz"))
