"""Independent reference values for the LTMC quantities of the delayed r-out-of-n policy.

Not used by the Julia code. It produces `test/reference/ltmc_reference.csv`, which
`test/runtests.jl` compares against `SimpleDelayed.policy_quantities`. It computes
the same quantities in two ways that do not share code with the Julia package:

  * closed forms of the paper (main.tex of Paper-LFMO-delayed): Lemma
    `lem:pbbs and exps` and Theorems `thm:LFMO` (exponential delay),
    `thm:determinstic` (deterministic delay) and `thm:LFMOsch` (inspection),
    evaluated with mpmath at 80 digits; and
  * Monte Carlo simulation of the actual network (structure function by graph
    search, LFMO failures, policy with delay), to check the signature
    decomposition and the policy semantics end to end.

Usage (from the repository root; needs python3 with mpmath and numpy):
    python3 test/reference/ltmc_reference.py            # writes the CSV
    python3 test/reference/ltmc_reference.py --mc       # also runs the Monte Carlo check

For every (system, model, beta, delay, r) the CSV gives, for one renewal cycle
under the policy (the paper's notation, Theorem thm:repairdelay):
r = 0..n for exponential and deterministic delays (r = 0: the crew is called at time 0,
with the paper's conventions T_{0:n} := 0 and P(N(T_{0:n}) = l) = 1{l = 0}), r = 1..n otherwise:
    p_fail = P(T_rep = T_fail)          (system failed before the repair started)
    e_nfail = E[N(T_rep)]               (number of failed components at the repair)
    e_trep = E[T_rep]                   (cycle length)
"""
import csv, math, os, random, sys
from fractions import Fraction
from mpmath import mp, mpf, binomial as B, exp, log, sqrt, log1p

mp.dps = 80
ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", ".."))

# ----------------------------------------------------------------------------- models
def psi_fun(model, a):
    a = mpf(a)
    f = {
        "Exp":     lambda x: mpf(x),
        "LFMO-C":  lambda x: mpf(x) / (1 + a * (x - 1)),                      # a = p
        "LFMO-G":  lambda x: log1p(x / a) / log1p(1 / a),                     # a = eta
        "LFMO-IG": lambda x: (sqrt(2 * x + a * a) - a) / (sqrt(2 + a * a) - a),  # a = eta
        "LFMO-S":  lambda x: mpf(x) ** a,                                     # a = alpha
        "BFR":     lambda x: (1 - (1 - a) ** x) / a,                          # a = p
        "LSh":     lambda x: (1 - a) * x + a,                                 # a = theta
    }[model]
    return lambda x: mpf(0) if x == 0 else f(x)

def calibrate(model, beta):
    """Parameter with 2 - Psi(2) = beta (Psi(1) = 1 for every family)."""
    beta = mpf(beta)
    if model == "Exp":
        return mpf(0)
    closed = {"LFMO-C": beta / (2 - beta), "BFR": beta, "LSh": beta,
              "LFMO-S": log(2 - beta) / log(2)}
    if model in closed:
        return closed[model]
    lo, hi = mpf("1e-30"), mpf("1e30")              # G, IG: beta decreasing in eta
    g = lambda eta: 2 - psi_fun(model, eta)(2) - beta
    if g(lo) < 0:
        raise ValueError(f"beta={beta} not attainable for {model}")
    for _ in range(400):
        mid = sqrt(lo * hi)
        if g(mid) > 0: lo = mid
        else: hi = mid
    return sqrt(lo * hi)

# ----------------------------------------------------------------------------- Lemma lem:pbbs and exps
class LFMO:
    def __init__(self, Psi, n):
        self.P = [Psi(x) for x in range(n + 1)]
    def Pw(self, w, f):                       # def:P, probability that f of w fail
        P = self.P
        return sum((-1) ** (l - (w - f) + 1) * B(w, f) * B(f, w - l) * P[l]
                   for l in range(w - f, w + 1)) / P[w]
    def _sumT(self, f, w, g):                  # common shape of eq:ETrn, Kclosed, EminDT, Ddet_Tkn
        return sum(B(w, l) * B(l - 1, w - f) * (-1) ** (l - (w - f + 1)) * g(self.P[l])
                   for l in range(w - f + 1, w + 1))
    def ET(self, f, w):        return mpf(0) if f == 0 else self._sumT(f, w, lambda p: 1 / p)
    def Kmin(self, f, w, t):   return self._sumT(f, w, lambda p: (1 - exp(-t * p)) / p)
    def PDlessT(self, f, w, g): return self._sumT(f, w, lambda p: g / (g + p))
    def PNwt(self, f, w, t):   # eq:PNwt
        return sum(B(w, l) * B(l, w - f) * (-1) ** (l - (w - f)) * exp(-t * self.P[l])
                   for l in range(w - f, w + 1))
    def H(self, f, w, t):      # def:H
        return sum(B(w, l) * B(l, w - f) * (-1) ** (l - (w - f)) / (1 - exp(-t * self.P[l]))
                   for l in range(w - f, w + 1))
    def _PNT(self, f1, f2, w, g):   # eq:PNTrn and its variants; convention T_{0:w} := 0
        if f1 == 0:
            return mpf(1) if f2 == 0 else mpf(0)
        if f1 > f2:
            return mpf(0)
        return sum(self.P[w - gg] * self.Pw(w - gg, f2 - gg) *
                   sum(B(w, l) * B(l, w - gg) * (-1) ** (l - (w - gg)) * g(self.P[l])
                       for l in range(w - gg, w + 1))
                   for gg in range(0, f1))
    def PNT(self, f1, f2, w):       return self._PNT(f1, f2, w, lambda p: 1 / p)
    def PNT_t(self, f1, f2, w, t):  return self._PNT(f1, f2, w, lambda p: (1 - exp(-t * p)) / p)
    def PNT_D(self, f1, f2, w, gm): return self._PNT(f1, f2, w, lambda p: 1 / (gm + p))
    def PNminD(self, f1, f2, w, gm):   # eq:PNcDdetT (f2 = 0 uses T_{0:w} := 0)
        base = self.PNT_D(min(f1, f2), f2, w, gm)
        return base * (gm / (gm + self.P[w - f2]) if f2 < f1 else 1)

# ----------------------------------------------------------------------------- Theorems, per k
def per_k_closed(L, n, r, k, dtype, d):
    """(P(T_k <= T_r + D), E[N(T_rep^k)], E[T_rep^k]) for the k-out-of-n system."""
    m = min(r, k)
    if dtype == "none" or r >= k:
        pf = sum(L.PNT(r, l, n) for l in range(k, n + 1)) if r < k else mpf(1)
        en = sum(j * L.PNT(m, j, n) for j in range(1, n + 1))
        return pf, en, L.ET(m, n)
    PNr = {l: L.PNT(r, l, n) for l in range(r, n + 1)}
    if dtype == "exp":
        gm = 1 / mpf(d)
        psurv = sum(PNr[l] * L.PDlessT(k - l, n - l, gm) for l in range(r, k))       # eqtheo:PTrn
        et = L.ET(r, n) + psurv / gm                                                  # eqtheo:ETrn
        en = mpf(0)
        for j in range(1, n + 1):                                                     # eqtheo:PNTrn
            pj = L.PNT(m, j, n) if j >= k else mpf(0)
            if r <= j:
                pj += sum(PNr[l] * L.PNminD(k - l, j - l, n - l, gm) for l in range(r, min(j, k - 1) + 1))
            en += j * pj
        return 1 - psurv, en, et
    if dtype == "det":
        t = mpf(d)
        psurv = sum(PNr[l] * sum(L.PNwt(i, n - l, t) for i in range(0, k - l)) for l in range(r, k))
        et = L.ET(r, n) + sum(PNr[l] * L.Kmin(k - l, n - l, t) for l in range(r, k))
        en = mpf(0)
        for j in range(1, n + 1):                                                     # eqtheo:PNTrn det
            pj = L.PNT(m, j, n) if j >= k else mpf(0)
            if r <= j < k:
                pj += sum(PNr[l] * L.PNwt(j - l, n - l, t) for l in range(r, j + 1))
            if k <= j:
                pj += sum(PNr[l] * L.PNT_t(k - l, j - l, n - l, t) for l in range(r, k))
            en += j * pj
        return 1 - psurv, en, et
    if dtype == "insp":
        t = mpf(d)
        Hr = {rp: L.H(rp, n, t) for rp in range(0, r)}
        psurv = sum(Hr[rp] * L.PNwt(l - rp, n - rp, t) for rp in range(r) for l in range(r, k))
        et = sum(Hr[rp] * L.Kmin(k - rp, n - rp, t) for rp in range(r))
        en = mpf(0)
        for j in range(r, n + 1):
            if j < k:
                pj = sum(Hr[rp] * L.PNwt(j - rp, n - rp, t) for rp in range(r))
            else:
                pj = sum(Hr[rp] * L.PNT_t(k - rp, j - rp, n - rp, t) for rp in range(r))
            en += j * pj
        return 1 - psurv, en, et
    raise ValueError(dtype)

def policy_closed(s, Psi, n, r, dtype, d):
    L = LFMO(Psi, n)
    pf = en = et = mpf(0)
    for k in range(1, n + 1):
        if s[k - 1] == 0:
            continue
        a, b, c = per_k_closed(L, n, r, k, dtype, d)
        sk = mpf(s[k - 1].numerator) / s[k - 1].denominator
        pf += sk * a; en += sk * b; et += sk * c
    return pf, en, et

# ----------------------------------------------------------------------------- data
def load_signature(name):
    s = []
    with open(os.path.join(ROOT, "data", "signatures", name + ".csv")) as fh:
        rows = [l for l in fh if not l.startswith("#")]
    for row in csv.DictReader(rows):
        if int(row["k"]) >= 1:
            s.append(Fraction(row["signature"]))
    return s

def load_edges(system):
    edges = []
    with open(os.path.join(ROOT, "data", "systems", "yeh2024_edges.csv")) as fh:
        for row in csv.DictReader(fh):
            if int(row["system"]) == system:
                edges.append((row["from"], row["to"]))
    return edges

# ----------------------------------------------------------------------------- Monte Carlo
def works(edges, alive):
    adj = {}
    for e, (u, v) in enumerate(edges):
        if alive[e]:
            adj.setdefault(u, []).append(v); adj.setdefault(v, []).append(u)
    seen, stack = {"Start"}, ["Start"]
    while stack:
        u = stack.pop()
        if u == "End":
            return True
        for v in adj.get(u, []):
            if v not in seen:
                seen.add(v); stack.append(v)
    return False

def monte_carlo(system, Psi, r, dtype, d, cycles, seed):
    edges = load_edges(system); n = len(edges)
    L = LFMO(Psi, n)
    rate = [float(L.P[w]) for w in range(n + 1)]
    batch = {w: [float(L.Pw(w, f)) for f in range(1, w + 1)] for w in range(1, n + 1)}
    rng = random.Random(seed)
    acc = [0.0] * 6       # sum pf, sum N, sum T, sum T^2, sum N^2, sum pf^2 (for SEs)
    for _ in range(cycles):
        alive = [True] * n; nfail = 0; t = 0.0; trep = math.inf; failed_sys = False
        while True:
            w = n - nfail
            dt = rng.expovariate(rate[w]) if w > 0 else math.inf
            if t + dt >= trep:                     # repair starts before the next failure
                t = trep; break
            t += dt
            f = rng.choices(range(1, w + 1), weights=batch[w])[0]
            for e in rng.sample([i for i in range(n) if alive[i]], f):
                alive[e] = False
            nfail += f
            if not works(edges, alive):            # system failure: repair immediately
                failed_sys = True; break
            if nfail >= r and trep == math.inf:    # trigger
                if dtype == "none": break
                if dtype == "exp": trep = t + rng.expovariate(1 / d)
                if dtype == "det": trep = t + d
                if dtype == "insp": trep = math.ceil(t / d) * d
        acc[0] += failed_sys; acc[1] += nfail; acc[2] += t
        acc[3] += t * t; acc[4] += nfail * nfail
    N = cycles
    mean = [acc[0] / N, acc[1] / N, acc[2] / N]
    se = [math.sqrt(mean[0] * (1 - mean[0]) / N),
          math.sqrt(max(acc[4] / N - mean[1] ** 2, 0) / N),
          math.sqrt(max(acc[3] / N - mean[2] ** 2, 0) / N)]
    return mean, se

# ----------------------------------------------------------------------------- main
CASES = []     # (system, model, beta, dtype, delay)
for model, beta in [("Exp", "0"), ("LFMO-C", "0.2"), ("LFMO-G", "0.05"), ("LFMO-IG", "0.5"),
                    ("LFMO-S", "0.2"), ("BFR", "0.5"), ("LSh", "0.2")]:
    for dtype, d in [("none", "0"), ("exp", "0.1"), ("det", "0.3"), ("insp", "0.1")]:
        CASES.append((1, model, beta, dtype, d))
for dtype, d in [("none", "0"), ("exp", "0.03"), ("det", "0.1"), ("insp", "0.3"), ("insp", "1")]:
    CASES.append((3, "BFR", "0.2", dtype, d))
    CASES.append((5, "LFMO-G", "0.2", dtype, d))
    CASES.append((5, "LFMO-S", "0.01", dtype, d))

def main():
    out = os.path.join(ROOT, "test", "reference", "ltmc_reference.csv")
    rows = []
    for system, model, beta, dtype, d in CASES:
        name = f"yeh2024_{system:02d}"
        s = load_signature(name); n = len(s)
        a = calibrate(model, beta); Psi = psi_fun(model, a)
        for r in range(0 if dtype in ("exp", "det") else 1, n + 1):   # r = 0: crew called at time 0
            pf, en, et = policy_closed(s, Psi, n, r, dtype, mpf(d))
            rows.append([system, model, beta, dtype, d, r,
                         mp.nstr(pf, 17), mp.nstr(en, 17), mp.nstr(et, 17)])
        print(f"{name} {model} beta={beta} {dtype} {d}: done", flush=True)
    with open(out, "w") as fh:
        fh.write("# Closed forms of the paper (thm:LFMO, thm:determinstic, thm:LFMOsch), mpmath 80 digits.\n")
        fh.write("# Generated by test/reference/ltmc_reference.py; do not edit by hand.\n")
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["system", "model", "beta", "delay_type", "delay", "r", "p_fail", "e_nfail", "e_trep"])
        w.writerows(rows)
    print("wrote", out)
    if "--mc" in sys.argv:
        print("Monte Carlo check (|z| = |closed form - MC| / SE):")
        for system, model, beta, dtype, d, r in [(1, "LFMO-C", "0.2", "exp", 0.1, 2),
                                                 (1, "BFR", "0.5", "insp", 0.1, 1),
                                                 (2, "LSh", "0.2", "det", 0.3, 2),
                                                 (2, "LFMO-S", "0.2", "none", 0, 3),
                                                 (4, "LFMO-IG", "0.5", "exp", 0.3, 3)]:
            s = load_signature(f"yeh2024_{system:02d}"); n = len(s)
            Psi = psi_fun(model, calibrate(model, beta))
            cf = policy_closed(s, Psi, n, r, dtype, mpf(d))
            mc, se = monte_carlo(system, Psi, r, dtype, d, 200000, seed=12345)
            z = [abs(float(c) - m) / s_ for c, m, s_ in zip(cf, mc, se)]
            print(f"  sys {system} {model} beta={beta} {dtype} {d} r={r}: "
                  f"closed={[round(float(c), 5) for c in cf]} mc={[round(m, 5) for m in mc]} "
                  f"z={[round(x, 2) for x in z]}", flush=True)

if __name__ == "__main__":
    main()
