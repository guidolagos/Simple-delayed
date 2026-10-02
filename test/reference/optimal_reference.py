"""Independent reference values for the optimal repair policies (Section sec:optimal).

Not used by the Julia code. It produces `test/reference/optimal_reference.csv`,
which `test/runtests.jl` compares with `SimpleDelayed.optimal_count` and
`SimpleDelayed.optimal_state`. It shares no code with the Julia package and
follows the paper as literally as possible, with mpmath arithmetic:

  * best r-out-of-n policy, r in {0,...,n} (r = 0 only with continuous inspection
    and a delay): closed forms thm:LFMO / thm:determinstic / thm:LFMOsch
    (via ltmc_reference.py);
  * count-based optimum: Dinkelbach (alg:dinkelbach) with the passes of the paper:
      exponential delay   eq:ccall exp - eq:lcall exp (recursion, NOT the signature form)
      deterministic delay prop:after call det, with the closed forms of lem:pbbs and exps
      fixed inspections   lem:period and eq:bellman fixed
      no delay            calling = repairing at once (the case D = 0 of the det formulas)
    and, for small n, brute force over all stopping sets;
  * state-based optimum, for small n: the same Bellman equations on the full state
    space, with the sums over sub-states enumerated directly (3^n pairs) and the
    per-state signatures counted by brute force; for the exponential delay the
    after-call quantities use the state-space version of eq:ccall exp;
  * Monte Carlo simulation of the actual networks under the computed optimal
    policies, to check the semantics end to end (python3 optimal_reference.py --mc).

Usage (from the repository root; needs mpmath):
    python3 test/reference/optimal_reference.py          # writes the CSV
    python3 test/reference/optimal_reference.py --mc     # also Monte Carlo checks
"""
import csv, itertools, math, os, random, sys
from fractions import Fraction
from mpmath import mp, mpf, binomial as B, exp

sys.path.insert(0, os.path.dirname(__file__))
import ltmc_reference as R

mp.dps = 40
ROOT = R.ROOT

# ----------------------------------------------------------------------------- data
def load_W(name):
    with open(os.path.join(ROOT, "data", "signatures", name + ".csv")) as fh:
        rows = [l for l in fh if not l.startswith("#")]
    return [int(row["W"]) for row in csv.DictReader(rows)]

class System:
    def __init__(self, sid):
        self.id = sid
        self.name = f"yeh2024_{sid:02d}"
        self.edges = R.load_edges(sid)
        self.n = len(self.edges)
        self.W = load_W(self.name)
        self.Sbar = [mpf(self.W[f]) / B(self.n, f) for f in range(self.n + 1)]
        self.s = R.load_signature(self.name)
        self._phi = {}
    def phi(self, x):                        # bit e of x = 1 iff component e+1 works
        if x not in self._phi:
            alive = [(x >> e) & 1 == 1 for e in range(self.n)]
            self._phi[x] = R.works(self.edges, alive)
        return self._phi[x]

class Costs:
    def __init__(self, n, rho, sigma, csys=1):
        self.n, self.rho, self.sigma = n, mpf(rho), mpf(sigma)
        self.csys = mpf(csys); self.cset = self.sigma * self.csys; self.ccmp = self.rho * self.csys / n
    def rep(self, f):  return self.cset + self.ccmp * f          # c_cmp(f)
    def fail(self, f): return self.rep(f) + self.csys             # c^fail(f)

def is_continuous(dtype): return dtype in ("none", "exp", "det")

# ----------------------------------------------------------------------------- r-out-of-n
def admissible_r(n, dtype):
    return range(0, n + 1) if dtype in ("exp", "det") else range(1, n + 1)

def best_r(sysm, Psi, dtype, d, costs):
    out = {}
    for r in admissible_r(sysm.n, dtype):
        pf, en, et = R.policy_closed(sysm.s, Psi, sysm.n, r, dtype, mpf(d))
        out[r] = (costs.cset + costs.ccmp * en + costs.csys * pf) / et
    r = min(out, key=lambda k: out[k])
    return r, out[r], out

# ----------------------------------------------------------------------------- after-call quantities
class Tables:
    """Closed forms of the paper for a residual LFMO system with w components."""
    def __init__(self, Psi, n, dtype, d):
        self.n, self.dtype = n, dtype
        self.L = R.LFMO(Psi, n)
        self.P = self.L.P
        self.d = mpf(d)
        self.Pw = {(w, j): self.L.Pw(w, j) for w in range(1, n + 1) for j in range(1, w + 1)}
        if dtype in ("det", "insp"):
            t = self.d
            self.pN = {(w, j): self.L.PNwt(j, w, t) for w in range(1, n + 1) for j in range(0, w + 1)}
            self.K = {(w, i): self.L.Kmin(i, w, t) for w in range(1, n + 1) for i in range(1, w + 1)}
            self.PT = {(w, i, j): self.L.PNT_t(i, j, w, t)
                       for w in range(1, n + 1) for i in range(1, w + 1) for j in range(i, w + 1)}

    def a_b(self, w, i, f, costs):
        """prop:after call det: (a^{(w)}_i, b^{(w)}_i) for a call at f = n - w failed (det delay)."""
        if self.dtype == "none":
            return costs.rep(f), mpf(0)
        t = self.d
        pj = [self.pN[(w, j)] if j < i else self.PT[(w, i, j)] for j in range(0, w + 1)]
        pfail = 1 - sum(self.pN[(w, j)] for j in range(0, i))
        a = sum(costs.rep(f + j) * pj[j] for j in range(0, w + 1)) + costs.csys * pfail
        return a, self.K[(w, i)]

def cond_sig(sysm, f):
    w = sysm.n - f
    return [(sysm.Sbar[f + i - 1] - sysm.Sbar[f + i]) / sysm.Sbar[f] for i in range(1, w + 1)]

# ----------------------------------------------------------------------------- count-based
def count_pass(sysm, T, costs, g, forced=None):
    """One pass of alg:count exp / its det, none and insp variants. forced(f) -> bool overrides."""
    n, Sb, P = sysm.n, sysm.Sbar, T.P
    c, l, act = {}, {}, {}
    # exponential delay: eq:ccall exp - eq:lcall exp, computed for all f first
    if T.dtype == "exp":
        gam = 1 / T.d
        cc, lc = {}, {}
        for f in range(n - 1, -1, -1):
            if Sb[f] == 0: continue
            w = n - f
            sc = sl = mpf(0)
            for j in range(1, w + 1):
                ps = Sb[f + j] / Sb[f]
                sc += T.Pw[(w, j)] * ((ps * cc[f + j] if ps > 0 else 0) + (1 - ps) * costs.fail(f + j))
                sl += T.Pw[(w, j)] * (ps * lc[f + j] if ps > 0 else 0)
            cc[f] = (gam * costs.rep(f) + P[w] * sc) / (gam + P[w])
            lc[f] = (1 + P[w] * sl) / (gam + P[w])
    for f in range(n - 1, -1, -1):
        if Sb[f] == 0: continue
        w = n - f
        if T.dtype == "insp":
            sig = cond_sig(sysm, f)
            q = [T.pN[(w, j)] * Sb[f + j] / Sb[f] for j in range(0, w)]
            cD = sum(sig[i - 1] * sum(costs.fail(f + j) * T.PT[(w, i, j)] for j in range(i, w + 1))
                     for i in range(1, w + 1))
            lD = sum(sig[i - 1] * T.K[(w, i)] for i in range(1, w + 1))
            cw = (cD + sum(q[j] * c[f + j] for j in range(1, w) if q[j] > 0)) / (1 - q[0])
            lw = (lD + sum(q[j] * l[f + j] for j in range(1, w) if q[j] > 0)) / (1 - q[0])
            if f >= 1:
                call = (forced(f) if forced else costs.rep(f) <= cw - g * lw)
            else:
                call = False
            act[f] = call
            c[f], l[f] = (costs.rep(f), mpf(0)) if call else (cw, lw)
        else:
            if T.dtype == "exp":
                ccall, lcall = cc[f], lc[f]
            else:
                sig = cond_sig(sysm, f)
                ab = [T.a_b(w, i, f, costs) for i in range(1, w + 1)]
                ccall = sum(si * a for si, (a, _) in zip(sig, ab))
                lcall = sum(si * b for si, (_, b) in zip(sig, ab))
            ccont = lcont = mpf(0)
            for j in range(1, w + 1):
                ps = Sb[f + j] / Sb[f]
                ccont += T.Pw[(w, j)] * ((ps * c[f + j] if ps > 0 else 0) + (1 - ps) * costs.fail(f + j))
                lcont += T.Pw[(w, j)] * (ps * l[f + j] if ps > 0 else 0)
            lcont += 1 / P[w]
            allowed = not (T.dtype == "none" and f == 0)
            if forced:
                call = allowed and forced(f)
            else:
                call = allowed and (ccall - g * lcall <= ccont - g * lcont)
            act[f] = call
            c[f], l[f] = (ccall, lcall) if call else (ccont, lcont)
    return act, c[0], l[0]

def dinkelbach(passf, g0, tol=mpf(10) ** -25, maxit=100):
    g = g0
    for it in range(1, maxit + 1):
        act, c, l = passf(g)
        if c - g * l >= -tol * max(1, abs(c)):
            return c / l, act, it
        g = c / l
    raise RuntimeError("Dinkelbach did not converge")

def policy_string(sysm, act):
    return "".join("-" if sysm.Sbar[f] == 0 else ("C" if act.get(f) else ".") for f in range(sysm.n))

def brute_force_count(sysm, T, costs):
    """min over all stopping sets of g = c/l (small n)."""
    fs = [f for f in range(sysm.n) if sysm.Sbar[f] > 0 and not (f == 0 and T.dtype in ("none", "insp"))]
    best = None
    for k in range(len(fs) + 1):
        for S in itertools.combinations(fs, k):
            S = set(S)
            _, c, l = count_pass(sysm, T, costs, 0, forced=lambda f: f in S)
            if l > 0 and (best is None or c / l < best):
                best = c / l
    return best

# ----------------------------------------------------------------------------- state-based
class StateSpace:
    """Full state space of a small system; sub-state sums enumerated directly."""
    def __init__(self, sysm):
        self.sysm = sysm
        n = sysm.n
        self.order = sorted(range(1 << n), key=lambda x: bin(x).count("1"))
        self.subs = {}
        for x in range(1 << n):
            lst, y = [], x
            while True:                      # all submasks of x
                y = (y - 1) & x
                lst.append(y)
                if y == 0: break
            self.subs[x] = [y for y in lst if y != x]
        self.Wx = {}
        for x in range(1 << n):
            w = bin(x).count("1")
            cnt = [0] * (w + 1)
            for y in self.subs[x] + [x]:
                if sysm.phi(y): cnt[w - bin(y).count("1")] += 1
            self.Wx[x] = cnt
    def sig(self, x):
        w = bin(x).count("1"); W = self.Wx[x]
        Sb = [mpf(W[i]) / B(w, i) for i in range(w + 1)]
        return [Sb[i - 1] - Sb[i] for i in range(1, w + 1)]

def state_pass(S, T, costs, g, forced=None):
    sysm, n, P = S.sysm, S.sysm.n, T.P
    phi = sysm.phi
    hc, hl, act = {}, {}, {}
    cc, lc = {}, {}
    gam = 1 / T.d if T.dtype == "exp" else None
    for x in S.order:
        w = bin(x).count("1"); f = n - w
        if not phi(x):
            hc[x] = costs.fail(f) if T.dtype != "insp" else mpf(0)
            hl[x] = mpf(0)
            continue
        if T.dtype == "insp":
            sig = S.sig(x)
            p = [T.pN[(w, j)] / B(w, j) for j in range(0, w + 1)]
            cD = sum(sig[i - 1] * sum(costs.fail(f + j) * T.PT[(w, i, j)] for j in range(i, w + 1))
                     for i in range(1, w + 1))
            lD = sum(sig[i - 1] * T.K[(w, i)] for i in range(1, w + 1))
            sc = sum(p[w - bin(y).count("1")] * hc[y] for y in S.subs[x])
            sl = sum(p[w - bin(y).count("1")] * hl[y] for y in S.subs[x])
            cw, lw = (cD + sc) / (1 - p[0]), (lD + sl) / (1 - p[0])
            call = f >= 1 and (forced(x, f) if forced else costs.rep(f) <= cw - g * lw)
            act[x] = call
            hc[x], hl[x] = (costs.rep(f), mpf(0)) if call else (cw, lw)
            continue
        # continuous inspection
        if T.dtype == "exp":                 # eq:ccall exp on the state space
            sc = sum(T.Pw[(w, w - bin(y).count("1"))] / B(w, w - bin(y).count("1")) *
                     (cc[y] if phi(y) else costs.fail(n - bin(y).count("1"))) for y in S.subs[x])
            sl = sum(T.Pw[(w, w - bin(y).count("1"))] / B(w, w - bin(y).count("1")) *
                     (lc[y] if phi(y) else 0) for y in S.subs[x])
            cc[x] = (gam * costs.rep(f) + P[w] * sc) / (gam + P[w])
            lc[x] = (1 + P[w] * sl) / (gam + P[w])
            ccall, lcall = cc[x], lc[x]
        else:
            sig = S.sig(x)
            ab = [T.a_b(w, i, f, costs) for i in range(1, w + 1)]
            ccall = sum(si * a for si, (a, _) in zip(sig, ab))
            lcall = sum(si * b for si, (_, b) in zip(sig, ab))
        ccont = sum(T.Pw[(w, w - bin(y).count("1"))] / B(w, w - bin(y).count("1")) * hc[y] for y in S.subs[x])
        lcont = 1 / P[w] + sum(T.Pw[(w, w - bin(y).count("1"))] / B(w, w - bin(y).count("1")) * hl[y]
                               for y in S.subs[x])
        allowed = not (T.dtype == "none" and f == 0)
        call = allowed and (forced(x, f) if forced else ccall - g * lcall <= ccont - g * lcont)
        act[x] = call
        hc[x], hl[x] = (ccall, lcall) if call else (ccont, lcont)
    full = (1 << n) - 1
    return act, hc[full], hl[full]

# ----------------------------------------------------------------------------- Monte Carlo
def simulate(sysm, Psi, dtype, d, costs, decide, cycles, seed):
    """Long-run cost rate of a policy on the actual network. decide(x, f) -> call/repair."""
    n = sysm.n
    L = R.LFMO(Psi, n)
    rate = [float(L.P[w]) for w in range(n + 1)]
    batch = {w: [float(L.Pw(w, j)) for j in range(1, w + 1)] for w in range(1, n + 1)}
    rng = random.Random(seed)
    tot_c = tot_t = 0.0
    cs, ts = [], []
    d = float(d)
    rep = lambda f: float(costs.rep(f)); csys = float(costs.csys)
    for _ in range(cycles):
        x = (1 << n) - 1; t = 0.0; arrive = math.inf
        if dtype in ("exp", "det") and decide(x, 0):
            arrive = rng.expovariate(1 / d) if dtype == "exp" else d
        next_insp = d if dtype == "insp" else math.inf
        while True:
            w = bin(x).count("1")
            dt = rng.expovariate(rate[w]) if w > 0 else math.inf
            if dtype == "insp":
                if t + dt >= next_insp:      # inspection first
                    t = next_insp; f = n - w
                    if f >= 1 and decide(x, f):
                        c = rep(f); break
                    next_insp += d
                    continue
            elif t + dt >= arrive:           # crew arrives first
                t = arrive; c = rep(n - w); break
            t += dt
            j = rng.choices(range(1, w + 1), weights=batch[w])[0]
            alive = [e for e in range(n) if (x >> e) & 1]
            for e in rng.sample(alive, j):
                x &= ~(1 << e)
            f = n - bin(x).count("1")
            if not sysm.phi(x):
                c = rep(f) + csys; break
            if dtype in ("none", "exp", "det") and arrive == math.inf and decide(x, f):
                if dtype == "none":
                    c = rep(f); break
                arrive = t + (rng.expovariate(1 / d) if dtype == "exp" else d)
        cs.append(c); ts.append(t)
    C, Tt = sum(cs) / cycles, sum(ts) / cycles
    g = C / Tt
    var = sum((ci - g * ti) ** 2 for ci, ti in zip(cs, ts)) / (cycles - 1)
    return g, math.sqrt(var / cycles) / Tt

# ----------------------------------------------------------------------------- cases
def run_case(sid, model, beta, dtype, d, rho, sigma, state=True, brute=False):
    sysm = System(sid)
    Psi = R.psi_fun(model, R.calibrate(model, beta))
    costs = Costs(sysm.n, rho, sigma)
    T = Tables(Psi, sysm.n, dtype, d)
    r0, g_r, _ = best_r(sysm, Psi, dtype, d, costs)
    g_c, act_c, it_c = dinkelbach(lambda g: count_pass(sysm, T, costs, g), g_r)
    res = dict(system=sid, n=sysm.n, model=model, beta=beta, delay_type=dtype, delay=d,
               ratio=rho, setup=sigma, best_r=r0, ltmc_best_r=g_r, ltmc_count=g_c,
               count_policy=policy_string(sysm, act_c), ltmc_state=None)
    if brute:
        res["ltmc_count_brute"] = brute_force_count(sysm, T, costs)
    if state:
        S = StateSpace(sysm)
        g_s, act_s, it_s = dinkelbach(lambda g: state_pass(S, T, costs, g), g_c)
        res["ltmc_state"] = g_s
        res["_state_act"] = act_s
    res["_count_act"] = act_c
    res["_objs"] = (sysm, Psi, T, costs)
    return res

CASES = []
for model, beta in [("Exp", "0"), ("LFMO-C", "0.2"), ("BFR", "0.5"), ("LSh", "0.2"), ("LFMO-S", "0.5")]:
    for dtype, d in [("none", "0"), ("exp", "0.1"), ("det", "0.3"), ("insp", "0.3")]:
        for rho, sigma in [("1", "0.1"), ("0.5", "0.05")]:
            CASES.append((1, model, beta, dtype, d, rho, sigma))
for model, beta in [("LFMO-G", "0.2"), ("BFR", "0.2"), ("Exp", "0")]:
    for dtype, d in [("none", "0"), ("exp", "0.3"), ("det", "0.1"), ("insp", "0.1")]:
        CASES.append((3, model, beta, dtype, d, "2", "0.2"))
        CASES.append((2, model, beta, dtype, d, "0.2", "0.05"))

def main():
    rows = []
    for case in CASES:
        res = run_case(*case, state=True, brute=True)
        gb = res.pop("ltmc_count_brute")
        assert abs(gb - res["ltmc_count"]) <= mpf(10) ** -20 * abs(gb), (case, gb, res["ltmc_count"])
        assert res["ltmc_state"] <= res["ltmc_count"] * (1 + mpf(10) ** -25)
        assert res["ltmc_count"] <= res["ltmc_best_r"] * (1 + mpf(10) ** -25)
        rows.append(res)
        print(f"{case}: best r={res['best_r']} {float(res['ltmc_best_r']):.10f} count "
              f"{float(res['ltmc_count']):.10f} [{res['count_policy']}] state {float(res['ltmc_state']):.10f}",
              flush=True)
    out = os.path.join(ROOT, "test", "reference", "optimal_reference.csv")
    with open(out, "w") as fh:
        fh.write("# Optimal policies (sec:optimal) by an independent mpmath implementation; see optimal_reference.py.\n")
        fh.write("# Count-based optimum also checked by brute force over all stopping sets.\n")
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["system", "model", "beta", "delay_type", "delay", "ratio", "setup",
                    "best_r", "ltmc_best_r", "ltmc_count", "count_policy", "ltmc_state"])
        for r in rows:
            w.writerow([r["system"], r["model"], r["beta"], r["delay_type"], r["delay"], r["ratio"], r["setup"],
                        r["best_r"], mp.nstr(r["ltmc_best_r"], 17), mp.nstr(r["ltmc_count"], 17),
                        r["count_policy"], mp.nstr(r["ltmc_state"], 17)])
    print("wrote", out)
    if "--mc" in sys.argv:
        print("Monte Carlo (|z| = |computed - simulated| / SE):")
        for case in [(1, "LFMO-C", "0.2", "exp", "0.1", "1", "0.1"), (1, "BFR", "0.5", "det", "0.3", "0.5", "0.05"),
                     (3, "LFMO-G", "0.2", "insp", "0.1", "2", "0.2"), (2, "BFR", "0.2", "none", "0", "0.2", "0.05"),
                     (3, "BFR", "0.2", "det", "0.1", "2", "0.2")]:
            res = run_case(*case, state=True)
            sysm, Psi, T, costs = res["_objs"]
            for label, act, g in [("count", res["_count_act"], res["ltmc_count"]),
                                  ("state", res["_state_act"], res["ltmc_state"])]:
                if label == "count":
                    decide = lambda x, f, a=act: a.get(f, False)
                else:
                    decide = lambda x, f, a=act: a.get(x, False)
                gs, se = simulate(sysm, Psi, case[3], case[4], costs, decide, 200000, seed=7)
                print(f"  {case} {label}: computed {float(g):.5f} simulated {gs:.5f} z={abs(float(g) - gs) / se:.2f}",
                      flush=True)

if __name__ == "__main__":
    main()
