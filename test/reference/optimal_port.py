"""Python port of src/optimal.jl (Float64), used to validate the algorithm design.

It mirrors the Julia code function by function:
  call_tables      after-call quantities a^{(w)}_i, b^{(w)}_i (prop:after call det) and the
                   one-period quantities of lem:period, for every w and i, from the
                   survivor-count chain with nonnegative-term methods (resolvent or
                   uniformization), i.e. without the alternating sums of lem:pbbs and exps;
  count_pass       alg:count exp and its det / none / insp variants (conditional signature);
  StateSpace       levels and ranks of the 2^n states;
  state_pre        per-state signature (lem:level sums with u = Phi) and per-state after-call
                   quantities, computed once per (system, model, delay);
  state_pass       alg:state exp / alg:state det / alg:state fixed (level-wise sweep, two
                   quantities, closed forms at failed states);
  dinkelbach       alg:dinkelbach.
Compared with optimal_reference.py (independent mpmath implementation) by
`python3 test/reference/optimal_port.py`.
"""
import math, os, sys
sys.path.insert(0, os.path.dirname(__file__))
import ltmc_reference as R
import ctmc_port as C
from math import comb

def rates_float(Psi, n):
    Rt, ex = C.rates(Psi, n)
    return [[float(v) for v in row] for row in Rt], [float(v) for v in ex]

# ----------------------------------------------------------------------------- tables
def _poisson(lam):
    M = int(math.ceil(lam + 10 * math.sqrt(lam) + 25))
    p = [math.exp(-lam)]
    for m in range(1, M + 1):
        p.append(p[-1] * lam / m)
    tail = [0.0] * (M + 1); acc = 0.0
    for m in range(M, -1, -1):
        tail[m] = acc; acc += p[m]
    return p, tail

def _row0(Rw, exw, i, kind, value):
    """Law E0 of the count at the end of the delay, and occupation F0 of counts < i,
    for the w-component chain frozen at counts >= i, started at count 0."""
    w = len(exw) - 1
    E0 = [0.0] * (w + 1); F0 = [0.0] * (w + 1)
    if kind == "exp":
        gam = 1.0 / value
        X = [0.0] * (w + 1)
        X[0] = 1.0 / (gam + exw[0])
        for c in range(1, w + 1):
            acc = 0.0
            for c1 in range(0, min(c, i)):
                acc += X[c1] * Rw[c1][c]
            X[c] = acc / (gam + (exw[c] if c < i else 0.0))
        for c in range(w + 1):
            E0[c] = gam * X[c]
            F0[c] = X[c] if c < i else 0.0
        return E0, F0
    # det / insp: uniformization with rate q = exw[0] = Psi(w)
    t = value; q = exw[0]
    p, tail = _poisson(q * t)
    v = [0.0] * (w + 1); v[0] = 1.0
    for m in range(len(p)):
        for c in range(w + 1):
            E0[c] += p[m] * v[c]
            F0[c] += tail[m] / q * v[c]
        nv = [0.0] * (w + 1)
        for c in range(w + 1):
            if v[c] == 0.0: continue
            if c < i:
                nv[c] += v[c] * (1.0 - exw[c] / q)
                for c2 in range(c + 1, w + 1):
                    nv[c2] += v[c] * Rw[c][c2] / q
            else:
                nv[c] += v[c]
        v = nv
    E0[0] = math.exp(-exw[0] * t)
    F0[0] = -math.expm1(-exw[0] * t) / exw[0]
    for c in range(i, w + 1):
        F0[c] = 0.0
    return E0, F0

def call_tables(Rt, ex, n, kind, value):
    """pf[w][i], enA[w][i], enF[w][i], ell[w][i] for w = 1..n, i = 1..w; pN[w][j] (insp)."""
    T = dict(pf={}, enA={}, enF={}, ell={}, pN={})
    for w in range(1, n + 1):
        f = n - w
        Rw = [[Rt[f + a][f + b] for b in range(w + 1)] for a in range(w + 1)]
        exw = [ex[f + a] for a in range(w + 1)]
        for i in range(1, w + 1):
            if kind == "none":
                T["pf"][w, i] = T["enA"][w, i] = T["enF"][w, i] = T["ell"][w, i] = 0.0
                continue
            E0, F0 = _row0(Rw, exw, i, kind, value)
            T["pf"][w, i] = sum(E0[c] for c in range(i, w + 1))
            T["enA"][w, i] = sum(c * E0[c] for c in range(w + 1))
            T["enF"][w, i] = sum(c * E0[c] for c in range(i, w + 1))
            T["ell"][w, i] = sum(F0[c] for c in range(i))
            if i == w and kind == "insp":
                T["pN"][w] = E0[:]
    return T

# ----------------------------------------------------------------------------- costs
class Costs:
    def __init__(self, n, rho, sigma, csys=1.0):
        self.csys = csys; self.cset = sigma * csys; self.ccmp = rho * csys / n
    def rep(self, f): return self.cset + self.ccmp * f
    def fail(self, f): return self.cset + self.ccmp * f + self.csys

# ----------------------------------------------------------------------------- count-based
class CountModel:
    def __init__(self, n, W, Rt, ex, kind, value):
        self.n, self.kind = n, kind
        self.Sb = [W[f] / comb(n, f) for f in range(n + 1)]
        self.Rt, self.ex = Rt, ex
        T = call_tables(Rt, ex, n, kind, value)
        self.T = T
        self.A = {}       # per f: (Pf, EN, L) after a call, or one-period (PfD, ENfD, LD) for insp
        self.q = {}
        for f in range(n):
            if self.Sb[f] == 0: continue
            w = n - f
            sig = [(self.Sb[f + i - 1] - self.Sb[f + i]) / self.Sb[f] for i in range(1, w + 1)]
            pf = sum(sig[i - 1] * T["pf"][w, i] for i in range(1, w + 1))
            en = sum(sig[i - 1] * (T["enF"][w, i] if kind == "insp" else T["enA"][w, i]) for i in range(1, w + 1))
            ell = sum(sig[i - 1] * T["ell"][w, i] for i in range(1, w + 1))
            self.A[f] = (pf, en, ell)
            if kind == "insp":
                self.q[f] = [T["pN"][w][j] * self.Sb[f + j] / self.Sb[f] for j in range(0, w)]

    def call_allowed(self, f):
        return not ((self.kind == "none" or self.kind == "insp") and f == 0)

    def pass_(self, g, costs, forced=None):
        n, Sb = self.n, self.Sb
        c = [0.0] * (n + 1); l = [0.0] * (n + 1); act = {}
        for f in range(n - 1, -1, -1):
            if Sb[f] == 0: continue
            w = n - f
            pf, en, ell = self.A[f]
            if self.kind == "insp":
                q = self.q[f]
                cD = (costs.cset + costs.csys + costs.ccmp * f) * pf + costs.ccmp * en
                cw = cD; lw = ell
                for j in range(1, w):
                    if q[j] > 0:
                        cw += q[j] * c[f + j]; lw += q[j] * l[f + j]
                cw /= (1 - q[0]); lw /= (1 - q[0])
                ccall, lcall = costs.rep(f), 0.0
                ccont, lcont = cw, lw
            else:
                ccall = costs.cset + costs.ccmp * (f + en) + costs.csys * pf
                lcall = ell
                ccont = 0.0; lcont = 1.0 / self.ex[f]
                for j in range(1, w + 1):
                    Pj = self.Rt[f][f + j] / self.ex[f]
                    ps = Sb[f + j] / Sb[f]
                    ccont += Pj * (ps * c[f + j] + (1 - ps) * costs.fail(f + j))
                    lcont += Pj * ps * l[f + j]
            if not self.call_allowed(f):
                call = False
            elif forced is not None:
                call = forced(f)
            else:
                call = ccall - g * lcall <= ccont - g * lcont
            act[f] = call
            c[f], l[f] = (ccall, lcall) if call else (ccont, lcont)
        return act, c[0], l[0]

def dinkelbach(passf, g0, maxit=100):
    g = g0
    for it in range(1, maxit + 1):
        act, c, l = passf(g)
        if c - g * l >= -1e-12 * max(1.0, abs(c)):
            return c / l, act, it
        g = c / l
    raise RuntimeError("no convergence")

# ----------------------------------------------------------------------------- state-based
class StateSpace:
    def __init__(self, n, phi):
        self.n = n
        self.phi = [bool(phi(x)) for x in range(1 << n)]
        self.levels = [[] for _ in range(n + 1)]
        self.rank = [0] * (1 << n)
        for x in range(1 << n):
            w = bin(x).count("1")
            self.rank[x] = len(self.levels[w]); self.levels[w].append(x)

def state_pre(S, cm):
    """Per-state signature by lem:level sums (u = Phi, integers) and after-call quantities."""
    n, T, kind = S.n, cm.T, cm.kind
    A = {}
    prevW = None
    for w in range(n + 1):
        curW = []
        for x in S.levels[w]:
            if not S.phi[x]:
                curW.append([0] * (w + 1)); continue
            U = [0] * (w + 1); U[0] = 1
            y = x
            while y:
                b = y & -y; y ^= b
                son = S.rank[x ^ b]
                col = prevW[son]
                for j in range(1, w + 1):
                    U[j] += col[j - 1]
            for j in range(1, w + 1):
                assert U[j] % j == 0
                U[j] //= j
            curW.append(U)
            f = n - w
            if f < n:
                Sx = [U[i] / comb(w, i) for i in range(w + 1)]
                sig = [Sx[i - 1] - Sx[i] for i in range(1, w + 1)]
                pf = sum(sig[i - 1] * T["pf"][w, i] for i in range(1, w + 1))
                en = sum(sig[i - 1] * (T["enF"][w, i] if kind == "insp" else T["enA"][w, i]) for i in range(1, w + 1))
                ell = sum(sig[i - 1] * T["ell"][w, i] for i in range(1, w + 1))
                A[x] = (pf, en, ell)
        prevW = curW
    return A

def state_pass(S, cm, A, g, costs, forced=None):
    n, kind = S.n, cm.kind
    prevC = prevL = None
    calls = [0] * (n + 1); working = [0] * (n + 1)
    for w in range(n + 1):
        f = n - w
        curC, curL = [], []
        for x in S.levels[w]:
            if not S.phi[x]:
                if kind == "insp":
                    curC.append([0.0] * (w + 1))
                else:
                    curC.append([comb(w, j) * costs.fail(f + j) for j in range(w + 1)])
                curL.append([0.0] * (w + 1))
                continue
            UC = [0.0] * (w + 1); UL = [0.0] * (w + 1)
            y = x
            while y:
                b = y & -y; y ^= b
                r = S.rank[x ^ b]
                cc, cl = prevC[r], prevL[r]
                for j in range(1, w + 1):
                    UC[j] += cc[j - 1]; UL[j] += cl[j - 1]
            for j in range(1, w + 1):
                UC[j] /= j; UL[j] /= j
            pf, en, ell = A[x]
            if kind == "insp":
                pN = cm.T["pN"][w]
                cw = (costs.cset + costs.csys + costs.ccmp * f) * pf + costs.ccmp * en
                lw = ell
                for j in range(1, w + 1):
                    pj = pN[j] / comb(w, j)
                    cw += pj * UC[j]; lw += pj * UL[j]
                cw /= (1 - pN[0]); lw /= (1 - pN[0])
                ccall, lcall, ccont, lcont = costs.rep(f), 0.0, cw, lw
            else:
                ccall = costs.cset + costs.ccmp * (f + en) + costs.csys * pf
                lcall = ell
                ccont = 0.0; lcont = 1.0 / cm.ex[f]
                for j in range(1, w + 1):
                    pj = cm.Rt[f][f + j] / cm.ex[f] / comb(w, j)
                    ccont += pj * UC[j]; lcont += pj * UL[j]
            if not cm.call_allowed(f):
                call = False
            elif forced is not None:
                call = forced(x, f)
            else:
                call = ccall - g * lcall <= ccont - g * lcont
            working[f] += 1; calls[f] += call
            v, t = (ccall, lcall) if call else (ccont, lcont)
            UC[0] = v; UL[0] = t
            curC.append(UC); curL.append(UL)
        prevC, prevL = curC, curL
    return (calls, working), prevC[0][0], prevL[0][0]

# ----------------------------------------------------------------------------- driver
def run_case(sid, model, beta, kind, value, rho, sigma, state=True):
    s = R.load_signature(f"yeh2024_{sid:02d}"); n = len(s)
    import optimal_reference as O
    W = O.load_W(f"yeh2024_{sid:02d}")
    Psi = R.psi_fun(model, R.calibrate(model, beta))
    Rt, ex = rates_float(Psi, n)
    costs = Costs(n, float(rho), float(sigma))
    cm = CountModel(n, W, Rt, ex, kind, float(value))
    # r-out-of-n with the CTMC engine (r = 0 included when admissible)
    r0 = kind in ("exp", "det")
    q = C.policy_ctmc(s, Psi, n, kind, value, r0=r0)
    q = q if r0 else [None] + q                     # q[r] for r = 0..n
    rs = list(range(0 if kind in ("exp", "det") else 1, n + 1))
    L = {r: (costs.cset + costs.ccmp * q[r][1] + costs.csys * q[r][0]) / q[r][2] for r in rs}
    rbest = min(L, key=lambda r: L[r])
    gc, actc, itc = dinkelbach(lambda g: cm.pass_(g, costs), L[rbest])
    out = dict(n=n, best_r=rbest, ltmc_best_r=L[rbest], ltmc_count=gc, act_count=actc, it_count=itc,
               L=L, cm=cm, costs=costs)
    if state:
        edges = R.load_edges(sid)
        phi = lambda x: R.works(edges, [(x >> e) & 1 == 1 for e in range(n)])
        S = StateSpace(n, phi)
        A = state_pre(S, cm)
        gs, acts, its = dinkelbach(lambda g: state_pass(S, cm, A, g, costs), gc)
        out.update(ltmc_state=gs, it_state=its, act_state=acts, S=S, A=A)
    return out

if __name__ == "__main__":
    import csv
    path = os.path.join(R.ROOT, "test", "reference", "optimal_reference.csv")
    rows = list(csv.DictReader(l for l in open(path) if not l.startswith("#")))
    worst = {"ltmc_best_r": 0, "ltmc_count": 0, "ltmc_state": 0}
    bad_policy = 0
    for r in rows:
        o = run_case(int(r["system"]), r["model"], r["beta"], r["delay_type"], r["delay"], r["ratio"], r["setup"])
        for k in worst:
            worst[k] = max(worst[k], abs(o[k] - float(r[k])) / float(r[k]))
        pol = "".join("-" if o["cm"].Sb[f] == 0 else ("C" if o["act_count"].get(f) else ".") for f in range(o["n"]))
        if pol != r["count_policy"] or o["best_r"] != int(r["best_r"]):
            bad_policy += 1
            print("policy mismatch", r, pol, o["best_r"])
    print(f"{len(rows)} reference cases: worst relative differences {worst}; policy/r mismatches {bad_policy}")
