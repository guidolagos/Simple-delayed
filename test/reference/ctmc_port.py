"""Python port of the Julia engine (src/ltmc.jl), generic in the number type.

Used only to validate the algorithm before it runs on the cluster:
  float CTMC  vs  mpmath closed forms of the paper   (systems with n <= 12)
  float CTMC  vs  mpmath CTMC                         (n = 30, float stability)
"""
import math, sys, os
from mpmath import mp, mpf
sys.path.insert(0, os.path.dirname(__file__))
import ltmc_reference as R

def rates(Psi, n):
    """Rt[c][c2] = rate from c to c2 > c failed components (exact, then converted)."""
    L = R.LFMO(Psi, n)
    Rt = [[mpf(0)] * (n + 1) for _ in range(n + 1)]
    for c in range(n):
        w = n - c
        for f in range(1, w + 1):
            Rt[c][c + f] = L.P[w] * L.Pw(w, f)
    return Rt, [L.P[n - c] for c in range(n + 1)]

def poisson_weights(lam, one, exp_):
    M = int(math.ceil(float(lam) + 10 * math.sqrt(float(lam)) + 25))
    p = [exp_(-lam)]
    for m in range(1, M + 1):
        p.append(p[-1] * lam / m)
    tail = [0 * one] * (M + 1)
    acc = 0 * one
    for m in range(M, -1, -1):
        tail[m] = acc          # P(N > m), ignoring mass beyond M
        acc += p[m]
    return p, tail

def unif(Q, rows, q, t, one, exp_, expm1_):
    """E = exp(Q t), F = int_0^t exp(Q s) ds, via uniformization. Q upper triangular."""
    N = len(Q)
    P = [[(1 if i == j else 0) * one + Q[i][j] / q for j in range(N)] for i in range(N)]
    p, tail = poisson_weights(q * t, one, exp_)
    Pm = [[(1 if i == j else 0) * one for j in range(N)] for i in range(N)]
    E = [[0 * one] * N for _ in range(N)]; F = [[0 * one] * N for _ in range(N)]
    for m in range(len(p)):
        for i in rows:
            for j in range(i, N):
                E[i][j] += p[m] * Pm[i][j]
                F[i][j] += tail[m] * Pm[i][j] / q
        Pm = [[sum(Pm[i][l] * P[l][j] for l in range(i, j + 1)) for j in range(N)] if i in rows
              else Pm[i] for i in range(N)]
    for i in rows:                      # exact diagonal (triangular generator)
        E[i][i] = exp_(Q[i][i] * t)
        F[i][i] = t * one if Q[i][i] == 0 else expm1_(Q[i][i] * t) / Q[i][i]
    return E, F

def policy_ctmc(s, Psi, n, dtype, d, num="float"):
    if num == "float":
        one, exp_, expm1_ = 1.0, math.exp, math.expm1
        conv = float
    else:
        one, exp_, expm1_ = mpf(1), mp.exp, mp.expm1
        conv = lambda x: x
    Rt, ex = rates(Psi, n)
    Rt = [[conv(x) for x in row] for row in Rt]; ex = [conv(x) for x in ex]
    d = conv(mpf(d))
    out = [[0 * one] * 3 for _ in range(n + 1)]     # out[r] = (A, B, C)
    for k in range(1, n + 1):
        sk = s[k - 1]
        if sk == 0:
            continue
        sk = conv(mpf(sk.numerator) / sk.denominator)
        N = n + 1
        Q = [[0 * one] * N for _ in range(N)]
        for c in range(k):
            Q[c][c] = -ex[c]
            for c2 in range(c + 1, N):
                Q[c][c2] = Rt[c][c2]
        rows = list(range(k))
        if dtype == "exp":
            g = 1 / d
            X = [[0 * one] * N for _ in range(N)]
            for c0 in range(k):
                X[c0][c0] = 1 / (g - Q[c0][c0])
                for c in range(c0 + 1, N):
                    X[c0][c] = sum(X[c0][c1] * Q[c1][c] for c1 in range(c0, min(c, k))) / (g - Q[c][c])
            Mst = [[g * x for x in row] for row in X]; F = X
        elif dtype in ("det", "insp"):
            Mst, F = unif(Q, rows, ex[0], d, one, exp_, expm1_)
        for r in range(1, n + 1):
            m = min(r, k)
            if dtype == "insp":
                G = [0 * one] * m
                for c in range(m):
                    acc = (1 if c == 0 else 0) * one + sum(G[c1] * Mst[c1][c] for c1 in range(c))
                    G[c] = acc / (-expm1_(Q[c][c] * d))
                A = sum(G[c1] * Mst[c1][c] for c1 in range(m) for c in range(k, N))
                Bv = sum(G[c1] * c * Mst[c1][c] for c1 in range(m) for c in range(m, N))
                C = sum(G[c1] * F[c1][c] for c1 in range(m) for c in range(k))
            else:
                v = [0 * one] * m; v[0] = one
                for c in range(1, m):
                    v[c] = sum(v[c1] * Rt[c1][c] / ex[c1] for c1 in range(c))
                h = [0 * one] * N
                for c in range(m, N):
                    h[c] = sum(v[c1] * Rt[c1][c] / ex[c1] for c1 in range(m))
                C = sum(v[c] / ex[c] for c in range(m))
                A = sum(h[c] for c in range(k, N))
                if dtype == "none" or r >= k:
                    Bv = sum(c * h[c] for c in range(m, N))
                else:
                    Bv = sum(c * h[c] for c in range(k, N))
                    for c0 in range(r, k):
                        A += h[c0] * sum(Mst[c0][c] for c in range(k, N))
                        Bv += h[c0] * sum(c * Mst[c0][c] for c in range(N))
                        C += h[c0] * sum(F[c0][c] for c in range(k))
            out[r][0] += sk * A; out[r][1] += sk * Bv; out[r][2] += sk * C
    return out[1:]

if __name__ == "__main__":
    import time
    worst = 0.0
    for system in (1, 3, 5):
        s = R.load_signature(f"yeh2024_{system:02d}"); n = len(s)
        for model, beta in [("Exp", "0"), ("LFMO-C", "0.2"), ("LFMO-G", "0.05"), ("LFMO-IG", "0.5"),
                            ("LFMO-S", "0.2"), ("BFR", "0.5"), ("LSh", "0.2")]:
            Psi = R.psi_fun(model, R.calibrate(model, beta))
            for dtype, d in [("none", "0"), ("exp", "0.1"), ("det", "0.3"), ("insp", "0.1"), ("insp", "1")]:
                ct = policy_ctmc(s, Psi, n, dtype, d)
                rs = range(1, n + 1) if n <= 8 else (1, 3, n // 2, n)
                for r in rs:
                    cf = R.policy_closed(s, Psi, n, r, dtype, mpf(d))
                    err = max(abs(float(a) - b) / max(1.0, abs(float(a))) for a, b in zip(cf, ct[r - 1]))
                    worst = max(worst, err)
            print(f"system {system} (n={n}) {model}: max rel err so far {worst:.2e}", flush=True)
    print("float CTMC vs mp closed forms, worst relative error:", worst)
