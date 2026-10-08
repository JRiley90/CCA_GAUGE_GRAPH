"""Core computations for CCA-GaugeGraph v1.9.2: graph invariants, logical masks, C2/C3, analytic
crossing quantities, Rung-1 microscopic energies, and the binary64 crossing ladder."""
import json, collections, itertools, math
import numpy as np, scipy.sparse as sp
from scipy.sparse.linalg import eigsh
from scipy.optimize import brentq
from graph_data import EDGES, BASES

pc = lambda x: bin(x).count('1')
NV = {'T8': 8, 'Q3': 8, 'W8': 8, 'T10': 10, 'P5': 10}
J, K, G = 1.0, 0.5, 0.15

# ---------- graph invariants ----------
def connected(k):
    adj = collections.defaultdict(set)
    for u, v in EDGES[k]: adj[u].add(v); adj[v].add(u)
    seen, st = {0}, [0]
    while st:
        x = st.pop()
        for y in adj[x]:
            if y not in seen: seen.add(y); st.append(y)
    return len(seen) == NV[k]
def simple(k):
    es = [tuple(sorted(e)) for e in EDGES[k]]
    return len(set(es)) == len(es) and all(u != v for u, v in es)
def is_simple_cycle(ids, k):
    deg = collections.Counter(); adj = collections.defaultdict(set)
    for e in ids:
        u, v = EDGES[k][e]; deg[u] += 1; deg[v] += 1; adj[u].add(v); adj[v].add(u)
    if any(d != 2 for d in deg.values()): return False
    s = next(iter(adj)); seen = {s}; st = [s]
    while st:
        x = st.pop()
        for y in adj[x]:
            if y not in seen: seen.add(y); st.append(y)
    return len(seen) == len(adj)
def rank_f2(vs):
    red = []
    for x in vs:
        for r in red: x = min(x, x ^ r)
        if x: red.append(x); red.sort(reverse=True)
    return len(red)
def ke_distribution(k):
    c = collections.Counter(sum(1 for C in BASES[k] if e in C) for e in range(len(EDGES[k])))
    return {int(a): n for a, n in sorted(c.items())}
def masks(k):
    n = len(EDGES[k]); c = collections.Counter(sum(1 << i for i, C in enumerate(BASES[k]) if e in C) for e in range(n))
    c.pop(0, None); return c

# ---------- effective functionals ----------
def a_of(m, tau): return m * math.exp(-tau * (m - 1))
def C2(k, tau): return sum(a_of(m, tau) ** 2 / pc(v) for v, m in masks(k).items())
def C3(k, tau):
    mk = masks(k); a = {v: a_of(m, tau) for v, m in mk.items()}; t = 0.0
    for u in a:
        for w in a:
            r = u ^ w
            if r in a: t += a[u] * a[w] * a[r] / (pc(u) * pc(r))
    return t
def xor_closures(k):
    ks = list(masks(k)); return [t for t in itertools.combinations(ks, 3) if t[0] ^ t[1] ^ t[2] == 0]
def analytic():
    uc = (math.sqrt(226) - 8) / 18; tc = -0.5 * math.log(uc)
    d2 = -16 * math.exp(-2 * tc) - 36 * math.exp(-4 * tc)
    c3t = C3('T10', tc); dC3 = c3t - C3('P5', 0)
    return {'u_c': uc, 'tau_c2': tc, 'C2prime': d2, 'C3_T10_tauc': c3t, 'C3_P5': C3('P5', 0), 'dC3': dC3, 'a1': -dC3 / d2}

# ---------- Rung-1 microscopic energies ----------
def rung1(k, basis=None):
    basis = BASES[k] if basis is None else basis
    N = len(EDGES[k]); D = 1 << N; s = np.arange(D)
    par = lambda x: np.array([bin(i).count('1') & 1 for i in x])
    d = -J * sum(1 - 2 * ((s >> e) & 1) for e in range(N)).astype(float)
    for C in basis: d -= K * (1 - 2 * par(s & sum(1 << e for e in C)))
    H = sp.diags(d)
    for e in range(N): H = H + sp.csr_matrix((np.full(D, -G), (s, s ^ (1 << e))), shape=(D, D))
    H = H.tocsr(); w, v = eigsh(H, k=1, which='SA', v0=np.ones(D) / np.sqrt(D))
    return float(w[0]), float(np.linalg.norm(H @ v[:, 0] - w[0] * v[:, 0]))

# ---------- binary64 effective ladder ----------
def _HA(k, tau, g):
    M = np.zeros((64, 64)); s = np.arange(64)
    for x in s: M[x, x] = -K * sum(1 - 2 * ((x >> i) & 1) for i in range(6))
    for v, m in masks(k).items():
        a = a_of(m, tau)
        for x in s: M[x, x ^ v] -= g * a
    return M
def binary64_root(g):
    f = lambda t: np.linalg.eigvalsh(_HA('T10', t, g))[0] - np.linalg.eigvalsh(_HA('P5', t, g))[0]
    return brentq(f, 0.4, 0.5, xtol=1e-14, rtol=4 * 2.220446049250313e-16)
