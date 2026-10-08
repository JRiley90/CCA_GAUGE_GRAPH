"""W8 basis sensitivity of the Rung-1 microscopic benchmark (Class B, supplementary).
Run: python w8_basis_sensitivity.py [--out w8_basis_fixture.json]"""
import json, sys, platform, collections
import numpy as np, scipy.sparse as sp
from scipy.sparse.linalg import eigsh
from graph_data import EDGES, BASES
J, K, G = 1.0, 0.5, 0.15
pc = lambda x: bin(x).count('1')

def is_simple_cycle(edge_ids, edges):
    deg = collections.Counter(); adj = collections.defaultdict(set)
    for e in edge_ids:
        u, v = edges[e]; deg[u] += 1; deg[v] += 1; adj[u].add(v); adj[v].add(u)
    if any(d != 2 for d in deg.values()): return False
    start = next(iter(adj)); seen = {start}; st = [start]
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

def cycle_space(basis):
    ev = [sum(1 << e for e in C) for C in basis]; n = len(ev); S = set()
    for m in range(1, 1 << n):
        x = 0
        for i in range(n):
            if (m >> i) & 1: x ^= ev[i]
        S.add(x)
    return S

def min_length_basis(basis):
    """Deterministic greedy (weight, integer order): matroid greedy => minimum total edge-length basis."""
    n = len(basis); els = sorted(cycle_space(basis), key=lambda x: (pc(x), x)); red, out = [], []
    for x in els:
        y = x
        for r in red: y = min(y, y ^ r)
        if y: red.append(y); red.sort(reverse=True); out.append(x)
        if len(out) == n: break
    return [[e for e in range(len(EDGES['W8'])) if (x >> e) & 1] for x in out]

def ground_energy(n_edges, basis):
    D = 1 << n_edges; s = np.arange(D)
    par = lambda x: np.array([bin(i).count('1') & 1 for i in x])
    d = -J * sum(1 - 2 * ((s >> e) & 1) for e in range(n_edges)).astype(float)
    for C in basis: d -= K * (1 - 2 * par(s & sum(1 << e for e in C)))
    H = sp.diags(d)
    for e in range(n_edges): H = H + sp.csr_matrix((np.full(D, -G), (s, s ^ (1 << e))), shape=(D, D))
    H = H.tocsr(); w, v = eigsh(H, k=1, which='SA', v0=np.ones(D) / np.sqrt(D))
    return float(w[0]), float(np.linalg.norm(H @ v[:, 0] - w[0] * v[:, 0]))

def run():
    can = BASES['W8']; mn = min_length_basis(can)
    checks = {'canonical_cycles_simple': all(is_simple_cycle(c, EDGES['W8']) for c in can),
              'min_basis_cycles_simple': all(is_simple_cycle(c, EDGES['W8']) for c in mn),
              'min_basis_rank': rank_f2([sum(1 << e for e in c) for c in mn]),
              'canonical_total_length': sum(len(c) for c in can), 'min_total_length': sum(len(c) for c in mn),
              'same_cycle_space': cycle_space(can) == cycle_space(mn)}
    e = {}
    for name, b in (('W8_canonical', can), ('W8_min_length', mn), ('T8', BASES['T8']), ('Q3', BASES['Q3'])):
        E, r = ground_energy(12, b); e[name] = {'E0': E, 'residual': r}
    order = lambda w8: e['T8']['E0'] < e['Q3']['E0'] < w8
    return {'schema': 'w8_basis_fixture/v1', 'class': 'B (supplementary)', 'min_length_basis': mn, 'checks': checks,
            'energies': e, 'shift_W8': e['W8_min_length']['E0'] - e['W8_canonical']['E0'],
            'ordering_T8<Q3<W8_canonical': order(e['W8_canonical']['E0']),
            'ordering_T8<Q3<W8_min_length': order(e['W8_min_length']['E0']),
            'note': 'minimum-length basis is not unique; this is the deterministic-greedy choice',
            'environment': {'python': platform.python_version(), 'numpy': np.__version__}}

def verify(fx, ref_canonical=-14.568567641496603, tol=1e-10):
    f = []; c = fx['checks']
    if not (c['canonical_cycles_simple'] and c['min_basis_cycles_simple']): f.append('non-simple cycle in basis')
    if c['min_basis_rank'] != 5 or not c['same_cycle_space']: f.append('min basis not a basis of the same cycle space')
    if c['min_total_length'] >= c['canonical_total_length']: f.append('min basis not shorter')
    if abs(fx['energies']['W8_canonical']['E0'] - ref_canonical) > tol: f.append('canonical W8 energy disagrees with CLM-003')
    if max(v['residual'] for v in fx['energies'].values()) > 1e-9: f.append('residual too large')
    if not fx['ordering_T8<Q3<W8_min_length']: f.append('ordering changed under min-length basis')
    return f

def claim(fx):
    e = fx['energies']
    return {'CLM-015': {'quantity': 'basis sensitivity of the Rung-1 microscopic benchmark (W8)', 'class': 'B', 'status': 'supplementary-canonical',
        'value': {'E0_W8_min_length_basis': e['W8_min_length']['E0'], 'E0_W8_canonical_basis': e['W8_canonical']['E0'],
                  'shift': fx['shift_W8'], 'canonical_total_length': fx['checks']['canonical_total_length'],
                  'min_total_length': fx['checks']['min_total_length'], 'ordering_preserved': fx['ordering_T8<Q3<W8_min_length'],
                  'residual': e['W8_min_length']['residual']},
        'statement': 'Hmicro depends on the chosen cycle basis; ordering T8<Q3<W8 is unchanged for the greedy minimum-length W8 basis',
        'fixture': 'w8_basis_fixture.json'}}

if __name__ == '__main__':
    out = sys.argv[sys.argv.index('--out') + 1] if '--out' in sys.argv else 'w8_basis_fixture.json'
    fx = run(); json.dump(fx, open(out, 'w'), indent=1); print(json.dumps(fx, indent=1))
    print('verify:', verify(fx)); json.dump(claim(fx), open('claim_clm015.json', 'w'), indent=1)
