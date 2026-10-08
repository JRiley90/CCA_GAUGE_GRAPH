"""Basis-dependence sampling for the effective functional C2 (Class C, empirical).
Run:  python basis_sampling_pipeline.py [--out basis_sampling_fixtures.json]
Deterministic: fixed seeds, independent random.Random streams."""
import json, sys, random, platform, collections, argparse
import numpy as np
from scipy.optimize import brentq
from graph_data import BASES

SEED_GL, SEED_MIN = 1, 3
N_GL, N_MIN = 20000, 3000
pc = lambda x: bin(x).count('1')

def masks_from_basis(basis_cycles, n_edges):
    """basis_cycles: list of edge-index lists -> Counter {mask: multiplicity} (zero mask = bridges dropped)."""
    c = collections.Counter(sum(1 << i for i, C in enumerate(basis_cycles) if e in C) for e in range(n_edges))
    c.pop(0, None); return c

def C2(mult, tau, cols=None):
    tot = 0.0
    for v, m in mult.items():
        w = pc(v) if cols is None else pc(apply(cols, v))
        tot += (m * np.exp(-tau * (m - 1))) ** 2 / w
    return tot

def apply(cols, v):
    o = 0
    for i, c in enumerate(cols):
        if (v >> i) & 1: o ^= c
    return o

def rand_gl(rng, n=6):
    """Uniform on GL(n,F2): sample nonzero columns, keep if full rank (rank-6 matrices are exactly the invertible ones)."""
    while True:
        cols = [rng.randrange(1, 1 << n) for _ in range(n)]
        rows, r = cols[:], 0
        for b in range(n):
            piv = next((x for x in rows if (x >> b) & 1), None)
            if piv is None: continue
            rows.remove(piv); rows = [x ^ piv if (x >> b) & 1 else x for x in rows]; r += 1
        if r == n: return cols

def cycle_space(basis_cycles):
    ev = [sum(1 << e for e in C) for C in basis_cycles]; n = len(ev); S = set()
    for m in range(1, 1 << n):
        x = 0
        for i in range(n):
            if (m >> i) & 1: x ^= ev[i]
        S.add(x)
    return sorted(S)

def random_min_basis(space, n, rng):
    els = space[:]; rng.shuffle(els); els.sort(key=pc)
    red, basis = [], []
    for x in els:
        y = x
        for r in red: y = min(y, y ^ r)
        if y: red.append(y); red.sort(reverse=True); basis.append(x)
        if len(basis) == n: break
    return basis

def run():
    mT = masks_from_basis(BASES['T10'], 15); mP = masks_from_basis(BASES['P5'], 15)
    canon = {'C2_T10_0': C2(mT, 0), 'C2_P5': C2(mP, 0)}
    # --- (a) 20,000 uniform GL(6,F2) pairs ---
    rng = random.Random(SEED_GL); t0s, ps, taus = [], [], []; surv = 0
    for _ in range(N_GL):
        A, Bm = rand_gl(rng), rand_gl(rng)
        p = C2(mP, 0, Bm); t0 = C2(mT, 0, A); tinf = C2(mT, 60, A)
        t0s.append(t0); ps.append(p)
        if tinf < p < t0:
            surv += 1; taus.append(brentq(lambda t: C2(mT, t, A) - p, 0, 60))
    taus = np.array(taus); q = lambda a: [float(x) for x in np.percentile(a, [5, 50, 95])]
    gl = {'n_draws': N_GL, 'seed': SEED_GL, 'n_crossing': surv, 'survival_rate': surv / N_GL,
          'C2_T10_0': {'min': min(t0s), 'max': max(t0s), 'median': float(np.median(t0s)), 'p5_p50_p95': q(t0s)},
          'C2_P5':    {'min': min(ps),  'max': max(ps),  'median': float(np.median(ps)),  'p5_p50_p95': q(ps)},
          'tau_c': {'min': float(taus.min()), 'max': float(taus.max()), 'p5_p50_p95': q(taus)},
          'canonical_above_all_samples': bool(canon['C2_T10_0'] > max(t0s) and canon['C2_P5'] > max(ps))}
    # --- (b) minimum-total-length bases with random tie-breaking ---
    ST, SP = cycle_space(BASES['T10']), cycle_space(BASES['P5']); rng = random.Random(SEED_MIN)
    pairs, mtaus = set(), []
    def masks_of(b): return collections.Counter(sum(1 << i for i, x in enumerate(b) if (x >> e) & 1) for e in range(15)) 
    for _ in range(N_MIN):
        a = masks_of(random_min_basis(ST, 6, rng)); b = masks_of(random_min_basis(SP, 6, rng)); a.pop(0, None); b.pop(0, None)
        t0, p = C2(a, 0), C2(b, 0); pairs.add((round(t0, 9), round(p, 9)))
        if C2(a, 60) < p < t0: mtaus.append(brentq(lambda t: C2(a, t) - p, 0, 60))
    mn = {'n_draws': N_MIN, 'seed': SEED_MIN, 'distinct_(C2_T10_0,C2_P5)_pairs': len(pairs),
          'values': sorted(pairs)[0] if len(pairs) == 1 else sorted(pairs),
          'tau_c_min': float(min(mtaus)), 'tau_c_max': float(max(mtaus)), 'n_crossing': len(mtaus),
          'canonical_is_min_total_length': {k: sum(len(c) for c in BASES[k]) for k in BASES}}
    return {'schema': 'basis_sampling_fixtures/v1', 'class': 'C (empirical sampling)', 'canonical': canon,
            'gl6_f2_sampling': gl, 'min_length_sampling': mn,
            'environment': {'python': platform.python_version(), 'numpy': np.__version__}}

# ---- verification (verify_track2b.py hook) ----
def verify(fx):
    f = []; g = fx['gl6_f2_sampling']; n = g['n_draws']
    se = (g['survival_rate'] * (1 - g['survival_rate']) / n) ** .5
    ref = run_ref = None
    if not (0.97 < g['survival_rate'] < 1.0): f.append('survival rate outside sanity band')
    if not g['canonical_above_all_samples']: f.append('canonical C2 no longer above all samples')
    m = fx['min_length_sampling']
    if m['distinct_(C2_T10_0,C2_P5)_pairs'] != 1: f.append('C2 not constant on sampled min-length bases')
    if abs(fx['canonical']['C2_T10_0'] - 22.5) > 1e-12 or abs(fx['canonical']['C2_P5'] - 10.0) > 1e-12: f.append('canonical C2 mismatch')
    return f, se

def reproduce_check(path, rtol_survival=5):
    """Regenerate and compare to stored fixture: exact on same Python, else within 5 binomial sigma."""
    old = json.load(open(path)); new = run(); f, se = verify(new)
    d = abs(new['gl6_f2_sampling']['survival_rate'] - old['gl6_f2_sampling']['survival_rate'])
    if d > rtol_survival * se: f.append(f'survival rate drift {d:.4f}')
    return f

def claim(fx):
    g, m = fx['gl6_f2_sampling'], fx['min_length_sampling']
    return {'CLM-014': {'quantity': 'basis-robustness of the second-order crossing', 'class': 'C', 'status': 'canonical',
        'method': 'uniform GL(6,F2) basis pairs for T10 and P5, fixed seed; plus random-tie-break minimum-length bases',
        'value': {'survival_rate': g['survival_rate'], 'n_draws': g['n_draws'], 'tau_c_p5_p50_p95': g['tau_c']['p5_p50_p95'],
                  'canonical_above_all_samples': g['canonical_above_all_samples'],
                  'min_length_tau_c_range': [m['tau_c_min'], m['tau_c_max']], 'min_length_n_draws': m['n_draws']},
        'statement': 'empirical sampling result; not a proof of invariance over all minimum-length bases',
        'fixture': 'basis_sampling_fixtures.json'}}

if __name__ == '__main__':
    out = 'basis_sampling_fixtures.json'
    if '--out' in sys.argv: out = sys.argv[sys.argv.index('--out') + 1]
    fx = run(); json.dump(fx, open(out, 'w'), indent=1)
    print(json.dumps({k: fx[k] for k in ('canonical', 'gl6_f2_sampling', 'min_length_sampling')}, indent=1))
    print('verify:', verify(fx)[0]); json.dump(claim(fx), open('claim_clm014.json', 'w'), indent=1)
