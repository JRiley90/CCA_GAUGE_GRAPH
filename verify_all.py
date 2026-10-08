"""Verification suite for CCA-GaugeGraph v1.9.2. Run: python verify_all.py   (exit code 0 = all pass)"""
import json, os, sys, hashlib, math
import numpy as np, mpmath as mp
import canonical_compute as cc, highprecision_pipeline as hp, basis_sampling_pipeline as bs, w8_basis_sensitivity as w8
import build_claim_manifest as bm
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')); FIX = os.path.join(ROOT, 'fixtures') + os.sep
cv = json.load(open(FIX + 'canonical_values.json')); results = []
def check(n, name, fails):
    results.append((n, name, not fails)); print(f"[{'PASS' if not fails else 'FAIL'}] {n:2d} {name}" + ('' if not fails else ' :: ' + '; '.join(fails)))

# 1 graphs connected & simple
f = []
for k in cc.EDGES:
    if not cc.connected(k): f.append(k + ' disconnected')
    if not cc.simple(k): f.append(k + ' not simple')
check(1, 'connectedness and simplicity', f)
# 2 cycle bases
f = []
for k in cc.EDGES:
    b1 = len(cc.EDGES[k]) - cc.NV[k] + 1
    if len(cc.BASES[k]) != b1: f.append(k + ' basis size')
    if cc.rank_f2([sum(1 << e for e in C) for C in cc.BASES[k]]) != b1: f.append(k + ' rank')
    if not all(cc.is_simple_cycle(C, k) for C in cc.BASES[k]): f.append(k + ' non-simple cycle')
check(2, 'cycle-basis size, rank, simple cycles', f)
# 3 Table 1
f = []
for k, s in cv['table1_sum_ke'].items():
    if sum(len(C) for C in cc.BASES[k]) != s: f.append(k + ' sum k_e')
exp = {'T8': {1: 8, 2: 4}, 'Q3': {1: 4, 2: 8}, 'W8': {1: 5, 2: 2, 3: 2, 4: 2, 5: 1}, 'T10': {1: 10, 2: 5}, 'P5': {1: 5, 2: 10}}
for k in exp:
    if cc.ke_distribution(k) != exp[k]: f.append(k + ' k_e distribution')
check(3, 'Table 1 k_e distributions and sums', f)
# 4 masks
f = []; mT = cc.masks('T10'); mult = sorted(mT.values())
if mult != [1] * 8 + [2] * 2 + [3]: f.append('T10 multiset')
if sorted(mT) != [1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48]: f.append('T10 mask set')
if sum(mT.values()) != 15: f.append('T10 total')
if set(cc.masks('P5').values()) != {1} or len(cc.masks('P5')) != 15: f.append('P5 distinct')
check(4, 'T10 mask multiset {1^8,2^2,3^1}; P5 all distinct', f)
# 5 XOR closures
f = []; XT = cc.xor_closures('T10'); XP = cc.xor_closures('P5')
if len(XT) != 5 or len(XP) != 10: f.append(f'counts {len(XT)},{len(XP)}')
if sorted(sorted(t) for t in XT) != [[1, 2, 3], [2, 4, 6], [4, 8, 12], [8, 16, 24], [16, 32, 48]]: f.append('T10 triples')
check(5, 'XOR closures (T10: 5, P5: 10)', f)
# 6 analytic
f = []; a = cc.analytic(); ref = cv['analytic']
for k in ref:
    if abs(a[k] - ref[k]) > 1e-12 * max(1, abs(ref[k])): f.append(f'{k}: {a[k]} vs {ref[k]}')
if abs(cc.C2('T10', 0) - 22.5) > 1e-12 or abs(cc.C2('P5', 0) - 10) > 1e-12: f.append('C2 at 0')
if abs(cc.C2('T10', a['tau_c2']) - 10) > 1e-10: f.append('crossing')
if abs(cc.C3('T10', 0) - 48) > 1e-12: f.append('C3(T10,0)')
check(6, 'analytic values u_c, tau_c, C2\', C3, dC3, a1', f)
# 7 Rung-1
f = []; r1 = json.load(open(FIX + 'rung1_fixtures.json'))['energies']
for k, c in cv['rung1'].items():
    if abs(r1[k]['E0'] - c['E0']) > 1e-10: f.append(f"{k} E0 differs {r1[k]['E0'] - c['E0']:.2e}")
for k, v in r1.items():
    if v['residual'] > 1e-9: f.append(k + ' residual')
    E, _ = cc.rung1(k)
    if abs(E - v['E0']) > 1e-10: f.append(k + ' fixture not reproduced')
if not (r1['T8']['E0'] < r1['Q3']['E0'] < r1['W8']['E0']): f.append('ordering')
check(7, 'Rung-1 energies vs frozen CLM-001..003, residuals, ordering', f)
# 8 T10/P5 supplementary
f = []
if not r1['T10']['E0'] < r1['P5']['E0']: f.append('T10<P5 ordering')
if abs(r1['T10']['E0'] + 18.10316) > 5e-6 or abs(r1['P5']['E0'] + 18.09386) > 5e-6: f.append('values vs manuscript')
check(8, 'T10/P5 microscopic energies (CLM-016)', f)
# 9 binary64 ladder vs HP within envelope
lad, rows = hp.load_ladder_fixture(FIX + 'ladder_highprecision.json'); f = []
fails, rep = hp.verify_ladder(lad); f += fails
for r in lad['rows']:
    g = r['g']; mine = cc.binary64_root(float(g)); d = abs(mine - float(mp.mpf(r['tau_c_hp40'])))
    if d > hp.noise_tol(g): f.append(f'g={g}: recomputed binary64 off by {d:.2e} > {hp.noise_tol(g):.2e}')
check(9, 'binary64 ladder (stored and recomputed) vs 40-digit reference within noise envelope (kappa=20)', f)
# 10 HP fixture consistency
f = []; ex = hp.a2_extrap(lad); a2 = ex['a2']; gs = ['0.0001', '0.0002', '0.0005', '0.001', '0.002']
vals = [float(a2[g]) for g in gs]
if not all(vals[i] > vals[i + 1] for i in range(4)): f.append('a2 not monotone')
for key in ('pair_1e-4_2e-4', 'pair_5e-4_1e-3'):
    if abs(float(ex[key]) + 8.404) > 2e-3: f.append(key)
if abs(float(ex['pair_1e-4_2e-4'] - ex['pair_5e-4_1e-3'])) > 1e-4: f.append('pairs disagree')
mp.mp.dps = 40; tc2 = -mp.log((mp.sqrt(226) - 8) / 18) / 2
s1 = (mp.mpf(rows['0.0005']['tau_c_hp40']) - tc2) / mp.mpf('0.0005'); s2 = (mp.mpf(rows['0.001']['tau_c_hp40']) - tc2) / mp.mpf('0.001')
if abs(float(2 * s1 - s2) - cv['analytic']['a1']) > 5e-5: f.append('a1 Richardson')
check(10, 'high-precision ladder: monotone a2, Richardson a2_extrap ~ -8.404, a1 cross-check', f)
# 11 basis sampling (rerun and compare)
f = bs.reproduce_check(FIX + 'basis_sampling_fixtures.json'); check(11, 'basis-sampling fixture reproduced and sane (CLM-014)', f)
# 12 W8
w8fx = json.load(open(FIX + 'w8_basis_fixture.json')); f = w8.verify(w8fx)
f += ['not reproduced'] if abs(w8.run()['shift_W8'] - w8fx['shift_W8']) > 1e-10 else []
check(12, 'W8 minimum-length basis sensitivity (CLM-015)', f)
# 13 manifest determinism and consistency
f = []; m1 = bm.dumps(bm.build()); m2 = bm.dumps(bm.build())
if m1 != m2: f.append('non-deterministic')
disk = open(os.path.join(ROOT, 'CLAIM_MANIFEST.json')).read()
if disk != m1: f.append('CLAIM_MANIFEST.json differs from regenerated manifest')
man = json.loads(m1)
if sorted(man['claims']) != [f'CLM-{i:03d}' for i in range(1, 17)]: f.append('claim ids')
if any(man['claims'][k]['class'] != 'A|D' for k in ('CLM-004', 'CLM-005', 'CLM-006')): f.append('A|D labels')
check(13, 'manifest deterministic, matches disk, CLM-001..016, A|D labels', f)
# 14 SHA256SUMS
f = []; sp_ = os.path.join(ROOT, 'SHA256SUMS')
if not os.path.exists(sp_): f.append('SHA256SUMS missing (run generate_hashes.py)')
else:
    for line in open(sp_).read().splitlines():
        h, p = line.split('  ', 1)
        if hashlib.sha256(open(os.path.join(ROOT, p), 'rb').read()).hexdigest() != h: f.append('hash mismatch ' + p)
check(14, 'SHA256SUMS verification', f)
n_ok = sum(r[2] for r in results); print(f'\n{n_ok}/{len(results)} checks passed'); sys.exit(0 if n_ok == len(results) else 1)
