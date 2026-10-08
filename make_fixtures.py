"""Regenerate rung1_fixtures.json, graphs.json, and ladder_highprecision.json (40 dps, confirmed at 60)."""
import json, sys, platform, mpmath as mp, numpy as np, scipy
import canonical_compute as cc, highprecision_pipeline as hp
from graph_data import EDGES, BASES
FIX = '../fixtures/'
def graphs():
    json.dump({'edges': EDGES, 'bases': BASES, 'table1': {k: {'V': cc.NV[k], 'E': len(EDGES[k]), 'beta1': len(BASES[k]),
        'sum_ke': sum(len(c) for c in BASES[k]), 'ke_dist': cc.ke_distribution(k)} for k in EDGES}}, open(FIX + 'graphs.json', 'w'), indent=1)
def rung1():
    d = {}
    for k in EDGES:
        E, r = cc.rung1(k); d[k] = {'E0': E, 'residual': r}
    json.dump({'schema': 'rung1_fixtures/v1', 'params': {'J': 1.0, 'K': 0.5, 'g': 0.15}, 'energies': d,
               'environment': env()}, open(FIX + 'rung1_fixtures.json', 'w'), indent=1)
def env(): return {'python': platform.python_version(), 'numpy': np.__version__, 'scipy': scipy.__version__, 'mpmath': mp.__version__}
def ladder():
    gs = ['0.002', '0.001', '0.0005', '0.0002', '0.0001']
    r40 = hp.compute_highprecision_ladder(cc.masks('T10'), cc.masks('P5'), gs, dps=40, tol_exp=30)
    r60 = hp.compute_highprecision_ladder(cc.masks('T10'), cc.masks('P5'), gs, dps=60, tol_exp=45)
    canon = json.load(open(FIX + 'canonical_values.json'))['ladder_binary64_canonical']
    rows = [{'g': g, 'tau_c_hp40': mp.nstr(r40[g], 25), 'tau_c_hp60': mp.nstr(r60[g], 25),
             'tau_c_binary64_canonical': canon[g]} for g in gs]
    json.dump({'mpmath_version': mp.__version__, 'working_precision_digits': 40, 'confirmation_precision_digits': 60,
        'solver': 'mpmath.findroot(anderson), tol 1e-30 (40 dps) / 1e-45 (60 dps); eigenvalues via mpmath.eigsy',
        'binary64_solver': 'scipy brentq xtol=1e-14 rtol=4*eps',
        'note': 'arbitrary-precision floating point, not interval enclosures', 'rows': rows, 'environment': env()},
        open(FIX + 'ladder_highprecision.json', 'w'), indent=1)
if __name__ == '__main__':
    for a in sys.argv[1:] or ['graphs', 'rung1', 'ladder']: globals()[a]()
