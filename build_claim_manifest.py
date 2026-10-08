"""Build CLAIM_MANIFEST.json (CLM-001..016) deterministically from fixtures."""
import json, sys, os
import highprecision_pipeline as hp, basis_sampling_pipeline as bs, w8_basis_sensitivity as w8
FIX = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'fixtures') + os.sep
def build():
    cv = json.load(open(FIX + 'canonical_values.json')); r1 = json.load(open(FIX + 'rung1_fixtures.json'))
    lad, _ = hp.load_ladder_fixture(FIX + 'ladder_highprecision.json')
    bsx = json.load(open(FIX + 'basis_sampling_fixtures.json')); w8x = json.load(open(FIX + 'w8_basis_fixture.json'))
    M = {}
    for k in ('T8', 'Q3', 'W8'):
        c = cv['rung1'][k]
        M[c['claim']] = {'quantity': f'E0({k})', 'class': 'B', 'status': 'canonical', 'value': c['E0'], 'residual_v1_9': c['residual'],
                         'reproduced_value': r1['energies'][k]['E0'], 'reproduced_residual': r1['energies'][k]['residual']}
    a = cv['analytic']
    M['CLM-004'] = {'quantity': 'u_c', 'class': 'A|D', 'status': 'canonical', 'value': a['u_c'], 'exact': '(sqrt(226)-8)/18'}
    M['CLM-005'] = {'quantity': 'tau_c^(2)', 'class': 'A|D', 'status': 'canonical', 'value': a['tau_c2'], 'exact': '-(1/2) ln u_c'}
    M['CLM-006'] = {'quantity': 'a1', 'class': 'A|D', 'status': 'canonical', 'value': a['a1'], 'C2prime': a['C2prime'], 'dC3': a['dC3']}
    M.update(hp.build_ladder_claims(lad))   # CLM-007,008,009,012,013 (B), CLM-010 (C), relabels 004-006
    for k in ('CLM-004', 'CLM-005', 'CLM-006'):
        M[k]['class'] = 'A|D'                  # keep value fields, enforce label
    M['CLM-004'].pop('update', None)
    M['CLM-011'] = {'quantity': 'attenuation ansatz', 'class': 'D', 'status': 'canonical', 'expression': 'm*exp(-tau*(m-1))',
                    'statement': 'postulated effective model, not a microscopic derivation'}
    M.update(bs.claim(bsx)); M.update(w8.claim(w8x))
    M['CLM-016'] = {'quantity': 'microscopic ground-state energies of T10 and P5', 'class': 'B', 'status': 'supplementary-canonical',
                    'value': {'E0_T10': r1['energies']['T10']['E0'], 'E0_P5': r1['energies']['P5']['E0'],
                              'ordering_T10<P5': r1['energies']['T10']['E0'] < r1['energies']['P5']['E0']},
                    'residuals': {'T10': r1['energies']['T10']['residual'], 'P5': r1['energies']['P5']['residual']},
                    'fixture': 'rung1_fixtures.json'}
    return {'schema': 'claim_manifest/v1', 'release': 'CCA-GaugeGraph-v1.9.2-frozen', 'n_claims': len(M),
            'claims': {k: M[k] for k in sorted(M)}}
def dumps(m): return json.dumps(m, indent=1, sort_keys=True, default=str) + '\n'
if __name__ == '__main__':
    out = os.path.join(FIX, '..', 'CLAIM_MANIFEST.json')
    m = build(); open(out, 'w').write(dumps(m)); print('wrote', out, m['n_claims'], 'claims')
