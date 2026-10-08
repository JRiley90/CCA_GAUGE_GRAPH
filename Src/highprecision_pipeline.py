"""v1.9.2 additions for CCA-GaugeGraph. Drop into src/ (merge into canonical_compute.py,
verify_track2b.py, build_claim_manifest.py, generate_hashes.py as marked)."""
import json, hashlib, io, os, subprocess, tarfile, sys
import mpmath as mp

C2PRIME_AT_TC2 = -11.74818099700186      # Class A|D, manuscript Eq. (a1)
KAPPA = 20.0
EPS_ABS = 1e-15                          # absolute eigenvalue error scale

# ---------- (1) canonical_compute.py: high-precision ladder as Class B reference ----------
def compute_highprecision_ladder(mult_T10, mult_P5, g_values, dps=40, tol_exp=30):
    """mult_*: {mask:int -> multiplicity}. Returns {g:str -> mp root}. K=1/2, 6 logical qubits."""
    mp.mp.dps = dps
    K = mp.mpf(1) / 2
    def H(mult, tau, g):
        M = mp.zeros(64, 64)
        for x in range(64):
            M[x, x] = -K * sum(1 - 2 * ((x >> i) & 1) for i in range(6))
        for v, m in mult.items():
            a = m * mp.e ** (-tau * (m - 1))
            for x in range(64):
                M[x, x ^ v] -= g * a
        return M
    e0 = lambda mult, t, g: min(mp.eigsy(H(mult, t, g), eigvals_only=True))
    tc2 = -mp.log((mp.sqrt(226) - 8) / 18) / 2
    out = {}
    for gs in g_values:
        g = mp.mpf(gs); guess = tc2 - mp.mpf('0.30954679241997696') * g
        f = lambda t: e0(mult_T10, t, g) - e0(mult_P5, t, g)
        out[gs] = mp.findroot(f, (guess - mp.mpf('1e-4'), guess + mp.mpf('1e-4')),
                              solver='anderson', tol=mp.mpf(10) ** (-tol_exp), maxsteps=60)
    return out

def load_ladder_fixture(path):
    d = json.load(open(path))
    return d, {r['g']: r for r in d['rows']}

# ---------- (2) verify_track2b.py: noise-scaled tolerance ----------
def noise_tol(g, kappa=KAPPA):
    return kappa * EPS_ABS / (float(g) ** 2 * abs(C2PRIME_AT_TC2))

def verify_ladder(fixture, strict_hp_agreement=1e-20):
    """Binary64 vs high precision within noise envelope; 40 vs 60 digit agreement."""
    failures, report = [], []
    for r in fixture['rows']:
        g = r['g']
        hp40, hp60, b64 = mp.mpf(r['tau_c_hp40']), mp.mpf(r['tau_c_hp60']), mp.mpf(r['tau_c_binary64_canonical'])
        d = abs(b64 - hp40); tol = noise_tol(g)
        agree = abs(hp40 - hp60)
        report.append((g, float(d), tol, float(agree)))
        if d > tol: failures.append(f"g={g}: |binary64-HP|={float(d):.2e} > tol={tol:.2e}")
        if agree > strict_hp_agreement: failures.append(f"g={g}: 40/60-digit disagreement {float(agree):.1e}")
    return failures, report

def a2_extrap(fixture):
    mp.mp.dps = 40
    tc2 = -mp.log((mp.sqrt(226) - 8) / 18) / 2
    a1 = mp.mpf('-0.30954679241997696')
    a2 = {r['g']: ((mp.mpf(r['tau_c_hp40']) - tc2) / mp.mpf(r['g']) - a1) / mp.mpf(r['g']) for r in fixture['rows']}
    r = lambda g1, g2: 2 * a2[g1] - a2[g2]          # Richardson for g2 = 2*g1
    return {'pair_1e-4_2e-4': r('0.0001', '0.0002'), 'pair_5e-4_1e-3': r('0.0005', '0.001'), 'a2': a2}

# ---------- (3) build_claim_manifest.py: CLM-007..013 ----------
def build_ladder_claims(fixture):
    ex = a2_extrap(fixture)
    rows = {r['g']: r for r in fixture['rows']}
    prov = {k: fixture[k] for k in ('mpmath_version', 'working_precision_digits', 'confirmation_precision_digits', 'solver', 'binary64_solver')}
    claims = {}
    for clm, g in (('CLM-007', '0.002'), ('CLM-008', '0.001'), ('CLM-009', '0.0005'), ('CLM-012', '0.0002'), ('CLM-013', '0.0001')):
        r = rows[g]
        claims[clm] = {'quantity': f'tau_c(g={g})', 'class': 'B', 'status': 'canonical',
                       'reference': 'high_precision', 'value_hp40': r['tau_c_hp40'],
                       'value_binary64': r['tau_c_binary64_canonical'],
                       'noise_envelope': noise_tol(g), 'solver': prov}
    avg = (ex['pair_1e-4_2e-4'] + ex['pair_5e-4_1e-3']) / 2
    claims['CLM-010'] = {'quantity': 'a2_extrap', 'class': 'C', 'status': 'canonical',
                         'value': mp.nstr(avg, 6), 'pairs': {k: mp.nstr(ex[k], 7) for k in ex if k.startswith('pair')},
                         'method': 'linear Richardson 2*a2(g)-a2(2g); empirical, not asymptotic'}
    for k in ('CLM-004', 'CLM-005', 'CLM-006'):       # relabel, values unchanged
        claims[k] = {'class': 'A|D', 'update': 'epistemic label only'}
    return claims

# ---------- (4) generate_hashes.py: deterministic archive + release_info.tex ----------
EXCLUDE = ('.git', 'manuscript', '__pycache__', 'SHA256SUMS', 'release_info.tex')
def deterministic_tar(root):
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode='w', format=tarfile.PAX_FORMAT) as t:
        for dp, dn, fn in os.walk(root):
            dn[:] = sorted(d for d in dn if d not in EXCLUDE)
            for f in sorted(fn):
                if f in EXCLUDE: continue
                p = os.path.join(dp, f); ti = t.gettarinfo(p, os.path.relpath(p, root))
                ti.mtime = 0; ti.uid = ti.gid = 0; ti.uname = ti.gname = ''
                ti.mode = 0o644
                with open(p, 'rb') as fh: t.addfile(ti, fh)
    return buf.getvalue()

def generate_hashes(root, release='CCA-GaugeGraph-v1.9.2-frozen', out_tex=None):
    sums = []
    for dp, dn, fn in os.walk(root):
        dn[:] = sorted(d for d in dn if d not in EXCLUDE)
        for f in sorted(fn):
            if f in EXCLUDE: continue
            p = os.path.join(dp, f)
            sums.append(f"{hashlib.sha256(open(p,'rb').read()).hexdigest()}  {os.path.relpath(p, root)}")
    open(os.path.join(root, 'SHA256SUMS'), 'w').write('\n'.join(sums) + '\n')
    archive = hashlib.sha256(deterministic_tar(root)).hexdigest()
    try: commit = subprocess.check_output(['git', '-C', root, 'rev-parse', 'HEAD'], text=True).strip()
    except Exception: commit = None
    if out_tex and commit:
        h = lambda s: s[:32] + r'\allowbreak ' + s[32:] if len(s) > 32 else s
        open(out_tex, 'w').write(
            f"\\newcommand{{\\relName}}{{{release}}}\n\\newcommand{{\\relCommit}}{{{h(commit)}}}\n"
            f"\\newcommand{{\\relSHA}}{{{h(archive)}}}\n")
    return archive, commit

if __name__ == '__main__':
    fx, _ = load_ladder_fixture(sys.argv[1] if len(sys.argv) > 1 else 'ladder_highprecision.json')
    fails, rep = verify_ladder(fx)
    for g, d, t, a in rep: print(f"g={g:>7}  |b64-HP|={d:.2e}  tol(kappa=20)={t:.2e}  HP40-60={a:.1e}")
    print('FAIL' if fails else 'PASS', fails)
    print(json.dumps(build_ladder_claims(fx), indent=1, default=str)[:1500])
