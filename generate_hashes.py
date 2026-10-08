"""Write SHA256SUMS, commit (if git repo), compute deterministic archive hash, write release_info.tex."""
import sys, os, subprocess
import highprecision_pipeline as hp
root = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
def main(out_tex=None):
    hp.EXCLUDE = ('.git', 'manuscript', '__pycache__', 'SHA256SUMS', 'release_info.tex')
    sums = []
    for dp, dn, fn in os.walk(root):
        dn[:] = sorted(d for d in dn if d not in ('.git', 'manuscript', '__pycache__'))
        for f in sorted(fn):
            if f in ('SHA256SUMS', 'release_info.tex') or f.endswith('.pyc'): continue
            p = os.path.join(dp, f)
            import hashlib; sums.append(f"{hashlib.sha256(open(p, 'rb').read()).hexdigest()}  {os.path.relpath(p, root)}")
    open(os.path.join(root, 'SHA256SUMS'), 'w').write('\n'.join(sums) + '\n')
    subprocess.check_call(['git', '-C', root, 'add', '-A']); 
    subprocess.check_call(['git', '-C', root, '-c', 'user.name=Jake Riley', '-c', 'user.email=245994353+JRiley90@users.noreply.github.com', 'commit', '-q', '-m', 'CCA-GaugeGraph v1.9.2 frozen release'])
    commit = subprocess.check_output(['git', '-C', root, 'rev-parse', 'HEAD'], text=True).strip()
    import hashlib, io, tarfile
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode='w', format=tarfile.PAX_FORMAT) as t:
        for dp, dn, fn in os.walk(root):
            dn[:] = sorted(d for d in dn if d not in ('.git', 'manuscript', '__pycache__'))
            for f in sorted(fn):
                if f == 'release_info.tex' or f.endswith('.pyc'): continue
                p = os.path.join(dp, f); ti = t.gettarinfo(p, os.path.relpath(p, root))
                ti.mtime = 0; ti.uid = ti.gid = 0; ti.uname = ti.gname = ''; ti.mode = 0o644
                with open(p, 'rb') as fh: t.addfile(ti, fh)
    arch = hashlib.sha256(buf.getvalue()).hexdigest()
    if out_tex:
        h = lambda s: s[:32] + r'\allowbreak ' + s[32:] if len(s) > 32 else s
        open(out_tex, 'w').write("\\newcommand{\\relName}{CCA-GaugeGraph-v1.9.2-frozen}\n\\newcommand{\\relCommit}{%s}\n\\newcommand{\\relSHA}{%s}\n" % (h(commit), h(arch)))
    print('commit', commit); print('archive_sha256', arch)
if __name__ == '__main__': main(sys.argv[1] if len(sys.argv) > 1 else None)
