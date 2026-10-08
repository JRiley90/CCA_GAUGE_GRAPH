CCA-GaugeGraph v1.9.3 package
manuscript/  gauge_cycle_multiplicity.tex + fig1-3.png; PDF is a PRE-FREEZE build (Sec. 16 / App. E show "[PRODUCTION BUILD]" until release_info.tex exists). src/         v192_pipeline.py  - high-precision ladder, noise-envelope check (kappa=20), CLM-007..013 claims, deterministic archive hash, release_info.tex writer basis_sampling_pipeline.py - CLM-014 fixture/claim/verify w8_basis_sensitivity.py    - CLM-015 fixture/claim/verify graph_data_v192.py         - embedded copy of Appendix A/B graph data fixtures/    ladder_highprecision.json, basis_sampling_fixtures.json, w8_basis_fixture.json (generated with Python 3.12.3 / NumPy 2.4.4, NOT your frozen 3.13.5 / 2.3.5: regenerate) claims/      claim entries to merge into CLAIM_MANIFEST.json crosschecks/ my original scripts (ladder.py, hp2.py, basis.py, minb.py, mcb.py, w8.py)
Before submission:
Merge src/ into your repo; regenerate fixtures in the frozen environment.
Confirm your pipeline's 40-digit roots match fixtures/ladder_highprecision.json to ~1e-20.
Add CLM-014/015 (and a claim ID for the T10/P5 microscopic energies, or keep the Sec. 8 caveat).
Commit, run generate_hashes (writes release_info.tex), recompile the manuscript.
