# Where each reported number is computed

Every quantity in the paper --- each table row, each figure, each number quoted in the text ---
with the driver that produced it and the file and field in this release that hold it. Paths are
relative to the repository root.

**The `verified` column is not decoration.** Rows marked *verified* were checked by machine: the
value the paper prints was reproduced from the named file and field. Rows marked *not
machine-verified* name the file on the basis of the run record, and the reader should treat the
mapping as a pointer rather than a checked identity. Eight of the twenty-one rows are verified;
the unverified ones are the appendix tables from the earlier campaign, whose values are spread
over several files, together with the schematic figure and the tables of costs, resolutions and
counts that carry no long decimal to probe.

| paper float | quantity reported | driver | file in this release | field | verified |
|---|---|---|---|---|---|
| `eq:headline` | reported enclosure of lambda_1(R^3) on Omega_2 | `run_cert_g.jl SL1` | `certificates/jcp2026/cert_SL1_neumann.json` | `CERTIFIED_ENCLOSURE.lambda1_lower / .lambda1_upper / .width` | verified 2/2 primary values |
| `tab:certified` | Omega_1, Omega_2: width, margin, admissibility pair | `run_cert_g.jl SL8 | SL1` | `certificates/jcp2026/cert_SL8_stagelg.json, cert_SL8_dir_encl.json, cert_SL1_neumann.json` | `Om_1 lower: lambda1_lower_certified; Om_1 upper: mu1_galerkin[1]; Om_2: CERTIFIED_ENCLOSURE` | verified: CERTIFIED_ENCLOSURE.width reproduces every row |
| `tab:certified` | three-centre rows D1, D3a, D3c (not bounds on lambda_1) | `run_cert.jl D1 | D3a | D3c` | `certificates/cert_cmatch_D1.json, cert_cmin_D3a.json, cert_cmin_D3c.json` | `CERTIFIED_ENCLOSURE` | verified: CERTIFIED_ENCLOSURE.width reproduces every row |
| `tab:configs` | box, basis and stage resolutions | `the same certificates` | `certificates/jcp2026/cert_SL8_stageA.json, cert_SL1_stageA.json, cert_SL8_stagemu2.json, cert_SL1_stagemu2.json` | `N; Nprime; N_stageA; Nprime_sep_configured` | not machine-verified |
| `tab:headline-tiers` | Omega_2 certificate: rho, A_2, B, nu_1, admissibility gate | `run_cert_g.jl SL1` | `certificates/jcp2026/cert_SL1_stagelg.json` | `separator_variants.certified_oee_mu2: rho, B, nu; A2; goerisch_variants.0.65` | verified 3/3 primary values |
| `tab:stagea` | certified Stage A on Omega_2: epsilon, eta, sigma, L_k | `run_cert_g.jl SL1` | `certificates/corrected/results_final.json, certificates/corrected/gates_record.json` | `stages.A` | pre-shift campaign; source spans several files, not machine-verified |
| `tab:ceps-cert` | certified coercivity shift C_eps at both boxes | `run_ceps.jl` | `certificates/jcp2026/ceps_certified.json` | `results.N_aux=<N>.eps=<e>.Ceps_cert` | verified 2/2 primary values |
| `tab:formfactor-all` | certified tail-refined form factor eta_V^2 | `run_formfactor_all.jl` | `certificates/sector_free/formfactor_eee_N48_sig1p8146568831188563.json, formfactor_eee_N56_sig1p991749904742849.json, formfactor_eoe_N48_*.json, formfactor_eoe_N56_*.json` | `eta_prime_sq_certified_upper (one file per configuration and block)` | per-configuration files identified; values are bracket ranges, not machine-verified |
| `tab:sepfree` | the three Theorem-18 candidates on Omega_2 and the clearance margin | `run_formfactor.jl eee | eoe` | `certificates/sector_free/assumption_free_separator.csv` | `row cfg=SL1: L_l2_eee; L_l1_eoe; rho_rec; margin_to_rho` | verified 2/2 primary values |
| `tab:dir` | two-sided Dirichlet brackets and the truncation bias | `run_cert_gd.jl` | `certificates/jcp2026/cert_SL8_dirLG.json, cert_SL1_dirLG.json` | `CERTIFIED_LOWER; lambda1D_upper` | verified 1/1 primary values |
| `tab:cost` | per-stage wall time and peak memory of the two certified runs | `the same certificates` | `certificates/jcp2026/ (every cert_SL8_* and cert_SL1_* file)` | `wall_seconds; peak_rss_gb` | not machine-verified |
| `tab:mu3-needed` | mu_1 lower bound with and without the intermediate mu_2 step | `run_cert_g.jl` | `certificates/corrected/gates_record.json, certificates/unconditional_mu2_separator.json` | `separator_variants.<source>.L1_LG` | pre-shift campaign; source spans several files, not machine-verified |
| `tab:audit` | exact/approximate audit of every pipeline step | `cert_selftest.jl` | `exit status` | `-` | not machine-verified |
| `tab:sigmaloc-cert` | the three admissible shifts compared end to end on Omega_2 | `run_cert_g.jl SL1` | `certificates/corrected/results_final.json, certificates/cert_etashift_SL1.json` | `goerisch_variants.<eps_prime>.L1_LG` | pre-shift campaign; source spans several files, not machine-verified |
| `tab:noassume` | cost of removing the sector identification | `run_cert_g.jl SL1` | `certificates/corrected/gates_record.json, certificates/dirichlet_brackets_all_domains.json` | `separator_variants.stageA_L2.L1_LG` | pre-shift campaign; source spans several files, not machine-verified |
| `tab:reduction` | what the parity reduction saves | `run_cert_g.jl / run_cert_gs.jl` | `certificates/jcp2026/cert_SL1_stagelg.json` | `D_aux; wall_seconds` | not machine-verified |
| `tab:symfree` | the certified chain without the parity reduction | `run_cert_gs.jl` | `certificates/symmetry_free_stageA.json` | `rows[*]: eta2, mu2N, L2, peak_gb, wall_s` | verified 2/2 primary values |
| `fig:framework` | schematic; no computed data | `-` | `-` | `-` | not machine-verified |
| `fig:twosided` | brackets for both boundary conditions at matched resolution | `figure_scripts/fig2_twosided.py` | `certificates/jcp2026/cert_SL1_stagelg.json, certificates/dirichlet_brackets_all_domains.json, data/approx_mode_rows.csv` | `named in the generator header` | not machine-verified |
| `fig:optimalbox` | enclosure width against box size | `run_sweep.jl then build_surface.py` | `data/optimal_box_surface.csv` | `enclosure_width` | not machine-verified |
| `fig:convergence` | two convergence diagnostics at N=32,48 | `figure_scripts/fig3_convergence.py` | `figure_scripts/figs.json, certificates/certified_stage2_final.json, certificates/upper_lower_truncation.json` | `fig3.stageA_width -> rows[*].bracket_w; fig3.stageB_width -> rows[*].bracket_width` | not machine-verified |

## One superseded record, retained and labelled

`certificates/jcp2026/cert_SL8_neumann.json` carried a combined enclosure block computed before the
Lehmann--Goerisch lower endpoint at that box was sharpened; its width, 1.033177e-02, contradicts the
9.146048e-03 the paper reports. The block has been renamed `CERTIFIED_ENCLOSURE_SUPERSEDED`, its
values left unaltered, and a note added pointing at the two files that hold the reported endpoints.
The file is retained because it is the only source for the Neumann $\mu_{1,N}$ at that box.

## Naming

`SL8` is the box $\Omega_1$ ($L=8$) and `SL1` is $\Omega_2$ ($L=10$); the tags predate the labels
used in the paper. Files under `certificates/jcp2026/` are the runs reported in the paper. Files
directly under `certificates/` are earlier runs retained for the box-size study, the three-centre
configurations and the appendix tables.

## Reading a certificate

`cert_*_stagelg.json` carries the certified lower endpoint in `lambda1_lower_certified` (and the
enclosing interval in `L1_LG`), the separator and its provenance (`separator_candidates`,
`separator_variants`), the certificate quantities `A2`, `B`, `nu`, the shift variants in
`goerisch_variants`, and `wall_seconds` / `peak_rss_gb`. The combined two-sided enclosure is in
`cert_*_neumann.json` under `CERTIFIED_ENCLOSURE`, with `lower_source` and `upper_source` naming
the files it was assembled from.
