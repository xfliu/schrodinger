# Certified enclosures for Coulomb-singular Schrodinger operators (JCP submission)

Release material for the manuscript *Certified two-sided enclosures for Coulomb-singular
Schrodinger operators on R^3*. Everything reported in the paper is produced by the drivers here
and recorded in the certificate files here.

## Units

The operator is H = -Laplacian + V, with no factor of one half on the Laplacian. Consequently

* energies in these files are in units of 1/2 hartree,
* lengths are in units of 2 bohr.

A paper energy w corresponds to w x 13605.693 meV. Every geometry in these files is in paper
units; converting a bond length to bohr means multiplying by two.

## Layout

    drivers/          the verified pipeline and the sweep tooling
    certificates/     one JSON per certified configuration, plus the auxiliary certificates
    figure_scripts/   figure generators, with the data files they read
    data/             the box-size x basis-size surface and the balance-law fits

## Which file backs which reported quantity

| reported quantity | driver | output |
|---|---|---|
| Certified results table, C1/C2/D1/C3a/C3c | `run_cert.jl <tag>` | `certificates/cert_cmin_C1.json`, `cert_cmatch_C2.json`, `cert_cmatch_D1.json`, `cert_cmin_C3a.json`, `cert_cmin_C3c.json`, field `CERTIFIED_ENCLOSURE` |
| the headline enclosure and the factor 1.96 | `run_cert.jl C2` | `certificates/cert_cmatch_C2.json` |
| Stage-A table (eps, eta, sigma, L_k) | `pipeline2.jl` | `certificates/certified_precision.json` |
| oee separator | `lg_oee.jl` | `certificates/lg_oee_certified_mu2_O2.json` |
| eoo/ooo Dirichlet rows | `dirichlet_lg.jl` | `certificates/dir_lambda2D_48_64.json` |
| (L,N) surface | `run_sweep.jl` | dense and LG tier CSVs, merged below |
| balance-law fits | `build_surface.py` | `data/optimal_box_surface.csv`, `data/balance_law.json` |
| per-stage cost | `run_cert.jl` | the same certificates, `stages.*.wall_seconds` and `peak_rss_gb` |
| moment-enclosure validation | `run_checks.jl` | `certificates/step3_checks.json` |
| assembly bit-for-bit gate | `cert_selftest.jl` | exit status |
| projection form factor, all sectors | `run_formfactor.jl`, `run_formfactor_h3_sect.jl` | per-run form-factor certificates |
| assumption-free second-eigenvalue separator | `run_formfactor_h3_sect.jl` | `certificates/unconditional_mu2_separator.json` |
| chain run without the parity reduction | `run_formfactor.jl` | `certificates/symmetry_free_stageA.json` |

The `cmin` / `cmatch` prefix on a certificate name records how the Goerisch shift was selected.

## Certified versus floating point

Certified quantities are interval-valued and carry both endpoints. Floating-point quantities are
labelled as estimates wherever they appear. In `data/optimal_box_surface.csv` the column
`LG_lower_status` marks each row: the Lehmann-Goerisch surrogate values for the three-centre system
are annotated `DO_NOT_QUOTE`, because they are float estimates that disagree with the certified
values at the 1e-4 level. The column `derived_from_LG` names the columns that inherit this.

## Reproducing the figures

    cd figure_scripts
    python fig_optimal_box.py     # box-size figure, reads the sweep and four certificates
    python fig_cost.py            # cost figure, reads timings and widths from the certificates
    python fig3_convergence.py    # convergence figure, reads figs.json
    python make_fig1_framework.py # schematic, no data

No figure script hardcodes a certified result; each reads the certificates by filename from its own
directory. The two-sided boundary-condition schematic in the paper is a drawn diagram and has no
generator.

## Note on file paths inside the certificates

The certificates record, for each imported constant, the file it was read from. Those fields were
rewritten to bare filenames for this release; the directory components referred to scratch
locations on the machine the runs were performed on and carry no information. All numeric content
is byte-identical to the certificates as produced.
