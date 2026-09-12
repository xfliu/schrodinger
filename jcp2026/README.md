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
    data_provenance.csv   every reported quantity -> its driver, file and field (21 rows)
    PROVENANCE.md         the same mapping in prose, with how to read a certificate

## Which file backs which reported quantity

`data_provenance.csv` is the authoritative mapping: one row per paper float or reported number,
naming the driver, the file in this release and the field inside it. Its `verified` column records
which rows were checked by matching the file's contents against the number printed in the paper,
rather than merely checking that the path exists. `PROVENANCE.md` carries the same mapping in prose
together with instructions for reading a certificate. No reported bound requires combining files.

The `cmin` / `cmatch` prefix on a certificate name records how the Goerisch shift was selected.

## Certified versus floating point

Certified quantities are interval-valued and carry both endpoints. Floating-point quantities are
labelled as estimates wherever they appear. In `data/optimal_box_surface.csv` the column
`LG_lower_status` marks each row: the Lehmann-Goerisch surrogate values for the three-centre system
are annotated `DO_NOT_QUOTE`, because they are float estimates that disagree with the certified
values at the 1e-4 level. The column `derived_from_LG` names the columns that inherit this.

## Configuration tags renamed after the runs

The two linear H$_3^{2+}$ configurations at $d=4$ were recorded during the runs as `C3a` and `C3c`
and are named `D3a` and `D3c` in the manuscript and throughout this release, so that the letter
records the system: `C` for H$_2^+$, `D` for the three-centre ion. The rename is textual only --- the
certificate filenames, the `config` field, the driver tag (`julia run_cert.jl D3a`) and the recorded
stage-file names were all updated together, and no numeric value in any certificate was touched. A
reader holding an older copy should read `C3a` as `D3a` and `C3c` as `D3c`.

The two H$_2^+$ boxes carry run tags that do not match the manuscript's symbols, because the
manuscript's box labels were reassigned after the runs. The mapping is:

| manuscript | half-sides (paper units) | run tag in this release |
|---|---|---|
| $\Omega_1$ | $[-8,8]\times[-6.4,6.4]^2$  | `SL8` |
| $\Omega_2$ | $[-10,10]\times[-8,8]^2$    | `SL1` |
| $\Omega_3$ (appendix only) | $[-20,20]\times[-16,16]^2$ | the `L=20` sweep rows |

`SL1` was the first box studied, hence the tag; it is the manuscript's $\Omega_2$ and carries the
headline enclosure. No numeric value depends on the labelling.

## A note on one field name

`certified_optimum.json` carries a field named `published_baseline`. The name is a misnomer kept for
provenance: the L=20, N=64 configuration it records is an earlier run of this same code by the same
authors, not a published result. The manuscript describes it as such throughout.

## Two tiers of certificate

`certificates/*.json` are the **Galerkin-tier** certificates: the Lehmann-Goerisch stage is realized on
the auxiliary Galerkin space, so they certify the auxiliary Galerkin eigenvalue, not lambda_1 (Section 6.1
of the paper). `certificates/corrected/` holds **operator-certified** runs produced
by `drivers/run_cert_g.jl` with the corrected Goerisch block `drivers/lg_goerisch_correct2.jl`. The runs
reported in the current manuscript are in `certificates/jcp2026/` (see below); the files directly under
`corrected/` are the earlier campaign at the configurations `C1` and `C2`, retained because the
box-size discussion in the appendix refers to them. `corrected/gate/`
holds the G1 reduction-gate runs (discretization term forced to zero; they reproduce the Galerkin-tier
certificates bit-for-bit), `corrected/gates_record.json` the four gate verdicts, and
`corrected/results_final.json` the per-configuration summary. `ceps_certified_*_N80.json` are the
admissibility-gate constants C_eps'. Path strings inside these files were reduced to bare filenames.

### `certificates/jcp2026/` --- the runs the manuscript reports

This directory holds the current campaign and is what `data_provenance.csv` points at. Per box
(`SL8` = $\Omega_1$, `SL1` = $\Omega_2$) there is one certificate per stage:
`*_stageA.json` the projection stage, `*_stagemu2.json` the separator stage,
`*_stagelg.json` the Lehmann--Goerisch stage carrying the reported lower endpoint and both
separator variants, `*_dirLG.json` and `*_neumann.json` the Dirichlet and Neumann brackets, and
`ceps_certified_*_N80.json` the admissibility constant. The `formfactor_*` files carry the
tail-refined form factor at each shift used. One file is superseded in part: the combined-enclosure
block of `cert_SL8_neumann.json` predates a sharpening and its width disagrees with the manuscript;
the block is marked as superseded in place and the current endpoints are in `cert_SL8_stagelg.json`
and `cert_SL8_dirLG.json`.

## The separator without a sector identification

`certificates/sector_free/` holds the two further certified sector bounds required by the
minimum-of-three theorem: `formfactor_eee_*` carries the eee block's second Ritz value and its
projection bound at k=2, `formfactor_eoe_*` the eoe block's ground bound at k=1, each with the form
factor re-certified on that block. `assumption_free_separator.csv` collects the candidates and shows
that rho = min is the oee value, so the separator carries no sector identification -- for three
configurations: `SL1` (the manuscript's $\Omega_2$) and the two departed configurations `C1` and `C2`.

**The coverage is not uniform across the two boxes reported in the manuscript, and the manuscript
says so.** On $\Omega_2$ all three candidates are certified. On $\Omega_1$ the eee candidate is
certified -- `certificates/jcp2026/formfactor_eee_N72_sig1p4701145303271432.json`, the run at that
box's shift -- and it clears the oee value, but **the eoe candidate was never certified there**; the
manuscript reports its approximate-mode margin and prices the unconditional fallback that removes
the dependence. There is therefore no $\Omega_1$ row in `assumption_free_separator.csv`, and its
absence is a fact about the computation, not an omission from the release.

The N=24 eoe/eeo pair is the transverse-symmetry smoke test.

## Running the drivers

The drivers live in `drivers/` and are run from that directory with the repository root as the
Julia project. The verified-eigenvalue package `Veigs` (Lehmann-Behnke, verified Cholesky, inertia)
is bundled at `lib/Veigs.jl` in the repository root and is wired in as a path dependency.

    # once, from the repository root (Julia >= 1.11)
    julia --project=. -e 'using Pkg; Pkg.instantiate()'

    # the hard gate: bit-for-bit thread-exactness of the interval assembly (exit 0 = pass)
    cd jcp2026/drivers
    VEIGS_SRC=../../lib/Veigs.jl/src julia --project=../.. -t 4 cert_selftest.jl

    # a certified configuration; tags: SL8 and SL1 are the manuscript's Omega_1 and Omega_2,
    # D1/D3a/D3c the three-centre ion, C1/C2 the departed configurations of the appendix
    VEIGS_SRC=../../lib/Veigs.jl/src julia --project=../.. -t <threads> run_cert.jl SL1
    # operator-certified (corrected realization); CERT_G_DISC=0 reproduces the G1 gate run
    VEIGS_SRC=../../lib/Veigs.jl/src julia --project=../.. -t <threads> run_cert_g.jl SL1

Every `include` in the drivers resolves inside `drivers/`; the chain is
`run_cert.jl -> cert_core.jl -> moments_verified.jl`, `lg_oee.jl / pipeline2.jl -> assembly_verified.jl`,
`dirichlet_lg.jl -> dirichlet_assembly.jl`, `run_sweep.jl / run_checks.jl -> sweep_core.jl -> float_core.jl`.
The production runs need large memory (the certificates record `peak_rss_gb` per stage); the self-test
runs in seconds on a laptop.

## Reproducing the figures

    cd figure_scripts
    python fig_optimal_box.py     # box-size figure, reads the sweep and four certificates
    python fig_cost.py            # cost figure, reads timings and widths from the certificates
    python fig3_convergence.py    # convergence figure, reads figs.json
    python make_fig1_framework.py # schematic, no data

No figure script hardcodes a certified result; each reads the certificates from `../certificates/`
and the sweep from `../data/`. The two-sided boundary-condition schematic in the paper is a drawn diagram and has no
generator.

## Note on file paths inside the certificates

The certificates record, for each imported constant, the file it was read from. Those fields were
rewritten to bare filenames for this release; the directory components referred to scratch
locations on the machine the runs were performed on and carry no information. All numeric content
is byte-identical to the certificates as produced.
