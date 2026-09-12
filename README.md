# Certified two-sided eigenvalue enclosures for Coulomb-singular Schrodinger operators

Reproduction material for certified interval-arithmetic enclosures of the ground-state eigenvalue
of H = -Delta + V on R^3 with attractive Coulomb singularities. Every bound here is a **rigorous
interval-arithmetic enclosure** computed with directed rounding, so the printed interval is
guaranteed to contain the true value.

## Start here

**[`jcp2026/`](jcp2026/) is the release for the current manuscript.** It carries the drivers, one
certificate per reported configuration, the figure generators, and
[`jcp2026/data_provenance.csv`](jcp2026/data_provenance.csv) -- a table mapping every reported
quantity to the file and field that holds it. No reported bound requires combining files. Read
[`jcp2026/README.md`](jcp2026/README.md) first.

The manuscript reports the hydrogen molecular ion H_2^+ at the equilibrium separation and, for
transferability, the linear three-proton ion H_3^(2+). The headline certified enclosure is

    lambda_1(R^3)  in  [ -0.5532907722271889 , -0.5507405804827743 ]      width 2.550e-3

**Units.** The operator carries no factor of one half on the Laplacian, so these energies are in
units of **1/2 hartree**, not hartree. Doubling gives the enclosure in hartree,
[-1.1065815, -1.1014812], which contains the H_2^+ ground state at R = 2 bohr, -1.1026342 hartree.
Lengths are in units of 2 bohr; a paper energy w equals w x 13605.693 meV.

## The bound chain

The operator lives on R^3; it is solved on a finite box Omega and the box computation is turned
into a two-sided enclosure of the whole-space eigenvalue by

    L_1_LG  <=  mu_1^Neu(Omega)  <=  lambda_1(R^3)  <=  lambda_1^Dir(Omega)  <=  U_1^Dir

- **Neumann -> lower bound.** Under the confinement condition sigma(Omega) = inf_{x not in Omega} V
  > lambda_1, the first Neumann eigenvalue of the box lies below lambda_1(R^3); a certified
  Lehmann-Goerisch step bounds it from below, seeded by a separator the method generates for itself.
- **Dirichlet -> upper bound.** By domain monotonicity lambda_1(R^3) <= lambda_1^Dir(Omega), and an
  interval Rayleigh-Ritz value on the Dirichlet form bounds that from above.

Boxes, in the manuscript's labelling (half-sides in paper units, aspect ratio L_y = L_z = 0.8 L_x):

| | half-sides | role | run tag |
|---|---|---|---|
| Omega_1 | [-8,8] x [-6.4,6.4]^2   | second reported box | `SL8` |
| Omega_2 | [-10,10] x [-8,8]^2     | headline box        | `SL1` |
| Omega_3 | [-20,20] x [-16,16]^2   | appendix comparison | `L=20` sweep rows |

## Repository layout

    Project.toml            top-level Julia environment (bundles Veigs as a path dependency)
    README.md               this file
    FILE_MANIFEST.md        source file -> the result it produces (earlier campaign)
    lib/Veigs.jl/           verified-eigenvalue package (lehmann_behnke, verified LDL^T / inertia)
    jcp2026/                *** the release for the current manuscript ***
    src/, results/, docs/   the earlier campaign (see below)

### The earlier campaign: `src/`, `results/`, `docs/`

These directories are the **first** campaign, at the boxes of half-length L = 10 and L = 20. They
are retained because the appendix box-size discussion refers to those runs, and because
`lib/Veigs.jl` is shared. Two warnings for anyone reading them:

- **They predate the manuscript's box labelling.** In `src/` and `results/` the directory and field
  names `omega1` and `omega2` mean the L = 10 and L = 20 boxes, which the manuscript now calls
  **Omega_2 and Omega_3**. The manuscript's Omega_1 (L = 8) has no directory here; it was computed
  with the `jcp2026/` drivers.
- **Their reported enclosure is superseded.** `results/RESULTS.md` records
  [-0.5514436010, -0.5509672618], from the Galerkin-tier realization of the Lehmann-Goerisch stage.
  The manuscript reports operator-certified bounds from the corrected realization; the two are not
  comparable, and the current numbers are the ones in `jcp2026/certificates/jcp2026/`.

The reproduction instructions below apply to that earlier campaign. For the manuscript's runs, use
`jcp2026/README.md` instead.

## Environment setup

**Prerequisites**
- **Julia 1.11 or newer** (developed on 1.12.6). Install via [`juliaup`](https://github.com/JuliaLang/juliaup):
  `curl -fsSL https://install.julialang.org | sh` then `juliaup add 1.12` .
- Rigour depends on **IntervalArithmetic.jl 1.0.8** with CRlibm correctly-rounded directed
  rounding (pinned in `Project.toml` and in `lib/Veigs.jl/Project.toml`).

**Instantiate the environment** (from the repository root):

    julia --project=. -e 'using Pkg; Pkg.instantiate()'

This resolves the top-level `Project.toml`, which pulls in IntervalArithmetic, KrylovKit,
IterativeSolvers, Arpack, SpecialFunctions, JSON, and the bundled **Veigs** package via the
`[sources] Veigs = {path = "lib/Veigs.jl"}` entry — no manual path setup is needed.

> **Julia < 1.11 fallback.** The `[sources]` path-dependency mechanism requires Julia ≥ 1.11. On
> an older Julia, instead run each script with `--project=lib/Veigs.jl` (whose own Project.toml
> supplies IntervalArithmetic/KrylovKit/Arpack) after `] add JSON SpecialFunctions IterativeSolvers`
> into that project.

**Verify the install** by running the moment-enclosure validation (fast, ~1 min):

    cd src/stage2/omega1_neumann
    julia --project=../../.. validate_A1.jl

Expect zero interval-containment failures and max interval width ≈ 1.9e−12.

## Hardware requirements and caveats

The interval matrices are dense: dimension D = ((N+1)/2)³ for the largest (auxiliary N′) basis, and
each `Interval{Float64}` entry is 16 bytes. Peak resident memory during **assembly** at the
largest configuration (N′=80, D=68,921) is **~213 GB**. Guidance:

- The N=48/N′=64 configurations (D≈35,000) peak ~35–50 GB — comfortable on a 128 GB machine.
- The **N=64/N′=80 configurations peak ~213 GB** — require a large-memory node (≥ 251 GB used here).
- **Run only ONE dense N ≥ 64 interval job at a time.** Two such jobs, or an N′=80 job alongside a
  large concurrent process, will exceed 251 GB and be OOM-killed mid-assembly.
- On systems with **`systemd-oomd`** active (userspace, PSI-pressure OOM killer, logs to journald
  not dmesg), a job hitting the memory ceiling dies **silently** — no Julia stacktrace, no dmesg
  entry. If a long run vanishes after assembly with no output JSON, suspect memory pressure.
- Use `julia -t <N>` to enable multithreading for the interval matrix-vector products in the LG step.

## Generating the results

All commands are run from inside the relevant `src/stage2/<box>_<sector>/` directory with
`--project` pointing at the repository root (`../../..`). Each runner has the basis sizes N/N′
baked in (no command-line arguments) and writes a JSON result file into the current directory. The
certified value to expect and the reference file to diff against are given for each.

Times below are wall-clock on a 96-core node (`julia -t 8`); they are dominated by the interval
assembly. Memory peaks are per the caveats above.

### 0. Validate the moment enclosure (prerequisite, ~1 min)

    cd src/stage2/omega1_neumann
    julia --project=../../.. validate_A1.jl
    # → 0 containment failures, max width ≈ 1.9e−12

### 1. Neumann Lehmann–Goerisch lower bounds (λ₁ lower)

**L=20 box, Stage-A separators** (needed as ρ for the LG step):

    cd src/stage2/omega2_neumann
    julia -t 8 --project=../../.. run_stageA.jl        # N=32,48  → stageA_dom2.json     (~50 min)
    julia -t 8 --project=../../.. run_stageA_N64.jl     # N=64     → stageA_dom2_N64.json  (~4.5 hr, ~90 GB)

**L=20 box, Neumann LG lower** (uses the Stage-A separator):

    julia -t 8 --project=../../.. run_lg_O2.jl          # N=32/48, aux N'=64 → lg_verified_O2.json      (~2.5 hr)
    julia -t 8 --project=../../.. run_lg_O2_N64.jl       # N=64,   aux N'=80 → lg_verified_O2_N64.json    (~8.5 hr, ~213 GB)
    # expect L1_LG(N=64) = -0.5520185977  (Stage-A separator)

**L=20 box, mu_2 separator** (oee-sector LG → sharper ρ that tightens the N=64 lower bound to −0.5514436010):

    julia -t 8 --project=../../.. run_lg_oee.jl          # → lg_oee_O2_48_64.json,  mu2(Omega2) >= -0.3395656252

**L=10 box, Neumann LG lower:**

    cd ../omega1_neumann
    julia -t 8 --project=../../.. run2b.jl               # Stage-A separators (Omega1)
    julia -t 8 --project=../../.. run_lg_shared.jl       # LG lower N=32/48, aux N'=64  → L1_LG(N=48) = -0.5520546046

### 2. Dirichlet upper bounds (λ₁ᴰ upper, Rayleigh–Ritz)

    cd src/stage2/omega1_neumann
    julia -t 8 --project=../../.. run_dirichlet.jl       # Omega1 N=32,48 → -0.5506773767 (N=48)
    cd ../omega2_neumann
    julia -t 8 --project=../../.. run_dirichlet.jl       # Omega2 N=32,48 → -0.5505084907 (N=48)
    julia -t 8 --project=../../.. run_dir_N64.jl         # Omega2 N=64    → -0.5509672618

### 3. Dirichlet-sector LG (two-sided bracket on λ₁ᴰ — the Dirichlet-limitation result)

Each run first certifies the λ₂ᴰ separator (eoo sector) then the λ₁ᴰ lower bound (ooo sector);
the ooo Galerkin ground printed alongside is the Dirichlet Ritz upper bound.

**L=20 box:**

    cd src/stage2/omega2_dirichlet
    julia -t 8 --project=../../.. run_dir_lambda2D.jl    # lambda_2^Dir(L=20) ≥ -0.3382401620  → dir_lambda2D_48_64.json  (~2.4 hr)
    julia -t 8 --project=../../.. run_dir_lambda1D.jl    # lambda_1^Dir(L=20) N=48 ≥ -0.5519149343 → dir_lambda1D_48.json  (~2.3 hr)
    julia -t 8 --project=../../.. run_dir_N64.jl         # lambda_1^Dir(L=20) N=64 ≥ -0.5514408191 → dir_lambda1D_64.json  (~8.5 hr, ~213 GB)

**L=10 box:**

    cd ../omega1_dirichlet
    julia -t 8 --project=../../.. run_dir1_all.jl        # lambda_2^Dir(L=10) ≥ -0.3312260682 + lambda_1^Dir(L=10) N=48 ≥ -0.5508692004  (~4.8 hr)
    julia -t 8 --project=../../.. run_dir1_N64.jl        # lambda_1^Dir(L=10) N=64 ≥ -0.5508176726 → dir1_lambda1D_64.json   (~8.6 hr, ~213 GB)

Each finished run prints `Bpos=true nu<1=true` — the two Lehmann–Goerisch validity conditions.
Compare your regenerated JSON against the bundled reference in `results/` (same filename); the
`LG_lower` field should match to all printed digits.

### (optional) Stage-1 double-precision reference

    cd src/stage1
    julia --project=. -e 'using Pkg; Pkg.instantiate()'
    julia --project=. lg_ref2_driver.jl                 # float LG reference (N'>N enriched)

## Reproducing the summary figure

`docs/dirichlet_brackets_all_domains.png` (and its `.json`/`.csv`) collect the four Dirichlet
brackets against the true-λ₁ enclosure — the visual statement of the Dirichlet limitation. The
underlying numbers are exactly the `results/dir*_lambda1D_*.json` files.

## Notes on rigour

- The auxiliary (A4) Lehmann–Goerisch condition is **not** enforced exactly (impossible with the
  Coulomb 1/|x−aᵢ| term in a finite basis). It is validated by the **Goerisch defect identity**:
  A₂ = 2⟨w̃,v⟩ − ⟨w̃,H′w̃⟩ + ⟨r,H′⁻¹r⟩ with the defect energy rigorously enclosed in [0, ‖r‖²/λ_min].
  See `docs/Validation_of_LG_method.html`.
- The cosine (Neumann) basis is L²-orthonormal, so the mass matrix is exactly the identity.
- Interval moments use composite verified Gauss–Legendre with a per-panel Bernstein remainder
  (total remainder ≤ 1e−13 even at N=96); validated by `validate_A1.jl`.

## License and citation

See `LICENSE`. If you use this code, please cite the accompanying paper. The `Veigs` package is by
X. Liu and Y. Yanagisawa.
