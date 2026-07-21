# Rigorous two-sided bounds for the H₂⁺ ground-state energy

Reproduction code for the certified interval-arithmetic computation of the H₂⁺ (hydrogen molecular
ion) electronic ground-state energy λ₁. Every bound in the paper is a **rigorous
interval-arithmetic enclosure** — computed with directed rounding so that the printed interval is
mathematically guaranteed to contain the true value.

## What is computed

The physical operator H = −Δ + V (with the two-centre Coulomb potential V = −1/|x−a₁| − 1/|x−a₂|,
nuclei at a = ±2 on the x-axis) lives on ℝ³. We solve it on a finite box Ω and turn the box
computation into a rigorous two-sided enclosure of the true infinite-domain eigenvalue via the
**bound chain**

    L₁_LG  ≤  μ₁ᴺ(Ω)  ≤  λ₁(ℝ³)  ≤  λ₁ᴰ(Ω)  ≤  λ₁,N^D

- **Neumann → lower bound.** Under a confinement condition σ(Ω) = inf_{x∉Ω} V > λ₁, the first
  Neumann eigenvalue μ₁ᴺ(Ω) is ≤ λ₁(ℝ³). The certified **Lehmann–Goerisch (LG)** lower bound
  L₁_LG rigorously bounds μ₁ᴺ from below.
- **Dirichlet → upper bound.** By domain monotonicity λ₁(ℝ³) ≤ λ₁ᴰ(Ω); a sine-Galerkin
  Rayleigh–Ritz value λ₁,N^D is a rigorous upper bound.

Two boxes are used: **Ω₁ = [−10,10]×[−8,8]²** (small) and **Ω₂ = [−20,20]×[−16,16]²** (doubled).
A separate Dirichlet-sector LG computation additionally brackets the Dirichlet eigenvalue λ₁ᴰ
two-sidedly, which rigorously demonstrates the *limitation* of Dirichlet truncation on a small box
(on Ω₁ the certified λ₁ᴰ lower bound provably exceeds the true ground energy).

**Best certified enclosure:** λ₁(ℝ³) ∈ [−0.5514436010, −0.5509672618] Hartree (reference ≈ −0.5513).
The full results table is in [`results/RESULTS.md`](results/RESULTS.md); the method is described in
[`docs/domain_truncation_methods.md`](docs/domain_truncation_methods.md) and the Lehmann–Goerisch
A4-condition validation in [`docs/Validation_of_LG_method.html`](docs/Validation_of_LG_method.html).

## Repository layout

    h2plus-rigorous-bounds/
    ├── Project.toml            top-level Julia environment (bundles Veigs as a path dep)
    ├── README.md               this file
    ├── FILE_MANIFEST.md        every source file → the paper result it produces
    ├── lib/Veigs.jl/           verified-eigenvalue package (lehmann_behnke, verified LDLᵀ/inertia)
    ├── src/
    │   ├── stage2/             RIGOROUS interval computation (the paper's bounds)
    │   │   ├── omega1_neumann/    Ω₁ Neumann Stage-A + LG lower + Dirichlet upper
    │   │   ├── omega1_dirichlet/  Ω₁ Dirichlet-sector LG (λ₂ᴰ separator + λ₁ᴰ lower)
    │   │   ├── omega2_neumann/    Ω₂ Neumann Stage-A + LG lower + μ₂ separator + Dirichlet upper
    │   │   └── omega2_dirichlet/  Ω₂ Dirichlet-sector LG (λ₂ᴰ separator + λ₁ᴰ lower)
    │   └── stage1/             double-precision reference + analytic moment-method seeds (secondary)
    ├── results/                certified result JSONs + RESULTS.md summary table
    └── docs/                   methods note, LG-validation report, limitation figure

Each `src/stage2/<box>_<sector>/` directory is **self-contained**: it carries its own copy of the
box-independent `moments_verified.jl` and the box-specific assembly module, with the box constants
(LX,LY,LZ) baked into the source and all `include()` paths relative. A directory runs as-is with no
path configuration. (This deliberate duplication keeps each domain independently runnable; three of the four
`moments_verified.jl` copies are byte-identical and the fourth — `omega1_neumann/` — differs only
in the `moment_smallt` default `npanel=48` vs `96`, which is never exercised because every call
site passes `npanel=npanel` explicitly. The assembly/dirichlet modules differ only in the
`const LX,LY,LZ` line.)

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

**Ω₂ Stage-A separators** (needed as ρ for the LG step):

    cd src/stage2/omega2_neumann
    julia -t 8 --project=../../.. run_stageA.jl        # N=32,48  → stageA_dom2.json     (~50 min)
    julia -t 8 --project=../../.. run_stageA_N64.jl     # N=64     → stageA_dom2_N64.json  (~4.5 hr, ~90 GB)

**Ω₂ Neumann LG lower** (uses the Stage-A separator):

    julia -t 8 --project=../../.. run_lg_O2.jl          # N=32/48, aux N'=64 → lg_verified_O2.json      (~2.5 hr)
    julia -t 8 --project=../../.. run_lg_O2_N64.jl       # N=64,   aux N'=80 → lg_verified_O2_N64.json    (~8.5 hr, ~213 GB)
    # expect L1_LG(N=64) = -0.5520185977  (Stage-A separator)

**Ω₂ μ₂ separator** (oee-sector LG → sharper ρ that tightens the N=64 lower bound to −0.5514436010):

    julia -t 8 --project=../../.. run_lg_oee.jl          # → lg_oee_O2_48_64.json,  mu2(Omega2) >= -0.3395656252

**Ω₁ Neumann LG lower:**

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

**Ω₂:**

    cd src/stage2/omega2_dirichlet
    julia -t 8 --project=../../.. run_dir_lambda2D.jl    # λ₂ᴰ(Ω₂) ≥ -0.3382401620  → dir_lambda2D_48_64.json  (~2.4 hr)
    julia -t 8 --project=../../.. run_dir_lambda1D.jl    # λ₁ᴰ(Ω₂) N=48 ≥ -0.5519149343 → dir_lambda1D_48.json  (~2.3 hr)
    julia -t 8 --project=../../.. run_dir_N64.jl         # λ₁ᴰ(Ω₂) N=64 ≥ -0.5514408191 → dir_lambda1D_64.json  (~8.5 hr, ~213 GB)

**Ω₁:**

    cd ../omega1_dirichlet
    julia -t 8 --project=../../.. run_dir1_all.jl        # λ₂ᴰ(Ω₁) ≥ -0.3312260682 + λ₁ᴰ(Ω₁) N=48 ≥ -0.5508692004  (~4.8 hr)
    julia -t 8 --project=../../.. run_dir1_N64.jl        # λ₁ᴰ(Ω₁) N=64 ≥ -0.5508176726 → dir1_lambda1D_64.json   (~8.6 hr, ~213 GB)

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

See `LICENSE`. If you use this code, please cite the accompanying paper (H₂⁺ rigorous ground-state
bounds). The `Veigs` package is by X. Liu and Y. Yanagisawa.
