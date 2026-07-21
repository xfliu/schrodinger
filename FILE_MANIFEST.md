# FILE_MANIFEST — H₂⁺ rigorous-bounds reproduction repository

This manifest maps every **kept** source file to the paper result it produces, and lists the
scratch/diagnostic files that were **excluded** from the reproduction set.

The code follows a **self-contained-per-directory** design: each computation directory carries its
own copy of the box-independent `moments_verified.jl` and the box-specific `assembly_verified.jl` /
Dirichlet modules, with the box constants (LX,LY,LZ) baked into the source and all `include()`
paths relative. A directory can therefore be run as-is without any path configuration. The two
domains are Ω₁ = [−10,10]×[−8,8]² (small box) and Ω₂ = [−20,20]×[−16,16]² (doubled box); nuclei at ±2.

---

## lib/ — interval-eigenvalue dependency
| file | role |
|------|------|
| `lib/Veigs.jl/` | The `Veigs` package: `lehmann_behnke` (rigorous verified eigenvalue enclosure), verified LDLᵀ / inertia / SPD checks. Required by every Stage-2 module. Project.toml + Manifest.toml pin IntervalArithmetic 1.0.8 (CRlibm rounding), KrylovKit, SpecialFunctions. |

## src/stage2/common primitives (present in every stage-2 directory)
| file | module | produces / role |
|------|--------|-----------------|
| `moments_verified.jl` | `MomentsVerified` | **A1** — verified interval Coulomb/kinetic moment enclosures (composite Gauss–Legendre + Bernstein remainder). Box-independent; validated by `validate_A1.jl`. |
| `assembly_verified.jl` | `AssemblyVerified` | **A2** — verified interval Hamiltonian assembly (cosine Neumann basis) via interval `kron`. Box constants LX,LY,LZ set here. |

## src/stage2/omega1_neumann/ — Ω₁ Neumann + Dirichlet-upper
| file | module | produces (paper result) |
|------|--------|--------------------------|
| `pipeline2b.jl` | `Pipeline2b` | Stage-A: certified Neumann Galerkin μ₁,μ₂ + separators (memory-lean interval pipeline). |
| `lg_verified.jl` | `LGVerified` | **Ω₁ Neumann LG lower bound** λ₁ ≥ −0.5520546046 (N=48/N′=64). |
| `dirichlet_verified.jl` | `DirichletVerified` | **Ω₁ Dirichlet upper bound** λ₁ᴰ ≤ −0.5506773767 (N=48). |
| `validate_A1.jl` | — | Runs the A1 moment-enclosure validation (0 containment failures, max width ~1.9e−12). |
| `run2b.jl` | — | Runner → Stage-A separators (Ω₁). |
| `run_lg_shared.jl` | — | Runner → Ω₁ Neumann LG lower (assembles aux N′=64 once, primaries 32/48). |
| `run_dirichlet.jl` | — | Runner → Ω₁ Dirichlet upper bounds. |

## src/stage2/omega1_dirichlet/ — Ω₁ Dirichlet-sector LG (lower bound on λ₁ᴰ)
| file | module | produces (paper result) |
|------|--------|--------------------------|
| `dirichlet_assembly.jl` | `DirichletAssembly` | Sine-basis Dirichlet interval assembly (LX=10). |
| `dirichlet_lg.jl` | `DirichletLG` | **Ω₁ λ₂ᴰ separator** ≥ −0.3312260682 (eoo) and **λ₁ᴰ lower** ≥ −0.5508692 (N=48), ≥ −0.5508177 (N=64). Lean single-D×D-interval-matrix assembler. |
| `run_dir1_all.jl` | — | Runner → Ω₁ Stage A (λ₂ᴰ separator) + Stage B N=48 in one job. |
| `run_dir1_N64.jl` | — | Runner → Ω₁ Stage B N=64/N′=80 (standalone; ~213 GB assembly peak). |

## src/stage2/omega2_neumann/ — Ω₂ Neumann + Dirichlet-upper + μ₂ separator
| file | module | produces (paper result) |
|------|--------|--------------------------|
| `pipeline2b.jl` | `Pipeline2b` | Stage-A (Ω₂, LX=20). |
| `lg_verified.jl` | `LGVerified` | **Ω₂ Neumann LG lower** λ₁ ≥ −0.5520185977 (Stage-A sep), sharpened to ≥ −0.5514436010 with the μ₂ separator. |
| `lg_oee.jl` | `LGOEE` | **Ω₂ certified μ₂ lower** ≥ −0.3395656252 (oee-sector single-vector LG) — the sharper separator ρ. |
| `dirichlet_verified.jl` | `DirichletVerified` | **Ω₂ Dirichlet upper** λ₁ᴰ ≤ −0.5509672618 (N=64). |
| `run_stageA.jl`, `run_stageA_N64.jl` | — | Runners → Ω₂ Stage-A separators (N=32/48, N=64). |
| `run_lg_O2.jl`, `run_lg_O2_N64.jl` | — | Runners → Ω₂ Neumann LG (N=32/48 shared aux; N=64/N′=80). |
| `run_lg_oee.jl`, `run_lg_oee48.jl` | — | Runners → Ω₂ μ₂ separator via oee LG. |
| `run_dirichlet.jl`, `run_dir_N64.jl` | — | Runners → Ω₂ Dirichlet upper (N≤48; N=64). |

## src/stage2/omega2_dirichlet/ — Ω₂ Dirichlet-sector LG (lower bound on λ₁ᴰ)
| file | module | produces (paper result) |
|------|--------|--------------------------|
| `dirichlet_assembly.jl` | `DirichletAssembly` | Sine-basis Dirichlet interval assembly (LX=20). |
| `dirichlet_lg.jl` | `DirichletLG` | **Ω₂ λ₂ᴰ separator** ≥ −0.3382401620 and **λ₁ᴰ lower** ≥ −0.5519149343 (N=48), ≥ −0.5514408191 (N=64). |
| `run_dir_lambda2D.jl` | — | Runner → Ω₂ λ₂ᴰ separator (eoo). |
| `run_dir_lambda1D.jl` | — | Runner → Ω₂ λ₁ᴰ lower (N=48). |
| `run_dir_N64.jl` | — | Runner → Ω₂ λ₁ᴰ lower (N=64/N′=80). |

## src/stage1/ — double-precision reference + moment-method seeds (secondary)
Self-contained project (own Project.toml/Manifest.toml; uses LinearAlgebra + SpecialFunctions +
IterativeSolvers, **not** interval arithmetic). Generates the double-precision LG reference bounds
and the analytic moment-method seeds that the Stage-2 rigorous computation is compared against.
| file | module | role |
|------|--------|------|
| `h2plus_bounds.jl` | `H2plusBounds` | Analytic moment-method H₂⁺ bounds (float). |
| `h2plus_lg_bounds.jl` | `H2plusLG` | Double-precision Lehmann–Goerisch reference (float). |
| `hydrogen_bounds.jl` | `HydrogenBounds` | Hydrogen validation of the method. |
| `lg_cert_driver.jl`, `lg_ref_driver.jl`, `lg_ref2_driver.jl` | — | Drivers → double-precision LG reference (N′=N and enriched N′>N). |
| `aux_eta_ceps.jl`, `aux_eta_driver.jl` | — | Cε auxiliary-quadrature constants for the sharp-η formula. |

---

## Excluded (scratch / diagnostic — NOT part of the reproduction)
Superseded pipeline versions (final path uses `pipeline2b.jl`):
`pipeline2.jl`, `pipeline_verified.jl`, `run2.jl`, `run_certified.jl`, `run_test.jl`

Development validation / diagnostics (one-off checks, not paper results):
`validate_A2.jl`, `validate_A3.jl`, `diag.jl`, `diag_gaps.jl`, `tbig.jl`, `tcomp.jl`,
`test_lb_direct.jl`, `float_pipeline.jl`

Ω₂ scratch: `aux_diag.jl`, `eee_sharp.jl`, `float_diag.jl`, `rho_scan.jl`, `rho_scan2.jl`,
`sep_est.jl`, `val_lean.jl`, `val_oee.jl`, `run_lg_O2.jl`-era diagnostics

Dirichlet-sector scratch: `dir_sep_est.jl`, `sector_diag.jl`, `sector_diag1.jl`, `eoo_spec1.jl`,
`val_B.jl`, `val_dir.jl`, `val_lean_dir.jl`, `val_o1A.jl`, `run_dir1_AB.jl` (superseded by `run_dir1_all.jl`)

Stage-1 scratch: `*.bak_deltar1`, `*.done`, `*.out`, `*.log`, `aux_eps_sweep.jl`,
`track_stage1.html`, `tab_3d_h2plus.tex`
