# RESULTS — certified bounds and reference values

> **Superseded, and at the earlier box labelling.** These are the FIRST campaign's results, at the
> Galerkin-tier realization of the Lehmann-Goerisch stage. Here `Omega_1` and `Omega_2` mean the
> boxes of half-length L = 10 and L = 20, which the current manuscript calls **Omega_2 and
> Omega_3**; the manuscript's Omega_1 (L = 8) does not appear below. The manuscript reports
> operator-certified bounds from the corrected realization, which are not comparable with these and
> are found in `jcp2026/certificates/jcp2026/`.

All values are **rigorous interval-arithmetic outputs** (Julia + IntervalArithmetic 1.0.8, CRlibm
directed rounding). Energies in Hartree. Every certificate satisfies the two Lehmann–Goerisch
validity conditions B > 0 and ν < 1 (verified in the JSON `B_positive` / `nu_lt_1` fields).

Domains: **Ω₁ = [−10,10]×[−8,8]²** (small box), **Ω₂ = [−20,20]×[−16,16]²** (doubled box); nuclei at ±2.
`N` = primary basis size per axis, `N′` = auxiliary (enriched) basis size for the LG lower bound.

## Best certified enclosure of the true ground energy λ₁(ℝ³)
Sharpest lower bound (Neumann-LG on Ω₂, sharpened separator) + tightest Dirichlet upper (Ω₂ N=64):

    λ₁(ℝ³) ∈ [ −0.5514436010 , −0.5509672618 ]     width 4.76e−4     (reference ≈ −0.5513)

The lower bound comes from `lg_verified_O2_N64.json` re-Möbius'd with the μ₂ separator
`lg_oee_O2_48_64.json`; the upper from the Ω₂ N=64 Dirichlet Ritz value (`dir_lambda1D_64.json` →
`ritz_primary_float`).

## Neumann Lehmann–Goerisch lower bounds (λ₁ lower)
| domain | N (N′) | LG lower bound | file |
|--------|--------|----------------|------|
| Ω₂ | 32 (64) | −0.7213431719 | lg_verified_O2.json[0] |
| Ω₂ | 48 (64) | −0.5839662573 | lg_verified_O2.json[1] |
| Ω₂ | 64 (80) | −0.5520185977 (Stage-A sep) | lg_verified_O2_N64.json |
| Ω₂ | 64 (80) | **−0.5514436010** (μ₂ sep) | + lg_oee_O2_48_64.json |

μ₂ separator (oee-sector single-vector LG): **μ₂(Ω₂) ≥ −0.3395656252** (`lg_oee_O2_48_64.json`).

## Neumann Stage-A Galerkin eigenvalues + separators
| domain | N | μ₁ᴺ (Ritz upper) | μ₂ᴺ | L₂ (λ₂ lower / separator) | file |
|--------|---|------------------|-----|----------------------------|------|
| Ω₂ | 32 | −0.5485270656 | −0.3311053099 | −0.6549971704 | stageA_dom2.json[0] |
| Ω₂ | 48 | −0.5504842149 | −0.3327902841 | −0.5446951347 | stageA_dom2.json[1] |
| Ω₂ | 64 | −0.5509505957 | −0.3333398828 | −0.4732191311 | stageA_dom2_N64.json |

## Dirichlet eigenvalue λ₁ᴰ — two-sided certified brackets
Lower = single-vector interval LG on the ooo sine sector (separator = certified λ₂ᴰ lower, eoo sector);
upper = ooo Galerkin Rayleigh–Ritz.
| domain | N (N′) | λ₁ᴰ lower (LG) | λ₁ᴰ upper (Ritz) | width | truncation error | files |
|--------|--------|----------------|-------------------|-------|------------------|-------|
| Ω₁ | 48 (64) | −0.5508692004 | −0.5506773767 | 1.918e−4 | [9.8e−5, 7.7e−4] **>0** | dir1_lambda1D_48.json |
| Ω₁ | 64 (80) | −0.5508176726 | −0.5507405805 | 7.71e−5 | [1.50e−4, 7.0e−4] **>0** | dir1_lambda1D_64.json |
| Ω₂ | 48 (64) | −0.5519149343 | −0.5505084907 | 1.406e−3 | [−9.5e−4, 9.4e−4] | dir_lambda1D_48.json |
| Ω₂ | 64 (80) | −0.5514408191 | −0.5509672618 | 4.736e−4 | [−4.7e−4, 4.8e−4] | dir_lambda1D_64.json |

Certified Dirichlet separators (λ₂ᴰ lower, eoo sector):
**λ₂ᴰ(Ω₁) ≥ −0.3312260682** (`dir1_lambda2D_48_64.json`),
**λ₂ᴰ(Ω₂) ≥ −0.3382401620** (`dir_lambda2D_48_64.json`).

**The Dirichlet limitation (key finding).** On the small box Ω₁ both certified Dirichlet lower
bounds lie strictly ABOVE the true λ₁(ℝ³) ≤ −0.5509673, so the Dirichlet eigenvalue **provably
overestimates** the infinite-domain ground energy (truncation error > 0, and its lower bound rises
with N). On the doubled box Ω₂ the true value sits inside the N=64 Dirichlet bracket.

## Precision / separator reference
`certified_precision.json` — Ω₁ certified L₂ (λ₂ lower) at N=32/48/64, used as Neumann separators.
`dirichlet_N64.json` — Ω₂ N=64 Dirichlet Ritz upper raw record.
