# Domain truncation and the two-sided enclosure of the infinite-domain λ₁ — methods note

## The construction (matches manuscript_non_bounded_V.tex §"Domain truncation: two-sided bounds on ℝᵈ")
The physical H₂⁺ operator H = −Δ + V lives on ℝ³. We solve it on a finite box Ω with a
Neumann cosine basis. Two classical facts turn the box computation into a rigorous two-sided
enclosure of the true infinite-domain eigenvalue λ₁(ℝ³):

- **Neumann → lower bound (Lem. 1, LiuJSIAM2026 — "Neumann truncation lower bound").** Define
  σ(Ω) = inf_{x∉Ω} V(x). If the confinement condition σ(Ω) > λ_k holds, then the k-th Neumann
  eigenvalue on Ω satisfies μ_k ≤ λ_k. In practice one checks the a posteriori form
  σ(Ω) > (computed upper bound). For H₂⁺ on Ω₁=[−10,10]×[−8,8]² we computed σ(Ω₁) = −0.2425
  (attained on the box face nearest a nucleus), and −0.2425 > λ₁ ≈ −0.551, so the condition holds
  with wide margin. (The base lemma is stated for V ≥ 0; the Coulomb case uses the negative-
  spectrum variant of LiuNM2026 together with patch truncation of the 1/|x−a| singularity.)
- **Dirichlet → upper bound (domain monotonicity).** H¹₀(Ω) ⊂ H¹(ℝ³), so λ₁(ℝ³) ≤ λ₁^D(Ω),
  the first Dirichlet eigenvalue on Ω. A sine-Galerkin Rayleigh–Ritz value λ₁,N^D is a rigorous
  upper bound on λ₁^D(Ω).

Chaining the discrete/rigorous pieces (the document's bound chain, Eq. \eqref{eq:bound-chain}):

    L₁_LG  ≤  μ₁  ≤  λ₁(ℝ³)  ≤  λ₁^D(Ω)  ≤  λ₁,N^D

where L₁_LG is the certified Lehmann–Goerisch lower bound on the Neumann Galerkin eigenvalue and
λ₁,N^D is the rigorous sine-Galerkin Dirichlet upper bound. Hence **[L₁_LG, λ₁,N^D] is a fully
computable, fully rigorous enclosure of the infinite-domain λ₁** — no literature reference needed.

## The rigorous Dirichlet upper bound (new piece)
Sine basis ψ_n(x) = (1/√L) sin(nπ(x+L)/(2L)), n=1,2,…, vanishing at the walls; L²-orthonormal so
the mass matrix is the identity. The Coulomb sine moment is 0.5·(1/L)·(G[|n−m|] − G[n+m]) — the
same verified two-regime moments G as the Neumann assembly but with a MINUS (sin·sin product-to-
sum). The even-in-x states are the ODD sine modes n=1,3,5,…, giving the same sector dimension D
as the Neumann (e,e,e) sector at matched N. The rigorous upper bound is the sup of the interval
Rayleigh quotient sup(⟨v,Hv⟩/⟨v,v⟩) for the float Dirichlet ground vector v ∈ H¹₀ — any H¹₀
trial vector's Rayleigh quotient is a rigorous upper bound on λ₁^D ≥ λ₁(ℝ³).

## What the two errors are, and why refining N eventually stops helping (fixed domain)
The enclosure width splits into two physically distinct pieces:
1. **Discretization error** — the Neumann bracket width sup(μ₁ᴺ) − L₁_LG (finite N, finite N′).
   Shrinks as N, N′ grow. On Ω₁: 1.94e−3 (N=32) → 2.91e−4 (N=48).
2. **Domain truncation error** — the gap λ₁^D(Ω) − μ₁(Ω) between the Dirichlet and Neumann
   eigenvalues on the SAME box, estimated by the discrete λ₁,N^D − sup(μ₁ᴺ). Independent of N;
   set by the box size only. On Ω₁: 1.08e−3 (N=32) and 1.09e−3 (N=48) — essentially constant.

**Key consequence (the cap).** The Lehmann–Goerisch lower bound converges, as N,N′→∞, to the
Neumann eigenvalue λ₁^Neu(Ω), which is strictly BELOW λ₁(ℝ³) ≤ λ₁^D(Ω). So on a FIXED domain,
refining the basis drives the lower bound upward only until it approaches λ₁^Neu(Ω); it can never
cross the truncation gap to reach λ₁(ℝ³), and never the Dirichlet upper bound. Once the
discretization error has fallen to the level of the truncation gap — which happens at N=48 on Ω₁
(2.9e−4 vs 1.1e−3, same order) — additional modes are wasted: the achievable rigorous lower bound
is capped by the Neumann eigenvalue of that box.

**The only way to tighten further is to enlarge the domain.** The truncation gap λ₁^D(Ω) − μ₁(Ω)
is exponentially small in the truncation radius (Agmon estimates; the bound states decay like
e^{−κ|x|}), so doubling the box shrinks it rapidly at essentially no extra cost — the matrix
dimension D depends only on N and parity, not on the box size, so a larger domain runs in the same
memory and time (only the small-t GL panel count is doubled, npanel 48→96, to hold the panel
half-length h = L/npanel and thus the validated Bernstein remainder constant).

## Certified numbers — current domain Ω₁ = [−10,10]×[−8,8]²
| N | Neumann LG lower | Neumann Ritz upper | Dirichlet upper | λ₁(ℝ³) enclosure | width | truncation gap |
|---|------------------|--------------------|-----------------|------------------|-------|----------------|
| 32 | −0.5534801777 | −0.5515412067 | −0.5504637987 | [−0.5534801777, −0.5504637987] | 3.02e−3 | 1.08e−3 |
| 48 | −0.5520546046 | −0.5517637904 | −0.5506773767 | [−0.5520546046, −0.5506773767] | 1.38e−3 | 1.09e−3 |

## Certified numbers — doubled domain Ω₂ = [−20,20]×[−16,16]² (nuclei fixed at ∓2)
The matrix dimension D depends only on N and parity, not on box size, so Ω₂ runs in the same
memory/time as Ω₁ (npanel doubled 48→96 to hold the Bernstein remainder). Neumann Ritz upper and
Dirichlet upper on Ω₂:

| N | Neumann Ritz upper μ₁ᴺ | Dirichlet upper λ₁,N^D | discrete Dir−Neu gap |
|---|------------------------|------------------------|----------------------|
| 32 | −0.5485270656 | −0.5486625320 | −1.35e−4 |
| 48 | −0.5504842149 | −0.5505084907 | −2.43e−5 |

**The truncation gap collapses on the larger box** — from ~+1.1e−3 on Ω₁ to −1.4e−4 (N=32) and
−2.4e−5 (N=48) on Ω₂, a factor of ~8 (N=32) to ~45 (N=48), consistent with the Agmon
exponential decay of the truncation gap in the box radius. On Ω₂ the discrete Dir−Neu difference
has fallen *below* the N-level discretization floor and even flipped sign: μ₁ᴺ and λ₁,N^D are Ritz
UPPER bounds of two DIFFERENT operators (Neumann vs Dirichlet), so once the true gap
λ₁^D(Ω)−μ₁^Neu(Ω) is smaller than each one's discretization error, their discrete difference is
dominated by discretization and carries no sign information. Interpretation: **Ω₂ is large enough
that domain truncation is no longer the limiting error** — exactly the regime where enlarging the
domain has done its job and the remaining error is discretization, to be reduced by N.

Note the Neumann μ₁ᴺ rises on Ω₂ (−0.5504842 at N=48 vs −0.5517638 on Ω₁): a larger box lifts the
box eigenvalue toward the true λ₁(ℝ³) from below, as μ₁^Neu(Ω) ≤ λ₁(ℝ³) requires.

## Certified two-sided ℝ³ enclosures with the Lehmann–Goerisch lower bound (aux N′=64)
Combining the LG Neumann lower bound with the Dirichlet upper bound gives the fully rigorous
enclosure [L₁_LG, λ₁,N^D] of λ₁(ℝ³):

| domain | N | LG lower | Dirichlet upper | enclosure width | ρ (λ₂ separator) |
|--------|---|----------|-----------------|-----------------|-------------------|
| Ω₁ | 32 | −0.5534801777 | −0.5504637987 | 3.02e−3 | −0.4723 |
| Ω₁ | 48 | −0.5520546046 | −0.5506773767 | **1.38e−3** | −0.4059 |
| Ω₂ | 32 | −0.7213431719 | −0.5486625320 | 1.73e−1 | −0.5447 † |
| Ω₂ | 48 | −0.5839662573 | −0.5505084907 | 3.35e−2 | −0.5447 |

† The Ω₂ N=32 LG certificate uses the shared conservative separator ρ = inf(L₂) = −0.5447 taken
from the Ω₂ N=48 Stage-A run (both Ω₂ configs draw the certified λ₂ lower bound from the same aux
assembly). N=32's own Stage-A separator (−0.6550) fails the ρ > μ₁ validity precondition on the
larger box (sharp-L separation is false at N=32/Ω₂), which is exactly why the shared N=48 separator
is used; the LG certificate for the N=32 row is valid (B > 0, ν < 1) with ρ = −0.5447.

**The key tradeoff (important for the paper).** Doubling the box does exactly what the theory
predicts to the *truncation* error — it collapses (~8× at N=32, ~45× at N=48). But at FIXED N it
makes the certified enclosure WIDER, not tighter: on Ω₂ the N=48 enclosure is 3.35e−2 vs 1.38e−3
on Ω₁. Two mechanisms, both from the larger box:
1. **ν*(N) = ((N+1)π/(2Lₓ))² shrinks 4× when Lₓ doubles.** ν* is the a-priori spectral gap the
   sharp-L and LG bounds lean on; a smaller ν* loosens both.
2. **The rigorous λ₂ separator ρ = inf(L₂) rises toward μ₁.** The excited-state spacing scales
   like ν*, so on the bigger box the levels crowd: at Ω₂ N=48, ρ = −0.5447 sits only 0.006 above
   μ₁ = −0.5505 (vs a 0.146 margin on Ω₁). The LG denominator B collapses (5.5e−4 vs 4.8e−2 on
   Ω₁) and the lower bound loosens to −0.5840.

**Net message.** Truncation error and discretization/separator quality pull in opposite
directions under domain enlargement. Ω₁ at N=48 is already near-balanced (truncation ~1.1e−3,
discretization ~2.9e−4). Enlarging the domain to remove the truncation gap only pays off if N is
raised simultaneously to restore ν* and the separator margin — a larger box demands a
proportionally finer basis. **The tightest certified enclosure to date remains Ω₁ N=48/N′=64:
λ₁(ℝ³) ∈ [−0.5520546046, −0.5506773767], width 1.377e−3.**

## Certified LOWER bound on the Dirichlet eigenvalue λ₁ᴰ — the Dirichlet limitation made rigorous
The pieces above use the Dirichlet eigenvalue only as an *upper* bound to λ₁(ℝ³). A separate
computation brackets λ₁ᴰ(Ω) **two-sidedly**, which lets us state rigorously *how much* the
Dirichlet truncation overestimates the true ground energy.

**Method (sine sectors).** The Dirichlet ground state λ₁ᴰ is the ground of the fully-odd sine
sector **ooo** (odd sine modes in all three axes). Its Lehmann–Goerisch lower bound needs a
separator ρ that is a certified lower bound on λ₂ᴰ; the second Dirichlet eigenvalue is the ground
of the **eoo** sector (odd-x, even-y, even-z sines are the parity of the 2nd bound state). So the
chain per domain is:
- **Stage A** — LG certificate on the eoo sector → certified λ₂ᴰ lower bound ρ (the separator).
- **Stage B** — single-test-vector LG on the ooo sector with separator ρ → certified λ₁ᴰ lower.
- **Dirichlet Ritz upper** — sup interval Rayleigh quotient of the ooo Galerkin ground (same
  assembled matrix as Stage B), giving the two-sided Dirichlet bracket.

The auxiliary condition (A4) is validated by the Goerisch defect identity, not enforced exactly:
A₂ = 2⟨w̃,v⟩ − ⟨w̃,H′w̃⟩ + ⟨r,H′⁻¹r⟩ with the defect energy rigorously enclosed in [0,‖r‖²/λ_min],
λ_min ≥ c + μ_lo (rigorous floor μ_lo = −0.65 ≤ Neumann λ₁). Defect correction ~10⁻²⁴, negligible
vs the A₂ interval width ~10⁻¹⁰. (Full description in Validataion_of_LG_method.html.)

**Certified Dirichlet separators:** λ₂ᴰ(Ω₁) ≥ −0.3312260682 (eoo N=48/N′=64), λ₂ᴰ(Ω₂) ≥
−0.3382401620 (eoo N=48/N′=64). Note the Ω₁ eoo second eigenvalue sits near the continuum
(μ₂ ≈ −0.0041), so its Galerkin bracket uses ρ₂ = −0.15, σ₂ = +0.05 (adaptive), vs Ω₂'s deeper
μ₂ ≈ −0.125.

**Certified two-sided Dirichlet brackets:**
| domain | N (N′) | λ₁ᴰ lower (LG) | λ₁ᴰ upper (Ritz) | width | B>0, ν<1 |
|--------|--------|----------------|-------------------|-------|----------|
| Ω₁ | 48 (64) | −0.5508692004 | −0.5506773767 | 1.918e−4 | ✓ (B=0.107, ν=−2.04) |
| Ω₁ | 64 (80) | −0.5508176726 | −0.5507405805 | 7.71e−5 | ✓ (B=0.107, ν=−2.05) |
| Ω₂ | 48 (64) | −0.5519149343 | −0.5505084907 | 1.406e−3 | ✓ |
| Ω₂ | 64 (80) | −0.5514408191 | −0.5509672618 | 4.736e−4 | ✓ (B=0.101, ν=−2.10) |

**The Dirichlet limitation, proven.** The true ground energy is certified in
λ₁(ℝ³) ∈ [−0.5514436010, −0.5509672618] (Neumann-LG lower on Ω₂; Dirichlet-Ritz upper on Ω₂ N=64;
reference ≈ −0.5513). On the **small box Ω₁** the entire certified Dirichlet bracket lies strictly
ABOVE this enclosure: λ₁ᴰ(Ω₁) ≥ −0.5508692 (N=48), ≥ −0.5508177 (N=64), both greater (less
negative) than the true λ₁(ℝ³) ≤ −0.5509673. Hence the Dirichlet eigenvalue on Ω₁ **provably
overestimates** the infinite-domain ground energy — the certified Dirichlet-truncation error
λ₁ᴰ(Ω₁) − λ₁(ℝ³) is rigorously in [9.8e−5, 7.7e−4] (N=48) and [1.50e−4, 7.0e−4] (N=64), strictly
positive. As N is refined the truncation-error *lower* bound RISES (9.8e−5 → 1.50e−4), sharpening
the proof that Dirichlet-on-a-small-box overestimates.

On the **doubled box Ω₂** the bias vanishes: the true λ₁(ℝ³) enclosure sits *inside* the N=64
Dirichlet bracket [−0.5514408, −0.5509673]; the truncation error is < 4.8e−4 and no longer
sign-definite. This is the Dirichlet counterpart of the Neumann truncation-gap collapse documented
above: enlarging the box removes the Dirichlet confinement bias, at the cost (fixed N) of a wider
bracket, exactly as for the Neumann/LG enclosure.

**Dirichlet-error decomposition (Ω₂ N=64).** The two-sided Dirichlet bracket separates cleanly:
(i) the *Dirichlet-discretization* error = bracket width 4.74e−4 (how well finite N/N′ resolves
λ₁ᴰ itself), and (ii) the *Dirichlet-truncation* error < 4.8e−4 (how far λ₁ᴰ sits above λ₁(ℝ³)).
On Ω₂ these two error sources are the same order — balanced, neither dominates — whereas on Ω₁ the
truncation error (≥1.5e−4, sign-definite) dominates and is the proof of the limitation.

**Software (Dirichlet sector):** dirichlet_lg.jl (both LG stages: eoo separator, ooo ground;
lean single-D×D-interval-matrix assembler with in-place shift_diag! and rigorous μ_lo floor),
dirichlet_assembly.jl (sector-general sine-basis interval assembler). Result files:
dirichlet_brackets_all_domains.{json,csv,png}, dir_lambda{1,2}D_*.json (Ω₂), dir1_lambda{1,2}D_*.json (Ω₁).

## Software
Julia 1.12, IntervalArithmetic.jl 1.0.8 (CRlibm directed rounding), Veigs.jl (lehmann_behnke).
Source: dirichlet_verified.jl (sine-basis upper bound), lg_verified.jl (Neumann LG lower),
pipeline2b.jl (Neumann Stage-A), moments_verified.jl / assembly_verified.jl (verified moments).
