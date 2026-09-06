# Removing the parity-sector identification from the separator

*Generated 20260906T074328Z. All numbers read from the completed per-run certificate files listed below;
none retyped.*

The separator `rho` fed to the corrected Lehmann–Goerisch stage was, until now, the oee-block
ground-state bound read as a bound on `mu_2(Omega)` through a parity-sector identification that was
assumed rather than certified. Theorem 18 makes the identification unnecessary:

    mu_2(Omega) = min{ lambda_2^eee, lambda_1^oee, lambda_1^eoe }

(`lambda_1^eeo = lambda_1^eoe` exactly, by the `x2 <-> x3` symmetry of `Omega` and `V` — `L_y = L_z`
and the nuclei lie on the `x1`-axis; the two- and three-odd sectors are dominated by a single-odd
sector and drop out, Lemma 17.) So

    rho_assumption_free := min{ rho_recorded_oee, L(lambda_2^eee), L(lambda_1^eoe) }

is an admissible separator with no identification anywhere, provided the two new terms are certified
lower bounds. This note reports them.

## Result

| configuration | N | L(lambda_2^eee) | L(lambda_1^eoe) | rho_recorded (oee) | rho_assumption_free | unchanged |
|---|---|---|---|---|---|---|
| C1 (L=12) | 48 | -0.339381770240 | -0.324950457017 | -0.387173683480 | -0.387173683480 | yes |
| C2 (L=14) | 56 | -0.355282727591 | -0.341268501446 | -0.418692731702 | -0.418692731702 | yes |
| Omega_1 (L=10) | 72 | -0.263131624573 | -0.270653466555 | -0.397309606637 | -0.397309606637 | yes |

**The identification is removed at every configuration, and at zero cost to the separator.** At all
three, both new candidates lie strictly above the recorded oee separator, so the minimum is still
attained by `rho_recorded_oee` and `rho_assumption_free == rho_recorded` to the last bit. No
downstream quantity moves: the Lehmann–Goerisch pencil is fed the same `rho`, so `L_1`, the reported
enclosures and the recorded certificate fields are untouched. What changes is the *status* of `rho` —
from "oee bound plus an assumed sector identification" to "minimum of three certified lower bounds".

Margins of the new candidates above the recorded separator:

| configuration | L(lambda_2^eee) − rho_recorded | L(lambda_1^eoe) − rho_recorded |
|---|---|---|
| C1 (L=12) | 0.047792 | 0.062223 |
| C2 (L=14) | 0.063410 | 0.077424 |
| Omega_1 (L=10) | 0.134178 | 0.126656 |

## Gates

All four gates pass at all three configurations, with numbers, before any of the above was concluded.

**G1 same-space.** The eee run must reproduce the *recorded* `mu_{1,N}` enclosure — this is what
proves the block was built on the same space at the same shift.

| configuration | run enclosure (form-factor driver) | recorded enclosure (Stage A) | contained | max endpoint deviation | recorded width |
|---|---|---|---|---|---|
| C1 | [-0.5511870476181780, -0.5511870476162194] | [-0.5511870476195102, -0.5511870476148756] | yes | 1.344e-12 | 4.635e-12 |
| C2 | [-0.5511231250043536, -0.5511231250020505] | [-0.5511231250069156, -0.5511231249994741] | yes | 2.576e-12 | 7.441e-12 |
| Omega_1 | [-0.5518514996673493, -0.5518514996634128] | [-0.5518514996719721, -0.5518514996588340] | yes | 4.623e-12 | 1.314e-11 |

Containment holds at all three, and the endpoints agree to within the recorded enclosure width in
every case. The two enclosures are *not* the same width: the form-factor run's is narrower by a
factor 2.4 (C1) / 3.2 (C2) / 3.3 (Omega_1), because the two routes bound the interval residual
differently — the form-factor driver uses the row-sum bound on `rad(P)` from its own interval
assembly, Stage A its verified Galerkin residual. Containment of the narrower in the wider, with
endpoint agreement inside the recorded width, is the meaningful reading and is what is reported.

A stronger check was available at Omega_1 and it passes exactly: the eee run reproduces the
pre-existing `formfactor_eee_N72_sig1p7656061569520733.json` — produced by the released two-centre
driver `run_formfactor.jl` — in **28 of 28 numeric fields bit-for-bit** (cost fields excluded),
including `eta'^2 = 0.55`, both `mu_1` endpoints, `mu_2`, `Chat^2` and both Theorem-8 rows. This
re-establishes on the two-centre eee sector the bit-identity regression that the sector-capable
driver was originally gated on.

**G2 direction.** Theory forces `L(lambda_k) <= (that block's k-th Ritz value lower endpoint)`.

| configuration | eee k=2: L / m_2 / slack | eoe k=1: L / m_1 / slack |
|---|---|---|
| C1 | -0.339381770 / -0.205118629 / 0.134263 | -0.324950457 / -0.226702808 / 0.098248 |
| C2 | -0.355282728 / -0.193773566 / 0.161509 | -0.341268501 / -0.219376962 / 0.121892 |
| Omega_1 | -0.263131625 / -0.220368144 / 0.042763 | -0.270653467 / -0.241005547 / 0.029648 |

**G3 validity floor.** Each candidate's `L` must exceed `U_1`, the certified upper bound on
`mu_1(Omega)` (Remark 19) — far weaker than clearing `rho`, and what makes a candidate admissible at
all.

| configuration | U_1 | margin of L(lambda_2^eee) | margin of L(lambda_1^eoe) |
|---|---|---|---|
| C1 | -0.5511870476148756 | 0.211805 | 0.226237 |
| C2 | -0.5511231249994741 | 0.195840 | 0.209855 |
| Omega_1 | -0.5518514996588340 | 0.288720 | 0.281198 |

No candidate fails G3; the resolution is adequate in both sectors at every configuration, and no
larger `N` is needed.

**G4 self-consistency, machine-checked.** For each candidate, fifteen field comparisons were
computed (not asserted) against the configuration's own record: `LX`; `LY = LZ = 0.8 LX`; `sigma`
equal to the configuration's `sigma` *and* to `Ceps^cert` at `eps = 0.4` read from
`cert_<tag>.json stages.A.Ceps['0.4']`; `N` equal to `N_stageA`; sector equal to the candidate's
sector; geometry equal to the two unit charges at `x = ±2.0`; `npanel = round(4.8 LX)`; `nt = 48`;
`t_star = 1.0`; `D` equal to the input file's `D_eee`/`D_eoe`; and `nu_*` computed from that same
`(N, LX)`. Every candidate's `L` is read from a **single** certificate file, so `eta_V^2`, `m_1`
(for `g`), `m_k`, `sigma` and `nu_*` are necessarily from one `(box, N, sector, sigma)`. All
comparisons are `true` at all three configurations; the per-field booleans are in
`assumption_free_separator.json` under `gates.G4_self_consistency`.

In particular `eta_V^2` was certified on each block's own space by that block's own ladder — the eee
value does not transfer to eoe, and, as the measurements below show, does not even transfer between
configurations.

**Symmetry smoke test (eoe vs eeo at N=24).** Since `L_y = L_z`, the eoe and eeo blocks are the same
operator up to relabelling the two transverse axes; a transverse-moment matrix copied from the wrong
axis — a defect a previous driver version actually had, under a comment asserting equal parity —
would show here. Run at the C1 geometry, `D = 2028`: **34 of 54 numeric fields bit-identical, all
others agreeing to 5.56e-15 absolute**, with the certified `eta'^2` rung identical (0.04) and no
non-numeric mismatch. It is *not* bit-exact, and it need not be: the driver accumulates each entry as
`c*Bx*By*Bz` left to right, so exchanging which transverse axis carries the odd-parity block
exchanges two factors of a Float64 product, which is not associative. A wrong-axis moment matrix
would change the block itself and register at the 1e-2 level, not 1e-15.

## Certified form factors measured on each block's own space

| configuration | eee: last failing / certified | eoe: last failing / certified |
|---|---|---|
| C1 | 0.45 / **0.5** | 0.03 / **0.04** |
| C2 | 0.45 / **0.5** | 0.03 / **0.04** |
| Omega_1 | 0.5 / **0.55** | 0.03 / **0.04** |

Two observations that matter for the manuscript:

* the eee form factor is **0.50** at C1 and C2, not the 0.55 certified at Omega_1. 0.55 is a property
  of `V_72^eee` at `L = 10`; each configuration's own space certifies its own value, and here the
  smaller boxes certify a smaller one.
* the eoe form factor is **0.04** at all three configurations, an order of magnitude below the eee
  value. This is consistent with the eoe trial functions being odd in `x2` and therefore vanishing on
  the plane `x2 = 0` that contains both nuclei, where `V` is singular — but that reading is an
  interpretation, not something measured here. What is measured is the bracket `(0.03, 0.04]`.

`nu_*` as used by the driver, `((N+1) pi / (2 L_x))^2`, was checked against the true smallest omitted
mode in each parity class by exact enumeration; it is strictly conservative in all six blocks (e.g.
41.140 used against 42.837 true at C1 eee, minimising triple (50,0,0)). Values per candidate are in
the JSON.

## Predictions versus measurement

The predicted `L` values supplied with the inputs used literature levels (`2sigma_g = -0.2204`,
`pi_u = -0.2144`) as stand-ins for the Ritz values, and `eta_V^2 = 0.55` throughout. The measured
values differ, and the measurement is what is reported:

| configuration | candidate | predicted | measured | measured − predicted |
|---|---|---|---|---|
| C1 | L(lambda_2^eee) | -0.355651 | -0.339382 | +0.016269 |
| C1 | L(lambda_1^eoe) | -0.350627 | -0.324950 | +0.025677 |
| C2 | L(lambda_2^eee) | -0.381047 | -0.355283 | +0.025764 |
| C2 | L(lambda_1^eoe) | -0.376088 | -0.341269 | +0.034819 |
| Omega_1 | L(lambda_2^eee) | -0.263162 | -0.263132 | +0.000030 |
| Omega_1 | L(lambda_1^eoe) | -0.257490 | -0.270653 | -0.013163 |

Five of the six measured values are **higher** (better) than predicted, by 3.0e-5 to 3.5e-2. The one
that is lower is `L(lambda_1^eoe)` at Omega_1, by 1.32e-2: the measured eoe ground Ritz value there
is -0.24101, materially below the `pi_u = -0.2144` stand-in, and that outweighs the gain from the
much smaller measured form factor. It still clears `rho_recorded` by 0.1267, so the conclusion is
unaffected.

Two sources of the discrepancy, both worth stating explicitly:

* the literature stand-ins are not the Ritz values of these blocks, and the gap is
  resolution-dominated rather than box-dominated. The eee second Ritz value measures -0.22037 at
  `L=10, N=72` (mesh `2L/N = 0.28`), -0.20512 at `L=12, N=48` (0.50) and -0.19377 at `L=14, N=56`
  (0.50); the finer Omega_1 space sits closest to the literature value, and the two coarser spaces
  sit above it. Same pattern in eoe: -0.24101, -0.22670, -0.21938 against `pi_u = -0.2144`.
* the prediction's `Chat^2` used the eee `mu_{1,N}` for the eoe candidate as well (the "foreign
  `m_1`" reading), whereas the driver uses each block's own `m_1`. Since `dL/dm_1 > 0`, the foreign
  reading is conservative, and the measured eoe values benefit from the own-block reading.

## Independent recomputation

Each candidate's `L` was recomputed from its three inputs (`eta_V^2`, `m_1`, `m_k`, plus `sigma`, `N`,
`L_x`) in 200-bit outward-rounded interval arithmetic, independently of the driver's own
`IntervalArithmetic` evaluation. Agreement is 1.1e-16 to 7.2e-16 across the six candidates. The
driver's value is the one quoted in every table above; the recomputation is a cross-check, not a
substitute.

## Provenance — which driver and which file for each number

Driver for all eight runs (six production, two smoke): `drivers/run_formfactor_h3_sect.jl`,
md5 `b6a36e79761dc9b09473d42f01cb5989`, the sector-capable form factor driver written for the
three-centre eoe computation. Its whole include chain is md5-identical to the released tree
`jcp2026/drivers/` on `ssh:hpc` (`cert_core.jl` `aab96d3d…`, `moments_verified.jl` `63a936e6…`), as is
the Julia project (`Project.toml` `8361d273…`) and every file of the bundled `lib/Veigs.jl/src`.
Run on `waseda-hpc` with `julia 1.12.7`, `-t 32`, `OPENBLAS_NUM_THREADS=32`,
`numactl --cpunodebind=k --preferred=k`, `FF_FLOATQ=1 FF_NOCACHE=1`, one job per NUMA node.
The Theorem-8 row used is `theorem8_evaluation.keeps_1_minus_eps_0p4` (`fac = 1-eps = 0.6`), which
`certificates/symmetry_free_stageA.json` records as the reading behind the manuscript's numbers.

| number | file | field |
|---|---|---|
| C1 L(lambda_2^eee) | `formfactor_eee_N48_sig1p8146568831188563.json` | `theorem8_evaluation.keeps_1_minus_eps_0p4.L2_lower` |
| C1 L(lambda_1^eoe) | `formfactor_eoe_N48_sig1p8146568831188563.json` | `theorem8_evaluation.keeps_1_minus_eps_0p4.L1_lower` |
| C1 eta_V^2 (eee / eoe) | `formfactor_eee_N48_sig1p8146568831188563.json` / `formfactor_eoe_N48_sig1p8146568831188563.json` | `eta_prime_sq_certified_upper` |
| C1 m_2^eee, m_1^eee | `formfactor_eee_N48_sig1p8146568831188563.json` | `galerkin.mu2_lower`, `galerkin.mu1_lower` |
| C1 m_1^eoe | `formfactor_eoe_N48_sig1p8146568831188563.json` | `galerkin.mu1_lower` |
| C1 sigma, U_1, rho_recorded | `certificates/corrected/cert_C1.json` | `stages.A.Ceps['0.4']`, `stages.A.mu1N_galerkin` sup, `stages.lg.separator_rho` |
| C2 L(lambda_2^eee) | `formfactor_eee_N56_sig1p991749904742849.json` | `theorem8_evaluation.keeps_1_minus_eps_0p4.L2_lower` |
| C2 L(lambda_1^eoe) | `formfactor_eoe_N56_sig1p991749904742849.json` | `theorem8_evaluation.keeps_1_minus_eps_0p4.L1_lower` |
| C2 eta_V^2 (eee / eoe) | `formfactor_eee_N56_sig1p991749904742849.json` / `formfactor_eoe_N56_sig1p991749904742849.json` | `eta_prime_sq_certified_upper` |
| C2 m_2^eee, m_1^eee | `formfactor_eee_N56_sig1p991749904742849.json` | `galerkin.mu2_lower`, `galerkin.mu1_lower` |
| C2 m_1^eoe | `formfactor_eoe_N56_sig1p991749904742849.json` | `galerkin.mu1_lower` |
| C2 sigma, U_1, rho_recorded | `certificates/corrected/cert_C2.json` | `stages.A.Ceps['0.4']`, `stages.A.mu1N_galerkin` sup, `stages.lg.separator_rho` |
| Omega_1 L(lambda_2^eee) | `formfactor_eee_N72_sig1p7656061569520733.json` | `theorem8_evaluation.keeps_1_minus_eps_0p4.L2_lower` |
| Omega_1 L(lambda_1^eoe) | `formfactor_eoe_N72_sig1p7656061569520733.json` | `theorem8_evaluation.keeps_1_minus_eps_0p4.L1_lower` |
| Omega_1 eta_V^2 (eee / eoe) | `formfactor_eee_N72_sig1p7656061569520733.json` / `formfactor_eoe_N72_sig1p7656061569520733.json` | `eta_prime_sq_certified_upper` |
| Omega_1 m_2^eee, m_1^eee | `formfactor_eee_N72_sig1p7656061569520733.json` | `galerkin.mu2_lower`, `galerkin.mu1_lower` |
| Omega_1 m_1^eoe | `formfactor_eoe_N72_sig1p7656061569520733.json` | `galerkin.mu1_lower` |
| Omega_1 sigma, U_1, rho_recorded | `certificates/corrected/cert_SL1.json` | `stages.A.Ceps['0.4']`, `stages.A.mu1N_galerkin` sup, `stages.lg.separator_rho` |

Each run also wrote an immutable timestamped copy alongside its live certificate; both are in
`afsep_runs_20260906T074328Z.tar.gz`, together with every run log and the launch scripts.

## Cost

Wall time and peak RSS below are read from the completed per-run certificate files
(`timings.total`, `peak_rss_gb`), never from a mid-run poll, and all six ran on the same host.

| configuration | block | D | wall (s) | peak RSS (GB) | predicted RSS (GB) |
|---|---|---|---|---|---|
| C1 | eee | 15625 | 838 | 15.66 | 14.6 |
| C1 | eoe | 15000 | 1005 | 14.61 | 13.4 |
| C2 | eee | 24389 | 1574 | 36.62 | 35.5 |
| C2 | eoe | 23548 | 1696 | 34.27 | 33.1 |
| Omega_1 | eee | 50653 | 5390 | 153.66 | 153.2 |
| Omega_1 | eoe | 49284 | 5249 | 143.75 | 145.0 |

Six production runs, 15,753 s of aggregate wall time, plus two N=24 smoke runs. No parameter was
adjusted to make any inequality hold.
