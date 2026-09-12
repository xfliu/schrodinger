"""Figure 2: two-sided brackets for both boundary conditions at L = 10 and L = 14.

Panel (a)  L = 10 (Omega_1) -- EVERY ENDPOINT CERTIFIED
    Neumann   lower  cert_SL1.json        /CERTIFIED_ENCLOSURE/lambda1_lower   (corrected LG, N=64/N'=80)
              upper  cert_SL1.json        /stages/A/mu1N_galerkin[1]           (interval Rayleigh-Ritz, N_stageA)
    Dirichlet lower  cert_SL1_dirLG.json  /CERTIFIED_LOWER/lower_certified     (corrected LG, N=48/N'=64)
              upper  cert_SL1.json        /CERTIFIED_ENCLOSURE/lambda1_upper   (interval Rayleigh-Ritz, N=64)

Panel (b)  L = 14 (C2) -- APPROXIMATE MODE, NOT BOUNDS
    one Float64 matrix-free run at N=64/N'=80: LG_lower, mu1N_Ritz, lambda1D_Ritz.
    The approximate mode produces no Dirichlet lower endpoint; that slot is marked absent
    rather than filled from another tier.

Every value is read from the files named above; nothing is typed in.  Energies are in paper
units (= 1/2 hartree).
"""
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt

REF = -0.55131710724745                      # exact lambda_1(R^3), paper units

# ------------------------------------------------------------------ inputs ---
SL1 = json.load(open("inputs/cert_SL1.json"))
DIR = json.load(open("inputs/cert_SL1_dirLG.json"))
APX = json.load(open("inputs/fig2_inputs.json"))["approx_L14_row"]

cert = dict(
    neu_lo=SL1["CERTIFIED_ENCLOSURE"]["lambda1_lower"],
    neu_up=SL1["stages"]["A"]["mu1N_galerkin"][1],
    dir_lo=DIR["CERTIFIED_LOWER"]["lower_certified"],
    dir_up=SL1["CERTIFIED_ENCLOSURE"]["lambda1_upper"],
    N=SL1["N"], Np=SL1["Nprime"], N_stageA=SL1["stages"]["A"]["N_stageA"],
    dir_N=DIR["stages"]["dirlg"]["N"], dir_Np=DIR["stages"]["dirlg"]["Nprime"],
    L=SL1["box"]["LX"], dir_L=DIR["box"]["LX"],
)
apx = dict(
    neu_lo=float(APX["LG_lower"]), neu_up=float(APX["mu1N_Ritz"]),
    dir_up=float(APX["lambda1D_Ritz"]),
    N=int(APX["N"]), Np=int(APX["Nprime"]), L=float(APX["L"]),
)

# sanity: the chain the figure draws must hold in the data
assert cert["L"] == cert["dir_L"], "the two certificates are not on the same box"
assert cert["neu_lo"] < cert["neu_up"] < cert["dir_up"], "certified Neumann bracket misordered"
assert cert["dir_lo"] < cert["dir_up"], "certified Dirichlet bracket misordered"
assert cert["neu_lo"] <= REF <= cert["dir_up"], "certified enclosure does not contain the reference"
assert apx["neu_lo"] <= REF <= apx["dir_up"], "approximate enclosure does not contain the reference"

# ------------------------------------------------------------------ style ----
NEU, DIR_C, GREY, ENC = "#1f6f8b", "#b5651d", "#8a8a8a", "#dfe7ee"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.linewidth": .8,
                     "xtick.direction": "out", "ytick.direction": "out",
                     "xtick.major.width": .8, "ytick.major.width": .8})

fig, axes = plt.subplots(1, 2, figsize=(7.1, 3.8), gridspec_kw=dict(wspace=0.42))


def bracket(ax, i, lo, up, colour, dashed=False, label=None):
    """Two-sided bracket: caps at both endpoints, thin stem between."""
    ls = (0, (3.2, 1.7)) if dashed else "-"
    ax.plot([i, i], [lo, up], color=colour, lw=1.1, ls=ls, zorder=2)
    for y in (lo, up):
        ax.plot([i - 0.20, i + 0.20], [y, y], color=colour, lw=2.7,
                ls="-" if not dashed else (0, (2.2, 1.3)), solid_capstyle="butt", zorder=3)
    if label:
        ax.annotate(label, (i + 0.26, 0.5 * (lo + up)), ha="left", va="center",
                    fontsize=7, color=colour)


# ------------------------------------------------- panel (a): all certified --
axA = axes[0]
axA.axhspan(cert["neu_lo"], cert["dir_up"], color=ENC, zorder=0)
bracket(axA, 0, cert["neu_lo"], cert["neu_up"], NEU,
        label=r"$w=%.2f\times10^{-3}$" % ((cert["neu_up"] - cert["neu_lo"]) * 1e3))
bracket(axA, 1, cert["dir_lo"], cert["dir_up"], DIR_C,
        label=r"$w=%.2f\times10^{-3}$" % ((cert["dir_up"] - cert["dir_lo"]) * 1e3))
encA = cert["dir_up"] - cert["neu_lo"]
axA.annotate(r"$\lambda_1(\mathbb{R}^3)$ enclosure: $w=%.2f\times10^{-3}$" "\n" r"(%.2f%% of $|\lambda_1|$)"
             % (encA * 1e3, encA / abs(REF) * 100),
             (2.03, cert["neu_lo"]), xytext=(-2, 4), textcoords="offset points",
             ha="right", va="bottom", fontsize=7, color="#2f4f63")
axA.set_xticks([0, 1])
axA.set_xticklabels(["Neumann EVP\n" r"$\mu_1(\Omega)$",
                     "Dirichlet EVP\n" r"$\lambda_1^{\rm Dir}(\Omega)$"])
axA.annotate(r"$N{=}%d/N'{=}%d$, except $\mu_{1,N}$ at $N{=}%d$" "\n"
             r"and the Dirichlet lower at $N{=}%d/N'{=}%d$"
             % (cert["N"], cert["Np"], cert["N_stageA"], cert["dir_N"], cert["dir_Np"]),
             (-0.48, cert["dir_up"]), xytext=(0, 6), textcoords="offset points",
             ha="left", va="bottom", fontsize=6.6, color="0.35")
axA.set_title("(a) $L=10$: every endpoint certified", loc="left", fontsize=8.5)
axA.set_ylabel(r"energy   (paper units, $=\frac{1}{2}$ Ha)", fontsize=8)

# ------------------------------------------------ panel (b): approximate -----
axB = axes[1]
axB.axhspan(apx["neu_lo"], apx["dir_up"], color=ENC, zorder=0)
bracket(axB, 0, apx["neu_lo"], apx["neu_up"], NEU, dashed=True,
        label=r"$w=%.2f\times10^{-4}$" % ((apx["neu_up"] - apx["neu_lo"]) * 1e4))
axB.plot([1 - 0.20, 1 + 0.20], [apx["dir_up"]] * 2, color=DIR_C, lw=2.7,
         ls=(0, (2.2, 1.3)), solid_capstyle="butt", zorder=3)
encB = apx["dir_up"] - apx["neu_lo"]
axB.annotate(r"$\lambda_1(\mathbb{R}^3)$ enclosure: $w=%.2f\times10^{-4}$" "\n" r"(%.3f%% of $|\lambda_1|$)"
             % (encB * 1e4, encB / abs(REF) * 100),
             (2.03, apx["neu_lo"]), xytext=(-2, 4), textcoords="offset points",
             ha="right", va="bottom", fontsize=7, color="#2f4f63")
axB.set_xticks([0, 1])
axB.set_xticklabels(["Neumann EVP\n" r"$\mu_1(\Omega)$",
                     "Dirichlet EVP\n" r"$\lambda_1^{\rm Dir}(\Omega)$"])
axB.annotate(r"$N{=}%d/N'{=}%d$ throughout, one run" % (apx["N"], apx["Np"]),
             (-0.48, apx["dir_up"]), xytext=(0, 6), textcoords="offset points",
             ha="left", va="bottom", fontsize=6.6, color="0.35")
axB.set_title("(b) $L=14$: approximate mode, not bounds", loc="left", fontsize=8.5)

# the absent Dirichlet lower is marked, not left as an empty slot
ylo_b = apx["neu_lo"] - 0.34 * encB
axB.plot([1 - 0.20, 1 + 0.20], [ylo_b] * 2, color=DIR_C, lw=1.6, ls=(0, (1, 1.7)), zorder=2)
axB.annotate("lower endpoint not produced\nin the approximate mode", (1 + 0.26, ylo_b),
             ha="left", va="center", fontsize=6.6, color=DIR_C, style="italic")

# ------------------------------------------------------------ both panels ----
for ax, lo, hi, fpad in ((axA, cert["neu_lo"], cert["dir_up"], 0.22),
                         (axB, ylo_b, apx["dir_up"], 0.12)):
    ax.axhline(REF, ls=":", lw=1.0, color=GREY, zorder=1)
    pad = fpad * (hi - lo)
    ax.set_xlim(-0.52, 2.05)
    ax.set_ylim(lo - pad, hi + pad)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=7)
    ax.ticklabel_format(axis="y", useOffset=False, style="plain")

axA.annotate(r"reference $\lambda_1$", (-0.48, REF), xytext=(0, 3), textcoords="offset points",
             ha="left", fontsize=6.8, color=GREY)

handles = [
    mpl.lines.Line2D([], [], color="0.35", lw=2.7, ls="-",
                     label="certified endpoint (interval arithmetic)"),
    mpl.lines.Line2D([], [], color="0.35", lw=2.7, ls=(0, (2.2, 1.3)),
                     label="approximate mode: Float64, not a bound"),
    mpl.patches.Patch(facecolor=ENC, edgecolor="none",
                      label=r"enclosure of $\lambda_1(\mathbb{R}^3)$: Neumann lower to Dirichlet upper"),
]
fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.155), frameon=False,
           fontsize=7, handlelength=2.1, labelspacing=0.45)

fig.savefig("fig2_L10_L14.pdf", bbox_inches="tight")
fig.savefig("fig2_L10_L14.png", dpi=300, bbox_inches="tight")
print("certified  (a) Neumann [%.10f, %.10f] w=%.4e" % (cert["neu_lo"], cert["neu_up"], cert["neu_up"] - cert["neu_lo"]))
print("certified  (a) Dirichlet [%.10f, %.10f] w=%.4e" % (cert["dir_lo"], cert["dir_up"], cert["dir_up"] - cert["dir_lo"]))
print("certified  (a) enclosure w=%.6e = %.3f%%" % (encA, encA / abs(REF) * 100))
print("approx     (b) Neumann [%.10f, %.10f] w=%.4e" % (apx["neu_lo"], apx["neu_up"], apx["neu_up"] - apx["neu_lo"]))
print("approx     (b) Dirichlet upper %.10f ; no lower in this mode" % apx["dir_up"])
print("approx     (b) enclosure w=%.6e = %.4f%%" % (encB, encB / abs(REF) * 100))
print("gap (a) dir_lo - neu_up = %+.4e  -> %s" % (cert["dir_lo"] - cert["neu_up"],
      "disjoint" if cert["dir_lo"] > cert["neu_up"] else "OVERLAP, no band drawn"))
print("gap (b) dir_up - neu_up = %+.4e" % (apx["dir_up"] - apx["neu_up"]))
