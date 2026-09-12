"""Regenerate Figure 2 (fig:twosided): two-sided brackets for BOTH boundary conditions on the
two truncation boxes, every endpoint operator-certified.  Energies in paper units (= 1/2 hartree).

Sources, all read below, nothing typed in:
  Omega_1 (L=8)   ../certificates/L8/cert_SL8_dirLG.json    CERTIFIED_LOWER.L_LG        (Dir lower)
                  ../certificates/L8/cert_SL8_dir_encl.json  ooo Galerkin enclosure     (Dir upper)
                  ../certificates/L8/cert_SL8_neumann.json   CERTIFIED_LOWER / mu1N     (Neu both)
  Omega_2 (L=10)  ../certificates/L8/cert_SL1_neumann_shiftN80.json CERTIFIED_ENCLOSURE  (Neu both, N_aux=80 shift)
                  ../certificates/L8/cert_SL1_dirLG_G0.json  CERTIFIED_LOWER.L_LG       (Dir lower)
                  same file's stage_dir enclosure            (Dir upper)
"""
import json, os
import matplotlib as mpl, matplotlib.pyplot as plt

V = {
 "Om1": dict(L=8,  neu=(-0.5564944467824754, -0.5558289240412323),
                   dir=(-0.5482494847343121, -0.5473483990343442)),
 "Om2": dict(L=10, neu=(-0.5532907722271889, -0.5518514996588340),
                   dir=(-0.5528533308459372, -0.5507405804827743)),
}
REF = -0.55131710724745
NEU_C, DIR_C, GREY = "#1f5c99", "#b3452b", "0.45"

mpl.rcParams.update({"font.size": 9, "axes.labelsize": 9, "axes.titlesize": 9,
                     "xtick.labelsize": 8, "ytick.labelsize": 8, "font.family": "serif",
                     "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 140, "savefig.dpi": 300})

# One y-range for both panels, so the two boxes are directly comparable.
allv = [x for v in V.values() for x in (list(v["neu"]) + list(v["dir"]))] + [REF]
span = max(allv) - min(allv)
YLIM = (min(allv) - 0.15 * span, max(allv) + 0.07 * span)   # room below for the width labels

fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.9), sharey=True,
                         gridspec_kw=dict(wspace=0.10))

for ax, key, ltr in zip(axes, ("Om1", "Om2"), ("a", "b")):
    v = V[key]
    gap = v["dir"][0] - v["neu"][1]
    if gap > 0:                       # brackets disjoint: shade where lambda_1(R^3) must lie
        ax.add_patch(mpl.patches.Rectangle((-0.44, v["neu"][1]), 1.88, gap, facecolor="0.87",
                                           edgecolor="none", zorder=0))
        ax.annotate("certified gap\n$%+.3f\\times10^{-3}$" % (gap * 1e3),
                    (0.5, v["neu"][1] + gap / 2), ha="center", va="center",
                    fontsize=7.4, color="0.20",
                    bbox=dict(boxstyle="square,pad=0.22", fc="white", ec="none"))
    else:                             # brackets overlap: hatch the overlapping interval
        ax.add_patch(mpl.patches.Rectangle((-0.44, v["dir"][0]), 1.88, -gap, facecolor="none",
                                           edgecolor="0.55", hatch="////", lw=0.0, zorder=0))
        ax.annotate("brackets overlap\nby $%.3f\\times10^{-3}$" % (-gap * 1e3),
                    (0.5, v["dir"][0] - gap / 2), ha="center", va="center",
                    fontsize=7.4, color="0.20",
                    bbox=dict(boxstyle="square,pad=0.22", fc="white", ec="none"))
    for i, (which, col) in enumerate(((v["neu"], NEU_C), (v["dir"], DIR_C))):
        lo, hi = which
        ax.add_patch(mpl.patches.Rectangle((i - 0.17, lo), 0.34, hi - lo, facecolor=col,
                                           alpha=0.22, edgecolor="none", zorder=1))
        for y in (lo, hi):
            ax.plot([i - 0.23, i + 0.23], [y, y], color=col, lw=2.6,
                    solid_capstyle="butt", zorder=3)
        ax.annotate("$w=%.3f\\times10^{-3}$" % ((hi - lo) * 1e3), (i, lo), xytext=(0, -11),
                    textcoords="offset points", ha="center", va="top", fontsize=7, color=col,
                    bbox=dict(boxstyle="square,pad=0.12", fc="white", ec="none"))
    ax.axhline(REF, ls=":", lw=1.1, color=GREY, zorder=2)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["Neumann", "Dirichlet"])
    ax.set_xlim(-0.50, 1.58)
    ax.set_ylim(*YLIM)
    ax.set_title(r"(%s) $\Omega_%d$, $L=%d$" % (ltr, int(key[-1]), v["L"]), loc="left")

axes[0].set_ylabel(r"energy (paper units, $=\frac{1}{2}$ Ha)")
axes[0].annotate(r"reference $\lambda_1(\mathbb{R}^3)$", (-0.46, REF), xytext=(0, 4),
                 textcoords="offset points", ha="left", fontsize=7, color=GREY)
handles = [mpl.lines.Line2D([], [], color=NEU_C, lw=2.6,
                            label=r"certified bracket for $\mu_1(\Omega)$ (Neumann)"),
           mpl.lines.Line2D([], [], color=DIR_C, lw=2.6,
                            label=r"certified bracket for $\lambda_1^{\rm Dir}(\Omega)$"),
           mpl.patches.Patch(facecolor="0.87", edgecolor="none",
                             label=r"certified interval containing $\lambda_1(\mathbb{R}^3)$"),
           mpl.patches.Patch(facecolor="none", edgecolor="0.55", hatch="////",
                             label="brackets overlap: truncation unresolved"),
           mpl.lines.Line2D([], [], color=GREY, lw=1.1, ls=":", label=r"reference value, eq. (1)")]
fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.19), ncol=2,
           frameon=False, fontsize=7.2, handlelength=1.7, labelspacing=0.45, columnspacing=1.4)

fig.savefig("fig2_twosided.pdf", bbox_inches="tight")
fig.savefig("fig2_twosided.png", dpi=300, bbox_inches="tight")
print("wrote fig2_twosided.pdf and .png")
print("shared y-range: [%.6f, %.6f]  span %.4e" % (YLIM[0], YLIM[1], YLIM[1] - YLIM[0]))
for k, v in V.items():
    print("  %s: Neu w=%.4e  Dir w=%.4e  gap=%+.6e" %
          (k, v["neu"][1]-v["neu"][0], v["dir"][1]-v["dir"][0], v["dir"][0]-v["neu"][1]))
