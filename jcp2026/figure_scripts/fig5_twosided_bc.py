"""Regenerate Figure 5 (fig5_twosided_bc): two-sided certified brackets for each boundary
condition, Galerkin tier.  All values are read from the certificate files named below; nothing is
typed in.  Energies are in paper units (= 1/2 hartree).

  panel (a)  Omega_1 (L=10), N=48
      Neumann  mu_1(Omega):  ../certificates/upper_lower_truncation.json  rows[N=48] lg_lower / ritz_upper
      Dirichlet lambda_1^D:  ../certificates/dirichlet_brackets_all_domains.json domains.Omega1.rows.48
  panel (b)  Omega_2 (L=20), N=64
      Neumann  mu_1(Omega):  ../certificates/R3_enclosure_Omega2_N64.json  two_sided_enclosure_N64.lower_certified_best
                             upper = Neumann Ritz value mu_{1,N}, ../certificates/handoff_b2_omega2.json  mu1N
      Dirichlet lambda_1^D:  ../certificates/dirichlet_brackets_all_domains.json domains.Omega2.rows.64
"""
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REF = -0.55131710724745          # exact lambda_1(R^3), paper units, eq. (lam1exact)
UL  = json.load(open("../certificates/upper_lower_truncation.json"))
DB  = json.load(open("../certificates/dirichlet_brackets_all_domains.json"))
R3  = json.load(open("../certificates/R3_enclosure_Omega2_N64.json"))
HB  = json.load(open("../certificates/handoff_b2_omega2.json"))

row48 = next(r for r in UL["rows"] if int(r["N"]) == 48)
neu1 = (row48["lg_lower"], row48["ritz_upper"])
d1   = DB["domains"]["Omega1"]["rows"]["48"]
dir1 = (d1["lambda1D_lower"], d1["lambda1D_upper"])
gap1 = dir1[0] - neu1[1]

neu2 = (R3["two_sided_enclosure_N64"]["lower_certified_best"], HB["mu1N"])
d2   = DB["domains"]["Omega2"]["rows"]["64"]
dir2 = (d2["lambda1D_lower"], d2["lambda1D_upper"])
enc2 = (neu2[0], dir2[1])

BLUE, RED, GREY = "#3b6fb6", "#c0392b", "#888888"

def bracket(ax, x, lo, hi, color, label):
    ax.plot([x, x], [lo, hi], color=color, lw=9, alpha=0.35, solid_capstyle="butt")
    for y in (lo, hi):
        ax.plot([x - 0.14, x + 0.14], [y, y], color=color, lw=2.2)
    ax.text(x, lo - 0.00004, label + f"  width={hi-lo:.1e}", ha="center", va="top", fontsize=7.2, color=color)

fig, (ax, bx) = plt.subplots(1, 2, figsize=(8.7, 2.7))
for a in (ax, bx):
    a.set_xticks([0, 1]); a.set_xticklabels(["Neumann EVP", "Dirichlet EVP"], fontsize=8)
    a.set_xlim(-0.5, 1.5); a.spines["top"].set_visible(False); a.spines["right"].set_visible(False)
    a.axhline(REF, ls=":", color=GREY, lw=1.1); a.tick_params(labelsize=7.5)

bracket(ax, 0, *neu1, BLUE, r"$\mu_1(\Omega)$")
bracket(ax, 1, *dir1, RED,  r"$\lambda_1^D(\Omega)$")
ax.axhspan(neu1[1], dir1[0], color="#dddddd", alpha=0.6, zorder=0)
ax.text(0.5, (neu1[1] + dir1[0]) / 2 + 0.00012, f"truncation gap ({gap1:.1e})", ha="center", va="center", fontsize=7.2)
ax.text(1.42, REF - 0.00004, f"ref {REF:.6f}", ha="right", va="top", fontsize=6.8, color=GREY)
ax.set_ylabel(r"box eigenvalue  (paper units, $=\frac{1}{2}$ Ha)", fontsize=8)
ax.set_title(r"(a) $\Omega_1$, $N=48$: truncation limits the bound", fontsize=8.5, loc="left")

bracket(bx, 0, *neu2, BLUE, r"$\mu_1(\Omega)$")
bracket(bx, 1, *dir2, RED,  r"$\lambda_1^D(\Omega)$")
bx.annotate("", xy=(0.5, enc2[1]), xytext=(0.5, enc2[0]), arrowprops=dict(arrowstyle="<->", lw=1.1))
bx.text(0.56, (enc2[0] + enc2[1]) / 2, "$\\mathbb{R}^3$ enclosure\n" + f"GAP={enc2[1]-enc2[0]:.2e}", fontsize=7.2, va="center")
bx.text(-0.45, REF + 0.00002, f"ref {REF:.6f}", fontsize=6.8, color=GREY)
bx.set_ylabel(r"$\lambda_1(\mathbb{R}^3)$  (paper units, $=\frac{1}{2}$ Ha)", fontsize=8)
bx.set_title(r"(b) $\Omega_2$, $N=64$: both two-sided; enclosure between them", fontsize=8.5, loc="left")

for a in (ax, bx):
    lo_, hi_ = a.get_ylim(); a.set_ylim(lo_ - 0.25*(hi_-lo_), hi_ + 0.05*(hi_-lo_))
fig.tight_layout()
fig.savefig("fig5_twosided_bc.pdf"); fig.savefig("fig5_twosided_bc.png", dpi=200)
print(f"(a) gap={gap1:.4e}  neu={neu1} dir={dir1}")
print(f"(b) enclosure={enc2} width={enc2[1]-enc2[0]:.6e}")
