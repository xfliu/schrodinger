"""Regenerate Figure 4 (fig_optimal_box): certified enclosure width against truncation box size.

Data sources (both artifacts of the certified pipeline):
  optimal_box_surface.csv  -- the sweep; column `enclosure_width` is the Lehmann-Goerisch
                              tier's FLOAT estimate, not a bound (see the file header).
  cert_*.json              -- CERTIFIED_ENCLOSURE.width for the certified configurations.
Curves are the float sweep at the certified grid; the overlaid stars are certified widths.
Widths and L are in paper units (energy = 1/2 hartree, length = 2 bohr).
"""
import json, numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt

SWEEP = "optimal_box_surface.csv"
CERTS = {"C1": "cert_cmin_C1.json", "C2": "cert_cmatch_C2.json",
         "C3a": "cert_cmin_C3a.json", "C3c": "cert_cmin_C3c.json"}
NGRID = [32, 48, 64]
GREYS = {32: "#b8b8b8", 48: "#6e6e6e", 64: "#1f3f6e"}
STAR, PRIOR = "#7b3294", "#c0392b"

def load():
    sw = pd.read_csv(SWEEP, comment="#")
    sw = sw[sw.tier.astype(str).str.startswith("lg")]
    cw = {}
    for lab, fn in CERTS.items():
        d = json.load(open(fn))
        cw[lab] = dict(width=d["CERTIFIED_ENCLOSURE"]["width"],
                       L=d["box"]["LX"], N=d["N"])
    return sw, cw

def panel(ax, sub, stars, cw, ref_lambda=None, prior=None):
    for N in NGRID:
        g = sub[sub.N == N].sort_values("L").dropna(subset=["enclosure_width"])
        if g.empty:
            continue
        ax.plot(g.L, g.enclosure_width, "-o", color=GREYS[N], ms=4, lw=1.2, zorder=2)
        ax.annotate(f"$N={N}$", (g.L.iloc[-1], g.enclosure_width.iloc[-1]),
                    xytext=(5, 0), textcoords="offset points", color=GREYS[N],
                    fontsize=6, va="center")
    for k, lab in enumerate(stars):
        c = cw[lab]
        ax.plot([c["L"]], [c["width"]], "*", color=STAR, ms=15, zorder=5,
                mec="white", mew=0.6,
                label="certified enclosure" if k == 0 else None)
        ax.annotate(lab, (c["L"], c["width"]), xytext=(9, -1),
                    textcoords="offset points", color=STAR, fontsize=6, va="center")
    if prior is not None:
        ax.plot([prior[0]], [prior[1]], "s", color=PRIOR, ms=7, mfc="none", mew=1.6,
                zorder=4, label="our earlier configuration")
    if ref_lambda is not None:
        thr = 0.05 * abs(ref_lambda)
        ax.axhline(thr, ls=":", color=PRIOR, lw=1.1, zorder=1)
        ax.annotate(f"5% of $|\\lambda_1|$ = {thr:.1e}", (ax.get_xlim()[0], thr),
                    xytext=(4, 4), textcoords="offset points", color=PRIOR, fontsize=6)
    ax.set_yscale("log")
    ax.set_xlabel("box half-length $L$   (paper units, $=2$ bohr)")
    ax.margins(0.06)

def build():
    sw, cw = load()
    fig, (axa, axb) = plt.subplots(1, 2, figsize=(7.2, 3.1))
    a = sw[sw.system == "H2plus"]
    prior = a[(a.L == 20) & (a.N == 64)].enclosure_width
    panel(axa, a, ["C1", "C2"], cw, prior=(20, float(prior.iloc[0])) if len(prior) else None)
    axa.set_ylabel("certified enclosure width\n(paper units, $=\\frac{1}{2}$ Ha)")
    axa.set_title("H$_2^+$: the balance point is interior in $L$", fontsize=8, loc="left")
    axa.legend(frameon=False, fontsize=6, loc="upper left")
    b = sw[(sw.system == "H3plus") & (sw.d_paper == 4.0)]
    panel(axb, b, ["C3a", "C3c"], cw, ref_lambda=float(b.lambda_ref.dropna().iloc[0]))
    axb.set_title("Linear H$_3^{2+}$, $d=4$: certified widths exceed the float estimates",
                  fontsize=8, loc="left")
    for ax, L in ((axa, "a"), (axb, "b")):
        ax.text(-0.16, 1.06, L, transform=ax.transAxes, fontweight="bold", fontsize=10)
    fig.tight_layout()
    return fig, cw

if __name__ == "__main__":
    fig, cw = build()
    for lab, c in cw.items():
        print(f"{lab}: certified width {c['width']:.6e} at L={c['L']:.0f}, N={c['N']:.0f}")
    fig.savefig("fig_optimal_box_regen.pdf")
    fig.savefig("fig_optimal_box_regen.png", dpi=300)
