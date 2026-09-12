"""Regenerate Figure 4 (fig_optimal_box): Galerkin-tier enclosure width against truncation box size.

Data sources (both artifacts of the certified pipeline):
  optimal_box_surface.csv  -- the sweep; column `enclosure_width` is the Lehmann-Goerisch
                              tier's FLOAT estimate, not a bound (see the file header).
  cert_*.json              -- CERTIFIED_ENCLOSURE.width for the certified configurations.
Curves are the float sweep at the certified grid; the overlaid stars are certified widths.
Widths and L are in paper units (energy = 1/2 hartree, length = 2 bohr).
"""
import json, numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt


def apply_figure_style(*, frame="open", font=None, sizes=(8, 7, 6), grid=False):
    import matplotlib as mpl
    if frame not in ("open", "boxed", "none"):
        raise ValueError(f"frame must be 'open'|'boxed'|'none', got {frame!r}")

    try:
        import os, sys, glob, matplotlib.font_manager as fm
        fdir = os.path.join(os.environ.get("CONDA_PREFIX") or sys.prefix, "fonts")
        if os.path.isdir(fdir):
            known = {f.fname for f in fm.fontManager.ttflist}
            for f in glob.glob(os.path.join(fdir, "*.ttf")):
                if f not in known:
                    fm.fontManager.addfont(f)
    except Exception:
        pass
    base, secondary, tick = sizes
    boxed = (frame == "boxed")
    rc = {
        "font.family": "sans-serif",
        "font.size": base,
        "axes.labelsize": base,
        "axes.titlesize": base,
        "legend.fontsize": secondary,
        "xtick.labelsize": tick,
        "ytick.labelsize": tick,
        "axes.linewidth": 0.6,
        "xtick.direction": "out", "ytick.direction": "out",
        "xtick.major.size": 3, "ytick.major.size": 3,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "axes.spines.top": boxed, "axes.spines.right": boxed,
        "axes.spines.left": frame != "none", "axes.spines.bottom": frame != "none",
        "axes.grid": bool(grid),
        "legend.frameon": False,
        "figure.dpi": 200,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "axes.titleweight": "normal",
        "axes.titlelocation": "left",
        "axes.labelweight": "normal",
        "lines.linewidth": 1.2,
        "patch.linewidth": 0.6,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    }
    if font:
        rc["font.sans-serif"] = [font, "DejaVu Sans"]
    mpl.rcParams.update(rc)



apply_figure_style(frame="open")

SWEEP = "../data/optimal_box_surface.csv"
CERTS = {"C1": "../certificates/cert_cmin_C1.json", "C2": "../certificates/cert_cmatch_C2.json",
         "D3a": "../certificates/cert_cmin_D3a.json", "D3c": "../certificates/cert_cmin_D3c.json"}
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

DISP = {"C1": "$L=12$", "C2": "$L=14$", "D3a": "D3a", "D3c": "D3c"}

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
                label="Galerkin-tier enclosure" if k == 0 else None)
        ax.annotate(DISP.get(lab, lab), (c["L"], c["width"]), xytext=(9, -1),
                    textcoords="offset points", color=STAR, fontsize=6, va="center")
    if prior is not None:
        ax.plot([prior[0]], [prior[1]], "s", color=PRIOR, ms=7, mfc="none", mew=1.6,
                zorder=4, label="$L=20$ comparison run")
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
    axa.set_ylabel("Galerkin-tier enclosure width\n(paper units, $=\\frac{1}{2}$ Ha)")
    axa.set_title("H$_2^+$: the balance point is interior in $L$", fontsize=8, loc="left")
    axa.legend(frameon=False, fontsize=6, loc="upper left")
    b = sw[(sw.system == "H3plus") & (sw.d_paper == 4.0)]
    panel(axb, b, ["D3a", "D3c"], cw, ref_lambda=float(b.lambda_ref.dropna().iloc[0]))
    axb.set_title("Linear H$_3^{2+}$, $d=4$: Galerkin-tier widths exceed the float estimates",
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
