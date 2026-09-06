import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import json

# skill:figure-style kernel.py (auto-injected on skill load)
META_GREY = "#888888"


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


def set_frame(ax, style="open"):
    show = {"open": (False, False, True, True),
            "boxed": (True, True, True, True),
            "none": (False, False, False, False)}[style]
    for side, vis in zip(("top", "right", "bottom", "left"), show):
        ax.spines[side].set_visible(vis)
        if vis:
            ax.spines[side].set_linewidth(0.6)
    ax.tick_params(direction="out", length=0 if style == "none" else 3, width=0.6)


def panel_letter(ax, letter, dx=-0.18, dy=1.02, case="lower", fontsize=None):
    import matplotlib.pyplot as plt
    if fontsize is None:
        fontsize = plt.rcParams.get("font.size", 8) + 1
    s = letter.lower() if case == "lower" else letter.upper()
    ax.text(dx, dy, s, transform=ax.transAxes,
            fontweight="bold", fontsize=fontsize, va="bottom", ha="left")


apply_figure_style(sizes=(9, 8, 7))

PUB_W, PUB_L, PUB_N = 4.763391629341e-04, 20.0, 64
PUB_T = 8.5 * 3600

# Widths and per-stage timings are read from the configuration certificates, so that this
# figure cannot drift from the tables. Do not hardcode widths here: the values that used to
# be inlined were the pre-certification ones, obtained before the coercivity shift was
# certified, and they are 1.02x-2.97x smaller than the certified widths.
CERT_FILES = {"C1": "../certificates/cert_cmin_C1.json", "C2": "../certificates/cert_cmatch_C2.json",
              "C3a": "../certificates/cert_cmin_C3a.json", "C3c": "../certificates/cert_cmin_C3c.json"}
GEOM = {"C1": dict(system="H2plus", d=2.0, L=12, N=48, Np=64),
        "C2": dict(system="H2plus", d=2.0, L=14, N=64, Np=80),
        "C3a": dict(system="H3plus", d=4.0, L=12, N=32, Np=48),
        "C3c": dict(system="H3plus", d=4.0, L=12, N=24, Np=40)}
CERT = {}
for _lab, _fn in CERT_FILES.items():
    _d = json.load(open(_fn))
    CERT[_lab] = dict(GEOM[_lab], w=_d["CERTIFIED_ENCLOSURE"]["width"], raw=_d)

C1, C2 = CERT["C1"]["raw"], CERT["C2"]["raw"]

import pandas as pd

stages = ["A", "mu2", "lg", "dirichlet"]
nice = {"A": "Stage A\n(projection)", "mu2": "separator\n(oee LG)", "lg": "Lehmann–Goerisch\n(eee)", "dirichlet": "Dirichlet\nupper"}
rows = []
for nm, C, D, Da in (("C1  N=48/N'=64", C1, 15625, 35937), ("C2  N=64/N'=80", C2, 35937, 68921)):
    for s in stages:
        rows.append(dict(config=nm, stage=nice[s], D=D, D_aux=Da,
                         wall_s=C["stages"][s]["wall_seconds"],
                         peak_gb=C["stages"][s]["peak_rss_gb"]))
cost = pd.DataFrame(rows)
cost["frac_of_total"] = cost.groupby("config").wall_s.transform(lambda x: x / x.sum())

tot1, tot2 = C1["wall_seconds_total"], C2["wall_seconds_total"]

order = [nice[s] for s in stages]
cfgs = ["C1  N=48/N'=64", "C2  N=64/N'=80"]
cols = ["#8fb3d1", "#1b3b5f"]
x = np.arange(len(order))
w = 0.36

nice2 = {"A": "Stage A\n(projection)", "mu2": "separator\n(oee LG)",
         "lg": "LG bound\n(eee)", "dirichlet": "Dirichlet\nupper"}
order2 = [nice2[s] for s in stages]

pts = [("published\n$L$=20, $N$=64", PUB_T, PUB_W, "#c1440e", "s", 8),
       ("C1  $L$=12, $N$=48", tot1, CERT["C1"]["w"], "#7b2d8e", "*", 15),
       ("C2  $L$=14, $N$=64", tot2, CERT["C2"]["w"], "#7b2d8e", "*", 15)]

fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.2))
ax = axes[0]
for i, (cfg_, col) in enumerate(zip(cfgs, cols)):
    sub = cost[cost.config == cfg_].set_index("stage").loc[order]
    ax.bar(x + (i - 0.5) * w, sub.wall_s, w, color=col, zorder=3, label=cfg_.replace("  ", "   "))
    for xi, v in zip(x + (i - 0.5) * w, sub.wall_s):
        ax.annotate(f"{v:.0f}", (xi, v), xytext=(0, 2), textcoords="offset points",
                    ha="center", va="bottom", fontsize=6.5, color=col)
ax.set_yscale("log"); ax.set_ylim(60, 2.0e4)
ax.set_xticks(x); ax.set_xticklabels(order2, fontsize=7)
ax.set_ylabel("wall time  (s)")
ax.set_title("Lehmann\u2013Goerisch dominates; Stage A and the\nseparator run at decoupled resolutions", loc="left")
ax.legend(frameon=False, fontsize=7, loc="upper left", handlelength=1.2, borderpad=0.2)
ax.annotate(f"{cost[cost.config == cfgs[1]].set_index('stage').loc[nice['lg']].frac_of_total:.0%} of C2",
            (2 + 0.5 * w, C2["stages"]["lg"]["wall_seconds"]), xytext=(7, 8), textcoords="offset points",
            fontsize=7, color="#1b3b5f", ha="left", va="bottom")
set_frame(ax); panel_letter(ax, "a")

ax = axes[1]
for lab, t, wv, col, mk, ms in pts:
    ax.plot([t / 3600], [wv], mk, color=col, ms=ms, mec="white", mew=0.8, zorder=6, clip_on=False)
ax.annotate("published\n$L$=20, $N$=64", (PUB_T / 3600, PUB_W), xytext=(9, 2),
            textcoords="offset points", ha="left", va="center", fontsize=7.5, color="#c1440e")
ax.annotate("C1  $L$=12, $N$=48\nsame width as published", (tot1 / 3600, CERT["C1"]["w"]),
            xytext=(7, -13), textcoords="offset points", ha="left", va="top",
            fontsize=7.5, color="#7b2d8e")
ax.annotate(f"C2  $L$=14, $N$=64\n{PUB_W / CERT['C2']['w']:.2f}$\\times$ tighter (box choice)",
            (tot2 / 3600, CERT["C2"]["w"]), xytext=(11, 2), textcoords="offset points",
            ha="left", va="center", fontsize=7.5, color="#7b2d8e")
ax.annotate("", xy=(tot2 / 3600, CERT["C2"]["w"]), xytext=(PUB_T / 3600, PUB_W),
            arrowprops=dict(arrowstyle="->", color="#6a6a6a", lw=1.0,
                            shrinkA=8, shrinkB=9, connectionstyle="arc3,rad=0.2"))
ax.annotate("better", (4.4, 3.0e-4), fontsize=7, color="#6a6a6a", ha="center",
            rotation=-34, va="center")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlim(0.22, 55); ax.set_ylim(1.95e-4, 7.2e-4)
ax.set_xlabel("certified end-to-end wall time  (h)")
ax.set_ylabel(r"certified enclosure width" "\n" r"(paper units, $=\frac{1}{2}$ Ha)")
ax.set_title("Choosing the box moves the\ncost–accuracy frontier", loc="left")
set_frame(ax); panel_letter(ax, "b")

fig.tight_layout(w_pad=2.6)
fig.savefig("fig_cost.png", dpi=300, bbox_inches="tight")

ax.annotate("horizontal axis mixes two implementations:\nthe published point predates the threaded assembly", (0.5, 0.02), xycoords="axes fraction", fontsize=6.0, color="#6a6a6a", ha="center", va="bottom")
fig.savefig("fig_cost.pdf", bbox_inches="tight")