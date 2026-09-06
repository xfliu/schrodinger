"""Regenerates figures/fig1_framework.pdf for the JCP manuscript.

Two colours only: ink (#1a1a1a) and one accent hue (#0b6b8f) with its light tint.
Drawn at 6.9 x 3.15 in and included at \\textwidth, so no downscaling shrinks the
labels below their nominal size.

The box labels state the conditions the pipeline actually verifies at run time:
a certified C_eps that is an UPPER bound on C_eps^opt (never the sharp value,
which lies on the inadmissible side), and the separator conditions
rho > Lambda_n and rho <= lambda_{m+1}.
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

INK, ACC, TINT = "#1a1a1a", "#0b6b8f", "#e8f1f5"
mpl.rcParams.update({"font.size": 9.0, "pdf.fonttype": 42, "ps.fonttype": 42,
                     "savefig.bbox": "tight"})

fig, ax = plt.subplots(figsize=(6.9, 3.35))
ax.set_xlim(0, 100)
ax.set_ylim(0, 44)
ax.axis("off")

# --- geometry, in axis units -------------------------------------------------
# The y-axis spans 44 units over the figure height, so one unit is a fixed
# number of points; text extents below are expressed in units on that basis.
# Placement is computed from the content rather than hard-coded offsets, so a
# one-line and a two-line body both keep the same clearance from the border.
FS_HEAD, FS_BODY = 8.4, 7.4
LINESP = 1.55
PT_PER_UNIT = 44.0 / (3.35 * 72.0)          # units per point
PAD_TOP, PAD_BOT, GAP = 1.9, 1.9, 1.7        # interior margins, in units


def _extent(fontsize, nlines):
    """Approximate vertical extent of a text block, in axis units."""
    per_line = fontsize * (LINESP if nlines > 1 else 1.0)
    return (per_line * (nlines - 1) + fontsize) * PT_PER_UNIT


def box(x, y, w, h, head, body, accent=False):
    ec = ACC if accent else INK
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                                boxstyle="round,pad=0.0,rounding_size=1.2",
                                linewidth=1.6 if accent else 0.95,
                                edgecolor=ec, facecolor=TINT if accent else "none"))
    nb = body.count("\n") + 1
    h_head = _extent(FS_HEAD, 1)
    h_body = _extent(FS_BODY, nb)
    # head hangs from the top pad; body is centred in what is left above the
    # bottom pad, so a single-line body does not float against the border.
    y_head = y + h - PAD_TOP - h_head / 2.0
    top_of_body_zone = y_head - h_head / 2.0 - GAP
    bot_of_body_zone = y + PAD_BOT
    y_body = 0.5 * (top_of_body_zone + bot_of_body_zone)
    assert y_body - h_body / 2.0 >= y + 0.9, (head, "body too close to bottom border")
    assert y_head + h_head / 2.0 <= y + h - 0.9, (head, "head too close to top border")
    ax.text(x + w / 2, y_head, head, ha="center", va="center", fontsize=FS_HEAD,
            color=ec, fontweight="bold" if accent else "normal")
    ax.text(x + w / 2, y_body, body, ha="center", va="center",
            fontsize=FS_BODY, color=INK, linespacing=LINESP)


def arrow(x0, y0, x1, y1, label=None, ly=2.0):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>",
                                 mutation_scale=8.5, linewidth=0.95, color=INK,
                                 shrinkA=0, shrinkB=0))
    if label:
        ax.text((x0 + x1) / 2, (y0 + y1) / 2 + ly, label, ha="center", va="center",
                fontsize=7.2, color=INK)


ax.text(50, 42.2,
        "Guaranteed two-sided enclosure of the true $\\mathbb{R}^3$ ground-state energy",
        ha="center", va="center", fontsize=9.0, color=INK)

ROW1_Y, ROW2_Y, BOX_H = 25.6, 6.4, 13.2
R1M, R2M = ROW1_Y + BOX_H / 2.0, ROW2_Y + BOX_H / 2.0

# Widths are NOT uniform: each box is sized to its own widest line so that the
# side clearance is comparable across the figure. "Stage B (Lehmann-Goerisch)"
# is the longest string in the diagram and needs the widest box; the
# "separation" box holds the shortest content and gives the space back.
# Inter-box gaps are wide enough to hold the arrow labels clear of both
# borders: the row-2 labels U^Dir and L^Neu_LG are the constraint, so those
# gaps are 6 units against row 1's 5.5.
box(2.5, ROW1_Y, 26.5, BOX_H, "Stage A  (projection)",
    "certified $C_\\epsilon \\geq C_\\epsilon^{\\mathrm{opt}}$\n$\\rightarrow\\ L_1,\\ L_2$")
box(34.5, ROW1_Y, 20, BOX_H, "separation",
    "$\\rho > \\Lambda_n$,\u2002$\\rho \\leq \\lambda_{m+1}$")
box(60.0, ROW1_Y, 37, BOX_H, "Stage B  (Lehmann\u2013Goerisch)",
    "seeded by $\\rho$\n$\\rightarrow$ sharpened lower bound")
arrow(29.0, R1M, 34.3, R1M)
arrow(54.5, R1M, 59.8, R1M, label="$\\rho$")

ax.text(50, 22.0, "lift from truncated box $\\Omega$ to $\\mathbb{R}^3$",
        ha="center", va="center", fontsize=7.6, color=INK, style="italic")
arrow(84.0, ROW1_Y, 84.0, ROW2_Y + BOX_H)

box(0.2, ROW2_Y, 30.3, BOX_H, "Dirichlet bound on $\\Omega$",
    "domain monotonicity $\\Omega\\subset\\mathbb{R}^3$\n"
    "$\\Rightarrow$ upper bound on $\\lambda_1(\\mathbb{R}^3)$")
box(36.5, ROW2_Y, 29.6, BOX_H, "certified $\\mathbb{R}^3$ enclosure",
    "$\\lambda_1(\\mathbb{R}^3)\\in[\\,L^{\\mathrm{Neu}}_{\\mathrm{LG}},\\ "
    "U^{\\mathrm{Dir}}\\,]$", accent=True)
box(72.1, ROW2_Y, 27.6, BOX_H, "Neumann bound on $\\Omega$",
    "confinement $\\sigma(\\Omega) > \\lambda_1$\n"
    "$\\Rightarrow$ lower bound on $\\lambda_1(\\mathbb{R}^3)$")
arrow(30.5, R2M, 36.2, R2M, label="$U^{\\mathrm{Dir}}$", ly=2.4)
arrow(72.1, R2M, 66.4, R2M, label="$L^{\\mathrm{Neu}}_{\\mathrm{LG}}$", ly=2.6)

ax.text(50, 2.6, "all steps in verified interval arithmetic",
        ha="center", va="center", fontsize=7.8, color=INK)

fig.savefig("fig1_framework.pdf")
fig.savefig("fig1_framework.png", dpi=300)
