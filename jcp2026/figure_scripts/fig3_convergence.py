"""Regenerate Figure 3 (fig3_convergence) from figs.json.

figs.json holds the plotted series; both trace to certified artifacts:
  fig3.stageA_width -> certified_stage2_final.json rows[*].bracket_w
  fig3.stageB_width -> upper_lower_truncation.json rows[*].bracket_width
Widths are in paper units (energy = 1/2 hartree).
"""
import json, matplotlib as mpl, matplotlib.pyplot as plt

G = json.load(open('figs.json'))
f3 = G["fig3"]; f4 = G["fig4"]

mpl.rcParams.update({
    "font.family":"sans-serif","font.size":10,
    "axes.labelsize":10,"axes.titlesize":10,
    "legend.fontsize":9,"xtick.labelsize":9,"ytick.labelsize":9,
    "axes.linewidth":0.8,"xtick.direction":"out","ytick.direction":"out",
    "xtick.major.size":3,"ytick.major.size":3,"xtick.major.width":0.8,"ytick.major.width":0.8,
    "axes.spines.top":False,"axes.spines.right":False,"legend.frameon":False,
    "figure.dpi":200,"savefig.dpi":300,"pdf.fonttype":42,"svg.fonttype":"none","lines.linewidth":1.4,
})
fig,(axA,axB)=plt.subplots(1,2,figsize=(6.8,2.28),constrained_layout=True)

cA="#8a8a8a"; cB="#2f6fb2"
axA.plot(f3["N"],f3["stageA_width"],"o-",color=cA,ms=6,label="Stage A (projection, sharp $\\eta$)")
axA.plot(f3["N"],f3["stageB_width"],"s-",color=cB,ms=6,label="Stage B (Lehmann–Goerisch)")
axA.set_yscale("log"); axA.set_xticks([32,48]); axA.set_xlim(28,52); axA.set_ylim(1.5e-4,1.6e-1)
axA.set_xlabel("spectral order $N$ (per axis)")
axA.set_ylabel(r"certified bracket width" "\n" r"(paper units, $=\frac{1}{2}$ Ha)")
rr=f3["stageA_width"][1]/f3["stageB_width"][1]
axA.annotate(f"$\\sim${rr:.0f}$\\times$ tighter\nat $N{{=}}48$",
             (48,(f3['stageA_width'][1]*f3['stageB_width'][1])**0.5),
             xytext=(-6,0),textcoords="offset points",ha="right",va="center",fontsize=8.5,color="#333")
axA.legend(loc="lower left",fontsize=8.5)
axA.set_title("(a) second stage tightens the bracket",loc="left",pad=6)

c1="#c0504d"; c2="#4f81bd"
axB.plot(f4["Omega1"]["N"],f4["Omega1"]["gap"],"o-",color=c1,ms=6,label="$\\Omega_1$ (small box)")
axB.plot(f4["Omega2"]["N"],f4["Omega2"]["gap"],"s-",color=c2,ms=6,label="$\\Omega_2$ (doubled box)")
axB.set_yscale("log"); axB.set_xticks([32,48]); axB.set_xlim(28,52); axB.set_ylim(1.5e-5,2.2e-3)
axB.set_xlabel("spectral order $N$ (per axis)")
axB.set_ylabel(r"Dirichlet–Neumann gap" "\n" r"(paper units, $=\frac{1}{2}$ Ha)")
axB.annotate("truncation floor ($N$-indep.)",(40,f4["Omega1"]["gap"][0]),
             xytext=(0,7),textcoords="offset points",ha="center",va="bottom",fontsize=8.5,color=c1)
axB.legend(loc="lower left",fontsize=8.5)
axB.set_title("(b) box size, not $N$, cuts truncation",loc="left",pad=6)

for ax in (axA,axB):
    for sp in ax.spines.values(): sp.set_linewidth(0.8)
    ax.tick_params(width=0.8)
fig.savefig("fig3_convergence.pdf")