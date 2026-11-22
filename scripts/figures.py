"""Figures that need no answers. Built 22 Nov 2026, before the OpenBind release."""
import csv, math, statistics as st
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

TEAL, CORAL, INK, MUTE = "#1D9E75", "#D85A30", "#2C2C2A", "#888780"
plt.rcParams.update({
    "figure.dpi": 160, "savefig.dpi": 160, "font.size": 9,
    "axes.edgecolor": "#D3D1C7", "axes.linewidth": 0.6,
    "axes.titlesize": 10, "axes.titleweight": "medium", "axes.titlelocation": "left",
    "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": MUTE, "ytick.color": MUTE,
    "xtick.labelcolor": INK, "ytick.labelcolor": INK,
    "xtick.major.size": 3, "ytick.major.size": 3,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.facecolor": "white", "axes.facecolor": "white",
})

sp = {r["id"]: r for r in csv.DictReader(open("analysis/spread.csv"))}
rows = []
for r in csv.DictReader(open("analysis/merged.csv")):
    if r["id"] in sp:
        r.update({k: v for k, v in sp[r["id"]].items() if k != "id"}); rows.append(r)
col = lambda k: np.array([float(r[k]) for r in rows])

def rank(x):
    o = np.argsort(x); rk = np.empty(len(x)); rk[o] = np.arange(1, len(x)+1); return rk
def spear(a, b):
    ra, rb = rank(a), rank(b)
    return float(np.corrcoef(ra, rb)[0, 1])

S, C, SIM, MW = col("spread_mean"), col("chain_pair_iptm_B_L"), col("max_sim_in_set"), col("mw")

# ---- Figure 1: one measure ranks the set, the other does not -----------------
fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.7))
ax[0].hist(C, bins=np.linspace(0.5, 1.0, 36), color=CORAL, edgecolor="white", linewidth=0.4)
ax[0].set_title("The model's own confidence")
ax[0].set_xlabel("chain_pair_iptm_B_L"); ax[0].set_ylabel("compounds")
ax[0].annotate(f"half the set lies\nbetween {np.percentile(C,25):.2f} and {np.percentile(C,75):.2f}",
               xy=(0.03, 0.93), xycoords="axes fraction", va="top", fontsize=8, color=MUTE)

ax[1].hist(S, bins=np.linspace(0, 7, 36), color=TEAL, edgecolor="white", linewidth=0.4)
ax[1].axvline(1.5, color=INK, lw=0.9, ls=(0, (4, 3)))
ax[1].set_title("Disagreement between its own 25 poses")
ax[1].set_xlabel("mean pairwise RMSD (Å)")
ax[1].annotate("1.5 Å, the challenge tolerance", xy=(1.6, ax[1].get_ylim()[1]*0.88),
               fontsize=8, color=INK)
ax[1].annotate(f"{(S>1.5).sum()} compounds ({100*(S>1.5).mean():.0f}%) sit\nabove the tolerance",
               xy=(0.42, 0.62), xycoords="axes fraction", va="top", fontsize=8, color=MUTE)
fig.tight_layout(); fig.savefig("figures/01_discrimination.png", bbox_inches="tight"); plt.close(fig)

# ---- Figure 2: the size artifact, flagged in advance -------------------------
fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.9), sharex=True)
for a, y, c, lab, ttl in ((ax[0], C, CORAL, "confidence", "Confidence rises with molecular size"),
                          (ax[1], S, TEAL, "mean pairwise RMSD (Å)", "Disagreement does not")):
    a.scatter(MW, y, s=9, color=c, alpha=0.45, linewidths=0)
    slope, intercept = np.polyfit(MW, y, 1)
    xs = np.linspace(MW.min(), MW.max(), 2)
    a.plot(xs, intercept + slope*xs, color=INK, lw=1.1)
    a.set_xlabel("molecular weight"); a.set_ylabel(lab); a.set_title(ttl)
    a.annotate(f"Spearman {spear(MW, y):+.2f}", xy=(0.03, 0.93), xycoords="axes fraction",
               va="top", fontsize=8.5, color=INK)
fig.tight_layout(); fig.savefig("figures/02_size_artifact.png", bbox_inches="tight"); plt.close(fig)

# ---- Figure 3: where the two measures conflict -------------------------------
conf_hi, agree = C > 0.95, S <= 1.5
grpA = conf_hi & ~agree          # confident yet self-contradicting
grpB = ~conf_hi & agree          # doubted yet self-consistent
rest = ~(grpA | grpB)
fig, ax = plt.subplots(figsize=(5.4, 4.0))
ax.scatter(C[rest], S[rest], s=11, color=MUTE, alpha=0.30, linewidths=0, label="the two agree")
ax.scatter(C[grpB], S[grpB], s=26, color=TEAL, alpha=0.85, linewidths=0,
           label=f"doubted yet self-consistent ({grpB.sum()})")
ax.scatter(C[grpA], S[grpA], s=46, color=CORAL, linewidths=0,
           label=f"confident yet self-contradicting ({grpA.sum()})")
ax.axhline(1.5, color=INK, lw=0.8, ls=(0, (4, 3)))
ax.axvline(0.95, color=INK, lw=0.8, ls=(0, (4, 3)))
named = ", ".join(sorted(r["id"].split("_")[-1] for r, g in zip(rows, grpA) if g))
ax.annotate("named in advance, predicted to fail:\n" + named,
            xy=(0.02, 0.13), xycoords="axes fraction", va="top",
            fontsize=8, color=INK)
ax.set_xlabel("the model's own confidence"); ax.set_ylabel("mean pairwise RMSD (Å)")
ax.set_title(f"Where the two measures disagree   (Spearman {spear(C, S):+.2f})")
ax.legend(frameon=False, fontsize=8, loc="upper right")
fig.tight_layout(); fig.savefig("figures/03_conflict.png", bbox_inches="tight"); plt.close(fig)

print("wrote 3 figures")
print(f"  conflict set: {grpA.sum()} + {grpB.sum()} = {grpA.sum()+grpB.sum()} compounds")
