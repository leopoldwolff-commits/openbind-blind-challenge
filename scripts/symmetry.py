"""How much does molecular symmetry inflate the disagreement measure?

    python3 scripts/symmetry.py

scripts/spread.py matches atoms by their order in the file. If a molecule has a
symmetry, the same pose written two different ways round looks like a large
movement even though nothing has changed chemically. The challenge's own metric,
BiSyRMSD, is symmetry-corrected for exactly this reason. Ours is not, so some of
the measured disagreement may be an artifact.

This counts each molecule's automorphisms (the number of ways it maps onto itself)
and asks whether high-disagreement compounds are disproportionately the symmetric
ones. If they are, the measure is inflated and that has to be said out loud.
"""
import csv, math, statistics as st
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")

def automorphisms(smiles):
    m = Chem.MolFromSmiles(smiles)
    if m is None: return None
    return len(m.GetSubstructMatches(m, uniquify=False, useChirality=False,
                                     maxMatches=10000))

def rank(xs):
    o = sorted(range(len(xs)), key=lambda i: xs[i]); rk=[0.0]*len(xs); i=0
    while i < len(o):
        j = i
        while j+1 < len(o) and xs[o[j+1]] == xs[o[i]]: j += 1
        a = (i+j)/2+1
        for k in range(i, j+1): rk[o[k]] = a
        i = j+1
    return rk

def pearson(a, b):
    ma, mb = st.mean(a), st.mean(b)
    num = sum((x-ma)*(y-mb) for x, y in zip(a, b))
    den = math.sqrt(sum((x-ma)**2 for x in a)*sum((y-mb)**2 for y in b))
    return num/den if den else 0.0

spearman = lambda a, b: pearson(rank(a), rank(b))

sp = {r["id"]: r for r in csv.DictReader(open("analysis/spread.csv"))}
rows = []
for r in csv.DictReader(open("analysis/merged.csv")):
    if r["id"] in sp:
        r.update({k: v for k, v in sp[r["id"]].items() if k != "id"})
        r["aut"] = automorphisms(r["smiles"])
        rows.append(r)
rows = [r for r in rows if r["aut"]]
print(f"{len(rows)} molecules\n")

sym  = [r for r in rows if r["aut"] > 1]
asym = [r for r in rows if r["aut"] == 1]
print(f"symmetric   {len(sym):3d} ({100*len(sym)/len(rows):.0f} percent)  "
      f"median disagreement {st.median([float(r['spread_mean']) for r in sym]):.2f} A")
print(f"asymmetric  {len(asym):3d} ({100*len(asym)/len(rows):.0f} percent)  "
      f"median disagreement {st.median([float(r['spread_mean']) for r in asym]):.2f} A")

S   = [float(r["spread_mean"]) for r in rows]
AUT = [float(r["aut"]) for r in rows]
print(f"\ndisagreement against number of symmetries  Spearman {spearman(S, AUT):+.3f}")

hi = [r for r in rows if float(r["spread_mean"]) > 1.5]
hs = sum(1 for r in hi if r["aut"] > 1)
base = len(sym) / len(rows)
print(f"\nof the {len(hi)} compounds above the 1.5 A tolerance, {hs} are symmetric "
      f"({100*hs/len(hi):.0f} percent)")
print(f"the set as a whole is {100*base:.0f} percent symmetric")
print("-> symmetric molecules are " +
      ("over" if hs/len(hi) > base else "under") +
      "-represented among the high-disagreement compounds")

worst = sorted(rows, key=lambda r: -r["aut"])[:5]
print("\nmost symmetric molecules in the set:")
for r in worst:
    print(f"  {r['id']}  {r['aut']:4d} symmetries  disagreement {float(r['spread_mean']):.2f} A")
