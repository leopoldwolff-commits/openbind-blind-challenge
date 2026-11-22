"""How much of the size-artifact result is signal and how much is sampling noise.

    python3 scripts/uncertainty.py

The 400 compounds are not independent: 234 scaffolds, the largest family holding
29 members. Resampling compounds one at a time would pretend there are 400
independent observations and give error bars that are too narrow. So the bootstrap
resamples whole scaffold families, which is what notes/PREREGISTRATION.md commits to.
"""
import csv, math, random, statistics as st
from collections import defaultdict

def rank(xs):
    o = sorted(range(len(xs)), key=lambda i: xs[i]); rk = [0.0]*len(xs); i = 0
    while i < len(o):
        j = i
        while j+1 < len(o) and xs[o[j+1]] == xs[o[i]]: j += 1
        a = (i+j)/2 + 1
        for k in range(i, j+1): rk[o[k]] = a
        i = j+1
    return rk

def pearson(a, b):
    ma, mb = st.mean(a), st.mean(b)
    num = sum((x-ma)*(y-mb) for x, y in zip(a, b))
    den = math.sqrt(sum((x-ma)**2 for x in a)*sum((y-mb)**2 for y in b))
    return num/den if den else 0.0

def spearman(a, b): return pearson(rank(a), rank(b))

sp = {r["id"]: r for r in csv.DictReader(open("analysis/spread.csv"))}
rows = []
for r in csv.DictReader(open("analysis/merged.csv")):
    if r["id"] in sp:
        r.update({k: v for k, v in sp[r["id"]].items() if k != "id"}); rows.append(r)

fam = defaultdict(list)
for r in rows: fam[r["scaffold"]].append(r)
families = list(fam.values())
sizes = sorted(len(f) for f in families)
print(f"{len(rows)} compounds in {len(families)} scaffold families")
print(f"  {sum(1 for s in sizes if s == 1)} families of one, largest holds {sizes[-1]}\n")

def stats_of(sample):
    mw   = [float(r["mw"]) for r in sample]
    conf = [float(r["chain_pair_iptm_B_L"]) for r in sample]
    spr  = [float(r["spread_mean"]) for r in sample]
    a, b = spearman(conf, mw), spearman(spr, mw)
    return a, b, a - b

obs = stats_of(rows)
rng = random.Random(42)
B = 4000
draws = []
for _ in range(B):
    s = []
    for _ in range(len(families)): s += rng.choice(families)
    draws.append(stats_of(s))

def ci(vals):
    v = sorted(vals); return v[int(0.025*len(v))], v[int(0.975*len(v))]

names = ["confidence vs molecular weight", "disagreement vs molecular weight",
         "the difference between them"]
print(f"bootstrap over scaffold families, {B} draws, 95 percent intervals\n")
for k, nm in enumerate(names):
    lo, hi = ci([d[k] for d in draws])
    crosses = lo <= 0 <= hi
    print(f"  {nm:34s} {obs[k]:+.3f}   [{lo:+.3f}, {hi:+.3f}]"
          f"{'   includes zero' if crosses else ''}")

frac = sum(1 for d in draws if d[2] <= 0) / B
print(f"\n  draws where disagreement tracks size at least as strongly "
      f"as confidence does: {frac*100:.1f} percent")
print(f"\n  for scale, confidence vs size explains about "
      f"{100*obs[0]**2:.0f} percent of the rank variation in confidence")
