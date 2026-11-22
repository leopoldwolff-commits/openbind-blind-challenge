"""Every number quoted in notes/PREREGISTRATION.md, recomputed from the csvs.

    python3 scripts/stats.py

Needs analysis/merged.csv (descriptors + confidence register) and
analysis/spread.csv (pose disagreement). No answers are used anywhere.
"""
import csv, math, statistics as st

def load():
    sp = {r["id"]: r for r in csv.DictReader(open("analysis/spread.csv"))}
    rows = []
    for r in csv.DictReader(open("analysis/merged.csv")):
        if r["id"] in sp:
            r.update({k: v for k, v in sp[r["id"]].items() if k != "id"})
            rows.append(r)
    return rows

def rank(xs):
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    rk = [0.0] * len(xs); i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and xs[order[j+1]] == xs[order[i]]: j += 1
        avg = (i + j) / 2 + 1
        for k in range(i, j + 1): rk[order[k]] = avg
        i = j + 1
    return rk

def pearson(a, b):
    ma, mb = st.mean(a), st.mean(b)
    num = sum((x - ma) * (y - mb) for x, y in zip(a, b))
    den = math.sqrt(sum((x - ma) ** 2 for x in a) * sum((y - mb) ** 2 for y in b))
    return num / den if den else 0.0

def spearman(a, b):
    return pearson(rank(a), rank(b))

def quartiles(xs):
    s = sorted(xs); n = len(s)
    return s[0], s[n // 4], s[n // 2], s[3 * n // 4], s[-1]

def main():
    rows = load()
    col = lambda k: [float(r[k]) for r in rows]
    S, C, SIM = col("spread_mean"), col("chain_pair_iptm_B_L"), col("max_sim_in_set")
    print(f"{len(rows)} compounds\n")

    print("-- distributions (min, q1, median, q3, max)")
    print("   disagreement  " + "  ".join(f"{v:.2f}" for v in quartiles(S)) + "  A")
    print("   confidence    " + "  ".join(f"{v:.3f}" for v in quartiles(C)))

    print("\n-- seed against sample")
    w, b = st.mean(col("within_seed")), st.mean(col("between_seed"))
    print(f"   within-seed {w:.2f} A, between-seed {b:.2f} A, "
          f"seed adds {100*(b-w)/w:.0f} percent")

    print("\n-- are the two measures the same signal")
    print(f"   disagreement vs confidence  Spearman {spearman(S, C):+.3f}")
    print(f"   disagreement vs similarity  Spearman {spearman(S, SIM):+.3f}")
    print(f"   confidence   vs similarity  Spearman {spearman(C, SIM):+.3f}")

    print("\n-- the size artifact flagged in advance")
    for v in ("mw", "heavy_atoms"):
        print(f"   confidence vs {v:12s} {spearman(C, col(v)):+.3f}"
              f"      disagreement vs {v:12s} {spearman(S, col(v)):+.3f}")

    print("\n-- the part of disagreement confidence does not explain")
    rs, rc = rank(S), rank(C)
    mc, ms = st.mean(rc), st.mean(rs)
    beta = sum((c - mc) * (s - ms) for c, s in zip(rc, rs)) / sum((c - mc) ** 2 for c in rc)
    resid = [s - (ms + beta * (c - mc)) for s, c in zip(rs, rc)]
    print(f"   residual vs similarity      Spearman {spearman(resid, SIM):+.3f}")

    print("\n-- compound class")
    frag = [float(r["spread_mean"]) for r in rows if r["fragment_like"] == "1"]
    lead = [float(r["spread_mean"]) for r in rows if r["fragment_like"] != "1"]
    print(f"   fragment-like n={len(frag):3d}  median disagreement {st.median(frag):.2f} A")
    print(f"   lead-like     n={len(lead):3d}  median disagreement {st.median(lead):.2f} A")

    print("\n-- self-inconsistent beyond the 1.5 A tolerance")
    n = sum(1 for x in S if x > 1.5)
    print(f"   {n} compounds ({100*n/len(S):.0f} percent)")

    print("\n-- the two named predictions")
    a = sorted(r["id"] for r in rows
               if float(r["chain_pair_iptm_B_L"]) > 0.95 and float(r["spread_mean"]) > 1.5)
    b_ = [r["id"] for r in rows
          if float(r["chain_pair_iptm_B_L"]) <= 0.95 and float(r["spread_mean"]) <= 1.5]
    print(f"   confident yet self-contradicting  n={len(a)}: "
          + ", ".join(x.split('_')[-1] for x in a))
    print(f"   doubted yet self-consistent       n={len(b_)}")
    print(f"   conflict set                      n={len(a)+len(b_)}")

if __name__ == "__main__":
    main()
