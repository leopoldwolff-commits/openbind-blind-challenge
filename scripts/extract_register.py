"""Build the prediction register: one row per compound, every confidence field.

Run where the predictions directory is readable (in Colab, with Drive mounted):

    python3 extract_register.py --preds /content/drive/MyDrive/openbind-blind-challenge/predictions \
                                --out register.csv

This must be done BEFORE the answers are released. The register is the record of
what the model said about itself while the outcome was still unknown.
"""
import argparse, csv, json, pathlib, sys


def flatten(d, prefix=""):
    """confidences_aggregated nests chain_ptm and chain_pair_iptm; spread them out."""
    flat = {}
    for k, v in d.items():
        key = f"{prefix}{k}"
        if isinstance(v, dict):
            for k2, v2 in v.items():
                clean = k2.replace("(", "").replace(")", "").replace(", ", "_").replace(" ", "")
                flat[f"{key}_{clean}"] = v2
        else:
            flat[key] = v
    return flat


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--preds", required=True)
    p.add_argument("--out", default="register.csv")
    a = p.parse_args()

    preds = pathlib.Path(a.preds)
    recs = []
    for d in sorted(preds.glob("blind_challenge_*")):
        if not d.is_dir():
            continue
        conf = sorted(d.rglob("*_confidences_aggregated.json"))
        pdbs = sorted(d.rglob("*_model.pdb"))
        if not conf or not pdbs:
            print(f"  incomplete: {d.name}", file=sys.stderr)
            continue
        rec = {"id": d.name, "n_samples": len(pdbs)}
        rec.update(flatten(json.load(open(conf[0]))))
        txt = pdbs[0].read_text()
        rec["ligand_atoms"] = sum(1 for l in txt.splitlines() if l.startswith("HETATM"))
        rec["pdb_lines"] = len(txt.splitlines())
        timing = d / pdbs[0].parent.name / "timing.json"
        if timing.exists():
            try:
                t = json.load(open(timing))
                rec["seconds"] = list(t.values())[0] if isinstance(t, dict) else t
            except Exception:
                pass
        recs.append(rec)

    if not recs:
        sys.exit("no predictions found")

    cols = list(recs[0])
    for r in recs:
        for k in r:
            if k not in cols:
                cols.append(k)

    out = pathlib.Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(recs)

    print(f"{len(recs)} compounds written to {out}")
    print(f"columns: {', '.join(cols)}")
    key = "chain_pair_iptm_B_L"
    if key in cols:
        vals = sorted(r[key] for r in recs if isinstance(r.get(key), (int, float)))
        n = len(vals)
        print(f"\n{key}: min {vals[0]:.3f}  q1 {vals[n//4]:.3f}  median {vals[n//2]:.3f} "
              f" q3 {vals[3*n//4]:.3f}  max {vals[-1]:.3f}")
    clash = sum(1 for r in recs if r.get("has_clash"))
    print(f"predictions flagged with a clash: {clash}")


if __name__ == "__main__":
    main()
