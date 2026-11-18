"""Check and package the 400 predictions into one submission ZIP.

Usage:
    python3 package_submission.py --preds predictions --smiles data/compounds.csv \
                                  --out submission/submission.zip

Refuses to build the ZIP if any compound is missing, so an incomplete submission
never goes out by accident. Every missing prediction scores as a failure.
"""
import argparse, csv, pathlib, sys, zipfile


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--preds", required=True, help="directory of per-compound PDB files")
    p.add_argument("--smiles", required=True)
    p.add_argument("--out", default="submission/submission.zip")
    p.add_argument("--force", action="store_true", help="zip even if some are missing")
    a = p.parse_args()

    rows = list(csv.DictReader(open(a.smiles, newline="")))
    cols = list(rows[0].keys())
    id_col = next((c for c in cols if c and "id" in c.lower()), cols[0])
    ids = [r[id_col].strip() for r in rows if r[id_col].strip()]
    print(f"{len(ids)} compounds expected, using id column {id_col!r}")

    preds = pathlib.Path(a.preds)
    found, missing, empty, multi = {}, [], [], []
    for cid in ids:
        # OpenFold3 names outputs <id>_seed_<n>_sample_<n>_model.pdb, so match
        # on the prefix rather than an exact filename.
        hits = sorted(preds.rglob(f"{cid}.pdb")) or sorted(preds.rglob(f"{cid}_*.pdb"))
        if not hits:
            missing.append(cid)
            continue
        if len(hits) > 1:
            multi.append(cid)
        pick = hits[0]
        if pick.stat().st_size < 500:
            empty.append(cid)
        found[cid] = pick

    print(f"found {len(found)}, missing {len(missing)}, suspiciously small {len(empty)}")
    if multi:
        print(f"{len(multi)} compounds have several samples; taking the first by name."
              " Rank them properly before the real submission.")
    if missing:
        print("missing:", ", ".join(missing[:15]), "..." if len(missing) > 15 else "")
    if empty:
        print("too small:", ", ".join(empty[:15]))
    if (missing or empty) and not a.force:
        sys.exit("\nnot packaging. Fix these, or re-run with --force.")

    out = pathlib.Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as z:
        for cid, path in sorted(found.items()):
            z.write(path, arcname=f"{cid}.pdb")
    mb = out.stat().st_size / 1e6
    print(f"\nwrote {out} with {len(found)} files, {mb:.1f} MB")


if __name__ == "__main__":
    main()
