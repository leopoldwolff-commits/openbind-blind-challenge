"""Is each submitted pose physically possible? No answers needed.

    python3 scripts/posebusters_check.py --predictions <dir> --out analysis/posebusters.csv

PoseBusters checks a pose against chemistry and against the protein: bond lengths
and angles, internal clashes, ring geometry, stereochemistry, and whether the
molecule overlaps the protein it is supposed to be sitting in. A pose that fails
any of these is wrong whatever its RMSD, so the challenge counts it as a failure.

The model writes protein and ligand into one PDB, and a PDB carries no bond
orders. So the ligand is rebuilt by matching the coordinates against the SMILES
the challenge supplied, which is where the bond orders come from.
"""
import argparse, csv, sys, tempfile, warnings
from pathlib import Path
warnings.filterwarnings("ignore")

def split_pdb(text):
    prot = [l for l in text.splitlines() if l.startswith(("ATOM", "TER"))]
    lig  = [l for l in text.splitlines() if l.startswith("HETATM")]
    return "\n".join(prot) + "\nEND\n", "\n".join(lig) + "\nEND\n"

def ligand_from(block, smiles):
    """Coordinates out of the PDB, bond orders out of the SMILES."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    raw = Chem.MolFromPDBBlock(block, sanitize=False, removeHs=False)
    if raw is None:
        raise ValueError("could not read the ligand block")
    ref = Chem.MolFromSmiles(smiles)
    if ref is None:
        raise ValueError("could not read the SMILES")
    return AllChem.AssignBondOrdersFromTemplate(ref, raw)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--predictions", required=True,
                    help="folder of blind_challenge_NNN/ directories")
    ap.add_argument("--smiles", default="data/openbind-blind-challenge-1-TEST-BLINDED.csv")
    ap.add_argument("--out", default="analysis/posebusters.csv")
    ap.add_argument("--limit", type=int, default=0, help="stop after N, for a trial run")
    a = ap.parse_args()

    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    from posebusters import PoseBusters

    smi = {}
    with open(a.smiles) as f:
        for i, r in enumerate(csv.DictReader(f), 1):
            key = r.get("Molecule_Name") or f"blind_challenge_{i:03d}"
            smi[key] = r["SMILES"]
    by_index = list(smi.values())

    buster = PoseBusters(config="dock")
    rows, failed = [], []
    dirs = sorted(Path(a.predictions).glob("blind_challenge_*"))
    if a.limit: dirs = dirs[:a.limit]
    print(f"{len(dirs)} compounds to check", flush=True)

    for n, d in enumerate(dirs, 1):
        cid = d.name
        pdbs = sorted(d.rglob("*_model.pdb"))
        if not pdbs:
            failed.append((cid, "no pose file")); continue
        try:
            idx = int(cid.split("_")[-1]) - 1
            s = smi.get(cid) or by_index[idx]
            prot, lig = split_pdb(pdbs[0].read_text())
            mol = ligand_from(lig, s)
            with tempfile.TemporaryDirectory() as tmp:
                pf = Path(tmp) / "protein.pdb"; pf.write_text(prot)
                lf = Path(tmp) / "ligand.sdf"
                w = Chem.SDWriter(str(lf)); w.write(mol); w.close()
                df = buster.bust([lf], None, pf)
            rec = {"id": cid}
            for k, v in df.iloc[0].items():
                rec[str(k).replace(" ", "_")] = bool(v) if isinstance(v, bool) else v
            checks = [v for k, v in rec.items() if k != "id" and isinstance(v, bool)]
            rec["n_checks"] = len(checks)
            rec["n_failed"] = sum(1 for v in checks if not v)
            rec["valid"] = rec["n_failed"] == 0
            rows.append(rec)
        except Exception as e:
            failed.append((cid, f"{type(e).__name__}: {e}"))
        if n % 25 == 0:
            print(f"  {n}/{len(dirs)}  {sum(1 for r in rows if r['valid'])} valid so far",
                  flush=True)

    if not rows:
        print("nothing checked"); sys.exit(1)
    keys = list(rows[0])
    Path(a.out).parent.mkdir(exist_ok=True)
    with open(a.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        w.writeheader(); w.writerows(rows)

    nv = sum(1 for r in rows if r["valid"])
    print(f"\n{len(rows)} checked, {nv} fully valid ({100*nv/len(rows):.0f} percent)")
    if failed:
        print(f"{len(failed)} could not be checked:")
        for cid, why in failed[:10]: print(f"  {cid}  {why}")
    bad = {}
    for r in rows:
        for k, v in r.items():
            if isinstance(v, bool) and k != "valid" and not v:
                bad[k] = bad.get(k, 0) + 1
    if bad:
        print("\nmost common failures:")
        for k, c in sorted(bad.items(), key=lambda x: -x[1])[:8]:
            print(f"  {c:4d}  {k}")
    print(f"\n-> {a.out}")

if __name__ == "__main__":
    main()
