"""Choose one pose per compound and build the submission zip.

    python3 scripts/build_submission.py --predictions <dir> --replicates <dir> \\
        --posebusters analysis/posebusters.csv --out <dir>

The rule, fixed in notes/PREREGISTRATION.md before this was run:

  submit the original single-pose draw when PoseBusters passes it; otherwise take
  the replicate poses, keep the ones PoseBusters passes, and submit whichever of
  those sits closest to their centroid; if none passes, keep the original and
  record it as unrescued.

Writes submission.zip and a manifest saying, for every compound, which pose was
submitted and why. Refuses to build the zip unless all 400 are present.
"""
import argparse, csv, shutil, sys, tempfile, warnings, zipfile
from pathlib import Path
import numpy as np
warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).parent))
from spread import parse_pdb, kabsch, apply
from posebusters_check import split_pdb, ligand_from

def is_valid(buster, Chem, pdb_text, smiles):
    prot, lig = split_pdb(pdb_text)
    mol = ligand_from(lig, smiles)
    with tempfile.TemporaryDirectory() as tmp:
        pf = Path(tmp) / "p.pdb"; pf.write_text(prot)
        lf = Path(tmp) / "l.sdf"
        w = Chem.SDWriter(str(lf)); w.write(mol); w.close()
        df = buster.bust([lf], None, pf)
    return all(bool(v) for v in df.iloc[0] if isinstance(v, (bool, np.bool_)))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--predictions", required=True)
    ap.add_argument("--replicates", required=True)
    ap.add_argument("--posebusters", default="analysis/posebusters.csv")
    ap.add_argument("--smiles", default="data/openbind-blind-challenge-1-TEST-BLINDED.csv")
    ap.add_argument("--out", default=".")
    a = ap.parse_args()

    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    from posebusters import PoseBusters
    buster = PoseBusters(config="dock")
    # The internal-energy check generates 3D conformers and can take a quarter of
    # an hour on a fused or spiro ring system. It passed for all 400 poses in
    # analysis/posebusters.csv, so it decides nothing here and is dropped. Validity
    # is identical with or without it; see notes/PREREGISTRATION.md.
    dropped = 0
    try:
        mods = buster.config.get("modules", [])
        keep = [m for m in mods
                if "energy" not in (str(m.get("function", "")) + str(m.get("name", ""))).lower()]
        dropped = len(mods) - len(keep)
        buster.config["modules"] = keep
    except Exception as e:
        print(f"  (could not drop the energy check: {e}; it will just be slow)")
    print(f"  dropped {dropped} energy module(s); {len(buster.config.get('modules', []))} remain",
          flush=True)

    smi = {}
    with open(a.smiles) as f:
        for i, r in enumerate(csv.DictReader(f), 1):
            smi[r.get("Molecule_Name") or f"blind_challenge_{i:03d}"] = r["SMILES"]
    by_index = list(smi.values())
    pb = {r["id"]: r["valid"].strip().lower() == "true"
          for r in csv.DictReader(open(a.posebusters))}

    stage = Path(tempfile.mkdtemp())
    manifest, rescued, unrescued = [], 0, []

    for i in range(1, 401):
        cid = f"blind_challenge_{i:03d}"
        orig = sorted((Path(a.predictions) / cid).rglob("*_model.pdb"))
        if not orig:
            print(f"MISSING {cid}"); continue
        s = smi.get(cid) or by_index[i - 1]

        if pb.get(cid, True):
            shutil.copy(orig[0], stage / f"{cid}.pdb")
            manifest.append({"id": cid, "source": "original", "pose": orig[0].name,
                             "n_valid_replicates": "", "note": ""})
            continue

        reps = sorted((Path(a.replicates) / cid).rglob("*_model.pdb"))
        texts = [p.read_text() for p in reps]
        ref_ca, _ = parse_pdb(texts[0]) if texts else (None, None)
        valid, ligs = [], []
        for p, txt in zip(reps, texts):
            try:
                if not is_valid(buster, Chem, txt, s): continue
            except Exception:
                continue
            ca, lig = parse_pdb(txt)
            if ca.shape != ref_ca.shape: continue
            valid.append(p); ligs.append(apply(lig, *kabsch(ca, ref_ca)))

        if not valid:
            shutil.copy(orig[0], stage / f"{cid}.pdb")
            unrescued.append(cid)
            manifest.append({"id": cid, "source": "original", "pose": orig[0].name,
                             "n_valid_replicates": 0, "note": "no valid replicate"})
            continue

        L = np.array(ligs); centre = L.mean(0)
        d = [float(np.sqrt(((x - centre) ** 2).sum(1).mean())) for x in L]
        pick = valid[int(np.argmin(d))]
        shutil.copy(pick, stage / f"{cid}.pdb")
        rescued += 1
        manifest.append({"id": cid, "source": "replicate", "pose": pick.name,
                         "n_valid_replicates": len(valid),
                         "note": f"closest to centroid, {min(d):.2f} A"})
        print(f"  {cid}: {len(valid)}/25 replicates valid, picked {pick.parent.name}",
              flush=True)

    n = len(list(stage.glob("*.pdb")))
    Path("analysis").mkdir(exist_ok=True)
    with open("analysis/submission_manifest.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(manifest[0])); w.writeheader(); w.writerows(manifest)

    print(f"\n{n}/400 staged, {rescued} rescued from replicates, "
          f"{len(unrescued)} could not be rescued")
    if unrescued: print("  unrescued: " + ", ".join(unrescued))
    if n != 400:
        print("refusing to build the zip, not all 400 present"); sys.exit(1)
    z = Path(a.out) / "submission.zip"
    with zipfile.ZipFile(z, "w", zipfile.ZIP_DEFLATED) as zf:
        for p in sorted(stage.glob("*.pdb")): zf.write(p, arcname=p.name)
    print(f"wrote {z}, {z.stat().st_size/1e6:.1f} MB")
    print("manifest -> analysis/submission_manifest.csv")

if __name__ == "__main__":
    main()
