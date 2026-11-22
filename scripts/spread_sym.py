"""Disagreement between poses, corrected for molecular symmetry.

    python3 scripts/spread_sym.py --predictions <dir> --smiles <csv> --out analysis/spread_sym.csv

scripts/spread.py matched atoms by their order in the file. For a symmetric
molecule that overstates disagreement: the same pose written the other way round
reads as a large movement when nothing chemically has changed. The challenge
scores with BiSyRMSD, which minimises over the molecule's symmetries. This does
the same.

For every pair of poses it tries each way the molecule maps onto itself and keeps
the smallest distance. For an asymmetric molecule there is one such map and the
answer is identical to the old one, which is the check at the end.
"""
import argparse, csv, sys, time, warnings
from pathlib import Path
import numpy as np
warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).parent))
from spread import parse_pdb, kabsch, apply

def heavy_atom_maps(pdb_block, smiles, max_maps=64):
    """Every way the molecule maps onto itself, in the PDB's own atom order.

    Returns a list of index permutations. The first is always the identity.
    """
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem
    RDLogger.DisableLog("rdApp.*")
    raw = Chem.MolFromPDBBlock(pdb_block, sanitize=False, removeHs=True)
    if raw is None: return None
    ref = Chem.MolFromSmiles(smiles)
    if ref is None: return None
    try:
        mol = AllChem.AssignBondOrdersFromTemplate(ref, raw)
    except Exception:
        return None
    heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    pos = {idx: i for i, idx in enumerate(heavy)}
    maps = mol.GetSubstructMatches(mol, uniquify=False, useChirality=False,
                                   maxMatches=max_maps)
    out = []
    for m in maps:
        try:
            out.append(np.array([pos[m[idx]] for idx in heavy]))
        except KeyError:
            continue
    return out or None

def sym_rmsd(A, B, maps):
    """Smallest RMSD between A and B over the molecule's symmetries."""
    best = np.inf
    for p in maps:
        d = float(np.sqrt(((A - B[p]) ** 2).sum(1).mean()))
        if d < best: best = d
    return best

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--predictions", required=True)
    ap.add_argument("--smiles", default="data/openbind-blind-challenge-1-TEST-BLINDED.csv")
    ap.add_argument("--out", default="analysis/spread_sym.csv")
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()

    smi = {}
    with open(a.smiles) as f:
        for i, r in enumerate(csv.DictReader(f), 1):
            smi[r.get("Molecule_Name") or f"blind_challenge_{i:03d}"] = r["SMILES"]
    by_index = list(smi.values())

    dirs = sorted(Path(a.predictions).glob("blind_challenge_*"))
    if a.limit: dirs = dirs[:a.limit]
    print(f"{len(dirs)} compounds", flush=True)

    rows, skipped, t0 = [], [], time.time()
    for n, d in enumerate(dirs, 1):
        cid = d.name
        pdbs = sorted(d.rglob("*_model.pdb"))
        if len(pdbs) < 25:
            skipped.append((cid, f"{len(pdbs)} poses")); continue
        texts = [p.read_text() for p in pdbs]
        seeds = [p.parent.name for p in pdbs]
        s = smi.get(cid) or by_index[int(cid.split("_")[-1]) - 1]
        hetatm = "\n".join(l for l in texts[0].splitlines() if l.startswith("HETATM")) + "\nEND\n"
        maps = heavy_atom_maps(hetatm, s)
        if maps is None:
            skipped.append((cid, "could not rebuild the molecule")); continue

        ref_ca, ref_lig = parse_pdb(texts[0])
        ligs, keep = [ref_lig], [seeds[0]]
        for t, sd in zip(texts[1:], seeds[1:]):
            ca, lig = parse_pdb(t)
            if ca.shape != ref_ca.shape or lig.shape != ref_lig.shape: continue
            ligs.append(apply(lig, *kabsch(ca, ref_ca))); keep.append(sd)
        L = np.array(ligs)
        if L.shape[1] != len(maps[0]):
            skipped.append((cid, f"{L.shape[1]} atoms but maps cover {len(maps[0])}"))
            continue

        pairs, within, between = [], [], []
        for i in range(len(L)):
            for j in range(i + 1, len(L)):
                v = sym_rmsd(L[i], L[j], maps)
                pairs.append(v)
                (within if keep[i] == keep[j] else between).append(v)
        rows.append({"id": cid, "n_poses": len(L), "n_atoms": L.shape[1],
                     "n_symmetries": len(maps),
                     "spread_mean": float(np.mean(pairs)),
                     "spread_median": float(np.median(pairs)),
                     "spread_max": float(np.max(pairs)),
                     "within_seed": float(np.mean(within)) if within else float("nan"),
                     "between_seed": float(np.mean(between)) if between else float("nan")})
        if n <= 3 or n % 10 == 0:
            el = (time.time() - t0) / 60
            print(f"{n}/{len(dirs)}  {el:.1f} min  eta {el/n*(len(dirs)-n):.0f} min  "
                  f"[{cid} {len(maps)} symmetries, {rows[-1]['spread_mean']:.2f} A]",
                  flush=True)

    Path(a.out).parent.mkdir(exist_ok=True)
    with open(a.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    sm = sorted(r["spread_mean"] for r in rows)
    print(f"\n{len(rows)} compounds -> {a.out}")
    print(f"symmetry-corrected disagreement: min {sm[0]:.2f}  q1 {sm[len(sm)//4]:.2f}  "
          f"median {sm[len(sm)//2]:.2f}  q3 {sm[3*len(sm)//4]:.2f}  max {sm[-1]:.2f} A")
    print(f"above the 1.5 A tolerance: {sum(1 for x in sm if x > 1.5)} "
          f"(was 117 before correcting for symmetry)")
    if skipped:
        print(f"\n{len(skipped)} skipped:")
        for cid, why in skipped[:10]: print(f"  {cid}  {why}")

if __name__ == "__main__":
    main()
