"""Molecular descriptors for the 400 challenge compounds.

Everything here comes from the SMILES alone, so it is known before any
prediction is made. These are the candidate predictors of success.

    python3 scripts/descriptors.py --smiles data/openbind-blind-challenge-1-TEST-BLINDED.csv \
                                   --out analysis/descriptors.csv
"""
import argparse, csv, pathlib, sys

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import rdFingerprintGenerator
from rdkit import DataStructs

RDLogger.DisableLog("rdApp.*")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--smiles", required=True)
    p.add_argument("--out", default="analysis/descriptors.csv")
    a = p.parse_args()

    rows = list(csv.DictReader(open(a.smiles, newline="")))
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

    recs, fps = [], []
    for r in rows:
        cid, smi = r["Molecule_Name"].strip(), r["SMILES"].strip()
        m = Chem.MolFromSmiles(smi)
        if m is None:
            print(f"  could not parse {cid}", file=sys.stderr)
            continue
        scaffold = MurckoScaffold.GetScaffoldForMol(m)
        recs.append(
            {
                "id": cid,
                "smiles": smi,
                "mw": round(Descriptors.MolWt(m), 2),
                "heavy_atoms": m.GetNumHeavyAtoms(),
                "rotatable_bonds": Lipinski.NumRotatableBonds(m),
                "hbd": Lipinski.NumHDonors(m),
                "hba": Lipinski.NumHAcceptors(m),
                "tpsa": round(Descriptors.TPSA(m), 2),
                "logp": round(Descriptors.MolLogP(m), 2),
                "rings": rdMolDescriptors.CalcNumRings(m),
                "aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(m),
                "fraction_csp3": round(rdMolDescriptors.CalcFractionCSP3(m), 3),
                "stereocentres": len(Chem.FindMolChiralCenters(m, includeUnassigned=True)),
                "scaffold": Chem.MolToSmiles(scaffold) if scaffold.GetNumAtoms() else "",
            }
        )
        fps.append(gen.GetFingerprint(m))

    # How distinctive is each compound within the set? A compound with a close
    # neighbour is less of a novelty test than one standing on its own.
    for i, rec in enumerate(recs):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps)
        sims[i] = -1.0
        rec["max_sim_in_set"] = round(max(sims), 3)
        rec["mean_sim_in_set"] = round(sum(s for s in sims if s >= 0) / (len(sims) - 1), 3)

    # Fragment-like by the rule of three, a common screening convention.
    for rec in recs:
        rec["fragment_like"] = int(
            rec["mw"] <= 300 and rec["hbd"] <= 3 and rec["hba"] <= 3
            and rec["logp"] <= 3 and rec["rotatable_bonds"] <= 3
        )

    out = pathlib.Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(recs[0]))
        w.writeheader()
        w.writerows(recs)

    n_frag = sum(r["fragment_like"] for r in recs)
    mws = sorted(r["mw"] for r in recs)
    scaffolds = {r["scaffold"] for r in recs}
    print(f"{len(recs)} compounds written to {out}")
    print(f"  molecular weight: {mws[0]:.0f} to {mws[-1]:.0f}, median {mws[len(mws)//2]:.0f}")
    print(f"  fragment-like (rule of three): {n_frag} ({100*n_frag/len(recs):.0f}%)")
    print(f"  distinct Murcko scaffolds: {len(scaffolds)}")


if __name__ == "__main__":
    main()
