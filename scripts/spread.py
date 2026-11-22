"""Pose disagreement: how far apart the model's repeated answers land.

Factored out of the Colab notebook so it can be tested. See test_spread.py.
"""
import numpy as np

def parse_pdb(text):
    """Protein CA coordinates and ligand heavy-atom coordinates, in file order."""
    ca, lig = [], []
    for l in text.splitlines():
        if l.startswith("ATOM") and l[12:16].strip() == "CA":
            ca.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
        elif l.startswith("HETATM") and l[76:78].strip() != "H":
            lig.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
    return np.array(ca), np.array(lig)

def kabsch(P, Q):
    """Rotation taking P onto Q, with both centroids. Returns (R, P_centre, Q_centre)."""
    pc, qc = P.mean(0), Q.mean(0)
    V, _, W = np.linalg.svd((P - pc).T @ (Q - qc))
    R = V @ np.diag([1, 1, np.sign(np.linalg.det(V @ W))]) @ W
    return R, pc, qc

def apply(X, R, pc, qc):
    return (X - pc) @ R + qc

def rmsd(A, B):
    """Root mean square deviation between two sets of matched coordinates."""
    return float(np.sqrt(((A - B) ** 2).sum(1).mean()))

def superpose_on_protein(poses):
    """poses: list of (ca, lig). Put every ligand in the first pose's frame.

    The protein is what gets aligned, so what survives is the ligand moving in
    the pocket rather than the whole complex sitting somewhere else in space.
    """
    ref_ca, ref_lig = poses[0]
    out = [ref_lig]
    for ca, lig in poses[1:]:
        if ca.shape != ref_ca.shape or lig.shape != ref_lig.shape:
            continue
        out.append(apply(lig, *kabsch(ca, ref_ca)))
    return np.array(out)

def disagreement(ligands, seeds=None):
    """Mean, median and max pairwise RMSD; split by seed when seeds are given."""
    n = len(ligands)
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    d = {(i, j): rmsd(ligands[i], ligands[j]) for i, j in pairs}
    allv = list(d.values())
    out = {"n_poses": n, "spread_mean": float(np.mean(allv)),
           "spread_median": float(np.median(allv)), "spread_max": float(np.max(allv))}
    if seeds is not None:
        within = [v for (i, j), v in d.items() if seeds[i] == seeds[j]]
        between = [v for (i, j), v in d.items() if seeds[i] != seeds[j]]
        out["within_seed"] = float(np.mean(within)) if within else float("nan")
        out["between_seed"] = float(np.mean(between)) if between else float("nan")
    return out
