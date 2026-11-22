"""Known-answer tests. The right answers here are true by construction, so a
failure means the code is wrong, not that the data is odd.

    python3 scripts/test_spread.py
"""
import numpy as np, sys
sys.path.insert(0, "scripts")
from spread import parse_pdb, kabsch, apply, rmsd, superpose_on_protein, disagreement

rng = np.random.default_rng(0)
fails = []
def check(name, cond, detail=""):
    print(f"  {'pass' if cond else 'FAIL'}  {name}{'  ' + detail if detail else ''}")
    if not cond: fails.append(name)

def rand_rotation():
    A = rng.normal(size=(3, 3))
    Q, R = np.linalg.qr(A)
    return Q * np.sign(np.diag(R)) * (1 if np.linalg.det(Q) > 0 else -1)

print("parsing")
PDB = (
"ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00  0.00           N\n"
"ATOM      2  CA  ALA A   1       4.000   5.000   6.000  1.00  0.00           C\n"
"ATOM      3  CA  GLY B   2       7.000   8.000   9.000  1.00  0.00           C\n"
"HETATM    4  C1  LIG L   1      10.000  11.000  12.000  1.00  0.00           C\n"
"HETATM    5  H1  LIG L   1      13.000  14.000  15.000  1.00  0.00           H\n")
ca, lig = parse_pdb(PDB)
check("picks out CA atoms only", ca.shape == (2, 3), f"got {ca.shape}")
check("first CA coordinates right", np.allclose(ca[0], [4, 5, 6]))
check("drops ligand hydrogens", lig.shape == (1, 3), f"got {lig.shape}")

print("\nrmsd")
X = rng.normal(size=(20, 3)) * 5
check("a pose against itself is zero", rmsd(X, X) == 0.0)
shift = np.array([2.0, 0, 0])
check("a pure 2 A shift reads 2.00 A", abs(rmsd(X, X + shift) - 2.0) < 1e-12,
      f"got {rmsd(X, X + shift):.6f}")

print("\nsuperposition")
prot = rng.normal(size=(60, 3)) * 10
lig0 = rng.normal(size=(18, 3)) * 3
R, t = rand_rotation(), rng.normal(size=3) * 20
moved = (prot @ R + t, lig0 @ R + t)          # whole complex rigidly moved
L = superpose_on_protein([(prot, lig0), moved])
check("rigid motion of the whole complex reads as zero disagreement",
      rmsd(L[0], L[1]) < 1e-8, f"got {rmsd(L[0], L[1]):.2e} A")

d = np.array([0.0, 0.0, 1.3])                  # ligand alone slides 1.3 A
L = superpose_on_protein([(prot, lig0), (prot, lig0 + d)])
check("ligand sliding 1.3 A while the protein holds still reads 1.3 A",
      abs(rmsd(L[0], L[1]) - 1.3) < 1e-8, f"got {rmsd(L[0], L[1]):.6f} A")

R2, t2 = rand_rotation(), rng.normal(size=3) * 30
# the ligand slides 1.3 A inside the pocket, then the whole complex is moved:
# protein and ligand must travel together or the input is not a real structure
L = superpose_on_protein([(prot, lig0),
                          (prot @ R2 + t2, (lig0 + d) @ R2 + t2)])
check("same slide, but the complex also moved, still reads 1.3 A",
      abs(rmsd(L[0], L[1]) - 1.3) < 1e-6, f"got {rmsd(L[0], L[1]):.6f} A")

print("\ndisagreement over a set")
ligs = np.array([lig0, lig0 + d, lig0 + 2 * d])
out = disagreement(ligs)
check("three poses give three pairs averaging 1.733 A",
      abs(out["spread_mean"] - (1.3 + 1.3 + 2.6) / 3) < 1e-9,
      f"got {out['spread_mean']:.4f}")
check("the widest pair is 2.6 A", abs(out["spread_max"] - 2.6) < 1e-9)

out = disagreement(np.array([lig0, lig0 + d, lig0 + 10 * d, lig0 + 11 * d]),
                   seeds=["s1", "s1", "s2", "s2"])
check("within-seed picks out the close pairs", abs(out["within_seed"] - 1.3) < 1e-9,
      f"got {out['within_seed']:.4f}")
check("between-seed picks out the far ones", abs(out["between_seed"] - 13.0) < 1e-9,
      f"got {out['between_seed']:.4f}")

print()
if fails:
    print(f"{len(fails)} FAILED: " + ", ".join(fails)); sys.exit(1)
print("all checks pass")

