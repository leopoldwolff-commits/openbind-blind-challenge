"""Build OpenFold3 query JSONs for the OpenBind blind challenge.

Schema verified against aqlaboratory/openfold-3 examples/example_inference_inputs.

Writes:
  data/queries/all_400.json        every query in one file (best MSA dedup)
  data/queries/batch_000.json ...  batches for checkpointed running

Usage:
  python3 scripts/make_queries.py --fasta data/rcsb_pdb_29HN.fasta \
      --smiles data/openbind-blind-challenge-1-TEST-BLINDED.csv --batch 25
"""
import argparse, csv, json, pathlib, sys

LIGAND_CHAIN = "L"


def read_fasta(path):
    seqs, name, buf = [], None, []
    for line in pathlib.Path(path).read_text().splitlines():
        line = line.strip()
        if line.startswith(">"):
            if name is not None:
                seqs.append((name, "".join(buf)))
            name, buf = line[1:], []
        elif line:
            buf.append(line)
    if name is not None:
        seqs.append((name, "".join(buf)))
    return seqs


def pick(fieldnames, wanted):
    for w in wanted:
        for f in fieldnames:
            if f and f.strip().lower() == w:
                return f
    for w in wanted:
        for f in fieldnames:
            if f and w in f.strip().lower():
                return f
    return None


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--fasta", required=True)
    p.add_argument("--smiles", required=True)
    p.add_argument("--outdir", default="data/queries")
    p.add_argument("--batch", type=int, default=25)
    a = p.parse_args()

    chains = read_fasta(a.fasta)
    if not chains:
        sys.exit("no sequences in the FASTA")
    letters = "ABCDEFGHIJK"
    protein_chains = []
    for i, (header, seq) in enumerate(chains):
        cid = letters[i]
        protein_chains.append(
            {"molecule_type": "protein", "chain_ids": cid, "sequence": seq}
        )
        print(f"  chain {cid}: {len(seq)} residues  ({header[:60]})")

    rows = list(csv.DictReader(open(a.smiles, newline="")))
    id_col = pick(rows[0].keys(), ["molecule_name", "id", "compound_id", "name"])
    smi_col = pick(rows[0].keys(), ["smiles", "canonical_smiles"])
    if not id_col or not smi_col:
        sys.exit(f"could not find id/SMILES columns in {list(rows[0])}")
    print(f"  id column {id_col!r}, SMILES column {smi_col!r}, {len(rows)} rows")

    queries = {}
    for row in rows:
        cid, smi = row[id_col].strip(), row[smi_col].strip()
        if not cid or not smi:
            continue
        queries[cid] = {
            "chains": protein_chains
            + [
                {
                    "molecule_type": "ligand",
                    "chain_ids": LIGAND_CHAIN,
                    "smiles": smi,
                }
            ]
        }

    out = pathlib.Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    (out / "all_400.json").write_text(json.dumps({"queries": queries}, indent=1))

    ids = list(queries)
    n = 0
    for i in range(0, len(ids), a.batch):
        chunk = {k: queries[k] for k in ids[i : i + a.batch]}
        (out / f"batch_{i // a.batch:03d}.json").write_text(
            json.dumps({"queries": chunk}, indent=1)
        )
        n += 1

    print(f"\n  wrote all_400.json ({len(queries)} queries) and {n} batch files to {out}")


if __name__ == "__main__":
    main()
