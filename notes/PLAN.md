# OpenBind Blind Challenge — Zika NS2B-NS3

Final submission: **16 December 2026, 23:59 UTC**
Task: predict a 3D protein-ligand complex for each of 400 compounds (one PDB file each, zipped).

## The angle

This is a prospective test of the hypothesis I already examined retrospectively in
`openbind0-generalisation`: prediction success falls as a compound becomes less
similar to the training data, and confidence scores do not fully warn us about it.

That earlier work used results the consortium had already published. Here the answers
are hidden, the compounds range from small fragments to lead-like molecules, and I
record my predicted success *before* the results are released. The report then says
which cases I expect to fail and why, not only where I placed.

## Fixed plan

| Step | What | When |
| --- | --- | --- |
| 1 | Join the HuggingFace space and Discord, download FASTA + 400 SMILES | this week |
| 2 | Decide compute route, run 2-3 test compounds end to end | by 29 Nov |
| 3 | ETH application submitted | 30 Nov |
| 4 | Compute the MSA once, run all 400 in batches with checkpointing | 1-10 Dec |
| 5 | EPFL application submitted | 12 Dec (3 days early) |
| 6 | Package, validate, submit, write the methods report | 13-15 Dec |
| 7 | Prediction register + analysis for the letter | 16-20 Dec |

Buffer: aim to have a complete submittable ZIP by **12 December**, even if some
predictions are poor. A complete mediocre submission beats an incomplete good one.

## The one optimisation that makes this feasible

The protein is the **same for all 400 compounds**. Only the ligand changes.
So the multiple sequence alignment is computed **once** and reused for every run.
Without this, 400 MSA searches dominate the runtime.

## Compute routes, in order of preference

1. **Google Colab Pro (~10 EUR/month), A100 or L4.** Run openfold-3 from the
   aqlaboratory repo with the OpenBind-0 weights. Checkpoint after every compound,
   so a dropped session loses one prediction and not the batch.
2. **NVIDIA build.nvidia.com** hosts OpenFold3 as a hosted endpoint with free
   starter credits. Check whether 400 calls fit the free tier.
3. **Tamarind Bio** runs OpenBind-0 with no code at all. Simplest, but commercial;
   check the price for 400 jobs before committing.

Fallback if co-folding proves too slow: dock all 400 into the provided apo structure
(PDB 29HN) with an open docking tool. Lower expected accuracy, but it guarantees a
complete submission, and the comparison between the two is itself a finding.

## Scoring (what counts as success)

A prediction counts only if **all three** hold:
- ligand BiSyRMSD < 1.5 A
- lDDT-PLI > 0.8
- passes every PoseBusters check

Missing or unevaluable predictions count as failures. So submit all 400, always.
PoseBusters is open source and can be run locally on my own predictions before
submitting, which catches chemically implausible poses early.

## Before I submit, write down (for the register)

For each compound, before seeing any result:
- its similarity to publicly available NS2B-NS3 structures
- the model's own confidence
- whether I predict it will pass

This is the part no leaderboard entry will have, and it is the part that belongs
in a statistics application.
