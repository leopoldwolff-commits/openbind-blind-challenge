# Pre-registration

Written 19 November 2026 and extended 20 November, while the experimental answers
are still sealed until after the 16 December deadline. The git history of this file
is the record that these predictions came before the results. The wording below was
made plainer on 20 November; the original is in the history and nothing about the
hypotheses changed.

## What has been done

All 400 compounds were predicted once each with OpenBind-0 (OpenFold3 architecture,
parameters `of3-ob-2025-06-30-174k`), one pose per compound, seed 42, alignments from
the ColabFold server, structural templates switched on, low-memory preset, on one L4
GPU. All 400 finished. None failed.

Two datasets exist, both built before any answer is known:

- `analysis/descriptors.csv` — what each molecule is like, worked out from its SMILES
- `analysis/register.csv` — what the model said about its own prediction

## What those already show

**The model is sure about most compounds and unsure about a few.** Its confidence in
the protein-ligand interface has a median of 0.957 and an upper quartile of 0.972, so
the top half of the set sits inside a range of 0.027. Everything interesting is in the
bottom quarter, which runs from 0.546 to 0.873.

**The ones it is unsure about are not random.** They are the small compounds and the
odd ones out. Fragments have a median confidence of 0.850 against 0.964 for larger
lead-like molecules, and 58 of the 95 fragments fall in the lowest confidence quarter.
Confidence tracks how much a compound resembles the rest of the set at a Spearman
correlation of 0.674, and molecular size at 0.33.

**No prediction was flagged for a steric clash.** So that field is the same for all
400 and cannot predict anything.

## Hypotheses

Success means what the challenge says it means: ligand BiSyRMSD under 1.5 A,
lDDT-PLI above 0.8, and every PoseBusters check passed.

**H1. Small fragments fail more often than larger lead-like compounds.** A small
molecule gives the model less shape to anchor against, and 1.5 A is a much bigger
error on a 13-atom molecule than on a 24-atom one.

**H2. The model's confidence predicts success, but only weakly at the top end.**
It is bunched so tightly above 0.95 that it cannot separate those compounds from one
another. Any real predictive power should come from the bottom quarter.

**H3. The confidence under-warns on unfamiliar compounds.** Put both confidence and
structural similarity into a logistic model of success, and similarity still matters.
In plain terms: knowing the model is confident is not enough, because unusual
compounds fail more often than their confidence suggests. This is the prospective
test of what the retrospective work in `openbind0-generalisation` found, and it is
the reason this project exists.

## Amendment, 20 November 2026

Added before any replicate was run, and before any answer is known.

### The idea

The hypotheses above use only what the model *says* about itself. There is a second
signal available: what it *does* when you run it again.

Run the same compound five times from different starting points and you get five
predicted positions. **Ensemble spread is how far apart those five answers are.** Ask
five people to point at where something is hidden: if they all point at the same spot
the group probably knows, and if they point in five directions they are guessing, even
if each one says they are sure.

### How it is measured

Each compound in the replicate set is predicted several times with different model
seeds. For each pair of replicates the structures are lined up on the protein
backbone, to put them in a common frame, and the distance between the ligand
positions is measured as an RMSD. The median of those pairwise distances is that
compound's ensemble spread, in angstroms.

### Additional hypotheses

**H4. Tight agreement means a better prediction.** Compounds whose replicates land in
the same place are more often correct than compounds whose replicates scatter. This
is not guaranteed: a model can be consistently wrong, landing confidently in the same
incorrect spot every time.

**H5. Agreement is a better warning sign than the model's own confidence.** Compared
by area under the ROC curve, which scores how well a number separates right from wrong
on a scale from 0.5 for useless to 1.0 for perfect. The retrospective work scored the
model's self-reported confidence at 0.83. Uncertainty on the comparison comes from a
bootstrap that resamples whole scaffold families, because compounds sharing a scaffold
are not independent.

**H6. Agreement does not have the same blind spot as self-reported confidence.** The
test: once the spread is known, does knowing that a compound is unusual still add bad
news? If it adds nothing, the spread already captured the problem and is the better
measure. If it still adds bad news, both signals are blind in the same way, which
would suggest the model has no internal access to its own unfamiliarity at all.
Either outcome is a result.

### The noise floor

Replicate disagreement also says how much of any single prediction is seed luck rather
than compound difficulty. No effect smaller than this noise floor will be claimed from
the single-pose run, and the floor will be reported whether or not it is convenient.

### Scope

If compute does not stretch to replicates for all 400, they will be run on a subset
spread across the range of molecular size and confidence, chosen and recorded here
before the replicates are run. The single-pose run of all 400 remains the submission
and remains the basis for H1 to H3.

## Problems known in advance

**Confidence and novelty are tangled together** at a correlation of 0.674. Separating
what each contributes is the central difficulty, and it will widen the error bars on
both. With 400 compounds it may not be possible to separate them cleanly. If it is
not, that is a reportable result and not a failure.

**The confidence measure may be biased by size for a boring reason.** It scores an
interface, and a smaller molecule simply makes fewer contacts, so lower confidence on
fragments could be an artifact of the measure rather than real doubt. Controlling for
the number of contacts would tell the two apart.

**The 400 are not independent.** There are 222 scaffolds, 191 of them appearing once,
while the largest family has 29 members. Errors will be clustered by scaffold, as
they were clustered by test system in the earlier work.

**Templates were switched on**, so the model was shown known NS2B-NS3 structures.
That is another route by which resembling existing knowledge helps, and this design
does not separate it from the training-data effect.

## What is deliberately not claimed

Similarity here is measured within the challenge set, which is a proxy. The real
question is similarity to what OpenBind-0 was trained on, which would need comparison
against protein-ligand structures deposited before 30 June 2025. That has not been
computed. Any claim about novelty rests on the proxy until it is.
