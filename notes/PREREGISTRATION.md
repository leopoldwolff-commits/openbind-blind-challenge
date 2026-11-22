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

### Scope, decided 20 November before the replicates were run

A timing test on three compounds settled the cost. The expensive step is the trunk,
the model building its internal picture of the protein and the compound. Drawing
several poses out of that picture afterwards is nearly free: one seed producing one
pose took 27 seconds, one seed producing five took 24.5. Cost therefore scales with
compounds times seeds, not with poses.

Replicates will be run on **all 400 compounds at five seeds, five poses per seed**,
giving 25 poses each, at a cost of roughly 14 hours. The five seeds are
2746317213, 1181241943, 958682846, 3163119785 and 1812140441, and the same five are
used for every compound, so a compound that scatters more is scattering because of the
compound and not because of which random numbers it drew. No subset is needed.

### The submitted poses change, and why

The original run drew one pose per compound. The replicate run draws 25, and the model
ranks them itself. The submission will therefore be the highest ranked of the 25 rather
than the single pose originally drawn.

This uses only the model's own ranking. No information about the answers enters the
choice, and the ranking was available before any answer was released.

The consequence is that H1 to H3 will be tested against the submitted poses, which are
now the ensemble-selected ones, not the original single draws. This is recorded here
before the replicate run begins. The original single-pose predictions are kept and
remain in the repository, so if the true structures are released rather than only a
score, both sets can be scored and the difference between them reported.

### Amendment, 20 November evening: an interruption, and one thing it revealed

The run stopped after seven batches, at 175 of 400, because a second Colab session had
been left connected and was draining compute. It was not a failure of the run. The
remaining 225 are being completed in order, unchanged: five seeds, 25 poses, same five
seed values. Nothing about the design changes.

One thing worth recording, found while checking whether a partial set would have been
usable. **The 400 are not in random order.** Compounds 1 to 175 are smaller (mean
weight 305 against 326), less similar to the rest of the set (0.599 against 0.658),
more fragment-like (29 percent against 20), and lower confidence on every measure.
Only 12 of 234 scaffolds appear in both halves. The file is ordered by chemical series.

This matters for anyone reading the results, and it would have mattered a great deal
had the run not been finished: the first 175 are the harder half, not a sample. It is
also a reason not to analyse the compounds in file order, and not to split them into
training and held-out groups by position in the file.

### Amendment, 22 November: what the replicates show, before any answer exists

All 400 now have 25 poses each, five seeds by five samples. The disagreement between
those poses was measured after superposing each pose on the protein backbone, so it is
the molecule moving in the pocket and not the complex drifting in space. Everything
below was computed from the model's own output with no access to any answer, and is
recorded here so none of it can be adjusted afterwards.

**The measure discriminates where confidence does not.** Mean pairwise disagreement
runs from 0.15 to 6.62 Angstrom, quartiles 0.34 and 1.61, median 0.75. Confidence over
the same 400 has quartiles 0.873 and 0.972 and sits against its ceiling. 117 compounds,
29 percent, have 25 poses that average further apart than the 1.5 Angstrom tolerance the
challenge allows. Those cannot all be within tolerance of one true structure.

**Seed identity turned out to contribute almost nothing.** Within-seed disagreement is
1.19 Angstrom, between-seed 1.23. The expectation when this was designed was that
rebuilding the model's internal picture from a new seed would move the answer more than
redrawing coordinates from a fixed picture. It does, by about three percent. Nearly all
the variation is in the diffusion step. One seed with 25 samples would have carried
almost the same information at a fifth of the cost. This is recorded as a result about
the model, not as an error to bury.

**H5 is weaker than hoped.** Disagreement and confidence correlate at Spearman -0.790.
They are substantially the same signal, so disagreement cannot beat confidence by much
overall. H5 now stands or falls on the compounds where the two disagree, which is a
sharper test than the one originally written and a harder one to pass.

**The size objection, which this document raised in advance, lands on confidence.
Whether it spares disagreement cannot be settled with 400 compounds.** The concern
recorded before any of this was run was that confidence scores an interface, so a
smaller molecule makes fewer contacts and scores lower for a mechanical reason rather
than real doubt. Confidence correlates with molecular weight at +0.331 and with heavy
atom count at +0.311. Disagreement correlates with the same two at +0.032 and +0.030.

Bootstrapping over scaffold families, 4000 draws, as committed to below:

| | estimate | 95 percent interval |
| --- | --- | --- |
| confidence against molecular weight | +0.331 | +0.113 to +0.505 |
| disagreement against molecular weight | +0.032 | -0.169 to +0.270 |
| the difference | +0.298 | -0.150 to +0.666 |

Confidence's dependence on size clears zero. Disagreement's is consistent with nothing.
**The difference between the two does not clear zero**, and in 10.4 percent of draws
disagreement tracked size at least as strongly as confidence did.

So the claim that disagreement escapes an artifact confidence carries is not supported,
only the weaker pair of statements above. The first reading of these numbers here
inferred a difference from one correlation being distinguishable from zero while the
other was not. That inference is invalid, and testing the difference directly is what
exposed it. Recorded because the mistake is instructive and because an uncorrected
version of it would have been the headline claim of this project.

The intervals are wide because the 400 compounds are 222 scaffold families, 191 of them
holding one member. Resampling compounds rather than families would have produced
narrow intervals and a confident wrong answer.

**The part of disagreement that confidence does not explain still tracks novelty.**
After removing everything confidence accounts for, the residual correlates with
within-set similarity at +0.150. Small, but in the direction H6 predicts.

#### Two named predictions, fixed now

Where confidence and disagreement conflict, one of them is wrong. Both groups are named
in full so neither can be revised later.

*Confident yet self-contradicting*, confidence above 0.95 while the 25 poses average
further apart than the tolerance. These are predicted to fail more often than their
confidence implies: blind_challenge_052, 106, 108, 139, 226.

*Doubted yet self-consistent*, confidence at or below 0.95 while the 25 poses agree
within tolerance. 69 of the 181 low-confidence compounds. These are predicted to succeed
more often than their confidence implies. The list is `analysis/spread.csv` joined to
`analysis/register.csv` under exactly those two cutoffs.

If the first group fails and the second succeeds at rates their confidence does not
predict, disagreement carries information confidence lacks. If both track their
confidence instead, it does not, and that will be reported as plainly as the alternative.

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
