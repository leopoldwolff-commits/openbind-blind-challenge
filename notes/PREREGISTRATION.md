# Pre-registration

Written 19 November 2026, with the experimental answers still withheld by the
OpenBind consortium until after the 16 December submission deadline. The git
history of this file is the record that these predictions predate the results.

## What has been done

All 400 compounds predicted with OpenBind-0 (OpenFold3 architecture, parameters
`of3-ob-2025-06-30-174k`), one diffusion sample per compound, seed 42, MSAs from
the ColabFold server, structural templates enabled, low-memory preset, on a single
L4 GPU. 400 of 400 completed, none failed.

Two datasets exist, both built before any answer is known:

- `analysis/descriptors.csv` — molecular properties computed from the SMILES alone
- `analysis/register.csv` — the model's own confidence metrics for each prediction

## What is already known from them

Interface confidence (`chain_pair_iptm_B_L`) is heavily right-skewed: median 0.957,
upper quartile 0.972, maximum 0.984, so the top half of the set spans a range of
0.027. All of the spread sits in the bottom quartile, from 0.546 to 0.873.

That confidence is strongly structured by molecular properties. Spearman
correlations: mean within-set similarity +0.674, nearest-neighbour similarity
+0.476, aromatic rings +0.417, molecular weight +0.331, heavy atoms +0.329.
Fragment-like compounds have median confidence 0.850 against 0.964 for lead-like
ones, and 58 of the 95 fragments fall in the lowest confidence quartile.

No prediction was flagged with a steric clash, so `has_clash` has no variance and
cannot serve as a predictor.

## Hypotheses

Stated before the answers are released. Success is the challenge's own definition:
ligand BiSyRMSD < 1.5 A, lDDT-PLI > 0.8, and all PoseBusters checks passed.

**H1. Fragments fail more often than lead-like compounds.** Small molecules offer
less shape to anchor a pose, and a 1.5 A threshold is proportionally a much larger
error on a 13-atom molecule than a 24-atom one.

**H2. Interface confidence predicts success,** but weakly in the upper range, because
it saturates. Any discriminating power should come from the bottom quartile.

**H3. Confidence under-warns on unfamiliar compounds.** In a logistic regression of
success on both confidence and structural similarity, the similarity term remains
significant. This is the prospective test of the retrospective finding in
`openbind0-generalisation`, and it is the hypothesis this project exists to test.

## Problems identified in advance

**Confidence and novelty are confounded** at rho = 0.674. Separating their
contributions is the central estimation difficulty, and collinearity will widen the
standard errors on both terms. The sample is 400, which may not be enough to
resolve them cleanly. If it is not, that is a reportable result rather than a failure.

**The confidence metric may be mechanically biased by size.** `iptm` measures an
interface, and a smaller ligand makes fewer contacts, so lower confidence on
fragments may be an artifact of the metric rather than genuine uncertainty.
Controlling for contact count would distinguish the two accounts.

**The 400 are not independent.** There are 222 Murcko scaffolds, 191 of them
appearing once, while the largest family has 29 members. Standard errors will be
clustered by scaffold, as they were clustered by test system in the earlier work.

**Only one diffusion sample per compound was generated,** so part of each outcome is
seed luck rather than difficulty. Without repeat sampling the noise floor is unknown,
and effects smaller than it cannot be claimed. A multi-seed run on a subset would
measure it, and remains possible before the deadline.

**Templates were enabled.** The model was shown known NS2B-NS3 structures. This is
another channel by which resemblance to existing knowledge helps, and it is not
separated from the training-data effect in this design.

## What is deliberately not claimed

Within-set similarity is a proxy. The hypothesis concerns similarity to OpenBind-0's
training data, which would require comparison against protein-ligand complexes
deposited before 30 June 2025. That measure has not been computed here, and any
finding about novelty rests on the proxy until it is.
