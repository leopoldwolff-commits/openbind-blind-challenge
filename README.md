# OpenBind Blind Challenge — Zika NS2B-NS3

My entry to the first [OpenBind blind challenge](https://openbind.uk/): predicting
the bound pose of 400 compounds in the Zika virus NS2B-NS3 protease, with the
experimental structures withheld until after the 16 December 2026 deadline.

## The question I am asking

This is a prospective test of a hypothesis I examined retrospectively in
[`openbind0-generalisation`](../openbind0-generalisation): that prediction success
falls as a compound becomes less similar to the training data, and that a model's
own confidence score does not fully warn us when it is out of its depth.

Here the answers are hidden. So before the results are released I record, for every
compound, the model's confidence, how unusual the compound is, and whether I expect
the prediction to pass. The leaderboard reports one number. This reports whether the
failure mode can be anticipated.

## Method

- Structures predicted with OpenBind-0 (OpenFold3 + the `of3-ob-2025-06-30-174k`
  parameters), run on a single GPU.
- One protein, 400 ligands, so the multiple sequence alignment is computed once and
  reused across the set.
- Scoring follows the challenge: ligand BiSyRMSD < 1.5 A, lDDT-PLI > 0.8, and all
  PoseBusters checks passed. Missing predictions count as failures.

## Repository

```
scripts/     input generation, the Colab run notebook, submission packaging
notes/       plan and running decisions
data/        challenge inputs (not redistributed here)
```

Challenge data is from the OpenBind Consortium and OpenADMET and is not included in
this repository; see the
[challenge dataset](https://huggingface.co/datasets/openadmet/openbind-blind-challenge-1).

## Running the analysis

`scripts/stats.py` needs nothing but Python. Everything else wants numpy:

```
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Then, from the repo root:

```
python3 scripts/test_spread.py   # known-answer tests of the disagreement code
python3 scripts/stats.py         # every number quoted in the pre-registration
python3 scripts/figures.py       # rebuilds figures/
```

Run them from the root, not from inside `scripts/`, or the imports will not resolve.
