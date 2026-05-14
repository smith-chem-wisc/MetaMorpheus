# Iterative PEP Training in MetaMorpheus — Technical Report

**Branch:** `pep-iterative-training`  **Date:** 2026-05-14
**Companion docs:** `PEP_Iterative_Training_Methodology.md` (lab notebook),
`PEP_Iterative_Training_Phase1_Diagnostic.md` (go/no-go diagnostic).

---

## Bottom line

Iterative PEP training was implemented, tested, and evaluated end-to-end on real
bottom-up data. It produces a **small but consistent positive gain — +86 target peptides
at PEP q-value < 0.01 (+0.66%)** — and measurably improves the underlying model
(LogLoss 0.0734 → 0.0602, AUC 0.99926 → 0.99950). It costs about **+45% search
wall-clock** on the test dataset. The feature is **opt-in and off by default**; the
single-pass pipeline is byte-for-byte unchanged unless a user sets
`IterativePepTraining = true`.

**Recommendation: merge as an opt-in feature.** The gain on this clean, near-ceiling
dataset is modest, but the implementation is correct, fully tested, and free of regression
risk while off. Two planned follow-up PRs — FastTree hyperparameter tuning and new PEP
features — are expected to compound with it, since iterative relabeling's value scales
with how much headroom the model and feature set have beyond raw score.

---

## 1. Problem

MetaMorpheus trains its PEP model (FastTree, Platt-calibrated, 4-fold cross-fit) on a
positive pool selected **once**, before training, by raw-search-score q-value. Raw score
is one of ~17 features the trained model uses; mid-confidence targets — strong in the full
feature space but unremarkable by raw score — are dropped from training entirely. The
model therefore never learns from the mid-score region where PEP discrimination is most
valuable.

**Iterative training** (the Percolator idea) re-selects the positive pool from the model's
own PEP-derived q-values after the first fit, retrains, and repeats until the pool
stabilizes.

## 2. What was built

A training-iteration loop around the existing 4-fold cross-fit fit in
`PepAnalysisEngine.ComputePEPValuesForAllPSMs()`:

- **Iteration 0** labels positives by raw-score q-value — identical to the historical
  pipeline.
- **Iterations ≥ 1** relabel positives from the previous iteration's **per-fold**,
  **cross-fit-isolated** PEP q-values.
- The loop stops at convergence (positive pool moves < 1%) or a configurable max
  (default 3).

Key design properties:

- **Cross-fit isolation is structural.** Fold *k*'s relabeling uses a PEP q-value map
  built only from fold *k*'s held-out predictions — i.e. from the model trained on folds
  ≠ *k*. The per-fold maps are pairwise disjoint by construction; there is no single
  global PEP that could leak a fold's data into its own labels. A dedicated test asserts
  this partition directly.
- **Ambiguity pruning is deferred to the final pass.** It mutates the PSM, so running it
  every iteration would compound; the PEP value is provably unaffected by deferring it.
- **Decoys never move** between iterations — only target labels are re-selected.
- **Off by default.** `CommonParameters.IterativePepTraining` (bool, default false) gates
  the feature end-to-end.

Per-iteration positive-pool sizes are recorded and printed in the PEP metrics block of
`results.txt`.

## 3. Test coverage

`Test/PepIterativeTrainingTest.cs` — five contracts, all passing:

1. The max iteration count is configurable and clamped to ≥ 1.
2. The loop never runs more than the max.
3. The positive pool size is logged once per iteration (and surfaced in metrics text).
4. Convergence triggers below 1% positive-pool change.
5. With iterative training off (or max = 1), the engine runs exactly one pass and never
   builds a relabel map — the regression safety net.
6. Cross-fit isolation: the per-fold relabel maps partition the PSMs exactly by fold.

**Full NUnit suite after the change: 1,453 passed, 0 failed, 1 skipped — zero
regressions.** No existing test required revision.

## 4. Evaluation

### 4.1 Dataset

Bottom-up human protein lysate, trypsin, two Velos `.raw` files,
`human_rev_GPTMDpruned.xml` database. ~47k PSMs / ~29k peptides. Baseline and iterative
runs differ only in the `IterativePepTraining` flag.

### 4.2 Headline — target peptides by PEP q-value, from `AllPeptides.psmtsv`

| Threshold | Baseline | Iterative | Δ |
|---|---:|---:|---:|
| **PEP q-value < 0.01** | **13,049** | **13,135** | **+86 (+0.66%)** |
| PEP q-value < 0.05 | 13,998 | 14,084 | +86 (+0.61%) |
| PEP q-value < 0.005 | 12,831 | 12,888 | +57 (+0.44%) |
| raw-score q-value < 0.01 (PEP-independent sanity) | 12,651 | 12,644 | −7 (noise) |

The raw-score count barely moves (as expected — raw q-value is not a function of the PEP
model); the PEP-derived counts all move up by a small, consistent margin.

### 4.3 Model-quality metrics

| Metric | Baseline | Iterative |
|---|---:|---:|
| LogLoss | 0.0734 | 0.0602 |
| Area under ROC | 0.99926 | 0.99950 |
| Area under PR curve | 0.99952 | 0.99965 |
| Accuracy | 0.9916 | 0.9946 |
| F1 | 0.9931 | 0.9957 |

The model improves more than the peptide count does — the count is a coarse, near-
saturated readout of an already-excellent model.

### 4.4 Iteration trajectory

Iterative training ran the full 3 iterations (did not converge):

| Iteration | Positive pool | Change |
|---|---:|---:|
| 0 (raw-score labels) | 23,300 | — |
| 1 | 24,959 | +7.12% |
| 2 | 25,222 | +1.05% |

Iteration 0's pool (23,300) is exactly the baseline's positive pool — confirming
iteration-0 labeling is identical to the single-pass pipeline. The iteration 1→2 change
(+1.05%) sits just above the 1% convergence threshold, so the third (final) iteration ran.

### 4.5 Time cost

| | Baseline | Iterative | Δ |
|---|---:|---:|---:|
| FDR analysis engine (contains PEP training) | 15.2 s | 37.9 s | +22.7 s (≈ 2.5×) |
| Whole search task, wall-clock | 49.7 s | 71.9 s | +22.2 s (+45%) |

PEP training does ~3 fitting passes instead of 1. The engine step grows ~2.5× rather than
a full 3× because q-value computation and the constructor's file-specific dictionaries are
fixed overhead. Marginal cost ≈ ~11 s per extra iteration on this dataset. The cost scales
with dataset size; for very large searches, running iterations ≥ 2 with fewer trees is a
reasonable optimization (not done here).

### 4.6 Calibration

Reliability tables — mean predicted PEP vs. observed decoy fraction, 10 bins.

**Baseline**

| PEP bin | n | mean predicted | decoy fraction |
|---|---:|---:|---:|
| [0.0, 0.1) | 12,828 | 0.0005 | 0.0040 |
| [0.1, 0.2) | 50 | 0.1446 | 0.1800 |
| [0.2, 0.3) | 37 | 0.2431 | 0.2162 |
| [0.3, 0.4) | 36 | 0.3581 | 0.2500 |
| [0.4, 0.5) | 26 | 0.4434 | 0.1923 |
| [0.5, 0.6) | 20 | 0.5493 | 0.1500 |
| [0.6, 0.7) | 18 | 0.6484 | 0.2778 |
| [0.7, 0.8) | 26 | 0.7622 | 0.1538 |
| [0.8, 0.9) | 45 | 0.8665 | 0.1333 |
| [0.9, 1.0) | 16,231 | 0.9994 | 0.4777 |

**Iterative**

| PEP bin | n | mean predicted | decoy fraction |
|---|---:|---:|---:|
| [0.0, 0.1) | 13,008 | 0.0003 | 0.0057 |
| [0.1, 0.2) | 26 | 0.1397 | 0.1923 |
| [0.2, 0.3) | 12 | 0.2339 | 0.2500 |
| [0.3, 0.4) | 14 | 0.3426 | 0.1429 |
| [0.4, 0.5) | 14 | 0.4598 | 0.2143 |
| [0.5, 0.6) | 12 | 0.5409 | 0.0833 |
| [0.6, 0.7) | 5 | 0.6524 | 0.0000 |
| [0.7, 0.8) | 18 | 0.7568 | 0.2778 |
| [0.8, 0.9) | 19 | 0.8658 | 0.1579 |
| [0.9, 1.0) | 16,188 | 0.9998 | 0.4792 |

The two are comparable. The `[0.0, 0.1)` bin — the one that drives the headline count —
stays well-calibrated in both; iterative pulls slightly more peptides into it (12,828 →
13,008) at a marginally higher decoy fraction (0.0040 → 0.0057). The mid bins (n ≈ 5–50)
are too sparse for either run to read into. The `[0.9, 1.0)` bin holds the rejected
matches — predominantly decoys and weak targets — and its ~0.48 decoy fraction at
predicted PEP ≈ 1.0 is the expected shape.

## 5. Why the gain is smaller than the diagnostic predicted

The Phase 1 diagnostic measured 7.54% positive-pool growth and, by the Percolator rule of
thumb, suggested a ~2–8% peptide gain. The measured gain is +0.66%. The diagnostic was
optimistic for a knowable reason:

The 907 "would-enter" peptides it counted all had **PEP q-value ≤ 0.005 already** — they
were already passing in the baseline results. Moving them into the *training* pool does
not move them into the *results*; they were never absent. The only available benefit is
indirect: a model trained on a larger, more diverse positive set ranks *other* borderline
peptides slightly better. The +86 peptides are exactly that second-order effect. The
diagnostic correctly predicted the training pool would grow ~7.5% (it grew 8.25% across
three iterations) but over-translated that into a results gain.

The model also began near its ceiling (baseline AUC 0.99926) on clean Orbitrap-class
bottom-up data, leaving little headroom for any change to exploit.

## 6. Recommendation and follow-ups

**Merge as an opt-in feature.** It is correct, tested, regression-free while off, and
positive — if modest — when on. The combination "small gain, real cost, zero risk while
off" is exactly what opt-in is for.

Expected to compound with this work, in separate planned PRs:

- **FastTree hyperparameter tuning** — the current settings (`NumberOfTrees = 400`,
  `LearningRate = 0.2`, `MinimumExampleCountPerLeaf = 10`) are untuned with no early
  stopping. A better-tuned model has more headroom for iterative relabeling to exploit.
- **New PEP features** — iterative training's value scales with how much signal the
  feature set can express beyond raw score; richer features give relabeling more to work
  with.

Other candidate follow-ups:

- Iterations ≥ 2 with fewer trees, to cut the time cost where the model is already mostly
  converged.
- A reliability-diagram in the metrics output, so calibration is monitored routinely.
- Evaluation on larger / lower-quality datasets, where the model has more headroom than it
  did here.

## 7. Reproduction

See `PEP_Iterative_Training_Methodology.md` §Reproduction for exact commands. In short:
run a search with the stock TOML (baseline) and with a TOML adding
`IterativePepTraining = true` (iterative), then compare target-peptide counts at
PEP q-value < 0.01 from each run's `Task1SearchTask/AllPeptides.psmtsv`.
