# PEP Iterative Training — Phase 1 Diagnostic

**Date:** 2026-05-14
**Branch:** `pep-iterative-training`
**Decision:** **GO** — proceed to Phase 2 (7.54% positive-pool growth, rule threshold ≥ 5%).

## Question

Does iterative training have meaningful upside on real MetaMorpheus data? Concretely:
how many peptides currently *excluded* from the PEP training pool (mid-confidence by raw
search score) would be pulled into the positive pool on a hypothetical iteration 2,
because the trained model's PEP-derived q-value already ranks them as confident?

## Method

1. Ran the current (single-pass) MetaMorpheus pipeline on the bottom-up human lysate test
   data — two Velos `.raw` files, `human_rev_GPTMDpruned.xml` database, the provided
   `Task-SearchTaskconfig.toml` (Classic search, trypsin). Output:
   `E:\TestData\MetaMorpheus\BottomUp\IterativePepResults\`.
2. Parsed `Task1SearchTask\AllPeptides.psmtsv` — peptide-level, which is the training
   grouping for this dataset (≥ 100 distinct full sequences ⇒ `UsePeptideLevelQValueForTraining = true`).
3. Counted peptide rows satisfying the **would-enter** condition:
   raw-score `QValue` > 0.005 **AND** `PEP_QValue` ≤ 0.005.
4. Expressed that count as a percentage of the **current positive pool**:
   target rows with raw-score `QValue` ≤ 0.005.

`QValueCutoff` = 0.005 confirmed from code (`PepAnalysisEngine` ctor:
`Math.Max(min QValueCutoffForPepCalculation, 0.005)` for "standard" search) and from the
search TOML (`QValueCutoffForPepCalculation = 0.005`).

The count was done with a standalone Python parser
(`E:\CodeReview\pep_phase1_diagnostic.py`) — no diagnostic code added to the MetaMorpheus
codebase.

## Results

| Metric | Value |
|---|---|
| Total peptide rows (`AllPeptides.psmtsv`) | 29,317 |
| — target / decoy | 21,464 / 7,853 (no contaminants) |
| Current positive pool (target, raw q ≤ 0.005) | 12,031 |
| Would-enter on iteration 2 (raw q > 0.005 AND PEP_q ≤ 0.005) | 907 |
| **Growth as % of positive pool** | **7.54%** |

Would-enter PEP distribution (model's own confidence in these "rescued" peptides):

| Stat | Value |
|---|---|
| min / median / max | 0.0000 / 0.0000 / 0.2365 |
| mean | 0.0110 |

| PEP bin | Count |
|---|---|
| [0.00, 0.05) | 843 |
| [0.05, 0.10) | 28 |
| [0.10, 0.20) | 31 |
| [0.20, 0.30) | 5 |
| ≥ 0.30 | 0 |

## Interpretation

7.54% clears the pre-agreed ≥ 5% threshold, so the work proceeds to Phase 2.

The would-enter set is not marginal noise: 843 of 907 peptides (93%) have PEP < 0.05, and
the maximum PEP in the whole set is 0.236. These are matches the trained model considers
genuinely confident, but which raw search score alone placed in the mid-confidence band
and the single-pass pipeline therefore dropped from training entirely. This is precisely
the population the iterative-training hypothesis predicts exists — peptides that look
target-like across the full ~17-feature space but had unremarkable raw scores.

A 7.54% positive-pool expansion, by the Percolator-literature rule of thumb cited in the
overview, corresponds to roughly a 2–8% peptide gain at 1% FDR — enough to be worth a
proper prototype.

### Note on the "23,300" figure in `results.txt`

`results.txt` reports "Targets Used for Training: 23,300", larger than the 12,031 positive
pool above. The two are measured at different granularities and are *not* in conflict:

- 12,031 = unique target peptide rows with q ≤ 0.005, counted from `AllPeptides.psmtsv`.
- 23,300 = training rows emitted by `PepAnalysisEngine.CreatePsmData`, which adds one
  `PsmData` row *per ambiguous bio-polymer hypothesis per peptide*. Ambiguous peptides
  contribute multiple rows.

The diagnostic deliberately keeps numerator and denominator at the **same** peptide-row
granularity (907 / 12,031). Dividing the unique-row numerator by the ambiguity-inflated
23,300 would mix granularities and understate the growth to ~3.9%.

## Reproduction

```
# baseline search
dotnet CMD\bin\Release\net8.0\CMD.dll ^
  -t E:\TestData\MetaMorpheus\BottomUp\Task-SearchTaskconfig.toml ^
  -d E:\TestData\MetaMorpheus\BottomUp\human_rev_GPTMDpruned.xml ^
  -s E:\TestData\MetaMorpheus\BottomUp\20100609_Velos1_TaGe_SA_293_3.raw ^
     E:\TestData\MetaMorpheus\BottomUp\20100609_Velos1_TaGe_SA_293_4.raw ^
  -o E:\TestData\MetaMorpheus\BottomUp\IterativePepResults

# count
python E:\CodeReview\pep_phase1_diagnostic.py
```

Baseline search wall-clock: ~50 s. Test suite state at branch point: 1,445/1,446 passing;
the single failure (`FragmentReanalysisRaceConditionTest`) is a pre-existing flaky MetaDraw
concurrency test, unrelated to FDR/PEP, and passes in isolation.

## Out-of-scope issues noted during code review

Per the overview's "top concerns beyond iterative training" — recorded here as candidate
follow-up work, **not** addressed in this PR unless they directly block the implementation:

- `BGDTreeOptions` hyperparameters (`NumberOfTrees = 400`, `LearningRate = 0.2`,
  `MinimumExampleCountPerLeaf = 10`) are untuned; no early stopping.
- No calibration diagnostics (reliability diagram) in `AggregateMetricsForOutput` — only
  accuracy/AUC/F1/LogLoss.
- The 3-of-4 fold concat in `ComputePEPValuesForAllPSMs` is hardcoded, blocking flexible
  fold counts.
- `AbsoluteProbabilityThatDistinguishesPeptides = 0.05` ambiguity pruning is an absolute
  probability difference, behaving very differently at low vs high PEP.
- `BuildFileSpecificDictionaries` runs globally before fold splitting — mild test-set
  leakage via file-specific statistics.
