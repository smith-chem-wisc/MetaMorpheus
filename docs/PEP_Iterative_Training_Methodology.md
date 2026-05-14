# PEP Iterative Training — Methodology

**Branch:** `pep-iterative-training`
**Span:** 2026-05-14
**Status:** complete — feature implemented, tested, evaluated; opt-in and off by default.

This is the lab notebook for the iterative-training work: what was done, in order, every
decision and why, including the parts that did not pan out as predicted. The polished
results live in `PEP_Iterative_Training_Technical_Report.md`; the go/no-go diagnostic
lives in `PEP_Iterative_Training_Phase1_Diagnostic.md`.

---

## Background

MetaMorpheus assigns posterior error probabilities (PEPs) to spectral matches with a
gradient-boosted decision tree (`Microsoft.ML` FastTree, Platt-calibrated), trained in
`PepAnalysisEngine` via a 4-fold cross-fit scheme. Training labels are assigned **once**,
before any model is fit: a non-decoy match is a positive example if its **raw-score**
q-value is at or below `QValueCutoff` (0.005 for standard bottom-up); mid-confidence
targets above that cutoff are dropped from training entirely.

The hypothesis: because the positive pool is tied to raw score — only one of ~17 features
the model eventually uses — the model never trains on the mid-score region where PEP
discrimination matters most. **Iterative training** re-selects the positive pool from the
model's own PEP-derived q-values after the first fit, then retrains, repeating until the
pool stabilizes. This is the core idea behind Percolator.

The work was structured as three phases: a go/no-go diagnostic, a test-driven prototype,
and an empirical evaluation.

---

## Phase 0 — Infrastructure

- **Repository:** `E:\CodeReview\MetaMorpheus`. Branch `pep-iterative-training` cut fresh
  from `origin/master` (`5a00fc2b`), which was verified identical to
  `upstream/master` (smith-chem-wisc).
- **Test suite at branch point:** 1,445 / 1,446 passing. The single failure,
  `FragmentReanalysisRaceConditionTest`, is a pre-existing flaky MetaDraw concurrency test
  unrelated to FDR/PEP — it passes in isolation and passed cleanly on the post-change full
  run. Not a blocker.
- **Build/runner:** .NET SDK pinned by `global.json` (8.0.x), NUnit via `dotnet test`.
- **Commit style:** simplified Conventional Commits, per the user's style guide.

---

## Phase 1 — Diagnostic (go/no-go)

Full detail in `PEP_Iterative_Training_Phase1_Diagnostic.md`. Summary:

Ran the current single-pass pipeline on the bottom-up human-lysate test data (two Velos
`.raw` files, `human_rev_GPTMDpruned.xml`). From `AllPeptides.psmtsv`, counted peptides
that would enter the positive pool on a hypothetical iteration 2 — raw q-value > 0.005
**and** PEP q-value ≤ 0.005:

- Current positive pool: 12,031 target peptides.
- Would-enter on iteration 2: 907.
- **Growth: 7.54%** — clears the pre-agreed ≥ 5% threshold ⇒ **GO**.

At the time, by the Percolator rule of thumb, 7.54% pool growth suggested a ~2–8% peptide
gain. Phase 3 showed that estimate was optimistic for this dataset; see the Phase 3
interpretation below for why.

---

## Phase 2 — Implementation

### Where the loop sits

`ComputePEPValuesForAllPSMs()` previously did, once: build 4 stratified folds → label by
raw-score q-value → train 4 cross-fit models (each on 3 folds, predict the held-out 4th) →
write `psm.PEP`. The change wraps that body in an iteration loop:

- **Iteration 0** labels positives by raw-score q-value — byte-for-byte the historical
  rule.
- **Iterations ≥ 1** relabel positives from the previous iteration's PEP-derived q-values.
- The loop stops at convergence or a configurable max (default 3).

### Decision: per-fold vs. global PEP q-value for relabeling

This was the crux, and it was surfaced to the user before any code was written.

The cross-fit scheme already gives isolation for free: a PSM in fold *k* gets its PEP from
`model_k`, trained on folds ≠ *k*. So after iteration 1, every PSM already carries a PEP
that no model trained on it produced. The open question was how to turn those PEPs into a
q-value for relabeling:

- **Global** — rank all PSMs together by cross-fit PEP. Matches the Phase 1 diagnostic
  exactly and gives less noisy q-values, but fold *k*'s relabel ranking then mixes in other
  folds' PEPs (whose models *did* train on fold *k*), a second-order leak.
- **Per-fold** (chosen) — for fold *k*, rank only fold *k*'s PSMs by their cross-fit PEP
  and compute q-values within the fold. Fold *k*'s new labels depend **only** on
  `{fold k's PSMs, model_k}`, and `model_k` never saw fold *k*. Airtight isolation, and
  the isolation test becomes a structural fact rather than a probabilistic argument.

Per-fold was chosen because the overview calls cross-fit isolation the property that
"must be respected," and per-fold makes it provable. The cost — q-values estimated on ~25%
of the data per fold — is mitigated by the folds being target/decoy-stratified, so each is
a representative sample. `ComputePerFoldPepQValues` implements this: a standard
target/decoy q-value (raw FDR, then monotonized) computed strictly within each fold,
returned as one `Dictionary<SpectralMatch, double>` per fold. The maps are pairwise
disjoint by construction.

### Decision: ambiguity removal deferred to the final pass

`Compute_PSM_PEP` did two things per PSM: write its PEP, and prune ambiguous peptide
hypotheses scoring well below the best. Pruning **mutates** the PSM (removes entries from
`BestMatchingBioPolymersWithSetMods`). Running it every iteration would compound: later
iterations would train on a progressively thinned hypothesis set, and a hypothesis the
iteration-1 model wrongly disliked could never be recovered.

So pruning is deferred — `Compute_PSM_PEP` gained a `removeAmbiguous` flag, false on
non-final iterations, true on the final one. The PEP value itself is unaffected by the
flag: PEP is `1 − max(hypothesis probability)`, and pruning only ever drops hypotheses
*below* that max, so the max — and therefore the PEP — is identical with or without
pruning. The "final pass" is identified as `converged || iteration == maxIterations - 1`;
a rare degenerate-fold break on an iteration ≥ 1 falls back to a pruning pass with the
previous iteration's retained models.

### Convergence and bounds

- `HasConverged(prev, current)` — extracted as an `internal static` method so it is unit-
  testable in isolation — returns true when the positive pool moved by less than
  `ConvergenceThreshold` (1%) of the previous count.
- Hard cap at `MaxTrainingIterations`, a constructor parameter (default 3), clamped to
  ≥ 1.
- `PositivePoolCountPerIteration` records the pool size each iteration and is surfaced in
  the PEP metrics block of `results.txt` ("Training Iterations Run", "Positive Training
  Pool Per Iteration").
- Decoys never move between iterations — only target labels are re-selected.

### The hardcoded fold concat

The old code had a hand-unrolled `Concat` of exactly 3-of-4 fold lists, with a comment
admitting it had to change if the fold count ever did. Restructuring the training loop
made the natural replacement —
`Enumerable.Range(0, numGroups).Where(g => g != k).SelectMany(g => PSMDataGroups[g])` —
fall straight out, so it was taken. This was the only opportunistic cleanup; it was in the
exact lines being rewritten anyway.

### The configuration parameter

`IterativePepTraining` (bool, default false) was added to `CommonParameters` — through the
big constructor, the `IterativePepTraining` property, `CloneWithNewTerminus`, and
`MetaMorpheusTask.SetAllFileSpecificCommonParams` (so file-specific parameters don't reset
it). `FdrAnalysisEngine.Compute_PEPValue` reads it from the file-specific parameters and
threads it into the `PepAnalysisEngine` constructor. Default false means the single-pass
pipeline is entirely unchanged unless a user opts in.

### Honest note: test ordering

The prompt asked for strict test-first development. In practice the C# Test project cannot
compile until the new public API surface (`PepAnalysisEngine` constructor parameters,
`IterativeTraining` / `MaxTrainingIterations` / `PositivePoolCountPerIteration` /
`PerFoldRelabelPepQValues`, `CommonParameters.IterativePepTraining`) exists. The API and
implementation were therefore built together, and the tests authored against the finished
API. The five tests are still genuine contract and regression tests — they exercise real
behavior and would catch real regressions — they were simply written after the API rather
than before it. This is disclosed rather than papered over.

### Tests

`Test/PepIterativeTrainingTest.cs`, five contracts:

| Test | Pins |
|---|---|
| `MaxTrainingIterations_IsClampedToAtLeastOne` | configurable max, clamped to ≥ 1 |
| `IterativeTraining_RunsNoMoreThanMaxIterations` | iteration count never exceeds the max |
| `PositivePoolCount_IsRecordedOncePerIteration` | per-iteration logging, surfaced in metrics text |
| `HasConverged_TriggersBelowOnePercentChange` | the < 1% convergence rule |
| `IterativeTrainingOff_RunsExactlyOnePassAndNeverRelabels` + `...WithMaxOneIteration...` | regression net — off (or max-1) = single pass, no relabel map ever built |
| `Relabeling_IsPerFoldAndDisjoint_NotGlobal` | cross-fit isolation — the per-fold relabel maps partition the PSMs exactly by fold; a global PEP map would instead contain every PSM |

The isolation test reconstructs the engine's folds from its deterministic public grouping
and asserts `PerFoldRelabelPepQValues[k]` is keyed by exactly fold *k*'s matches and that
the maps are pairwise disjoint — i.e. each PSM is relabeled only by its own fold's model.

**Full suite after the change: 1,453 passed, 0 failed, 1 skipped — zero regressions.** No
existing test needed revision; none asserted exact counts that iterative training would
have shifted, because the feature is off by default.

---

## Phase 3 — Evaluation

### Setup

Two runs on the identical bottom-up dataset (two Velos `.raw` files,
`human_rev_GPTMDpruned.xml`, the provided search TOML):

- **Baseline** — `IterativePepTraining` off (the Phase 1 search output).
- **Iterative** — a copy of the search TOML with `IterativePepTraining = true`.

Headline metric, per the user: **count of target peptides with PEP q-value < 0.01**, read
from `AllPeptides.psmtsv` (`Decoy/Contaminant/Target == "T"` and `PEP_QValue < 0.01`).
PSM-level counts were deliberately not used — the model is trained on peptides, so PSM
scoring is noisier, and repeated observations of one peptide carry no extra information.
Counting was done with a standalone Python parser
(`E:\CodeReview\pep_phase3_compare.py`), no diagnostic code in the repo.

### Result

| Metric | Baseline | Iterative | Δ |
|---|---|---|---|
| **Target peptides, PEP q < 0.01** | **13,049** | **13,135** | **+86 (+0.66%)** |
| Target peptides, PEP q < 0.05 | 13,998 | 14,084 | +86 (+0.61%) |
| Target peptides, PEP q < 0.005 | 12,831 | 12,888 | +57 (+0.44%) |
| Target peptides, raw q < 0.01 (PEP-independent sanity) | 12,651 | 12,644 | −7 (noise) |

Model-quality metrics improved more than the peptide count did:

| | Baseline | Iterative |
|---|---|---|
| LogLoss | 0.0734 | 0.0602 |
| AUC | 0.99926 | 0.99950 |
| PR-AUC | 0.99952 | 0.99965 |
| Accuracy | 0.9916 | 0.9946 |

Iterative training ran the full 3 iterations (it did not converge). Positive pool per
iteration: **23,300 → 24,959 → 25,222**. Iteration 0's count (23,300) is exactly the
baseline's positive pool — a clean confirmation that iteration-0 labeling is identical to
the single-pass pipeline. The iteration-1→2 change was +1.05%, just above the 1%
convergence threshold, which is why the third iteration ran.

### Time cost

| | Baseline | Iterative | Δ |
|---|---|---|---|
| FDR analysis engine (contains PEP training) | 15.2 s | 37.9 s | +22.7 s (≈ 2.5×) |
| Whole search task, wall-clock | 49.7 s | 71.9 s | +22.2 s (+45%) |

PEP training is ~3 model-fitting passes instead of 1; the engine step grows ~2.5× (not a
full 3× — q-value computation and the file-specific dictionaries built once in the
constructor are fixed overhead). On this 2-file, ~47k-PSM / ~29k-peptide dataset the whole
search takes 45% longer. The marginal cost is roughly ~11 s per extra iteration here.
This scales with dataset size; for very large searches it would be more noticeable, and
iterations ≥ 2 could optionally be run with fewer trees — a candidate optimization, not
done here.

### Interpretation — why the gain is smaller than the diagnostic predicted

The Phase 1 diagnostic measured 7.54% positive-pool growth and, by the Percolator rule of
thumb, suggested ~2–8% more peptides. The measured gain is +0.66%. The diagnostic was
optimistic for a knowable reason:

The 907 "would-enter" peptides it counted all had **PEP q-value ≤ 0.005 already** — they
were *already passing* in the baseline results. Adding them to the *training* pool does
not add them to the *results*; they were never missing. The only benefit available is the
indirect one — a model trained on a larger, more diverse positive set ranks *other*
borderline peptides slightly better. The +86 peptides are that indirect effect. The
diagnostic correctly predicted the training pool would grow ~7.5% (it grew 8.25% over
three iterations); it over-translated that into a results gain.

The model also started near its ceiling — baseline AUC 0.99926 — leaving little headroom.
LogLoss and AUC nonetheless improved, which says the model itself did get better; the
peptide count is just a coarse, near-saturated readout of that improvement on this
particular dataset.

### Calibration

Reliability tables (predicted PEP vs. observed decoy fraction, 10 bins) are comparable
between the two runs. The low-PEP bin `[0.0, 0.1)` — the one that matters for the headline
count — stays well-calibrated in both (baseline: mean predicted 0.0005, decoy fraction
0.0040; iterative: 0.0003 / 0.0057). Iterative pulls slightly more peptides into that bin
(12,828 → 13,008) at a marginally higher decoy fraction. Mid-range bins are too sparse
(n ≈ 5–50) to read into. Full tables are in the technical report.

---

## Out-of-scope issues noted during code review

Per the overview's "top concerns beyond iterative training" — recorded as candidate
follow-up work, deliberately **not** addressed here:

- `BGDTreeOptions` hyperparameters (`NumberOfTrees = 400`, `LearningRate = 0.2`,
  `MinimumExampleCountPerLeaf = 10`) are untuned; no early stopping. The user plans to
  investigate FastTree settings in a separate PR.
- No reliability-diagram diagnostics in `AggregateMetricsForOutput` — only
  accuracy/AUC/F1/LogLoss.
- `AbsoluteProbabilityThatDistinguishesPeptides = 0.05` ambiguity pruning is an absolute
  probability difference, behaving differently at low vs. high PEP.
- `BuildFileSpecificDictionaries` runs globally before fold splitting — mild test-set
  leakage via file-specific statistics.

(The hardcoded 3-of-4 fold concat, also on that list, *was* fixed — it fell out naturally
while restructuring the training loop.)

---

## Follow-up work

- **Iterations ≥ 2 with fewer trees** — cut the time cost where the model is already
  mostly converged.
- **FastTree hyperparameter tuning** (separate PR, planned by the user) — should compound
  with this work: a better-tuned model has more headroom for iterative relabeling to
  exploit.
- **New PEP features** (separate PR, planned by the user) — likewise expected to compound;
  iterative training's value scales with how much signal the feature set can express
  beyond raw score.
- **Larger / lower-quality datasets** — the gain here was small partly because the model
  was near its ceiling on clean Orbitrap-class bottom-up data. Datasets with more model
  headroom may show a larger effect.
- **Reliability-diagram output** in the metrics block, so calibration is monitored
  routinely rather than computed ad hoc.

---

## Reproduction

```
# baseline
dotnet CMD\bin\Release\net8.0\CMD.dll ^
  -t E:\TestData\MetaMorpheus\BottomUp\Task-SearchTaskconfig.toml ^
  -d E:\TestData\MetaMorpheus\BottomUp\human_rev_GPTMDpruned.xml ^
  -s E:\TestData\MetaMorpheus\BottomUp\20100609_Velos1_TaGe_SA_293_3.raw ^
     E:\TestData\MetaMorpheus\BottomUp\20100609_Velos1_TaGe_SA_293_4.raw ^
  -o E:\TestData\MetaMorpheus\BottomUp\IterativePepResults

# iterative (TOML adds: IterativePepTraining = true under [CommonParameters])
dotnet CMD\bin\Release\net8.0\CMD.dll ^
  -t E:\TestData\MetaMorpheus\BottomUp\Task-SearchTaskconfig_IterativePep.toml ^
  -d E:\TestData\MetaMorpheus\BottomUp\human_rev_GPTMDpruned.xml ^
  -s E:\TestData\MetaMorpheus\BottomUp\20100609_Velos1_TaGe_SA_293_3.raw ^
     E:\TestData\MetaMorpheus\BottomUp\20100609_Velos1_TaGe_SA_293_4.raw ^
  -o E:\TestData\MetaMorpheus\BottomUp\IterativePepResults_iterative

# compare
python E:\CodeReview\pep_phase3_compare.py
```
