# PR Summary — Tune the FastTree PEP learning rate (0.2 → 0.05)

## Summary

The PEP model's FastTree hyper-parameters have never been tuned and were hardcoded. A
three-round hyper-parameter sweep across three dataset sizes shows the learning rate is the
one knob that matters — and the current default of **0.2 is too aggressive**: it leaves the
model overconfident (it can produce infinite LogLoss) and, worst of all, **iterative PEP
training at LR 0.2 actively degrades calibration**.

This PR lowers the default learning rate **0.2 → 0.05**. It is a **calibration fix**, not
an identification-count change: on realistically sized data the peptide count at PEP
q‑value < 0.01 is unchanged-to-slightly-better, while LogLoss drops by roughly half to a
tenth. It also makes the FastTree options injectable so they can be tuned and tested.

## The change

`EngineLayer/FdrAnalysis/PEPAnalysisEngine.cs`:
- `BGDTreeOptions` default `LearningRate`: **0.2 → 0.05** (`NumberOfTrees` 400, `NumberOfLeaves`
  20, `MinimumExampleCountPerLeaf` 10 unchanged — the sweep showed those flat).
- `BGDTreeOptions` is now backed by an optional `customTreeOptions` constructor parameter
  (same pattern as the existing optional `rtPredictor`), so the hyper-parameters can be
  injected for tuning and tests; the fixed plumbing fields (column names, threads, seeds)
  are always enforced regardless.

## Evidence

### The decisive picture — learning rate across three dataset sizes

`lr=0.05` vs the old `lr=0.2`, target peptides at PEP q < 0.01 and model LogLoss:

| dataset | size | count effect | calibration (LogLoss) |
|---|---|---|---|
| tiny test snippet | 294 MS2 | −21 to −32% (underfits) | weak model regardless |
| **BottomUp 2‑min snips** | **~1,000 MS2** | **+1.6% (best of configs)** | **0.067 → 0.011** |
| **full BottomUp (4 raws)** | **~32,000 MS2** | **neutral (−0.1%)** | **0.048 → 0.033 (off), 0.068 → 0.023 (on)** |

The peptide *count* is essentially insensitive to learning rate on any realistically sized
search — the only place `lr=0.05` loses count is a 294‑scan **toy** dataset, where the
gentler rate underfits. No real search is that small. Meanwhile **calibration improves
substantially and consistently** from the medium dataset upward.

### Full-data validation (4 raw files, realistic search with modifications)

| run | target peptides, PEP q<0.01 | LogLoss | AUC |
|---|---:|---:|---:|
| LR 0.20, iterative off | 18,301 | 0.0482 | 0.99942 |
| LR 0.20, iterative on  | 18,355 | 0.0680 | 0.99748 |
| LR 0.05, iterative off | 18,285 | 0.0333 | 0.99955 |
| LR 0.05, iterative on  | 18,334 | 0.0233 | 0.99970 |

Count difference is noise (±0.1%); calibration improvement is large and real. Note the
worst cell in the table is **LR 0.2 + iterative training** (LogLoss 0.068) — the old
default fights the iterative-training feature. LR 0.05 + iterative is the best (0.023).

### Why the sweep was done three times

- **Round 1** (full data, one-factor-at-a-time): peptide count is flat for every FastTree
  knob; only learning rate moves calibration.
- **Round 2** (full data, LR × iterative on/off): `lr=0.05` is the best-calibrated and the
  old `lr=0.2` degrades under iterative retraining.
- **Round 3** (medium snips, after the regression suite revealed `lr=0.05` underfits the
  294‑scan toy test data): confirms the underfitting is confined to toy-sized data — at
  ~1,000 MS2 scans `lr=0.05` is already the best config on both count and calibration.

Full lab notebook: `docs/PEP_FastTree_Options_Methodology.md`.

## Tests

**New, in CI:**
- `Test/FastTreeLearningRateTest.cs` — integration test on the four BottomUp 2‑minute snips
  (newly added, see below). Runs PEP at LR 0.2 and LR 0.05 and asserts the tuned rate is
  better calibrated (lower LogLoss) at no cost in target peptides. Demonstrates the change
  with the high-vs-low learning rate result baked into CI.
- `Test/FastTreeOptionsTest.cs` — pins the tuned default (`LearningRate == 0.05`) and the
  custom-options injection seam.

**New test data:** `Test/TestData/FastTreeSnips/` — four ~2‑minute `.mzML` snips cut from
the retention-time heart of the BottomUp runs (~13 MB total) plus a 0.7 MB protein
database pruned to just the proteins those snips identify. This medium-sized dataset
(~1,000 MS2 scans) is the size point between the existing tiny test snippet and a full run.

**Updated brittle tests:** `Test/PostSearchAnalysisTaskTests.cs` — the two
`AllResultsAndResultsTxtContainsCorrectValues_*_BottomUp` tests assert exact PSM/peptide
counts on the 294‑scan toy test dataset. Because that dataset is small enough for the
gentler learning rate to underfit, their numbers shifted (e.g. pep‑q‑value PSMs 382 → 300).
The expected values were updated to the retuned model's output on that toy dataset, with a
comment in each test explaining why and pointing at the methodology doc. This is the
expected, documented consequence of the model change on toy-sized data — realistically
sized searches are unaffected (see the evidence above).

## Test plan

- [x] `FastTreeLearningRateTest`, `FastTreeOptionsTest` — pass
- [x] `PostSearchAnalysisTaskTests` (the two updated tests) — pass
- [x] Full NUnit suite — no regressions attributable to this change
- [x] End-to-end CMD searches on the full 4-raw BottomUp dataset, LR 0.2 vs 0.05, iterative
      off and on

## Caveats

- Single instrument class (Velos bottom-up) and a fixed training seed. The count deltas
  are within seed jitter; the LogLoss improvement is large, monotonic in learning rate,
  and reproduced across every round and the full-data CMD validation.
- `lr=0.05` underfits genuinely tiny (sub‑300‑MS2) datasets. This is below any realistic
  search size and is the reason the two `PostSearchAnalysisTaskTests` numbers moved; it is
  documented rather than hidden.

## Files changed

- `EngineLayer/FdrAnalysis/PEPAnalysisEngine.cs` — LR default 0.2 → 0.05; injectable `BGDTreeOptions`
- `Test/FastTreeLearningRateTest.cs` — new integration test (high vs low LR on the snips)
- `Test/FastTreeOptionsTest.cs` — new unit tests (default value, injection seam)
- `Test/PostSearchAnalysisTaskTests.cs` — updated expected values on the toy dataset + explanatory comments
- `Test/TestData/FastTreeSnips/` — four 2‑minute snips + a pruned protein database
- `Test/Test.csproj` — registers the new test data
- `docs/PEP_FastTree_Options_Methodology.md`, `docs/PEP_FastTree_Options_PR_Summary.md` — write-up
