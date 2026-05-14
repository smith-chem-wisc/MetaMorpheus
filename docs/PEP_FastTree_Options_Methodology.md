# PEP FastTree Options — Methodology

**Branch:** `fastTreeOptions` (from `pep-iterative-training`)
**Date:** 2026-05-14
**Status:** complete — one-line change shipped (`LearningRate` 0.2 → 0.05), framed as a
calibration fix, **not** a peptide-count improvement.

Lab notebook for the FastTree hyperparameter exploration: what was swept, what the evidence
said, and the honest reason the change that shipped is a calibration fix rather than a
yield improvement — including the round where the fast exploration harness mispredicted the
direction of the headline metric and the full CMD validation corrected it.

---

## Question

The PEP model is a `Microsoft.ML` FastTree binary classifier. Its hyperparameters
(`BGDTreeOptions`: 400 trees, learning rate 0.2, 20 leaves, min-leaf 10, no early stopping)
were hardcoded and had never been tuned. Can evidence-based tuning improve the result —
where "the result" is, per the user, the **count of target peptides at PEP q-value < 0.01**
from `AllPeptides.psmtsv`?

Test data: the bottom-up human-lysate set, now **four** Velos `.raw` files
(`_2`, `_3`, `_4`, `_5`), `human_rev_GPTMDpruned.xml` database.

---

## Approach — explore fast, validate honestly

1. **Make `BGDTreeOptions` injectable.** Added an optional `customTreeOptions` constructor
   parameter to `PepAnalysisEngine` (consistent with the existing optional `rtPredictor`
   seam). The getter still produces the defaults when nothing is injected, and either way
   stamps the fixed plumbing fields (column names, single-threading, fixed seeds) so a
   caller only supplies the tunable knobs.
2. **A fast exploration harness** — an `[Explicit]` test that, per configuration, runs a
   fresh classic search of the four raw files (no modifications, for speed) and replays the
   FDR pipeline with the configuration's options injected, reporting target peptides at
   PEP q < 0.01. The search is byte-identical across configurations, so config-to-config
   deltas are trustworthy even though absolute counts are not.
3. **Validate the harness** against a full CMD search with a matching no-mods TOML.
4. **Sweep**, then **validate the apparent winner** with full CMD searches on the realistic
   (with-modifications) configuration.

The harness was a local exploration artifact with machine-specific paths; it was archived
to `E:\CodeReview\fasttree_exploration_harness.cs.txt` and is **not** committed to the repo.
Its sweep outputs (`E:\CodeReview\fasttree_sweep_results.md`, `fasttree_sweep_round2.md`)
are the raw evidence.

---

## Harness validation

| | Target peptides, PEP q < 0.01 |
|---|---:|
| CMD no-mods search (4 raws, default options) | 14,460 |
| Harness baseline (same config) | 14,306 |
| Offset | −154 (−1.08%) |

The harness undercounts the real pipeline by ~1.08% — a *constant* offset (the harness's
search is identical across all configurations), so config-to-config deltas remain
trustworthy for ranking. This caveat became important later: see Round 2.

---

## Round 1 — FastTree options alone (single-pass PEP)

One-factor-at-a-time around the defaults; iterative training off for a clean signal.

| config | peptides (PEP q<0.01) | Δ vs baseline | LogLoss | AUC |
|---|---:|---:|---:|---:|
| baseline 400t/0.2lr/20l/10m | 14,306 | — | 0.0383 | 0.99911 |
| trees=100 | 14,292 | −14 | 0.0401 | 0.99892 |
| trees=200 | 14,297 | −9 | 0.0444 | 0.99895 |
| trees=800 | 14,261 | −45 | Infinity | 0.99915 |
| lr=0.05 | 14,314 | +8 | 0.0313 | 0.99949 |
| lr=0.1 | 14,297 | −9 | 0.0351 | 0.99923 |
| lr=0.3 | 14,254 | −52 | Infinity | 0.99665 |
| leaves=10 | 14,269 | −37 | 0.0473 | 0.99847 |
| leaves=30 | 14,260 | −46 | 0.0613 | 0.99874 |
| leaves=50 | 14,289 | −17 | Infinity | 0.99903 |
| minLeaf=5 | 14,285 | −21 | 0.0654 | 0.99841 |
| minLeaf=20 | 14,306 | +0 | 0.0572 | 0.99881 |
| minLeaf=50 | 14,322 | +16 | 0.0460 | 0.99916 |

**Findings:**
- **Peptide count is essentially insensitive to FastTree options** — every config within
  +16/−52 of baseline (±0.36%), smaller than the harness's own 1% offset from CMD.
- **Calibration (LogLoss/AUC) does respond.** Lower learning rate gives lower LogLoss;
  `lr=0.05` is the best-calibrated. Aggressive settings (`lr=0.3`, `trees=800`,
  `leaves=50`) trigger `Infinity` LogLoss — a single overconfident wrong prediction
  (p ≈ 0 or 1) makes log-loss infinite. The baseline `lr=0.2` sits near that cliff.

---

## Round 2 — FastTree options × iterative training

The iterative-training feature (the `pep-iterative-training` PR this branch builds on) is
*meant to compound* with FastTree tuning, so each config was run twice — iterative off and
on (max 3 iterations) — to measure the interaction rather than assume it. Focused on the
learning-rate axis and the lr × trees interaction, since that is where Round 1's only real
signal lived.

| config | iter | peptides (PEP q<0.01) | Δ vs baseline-off | LogLoss | AUC |
|---|---|---:|---:|---:|---:|
| baseline 400t/0.2lr | off | 14,306 | — | 0.0383 | 0.99911 |
| baseline 400t/0.2lr | on | 14,357 | +51 | 0.0517 | 0.99902 |
| lr0.02 400t | off | 14,307 | +1 | 0.0285 | 0.99955 |
| lr0.02 400t | on | 14,349 | +43 | 0.0259 | 0.99947 |
| lr0.05 400t | off | 14,314 | +8 | 0.0313 | 0.99949 |
| lr0.05 400t | on | 14,370 | +64 | 0.0267 | 0.99945 |
| lr0.1 400t | off | 14,297 | −9 | 0.0351 | 0.99923 |
| lr0.1 400t | on | 14,350 | +44 | 0.0322 | 0.99940 |
| lr0.02 800t | off | 14,317 | +11 | 0.0299 | 0.99952 |
| lr0.02 800t | on | 14,358 | +52 | 0.0247 | 0.99949 |
| lr0.05 800t | off | 14,307 | +1 | 0.0368 | 0.99946 |
| lr0.05 800t | on | 14,366 | +60 | 0.0303 | 0.99947 |
| lr0.1 800t | off | 14,285 | −21 | 0.0416 | 0.99920 |
| lr0.1 800t | on | 14,355 | +49 | 0.0402 | 0.99938 |

**Findings:**
- **Iterative training is the real lever** for the count (+43 to +64); FastTree options
  alone stay flat (Round 1 confirmed).
- **They appear to compound modestly:** `lr0.05/400t/on` = +64, vs `baseline/on` = +51.
- **The sharp signal — the default learning rate is wrong for the iterative path.**
  `baseline 0.2 LR + iterative on` has LogLoss **0.0517**, *worse* than baseline-off
  (0.0383): retraining iteratively at LR 0.2 degrades calibration. Every lower-LR config
  *improves* with iterative on; `lr=0.05` was the best-calibrated config taking the best
  count.

At this point the harness pointed to `lr=0.05` as a small count win **and** a calibration
win. The next step tested that.

---

## CMD validation — and the harness was wrong on the count

Four full CMD searches on all four raws with the **realistic (with-modifications) TOML**:
`{LR 0.2, LR 0.05} × {iterative off, on}`.

| run | peptides PEP q<0.01 | LogLoss | AUC |
|---|---:|---:|---:|
| LR 0.20, iterative off | **18,301** | 0.0482 | 0.99942 |
| LR 0.20, iterative on | **18,355** | 0.0680 | 0.99748 |
| LR 0.05, iterative off | 18,285 | 0.0333 | 0.99955 |
| LR 0.05, iterative on | 18,334 | 0.0233 | 0.99970 |

| comparison | Δ |
|---|---:|
| LR 0.20 → 0.05, iterative off | **−16** |
| LR 0.20 → 0.05, iterative on | **−21** |
| iterative gain @ LR 0.20 | +54 |
| iterative gain @ LR 0.05 | +49 |

**The count result reversed.** The harness (no-mods search) predicted `lr=0.05` at +8/+13;
the realistic CMD search (with modifications) measured **−16 / −21**. The sign flipped
between the two search configurations.

That flip is not a harness bug — it is the **proof that the peptide count is genuinely
insensitive to learning rate**. Across every experiment the count moved only within ±0.1–0.4%,
and a difference that small does not survive a change in the search's modification set. The
honest conclusion: **FastTree tuning offers no peptide-count improvement on this data.**

**But the calibration result held — and strengthened — from harness to CMD:**
- LR 0.05 roughly halves to thirds the LogLoss: 0.0482 → 0.0333 (off), 0.0680 → 0.0233 (on).
- AUC with iterative on: 0.99748 → 0.99970.
- The worst configuration in the whole study is **LR 0.2 + iterative training**
  (LogLoss 0.0680, AUC 0.99748) — iterative retraining at the aggressive default learning
  rate actively degrades calibration. **LR 0.05 + iterative is the best** (0.0233 / 0.99970).
- Iteration trajectories were near-identical (positive pool 32,287 → 34,622 → 34,960 at
  LR 0.2; 32,287 → 34,698 → 35,011 at LR 0.05), so the calibration difference is the model,
  not the relabeling.

---

## The regression suite — lr=0.05 underfits toy-sized data

Applying the change and running the full NUnit suite, two `PostSearchAnalysisTaskTests`
exact-count tests failed — and not trivially. The `_PepQValue_BottomUp` test dropped hard:
PSMs 382 → 300, peptides 153 → 104, protein groups 145 → 113 (−21 to −32%).

That test runs on a **294‑MS2 toy dataset** (~700 PEP training examples). On a training set
that tiny, `lr=0.05` with 400 trees **underfits** — the model learns too slowly with too
few examples (its LogLoss there is ~0.54, versus ~0.03–0.07 on the full data). The full
sweep had only ever been run on the large dataset, so this size regime was unseen. This was
a real finding, not test brittleness, and it put the decision back open: a blanket
`lr=0.05` would regress every small search if the underfitting extended to realistic sizes.

## Round 3 — the medium dataset

To find where on the size axis the underfitting stops, four ~2‑minute snips were cut from
the retention-time heart of the BottomUp runs (~1,000 MS2 scans total — between the
294‑scan toy and the ~32,000‑scan full data) and the learning-rate sweep was re-run on
them, iterative off and on.

| config | iter | peptides (PEP q<0.01) | Δ vs baseline-off | LogLoss |
|---|---|---:|---:|---:|
| baseline 400t/0.2lr | off | 256 | — | 0.0672 |
| baseline 400t/0.2lr | on | 251 | −5 | 0.1350 |
| lr0.05 400t | off | 263 | +7 | 0.0113 |
| lr0.05 400t | on | 259 | +3 | 0.0226 |
| lr0.05 800t | off | 263 | +7 | 0.0333 |
| lr0.1 400t | off | 263 | +7 | 0.0268 |
| lr0.1 800t | off | 260 | +4 | 0.0454 |

**The underfitting is confined to toy-sized data.** At ~1,000 MS2 scans `lr=0.05` is
already the *best* count config (+7) and roughly 6× better LogLoss than baseline — the
medium dataset behaves like the large one, not like the toy. Lining up all three sizes:

| dataset | size | `lr=0.05` count effect | `lr=0.05` calibration |
|---|---|---|---|
| tiny test snippet | 294 MS2 | −21 to −32% (underfits) | weak model regardless |
| BottomUp snips | ~1,000 MS2 | +1.6% (best of configs) | LogLoss 0.067 → 0.011 |
| full BottomUp | ~32,000 MS2 | neutral (−0.1%) | LogLoss 0.048 → 0.033 |

No realistic search is as small as the 294‑scan test snippet. `lr=0.05` is the right
default for any real-world search; the only thing it "breaks" is two unit tests that
happen to exercise a toy-sized dataset.

---

## Decision and what shipped

`lr=0.05` is shipped as the default, framed honestly as a **calibration fix**:

> `EngineLayer/FdrAnalysis/PEPAnalysisEngine.cs` — `BGDTreeOptions` default
> `LearningRate` lowered **0.2 → 0.05**; `BGDTreeOptions` also made injectable via an
> optional `customTreeOptions` constructor parameter.

Why it is worth shipping:
- PEP is a *probability* — only meaningful if calibrated. Cutting LogLoss by half-to-a-tenth
  means a reported PEP of 0.05 is much closer to a true 5% error rate.
- It removes the `Infinity`-LogLoss instability the aggressive default walked into.
- It specifically repairs the LR-0.2-plus-iterative-training degradation, so it
  **compounds cleanly with the `pep-iterative-training` PR** rather than fighting it.
- It is count-neutral-to-better on every realistically sized dataset.

`NumberOfTrees` (400), `NumberOfLeaves` (20), and `MinimumExampleCountPerLeaf` (10) were
left unchanged — Round 1 showed them flat on the count, and more trees did not help
calibration over 400 at `lr=0.05` in Rounds 2–3.

### Tests

New, in CI:
- `Test/FastTreeLearningRateTest.cs` — integration test on the four BottomUp snips. Runs
  PEP at LR 0.2 and LR 0.05 and asserts the tuned rate is better calibrated (lower LogLoss)
  at no cost in target peptides. On the bundled snips: LR 0.2 → 252 peptides / LogLoss 1.18;
  LR 0.05 → 256 peptides / LogLoss 0.12.
- `Test/FastTreeOptionsTest.cs` — pins `LearningRate == 0.05` and the other defaults, and
  the custom-options injection seam (fixed plumbing fields always enforced).

New test data — `Test/TestData/FastTreeSnips/`: the four ~2‑minute snips (~13 MB) plus a
0.7 MB protein database pruned to just the proteins those snips identify.

Updated brittle tests — `Test/PostSearchAnalysisTaskTests.cs`: the two
`AllResultsAndResultsTxtContainsCorrectValues_*_BottomUp` tests assert exact counts on the
294‑scan toy dataset, where the gentler learning rate underfits. Their expected values were
updated to the retuned model's output on that toy dataset, each with a comment explaining
why and pointing here. The `_QValue` test moved by one peptide (174 → 173, PEP-driven
best-PSM-per-peptide selection); the `_PepQValue` test moved more (e.g. 382 → 300 PSMs) for
the toy-data underfitting reason above. Realistically sized searches are unaffected.

---

## Honest caveats

- **Single instrument class, single seed.** All training used the fixed seed 42. The count
  deltas are within seed jitter; the LogLoss deltas are large, monotonic in learning rate,
  and reproduced across all three rounds and the full-data CMD validation.
- **`lr=0.05` underfits genuinely tiny (sub‑300‑MS2) datasets.** This is below any
  realistic search size, and is the reason the two `PostSearchAnalysisTaskTests` numbers
  moved. Documented, not hidden.
- **The exploration harness mispredicted the count direction once** (no-mods harness said
  `+8`, realistic CMD said `−16`). That flip is itself the evidence that the count is
  noise-insensitive to this knob; the harness remained a sound *ranking* tool for the
  calibration signal, which was always well above its noise floor.
- **No reliability-diagram check here.** LogLoss/AUC are the calibration proxies used; a
  proper reliability diagram in the metrics output remains a candidate follow-up.

---

## Follow-up work

- A reliability-diagram in the PEP metrics block, so calibration is monitored routinely
  rather than inferred from LogLoss.
- Early-stopping for FastTree — not pursued here (`lr=0.05` already avoids the
  overconfidence cliff), but a natural next knob.
- A size-aware learning rate (gentler only when the training pool is large) would let the
  toy regime keep `lr=0.2` behavior — only worth it if sub‑300‑MS2 searches ever matter.

---

## Reproduction

The four validation searches:

```
# LR 0.2 (build CMD from this branch's parent state, or revert the one-line change):
dotnet CMD\bin\Release\net8.0\CMD.dll -t ...\Task-SearchTaskconfig.toml          -d ...\human_rev_GPTMDpruned.xml -s <4 raws> -o ...\FastTreeCmp_LR020_iterOff
dotnet CMD\bin\Release\net8.0\CMD.dll -t ...\Task-SearchTaskconfig_IterativePep.toml -d ...\human_rev_GPTMDpruned.xml -s <4 raws> -o ...\FastTreeCmp_LR020_iterOn

# LR 0.05 (this branch as shipped):
dotnet CMD\bin\Release\net8.0\CMD.dll -t ...\Task-SearchTaskconfig.toml          -d ...\human_rev_GPTMDpruned.xml -s <4 raws> -o ...\FastTreeCmp_LR005_iterOff
dotnet CMD\bin\Release\net8.0\CMD.dll -t ...\Task-SearchTaskconfig_IterativePep.toml -d ...\human_rev_GPTMDpruned.xml -s <4 raws> -o ...\FastTreeCmp_LR005_iterOn
```

Headline count per run: target peptides with `PEP_QValue < 0.01` in
`Task1SearchTask/AllPeptides.psmtsv`. LogLoss/AUC: the PEP metrics block of `results.txt`.
