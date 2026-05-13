# PR Review Triage for smith-chem-wisc/MetaMorpheus#2649

This document summarizes my assessment of the review comments on the PR associated with the current branch, `parsimon_fixes`.

## Overall Assessment

Most of the major review findings appear valid.

The strongest remaining issues are:

- `ProteinScoringAndFdrEngine` still re-filters already-filtered inputs
- weak new tests that do not verify the intended behavior

The `magic booleans` comment is valid as a maintainability issue, but it is not a correctness bug by itself.

## Issue Triage

### 1. `PostSearchAnalysisTask` localization regression

- File: `MetaMorpheus/TaskLayer/SearchTask/PostSearchAnalysisTask.cs`
- Review assessment: Valid
- Severity: High
- Status: Finished

Base behavior enabled `DoLocalizationAnalysis` for `ModOpen`, `Open`, and `Custom` searches. That regression has been fixed and committed by restoring the previous automatic localization behavior while preserving the current handling for custom oligo searches.

Resolution:

- restore the previous unconditional enablement for `ModOpen` and `Open`
- preserve the intended `Custom` search exception for oligo workflows

### 2. `ProteinParsimonyEngine` filtered-overload can break revert intent

- File: `MetaMorpheus/EngineLayer/ProteinParsimony/ProteinParsimonyEngine.cs`
- Review assessment: Valid
- Severity: High
- Status: Finished

The new overload taking `FilteredPsms` forwarded `FilteredPsmsList` into the list-based constructor. If a caller used that overload, `_allPsms` became the already-filtered set rather than the full PSM set. That conflicted with the PR goal of using all PSMs for parsimony/disambiguation. The overload has now been removed.

Resolution:

- remove the `FilteredPsms` overload entirely and only accept the full unfiltered PSM list
- update production and test call sites to pass an explicit `List<SpectralMatch>`

Notes:

- `PostGlycoSearchAnalysisTask` now passes the full PSM list into `ProteinParsimonyEngine`
- existing tests were updated to pass `FilteredPsms.FilteredPsmsList` explicitly where they only need the filtered test data

### 3. `ProteinScoringAndFdrEngine` constructor re-derives filter semantics

- File: `MetaMorpheus/EngineLayer/ProteinScoringAndFdr/ProteinScoringAndFdrEngine.cs`
- Review assessment: Valid
- Severity: High
- Status: Finished

The older `FilteredPsms`-based path used `filteredPsms.FilterType` and `filteredPsms.FilterThreshold` directly. The new overload re-derived those values from `CommonParameters`, which could disagree with the actual filtered set that was passed in. That constructor path has now been removed.

Resolution:

- remove the `FilteredPsms` overload entirely
- update production and test call sites to pass `FilteredPsms.FilteredPsmsList` explicitly when they intend to score an already-filtered set

### 4. `ProteinScoringAndFdrEngine` re-filters already-filtered inputs

- File: `MetaMorpheus/EngineLayer/ProteinScoringAndFdr/ProteinScoringAndFdrEngine.cs`
- Review assessment: Valid
- Severity: High

Callers can now pass an already-filtered list explicitly via `FilteredPsms.FilteredPsmsList`, but `ScoreProteinGroups` still filters that list again. That is unnecessary work at minimum, and it can still become a correctness problem if the engine's `CommonParameters`-derived settings diverge from the list's original filtering semantics.

Potential solution:

- for the `FilteredPsms` overload, use the provided filtered list directly
- for the raw-list overload, filter once internally and carry that result through scoring

### 5. Tie behavior differs from `FilteredPsms`

- File: `MetaMorpheus/EngineLayer/ProteinScoringAndFdr/ProteinScoringAndFdrEngine.cs`
- Review assessment: Valid
- Severity: Medium-High

The current logic picks `PepQValue` when thresholds are equal because it uses `if (QValueThreshold < PepQValueThreshold) ... else PepQValue`. Existing `FilteredPsms` behavior keeps ties on `QValue`.

Potential solution:

- align the engine logic with `FilteredPsms`
- better yet, avoid re-deriving this decision where an existing `FilteredPsms` instance already captures it

### 6. Parsimony reporting test is too weak

- File: `MetaMorpheus/Test/StefanParsimonyTest.cs`
- Review assessment: Valid
- Severity: Medium

The new test only checks that `HypothesesAdded` and `HypothesesRemoved` are non-negative and that `ToString()` contains the field labels. That would still pass if the counts were always zero or swapped.

Potential solution:

- use a deterministic fixture and assert exact expected count values
- keep the `ToString()` assertions only as a secondary check

### 7. `QValue` filter-selection test is vacuous

- File: `MetaMorpheus/Test/StefanParsimonyTest.cs`
- Review assessment: Valid
- Severity: Medium

The assertion `results.SortedAndScoredProteinGroups.Count >= 0` is always true and does not verify branch behavior.

Potential solution:

- build a case where `QValue` and `PepQValue` produce different outcomes
- assert the exact output that proves the `QValue` branch was selected

### 8. `PepQValue` filter-selection test is vacuous

- File: `MetaMorpheus/Test/StefanParsimonyTest.cs`
- Review assessment: Valid
- Severity: Medium

This test has the same flaw as the `QValue` test. It does not validate the intended behavior.

Potential solution:

- use a crafted dataset where peptide-level filtering changes the result
- assert the exact observable outcome of selecting `PepQValue`

### 9. `FilterByQValue` uses positional booleans

- File: `MetaMorpheus/EngineLayer/ProteinScoringAndFdr/ProteinScoringAndFdrEngine.cs`
- Review assessment: Valid, but low priority
- Severity: Low

The current call is readable only if the reader knows the full method signature. The arguments appear to be correct, but the code is fragile and easy to misread.

Potential solution:

- switch the call to named arguments for clarity

## Recommended Fix Order

1. Remove double-filtering on already-filtered inputs.
2. Strengthen the new tests so they assert exact behavior rather than tautologies.
3. Clean up the positional-boolean call with named arguments.

## Notes

- GitHub CLI authentication was unavailable in this environment, so PR comments were gathered through public GitHub API endpoints.
- This assessment is based on the current branch contents and comparison against `upstream/master`.
