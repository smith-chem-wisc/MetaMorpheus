# Truncation search

`TruncationSearchTask` is a follow-on MetaMorpheus task that discovers **N- and C-terminally
truncated proteoforms** of the parents identified by an upstream top-down `SearchTask`. It is a
normal task: it has its own TOML schema, it runs from the GUI or from CMD in a run list such as
`[Search, Truncation]`, and it writes results using the standard `.psmtsv` columns.

This document is the design contract the `EngineLayer/Truncation` and
`TaskLayer/TruncationSearchTask` code refers to. Numbered decisions (`#1`-`#21`) are cited from
doc-comments throughout that code; the numbering is stable, and gaps in it (there is no `#11`) are
historical.

---

## Overview

A truncated proteoform is *not* enumerated in advance. It is **deduced per scan**:

1. build a fragment index from the full-length parents an upstream search already identified;
2. work out, from the asymmetry of the matched fragment ions, which terminus a scan retains;
3. chop residues off the *other* terminus until the deduced mass matches the observed precursor.

That keeps the search space the size of the parent set rather than the size of the parent set times
every possible cut point.

## Three passes

- **Pass 1** - a standard top-down `SearchTask` runs to completion, including its full
  `PostSearchAnalysisTask`. It produces the proteoform-level `SpectralMatch` set (the same set that
  is written to `AllProteoforms.psmtsv`).
- **Pass 2** - `TruncationSearchTask` builds a fragment index over the parent proteoforms, then for
  each MS2 scan scores every mass-eligible parent twice: once against the N-terminal ion series
  only, once against the C-terminal ion series only. The winning `(parent, series)` pair names the
  parent and the retained terminus.
- **Pass 3** - residues (with their PTMs) are chopped from the terminus opposite the winning series
  until the precursor mass matches at an allowed notch; the deduced truncated form is then scored
  against the scan with standard both-series scoring and emitted as a `SpectralMatch`.

Pass 1 intact matches are inherited as `full-length`, pooled with the Pass 3 truncations, and the
whole pool gets one fresh per-class FDR/PEP analysis.

---

## The locked decisions

### #1 - Pass 1 source and task structure

- `TruncationSearchTask : MetaMorpheusTask`, parallel to `SearchTask` / `CalibrationTask` /
  `GptmdTask`.
- Pass 1 results reach Pass 2 **in memory**, not through a file round-trip.
  `EverythingRunnerEngine` owns a `TaskChainContext` - a typed deposit/retrieve slot keyed by task
  id - which a finishing `SearchTask` deposits into and `TruncationSearchTask` reads from. The
  runner hands out the context only when the run list actually contains a consumer, and clears it
  once that consumer has run, so nothing is pinned for chains that never read it.
- Fallback for standalone runs: parse `AllProteoforms.psmtsv` from a previous output folder
  (`Pass1ProteoformsFilePath`).
- The payload deposited on the context is the search's PSM-level `SpectralMatch` set as-searched.
  The consumer dedups it to proteoform level itself (best PSM per `FullSequence`).

### #2 - Pipe-ambiguous proteoforms

- A Pass 1 row whose `FullSequence` carries pipe-separated alternatives
  (`PEPTIDEK|PEPTIDER|PEPTIDEM`) contributes each alternative as its own parent.
- All alternatives are tagged back to the originating Pass 1 row so downstream counting does not
  double-count.

### #3 - Parent inclusion filter

- Targets and decoys live in one combined list, filtered identically.
- Filter: **PEP q-value <= `ParentQValueThreshold`** when PEP was computed, otherwise **notch
  q-value <= `ParentQValueThreshold`**.
- The default threshold (0.10) is deliberately permissive: cast a wide net of parents and let the
  truncation search's own FDR police the final output.

### #4 - Pass 2 index content

- The index holds fragments of the **parent** proteoforms only - b/y, c/z, etc. per dissociation
  type. No pre-expansion of truncated subsequences.

### #4a - Which MS2 scans enter Pass 2

- Every MS2 scan, including ones Pass 1 already identified.
- A scan is skipped only when Pass 1 has a PSM for it **and** that parent's theoretical mass equals
  the scan's precursor mass within `CommonParameters.PrecursorMassTolerance`. Those scans inherit
  their Pass 1 PSM into the Pass 3 output with `Description = full-length`.

### #5 - Fragment ion types

- The index follows `CommonParameters.DissociationType`, exactly as the existing `IndexingEngine`
  does. Terminus deduction generalizes to the **N-terminal series (a, b, c)** versus the
  **C-terminal series (x, y, z)** - it is not hard-coded to b/y.

### #6 - Maximum parent mass

- Parents heavier than `MaxParentMass` (default 30 kDa, the existing
  `SearchParameters.DefaultMaxFragmentSize`) are excluded from the index, reported as a single
  summary line: `"{N} parent proteoforms exceeded MaxFragmentSize and were excluded."`
- Raise this for database-seeded runs, where a large protein's intact form may never be observed
  but its truncations are.

### #7 - Pass 2 precursor acceptor

- `TruncationAcceptor`, a `MassDiffAcceptor` whose rule is
  `0 < M_obs <= M_theo + precursor tolerance` - the parent must be at least as heavy as the scan,
  and may be exceeded only by the tolerance (which absorbs measurement error and the
  intact-equality case).
- Single category, single notch. Notches enter at Pass 3, not here.

### #8 - Pass 2 scoring and terminus deduction

- Each scan x parent gets **two** scores: N-series-only and C-series-only.
- The highest-scoring `(parent, series)` pair across all parents wins the scan.
- Opposing-terminus ions are ignored in the chosen direction - no purity threshold, no rejection.

### #9 - Chopping algorithm

- One residue at a time from the indicated terminus, subtracting the residue's monoisotopic mass
  **plus any PTM locked to that residue** in one step; an N-terminal acetyl-Ala leaves as a single
  step.
- Notch behaviour is inherited from the task's `MassDiffAcceptorType` (see #15) and applies to
  `M_obs - M_truncated_theo`. Tolerance is `CommonParameters.PrecursorMassTolerance`.
- Stop at the first match, or when no residues remain. Because chopping is monotonic in mass, the
  loop also stops as soon as the truncated form falls below the lowest mass the acceptor could
  accept for that scan - the acceptor reports that bound itself, via
  `GetAllowedPrecursorMassIntervalsFromObservedMass`.
- **No mass-shift absorption.** If no integer-residue chop matches at any allowed notch, the scan
  is dropped - no Pass 3 PSM, no open-search-style residual annotation.
- Only normal backbone amide bonds are cut.

### #10 - PTM bookkeeping and parent-ambiguity collapse

- PTMs lock to residues. Chopping a residue removes its mod; a terminal mod survives only while its
  terminus is intact.
- Parents that chop down to the same truncated `FullSequence` for the same scan **collapse into one
  PSM**; protein-accession ambiguity is retained as pipe-separated accessions, per MetaMorpheus
  convention.
- A clean removal of exactly the initiator methionine is labelled *N-terminal Met excision*, not a
  1-residue truncation, so NME is not miscounted.

### #12 - Pass 3 scoring

- Standard MetaMorpheus scoring: fragment the deduced truncated `FullSequence` with the task's
  dissociation type, match against **both** series (the single-series restriction was only for
  parent selection), score through `CalculatePeptideScore`.
- The Pass 3 score is what feeds FDR. Pass 2 scores are internal.

### #13 - Output schema

- The truncation type rides the existing **`Description`** column - no new column - mirroring the
  proteolysis-product descriptors already carried there (`chain(2-121)`). Values are
  `N-terminal truncation(<start>-<end>)`, `C-terminal truncation(<start>-<end>)`,
  `N-terminal Met excision(...)`, or `full-length`.
- **Protein accession** is the parent's, unchanged. **Start/end residue** are the truncated
  coordinates within the parent protein, which makes downstream cleavage-motif analysis trivial.

### #14 - Decoys

- Decoys come from two sources, pooled into the same index: Pass 1 decoy proteoforms that passed
  the #3 filter, plus fresh reverse decoys generated for every target parent that has no
  reverse-decoy counterpart.
- Generation reuses the existing `DecoyType.Reverse` infrastructure; PTMs lock to residues during
  reversal.
- Pass 2 scoring is target/decoy-blind - scans compete naturally.
- Diagnostic side files (`CandidateRanks.tsv`) and the perf log are best-effort: a failure to write
  one never aborts an otherwise successful task.

### #15 - Notches and FDR categories

- Notches are the standard monoisotopic offsets on `M_obs - M_truncated_theo` (not on parent mass),
  inherited from the task's `MassDiffAcceptorType`. The default `ThreeMM` keeps a precursor isotope
  error from blocking a chop. Truncation length does not affect notch assignment.
- N-terminal and C-terminal truncations are **pooled** into one FDR analysis, for statistical power.

### #16 / #17 - Output files and filters

- `AllTruncatedPSMs.psmtsv` - one row per scan x candidate match.
- `AllTruncatedProteoforms.psmtsv` - deduped by truncated `FullSequence`.
- **No q-value or PEP cutoff at write time.** Both files carry the full FDR columns (q-value, notch
  q-value, PEP, PEP q-value) and the user filters downstream, exactly as with `AllPSMs.psmtsv`.
- `WriteDecoys` and `WriteContaminants` default to `true`; either can be turned off in the TOML.

### #18 - Pooled FDR universe

- The table that feeds FDR + PEP is the union of the Pass 3 truncation PSMs and the Pass 1 intact
  matches inherited per #4a.
- All q-values and PEPs in the two output files are computed **fresh from that pooled
  distribution**. Pass 1's own `AllPSMs.psmtsv` keeps its independently computed values; the two
  output families are separate FDR universes by design.

### #19 - Out of scope for v1

- **No internal truncations** (both termini cut in the same proteoform). v1 is strictly terminal.
- No mass-shift absorption for unmatched chops. No charge-state analysis beyond existing code.

### #20 - Task wiring

- Own TOML schema; most parameters follow Pass 1 conventions (dissociation type, precursor and
  product tolerances, notch type, decoy type).
- Plumbed into the existing runner, so `[Search, Truncation]` works end-to-end with in-memory
  result passing, and into `CMD/Program.cs` so the task dispatches from a TOML on the command line.

### #21 - Test and benchmark data

- Primary benchmark: **PXD049969** (Kaulich et al. 2024, *Nature Methods*).
- Secondary validation: an in-house 20-file Jurkat top-down set.
- Unit-test fixtures are snippets carved from the middle of selected raw files with the existing
  mzLib snip utility.

---

## Optional: database-seeded parents

`SeedParentsFromDatabase` generates parents by digesting the database — theoretical full-length
proteoforms built the way a top-down search builds them, so annotated PTMs ride along — instead of
taking only the proteoforms an upstream search identified. That removes the requirement that a
truncation's parent be observed intact, at the cost of a much larger Pass 2 search; raise
`MaxParentMass` alongside it, since a large protein's intact form may never be observed even when
its truncations are.

---

## Benchmarking hook

Performance instrumentation is deliberately **not** part of the task's TOML schema. Set the
environment variable

```
MetaMorpheusTruncationPerfLog=<path to perf_log.tsv>
```

and the task appends one metrics-and-timings row per run to that rolling, append-only TSV (header
written once, never overwritten), and additionally writes a `CandidateRanks.tsv` side file
recording each scan's winning-parent rank in the index match-count ordering. Leave the variable
unset - the default - and nothing is logged; ordinary runs and tests are unaffected.

A row carries the run's PSM and proteoform counts at 1% and 5% FDR, truncation counts by class,
parent and decoy counts, and per-stage wall-clock (index build, Pass 2 scoring, Pass 3 chopping,
Pass 3 scoring, FDR/PEP). Run metadata is parsed from the parent run-folder name when it follows
the `<date>_<phase>_<datasetTag>_<runLabel>` convention.
