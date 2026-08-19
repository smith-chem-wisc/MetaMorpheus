## Summary

- Glyco search's post-search PEP/FDR analysis always used the Chronologer retention-time (RT) predictor, a heavier deep-learning model, which noticeably increased search time. This PR makes the RT predictor configurable instead of hardcoded.
- Adds a "Retention Time Predictor" option (Chronologer / SSRCalc / Prosit2019iRT / Prosit2020iRTTMT) to the Post-Search Analysis section of both `SearchTaskWindow` and `GlycoSearchTaskWindow`, backed by the existing `CommonParameters.RTPredictorName` field and `FdrAnalysisEngine.GetRTPredictor` switch (Prosit options already existed engine-side; SSRCalc is newly added).
- Default stays **Chronologer** for both task types — fully backward compatible, purely opt-in.
- The picker is hidden in **RNA mode** on `SearchTaskWindow`, since `FdrAnalysisEngine.GetRTPredictor` short-circuits to `null` before ever reading the name whenever `searchType != "standard"` — the choice is genuinely inert there, unlike CZE (see below). `GlycoSearchTaskWindow` has no RNA mode, so no change was needed there.
- The picker is **not** hidden for CZE separation type. An earlier version of this PR hid it there too, on the assumption CZE never uses `IRetentionTimePredictor`. That assumption was wrong: `PEPAnalysisEngine.BuildFileSpecificDictionaries` calls `ComputeRetentionTimeEquivalentValues(..., RetentionTimePredictor)` unconditionally when building the training dictionaries, so the selected predictor is invoked for CZE runs too — only the later per-PSM z-score step branches to `GetMobilityZScore` instead of using that result. Hiding the picker for CZE would have removed the user's only way to avoid paying for the expensive Chronologer path on CZE data, while doing nothing to reduce that cost. The hide was reverted.
- Fixed a separate, pre-existing bug this PR's new setting made observable: `MetaMorpheusTask.SetAllFileSpecificCommonParams` rebuilds `CommonParameters` from scratch whenever a file-specific `<basename>.toml` exists (e.g. the `<file>-calib.toml` that `CalibrationTask` writes automatically), and did not forward `RTPredictorName` into the rebuilt object — so any pipeline that includes a Calibration step before Search/GlycoSearch silently reverted the user's predictor choice back to Chronologer. Now fixed by forwarding `rtPredictorName` through that rebuild.
- The four predictor-name strings (`"Chronologer"`, `"SSRCalc"`, `"Prosit2019iRT"`, `"Prosit2020iRTTMT"`) were previously duplicated as literals in six places across two assemblies (both GUI windows and the engine switch). Consolidated into a single `EngineLayer.FdrAnalysis.RTPredictorNames` constants class referenced everywhere, so a typo is now a compile error instead of a silent fallback to no predictor.

## Changes

- `EngineLayer/FdrAnalysis/RTPredictorNames.cs` (new): shared `const string` names for the four predictors.
- `EngineLayer/CommonParameters.cs`: default `rtPredictorName` now references `RTPredictorNames.Chronologer`.
- `EngineLayer/FdrAnalysis/FdrAnalysisEngine.cs`: `GetRTPredictor`'s switch gains `RTPredictorNames.SSRCalc => new SSRCalc3RetentionTimePredictor()`; all four arms use the shared constants.
- `TaskLayer/MetaMorpheusTask.cs`: `SetAllFileSpecificCommonParams` now forwards `rtPredictorName: commonParams.RTPredictorName` into the rebuilt `CommonParameters` — fixes the file-specific-toml bug described above.
- `GUI/TaskWindows/SearchTaskWindow.xaml(.cs)`: new "Retention Time Predictor" `GroupBox` (4 radio buttons); save/restore wiring to `CommonParameters.RTPredictorName` via the shared constants; `Visibility` bound to the existing `CollapseOnRnaModeConverter` (already used 7 other places in this file) so the group collapses in RNA mode only.
- `GUI/TaskWindows/GlycoSearchTaskWindow.xaml(.cs)`: same 4-option predictor `GroupBox` and save/restore wiring — no visibility toggle, not applicable to this window.
- `Test/FdrTest.cs`:
  - `FdrAnalysisEngine_GetRTPredictor_ReturnsExpectedPrositPredictor` — covers `Prosit2019iRT`/`Prosit2020iRTTMT` (network-dependent, sits behind a currently-inactive `Explicit`/`Category("Koina")` gate).
  - `FdrAnalysisEngine_GetRTPredictor_ReturnsExpectedLocalPredictor` (new) — covers `Chronologer`/`SSRCalc`, deliberately kept in a separate, ungated method so this coverage can't be silently excluded if the Koina gate above is ever turned on.
  - `SearchTask_DefaultRTPredictor_ResolvedToChronologer` — proves both `SearchTask` and `GlycoSearchTask` default to Chronologer; uses an explicit `case`/failing `default` so a typo'd `[TestCase]` label fails loudly instead of silently testing the wrong task type twice.
  - `SetAllFileSpecificCommonParams_PreservesRTPredictorName` (new) — regression test for the `MetaMorpheusTask.cs` fix above.

## Test plan

- [x] `dotnet build` on the full solution — succeeds with 0 errors.
- [x] `Test/FdrTest.cs`, full class (24 tests) — all passing. Coverage of the four predictor names is split across two tests by network dependency: `Chronologer`/`SSRCalc` in `FdrAnalysisEngine_GetRTPredictor_ReturnsExpectedLocalPredictor`, `Prosit2019iRT`/`Prosit2020iRTTMT` in `FdrAnalysisEngine_GetRTPredictor_ReturnsExpectedPrositPredictor`; Chronologer's role as the *default* (as opposed to just a valid selection) is separately covered by `SearchTask_DefaultRTPredictor_ResolvedToChronologer`.
- [x] Manual GUI verification (not covered by automated tests — `Test.csproj` doesn't reference `GUI.csproj`): opened both task windows, confirmed each of the 4 radio options saves and restores correctly on reopen for both `SearchTask` and `GlycoSearchTask`; confirmed the predictor group box disappears in `SearchTaskWindow` when RNA mode is active and stays visible for CZE.
- [ ] No behavior change for existing saved `.toml` task files without an explicit `RTPredictorName` — they continue to default to Chronologer.
