# BioPolymerCoverage Color Mapping - Implementation Plan

Source plan: refined execution plan with Strategy + Factory Method + Template
Method design. Implements a class-based and numeric-based color mapping
abstraction with `Viridis` as the first gradient.

## 1. Gradient types

- [x] Create `ColorMapping/ColorGradientType.cs` with `Viridis` enum value.
- [x] Create `ColorMapping/ColorGradient.cs` abstract base with:
  - `ColorGradientType GradientType { get; }`
  - `int BinCount { get; }`
  - `SolidColorBrush GetBrush(double normalizedValue)`
  - `IReadOnlyList<SolidColorBrush> GetBrushes()`
- [x] Create `ColorMapping/ViridisColorGradient.cs` with the current 20-bin palette.
- [x] Create `ColorMapping/ColorGradientFactory.cs` static `Create(ColorGradientType)`.

## 2. Scale value object

- [x] Create `ColorMapping/BioPolymerCoverageColorScale.cs` with:
  - `MinValue`, `MaxValue`, `Range`, `Gradient` properties
  - `Normalize(double value)` helper

## 3. Mapper hierarchy

- [x] Create `ColorMapping/BioPolymerCoverageColorMapper.cs` abstract base with:
  - `ColorResultsBy ColorBy { get; }`
  - `bool IsNumeric { get; }`
  - `SolidColorBrush GetBrush(BioPolymerCoverageResultModel, BioPolymerCoverageColorScale?)`
  - `virtual string GetLegendTitle(BioPolymerCoverageColorScale?)`
- [x] Create `ColorMapping/CategoricalBioPolymerCoverageColorMapper.cs`
  - constructor takes `ColorResultsBy` and `Func<BioPolymerCoverageResultModel, SolidColorBrush>`
- [x] Create `ColorMapping/NumericBioPolymerCoverageColorMapper.cs` abstract base
  - `abstract double? GetNumericValue(BioPolymerCoverageResultModel)`
  - base `GetBrush` does normalize + gradient lookup
- [x] Create `ColorMapping/PrecursorIntensityColorMapper.cs`
- [x] Create `ColorMapping/ScoreColorMapper.cs`

## 4. Mapper factory

- [x] Create `ColorMapping/BioPolymerCoverageColorMapperFactory.cs`
  - `Create(ColorResultsBy, Func<string, SolidColorBrush> identifierBrushResolver, ColorGradientType gradientType = Viridis)`
  - registrations:
    - `None` -> gray categorical mapper
    - `CoverageType` -> categorical mapper using `MetaDrawSettings.BioPolymerCoverageColors`
    - `FileOrigin` -> categorical mapper using identifier resolver
    - `PrecursorIntensity` -> `PrecursorIntensityColorMapper`
    - `Score` -> `ScoreColorMapper`

## 5. View model refactor (`BioPolymerCoverageMapViewModel.cs`)

- [x] Add `CreateColorMapper()` helper.
- [x] Add `CreateNumericScale(BioPolymerCoverageColorMapper, List<BioPolymerCoverageResultModel>)` helper.
- [x] Pass `GetColorByIdentifier` to factory as identifier resolver.
- [x] Replace intensity-specific normalization block in `Redraw()` with mapper + scale flow.
- [x] Update drawing loop to use `mapper.GetBrush(result, scale)`.
- [x] Update `CreateLegendItems` to keep returning empty for numeric modes (gradient bar draws via new path).
- [x] Extract gradient-bar drawing into `DrawNumericLegend(...)` returning consumed height.
- [x] Remove `ViridisColors`, `InitViridisPalette`, `GetIntensityBrush`, and the direct intensity ternary.

## 6. Tests

- [x] Add `ViridisColorGradient` tests (count, clamp, end colors) - 12 tests
- [x] Add `BioPolymerCoverageColorScale` tests (normalization, clamping, null guard) - 9 tests
- [x] Add numeric mapper tests (extract value, brush usage) - 19 tests
- [x] Add mapper factory tests (categorical vs numeric resolution) - 10 tests
- [x] Update view-model drawing tests to use new abstraction paths.
- [x] Add `CreateColorMapper` test for mapper resolution by `ColorBy`.
- [x] Add `CreateNumericScale` tests for non-numeric returns null and numeric derives min/max.
- [x] Add fallback test: `PrecursorIntensity` with all null falls back to score through new abstraction.

## 7. Build + verify

- [x] Build `GuiFunctions` project: 0 errors.
- [x] Build `Test` project: 0 errors.
- [x] Run all new and existing tests in `BioPolymerCoverageMapViewModelTests` and gradient/mapper tests.
- [x] Run full MetaDraw test suite to confirm no regressions (451 tests passing).

## 8. Commit + push

- [ ] Commit (only when explicitly approved).
- [ ] Push branch to origin.
- [ ] Update PR #2684 description with new architecture notes.
