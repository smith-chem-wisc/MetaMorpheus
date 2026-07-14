# BioPolymerCoverage Color Mapping - Implementation Plan

Source plan: refined execution plan with Strategy + Factory Method + Template
Method design. Implements a class-based and numeric-based color mapping
abstraction with `Viridis` as the first gradient.

## 1. Gradient types

- [ ] Create `ColorMapping/ColorGradientType.cs` with `Viridis` enum value.
- [ ] Create `ColorMapping/ColorGradient.cs` abstract base with:
  - `ColorGradientType GradientType { get; }`
  - `int BinCount { get; }`
  - `SolidColorBrush GetBrush(double normalizedValue)`
  - `IReadOnlyList<SolidColorBrush> GetBrushes()`
- [ ] Create `ColorMapping/ViridisColorGradient.cs` with the current 20-bin palette.
- [ ] Create `ColorMapping/ColorGradientFactory.cs` static `Create(ColorGradientType)`.

## 2. Scale value object

- [ ] Create `ColorMapping/BioPolymerCoverageColorScale.cs` with:
  - `MinValue`, `MaxValue`, `Range`, `Gradient` properties
  - `Normalize(double value)` helper

## 3. Mapper hierarchy

- [ ] Create `ColorMapping/BioPolymerCoverageColorMapper.cs` abstract base with:
  - `ColorResultsBy ColorBy { get; }`
  - `bool IsNumeric { get; }`
  - `SolidColorBrush GetBrush(BioPolymerCoverageResultModel, BioPolymerCoverageColorScale?)`
  - `virtual string GetLegendTitle(BioPolymerCoverageColorScale?)`
- [ ] Create `ColorMapping/CategoricalBioPolymerCoverageColorMapper.cs`
  - constructor takes `ColorResultsBy` and `Func<BioPolymerCoverageResultModel, SolidColorBrush>`
- [ ] Create `ColorMapping/NumericBioPolymerCoverageColorMapper.cs` abstract base
  - `abstract double? GetNumericValue(BioPolymerCoverageResultModel)`
  - base `GetBrush` does normalize + gradient lookup
- [ ] Create `ColorMapping/PrecursorIntensityColorMapper.cs`
- [ ] Create `ColorMapping/ScoreColorMapper.cs`

## 4. Mapper factory

- [ ] Create `ColorMapping/BioPolymerCoverageColorMapperFactory.cs`
  - `Create(ColorResultsBy, Func<string, SolidColorBrush> identifierBrushResolver, ColorGradientType gradientType = Viridis)`
  - registrations:
    - `None` -> gray categorical mapper
    - `CoverageType` -> categorical mapper using `MetaDrawSettings.BioPolymerCoverageColors`
    - `FileOrigin` -> categorical mapper using identifier resolver
    - `PrecursorIntensity` -> `PrecursorIntensityColorMapper`
    - `Score` -> `ScoreColorMapper`

## 5. View model refactor (`BioPolymerCoverageMapViewModel.cs`)

- [ ] Add `CreateColorMapper()` helper.
- [ ] Add `CreateNumericScale(BioPolymerCoverageColorMapper, List<BioPolymerCoverageResultModel>)` helper.
- [ ] Rename `GetColorByIdentifier` to `GetIdentifierBrush` and pass it to factory.
- [ ] Replace intensity-specific normalization block in `Redraw()` with mapper + scale flow.
- [ ] Update drawing loop to use `mapper.GetBrush(result, scale)`.
- [ ] Update `CreateLegendItems` to accept mapper + scale and return empty for numeric modes.
- [ ] Extract gradient-bar drawing into `DrawNumericLegend(...)` returning consumed height.
- [ ] Remove `ViridisColors`, `InitViridisPalette`, `GetIntensityBrush`, and the direct intensity ternary.

## 6. Tests

- [ ] Add `ViridisColorGradient` tests (count, clamp, end colors).
- [ ] Add mapper factory tests (categorical vs numeric resolution).
- [ ] Add numeric mapper tests (extract value, brush usage).
- [ ] Update view-model drawing tests to use new abstraction paths.
- [ ] Add `CreateNumericScale` test for min/max derivation.
- [ ] Add fallback test: `PrecursorIntensity` with all null falls back to score through new abstraction.

## 7. Build + verify

- [ ] Build `GuiFunctions` project: 0 errors.
- [ ] Build `Test` project: 0 errors.
- [ ] Run all new and existing tests in `BioPolymerCoverageMapViewModelTests` and gradient/mapper tests.
- [ ] Run full test suite to confirm no regressions.

## 8. Commit + push

- [ ] Commit (only when explicitly approved).
- [ ] Push branch to origin.
- [ ] Update PR #2684 description with new architecture notes.
