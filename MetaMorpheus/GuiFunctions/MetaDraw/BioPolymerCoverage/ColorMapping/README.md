# BioPolymerCoverage Color Mapping

## Intent

Provide an extensible, object-oriented way to color peptide rectangles in the
BioPolymer coverage view. The current implementation supports two families of
coloration:

- **Class-based**: a discrete brush is selected per result based on a category
  (e.g. `CoverageType`, `FileOrigin`).
- **Numeric-based**: a numeric value is extracted from the result and mapped to
  a gradient using a prepared `[min, max]` scale (e.g. `PrecursorIntensity`,
  `Score`).

The numeric path does **not** own min/max discovery, log behavior, or filtering.
Those concerns stay in the view model. Mappers are also independent of any
specific gradient implementation, so new gradients can be added without
touching mapper code.

## Structure

```
ColorMapping/
  Gradient/
    ColorGradient.cs               # abstract gradient base
    ColorGradientType.cs           # enum of implemented gradients
    ColorGradientFactory.cs        # resolves a gradient by enum
    ViridisColorGradient.cs        # 20-bin viridis palette
    PlasmaColorGradient.cs         # 20-bin plasma palette
    InfernoColorGradient.cs        # 20-bin inferno palette
    TurboColorGradient.cs          # 20-bin turbo palette
    GrayscaleColorGradient.cs      # 20-bin grayscale palette
  BioPolymerCoverageColorScale.cs  # value object holding min/max + gradient
  BioPolymerCoverageColorMapper.cs # abstract mapper base
  CategoricalBioPolymerCoverageColorMapper.cs # categorical strategy
  NumericBioPolymerCoverageColorMapper.cs     # numeric base (Template Method)
  FileOriginColorMapper.cs         # stateful categorical mapper for FileOrigin
  PrecursorIntensityColorMapper.cs
  ScoreColorMapper.cs
  BioPolymerCoverageColorMapperFactory.cs     # resolves a mapper by enum
```

### Roles

- **`ColorResultsBy`** — user-facing mode (`None`, `CoverageType`, `FileOrigin`,
  `PrecursorIntensity`, `Score`, ...).
- **`ColorGradientType`** — enum of implemented gradients (`Viridis`,
  `Plasma`, `Inferno`, `Turbo`, `Grayscale`, ...).
- **`ColorGradient`** — abstract gradient with `BinCount`, `GetBrush(double
  normalizedValue)`, and `GetBrushes()` for legend rendering.
- **`ColorGradientFactory`** — resolves a gradient by `ColorGradientType`.
- **`BioPolymerCoverageColorScale`** — value object holding
  `MinValue` / `MaxValue` / `Range` / `Gradient` plus a `Normalize(double)`
  helper. Numeric mappers consume this; they do not compute it.
- **`BioPolymerCoverageColorMapper`** — abstract base. Exposes
  `ColorBy`, `IsNumeric`, `SupportsGradientSelection`, `SupportsLogScale`,
  `DefaultUseLogScale`, `Prepare(...)`, `GetBrush(result)`, and legend
  metadata properties (`LegendItems`, `GradientLegendTitle`, `GradientBrushes`,
  `GradientMinValue`, `GradientMaxValue`).
- **`CategoricalBioPolymerCoverageColorMapper`** — strategy that delegates
  brush selection to a `Func<result, brush>`. Used for `None` and `CoverageType`.
- **`FileOriginColorMapper`** — stateful categorical mapper that owns a
  color queue and dictionary to assign stable per-file colors. Used for
  `FileOrigin` mode.
- **`NumericBioPolymerCoverageColorMapper`** — abstract numeric base. Subclasses
  provide `GetNumericValue(result)`, `DisplayName`, and optionally
  `DefaultUseLogScale`. The base handles normalize, gradient lookup, and log
  transform.
- **`PrecursorIntensityColorMapper` — concrete numeric strategy with
  `DefaultUseLogScale = true`.
- **`ScoreColorMapper`** — concrete numeric strategy with
  `DefaultUseLogScale = false`.
- **`BioPolymerCoverageColorMapperFactory`** — resolves a mapper for a
  `ColorResultsBy`. No parameters needed — all dependency injection is
  internal to the mappers.

### How the view model uses it

The view model creates a mapper from the factory, prepares it with the
current filtered results and settings, then caches it. Drawing queries the
cached mapper for each result.

Typical flow inside `Redraw`:

```csharp
var mapper = BioPolymerCoverageColorMapperFactory.Create(_colorBy);
mapper.Prepare(filteredResults, selectedGradientType, useLogColorScale);

foreach (var result in filteredResults)
{
    var brush = mapper.GetBrush(result);
    // draw using brush
}
```

Legend data is also provided by the prepared mapper:
```csharp
var legendTitle = mapper.GradientLegendTitle;
var gradientBrushes = mapper.GradientBrushes;
var minValue = mapper.GradientMinValue;
var maxValue = mapper.GradientMaxValue;
```

## How to extend

### Add a new categorical mode (e.g. `DatabaseOrigin`)

1. Add the enum value to `ColorResultsBy`.
2. Register the mapping in `BioPolymerCoverageColorMapperFactory.Create` using
   `CategoricalBioPolymerCoverageColorMapper` and a brush selector delegate.
3. Optionally provide a legend builder delegate.
4. No changes to mappers, gradients, or numeric logic.

### Add a new numeric mode (e.g. `QValue`)

1. Add the enum value to `ColorResultsBy`.
2. Create a new `QValueColorMapper : NumericBioPolymerCoverageColorMapper`
   that implements `GetNumericValue(result)`, `DisplayName`, and optionally
   `DefaultUseLogScale` and `GetFallbackMapper()`.
3. Register it in the factory.
4. No changes to gradient code, normalization pipeline, or categorical
   mappers.

### Add a new gradient (e.g. `Plasma`)

1. Add the enum value to `ColorGradientType`.
2. Create `PlasmaColorGradient : ColorGradient`.
3. Update `ColorGradientFactory.Create` to return it for the new enum value.
4. No changes to mappers, mapper factory, or view model.

## Design notes

- **Strategy**: each coloration mode is a pluggable mapper.
- **Factory Method / Registry**: mappers are constructed by
  `BioPolymerCoverageColorMapperFactory` from a `ColorResultsBy` selection.
- **Template Method**: numeric mappers share a normalize-and-lookup flow and
  vary only by how they extract a numeric value.
- **Bridge-style separation**: gradient implementations are independent from
  numeric extraction.
