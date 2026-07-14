# BioPolymerCoverage Color Mapping

## Intent

Provide an extensible, object-oriented way to color peptide rectangles in the
BioPolymer coverage view. The current implementation supports two families of
coloration:

- **Class-based**: a discrete brush is selected per result based on a category
  (e.g. `CoverageType`, `FileOrigin`).
- **Numeric-based**: a numeric value is extracted from the result and mapped to
  a gradient using a provided `[min, max]` scale (e.g. `PrecursorIntensity`,
  `Score`).

The numeric path does **not** own min/max discovery, log behavior, or filtering.
Those concerns stay in the view model. Mappers are also independent of any
specific gradient implementation, so new gradients can be added without
touching mapper code.

## Structure

```
ColorMapping/
  ColorResultsBy.cs              # user-facing mode selector
  ColorGradientType.cs           # enum of implemented gradients
  ColorGradient.cs               # abstract gradient base
  ViridisColorGradient.cs        # first concrete gradient
  ColorGradientFactory.cs        # resolves a gradient by enum
  BioPolymerCoverageColorScale.cs # value object holding min/max + gradient
  BioPolymerCoverageColorMapper.cs # abstract mapper base
  CategoricalBioPolymerCoverageColorMapper.cs # categorical strategy
  NumericBioPolymerCoverageColorMapper.cs      # numeric base (Template Method)
  PrecursorIntensityColorMapper.cs
  ScoreColorMapper.cs
  BioPolymerCoverageColorMapperFactory.cs      # resolves a mapper by enum
```

### Roles

- **`ColorResultsBy`** — user-facing mode (`None`, `CoverageType`, `FileOrigin`,
  `PrecursorIntensity`, `Score`, ...).
- **`ColorGradientType`** — enum of implemented gradients (`Viridis`, ...).
- **`ColorGradient`** — abstract gradient with `BinCount`, `GetBrush(double
  normalizedValue)`, and `GetBrushes()` for legend rendering.
- **`ViridisColorGradient`** — current 20-bin viridis palette.
- **`ColorGradientFactory`** — resolves a gradient by `ColorGradientType`.
- **`BioPolymerCoverageColorScale`** — value object holding
  `MinValue` / `MaxValue` / `Range` / `Gradient` plus a `Normalize(double)`
  helper. Numeric mappers consume this; they do not compute it.
- **`BioPolymerCoverageColorMapper`** — abstract base. Exposes
  `ColorBy`, `IsNumeric`, and `GetBrush(result, scale?)`.
- **`CategoricalBioPolymerCoverageColorMapper`** — strategy that delegates
  brush selection to a `Func<result, brush>`.
- **`NumericBioPolymerCoverageColorMapper`** — abstract numeric base. Subclasses
  provide `GetNumericValue(result)`; the base handles normalize + gradient
  lookup.
- **`PrecursorIntensityColorMapper` / `ScoreColorMapper`** — concrete numeric
  strategies.
- **`BioPolymerCoverageColorMapperFactory`** — resolves a mapper for a
  `ColorResultsBy`, injecting required dependencies (e.g. identifier brush
  resolver, gradient type).

### How the view model uses it

The view model is responsible for:

- filtering
- min/max derivation (only for numeric mappers)
- legend layout
- identifier brush caching (currently `FileOrigin`)

The view model does **not** know how a specific mode produces a brush. It asks
the factory for a mapper and asks the mapper for a brush.

Typical flow inside `Redraw`:

```csharp
var mapper = BioPolymerCoverageColorMapperFactory.Create(
    ColorBy,
    identifierBrushResolver: GetIdentifierBrush,
    gradientType: ColorGradientType.Viridis);

var scale = BuildNumericScale(mapper, filteredResults);

foreach (var result in filteredResults)
{
    var brush = mapper.GetBrush(result, scale);
    // draw using brush
}
```

## How to extend

### Add a new categorical mode (e.g. `DatabaseOrigin`)

1. Add the enum value to `ColorResultsBy`.
2. Register the mapping in `BioPolymerCoverageColorMapperFactory.Create` using
   `CategoricalBioPolymerCoverageColorMapper` and a brush selector delegate.
3. No changes to mappers, gradients, or numeric logic.

### Add a new numeric mode (e.g. `QValue`)

1. Add the enum value to `ColorResultsBy`.
2. Create a new `QValueColorMapper : NumericBioPolymerCoverageColorMapper`
   that implements `GetNumericValue(result)`.
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
