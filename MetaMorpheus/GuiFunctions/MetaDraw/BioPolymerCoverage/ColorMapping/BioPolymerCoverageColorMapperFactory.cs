using System;
using System.Windows.Media;

namespace GuiFunctions.MetaDraw;

public static class BioPolymerCoverageColorMapperFactory
{
    public static BioPolymerCoverageColorMapper Create(
        ColorResultsBy colorBy,
        Func<string, SolidColorBrush>? identifierBrushResolver = null,
        ColorGradientType gradientType = ColorGradientType.Viridis)
    {
        switch (colorBy)
        {
            case ColorResultsBy.None:
                return new CategoricalBioPolymerCoverageColorMapper(
                    ColorResultsBy.None,
                    _ => new SolidColorBrush(Colors.Gray));

            case ColorResultsBy.CoverageType:
                return new CategoricalBioPolymerCoverageColorMapper(
                    ColorResultsBy.CoverageType,
                    r => MetaDrawSettings.BioPolymerCoverageColors[r.CoverageType]);

            case ColorResultsBy.FileOrigin:
                return new CategoricalBioPolymerCoverageColorMapper(
                    ColorResultsBy.FileOrigin,
                    r => identifierBrushResolver?.Invoke(r.Match.FileName) ?? new SolidColorBrush(Colors.Gray));

            case ColorResultsBy.PrecursorIntensity:
                return new PrecursorIntensityColorMapper();

            case ColorResultsBy.Score:
                return new ScoreColorMapper();

            default:
                throw new ArgumentOutOfRangeException(nameof(colorBy), colorBy, "Unsupported ColorResultsBy");
        }
    }
}
