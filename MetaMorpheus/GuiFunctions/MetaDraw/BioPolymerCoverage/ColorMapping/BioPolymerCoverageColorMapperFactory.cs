using System;
using System.Linq;
using System.Windows.Media;

namespace GuiFunctions.MetaDraw;

public static class BioPolymerCoverageColorMapperFactory
{
    private static readonly BioPolymerCoverageType[] LegendOrder =
    {
        BioPolymerCoverageType.Unique,
        BioPolymerCoverageType.UniqueMissedCleavage,
        BioPolymerCoverageType.TandemRepeat,
        BioPolymerCoverageType.TandemRepeatMissedCleavage,
        BioPolymerCoverageType.Shared,
        BioPolymerCoverageType.SharedMissedCleavage
    };

    private static string AddSpaces(string value)
    {
        if (string.IsNullOrEmpty(value))
            return value;

        var chars = new System.Text.StringBuilder();
        chars.Append(value[0]);
        for (int i = 1; i < value.Length; i++)
        {
            if (char.IsUpper(value[i]) && !char.IsWhiteSpace(value[i - 1]))
                chars.Append(' ');
            chars.Append(value[i]);
        }
        return chars.ToString();
    }

    public static BioPolymerCoverageColorMapper Create(
        ColorResultsBy colorBy)
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
                    r => DrawnSequence.ParseColorBrushFromOxyColor(MetaDrawSettings.BioPolymerCoverageColors[r.CoverageType]),
                    filteredResults => LegendOrder
                        .Where(t => filteredResults.Any(r => r.CoverageType == t))
                        .Select(t =>
                        {
                            var count = filteredResults.Count(r => r.CoverageType == t);
                            var unitLabel = filteredResults.FirstOrDefault()?.Match.GetDigestionProductLabel();
                            return (DrawnSequence.ParseColorBrushFromOxyColor(MetaDrawSettings.BioPolymerCoverageColors[t]), $"{AddSpaces(t.ToString())} {unitLabel}s: {count}");
                        })
                        .ToList());

            case ColorResultsBy.FileOrigin:
                return new FileOriginColorMapper();

            case ColorResultsBy.PrecursorIntensity:
                return new PrecursorIntensityColorMapper();

            case ColorResultsBy.Score:
                return new ScoreColorMapper();

            default:
                throw new ArgumentOutOfRangeException(nameof(colorBy), colorBy, "Unsupported ColorResultsBy");
        }
    }
}
