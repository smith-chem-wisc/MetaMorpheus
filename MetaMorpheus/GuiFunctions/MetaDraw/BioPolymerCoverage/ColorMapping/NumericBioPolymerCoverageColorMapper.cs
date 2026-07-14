using System;
using System.Collections.Generic;
using System.Linq;
using System.Windows.Media;

namespace GuiFunctions.MetaDraw;

public abstract class NumericBioPolymerCoverageColorMapper : BioPolymerCoverageColorMapper
{
    public override bool IsNumeric => true;
    public override bool SupportsGradientSelection => true;
    public override bool SupportsLogScale => true;

    public abstract double? GetNumericValue(BioPolymerCoverageResultModel result);
    public abstract string DisplayName { get; }
    protected virtual NumericBioPolymerCoverageColorMapper? GetFallbackMapper() => null;

    public BioPolymerCoverageNumericRenderContext CreateRenderContext(
        IReadOnlyList<BioPolymerCoverageResultModel> filteredResults,
        ColorGradient gradient,
        bool useLogScale)
    {
        var rawValues = filteredResults
            .Select(GetNumericValue)
            .Where(v => v.HasValue)
            .Select(v => v.Value)
            .ToList();

        NumericBioPolymerCoverageColorMapper effectiveMapper = this;
        var usable = FilterUsableValues(rawValues, useLogScale);

        if (usable.Count == 0)
        {
            var fallback = GetFallbackMapper();
            if (fallback is not null)
            {
                rawValues = filteredResults
                    .Select(fallback.GetNumericValue)
                    .Where(v => v.HasValue)
                    .Select(v => v.Value)
                    .ToList();
                effectiveMapper = fallback;
                usable = FilterUsableValues(rawValues, useLogScale);
            }
        }

        var transformed = usable.Select(v => TransformValue(v, useLogScale)).ToList();
        BioPolymerCoverageColorScale scale;
        if (transformed.Count == 0)
        {
            scale = new BioPolymerCoverageColorScale(0, 0, gradient);
        }
        else
        {
            double minVal = transformed.Min();
            double maxVal = transformed.Max();
            if (maxVal - minVal < double.Epsilon) maxVal = minVal + 1;
            scale = new BioPolymerCoverageColorScale(minVal, maxVal, gradient);
        }

        string title = effectiveMapper.DisplayName;
        if (useLogScale) title += " (log10)";

        return new BioPolymerCoverageNumericRenderContext(effectiveMapper, scale, useLogScale, title);
    }

    public override SolidColorBrush GetBrush(
        BioPolymerCoverageResultModel result,
        BioPolymerCoverageColorScale? scale)
    {
        if (scale is null)
            return new SolidColorBrush(Colors.Gray);
        var value = GetNumericValue(result);
        if (!value.HasValue)
            return new SolidColorBrush(Colors.Gray);
        return scale.Gradient.GetBrush(scale.Normalize(value.Value));
    }

    public override string GetLegendTitle(BioPolymerCoverageColorScale? scale) => DisplayName;

    private static List<double> FilterUsableValues(IEnumerable<double> values, bool useLog)
    {
        if (!useLog) return values.ToList();
        return values.Where(v => v > 0).ToList();
    }

    private static double TransformValue(double value, bool useLog)
    {
        return useLog && value > 0 ? Math.Log10(value) : value;
    }

    // Used by callers when they need to apply the same log transform to a single result
    public static double TransformForRendering(double value, bool useLogScale)
    {
        return TransformValue(value, useLogScale);
    }

    // Validates that a value can be rendered under the given log mode
    public static bool IsRenderable(double? value, bool useLogScale)
    {
        return value.HasValue && (!useLogScale || value.Value > 0);
    }
}
