using GuiFunctions.MetaDraw.BioPolymerCoverage.ColorMapping.Gradient;
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

    private BioPolymerCoverageColorScale? _preparedScale;
    private bool _preparedUseLogScale;
    private string? _preparedLegendTitle;

    public override void Prepare(
        IReadOnlyList<BioPolymerCoverageResultModel> filteredResults,
        ColorGradientType gradientType,
        bool useLogScale)
    {
        var gradient = ColorGradientFactory.Create(gradientType);
        var values = filteredResults
            .Select(GetNumericValue)
            .Where(v => v.HasValue)
            .Select(v => v.Value)
            .ToList();

        var usable = FilterUsableValues(values, useLogScale);
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

        string title = DisplayName;
        if (useLogScale) title += " (log10)";

        _preparedScale = scale;
        _preparedUseLogScale = useLogScale;
        _preparedLegendTitle = title;
    }

    public override SolidColorBrush GetBrush(BioPolymerCoverageResultModel result)
    {
        if (_preparedScale is null)
            return new SolidColorBrush(Colors.Gray);
        var value = GetNumericValue(result);
        if (!IsRenderable(value, _preparedUseLogScale))
            return new SolidColorBrush(Colors.Gray);
        var transformed = TransformForRendering(value!.Value, _preparedUseLogScale);
        return _preparedScale.Gradient.GetBrush(_preparedScale.Normalize(transformed));
    }

    public override string? GradientLegendTitle => _preparedLegendTitle;
    public override IReadOnlyList<SolidColorBrush>? GradientBrushes => _preparedScale?.Gradient.GetBrushes();
    public override double? GradientMinValue => _preparedScale?.MinValue;
    public override double? GradientMaxValue => _preparedScale?.MaxValue;

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
