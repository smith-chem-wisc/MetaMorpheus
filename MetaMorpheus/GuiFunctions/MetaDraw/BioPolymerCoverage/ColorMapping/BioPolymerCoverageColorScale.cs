using System;

namespace GuiFunctions.MetaDraw;

public sealed class BioPolymerCoverageColorScale
{
    public BioPolymerCoverageColorScale(double minValue, double maxValue, ColorGradient gradient)
    {
        MinValue = minValue;
        MaxValue = maxValue;
        Gradient = gradient ?? throw new ArgumentNullException(nameof(gradient));
    }

    public double MinValue { get; }
    public double MaxValue { get; }
    public double Range => MaxValue - MinValue;
    public ColorGradient Gradient { get; }

    public double Normalize(double value)
    {
        if (Range <= 0)
            return 0;
        var normalized = (value - MinValue) / Range;
        if (normalized < 0) return 0;
        if (normalized > 1) return 1;
        return normalized;
    }
}
