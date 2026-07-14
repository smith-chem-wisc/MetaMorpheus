using System;
using System.Collections.Generic;
using System.Linq;
using System.Windows.Media;

namespace GuiFunctions.MetaDraw;

public sealed class ViridisColorGradient : ColorGradient
{
    public override ColorGradientType GradientType => ColorGradientType.Viridis;
    public override int BinCount => 20;

    private static readonly IReadOnlyList<SolidColorBrush> Brushes = BuildBrushes();

    public override SolidColorBrush GetBrush(double normalizedValue)
    {
        int bin = (int)Math.Clamp(normalizedValue * (BinCount - 1), 0, BinCount - 1);
        return Brushes[bin];
    }

    public override IReadOnlyList<SolidColorBrush> GetBrushes() => Brushes;

    private static SolidColorBrush[] BuildBrushes()
    {
        var hex = new[]
        {
            "440154", "481467", "482576", "453781", "3f4d8a",
            "39568c", "2e6f8e", "238a8d", "20a386", "34b679",
            "60c85d", "93d443", "c7df2b", "eddf2b", "f8e51d",
            "fce61b", "fbdf0e", "f6d300", "edc600", "fde725"
        };
        return hex.Select(h => new SolidColorBrush((Color)ColorConverter.ConvertFromString("#" + h))).ToArray();
    }
}
