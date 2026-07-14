using System.Collections.Generic;
using System.Linq;
using System.Windows.Media;

namespace GuiFunctions.MetaDraw;

public sealed class TurboColorGradient : ColorGradient
{
    public override ColorGradientType GradientType => ColorGradientType.Turbo;

    private static readonly IReadOnlyList<SolidColorBrush> Brushes = BuildBrushes();

    public override IReadOnlyList<SolidColorBrush> GetBrushes() => Brushes;

    private static SolidColorBrush[] BuildBrushes()
    {
        var hex = new[]
        {
            "30123b", "4145ab", "4870d5", "4194e1", "3cb6e4",
            "41d4d8", "51edc1", "6dfeb6", "8fffa9", "affe98",
            "c9fc85", "dff66f", "f0ec5c", "fbdd4b", "fdca3b",
            "feb72d", "fca02a", "f78532", "eb663c", "d93f4c"
        };
        return hex.Select(h => new SolidColorBrush((Color)ColorConverter.ConvertFromString("#" + h))).ToArray();
    }
}
