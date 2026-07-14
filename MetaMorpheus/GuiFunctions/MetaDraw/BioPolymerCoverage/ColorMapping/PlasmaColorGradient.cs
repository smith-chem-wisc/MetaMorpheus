using System.Collections.Generic;
using System.Linq;
using System.Windows.Media;

namespace GuiFunctions.MetaDraw;

public sealed class PlasmaColorGradient : ColorGradient
{
    public override ColorGradientType GradientType => ColorGradientType.Plasma;

    private static readonly IReadOnlyList<SolidColorBrush> Brushes = BuildBrushes();

    public override IReadOnlyList<SolidColorBrush> GetBrushes() => Brushes;

    private static SolidColorBrush[] BuildBrushes()
    {
        var hex = new[]
        {
            "0d0887", "2c0594", "41049d", "5502a3", "6600a7",
            "7801a8", "8b0aa5", "9c179e", "ad2494", "bc3587",
            "ca457a", "d7566c", "e26760", "ea7b5a", "f18f57",
            "f6a45a", "f9ba5f", "fbd06a", "fbe67b", "f0f921"
        };
        return hex.Select(h => new SolidColorBrush((Color)ColorConverter.ConvertFromString("#" + h))).ToArray();
    }
}
