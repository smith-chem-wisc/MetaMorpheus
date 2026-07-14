using System.Collections.Generic;
using System.Linq;
using System.Windows.Media;

namespace GuiFunctions.MetaDraw;

public sealed class GrayscaleColorGradient : ColorGradient
{
    public override ColorGradientType GradientType => ColorGradientType.Grayscale;

    private static readonly IReadOnlyList<SolidColorBrush> Brushes = BuildBrushes();

    public override IReadOnlyList<SolidColorBrush> GetBrushes() => Brushes;

    private static SolidColorBrush[] BuildBrushes()
    {
        var hex = new[]
        {
            "ffffff", "f2f2f2", "e6e6e6", "dadada", "cecece",
            "c2c2c2", "b6b6b6", "aaaaaa", "9e9e9e", "929292",
            "868686", "7a7a7a", "6e6e6e", "626262", "565656",
            "4a4a4a", "3e3e3e", "323232", "262626", "000000"
        };
        return hex.Select(h => new SolidColorBrush((Color)ColorConverter.ConvertFromString("#" + h))).ToArray();
    }
}
