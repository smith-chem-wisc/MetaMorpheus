using System.Collections.Generic;
using System.Linq;
using System.Windows.Media;

namespace GuiFunctions.MetaDraw;

public sealed class InfernoColorGradient : ColorGradient
{
    public override ColorGradientType GradientType => ColorGradientType.Inferno;

    private static readonly IReadOnlyList<SolidColorBrush> Brushes = BuildBrushes();

    public override IReadOnlyList<SolidColorBrush> GetBrushes() => Brushes;

    private static SolidColorBrush[] BuildBrushes()
    {
        var hex = new[]
        {
            "000004", "0b0724", "1f0c48", "33106b", "4a0c6b",
            "61136e", "781c6d", "8f2469", "a52c60", "bc3754",
            "d24742", "e2512a", "f06511", "f77b08", "fb9107",
            "fca50a", "fcbb1f", "f9d23c", "f2ea68", "fcffa4"
        };
        return hex.Select(h => new SolidColorBrush((Color)ColorConverter.ConvertFromString("#" + h))).ToArray();
    }
}
