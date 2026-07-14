using System.Collections.Generic;
using System.Windows.Media;

namespace GuiFunctions.MetaDraw;

public abstract class ColorGradient
{
    public abstract ColorGradientType GradientType { get; }
    public abstract int BinCount { get; }
    public abstract SolidColorBrush GetBrush(double normalizedValue);
    public abstract IReadOnlyList<SolidColorBrush> GetBrushes();
}
