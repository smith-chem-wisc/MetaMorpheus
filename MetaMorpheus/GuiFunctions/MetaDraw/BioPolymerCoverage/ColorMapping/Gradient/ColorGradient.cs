using System;
using System.Collections.Generic;
using System.Windows.Media;

namespace GuiFunctions.MetaDraw.BioPolymerCoverage.ColorMapping.Gradient;

public abstract class ColorGradient
{
    public abstract ColorGradientType GradientType { get; }
    public virtual int BinCount => 20;

    public abstract IReadOnlyList<SolidColorBrush> GetBrushes();

    public virtual SolidColorBrush GetBrush(double normalizedValue)
    {
        var brushes = GetBrushes();
        int bin = (int)Math.Clamp(normalizedValue * (BinCount - 1), 0, BinCount - 1);
        return brushes[bin];
    }
}
