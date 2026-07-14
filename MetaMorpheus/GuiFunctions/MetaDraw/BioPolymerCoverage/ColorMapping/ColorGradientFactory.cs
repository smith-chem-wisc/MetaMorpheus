using System;

namespace GuiFunctions.MetaDraw;

public static class ColorGradientFactory
{
    public static ColorGradient Create(ColorGradientType gradientType)
    {
        return gradientType switch
        {
            ColorGradientType.Viridis => new ViridisColorGradient(),
            _ => throw new ArgumentOutOfRangeException(nameof(gradientType), gradientType, "Unsupported gradient type")
        };
    }
}
