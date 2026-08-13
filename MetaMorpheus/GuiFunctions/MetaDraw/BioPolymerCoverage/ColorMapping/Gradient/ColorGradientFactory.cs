using System;

namespace GuiFunctions.MetaDraw.BioPolymerCoverage.ColorMapping.Gradient;

public static class ColorGradientFactory
{
    public static ColorGradient Create(ColorGradientType gradientType)
    {
        return gradientType switch
        {
            ColorGradientType.Viridis => new ViridisColorGradient(),
            ColorGradientType.Plasma => new PlasmaColorGradient(),
            ColorGradientType.Inferno => new InfernoColorGradient(),
            ColorGradientType.Turbo => new TurboColorGradient(),
            ColorGradientType.Grayscale => new GrayscaleColorGradient(),
            _ => throw new ArgumentOutOfRangeException(nameof(gradientType), gradientType, "Unsupported gradient type")
        };
    }
}
