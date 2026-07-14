using System.Collections.Generic;
using System.Windows.Media;

namespace GuiFunctions.MetaDraw;

public abstract class BioPolymerCoverageColorMapper
{
    public abstract ColorResultsBy ColorBy { get; }
    public abstract bool IsNumeric { get; }
    public abstract bool SupportsGradientSelection { get; }
    public abstract bool SupportsLogScale { get; }
    public virtual bool DefaultUseLogScale => false;

    public abstract void Prepare(
        IReadOnlyList<BioPolymerCoverageResultModel> filteredResults,
        ColorGradientType gradientType,
        bool useLogScale);

    public abstract SolidColorBrush GetBrush(BioPolymerCoverageResultModel result);

    public virtual IReadOnlyList<(SolidColorBrush Brush, string Label)> LegendItems => [];
    public virtual string? GradientLegendTitle => null;
    public virtual IReadOnlyList<SolidColorBrush>? GradientBrushes => null;
    public virtual double? GradientMinValue => null;
    public virtual double? GradientMaxValue => null;
}
