using System.Windows.Media;

namespace GuiFunctions.MetaDraw;

public abstract class BioPolymerCoverageColorMapper
{
    public abstract ColorResultsBy ColorBy { get; }
    public abstract bool IsNumeric { get; }

    public abstract SolidColorBrush GetBrush(
        BioPolymerCoverageResultModel result,
        BioPolymerCoverageColorScale? scale);

    public virtual string GetLegendTitle(BioPolymerCoverageColorScale? scale) => ColorBy.ToString();
}
