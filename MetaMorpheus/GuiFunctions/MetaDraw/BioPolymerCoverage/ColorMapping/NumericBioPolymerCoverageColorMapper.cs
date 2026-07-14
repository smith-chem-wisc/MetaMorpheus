using System.Windows.Media;

namespace GuiFunctions.MetaDraw;

public abstract class NumericBioPolymerCoverageColorMapper : BioPolymerCoverageColorMapper
{
    public override bool IsNumeric => true;

    public abstract double? GetNumericValue(BioPolymerCoverageResultModel result);

    public override SolidColorBrush GetBrush(
        BioPolymerCoverageResultModel result,
        BioPolymerCoverageColorScale? scale)
    {
        var value = GetNumericValue(result);
        if (!value.HasValue)
            return new SolidColorBrush(Colors.Gray);
        if (scale == null)
            return new SolidColorBrush(Colors.Gray);
        return scale.Gradient.GetBrush(scale.Normalize(value.Value));
    }

    public override string GetLegendTitle(BioPolymerCoverageColorScale? scale) => ColorBy.ToString();
}
