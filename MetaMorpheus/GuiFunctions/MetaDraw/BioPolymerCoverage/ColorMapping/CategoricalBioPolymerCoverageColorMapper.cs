using System;
using System.Windows.Media;

namespace GuiFunctions.MetaDraw;

public sealed class CategoricalBioPolymerCoverageColorMapper : BioPolymerCoverageColorMapper
{
    private readonly Func<BioPolymerCoverageResultModel, SolidColorBrush> _brushSelector;

    public CategoricalBioPolymerCoverageColorMapper(
        ColorResultsBy colorBy,
        Func<BioPolymerCoverageResultModel, SolidColorBrush> brushSelector)
    {
        _colorBy = colorBy;
        _brushSelector = brushSelector ?? throw new ArgumentNullException(nameof(brushSelector));
    }

    private readonly ColorResultsBy _colorBy;
    public override ColorResultsBy ColorBy => _colorBy;
    public override bool IsNumeric => false;

    public override SolidColorBrush GetBrush(
        BioPolymerCoverageResultModel result,
        BioPolymerCoverageColorScale? scale)
    {
        return _brushSelector(result) ?? new SolidColorBrush(Colors.Gray);
    }
}
