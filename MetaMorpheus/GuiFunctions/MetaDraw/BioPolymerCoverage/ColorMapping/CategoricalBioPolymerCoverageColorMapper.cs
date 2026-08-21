using GuiFunctions.MetaDraw.BioPolymerCoverage.ColorMapping.Gradient;
using System;
using System.Collections.Generic;
using System.Windows.Media;

namespace GuiFunctions.MetaDraw;

public sealed class CategoricalBioPolymerCoverageColorMapper : BioPolymerCoverageColorMapper
{
    private readonly Func<BioPolymerCoverageResultModel, SolidColorBrush> _brushSelector;
    private readonly Func<IReadOnlyList<BioPolymerCoverageResultModel>, IReadOnlyList<(SolidColorBrush Brush, string Label)>> _legendBuilder;
    private IReadOnlyList<(SolidColorBrush Brush, string Label)> _preparedLegendItems = [];

    public CategoricalBioPolymerCoverageColorMapper(
        ColorResultsBy colorBy,
        Func<BioPolymerCoverageResultModel, SolidColorBrush> brushSelector,
        Func<IReadOnlyList<BioPolymerCoverageResultModel>, IReadOnlyList<(SolidColorBrush Brush, string Label)>>? legendBuilder = null)
    {
        _colorBy = colorBy;
        _brushSelector = brushSelector ?? throw new ArgumentNullException(nameof(brushSelector));
        _legendBuilder = legendBuilder ?? (_ => []);
    }

    private readonly ColorResultsBy _colorBy;
    public override ColorResultsBy ColorBy => _colorBy;
    public override bool IsNumeric => false;
    public override bool SupportsGradientSelection => false;
    public override bool SupportsLogScale => false;

    public override void Prepare(
        IReadOnlyList<BioPolymerCoverageResultModel> filteredResults,
        ColorGradientType gradientType,
        bool useLogScale)
    {
        _preparedLegendItems = _legendBuilder(filteredResults);
    }

    public override SolidColorBrush GetBrush(BioPolymerCoverageResultModel result)
    {
        return _brushSelector(result) ?? new SolidColorBrush(Colors.Gray);
    }

    public override IReadOnlyList<(SolidColorBrush Brush, string Label)> LegendItems => _preparedLegendItems;
}
