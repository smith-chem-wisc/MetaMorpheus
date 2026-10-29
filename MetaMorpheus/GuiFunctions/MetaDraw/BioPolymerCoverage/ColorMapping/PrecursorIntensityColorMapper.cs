namespace GuiFunctions.MetaDraw;

public sealed class PrecursorIntensityColorMapper : NumericBioPolymerCoverageColorMapper
{
    public override ColorResultsBy ColorBy => ColorResultsBy.PrecursorIntensity;
    public override string DisplayName => "Precursor Intensity";
    public override bool DefaultUseLogScale => true;

    public override double? GetNumericValue(BioPolymerCoverageResultModel result)
    {
        return result?.Match?.PrecursorIntensity;
    }
}
