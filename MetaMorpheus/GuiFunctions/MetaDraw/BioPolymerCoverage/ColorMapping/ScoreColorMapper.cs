namespace GuiFunctions.MetaDraw;

public sealed class ScoreColorMapper : NumericBioPolymerCoverageColorMapper
{
    public override ColorResultsBy ColorBy => ColorResultsBy.Score;
    public override string DisplayName => "Score";

    public override double? GetNumericValue(BioPolymerCoverageResultModel result)
    {
        return result?.Match?.Score;
    }
}
