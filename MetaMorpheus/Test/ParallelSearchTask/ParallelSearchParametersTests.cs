using EngineLayer;
using NUnit.Framework;
using TaskLayer;
using TaskLayer.ParallelSearch;
using UsefulProteomicsDatabases;

namespace Test.ParallelSearchTask;

[TestFixture]
public class ParallelSearchParametersTests
{
    [Test]
    public void Constructor_SetsParallelSearchDefaults()
    {
        var parameters = new ParallelSearchParameters();

        Assert.Multiple(() =>
        {
            Assert.That(parameters.OverwriteTransientSearchOutputs, Is.True);
            Assert.That(parameters.MaxSearchesInParallel, Is.EqualTo(4));
            Assert.That(parameters.WriteTransientResultsOnly, Is.True);
            Assert.That(parameters.WriteTransientSpectralLibrary, Is.False);
            Assert.That(parameters.DoParsimony, Is.True);
            Assert.That(parameters.NoOneHitWonders, Is.True);
            Assert.That(parameters.MassDiffAcceptorType, Is.EqualTo(MassDiffAcceptorType.Exact));
            Assert.That(parameters.SearchType, Is.EqualTo(SearchType.Classic));
            Assert.That(parameters.DoLabelFreeQuantification, Is.False);
            Assert.That(parameters.DoMultiplexQuantification, Is.False);
        });
    }

    [Test]
    public void CopyConstructor_CopiesBaseSearchParameterValues()
    {
        var source = new SearchParameters
        {
            DoParsimony = false,
            NoOneHitWonders = false,
            SearchTarget = false,
            DecoyType = DecoyType.None,
            MassDiffAcceptorType = MassDiffAcceptorType.OneMM,
            SearchType = SearchType.Modern,
            DoLocalizationAnalysis = false,
            WriteDecoys = false,
            WriteContaminants = false,
        };

        var parameters = new ParallelSearchParameters(source);

        Assert.Multiple(() =>
        {
            Assert.That(parameters.DoParsimony, Is.True);
            Assert.That(parameters.NoOneHitWonders, Is.True);
            Assert.That(parameters.SearchTarget, Is.False);
            Assert.That(parameters.DecoyType, Is.EqualTo(DecoyType.None));
            Assert.That(parameters.MassDiffAcceptorType, Is.EqualTo(MassDiffAcceptorType.Exact));
            Assert.That(parameters.SearchType, Is.EqualTo(SearchType.Classic));
            Assert.That(parameters.DoLocalizationAnalysis, Is.False);
            Assert.That(parameters.WriteDecoys, Is.False);
            Assert.That(parameters.WriteContaminants, Is.False);
            Assert.That(parameters.WriteTransientResultsOnly, Is.True);
        });
    }
}
