using System.Collections.Generic;
using System.Linq;
using EngineLayer;
using EngineLayer.ParallelSearch;
using MassSpectrometry;
using NUnit.Framework;
using Omics;
using Omics.Digestion;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using TaskLayer;

namespace Test.ParallelSearchTask.Engine;

[TestFixture]
public class TransientClassicSearchEngineTests
{
    private const string ProteaseName = "TransientClassicSearchEngineTests protease";

    private static CommonParameters CreateCommonParameters(int threads)
    {
        if (!ProteaseDictionary.Dictionary.ContainsKey(ProteaseName))
        {
            ProteaseDictionary.Dictionary.Add(ProteaseName, new Protease(ProteaseName, CleavageSpecificity.Full, null, null,
                new List<DigestionMotif> { new("K", null, 1, "") }));
        }

        return new CommonParameters(
            dissociationType: DissociationType.HCD,
            digestionParams: new DigestionParams(protease: ProteaseName, minPeptideLength: 1),
            scoreCutoff: 1)
        {
            MaxThreadsToUsePerFile = threads
        };
    }

    private static (SpectralMatch[] Baseline, SpectralMatch[] AfterTransient, TransientSearchEngineResults Results) Search(int threads)
    {
        var commonParameters = CreateCommonParameters(threads);
        var fileSpecificParameters = new List<(string FileName, CommonParameters Parameters)> { ("", commonParameters) };
        var scans = MetaMorpheusTask.GetMs2Scans(new TestDataFile(), null, new CommonParameters())
            .OrderBy(b => b.PrecursorMass).ToArray();
        var searchMode = new SinglePpmAroundZeroSearchMode(5);
        var noMods = new List<Modification>();

        // baseline: the scan is identified from the base database
        SpectralMatch[] baseline = new SpectralMatch[scans.Length];
        new TransientClassicSearchEngine(baseline, scans, noMods, noMods,
            new List<IBioPolymer> { new Protein("MNNNKQQQ", "Base") }, searchMode, commonParameters,
            fileSpecificParameters, new List<string>()).Run();

        // transient: a second protein with the same peptide, searched on a copy of the baseline array
        SpectralMatch[] afterTransient = baseline.ToArray();
        var results = (TransientSearchEngineResults)new TransientClassicSearchEngine(afterTransient, scans, noMods, noMods,
            new List<IBioPolymer> { new Protein("MNNKNKNKQQQ", "Transient") }, searchMode, commonParameters,
            fileSpecificParameters, new List<string>(), copyOnWriteEnabled: true).Run();

        return (baseline, afterTransient, results);
    }

    [TestCase(1)]
    [TestCase(2)]
    public void CopyOnWrite_AddsTheTransientMatchWithoutTouchingTheBaselinePsm(int threads)
    {
        var (baseline, afterTransient, results) = Search(threads);

        Assert.Multiple(() =>
        {
            Assert.That(baseline[0], Is.Not.Null);
            Assert.That(baseline[0].BaseSequence, Is.EqualTo("QQQ"));
            Assert.That(baseline[0].BestMatchingBioPolymersWithSetMods.Select(m => m.SpecificBioPolymer.Parent.Accession),
                Is.EquivalentTo(new[] { "Base" }));

            // the transient search wrote to a clone, so the shared baseline PSM is unchanged
            Assert.That(afterTransient[0], Is.Not.SameAs(baseline[0]));
            Assert.That(afterTransient[0].BestMatchingBioPolymersWithSetMods.Select(m => m.SpecificBioPolymer.Parent.Accession),
                Does.Contain("Transient"));
            Assert.That(results.UpdatedSpectralMatchIndexes, Is.EquivalentTo(new[] { 0 }));
            Assert.That(results.PeptidesSearched, Is.GreaterThan(0));
            Assert.That(results.ToString(), Does.Contain("Number of PSMs updated: 1"));
        });
    }

    [Test]
    public void SingleAndMultiThreadedSearches_Agree()
    {
        var single = Search(1);
        var multi = Search(2);

        Assert.Multiple(() =>
        {
            Assert.That(single.AfterTransient[0].Score, Is.EqualTo(multi.AfterTransient[0].Score).Within(1e-9));
            Assert.That(single.AfterTransient[0].BaseSequence, Is.EqualTo(multi.AfterTransient[0].BaseSequence));
            Assert.That(single.Results.PeptidesSearched, Is.EqualTo(multi.Results.PeptidesSearched));
            Assert.That(single.Results.UpdatedSpectralMatchIndexes, Is.EquivalentTo(multi.Results.UpdatedSpectralMatchIndexes));
        });
    }
}
