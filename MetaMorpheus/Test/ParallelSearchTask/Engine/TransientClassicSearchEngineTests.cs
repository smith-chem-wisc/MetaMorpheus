using System.Collections.Generic;
using System.Linq;
using EngineLayer;
using EngineLayer.FdrAnalysis;
using EngineLayer.ParallelSearch;
using EngineLayer.ParallelSearch.FdrAlignment;
using EngineLayer.SpectrumMatch;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Omics;
using Omics.Digestion;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Readers;
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
            scoreCutoff: 1,
            maxThreadsToUsePerFile: threads);
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

    [Test]
    public void PeptideWithAnUnknownResidue_MatchesNoScan()
    {
        // XQQQ has a NaN mass (X has none) but its y ions are those of the scan's QQQ, so it would score if it reached the scan
        var commonParameters = CreateCommonParameters(1);
        var scans = MetaMorpheusTask.GetMs2Scans(new TestDataFile(), null, new CommonParameters())
            .OrderBy(b => b.PrecursorMass).ToArray();
        var noMods = new List<Modification>();
        var psms = new SpectralMatch[scans.Length];

        var results = (TransientSearchEngineResults)new TransientClassicSearchEngine(psms, scans, noMods, noMods,
            new List<IBioPolymer> { new Protein("MNNNKXQQQ", "Unknown") }, new SinglePpmAroundZeroSearchMode(5),
            commonParameters, new List<(string FileName, CommonParameters Parameters)> { ("", commonParameters) },
            new List<string>()).Run();

        Assert.Multiple(() =>
        {
            Assert.That(psms, Has.All.Null);
            Assert.That(results.UpdatedSpectralMatchIndexes, Is.Empty);
        });
    }

    private static SpectralMatch[] SearchWithFilePrecursorTolerance(double fileTolerancePpm, double otherFileTolerancePpm, out PerFileTolerances perFile)
    {
        var commonParameters = CreateCommonParameters(1);
        var scans = MetaMorpheusTask.GetMs2Scans(new TestDataFile(), null, new CommonParameters())
            .OrderBy(b => b.PrecursorMass).ToArray();
        CommonParameters FileParams(double ppm) => new(precursorMassTolerance: new PpmTolerance(ppm),
            dissociationType: DissociationType.HCD, digestionParams: commonParameters.DigestionParams, scoreCutoff: 1);
        var perFileParameters = new List<(string FileName, CommonParameters Parameters)>
        {
            ("", FileParams(fileTolerancePpm)), ("other", FileParams(otherFileTolerancePpm))
        };
        perFile = PerFileTolerances.Create(scans, perFileParameters, commonParameters, out _);

        var psms = new SpectralMatch[scans.Length];
        var noMods = new List<Modification>();
        // the acceptor is built from the widest tolerance, as the task does
        new TransientClassicSearchEngine(psms, scans, noMods, noMods, new List<IBioPolymer> { new Protein("MNNNKQQQ", "Base") },
            new SinglePpmAroundZeroSearchMode(perFile.WidestPrecursor.Value), commonParameters, perFileParameters,
            new List<string>(), perFileTolerances: perFile).Run();
        return psms;
    }

    [Test]
    public void PerFilePrecursorTolerance_NarrowsTheWidestWindowForEachScan()
    {
        // the scan's file allows 5 ppm even though another file widens the window to 10 ppm: still matched
        Assert.That(SearchWithFilePrecursorTolerance(5, 10, out var loose)[0], Is.Not.Null);
        Assert.That(loose.PrecursorVaries, Is.True);

        // the scan's file is effectively exact, so the window widened by the other file must not admit it
        Assert.That(SearchWithFilePrecursorTolerance(1e-9, 10, out _)[0], Is.Null);
    }

    [Test]
    public void PerFileTolerances_AreNotCreatedWhenEveryFileMatchesTheTask()
    {
        var commonParameters = CreateCommonParameters(1);
        var scans = MetaMorpheusTask.GetMs2Scans(new TestDataFile(), null, new CommonParameters()).ToArray();
        var fileParameters = new List<(string FileName, CommonParameters Parameters)> { ("", commonParameters) };

        Assert.That(PerFileTolerances.Create(scans, fileParameters, commonParameters, out var warning), Is.Null);
        Assert.That(warning, Is.Null);
    }

    [Test]
    public void PerFileTolerances_IgnoreTheTaskUnitsWhenEveryFileHasItsOwnTolerance()
    {
        var taskParameters = new CommonParameters(precursorMassTolerance: new AbsoluteTolerance(0.01));
        var scans = MetaMorpheusTask.GetMs2Scans(new TestDataFile(), null, new CommonParameters()).ToArray();
        var fileParameters = new List<(string FileName, CommonParameters Parameters)>
        {
            ("", new CommonParameters(precursorMassTolerance: new PpmTolerance(3))),
            ("other", new CommonParameters(precursorMassTolerance: new PpmTolerance(4)))
        };

        var perFile = PerFileTolerances.Create(scans, fileParameters, taskParameters, out var warning);

        Assert.Multiple(() =>
        {
            Assert.That(warning, Is.Null);
            Assert.That(perFile, Is.Not.Null);
            Assert.That(perFile!.PrecursorIsPpm, Is.True);
            Assert.That(perFile.ScanPrecursorHalfWidth, Is.EqualTo(new[] { 3e-6 }).Within(1e-15));
        });
    }

    /// <summary>
    /// Regression: with a notch-tolerant acceptor, a baseline PSM that FDR left notch-ambiguous carries notch q-values on
    /// each hypothesis. A transient search that ties it adds a hypothesis without them, and ResolveAllAmbiguities used to
    /// read the missing value ("Nullable object must have a value") and stop the whole ParallelSearch.
    /// </summary>
    [Test]
    public void OneMmTransientTie_OnANotchAmbiguousBaselinePsmAfterFdr_DoesNotThrow()
    {
        var (baselinePsm, afterTransient, _, _, _) = OneMmTieSearch();

        Assert.Multiple(() =>
        {
            Assert.That(afterTransient[0], Is.Not.SameAs(baselinePsm));
            Assert.That(afterTransient[0].Notch, Is.Null);
            Assert.That(afterTransient[0].BestMatchingBioPolymersWithSetMods.Select(h => h.SpecificBioPolymer.Parent.Accession),
                Is.EquivalentTo(new[] { "BaseE", "BaseQ", "Transient" }));
            Assert.That(baselinePsm.BestMatchingBioPolymersWithSetMods.Select(h => h.SpecificBioPolymer.Parent.Accession),
                Is.EquivalentTo(new[] { "BaseE", "BaseQ" }));
        });
    }

    /// <summary>
    /// The notch q-value decides whether a match counts as confident (ConfidenceCutoff, FilteredPsms). A match whose tied
    /// hypotheses sit at different notches must get its best notch's value, and the value written must be the one counted.
    /// </summary>
    [Test]
    public void OneMmNotchAmbiguousMatches_GetTheirBestNotchQValue_InTheBaselineAndTheTransientCopy()
    {
        var (baselinePsm, afterTransient, baselineList, before, _) = OneMmTieSearch();
        // FDR alone leaves the match's own notch q-value unset (2), so it could never pass a q-value cutoff
        Assert.That(before, Is.EqualTo(2));

        BaselineFdr.UseBestNotch(baselineList);
        var expected = baselinePsm.BestMatchingBioPolymersWithSetMods.Min(h => h.QValueNotch!.Value);
        Assert.That(baselinePsm.PsmFdrInfo.QValueNotch, Is.EqualTo(expected));
        Assert.That(expected, Is.LessThanOrEqualTo(1));

        // the transient copy is aligned the way ParallelSearchTask aligns it
        var transient = afterTransient[0];
        BaselineFdr.AlignPsms(BaselineFdr.CreatePsmLookup(baselineList), new[] { transient },
            BaselineFdr.CreateNotchLookups(baselineList, peptideLevel: false));

        var written = new Dictionary<string, string>();
        PsmTsvWriter.AddMatchScoreData(written, transient, writePeptideLevelFdr: false);

        Assert.Multiple(() =>
        {
            Assert.That(transient.PsmFdrInfo.QValueNotch, Is.EqualTo(expected));
            Assert.That(transient.PsmFdrInfo, Is.Not.SameAs(baselinePsm.PsmFdrInfo));
            Assert.That(written[SpectrumMatchFromTsvHeader.QValueNotch],
                Is.EqualTo(transient.PsmFdrInfo.QValueNotch.ToString("F6", System.Globalization.CultureInfo.InvariantCulture)));
            // the hypotheses the transient copy shares with the baseline match are not written to
            Assert.That(baselinePsm.BestMatchingBioPolymersWithSetMods.Select(h => h.QValueNotch), Has.All.Not.Null);
        });
    }

    /// <summary>
    /// QQE (403.1703) and QQQ (402.1863, ~47 ppm off OneMM's +1.0029 Da notch) tie on one scan in the baseline, so the
    /// baseline match is notch-ambiguous after FDR; a transient copy of QQQ then ties it again.
    /// </summary>
    private static (SpectralMatch BaselinePsm, SpectralMatch[] AfterTransient, List<SpectralMatch> BaselineList,
        double QValueNotchAfterFdr, MassDiffAcceptor Acceptor) OneMmTieSearch()
    {
        var commonParameters = CreateCommonParameters(1);
        var fileSpecificParameters = new List<(string FileName, CommonParameters Parameters)> { ("", commonParameters) };
        // OneMM's second notch is +1.0029 Da, so QQQ (402.1863) sits ~47 ppm off notch 1 of a 403.1703 (QQE) precursor
        var acceptor = SearchTask.GetMassDiffAcceptor(new PpmTolerance(100), MassDiffAcceptorType.OneMM, null);
        var noMods = new List<Modification>();

        // y1, b2 (shared) and y2 of QQQ and of QQE at equal intensity, so the two peptides tie on score
        double[] mz = { 147.07642, 148.06043, 257.12444, 275.13500, 276.11901 };
        double[] intensities = { 100, 100, 100, 100, 100 };
        var scans = MetaMorpheusTask.GetMs2Scans(new TestDataFile(mz, intensities, 403.17031, 1), null, commonParameters)
            .OrderBy(b => b.PrecursorMass).ToArray();

        SpectralMatch[] baseline = new SpectralMatch[scans.Length];
        new TransientClassicSearchEngine(baseline, scans, noMods, noMods,
            new List<IBioPolymer> { new Protein("MNNNKQQE", "BaseE"), new Protein("MNNNKQQQ", "BaseQ") }, acceptor,
            commonParameters, fileSpecificParameters, new List<string>()).Run();
        var baselineList = baseline.Where(p => p != null).ToList();
        new FdrAnalysisEngine(baselineList, acceptor.NumNotches, commonParameters,
            fileSpecificParameters, new List<string>(), "PSM", false).Run();

        SpectralMatch baselinePsm = baseline[0];
        Assert.That(baselinePsm, Is.Not.Null);
        Assert.That(baselinePsm.Notch, Is.Null, "the baseline PSM must be notch-ambiguous for this test to mean anything");
        Assert.That(baselinePsm.BestMatchingBioPolymersWithSetMods.Select(h => h.QValueNotch), Has.All.Not.Null);
        double qValueNotchAfterFdr = baselinePsm.PsmFdrInfo.QValueNotch;

        // the transient hypothesis has no notch q-value; resolving the PSM used to read it and throw
        SpectralMatch[] afterTransient = baseline.ToArray();
        Assert.DoesNotThrow(() => new TransientClassicSearchEngine(afterTransient, scans, noMods, noMods,
            new List<IBioPolymer> { new Protein("MNNNKQQQ", "Transient") }, acceptor, commonParameters,
            fileSpecificParameters, new List<string>(), copyOnWriteEnabled: true).Run());

        return (baselinePsm, afterTransient, baselineList, qValueNotchAfterFdr, acceptor);
    }
}
