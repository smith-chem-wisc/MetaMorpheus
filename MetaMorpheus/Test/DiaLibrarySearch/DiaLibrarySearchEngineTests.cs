#nullable enable
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Security.Cryptography;
using Chromatography.RetentionTimeCalibration;
using EngineLayer;
using EngineLayer.DiaLibrarySearch;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Omics.SpectralMatch.MslSpectralLibrary;
using Readers.SpectralLibrary;

namespace Test.DiaLibrarySearch;

/// <summary>
/// The DIA library search's permanent contract, run on <see cref="SyntheticDiaRun"/>, where the answer is known. Every
/// later milestone deepens the engine; none may break these.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class DiaLibrarySearchEngineTests
{
    private string _directory = "";

    /// <summary>
    /// A linear run-RT to iRT map through the true curve's endpoints (iRT -20 at 0.5 min, iRT 120 at 9.936 min). The
    /// true curve bows away from it by up to about 12 iRT, inside the engine's default iRT window.
    /// </summary>
    private static readonly IrtCalibrationModel EndpointMap = IrtCalibration.Line(
        (new RtMinutes(SyntheticDiaRun.TrueRtMinutes(-20)), new Irt(-20)),
        (new RtMinutes(SyntheticDiaRun.TrueRtMinutes(120)), new Irt(120)));

    [OneTimeSetUp]
    public void SetUp()
    {
        _directory = Path.Combine(TestContext.CurrentContext.TestDirectory, "DiaLibrarySearchEngineTests");
        Directory.CreateDirectory(_directory);
    }

    [OneTimeTearDown]
    public void TearDown()
    {
        Directory.Delete(_directory, true);
    }

    private DiaLibrarySearchResults Search(SyntheticDiaRun run, IrtCalibrationModel map, out string libraryPath)
    {
        libraryPath = run.WriteLibrary(_directory);
        using var library = MslLibrary.Load(libraryPath);
        var engine = new DiaLibrarySearchEngine(run.Scans, library, map, new DiaLibrarySearchParameters(),
            new CommonParameters(), [], []);
        return (DiaLibrarySearchResults)engine.Run();
    }

    private DiaLibrarySearchResults Search(SyntheticDiaRun run) => Search(run, EndpointMap, out _);

    private static List<DiaPrecursorMatch> Search(MslLibrary library, IEnumerable<MsDataScan> scans) =>
        ((DiaLibrarySearchResults)new DiaLibrarySearchEngine(scans.ToArray(), library, EndpointMap, new DiaLibrarySearchParameters(),
            new CommonParameters(), [], []).Run()).Matches.OrderBy(m => m.PrecursorIndex).ToList();

    /// <summary>Three of every four targets elute. At least 90% of those are reported at 1% FDR.</summary>
    [Test]
    public void PlantedPrecursorsAreFound()
    {
        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0);

        var results = Search(run);

        var found = results.Matches.Where(m => !m.IsDecoy && m.QValue <= 0.01).Select(m => m.FullSequence).ToHashSet();
        int planted = run.PlantedSequences.Count;
        Assert.That(planted, Is.GreaterThan(100), "the fixture must plant enough precursors to mean something");
        Assert.That(found.Count(run.PlantedSequences.Contains), Is.GreaterThanOrEqualTo((int)Math.Ceiling(0.9 * planted)));
        Assert.That(found.Count(sequence => !run.PlantedSequences.Contains(sequence)), Is.LessThanOrEqualTo(Math.Max(1, found.Count / 100)),
            "a target that never eluted is a false discovery; at 1% FDR there can be at most about 1% of them");
    }

    /// <summary>A run of noise alone has nothing to find.</summary>
    [Test]
    public void NullRunReportsNothing()
    {
        var run = SyntheticDiaRun.Build(200, _ => false);

        var results = Search(run);

        Assert.That(results.Matches.Count(m => !m.IsDecoy && m.QValue <= 0.01), Is.LessThanOrEqualTo(Math.Max(1, results.Matches.Count / 100)));
    }

    /// <summary>When only the decoys elute, no target may pass: every real signal belongs to a decoy.</summary>
    [Test]
    public void ARunOfDecoysReportsNoTargets()
    {
        var run = SyntheticDiaRun.Build(200, entry => entry.IsDecoy);

        var results = Search(run);

        Assert.That(results.Matches.Count(m => !m.IsDecoy && m.QValue <= 0.01), Is.Zero);
    }

    /// <summary>With no decoys, target-decoy q-values are all zero, which reads as "every target is right".</summary>
    [Test]
    public void ALibraryWithoutDecoysThrows()
    {
        var run = SyntheticDiaRun.Build(50, entry => true, withDecoys: false);

        var e = Assert.Throws<MetaMorpheusException>(() => Search(run));
        Assert.That(e!.Message, Does.Contain("decoy"));
    }

    [Test]
    public void TargetQValuesNeverDecreaseAsScoreFalls()
    {
        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 2) == 0);

        var targets = Search(run).Matches.Where(m => !m.IsDecoy).OrderByDescending(m => m.Score).ToList();

        Assert.That(targets, Is.Not.Empty);
        for (int i = 1; i < targets.Count; i++)
            Assert.That(targets[i].QValue, Is.GreaterThanOrEqualTo(targets[i - 1].QValue), $"at rank {i}");
        Assert.That(targets.All(m => m.QValue is >= 0 and <= 1));
    }

    [Test]
    public void TargetAndDecoyCountsAreReported()
    {
        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0);

        var results = Search(run);

        int targets = results.Matches.Count(m => !m.IsDecoy);
        int decoys = results.Matches.Count(m => m.IsDecoy);
        int passing = results.Matches.Count(m => !m.IsDecoy && m.QValue <= 0.01);
        Assert.That(results.TargetCount, Is.EqualTo(targets));
        Assert.That(results.DecoyCount, Is.EqualTo(decoys));
        Assert.That(results.ToString(), Does.Contain($"Target precursors: {targets}"));
        Assert.That(results.ToString(), Does.Contain($"Decoy precursors: {decoys}"));
        Assert.That(results.ToString(), Does.Contain($"Target precursors with q-value <= 0.01: {passing}"));
    }

    /// <summary>
    /// iRT and run minutes are different quantities. The types keep them apart, with no conversion in either
    /// direction, not even to double, so the only way across is a map. A search told the identity map (iRT = minutes)
    /// looks in the wrong place and finds almost nothing.
    /// </summary>
    [Test]
    public void IrtAndMinutesDoNotMix()
    {
        foreach (var type in new[] { typeof(Irt), typeof(RtMinutes) })
        {
            var conversions = type.GetMethods(BindingFlags.Public | BindingFlags.Static)
                .Where(method => method.Name is "op_Implicit" or "op_Explicit");
            Assert.That(conversions, Is.Empty, $"{type.Name} must not convert implicitly or explicitly");
        }

        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0);
        var identity = IrtCalibration.Line((new RtMinutes(0), new Irt(0)), (new RtMinutes(1), new Irt(1)));

        var results = Search(run, identity, out _);

        int found = results.Matches.Count(m => !m.IsDecoy && m.QValue <= 0.01 && run.PlantedSequences.Contains(m.FullSequence));
        Assert.That(found, Is.LessThan(run.PlantedSequences.Count / 2));
    }

    /// <summary>The run-to-iRT map is applied to the run. The library is never rewritten in minutes.</summary>
    [Test]
    public void TheLibraryFileIsUnchangedBySearch()
    {
        var run = SyntheticDiaRun.Build(50, entry => !entry.IsDecoy);
        string path = run.WriteLibrary(_directory);
        byte[] before = SHA256.HashData(File.ReadAllBytes(path));

        using (var library = MslLibrary.Load(path))
            new DiaLibrarySearchEngine(run.Scans, library, EndpointMap, new DiaLibrarySearchParameters(), new CommonParameters(), [], []).Run();

        Assert.That(SHA256.HashData(File.ReadAllBytes(path)), Is.EqualTo(before));
    }

    [Test]
    public void ReportedApexIsInIrtAndInMinutes()
    {
        var run = SyntheticDiaRun.Build(100, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0);

        var best = Search(run).Matches.Where(m => !m.IsDecoy && run.PlantedSequences.Contains(m.FullSequence))
            .OrderByDescending(m => m.Score).First();

        var entry = run.Library.Single(e => e.FullSequence == best.FullSequence);
        Assert.That(best.LibraryIrt.Value, Is.EqualTo(entry.RetentionTime).Within(1e-3));
        Assert.That(best.ApexRt.Value, Is.EqualTo(SyntheticDiaRun.TrueRtMinutes(entry.RetentionTime)).Within(SyntheticDiaRun.CycleMinutes));
        Assert.That(best.ApexIrt, Is.EqualTo(EndpointMap.ToIrt(best.ApexRt)));
    }

    [Test]
    public void EachMatchCarriesItsLibraryPrecursor()
    {
        var run = SyntheticDiaRun.Build(100, entry => SyntheticDiaRun.Bucket(entry, 4) != 0);

        var matches = Search(run).Matches;

        Assert.That(matches, Is.Not.Empty);
        foreach (var match in matches)
        {
            var entry = run.Library.Single(e => e.FullSequence == match.FullSequence);
            Assert.That(match.Charge, Is.EqualTo(entry.ChargeState), match.FullSequence);
            Assert.That(match.PrecursorMz, Is.EqualTo(entry.PrecursorMz).Within(1e-3), match.FullSequence);
            Assert.That(match.IsDecoy, Is.EqualTo(entry.IsDecoy), match.FullSequence);
            Assert.That(match.LibraryIrt.Value, Is.EqualTo(entry.RetentionTime).Within(1e-3), match.FullSequence);
        }
    }

    /// <summary>
    /// A precursor isolated by two windows is reported once, with the better of its two scores. The second window
    /// covers the first shifted by 0.25 m/z and holds the same scans with every intensity squared, so a precursor's
    /// fragment pattern, and its score, differ between them.
    /// </summary>
    [Test]
    public void APrecursorIsolatedByTwoWindowsKeepsItsBestScore()
    {
        var run = SyntheticDiaRun.Build(200, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0);
        double firstWindowCentre = SyntheticDiaRun.FirstWindowLowMz + SyntheticDiaRun.WindowWidth / 2;
        int scanNumber = run.Scans.Length + 1;
        var shifted = run.Scans.Where(scan => scan.IsolationMz == firstWindowCentre)
            .Select(scan => SyntheticDiaRun.Scan(
                new MzSpectrum(scan.MassSpectrum.XArray.ToArray(), scan.MassSpectrum.YArray.Select(y => y * y).ToArray(), false),
                scanNumber++, scan.RetentionTime, firstWindowCentre + 0.25))
            .ToList();
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));

        var first = Search(library, run.Scans).ToDictionary(m => m.PrecursorIndex, m => m.Score);
        var second = Search(library, shifted).ToDictionary(m => m.PrecursorIndex, m => m.Score);
        var both = Search(library, run.Scans.Concat(shifted));

        var inBoth = second.Keys.Where(first.ContainsKey).ToList();
        Assert.That(inBoth.Count(i => second[i] > first[i]), Is.Positive, "some precursors must score better in the second window");
        Assert.That(inBoth.Count(i => second[i] < first[i]), Is.Positive, "and some worse");
        Assert.That(both.Select(m => m.PrecursorIndex), Is.Unique);
        var expected = first.Keys.Union(second.Keys).ToDictionary(i => i,
            i => Math.Max(first.GetValueOrDefault(i, double.NegativeInfinity), second.GetValueOrDefault(i, double.NegativeInfinity)));
        Assert.That(both.ToDictionary(m => m.PrecursorIndex, m => m.Score), Is.EquivalentTo(expected));
    }

    /// <summary>
    /// Only MS2 scans with an isolation range are searched, and a scan with no peaks contributes nothing. Adding an MS1
    /// scan and an MS2 scan without an isolation range, both full of library fragments, and an empty MS2 scan inside a
    /// window changes no match.
    /// </summary>
    [Test]
    public void ScansThatAreNotDiaWindowsChangeNothing()
    {
        var run = SyntheticDiaRun.Build(100, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0);
        var richest = run.Scans.MaxBy(scan => scan.MassSpectrum.SumOfAllY)!;
        int scanNumber = run.Scans.Length + 1;
        var extra = new[]
        {
            SyntheticDiaRun.Scan(richest.MassSpectrum, scanNumber++, richest.RetentionTime, isolationMz: null, msnOrder: 1),
            SyntheticDiaRun.Scan(richest.MassSpectrum, scanNumber++, richest.RetentionTime, isolationMz: null),
            SyntheticDiaRun.Scan(new MzSpectrum([], [], false), scanNumber, richest.RetentionTime, richest.IsolationMz),
        };
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));

        var plain = Search(library, run.Scans);
        var withExtra = Search(library, run.Scans.Concat(extra));

        Assert.That(plain, Is.Not.Empty);
        Assert.That(withExtra, Is.EqualTo(plain));
    }

    /// <summary>
    /// An entry with no fragments has nothing to score and is skipped. mzLib's validator refuses such an entry, but
    /// <see cref="MslLibrary.Save"/> does not run it, so a library written elsewhere can hold one.
    /// </summary>
    [Test]
    public void AnEntryWithoutFragmentsIsNotScored()
    {
        var run = SyntheticDiaRun.Build(50, entry => !entry.IsDecoy);
        var template = run.Library.First(entry => !entry.IsDecoy);
        var bare = new MslLibraryEntry
        {
            FullSequence = "PEPTIDEPEPTIDEK",
            BaseSequence = "PEPTIDEPEPTIDEK",
            PrecursorMz = template.PrecursorMz,
            ChargeState = template.ChargeState,
            RetentionTime = template.RetentionTime,
            IsDecoy = false,
        };
        string path = Path.Combine(_directory, $"bare_{Guid.NewGuid():N}.msl");
        MslLibrary.Save(path, run.Library.Append(bare).ToList());
        using var library = MslLibrary.Load(path);

        var sequences = Search(library, run.Scans).Select(m => m.FullSequence).ToList();

        Assert.That(sequences, Does.Contain(template.FullSequence), "its neighbour, at the same m/z and iRT, is still scored");
        Assert.That(sequences, Does.Not.Contain(bare.FullSequence));
    }

    /// <summary>Stop, pressed before the search reaches a window, leaves every window unsearched.</summary>
    [Test]
    [NonParallelizable]
    public void AStoppedSearchSearchesNoWindow()
    {
        var run = SyntheticDiaRun.Build(50, entry => !entry.IsDecoy);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));
        bool stopLoops = GlobalVariables.StopLoops;
        List<DiaPrecursorMatch> matches;
        try
        {
            GlobalVariables.StopLoops = true;
            matches = Search(library, run.Scans);
        }
        finally
        {
            GlobalVariables.StopLoops = stopLoops;
        }

        Assert.That(matches, Is.Empty);
        Assert.That(Search(library, run.Scans), Is.Not.Empty, "the same search, not stopped, finds precursors");
    }

    [TestCase("scans")]
    [TestCase("library")]
    [TestCase("irtMap")]
    [TestCase("parameters")]
    public void AMissingArgumentIsRefusedByName(string argument)
    {
        var run = SyntheticDiaRun.Build(5, _ => false);
        using var library = MslLibrary.Load(run.WriteLibrary(_directory));

        var e = Assert.Throws<ArgumentNullException>(() => new DiaLibrarySearchEngine(
            argument == "scans" ? null! : run.Scans,
            argument == "library" ? null! : library,
            argument == "irtMap" ? null! : EndpointMap,
            argument == "parameters" ? null! : new DiaLibrarySearchParameters(),
            new CommonParameters(), [], []));
        Assert.That(e!.ParamName, Is.EqualTo(argument));
    }
}
