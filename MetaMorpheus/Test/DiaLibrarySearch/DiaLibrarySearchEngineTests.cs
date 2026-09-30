#nullable enable
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Security.Cryptography;
using EngineLayer;
using EngineLayer.DiaLibrarySearch;
using NUnit.Framework;
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
    private static readonly LinearIrtMap EndpointMap = LinearIrtMap.ThroughPoints(
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

    private DiaLibrarySearchResults Search(SyntheticDiaRun run, IIrtMap map, out string libraryPath)
    {
        libraryPath = run.WriteLibrary(_directory);
        using var library = MslLibrary.Load(libraryPath);
        var engine = new DiaLibrarySearchEngine(run.Scans, library, map, new DiaLibrarySearchParameters(),
            new CommonParameters(), [], []);
        return (DiaLibrarySearchResults)engine.Run();
    }

    private DiaLibrarySearchResults Search(SyntheticDiaRun run) => Search(run, EndpointMap, out _);

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
        var identity = new LinearIrtMap(slope: 1, intercept: 0);

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
}
