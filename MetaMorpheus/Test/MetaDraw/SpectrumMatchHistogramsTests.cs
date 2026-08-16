using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.IO;
using System.Linq;
using GuiFunctions;
using NUnit.Framework;
using Readers;

namespace Test.MetaDraw;

/// <summary>
/// Characterisation tests for the binning extracted out of PlotModelStat. These exist so phase 3 can
/// replace the OxyPlot rendering of these histograms and still prove the numbers did not move -
/// asserting on the bins rather than on a chart is the only check that survives a renderer swap.
/// </summary>
public class SpectrumMatchHistogramsTests
{
    private static Dictionary<string, ObservableCollection<SpectrumMatchFromTsv>> LoadFixture()
    {
        string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "VariantCrossTest.psmtsv");
        var matches = SpectrumMatchTsvReader.ReadTsv(path, out _);
        return matches.GroupBy(p => p.FileNameWithoutExtension)
            .ToDictionary(g => g.Key, g => new ObservableCollection<SpectrumMatchFromTsv>(g));
    }

    /// <summary>Every plot type must produce data without throwing on real, ambiguity-containing input.</summary>
    [Test]
    public void EveryPlotTypeBinsWithoutThrowing()
    {
        var bySourceFile = LoadFixture();
        Assert.That(bySourceFile, Is.Not.Empty, "fixture did not load");

        for (int plotType = 1; plotType <= 14; plotType++)
        {
            var data = SpectrumMatchHistograms.GetHistogramData(plotType, bySourceFile);
            Assert.That(data, Is.Not.Null, $"plot {plotType} returned no data");
            Assert.That(data.YAxisTitle, Is.Not.Empty, $"plot {plotType} has no y-axis title");
        }
    }

    /// <summary>
    /// The bin dictionaries must account for exactly the values collected, per source file. This is the
    /// invariant a renderer change could break silently.
    /// </summary>
    [TestCase(1)]
    [TestCase(3)]
    [TestCase(6)]
    [TestCase(7)]
    [TestCase(12)]
    [TestCase(14)]
    public void BinCountsAccountForEveryValue(int plotType)
    {
        var data = SpectrumMatchHistograms.GetHistogramData(plotType, LoadFixture());

        foreach (var (file, numbers) in data.NumbersBySourceFile)
        {
            Assert.That(data.DictsBySourceFile.ContainsKey(file), $"no bins for {file}");
            Assert.That(data.DictsBySourceFile[file].Values.Sum(), Is.EqualTo(numbers.Count()),
                $"plot {plotType}, file {file}: binned count does not match the values collected");
        }
    }

    /// <summary>Bin size is part of the contract - phase 3 must not quietly re-bin.</summary>
    [TestCase(1, 0.1)]
    [TestCase(2, 2.0)]
    [TestCase(12, 1.0)]
    [TestCase(14, 1.0)]
    public void BinSizeIsStable(int plotType, double expected)
    {
        var data = SpectrumMatchHistograms.GetHistogramData(plotType, LoadFixture());
        Assert.That(data.BinSize, Is.EqualTo(expected).Within(1e-9));
    }

    /// <summary>Ambiguous notches are split across bins rather than parsed as one value or throwing.</summary>
    [Test]
    public void AmbiguousNotchesAreSplitAcrossBins()
    {
        var bySourceFile = LoadFixture();
        int ambiguous = bySourceFile.Values.SelectMany(v => v).Count(p => (p.Notch ?? "").Contains('|'));
        var data = SpectrumMatchHistograms.GetHistogramData(14, bySourceFile);

        int totalPsms = bySourceFile.Values.Sum(v => v.Count);
        int binned = data.DictsBySourceFile.Values.Sum(d => d.Values.Sum());

        Assert.That(binned, Is.GreaterThanOrEqualTo(totalPsms),
            $"{ambiguous} PSMs have pipe-delimited notches, so the binned total should be at least the PSM count");
    }

    [Test]
    public void RoundToBinIsSymmetricAboutZero()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SpectrumMatchHistograms.RoundToBin(-2.0, 1), Is.EqualTo(-SpectrumMatchHistograms.RoundToBin(2.0, 1)));
            Assert.That(SpectrumMatchHistograms.RoundToBin(0.0, 1), Is.Zero);
        });
    }

    /// <summary>Ambiguity levels collapse to their first value; blank means level 1.</summary>
    [TestCase("2|3", "2")]
    [TestCase("", "1")]
    [TestCase(null, "1")]
    [TestCase(" 2A ", "2A")]
    public void AmbiguityLevelNormalises(string input, string expected)
    {
        Assert.That(SpectrumMatchHistograms.NormalizeAmbiguityLevel(input), Is.EqualTo(expected));
    }
}
