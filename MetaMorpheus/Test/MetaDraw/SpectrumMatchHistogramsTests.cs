using System;
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
    [TestCase(14)]
    public void BinCountsAccountForEveryValue(int plotType)
    {
        var bySourceFile = LoadFixture();
        var data = SpectrumMatchHistograms.GetHistogramData(plotType, bySourceFile);

        Assert.That(data.NumbersBySourceFile.Keys, Is.EquivalentTo(bySourceFile.Keys));
        foreach (var (file, numbers) in data.NumbersBySourceFile)
        {
            Assert.That(data.DictsBySourceFile.ContainsKey(file), $"no bins for {file}");
            Assert.That(data.DictsBySourceFile[file].Values.Sum(), Is.EqualTo(numbers.Count()),
                $"plot {plotType}, file {file}: binned count does not match the values collected");
        }
    }

    /// <summary>
    /// Plots 12 and 13 bin PSMs directly and leave NumbersBySourceFile empty, so the check above cannot
    /// reach them. Each PSM lands in exactly one bin, so the bins must sum to the PSM count per file.
    /// </summary>
    [TestCase(12)]
    [TestCase(13)]
    public void BinCountsAccountForEveryPsm(int plotType)
    {
        var bySourceFile = LoadFixture();
        var data = SpectrumMatchHistograms.GetHistogramData(plotType, bySourceFile);

        Assert.That(data.DictsBySourceFile.Keys, Is.EquivalentTo(bySourceFile.Keys));
        foreach (var (file, psms) in bySourceFile)
        {
            Assert.That(data.DictsBySourceFile[file].Values.Sum(), Is.EqualTo(psms.Count),
                $"plot {plotType}, file {file}: binned count does not match the PSM count");
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

    /// <summary>
    /// Ambiguous notches are split across bins rather than parsed as one value or dropped. Every notch in
    /// the fixture is 0, so the first few rows are rewritten to "0|1": each must add one value to bin 1.
    /// </summary>
    [Test]
    public void AmbiguousNotchesAreSplitAcrossBins()
    {
        const int ambiguousRows = 3;
        string source = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "VariantCrossTest.psmtsv");
        string[] lines = File.ReadAllLines(source);
        int notchColumn = Array.IndexOf(lines[0].Split('\t'), SpectrumMatchFromTsvHeader.Notch);
        Assert.That(notchColumn, Is.GreaterThanOrEqualTo(0), "fixture has no notch column");

        string path = Path.Combine(TestContext.CurrentContext.TestDirectory, $"ambiguousNotch_{Guid.NewGuid():N}.psmtsv");
        File.WriteAllLines(path, lines.Select((line, i) =>
        {
            if (i == 0 || i > ambiguousRows) return line;
            string[] fields = line.Split('\t');
            Assert.That(fields[notchColumn], Is.EqualTo("0"), "fixture notch changed; this test assumes 0");
            fields[notchColumn] = "0|1";
            return string.Join('\t', fields);
        }));

        try
        {
            var bySourceFile = SpectrumMatchTsvReader.ReadTsv(path, out _)
                .GroupBy(p => p.FileNameWithoutExtension)
                .ToDictionary(g => g.Key, g => new ObservableCollection<SpectrumMatchFromTsv>(g));
            var data = SpectrumMatchHistograms.GetHistogramData(14, bySourceFile);

            int totalPsms = bySourceFile.Values.Sum(v => v.Count);
            int binned = data.DictsBySourceFile.Values.Sum(d => d.Values.Sum());
            int inBinOne = data.DictsBySourceFile.Values.Sum(d => d.GetValueOrDefault("1"));

            Assert.That(binned, Is.EqualTo(totalPsms + ambiguousRows));
            Assert.That(inBinOne, Is.EqualTo(ambiguousRows));
        }
        finally
        {
            File.Delete(path);
        }
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

    /// <summary>
    /// The exact semantics MetaDrawTest.TestPlotModelStat_Notch_PipeOnlyAndEmptyParts pins, asserted
    /// here too because that test is Windows-only. Dropping the DefaultIfEmpty("0") when this helper
    /// was moved made pipe-only notches vanish instead of binning at zero, and only Windows caught it.
    /// </summary>
    [TestCase("|", new[] { 0.0 })]
    [TestCase("||", new[] { 0.0 })]
    [TestCase("|1.003|", new[] { 1.003 })]
    [TestCase("", new[] { 0.0 })]
    [TestCase(null, new[] { 0.0 })]
    [TestCase("1|2", new[] { 1.0, 2.0 })]
    [TestCase(" 1 | 2 ", new[] { 1.0, 2.0 })]
    public void AmbiguousNotchesParseToExpectedValues(string notch, double[] expected)
    {
        Assert.That(SpectrumMatchHistograms.ParseAmbiguousNotchValues(notch).ToArray(),
            Is.EqualTo(expected).Within(1e-9));
    }
}
