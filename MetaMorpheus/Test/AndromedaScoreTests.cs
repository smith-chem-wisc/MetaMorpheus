using Chemistry;
using EngineLayer;
using EngineLayer.DatabaseLoading;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Omics.Fragmentation;
using Readers;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;

namespace Test
{
    [TestFixture]
    public static class AndromedaScoreTests
    {
        // Hand-computed from Cox et al. 2011: S = -10 log10 sum_{j=k..n} C(n,j) p^j (1-p)^(n-j), p = q/100.
        // n=4, k=2, q=10: P = 6(.01)(.81) + 4(.001)(.9) + .0001 = .0523, S = -10 log10(.0523) = 12.81498...
        [Test]
        [TestCase(4, 2, 10, 12.814983111327257)]
        // n=k=3, q=10: P = .1^3, S = 30 exactly
        [TestCase(3, 3, 10, 30.0)]
        // n=20, k=10, q=10: the value the review computed for the published form (the old code gave 104.1)
        [TestCase(20, 10, 10, 51.45639050989115)]
        // nothing matched: P(X >= 0) = 1
        [TestCase(20, 0, 10, 0.0)]
        public static void BinomialScoreMatchesHandComputedValue(int n, int k, int q, double expected)
        {
            Assert.That(AndromedaScoring.BinomialScore(n, k, q), Is.EqualTo(expected).Within(1e-9));
        }

        [Test]
        public static void BinomialScoreAgreesWithMathNetBinomialTail()
        {
            foreach (int n in new[] { 1, 5, 12, 30, 60 })
                for (int q = 1; q <= 12; q++)
                    for (int k = 1; k <= n; k++)
                    {
                        double tail = 1 - MathNet.Numerics.Distributions.Binomial.CDF(q / 100.0, n, k - 1);
                        if (tail < 1e-8) continue; // 1 - CDF loses precision there; the log-space sum does not
                        Assert.That(AndromedaScoring.BinomialScore(n, k, q), Is.EqualTo(-10 * Math.Log10(tail)).Within(1e-6), $"n={n} k={k} q={q}");
                    }
        }

        [Test]
        public static void BinomialScoreStaysFiniteForLongPeptides()
        {
            // C(300,150) overflows a double; the tail must still come out finite and ordered
            double all = AndromedaScoring.BinomialScore(300, 300, 10);
            Assert.That(all, Is.EqualTo(3000).Within(1e-6));
            Assert.That(AndromedaScoring.BinomialScore(300, 150, 10), Is.GreaterThan(0).And.LessThan(all));
        }

        [Test]
        public static void BinomialScoreRejectsImpossibleInputs()
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => AndromedaScoring.BinomialScore(4, 5, 10));
            Assert.Throws<ArgumentOutOfRangeException>(() => AndromedaScoring.BinomialScore(4, 2, 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => AndromedaScoring.BinomialScore(4, 2, 100));
        }

        /// <summary>
        /// Four theoretical fragments against a spectrum where one sits at intensity rank 0 of the
        /// 100-200 Th window, one at rank 2 of that window, one at rank 0 of 300-400 Th, and one is
        /// absent. Matches by depth: q=1,2 -> k=2; q>=3 -> k=3. By hand:
        ///   q=1 32.2766, q=2 26.3144, q=3 39.7646, q=4 36.0499 ... q=10 24.3180
        /// so the score is the q=3 value, which a single fixed q (or q=1) would miss.
        /// </summary>
        [Test]
        public static void SpectrumScoreKeepsTopQPerWindowAndMaximisesOverQ()
        {
            double[] fragmentMz = { 150.0, 175.0, 350.0, 520.0 };
            var products = fragmentMz
                .Select((mz, i) => new Product(ProductType.b, FragmentationTerminus.N, mz.ToMass(1), i + 1, i + 1, 0))
                .ToList();

            var peaks = new List<(double mz, double intensity)>
            {
                (120.0, 50), (150.0, 100), (160.0, 80), (175.0, 60), (190.0, 10), // 100-200: 150 is rank 0, 175 is rank 2
                (350.0, 30), (380.0, 5),                                           // 300-400: 350 is rank 0
                (510.0, 40), (530.0, 40),                                          // 500-600: nothing at 520
            };
            var spectrum = new MzSpectrum(peaks.Select(p => p.mz).ToArray(), peaks.Select(p => p.intensity).ToArray(), false);

            double score = AndromedaScoring.Score(spectrum, products, new PpmTolerance(10));

            Assert.That(score, Is.EqualTo(AndromedaScoring.BinomialScore(4, 3, 3)).Within(1e-12));
            Assert.That(score, Is.EqualTo(39.7646).Within(1e-4));
        }

        /// <summary>
        /// The old code used every peak's density as p, so a dense spectrum gave p > 1 and a negative
        /// score. Top-q filtering keeps p &lt;= q/100 whatever the peak count.
        /// </summary>
        [Test]
        public static void DenseSpectrumNeverScoresNegative()
        {
            var mz = Enumerable.Range(0, 2000).Select(i => 300 + i * 0.75).ToArray();
            var intensity = mz.Select((_, i) => (double)(i % 97 + 1)).ToArray();
            var spectrum = new MzSpectrum(mz, intensity, false);
            var products = Enumerable.Range(0, 20)
                .Select(i => new Product(ProductType.y, FragmentationTerminus.C, mz[i * 90].ToMass(1), i + 1, i + 1, 0))
                .ToList();

            Assert.That(AndromedaScoring.Score(spectrum, products, new PpmTolerance(10)), Is.GreaterThanOrEqualTo(0));
        }

        [Test]
        [TestCase(SearchType.Classic)]
        [TestCase(SearchType.Modern)]
        public static void SearchWritesAndromedaScoreForTargetsAndDecoysWhenAskedOnly(SearchType searchType)
        {
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestAndromedaScore_" + searchType);
            string mzml = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\TaGe_SA_HeLa_04_subset_longestSeq.mzML");
            string fasta = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\hela_snip_for_unitTest.fasta");

            try
            {
                foreach (bool on in new[] { false, true })
                {
                    string folder = Path.Combine(outputFolder, on ? "on" : "off");
                    Directory.CreateDirectory(folder);
                    var task = new SearchTask
                    {
                        SearchParameters = new SearchParameters
                        {
                            SearchType = searchType,
                            WriteAndromedaScore = on,
                            WriteDecoys = true,
                            WriteHighQValuePsms = true,
                            DoLabelFreeQuantification = false,
                        }
                    };
                    task.RunTask(folder, new List<DbForTask> { new DbForTask(fasta, false) }, new List<string> { mzml }, "andromeda");

                    string[] lines = File.ReadAllLines(Path.Combine(folder, "AllPSMs.psmtsv"));
                    string[] header = lines[0].Split('\t');
                    int column = Array.IndexOf(header, PsmTsvWriter.AndromedaScoreHeader);
                    if (!on)
                    {
                        Assert.That(column, Is.EqualTo(-1), "default output must not change");
                        continue;
                    }

                    Assert.That(column, Is.GreaterThan(-1));
                    int dct = Array.IndexOf(header, SpectrumMatchFromTsvHeader.DecoyContaminantTarget);
                    var rows = lines.Skip(1).Select(l => l.Split('\t')).ToList();
                    Assert.That(rows, Is.Not.Empty);
                    foreach (var row in rows)
                    {
                        Assert.That(double.TryParse(row[column], System.Globalization.NumberStyles.Float, System.Globalization.CultureInfo.InvariantCulture, out double s), Is.True, row[column]);
                        Assert.That(s, Is.GreaterThanOrEqualTo(0));
                    }
                    var decoys = rows.Where(r => r[dct].Contains('D')).ToList();
                    var targets = rows.Where(r => r[dct] == "T").ToList();
                    Assert.That(decoys, Is.Not.Empty, "fixture should produce decoy PSMs");
                    // decoys are scored on the same scale as targets, not given a 0 placeholder
                    Assert.That(decoys.Any(r => double.Parse(r[column], System.Globalization.CultureInfo.InvariantCulture) > 0), Is.True);
                    Assert.That(targets.Any(r => double.Parse(r[column], System.Globalization.CultureInfo.InvariantCulture) > 0), Is.True);
                }
            }
            finally
            {
                if (Directory.Exists(outputFolder))
                    Directory.Delete(outputFolder, true);
            }
        }
    }
}
