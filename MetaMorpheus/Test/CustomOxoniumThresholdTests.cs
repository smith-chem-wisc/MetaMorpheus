using EngineLayer;
using EngineLayer.GlycoSearch;
using Nett;
using NUnit.Framework;
using System.IO;
using TaskLayer;

namespace Test
{
    /// <summary>
    /// The user-set CustomOxoniumRelativeIntensityThreshold: it decides when a custom-monosaccharide diagnostic ion counts as
    /// observed, leaves the built-in NeuAc and HexNAc(1)Hex(1) rules at their fixed 0.02, and rejects text and values that are
    /// not a fraction in [0, 1]. Tests that register custom monosaccharides restore startup state in a finally block.
    /// </summary>
    [TestFixture]
    [NonParallelizable] // mutates process-wide Glycan monosaccharide statics
    public class CustomOxoniumThresholdTests
    {
        private const int Ion512 = 51219700; // 512.197 m/z * 1e5, not a built-in oxonium ion

        private static void RestoreStartupMonosaccharides()
        {
            Glycan.ResetCustomMonosaccharides();
            string shipped = GlobalVariables.CustomMonosaccharidePath;
            if (File.Exists(shipped))
            {
                GlycanDatabase.LoadCustomMonosaccharides(shipped);
            }
        }

        [Test]
        public void CheckOxoniumPresence_CountsOnlyIonsAboveTheGivenThreshold()
        {
            // reference = 100, so the ratios are 0.01, 0.05, 0.10 and 0.20.
            double[] intensities = { 1, 5, 10, 20 };
            Assert.Multiple(() =>
            {
                Assert.That(GlycoPeptides.CheckOxoniumPresence(intensities, 0, 100, 0.05), Is.False, "below");
                Assert.That(GlycoPeptides.CheckOxoniumPresence(intensities, 1, 100, 0.05), Is.False, "exactly at the threshold is not above it");
                Assert.That(GlycoPeptides.CheckOxoniumPresence(intensities, 2, 100, 0.05), Is.True, "above");
                Assert.That(GlycoPeptides.CheckOxoniumPresence(intensities, 2, 100, 0.15), Is.False, "raising the threshold drops the 0.10 ion");
                Assert.That(GlycoPeptides.CheckOxoniumPresence(intensities, 3, 100, 0.15), Is.True);
                Assert.That(GlycoPeptides.CheckOxoniumPresence(intensities, 0, 100, 0), Is.True, "threshold 0 counts any signal");
            });
        }

        [TestCase(0.02)]
        [TestCase(0.05)]
        [TestCase(0.3)]
        public void CustomOxoniumFilter_CandidateCarryingTheMono_IsCutWheneverItsIonIsAtOrBelowTheThreshold(double threshold)
        {
            // The candidate carries the custom monosaccharide, so the strict rule keeps it only when the ion is observed:
            // every intensity at or below threshold * 138.055-reference must cut it, every intensity above must keep it.
            try
            {
                Glycan.ResetCustomMonosaccharides();
                Glycan.RegisterCustomMonosaccharide("SugarU", 'U', 17603209, new[] { Ion512 });
                int customIndex = Glycan.NameCharDic["SugarU"].Item2;
                int customSlot = Glycan.AllOxoniumIons.Length;

                byte[] kind = new byte[Glycan.KindCapacity];
                kind[customIndex] = 1;

                const double reference = 100;
                Assert.Multiple(() =>
                {
                    for (int intensity = 0; intensity <= 50; intensity++)
                    {
                        double[] intensities = new double[Glycan.AllOxoniumIonsIncludingCustoms.Length];
                        intensities[OxoniumIonReservedIndices.R138] = reference;
                        intensities[customSlot] = intensity;

                        bool expectedKept = intensity / reference > threshold;
                        Assert.That(GlycoPeptides.CustomOxoniumFilter(intensities, kind, threshold), Is.EqualTo(expectedKept),
                            $"N-search path, intensity {intensity}, threshold {threshold}");
                        Assert.That(GlycoPeptides.DiagonsticFilter(intensities, new GlycanBox(kind), threshold), Is.EqualTo(expectedKept),
                            $"O / N+O path, intensity {intensity}, threshold {threshold}");
                    }
                });
            }
            finally
            {
                RestoreStartupMonosaccharides();
            }
        }

        [TestCase(0.0)]
        [TestCase(0.5)]
        [TestCase(1.0)]
        public void DiagonsticFilter_BuiltInRulesKeepTheirFixedTwoPercent_WhateverTheCustomThreshold(double customThreshold)
        {
            try
            {
                Glycan.ResetCustomMonosaccharides();
                byte[] noNeuAc = new byte[Glycan.KindCapacity];
                byte[] withNeuAc = new byte[Glycan.KindCapacity];
                withNeuAc[2] = 1;

                // NeuAc at 3% of 138.055: above the built-in 2%, so a box without NeuAc is rejected and one with it is kept.
                double[] neuAcAtThreePercent = new double[Glycan.AllOxoniumIons.Length];
                neuAcAtThreePercent[OxoniumIonReservedIndices.R138] = 100;
                neuAcAtThreePercent[OxoniumIonReservedIndices.NeuAc274] = 3;
                neuAcAtThreePercent[OxoniumIonReservedIndices.NeuAc292] = 3;

                // NeuAc at 1%: below the built-in 2%, so a box with NeuAc is rejected.
                double[] neuAcAtOnePercent = new double[Glycan.AllOxoniumIons.Length];
                neuAcAtOnePercent[OxoniumIonReservedIndices.R138] = 100;
                neuAcAtOnePercent[OxoniumIonReservedIndices.NeuAc274] = 1;
                neuAcAtOnePercent[OxoniumIonReservedIndices.NeuAc292] = 1;

                // HexNAc(1)Hex(1) 366 at 3%: a box with neither Hex nor HexNAc is rejected. DiagonsticFilter reaches the 366
                // rule only when the two NeuAc ratios fall on opposite sides of 2%, so one NeuAc ion is set above it.
                double[] hexHexNAcAtThreePercent = new double[Glycan.AllOxoniumIons.Length];
                hexHexNAcAtThreePercent[OxoniumIonReservedIndices.R138] = 100;
                hexHexNAcAtThreePercent[OxoniumIonReservedIndices.NeuAc274] = 3;
                hexHexNAcAtThreePercent[OxoniumIonReservedIndices.HexHexNAc366] = 3;

                Assert.Multiple(() =>
                {
                    Assert.That(GlycoPeptides.DiagonsticFilter(neuAcAtThreePercent, new GlycanBox(noNeuAc), customThreshold), Is.False);
                    Assert.That(GlycoPeptides.DiagonsticFilter(neuAcAtThreePercent, new GlycanBox(withNeuAc), customThreshold), Is.True);
                    Assert.That(GlycoPeptides.DiagonsticFilter(neuAcAtOnePercent, new GlycanBox(withNeuAc), customThreshold), Is.False);
                    Assert.That(GlycoPeptides.DiagonsticFilter(hexHexNAcAtThreePercent, new GlycanBox(noNeuAc), customThreshold), Is.False);
                });
            }
            finally
            {
                RestoreStartupMonosaccharides();
            }
        }

        [TestCase("0.02", 0.02)]
        [TestCase("0", 0.0)]
        [TestCase("1", 1.0)]
        [TestCase(" 0.5 ", 0.5)]
        [TestCase(".05", 0.05)]
        public void TryParseCustomOxoniumThreshold_AcceptsAFractionBetweenZeroAndOne(string text, double expected)
        {
            Assert.That(GlycoPeptides.TryParseCustomOxoniumRelativeIntensityThreshold(text, out double threshold, out string error), Is.True);
            Assert.That(threshold, Is.EqualTo(expected));
            Assert.That(error, Is.Null);
        }

        [TestCase("", "not a number")]
        [TestCase(null, "not a number")]
        [TestCase("abc", "not a number")]
        [TestCase("2%", "not a number")]
        [TestCase("0,02", "not a number")]
        [TestCase("0.0.2", "not a number")]
        [TestCase("1e", "not a number")]
        [TestCase("NaN", "not a number")]
        [TestCase("Infinity", "not a number")]
        [TestCase("-0.01", "cannot be negative")]
        [TestCase("1.01", "cannot be above 1")]
        [TestCase("2", "cannot be above 1")]
        public void TryParseCustomOxoniumThreshold_RejectsWithAMessageTheUserCanActOn(string text, string expectedMessage)
        {
            Assert.That(GlycoPeptides.TryParseCustomOxoniumRelativeIntensityThreshold(text, out _, out string error), Is.False);
            Assert.That(error, Does.Contain(expectedMessage));
        }

        [Test]
        public void GlycoSearchParameters_ThresholdDefaultsToTwoPercentAndRejectsOutOfRange()
        {
            var parameters = new GlycoSearchParameters();
            Assert.That(parameters.CustomOxoniumRelativeIntensityThreshold, Is.EqualTo(0.02));

            parameters.CustomOxoniumRelativeIntensityThreshold = 0;
            parameters.CustomOxoniumRelativeIntensityThreshold = 1;
            Assert.Multiple(() =>
            {
                Assert.Throws<MetaMorpheusException>(() => parameters.CustomOxoniumRelativeIntensityThreshold = -0.01);
                Assert.Throws<MetaMorpheusException>(() => parameters.CustomOxoniumRelativeIntensityThreshold = 1.01);
                Assert.Throws<MetaMorpheusException>(() => parameters.CustomOxoniumRelativeIntensityThreshold = double.NaN);
            });
            Assert.That(parameters.CustomOxoniumRelativeIntensityThreshold, Is.EqualTo(1), "a rejected value leaves the last valid one");
        }

        [Test]
        public void GlycoSearchTaskToml_OldFileGetsTheDefault_AndASetThresholdRoundTrips()
        {
            string oldToml = Path.Combine(TestContext.CurrentContext.TestDirectory, "GlycoTestData", "GlycoSnip.toml");
            Assert.That(File.ReadAllText(oldToml), Does.Not.Contain("CustomOxoniumRelativeIntensityThreshold"));
            var oldTask = Toml.ReadFile<GlycoSearchTask>(oldToml, MetaMorpheusTask.tomlConfig);
            Assert.That(oldTask._glycoSearchParameters.CustomOxoniumRelativeIntensityThreshold, Is.EqualTo(0.02));

            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, "CustomOxoniumThresholdTomlTest");
            Directory.CreateDirectory(outputFolder);
            try
            {
                string path = Path.Combine(outputFolder, "task.toml");
                oldTask._glycoSearchParameters.CustomOxoniumRelativeIntensityThreshold = 0.07;
                Toml.WriteFile(oldTask, path, MetaMorpheusTask.tomlConfig);

                var reloaded = Toml.ReadFile<GlycoSearchTask>(path, MetaMorpheusTask.tomlConfig);
                Assert.That(reloaded._glycoSearchParameters.CustomOxoniumRelativeIntensityThreshold, Is.EqualTo(0.07));
            }
            finally
            {
                Directory.Delete(outputFolder, true);
            }
        }
    }
}
