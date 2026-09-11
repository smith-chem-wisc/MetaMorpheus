using System;
using System.IO;
using System.Linq;
using System.Reflection;
using EngineLayer;
using NUnit.Framework;
using Readers;
using TaskLayer;
using Mzml = IO.MzML.Mzml;

namespace Test
{
    /// <summary>
    /// Tests for the external MS1 feature ("FromFile") precursor source wired in by PR #2650:
    /// adjacent-file auto-discovery, and the additive precursor source + dedup in _GetMs2Scans.
    /// </summary>
    [TestFixture]
    public static class ExternalMs1FeatureTests
    {
        // fix_008 — TryFindAdjacentMs1FeatureFile auto-discovery. The method is private static, so it is
        // exercised via reflection (no production visibility change just for a test).
        [Test]
        public static void TryFindAdjacentMs1FeatureFile_FindsSiblingThenReturnsNullWhenAbsent()
        {
            string dir = Path.Combine(Path.GetTempPath(), "mm_ms1feature_" + Guid.NewGuid().ToString("N"));
            Directory.CreateDirectory(dir);
            try
            {
                string raw = Path.Combine(dir, "myRun.mzML");
                File.WriteAllText(raw, "");                                  // contents irrelevant to discovery
                string expectedFeature = Path.Combine(dir, "myRun_ms1.feature");

                MethodInfo finder = typeof(MetaMorpheusTask).GetMethod(
                    "TryFindAdjacentMs1FeatureFile", BindingFlags.NonPublic | BindingFlags.Static);
                Assert.That(finder, Is.Not.Null, "TryFindAdjacentMs1FeatureFile not found via reflection");

                // No sibling yet -> null
                Assert.That((string)finder.Invoke(null, new object[] { raw }), Is.Null);

                // Sibling present -> returns its path
                File.WriteAllText(expectedFeature, "");
                Assert.That((string)finder.Invoke(null, new object[] { raw }), Is.EqualTo(expectedFeature));
            }
            finally
            {
                Directory.Delete(dir, true);
            }
        }

        // fix_009 — additive FromFile precursor source and dedup through _GetMs2Scans (via the public
        // GetMs2Scans). The _ms1.feature fixture was generated from this same mzML's MS1 scans (mzLib
        // PR #1069 consensus pipeline) so its features align in RT/mass with the file's MS2 precursors.
        [Test]
        public static void GetMs2Scans_FromFileAdditiveSource_InjectsScoredPrecursorsAndDedups()
        {
            string dataDir = Path.Combine(TestContext.CurrentContext.TestDirectory, "TopDownTestData");
            string mzml = Path.Combine(dataDir, "TDGPTMDSearchSingleSpectra.mzML");
            string featureFile = Path.Combine(dataDir, "TDGPTMDSearchSingleSpectra_ms1.feature");
            Assume.That(File.Exists(mzml), $"missing {mzml}");
            Assume.That(File.Exists(featureFile), $"missing {featureFile}");

            var dataFile = Mzml.LoadAllStaticData(mzml);

            // Classic precursor decon only (charge capped to keep the test light; the fixture's features
            // top out near charge 9).
            var classicOnly = new CommonParameters(deconvolutionMaxAssumedChargeState: 20);

            // FromFile only: disable classic decon and scan-header so the additive source is the *only*
            // contributor, isolating its behavior.
            var fromFileParams = new FromFileDeconvolutionParameters(featureFile, 1, 20);
            var fromFileOnly = new CommonParameters(
                doPrecursorDeconvolution: false,
                useProvidedPrecursorInfo: false,
                deconvolutionMaxAssumedChargeState: 20,
                additionalPrecursorDeconParams: fromFileParams);

            // Both sources -> exercises the shared PrecursorSet dedup across sources.
            var combinedParams = new CommonParameters(
                deconvolutionMaxAssumedChargeState: 20,
                additionalPrecursorDeconParams: new FromFileDeconvolutionParameters(featureFile, 1, 20));

            var classicScans = MetaMorpheusTask.GetMs2Scans(dataFile, mzml, classicOnly).ToList();
            var fromFileScans = MetaMorpheusTask.GetMs2Scans(dataFile, mzml, fromFileOnly).ToList();
            var combinedScans = MetaMorpheusTask.GetMs2Scans(dataFile, mzml, combinedParams).ToList();

            // The additive source actually injects precursors.
            Assert.That(fromFileScans, Is.Not.Empty, "FromFile source produced no precursors");

            // fix_003: FromFile precursors carry the computed generic decon score, not the default 0.
            Assert.That(fromFileScans.Any(s => s.PrecursorDeconvolutionScore != 0),
                "FromFile precursors all have a default (0) DeconvolutionScore");

            // Additive: combining sources never loses the classic precursors.
            Assert.That(combinedScans.Count, Is.GreaterThanOrEqualTo(classicScans.Count));

            // Dedup: combining can only remove overlaps, never exceed the union of the two sources.
            Assert.That(combinedScans.Count, Is.LessThanOrEqualTo(classicScans.Count + fromFileScans.Count));
        }

        // fix_005 (FromFile resolution leg) — a configured, existing Ms1FeatureFilePath is resolved by
        // SetAllFileSpecificCommonParams into a FromFileDeconvolutionParameters additional source.
        [Test]
        public static void SetAllFileSpecificCommonParams_ResolvesMs1FeatureFilePathToFromFileSource()
        {
            string featureFile = Path.Combine(
                TestContext.CurrentContext.TestDirectory, "TopDownTestData", "TDGPTMDSearchSingleSpectra_ms1.feature");
            Assume.That(File.Exists(featureFile), $"missing {featureFile}");

            var common = new CommonParameters(deconvolutionMaxAssumedChargeState: 20);
            Assert.That(common.AdditionalPrecursorDeconvolutionParameters, Is.Null); // precondition

            var fsp = new FileSpecificParameters { Ms1FeatureFilePath = featureFile };
            CommonParameters resolved = MetaMorpheusTask.SetAllFileSpecificCommonParams(common, fsp);

            Assert.That(resolved.AdditionalPrecursorDeconvolutionParameters,
                Is.TypeOf<FromFileDeconvolutionParameters>());
        }
    }
}
