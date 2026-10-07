using System;
using System.IO;
using System.Linq;
using System.Reflection;
using EngineLayer;
using EngineLayer.DatabaseLoading;
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

            // Dedup across sources. The count bounds above hold even with no dedup at all, so check the
            // property itself: some FromFile precursors sit within the dedup tolerance of a classic one
            // in the same scan (otherwise there is nothing to collapse), and after combining, no scan
            // holds two precursors of the same charge within that tolerance.
            var tolerance = combinedParams.DeconvolutionMassTolerance;
            static bool Same(Ms2ScanWithSpecificMass a, Ms2ScanWithSpecificMass b, MzLibUtil.Tolerance t) =>
                a.TheScan.OneBasedScanNumber == b.TheScan.OneBasedScanNumber
                && a.PrecursorCharge == b.PrecursorCharge
                && t.Within(a.PrecursorMonoisotopicPeakMz, b.PrecursorMonoisotopicPeakMz);

            var classicByScan = classicScans.ToLookup(s => s.TheScan.OneBasedScanNumber);
            int overlapping = fromFileScans.Count(f => classicByScan[f.TheScan.OneBasedScanNumber].Any(c => Same(f, c, tolerance)));
            Assert.That(overlapping, Is.GreaterThan(0), "fixture has no FromFile precursor that duplicates a classic one");

            foreach (var scan in combinedScans.GroupBy(s => s.TheScan.OneBasedScanNumber))
            {
                var list = scan.ToList();
                for (int i = 0; i < list.Count; i++)
                    for (int j = i + 1; j < list.Count; j++)
                        Assert.That(Same(list[i], list[j], tolerance), Is.False,
                            $"scan {scan.Key}: z={list[i].PrecursorCharge} precursors at m/z {list[i].PrecursorMonoisotopicPeakMz} and {list[j].PrecursorMonoisotopicPeakMz} were not deduplicated");
            }
            Assert.That(combinedScans.Count, Is.LessThan(classicScans.Count + fromFileScans.Count),
                "overlapping precursors exist, so the combined set must be smaller than the sum");
        }

        // Negative mode: mzLib's FromFile source expands neutral masses as positive ions, so it would yield
        // nothing. It must be skipped with a warning rather than "enabled" and silently empty.
        [Test]
        public static void SetAllFileSpecificCommonParams_NegativeMode_SkipsFeatureFileWithWarning()
        {
            string featureFile = Path.Combine(
                TestContext.CurrentContext.TestDirectory, "TopDownTestData", "TDGPTMDSearchSingleSpectra_ms1.feature");
            Assume.That(File.Exists(featureFile), $"missing {featureFile}");

            var warnings = new System.Collections.Generic.List<string>();
            EventHandler<StringEventArgs> onWarn = (_, e) => warnings.Add(e.S);
            MetaMorpheusTask.WarnHandler += onWarn;
            try
            {
                var negative = new CommonParameters(deconvolutionMaxAssumedChargeState: -20);
                var fsp = new FileSpecificParameters { Ms1FeatureFilePath = featureFile };
                CommonParameters resolved = MetaMorpheusTask.SetAllFileSpecificCommonParams(negative, fsp);

                Assert.That(resolved.AdditionalPrecursorDeconvolutionParameters, Is.Null);
                Assert.That(warnings.Any(w => w.Contains(featureFile) && w.Contains("negative mode")), Is.True);

                // Positive mode still resolves it, so the skip is about polarity and nothing else.
                var positive = MetaMorpheusTask.SetAllFileSpecificCommonParams(
                    new CommonParameters(deconvolutionMaxAssumedChargeState: 20), fsp);
                Assert.That(positive.AdditionalPrecursorDeconvolutionParameters, Is.InstanceOf<FromFileDeconvolutionParameters>());
            }
            finally
            {
                MetaMorpheusTask.WarnHandler -= onWarn;
            }
        }

        [Test]
        public static void Ms1FeatureFileCalibrator_CombinesRoundsAndCorrectsAtTheApex()
        {
            // Two rounds over three MS1 scans at 10, 20 and 30 min.
            var round1 = new System.Collections.Generic.List<(double, double)> { (10, 2e-6), (20, 4e-6), (30, -1e-6) };
            var round2 = new System.Collections.Generic.List<(double, double)> { (10, 1e-6), (20, 0), (30, 0) };
            var (rts, factors) = Ms1FeatureFileCalibrator.Combine(new[] { round1, round2 });
            Assert.That(factors[0], Is.EqualTo((1 - 2e-6) * (1 - 1e-6)).Within(1e-15));
            Assert.That(factors[1], Is.EqualTo(1 - 4e-6).Within(1e-15));

            Assert.That(Ms1FeatureFileCalibrator.FactorAt(rts, factors, 14.9), Is.EqualTo(factors[0]));
            Assert.That(Ms1FeatureFileCalibrator.FactorAt(rts, factors, 15.1), Is.EqualTo(factors[1]));
            Assert.That(Ms1FeatureFileCalibrator.FactorAt(rts, factors, 99), Is.EqualTo(factors[2]));
            Assert.That(Ms1FeatureFileCalibrator.FactorAt(rts, factors, 0), Is.EqualTo(factors[0]));

            Assert.That(Ms1FeatureFileCalibrator.CanCalibrate("a/run_ms1.feature"), Is.True);
            Assert.That(Ms1FeatureFileCalibrator.CanCalibrate("a/run.feature.tsv"), Is.False);

            // Seconds-based file (Time_end > 500): apex 1210 s = 20.17 min must use the 20-min scan.
            string dir = Path.Combine(Path.GetTempPath(), "mm_ms1featurecal_" + Guid.NewGuid().ToString("N"));
            Directory.CreateDirectory(dir);
            try
            {
                string src = Path.Combine(dir, "run_ms1.feature");
                string dst = Path.Combine(dir, "run-calib_ms1.feature");
                new Ms1FeatureFile
                {
                    Results = new System.Collections.Generic.List<Ms1Feature>
                    {
                        new() { Id = 1, Mass = 10000.0, Intensity = 1e6, RetentionTimeBegin = 1190, RetentionTimeEnd = 1230,
                                RetentionTimeApex = 1210, IntensityApex = 1e5, ChargeStateMin = 5, ChargeStateMax = 9 },
                    }
                }.WriteResults(src);

                int n = Ms1FeatureFileCalibrator.WriteCalibratedCopy(src, dst, new[] { round1, round2 }, rtInSeconds: true);
                Assert.That(n, Is.EqualTo(1));
                var row = new Ms1FeatureFile(dst).Results.Single();
                Assert.That(row.Mass, Is.EqualTo(10000.0 * (1 - 4e-6)).Within(1e-9));
                Assert.That(row.RetentionTimeApex, Is.EqualTo(1210), "units are left as they were");
                Assert.That(row.ChargeStateMin, Is.EqualTo(5));
                Assert.That(row.ChargeStateMax, Is.EqualTo(9));
            }
            finally
            {
                Directory.Delete(dir, true);
            }
        }

        // The review finding on #2650: after calibration, the feature masses were never calibrated, yet the
        // path was carried into the -calib.toml. Build a feature file whose masses are exactly this file's
        // classic precursor masses, calibrate, and check that the calibrated feature masses track the
        // calibrated classic masses rather than keeping the pre-calibration offset.
        [Test]
        public static void CalibrationTask_CalibratesExternalFeatureMassesAndPointsTheCalibToml()
        {
            string testDir = TestContext.CurrentContext.TestDirectory;
            string unitTestFolder = Path.Combine(testDir, "CalibrationTask_ExternalFeatures");
            string outputFolder = Path.Combine(unitTestFolder, "TaskOutput");
            if (Directory.Exists(unitTestFolder)) Directory.Delete(unitTestFolder, true);
            Directory.CreateDirectory(outputFolder);
            try
            {
                string mzml = Path.Combine(unitTestFolder, "calfeat.mzML");
                File.Copy(Path.Combine(testDir, "TestData", "SmallCalibratible_Yeast.mzML"), mzml, true);
                string db = Path.Combine(testDir, "TestData", "smalldb.fasta");

                // One feature per classic precursor, keyed back to its scan and charge.
                var classicParams = new CommonParameters();
                var uncalibrated = MetaMorpheusTask.GetMs2Scans(Mzml.LoadAllStaticData(mzml), mzml, classicParams).ToList();
                var rows = new System.Collections.Generic.List<Ms1Feature>();
                var keys = new System.Collections.Generic.List<(int Scan, int Charge, double Mass)>();
                foreach (var s in uncalibrated)
                {
                    double rt = s.TheScan.RetentionTime;
                    rows.Add(new Ms1Feature
                    {
                        Id = rows.Count, Mass = s.PrecursorMass, Intensity = 1e6, IntensityApex = 1e5,
                        RetentionTimeBegin = rt - 0.05, RetentionTimeEnd = rt + 0.05, RetentionTimeApex = rt,
                        ChargeStateMin = s.PrecursorCharge, ChargeStateMax = s.PrecursorCharge,
                    });
                    keys.Add((s.TheScan.OneBasedScanNumber, s.PrecursorCharge, s.PrecursorMass));
                }
                string feature = Path.Combine(unitTestFolder, "calfeat_ms1.feature");
                new Ms1FeatureFile { Results = rows }.WriteResults(feature);

                new CalibrationTask().RunTask(outputFolder,
                    new System.Collections.Generic.List<DbForTask> { new DbForTask(db, false) },
                    new System.Collections.Generic.List<string> { mzml }, "test");

                string calibratedMzml = Path.Combine(outputFolder, "calfeat-calib.mzML");
                string calibratedFeature = Path.Combine(outputFolder, "calfeat-calib_ms1.feature");
                string calibratedToml = Path.Combine(outputFolder, "calfeat-calib.toml");
                Assert.That(File.Exists(calibratedMzml), "fixture did not calibrate");
                Assert.That(File.Exists(calibratedFeature), "no calibrated feature file was written next to the calibrated mzML");

                var tomlParams = new FileSpecificParameters(Nett.Toml.ReadFile(calibratedToml, MetaMorpheusTask.tomlConfig));
                Assert.That(Path.GetFullPath(tomlParams.Ms1FeatureFilePath), Is.EqualTo(Path.GetFullPath(calibratedFeature)),
                    "the -calib.toml must point at the calibrated features, not the uncalibrated source");

                // Compare each feature with the classic precursor it was made from, both after calibration.
                var calibratedRows = new Ms1FeatureFile(calibratedFeature).Results;
                Assert.That(calibratedRows.Count, Is.EqualTo(rows.Count));
                var calibratedClassic = MetaMorpheusTask.GetMs2Scans(Mzml.LoadAllStaticData(calibratedMzml), calibratedMzml, classicParams)
                    .ToLookup(s => (s.TheScan.OneBasedScanNumber, s.PrecursorCharge));

                var before = new System.Collections.Generic.List<double>();
                var after = new System.Collections.Generic.List<double>();
                for (int i = 0; i < keys.Count; i++)
                {
                    var (scan, z, original) = keys[i];
                    var match = calibratedClassic[(scan, z)]
                        .OrderBy(c => Math.Abs(c.PrecursorMass - original)).FirstOrDefault();
                    if (match == null || Math.Abs(match.PrecursorMass - original) / original * 1e6 > 50)
                        continue;
                    before.Add(Math.Abs(original - match.PrecursorMass) / match.PrecursorMass * 1e6);
                    after.Add(Math.Abs(calibratedRows[i].Mass - match.PrecursorMass) / match.PrecursorMass * 1e6);
                }

                Assert.That(after.Count, Is.GreaterThan(50), "too few feature/precursor pairs to compare");
                double medianBefore = before.OrderBy(x => x).ElementAt(before.Count / 2);
                double medianAfter = after.OrderBy(x => x).ElementAt(after.Count / 2);
                TestContext.WriteLine($"pairs {after.Count}; median |feature - calibrated classic|: uncalibrated {medianBefore:F3} ppm, calibrated {medianAfter:F3} ppm");
                Assert.That(medianBefore, Is.GreaterThan(0.3), "calibration did not shift this fixture enough to test anything");
                Assert.That(medianAfter, Is.LessThan(0.1 * medianBefore),
                    "calibrated feature masses should track the calibrated classic precursors");
            }
            finally
            {
                if (Directory.Exists(unitTestFolder)) Directory.Delete(unitTestFolder, true);
            }
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

        // mzLib reads the feature file lazily, on the first MS2 scan that queries it, which is inside
        // GetMs2Scans and outside every per-file guard. An unreadable feature file (e.g. a stale or
        // truncated _ms1.feature auto-discovered next to the raw file) must instead disable the external
        // source with a warning, exactly as a missing one does, and the search carry on without it.
        [Test]
        public static void SetAllFileSpecificCommonParams_UnreadableMs1FeatureFile_DisablesSourceInsteadOfAbortingSearch()
        {
            string mzml = Path.Combine(TestContext.CurrentContext.TestDirectory, "TopDownTestData", "TDGPTMDSearchSingleSpectra.mzML");
            string dir = Path.Combine(Path.GetTempPath(), "mm_ms1feature_" + Guid.NewGuid().ToString("N"));
            Directory.CreateDirectory(dir);
            string badFeatureFile = Path.Combine(dir, "myRun_ms1.feature");
            File.WriteAllText(badFeatureFile, "this is not\ta feature file\n\u0001\u0002\n");

            var warnings = new System.Collections.Generic.List<string>();
            EventHandler<StringEventArgs> onWarn = (_, e) => warnings.Add(e.S);
            MetaMorpheusTask.WarnHandler += onWarn;
            try
            {
                var common = new CommonParameters(deconvolutionMaxAssumedChargeState: 20);
                var fsp = new FileSpecificParameters { Ms1FeatureFilePath = badFeatureFile };

                CommonParameters resolved = null;
                Assert.DoesNotThrow(() => resolved = MetaMorpheusTask.SetAllFileSpecificCommonParams(common, fsp));
                Assert.That(resolved.AdditionalPrecursorDeconvolutionParameters, Is.Null);
                Assert.That(warnings.Any(w => w.Contains(badFeatureFile)), Is.True,
                    "the user should be told their feature file was not used");

                var dataFile = Mzml.LoadAllStaticData(mzml);
                Assert.DoesNotThrow(() => MetaMorpheusTask.GetMs2Scans(dataFile, mzml, resolved).ToList());
            }
            finally
            {
                MetaMorpheusTask.WarnHandler -= onWarn;
                Directory.Delete(dir, true);
            }
        }

        // End to end through RunTask: a <basename>_ms1.feature beside the raw file is picked up without
        // any file-specific toml and announced as a warning; an unreadable one is announced, disabled,
        // and the search still completes.
        [Test]
        [TestCase(true)]
        [TestCase(false)]
        public static void RunTask_AdjacentMs1FeatureFile_IsAutoDiscoveredAndAnnounced(bool readable)
        {
            string testDir = TestContext.CurrentContext.TestDirectory;
            string outputFolder = Path.Combine(testDir, "TestAdjacentMs1Feature_" + readable);
            string inputFolder = Path.Combine(outputFolder, "inputs");
            Directory.CreateDirectory(inputFolder);
            string mzml = Path.Combine(inputFolder, "TDGPTMDSearchSingleSpectra.mzML");
            string feature = Path.Combine(inputFolder, "TDGPTMDSearchSingleSpectra_ms1.feature");
            string fasta = Path.Combine(inputFolder, "ThreeHumanHistone.fasta");
            File.Copy(Path.Combine(testDir, "TopDownTestData", "TDGPTMDSearchSingleSpectra.mzML"), mzml, true);
            File.Copy(Path.Combine(testDir, "TopDownTestData", "ThreeHumanHistone.fasta"), fasta, true);
            if (readable)
                File.Copy(Path.Combine(testDir, "TopDownTestData", "TDGPTMDSearchSingleSpectra_ms1.feature"), feature, true);
            else
                File.WriteAllText(feature, "this is not\ta feature file\n");

            var warnings = new System.Collections.Generic.List<string>();
            EventHandler<StringEventArgs> onWarn = (_, e) => warnings.Add(e.S);
            MetaMorpheusTask.WarnHandler += onWarn;
            try
            {
                var searchTask = new SearchTask();
                Assert.DoesNotThrow(() => searchTask.RunTask(outputFolder,
                    new System.Collections.Generic.List<DbForTask> { new DbForTask(fasta, false) },
                    new System.Collections.Generic.List<string> { mzml }, "normal"));

                Assert.That(warnings.Any(w => w.Contains("Found adjacent MS1 feature file") && w.Contains(feature)), Is.True);
                Assert.That(warnings.Any(w => w.Contains("Could not read Ms1FeatureFilePath")), Is.EqualTo(!readable));
                Assert.That(File.Exists(Path.Combine(outputFolder, "AllPSMs.psmtsv")), Is.True);
            }
            finally
            {
                MetaMorpheusTask.WarnHandler -= onWarn;
                Directory.Delete(outputFolder, true);
            }
        }
    }
}
