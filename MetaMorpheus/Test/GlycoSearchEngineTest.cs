using Chemistry;
using Easy.Common.Extensions;
using EngineLayer;
using EngineLayer.DatabaseLoading;
using EngineLayer.GlycoSearch;
using MassSpectrometry;
using MzLibUtil;
using Nett;
using NUnit.Framework;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using TaskLayer;

namespace Test
{
    [TestFixture]
    public class GlycoSearchEngineTest
    {
        [Test]
        public void CreateGsm_WithWideProductTolerance_ScanInfo_p_IsCappedToOne() 
        {
            // Use a very wide product mass tolerance to force p > 1 so ScanInfo_p is capped at 1.
            var commonParameters = new CommonParameters(productMassTolerance: new AbsoluteTolerance(1_000_000), dissociationType: DissociationType.HCD, trimMsMsPeaks: false);

            // ensure glycan DB paths used by GlycoSearchEngine ctor are registered (filenames must match ctor arguments)
            string oglycanPath = "OGlycan.gdb";
            string nglycanPath = "NGlycan_ForNoSearch.gdb";
            if (!GlobalVariables.OGlycanDatabasePaths.Contains(oglycanPath)) GlobalVariables.OGlycanDatabasePaths.Add(oglycanPath);
            if (!GlobalVariables.NGlycanDatabasePaths.Contains(nglycanPath)) GlobalVariables.NGlycanDatabasePaths.Add(nglycanPath);

            // Load the test MGF file and get MS2 scans
            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\2019_09_16_StcEmix_35trig_EThcD25_rep1_9906.mgf");
            var myFileManager = new MyFileManager(true);
            var msDataFile = myFileManager.LoadFile(spectraFile,commonParameters);
            var ms2ScanWithSpecificMass = MetaMorpheusTask.GetMs2Scans(msDataFile, spectraFile, commonParameters).ToArray();
            Assert.That(ms2ScanWithSpecificMass.Length, Is.GreaterThan(0), "No MS2 scans found in test MGF.");

            // Create a peptide instance to satisfy CreateGsm inputs.
            var protein = new Protein("AATVGSLAGQPLQER", "accession");
            var peptides = protein.Digest(new DigestionParams(minPeptideLength: 1), new List<Modification>(), new List<Modification>()).ToList();
            var peptide = peptides.First();
            var peptideIndex = new List<PeptideWithSetModifications> { peptide };
            var globalGsms = new List<GlycoSpectralMatch>[ms2ScanWithSpecificMass.Length];

            // Create a GlycoSearchEngine instance using the wide product mass tolerance.
            var engine = new GlycoSearchEngine(globalGsms, ms2ScanWithSpecificMass, peptideIndex, null, null, 0, commonParameters, null, oglycanPath, nglycanPath, glycoSearchType: GlycoSearchType.OGlycanSearch, 30, 3, false, null);

            // Create minimal inputs for CreateGsm.
            var route = new Route();
            route.AddPos(2, 2, true);
            var oxoniumIonIntensities = new double[Glycan.AllOxoniumIons.Length];
            var localizationGraphs = new List<LocalizationGraph>();

            // Access the private CreateGsm method via reflection.
            var createGsmMethod = typeof(GlycoSearchEngine).GetMethod("CreateGsm", BindingFlags.NonPublic | BindingFlags.Instance);
            Assert.That(createGsmMethod, Is.Not.Null, "Unable to find private CreateGsm method via reflection.");

            // Compute expected p and verify it exceeds 1 so capping behavior is exercised.
            var scan = ms2ScanWithSpecificMass.First();
            double spectrumSize = scan.TheScan.MassSpectrum.Size;
            double rangeWidth = scan.TheScan.MassSpectrum.Range.Width;
            double tolWidth = commonParameters.ProductMassTolerance.GetRange(1000).Width;
            double pEstimate = rangeWidth > 0 ? spectrumSize * tolWidth / rangeWidth : double.PositiveInfinity;
            var result = (GlycoSpectralMatch)createGsmMethod.Invoke(engine, new object[] { ms2ScanWithSpecificMass.First(), 0, 0, peptide, route, oxoniumIonIntensities, localizationGraphs });

            // Compute expected p and verify it exceeds 1 so capping behavior is exercised.
            Assert.That(pEstimate, Is.GreaterThan(1d), "Test must drive p > 1 to validate the cap.");
            Assert.That(result.ScanInfo_p, Is.EqualTo(1), "P value should be capped at 1 when product mass tolerance is very wide.");
        }

        [Test]
        public static void TestLowResToleranceConstruction_Default() 
        {
            var commonParameters = new CommonParameters();
            Assert.That(commonParameters.ProductMassTolerance_LowRes, Is.Not.Null, "ChildScanMassTolerance should be initialized by default.");

            // Default product tolerance is a AbsoluteTolerance of 0.35 Da; verify type and numeric width equivalence
            Assert.That(commonParameters.ProductMassTolerance_LowRes, Is.TypeOf<AbsoluteTolerance>());
            Assert.That(commonParameters.ProductMassTolerance_LowRes.GetRange(1000).Width, Is.EqualTo(new AbsoluteTolerance(0.35).GetRange(1000).Width));
        }

        [Test]
        public static void TestLowResToleranceConstruction_Default2() 
        {
            // In this toml settng, we omit the productMassTolerance_LowRes. The test confirms that the default value is 0.35 Da.
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, "TESTGlycoData");
            Directory.CreateDirectory(outputFolder);

            try
            {
                var inputTask = Toml.ReadFile<GlycoSearchTask>(
                    Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\GlycoSearchTaskconfigOGlycoTest_ToleranceTracking.toml"),
                    MetaMorpheusTask.tomlConfig);

                // Confirm read-in behavior (tracking expected when low-res omitted)
                Assert.That(inputTask.CommonParameters.ProductMassTolerance, Is.TypeOf<PpmTolerance>());
                Assert.That(inputTask.CommonParameters.ProductMassTolerance.Value, Is.EqualTo(40).Within(1e-9));
                Assert.That(inputTask.CommonParameters.ProductMassTolerance_LowRes, Is.TypeOf<AbsoluteTolerance>());
                Assert.That(inputTask.CommonParameters.ProductMassTolerance_LowRes.Value, Is.EqualTo(0.35));

                DbForTask db = new(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\P16150.fasta"), false);
                string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\2019_09_16_StcEmix_35trig_EThcD25_rep1_9906.mgf");

                new EverythingRunnerEngine(
                    new List<(string, MetaMorpheusTask)> { ("Task", inputTask) },
                    new List<string> { spectraFile },
                    new List<DbForTask> { db },
                    outputFolder).Run();

                // Confirm written TOML exists
                string writtenTaskToml = Path.Combine(outputFolder, "Task Settings", "Taskconfig.toml");
                Assert.That(File.Exists(writtenTaskToml), Is.True, "Expected task config TOML was not written.");

                // Confirm written TOML preserves tracked values
                var writtenTask = Toml.ReadFile<GlycoSearchTask>(writtenTaskToml, MetaMorpheusTask.tomlConfig);
                Assert.That(writtenTask.CommonParameters.ProductMassTolerance, Is.TypeOf<PpmTolerance>());
                Assert.That(writtenTask.CommonParameters.ProductMassTolerance.Value, Is.EqualTo(40).Within(1e-9));
                Assert.That(writtenTask.CommonParameters.ProductMassTolerance_LowRes, Is.TypeOf<AbsoluteTolerance>());
                Assert.That(writtenTask.CommonParameters.ProductMassTolerance_LowRes.Value, Is.EqualTo(0.35));
            }
            finally
            {
                if (Directory.Exists(outputFolder))
                {
                    Directory.Delete(outputFolder, true);
                }
            }
        }

        [Test]
        public static void TestLowResTolerance_Da() 
        {
            // Purpose:
            // - Validate that fragment matching use `ChildScanMassTolerance` for child scans (isChildScan flag).
            // - Confirm initial tolerances are set as expected: product, child, precursor.
            // - For matched fragment ions, verify the observed m/z (converted to mass) is within the correct tolerance.

            // --- initialize tolerances for test
            var lowResTolerance = new AbsoluteTolerance(0.5);
            var highResTolerance = new PpmTolerance(20);
            // Create CommonParameters with explicit tolerances
            var commonParameters = new CommonParameters(
                productMassTolerance: highResTolerance,
                productMassTolerance_LowRes: lowResTolerance,
                dissociationType: DissociationType.ETD,
                trimMsMsPeaks: false);
            Assert.That(commonParameters.ProductMassTolerance_LowRes.GetRange(1000).Width, Is.EqualTo(lowResTolerance.GetRange(1000).Width));
            Assert.That(commonParameters.ProductMassTolerance.GetRange(1000).Width, Is.EqualTo(highResTolerance.GetRange(1000).Width));

            // Sanity-check the parameter values were set
            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\2019_09_16_StcEmix_35trig_EThcD25_rep1_4565.mgf");
            var file = new MyFileManager(true).LoadFile(spectraFile, commonParameters);
            var scan = MetaMorpheusTask.GetMs2Scans(file, spectraFile, commonParameters).First();

            // Load a single MS2 scan from test data
            var peptide = new PeptideWithSetModifications("TTGSLEPSSGASGPQVSSVK");
            var glycoPeptide = new PeptideWithSetModifications("T[O-linked glycosylation:H1N1 on T]T[O-linked glycosylation:H1N1 on T]GSLEPSS[O-linked glycosylation:N1 on S]GASGPQVSSVK", GlobalVariables.AllModsKnownDictionary);
            var fragmentsForEachGlycoPeptide = GlycoPeptides.OGlyGetTheoreticalFragments(commonParameters.DissociationType, null, peptide, glycoPeptide);

            // Match fragments treating the scan as a child scan (isChildScan = true).
            // Expect matches to satisfy child scan tolerance (potentially wider).
            var childMatchedIons = MetaMorpheusEngine.MatchFragmentIons(scan, fragmentsForEachGlycoPeptide, commonParameters, isLowRes: true);
            Assert.That(childMatchedIons.Count, Is.GreaterThan(0));
            foreach (var ion in childMatchedIons)
            {
                double matchedMass = ion.Mz.ToMass(ion.Charge);
                double theoreticalMass = ion.NeutralTheoreticalProduct.NeutralMass;
                // Verify matched ion is within child-scan tolerance
                Assert.That(commonParameters.ProductMassTolerance_LowRes.Within(matchedMass, theoreticalMass));
                Assert.That(lowResTolerance.Within(matchedMass, theoreticalMass));
            }

            var childMatchedIonsHighRes = MetaMorpheusEngine.MatchFragmentIons(scan, fragmentsForEachGlycoPeptide, commonParameters, isLowRes: false);
            Assert.That(childMatchedIonsHighRes.Count, Is.GreaterThan(0));
            foreach (var ion in childMatchedIonsHighRes) 
            {
                double matchedMass = ion.Mz.ToMass(ion.Charge);
                double theoreticalMass = ion.NeutralTheoreticalProduct.NeutralMass;
                // Verify matched ion is within high-res tolerance
                Assert.That(commonParameters.ProductMassTolerance.Within(matchedMass, theoreticalMass));
                Assert.That(highResTolerance.Within(matchedMass, theoreticalMass));
            }
            // With a wider tolerance, we expect to match the same or more ions than with a narrower tolerance.
            Assert.That(childMatchedIons.Count >= childMatchedIonsHighRes.Count);
        }

        [Test]
        public static void TestLowResTolerance_ppm()
        {
            // Purpose:
            // - Validate that fragment matching use `ChildScanMassTolerance` for child scans (isChildScan flag).
            // - Confirm initial tolerances are set as expected: product, child, precursor.
            // - For matched fragment ions, verify the observed m/z (converted to mass) is within the correct tolerance.

            // --- initialize tolerances for test
            var lowResTolerance = new PpmTolerance(200);
            var highResTolerance = new PpmTolerance(20);

            // Create CommonParameters with explicit tolerances
            var commonParameters = new CommonParameters(
                productMassTolerance: highResTolerance,
                productMassTolerance_LowRes: lowResTolerance,
                dissociationType: DissociationType.ETD,
                trimMsMsPeaks: false);
            Assert.That(commonParameters.ProductMassTolerance_LowRes.GetRange(1000).Width, Is.EqualTo(lowResTolerance.GetRange(1000).Width));

            // Sanity-check the parameter values were set
            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\2019_09_16_StcEmix_35trig_EThcD25_rep1_4565.mgf");
            var file = new MyFileManager(true).LoadFile(spectraFile, commonParameters);
            var scan = MetaMorpheusTask.GetMs2Scans(file, spectraFile, commonParameters).First();

            // Load a single MS2 scan from test data
            var peptide = new PeptideWithSetModifications("TTGSLEPSSGASGPQVSSVK");
            var glycoPeptide = new PeptideWithSetModifications("T[O-linked glycosylation:H1N1 on T]T[O-linked glycosylation:H1N1 on T]GSLEPSS[O-linked glycosylation:N1 on S]GASGPQVSSVK", GlobalVariables.AllModsKnownDictionary);
            var fragmentsForEachGlycoPeptide = GlycoPeptides.OGlyGetTheoreticalFragments(commonParameters.DissociationType, null, peptide, glycoPeptide);

            // Match fragments treating the scan as a child scan (isChildScan = true).
            // Expect matches to satisfy child scan tolerance (potentially wider).
            var childMatchedIons = MetaMorpheusEngine.MatchFragmentIons(scan, fragmentsForEachGlycoPeptide, commonParameters, isLowRes: true);
            Assert.That(childMatchedIons.Count, Is.GreaterThan(0));
            foreach (var ion in childMatchedIons)
            {
                double matchedMass = ion.Mz.ToMass(ion.Charge);
                double theoreticalMass = ion.NeutralTheoreticalProduct.NeutralMass;
                // Verify matched ion is within child-scan tolerance
                Assert.That(commonParameters.ProductMassTolerance_LowRes.Within(matchedMass, theoreticalMass));
                Assert.That(lowResTolerance.Within(matchedMass, theoreticalMass));
            }

            // Match fragments with narrow tolerance, and the expected matches should be satisfied high-res tolerance.
            var childMatchedIonsHighRes = MetaMorpheusEngine.MatchFragmentIons(scan, fragmentsForEachGlycoPeptide, commonParameters, isLowRes: false);
            Assert.That(childMatchedIonsHighRes.Count, Is.GreaterThan(0));
            foreach (var ion in childMatchedIonsHighRes)
            {
                double matchedMass = ion.Mz.ToMass(ion.Charge);
                double theoreticalMass = ion.NeutralTheoreticalProduct.NeutralMass;
                // Verify matched ion is within child-scan tolerance
                Assert.That(commonParameters.ProductMassTolerance.Within(matchedMass, theoreticalMass));
                Assert.That(highResTolerance.Within(matchedMass, theoreticalMass));
            }

            // With a wider tolerance, we expect to match the same or more ions than with a narrower tolerance.
            Assert.That(childMatchedIons.Count >= childMatchedIonsHighRes.Count);
        }


        [Test]
        public static void LowRes_CalibrateAndSearch() 
        {
            // This test verifies that calibration-updated precursor/product tolerances are propagated to the subsequent GlycoSearch via the calibration-generated file-specific TOML, and confirms the GlycoSearch task config TOML is written with expected settings.
            string outputRoot = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestCalibThenGlyco");
            Directory.CreateDirectory(outputRoot);

            try 
            {
                // GlycoSearchEngine expects these DB names to exist in global path lists.
                string oglycanPath = "OGlycan.gdb";
                string nglycanPath = "NGlycan_ForNoSearch.gdb";
                if (!GlobalVariables.OGlycanDatabasePaths.Contains(oglycanPath))
                {
                    GlobalVariables.OGlycanDatabasePaths.Add(oglycanPath);
                }

                if (!GlobalVariables.NGlycanDatabasePaths.Contains(nglycanPath))
                {
                    GlobalVariables.NGlycanDatabasePaths.Add(nglycanPath);
                }

                // Use a calibratable mzML for calibration.
                string rawFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SmallCalibratible_Yeast.mzML");
                string database = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\smalldb.fasta");

                Tolerance originalPrecursorTolerance = new PpmTolerance(5);
                Tolerance originalProductTolerance = new PpmTolerance(20);

                // 1) Calibration task
                var calibrationTask = new CalibrationTask();
                calibrationTask.CommonParameters = new CommonParameters(
                    precursorMassTolerance: originalPrecursorTolerance,
                    productMassTolerance: originalProductTolerance);

                // 2) Glyco task (will consume calibrated file-specific toml generated by calibration)
                var glycoTask = new GlycoSearchTask();
                glycoTask.CommonParameters = new CommonParameters(
                    precursorMassTolerance: originalPrecursorTolerance,
                    productMassTolerance: originalProductTolerance);

                var runner = new EverythingRunnerEngine(
                    new List<(string, MetaMorpheusTask)>
                    {
                ("Calibration", calibrationTask),
                ("Glyco", glycoTask)
                    },
                    new List<string> { rawFile },
                    new List<DbForTask> { new DbForTask(database, false) },
                    outputRoot);

                runner.Run();

                // Calibration writes calibrated mzML and a file-specific toml next to it.
                string calibratedTomlPath = Path.Combine(outputRoot, "Calibration", "SmallCalibratible_Yeast-calib.toml");
                Assert.That(File.Exists(calibratedTomlPath), Is.True, "Calibration file-specific toml was not written.");

                // Parse calibrated tolerances from the calibration-generated file-specific toml.
                TomlTable calibratedTable = Toml.ReadFile(calibratedTomlPath, MetaMorpheusTask.tomlConfig);
                var calibratedFileSpecific = new FileSpecificParameters(calibratedTable);
                Assert.That(calibratedFileSpecific.PrecursorMassTolerance, Is.Not.Null);
                Assert.That(calibratedFileSpecific.ProductMassTolerance, Is.Not.Null);

                // Verify Glyco task consumed file-specific params from the calibrated file.
                Assert.That(glycoTask.FileSpecificParameters, Is.Not.Null);
                Assert.That(glycoTask.FileSpecificParameters.Count, Is.EqualTo(1));

                // Verify Glyco task consumed file-specific params from the calibrated file.
                Assert.That(glycoTask.FileSpecificParameters, Is.Not.Null);
                Assert.That(glycoTask.FileSpecificParameters.Count, Is.EqualTo(1));

                CommonParameters glycoFileSpecificParams = glycoTask.FileSpecificParameters[0].Parameters;
                Assert.That(glycoFileSpecificParams.PrecursorMassTolerance.GetType(), Is.EqualTo(calibratedFileSpecific.PrecursorMassTolerance.GetType()));
                Assert.That(glycoFileSpecificParams.ProductMassTolerance.GetType(), Is.EqualTo(calibratedFileSpecific.ProductMassTolerance.GetType()));
                Assert.That(glycoFileSpecificParams.PrecursorMassTolerance.Value, Is.EqualTo(calibratedFileSpecific.PrecursorMassTolerance.Value).Within(1e-9));
                Assert.That(glycoFileSpecificParams.ProductMassTolerance.Value, Is.EqualTo(calibratedFileSpecific.ProductMassTolerance.Value).Within(1e-9));

                // Low-res tolerance remains independent (default 0.35 Da unless explicitly set).
                Assert.That(glycoFileSpecificParams.ProductMassTolerance_LowRes, Is.TypeOf<AbsoluteTolerance>());
                Assert.That(glycoFileSpecificParams.ProductMassTolerance_LowRes.Value, Is.EqualTo(0.35).Within(1e-9));

                // Confirm Glyco task TOML is written.
                string glycoTaskTomlPath = Path.Combine(outputRoot, "Task Settings", "Glycoconfig.toml");
                Assert.That(File.Exists(glycoTaskTomlPath), Is.True, "Glyco task config toml was not written.");
            }
            finally
            {
                if (Directory.Exists(outputRoot))
                {
                    Directory.Delete(outputRoot, true);
                }
            }

        }

    }
}
