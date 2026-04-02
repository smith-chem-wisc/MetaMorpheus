using Chemistry;
using EngineLayer;
using EngineLayer.GlycoSearch;
using MassSpectrometry;
using MzLibUtil;
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
        public static void TestChildScanToleranceConstruction_Default() 
        {
            var commonParameters = new CommonParameters();
            Assert.That(commonParameters.ProductMassTolerance_LowRes, Is.Not.Null, "ChildScanMassTolerance should be initialized by default.");

            // Default product tolerance is a PpmTolerance of 20 ppm; verify type and numeric width equivalence
            Assert.That(commonParameters.ProductMassTolerance, Is.TypeOf<PpmTolerance>());
            Assert.That(commonParameters.ProductMassTolerance.GetRange(1000).Width, Is.EqualTo(new PpmTolerance(20).GetRange(1000).Width));

            // By default child tolerance should equal product tolerance
            Assert.That(commonParameters.ProductMassTolerance_LowRes.GetRange(1000).Width, Is.EqualTo(commonParameters.ProductMassTolerance.GetRange(1000).Width));
        }

        [Test]
        public static void TestChildScanTolerance_Da() 
        {
            // Purpose:
            // - Validate that fragment matching use `ChildScanMassTolerance` for child scans (isChildScan flag).
            // - Confirm initial tolerances are set as expected: product, child, precursor.
            // - For matched fragment ions, verify the observed m/z (converted to mass) is within the correct tolerance.

            // --- initialize tolerances for test
            var childScanTolerance = new AbsoluteTolerance(0.5);

            // Create CommonParameters with explicit tolerances
            var commonParameters = new CommonParameters(
                productMassTolerance_LowRes: childScanTolerance,
                dissociationType: DissociationType.ETD,
                trimMsMsPeaks: false);
            Assert.That(commonParameters.ProductMassTolerance_LowRes.GetRange(1000).Width, Is.EqualTo(childScanTolerance.GetRange(1000).Width));

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
                Assert.That(childScanTolerance.Within(matchedMass, theoreticalMass));
            }
        }

        [Test]
        public static void TestChildScanTolerance_ppm()
        {
            // Purpose:
            // - Validate that fragment matching use `ChildScanMassTolerance` for child scans (isChildScan flag).
            // - Confirm initial tolerances are set as expected: product, child, precursor.
            // - For matched fragment ions, verify the observed m/z (converted to mass) is within the correct tolerance.

            // --- initialize tolerances for test
            var childScanTolerance = new PpmTolerance(200);

            // Create CommonParameters with explicit tolerances
            var commonParameters = new CommonParameters(
                productMassTolerance_LowRes: childScanTolerance,
                dissociationType: DissociationType.ETD,
                trimMsMsPeaks: false);
            Assert.That(commonParameters.ProductMassTolerance_LowRes.GetRange(1000).Width, Is.EqualTo(childScanTolerance.GetRange(1000).Width));

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
                Assert.That(childScanTolerance.Within(matchedMass, theoreticalMass));
            }
        }
    }
}
