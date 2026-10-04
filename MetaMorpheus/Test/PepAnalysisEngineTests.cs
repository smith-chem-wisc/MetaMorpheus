using EngineLayer;
using MassSpectrometry;
using System.Collections.Generic;
using NUnit.Framework;
using Omics.Fragmentation;
using System.Linq;
using MzLibUtil;
using Proteomics.ProteolyticDigestion;
using Omics.Modifications;
using Omics.Digestion;
using Proteomics;
using System;
using Chemistry;
using Omics;

namespace Test
{
    public class PepAnalysisEngineTests
    {
        [TestCase(5, 2.07918119f)]
        [Test]
        public static void GetLog10Factorial_ReturnsCorrectValue(int n, float? expected)
        {
            // Act
            float? result = PepAnalysisEngine.GetLog10Factorial(n);

            // Assert
            Assert.That(expected, Is.EqualTo(result));
        }

        [Test]
        public static void TestGetFraggerHyperScore()
        {
            MassDiffAcceptor searchModes = new DotMassDiffAcceptor(null, new List<double> { 0, 1.0029 }, new PpmTolerance(5));

            var p = new Protein("PEPTIDE", "accession");
            var d = p.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>()).ToList();
            PeptideWithSetModifications pep = d.First();

            CommonParameters commonParameters = new CommonParameters();

            var digested = p.Digest(commonParameters.DigestionParams, new List<Modification>(), new List<Modification>()).ToList();

            TestDataFile t = new TestDataFile(new List<IBioPolymerWithSetMods> { pep });

            MsDataScan mzLibScan1 = t.GetOneBasedScan(2);
            Ms2ScanWithSpecificMass scan1 = new Ms2ScanWithSpecificMass(mzLibScan1, pep.MonoisotopicMass.ToMz(1), 1, null, new CommonParameters());

            var peptideFragmentIons = new List<MatchedFragmentIon>
            {
                new MatchedFragmentIon(new Product(ProductType.b, FragmentationTerminus.N, 100, 1, 1, 0), 100, 100, 1),
                new MatchedFragmentIon(new Product(ProductType.b, FragmentationTerminus.N, 200, 2, 2, 0), 200, 200, 2),
                new MatchedFragmentIon(new Product(ProductType.b, FragmentationTerminus.N, 300, 3, 3, 0), 300, 300, 3),
                new MatchedFragmentIon(new Product(ProductType.y, FragmentationTerminus.C, 100, 1, 1, 0), 100, 100, 1),
                new MatchedFragmentIon(new Product(ProductType.y, FragmentationTerminus.C, 200, 2, 2, 0), 200, 200, 1),
                new MatchedFragmentIon(new Product(ProductType.y, FragmentationTerminus.C, 300, 3, 3, 0), 300, 300, 1)
            };

            SpectralMatch psm1 = new PeptideSpectralMatch(pep, 0, 3, 0, scan1, commonParameters, peptideFragmentIons);

            psm1.ResolveAllAmbiguities();

            // Act

            float hyperScore = PepAnalysisEngine.GetFraggerHyperScore(psm1, psm1.BestMatchingBioPolymersWithSetMods.First().SpecificBioPolymer);


            // Assert. Kong et al. 2017 define hyperscore = log(Nb! * Ny! * sum(I_b) * sum(I_y)), the product of
            // the per-terminus sums: log10(3!) + log10(3!) + log10(600 * 600). The single-sum form would give 4.635484.
            Assert.That(7.112605f, Is.EqualTo(hyperScore).Within(0.000001f));

            // A one-sided match scores on the matched side only: log10(3!) + log10(600), with no arbitrary
            // constant standing in for the missing y-ion intensity term.
            var bOnly = peptideFragmentIons.Where(i => i.NeutralTheoreticalProduct.Terminus == FragmentationTerminus.N).ToList();
            SpectralMatch bOnlyPsm = new PeptideSpectralMatch(pep, 0, 3, 0, scan1, commonParameters, bOnly);
            bOnlyPsm.ResolveAllAmbiguities();
            float bOnlyScore = PepAnalysisEngine.GetFraggerHyperScore(bOnlyPsm, bOnlyPsm.BestMatchingBioPolymersWithSetMods.First().SpecificBioPolymer);
            Assert.That(bOnlyScore, Is.EqualTo((float)(Math.Log10(6) + Math.Log10(600))).Within(0.000001f));
        }

        [Test]
        [TestCase(1, 0.0)]
        [TestCase(2, 0.30103)]
        [TestCase(3, 0.7781513)]
        [TestCase(4, 1.38021123)]
        [TestCase(5, 2.07918119)]
        [TestCase(6, 2.85733247)]
        [TestCase(7, 3.70243049)]
        [TestCase(8, 4.60552073)]
        [TestCase(9, 5.559763)]
        [TestCase(10, 6.559763)]
        public static void GetLog10Factorial_ReturnsExpectedResult(int n, double? expected)
        {
            float? result = PepAnalysisEngine.GetLog10Factorial(n);
            Assert.That((float)expected, Is.EqualTo(result).Within(0.0001));
        }

        [Test]
        public static void GetLog10Factorial_NegativeInputReturnsNull()
        {
            Assert.That(PepAnalysisEngine.GetLog10Factorial(-1), Is.Null);
        }

        [Test]
        public static void GetLog10Factorial_LargeInput_ReturnsExpectedResult()
        {
            int n = 20;
            float? result = PepAnalysisEngine.GetLog10Factorial(n);
            double? expected = 0.0;
            for (int i = 1; i <= n; i++)
            {
                expected += Math.Log10(i);
            }
            Assert.That((float)expected, Is.EqualTo(result).Within(0.0001)); // Allowing a small tolerance for floating-point comparison
        }

        [Test]
        public void Xcorr_ValidInput_ReturnsExpectedResult()
        {
            // Arrange
            var xArray = new double[] { 100, 150, 200, 250, 300 };
            var yArray = new double[] { 10, 20, 30, 40, 50 };

            var fragments = new List<MatchedFragmentIon>
            {
                new MatchedFragmentIon(new Product(ProductType.b, Omics.Fragmentation.FragmentationTerminus.N, 150, 1,1,0), 150, 20, 1),
                new MatchedFragmentIon(new Product(ProductType.b, Omics.Fragmentation.FragmentationTerminus.N, 250, 1,1,0), 250, 40, 1),
            };

            var psm = CreateSpectralMatch(xArray, yArray, [150, 250], [20, 40], fragments);

            var selectedPeptide = psm.BestMatchingBioPolymersWithSetMods.First().SpecificBioPolymer;

            // Act
            float result = EngineLayer.PepAnalysisEngine.Xcorr(psm, selectedPeptide, new MzSpectrum(xArray, yArray, true));

            // Assert
            Assert.That(58.8, Is.EqualTo(result).Within(1)); // Allowing a small tolerance for floating-point comparison
        }

        [Test]
        public void Xcorr_EmptyFragments_ReturnsZero()
        {
            // Arrange
            var xArray = new double[] { 100, 150, 200, 250, 300 };
            var yArray = new double[] { 10, 20, 30, 40, 50 };

            var fragments = new List<MatchedFragmentIon>
            {
            };

            var psm = CreateSpectralMatch(xArray, yArray, new double[0], new double[0], fragments);

            var selectedPeptide = psm.BestMatchingBioPolymersWithSetMods.First().SpecificBioPolymer;

            // Act
            float result = EngineLayer.PepAnalysisEngine.Xcorr(psm, selectedPeptide, new MzSpectrum(xArray, yArray, true));

            // Assert
            Assert.That(0, Is.EqualTo(result));
        }

        [Test]
        public void Xcorr_IonOutsideTheSpectrumWindow_DoesNotThrow()
        {
            // No peak lies within +/-75 of this ion, so there is no background window to read.
            var xArray = new double[] { 100, 150, 200 };
            var yArray = new double[] { 10, 20, 30 };
            var fragments = new List<MatchedFragmentIon>
            {
                new MatchedFragmentIon(new Product(ProductType.b, FragmentationTerminus.N, 500, 1, 1, 0), 500, 40, 1),
            };
            var psm = CreateSpectralMatch(xArray, yArray, [500], [40], fragments);
            var selectedPeptide = psm.BestMatchingBioPolymersWithSetMods.First().SpecificBioPolymer;

            float result = 0;
            Assert.DoesNotThrow(() => result = PepAnalysisEngine.Xcorr(psm, selectedPeptide, new MzSpectrum(xArray, yArray, true)));
            Assert.That(result, Is.EqualTo(40).Within(1e-6));
        }

        private SpectralMatch CreateSpectralMatch(double[] xArray, double[] yArray, double[] fragmentMz, double[] fragmentIntensity, List<MatchedFragmentIon> matchedFragmentIons)
        {
            PeptideWithSetModifications pwsm = new PeptideWithSetModifications(new Protein("PEPTIDE", "ACCESSION", "ORGANISM"), new DigestionParams(), 1, 2, CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0);
            int notch = 0;
            double score = 0;
            int scanIndex = 1;
            Ms2ScanWithSpecificMass scan = CreateMs2ScanWithSpecificMass(xArray, yArray);
            CommonParameters commonParameters = new CommonParameters();

            return new PeptideSpectralMatch(pwsm, notch, score, scanIndex, scan, commonParameters, matchedFragmentIons);
        }

        private Ms2ScanWithSpecificMass CreateMs2ScanWithSpecificMass(double[] xArray, double[] yArray)
        {
            MsDataScan scan = CreateMsDataScan(xArray, yArray);
            double precursorMonoisotopicPeakMz = 1;
            int precursorCharge = 1;
            string fullFilePath = "";
            CommonParameters commonParam = new CommonParameters();

            return new Ms2ScanWithSpecificMass(scan, precursorMonoisotopicPeakMz, precursorCharge, fullFilePath, commonParam);
        }
        private MsDataScan CreateMsDataScan(double[] xArray, double[] yArray)
        {
            MzSpectrum massSpectrum = CreateMzSpectrum(xArray, yArray);
            int oneBasedScanNumber = 1;
            int msnOrder = 1;
            bool isCentroid = true;
            MassSpectrometry.Polarity polarity = MassSpectrometry.Polarity.Positive;
            double retentionTime = 1.0;
            MzRange scanWindowRange = new MzRange(1, 500);
            string scanFilter = "";
            MZAnalyzerType mzAnalyzer = MZAnalyzerType.Orbitrap;
            double totalIonCurrent = 1.0;
            double? injectionTime = 1.0;
            double[,] noiseData = new double[1, 1];
            string nativeId = "";

            return new MsDataScan(massSpectrum, oneBasedScanNumber, msnOrder, isCentroid, polarity, retentionTime, scanWindowRange, scanFilter, mzAnalyzer, totalIonCurrent, injectionTime, noiseData, nativeId);

        }
        private MzSpectrum CreateMzSpectrum(double[] xArray, double[] yArray)
        {
            double[] mz = xArray;
            double[] intensities = yArray;
            bool shouldCopy = true;
            return new MzSpectrum(mz, intensities, shouldCopy);
        }
    }
}
