using EngineLayer;
using MassSpectrometry;
using System.Collections.Generic;
using NUnit.Framework;
using Omics;
using Omics.Fragmentation;
using ThermoFisher.CommonCore.Data.Business;
using System.Linq;
using MzLibUtil;
using Proteomics.ProteolyticDigestion;
using Omics.Modifications;
using Omics.Digestion;
using Proteomics;

namespace Test
{
    public class PepAnalysisEngineTests
    {
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

            var selectedPeptide = psm.BestMatchingBioPolymersWithSetMods.First().Peptide;

            // Act
            float result = EngineLayer.PepAnalysisEngine.Xcorr(psm, selectedPeptide);

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

            var selectedPeptide = psm.BestMatchingBioPolymersWithSetMods.First().Peptide;

            // Act
            float result = EngineLayer.PepAnalysisEngine.Xcorr(psm, selectedPeptide);

            // Assert
            Assert.That(0, Is.EqualTo(result));
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
            MzRange scanWindowRange = new MzRange(1,500);
            string scanFilter = "";
            MZAnalyzerType mzAnalyzer = MZAnalyzerType.Orbitrap;
            double totalIonCurrent = 1.0;
            double? injectionTime = 1.0;
            double[,] noiseData = new double[1,1];
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
