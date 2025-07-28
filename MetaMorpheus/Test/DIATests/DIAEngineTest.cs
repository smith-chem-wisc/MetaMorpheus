using EngineLayer.DIA;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using EngineLayer;
using TaskLayer;
using System.IO;
using Chemistry;

namespace Test.DIATests
{
    public class DIAEngineTest
    {
        public static MsDataScan[] FakeMs1Scans { get; set; }
        public static MsDataScan[] FakeMs2Scans { get; set; }
        public static IsotopicDistribution PrecursorDist { get; set; }
        public static IsotopicDistribution FragDist { get; set; }
        public static double Intensity { get; set; }

        [OneTimeSetUp]
        public void OneTimeSetup()
        {
            string peptide = "PEPTIDE";
            Intensity = 1e6;
            ChemicalFormula cf = new Proteomics.AminoAcidPolymer.Peptide(peptide).GetChemicalFormula();
            PrecursorDist = IsotopicDistribution.GetDistribution(cf, 0.125, 1e-8);

            string fakeFrag1 = "PEP";
            ChemicalFormula cf2 = new Proteomics.AminoAcidPolymer.Peptide(fakeFrag1).GetChemicalFormula();
            FragDist = IsotopicDistribution.GetDistribution(cf2, 0.125, 0.0001);

            //Create MS1 scans with two precursors that do not have any overlap in their retention time; one of the precursors has an apex RT at the third scan, the other at the eighth scan
            FakeMs1Scans = new MsDataScan[10];
            double[] intensityMultipliers = { 1, 2, 3, 2, 1, 1, 2, 3, 2, 1 };
            var massArray1 = PrecursorDist.Masses.SelectMany(v => new List<double> { v.ToMz(1) }).ToArray();
            var massArray2 = PrecursorDist.Masses.SelectMany(v => new List<double> { (v + 1).ToMz(1) }).ToArray();
            double[][] masses = Enumerable.Repeat(massArray1, 5).Concat(Enumerable.Repeat(massArray2, 5)).ToArray();

            for (int s = 0; s < FakeMs1Scans.Length; s++)
            {
                double[] intensities = PrecursorDist.Intensities.SelectMany(v => new List<double> { v * Intensity * intensityMultipliers[s] }).ToArray();
                FakeMs1Scans[s] = new MsDataScan(massSpectrum: new MzSpectrum(masses[s], intensities, false), oneBasedScanNumber: s + 1, msnOrder: 1, isCentroid: true,
                    polarity: Polarity.Positive, retentionTime: 1.0 + s / 10.0, scanWindowRange: new MzRange(400, 1600), scanFilter: "f",
                    mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: intensities.Sum(), injectionTime: 1.0, noiseData: null, nativeId: "scan=" + (s + 1));
            }

            //Create MS2 scans with two sets of fragments each aligning perfectly with one of the precursors in FakeMs1Scans
            FakeMs2Scans = new MsDataScan[10];
            double[] intensityMultipliers2 = { 0.1, 0.2, 0.3, 0.2, 0.1, 0.1, 0.2, 0.3, 0.2, 0.1 };
            var fragArray1 = FragDist.Masses.SelectMany(v => new List<double> { v.ToMz(1) }).ToArray();
            var fragArray2 = FragDist.Masses.SelectMany(v => new List<double> { (v + 1).ToMz(1) }).ToArray();
            double[][] frags = Enumerable.Repeat(fragArray1, 5).Concat(Enumerable.Repeat(fragArray2, 5)).ToArray();
            for (int s = 0; s < FakeMs2Scans.Length; s++)
            {
                double[] intensities = FragDist.Intensities.SelectMany(v => new List<double> { v * Intensity * intensityMultipliers2[s] }).ToArray();
                FakeMs2Scans[s] = new MsDataScan(massSpectrum: new MzSpectrum(frags[s], intensities, false), oneBasedScanNumber: s + 1, msnOrder: 2, isCentroid: true,
                    polarity: Polarity.Positive, retentionTime: 1.01 + s / 10.0, scanWindowRange: new MzRange(400, 1600), scanFilter: "f",
                    mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: intensities.Sum(), injectionTime: 1.0, noiseData: null, nativeId: "scan=" + (s + 1));
            }
        }
        [Test]
        public static void DIAScanWindowMapTest()
        {
            string dataPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData\\DIA\\18300_REP2_500ng_HumanLysate_SWATH_1_RT25.63-25.81.mzML");
            var myFileManager = new MyFileManager(true);
            var dataFile = myFileManager.LoadFile(dataPath, new CommonParameters());
            var allMs2Scans = dataFile.GetAllScansList().Where(s => s.MsnOrder == 2).ToArray();
            var diaScanWindowMap = DIAEngine.ConstructMs2Groups(allMs2Scans);
            Assert.That(diaScanWindowMap.Count, Is.EqualTo(34));
            foreach(var group in diaScanWindowMap)
            {
                Assert.That(group.Value.Count, Is.EqualTo(3));
            }
        }

        [Test]
        public static void TestXicConstructor()
        {

        }

        [Test]
        public static void PseudoScanConversionTest()
        {
            var deconParameters = new ClassicDeconvolutionParameters(1, 20, 4, 3);
            var ms1XicConstructor = new NeutralMassXicConstructor(new PpmTolerance(20), 2, 1, 3, deconParameters);
            var ms1Xics = ms1XicConstructor.GetAllXics(FakeMs1Scans);

            //Test with PseudoMs2ConstructionType.MzPeak when we use mzPeak indexing on Ms2 scans
            var ms2XicConstructor = new MzPeakXicConstructor(new PpmTolerance(5), 2, 1, 3);
            var ms2Xics = ms2XicConstructor.GetAllXics(FakeMs2Scans);
            var xicGroupingEngine = new XicGrouping(0.1f, 0.5, 0.5);
            var allGroups = xicGroupingEngine.PrecursorFragmentGrouping(ms1Xics, ms2Xics);
            var pseudoScan1 = PrecursorFragmentsGroup.GetPseudoMs2ScanFromPfGroup(allGroups[0], PseudoMs2ConstructionType.MzPeak, new CommonParameters(), "test");
            //The precursor information should match with PrecursorDist
            Assert.That(pseudoScan1.PrecursorCharge, Is.EqualTo(1));
            Assert.That(pseudoScan1.PrecursorMass, Is.EqualTo(PrecursorDist.Masses.First()).Within(0.001));
            Assert.That(pseudoScan1.PrecursorMonoisotopicPeakMz, Is.EqualTo(PrecursorDist.Masses.First().ToMz(1)).Within(0.01));
            //There are 5 mz peaks that should be paired with the first precursor
            Assert.That(pseudoScan1.TheScan.MassSpectrum.XArray.Length, Is.EqualTo(5));
            //The 5 peaks that belong to the same isotopic distribution should give one deconvolution result, and the deconvoluted mass should agree with FragDist from which it is generated
            Assert.That(pseudoScan1.ExperimentalFragments.Count(), Is.EqualTo(1));
            Assert.That(pseudoScan1.ExperimentalFragments.First().MonoisotopicMass, Is.EqualTo(FragDist.Masses.First()).Within(0.001));

            //Test with PseudoMs2ConstructionType.Mass when we use mass indexing on Ms2 scans
            ms2Xics = ms1XicConstructor.GetAllXics(FakeMs2Scans);
            allGroups = xicGroupingEngine.PrecursorFragmentGrouping(ms1Xics, ms2Xics);
            var pseudoScan2 = PrecursorFragmentsGroup.GetPseudoMs2ScanFromPfGroup(allGroups[0], PseudoMs2ConstructionType.Mass, new CommonParameters(), "test");
            //The precursor information should still match
            Assert.That(pseudoScan2.PrecursorCharge, Is.EqualTo(1));
            Assert.That(pseudoScan2.PrecursorMass, Is.EqualTo(PrecursorDist.Masses.First()).Within(0.001));
            Assert.That(pseudoScan2.PrecursorMonoisotopicPeakMz, Is.EqualTo(PrecursorDist.Masses.First().ToMz(1)).Within(0.01));
            //Because ms2Xics are now in neutral mass space, there should be only one fragment paired with the precursor
            Assert.That(pseudoScan2.TheScan.MassSpectrum.XArray.Length, Is.EqualTo(1));
            Assert.That(pseudoScan2.ExperimentalFragments.Count(), Is.EqualTo(1));
            Assert.That(pseudoScan1.ExperimentalFragments.First().MonoisotopicMass, Is.EqualTo(FragDist.Masses.First()).Within(0.001));
        }
    }
}