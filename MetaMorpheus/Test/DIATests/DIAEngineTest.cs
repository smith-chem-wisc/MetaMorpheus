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
using Readers;
using MathNet.Numerics.Interpolation;

namespace Test.DIATests
{
    public class DIAEngineTest
    {
        public static MsDataScan[] GetSimpleFakeScans(string peptideSequence,  double[] intensityMultipliers, double intensity, double retentionTimeStart, int msLevel, out IsotopicDistribution isotopicDistribution, int charge = 1)
        {
            ChemicalFormula cf = new Proteomics.AminoAcidPolymer.Peptide(peptideSequence).GetChemicalFormula();
            isotopicDistribution = IsotopicDistribution.GetDistribution(cf, 0.125, 1e-4);
            var fakeScans = new MsDataScan[intensityMultipliers.Length];
            for (int s = 0; s < fakeScans.Length; s++)
            {
                double[] mzs = isotopicDistribution.Masses.SelectMany(v => new List<double> { v.ToMz(charge) }).ToArray();
                double[] intensities = isotopicDistribution.Intensities.SelectMany(v => new List<double> { v * intensity * intensityMultipliers[s] }).ToArray();
                fakeScans[s] = new MsDataScan(massSpectrum: new MzSpectrum(mzs, intensities, false), oneBasedScanNumber: s + 1, msnOrder: msLevel, isCentroid: true,
                    polarity: Polarity.Positive, retentionTime: retentionTimeStart + s / 10.0, scanWindowRange: new MzRange(400, 1600), scanFilter: "f",
                    mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: intensities.Sum(), injectionTime: 1.0, noiseData: null, nativeId: "scan=" + (s + 1), dissociationType: DissociationType.Unknown);
            }
            return fakeScans;
        }

        [Test]
        public static void DIAScanWindowMapTest()
        {
            //This is a snip file from one of the published DIA datasets; the snip file has 34 different precursor isolation windows, each containing 3 MS2 scans
            string dataPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData\\DIA\\18300_REP2_500ng_HumanLysate_SWATH_1_RT25.63-25.81.mzML");
            var myFileManager = new MyFileManager(true);
            var dataFile = myFileManager.LoadFile(dataPath, new CommonParameters());
            var allMs2Scans = dataFile.GetAllScansList().Where(s => s.MsnOrder == 2).ToArray();
            var diaScanWindowMap = DIAEngine.ConstructMs2Groups(allMs2Scans);
            //There should be 34 different DIA windows
            Assert.That(diaScanWindowMap.Count, Is.EqualTo(34));
            foreach(var group in diaScanWindowMap)
            {
                //each window has 3 MS2 scans
                Assert.That(group.Value.Count, Is.EqualTo(3));
            }
        }

        [Test]
        public static void TestXicConstructor()
        {
            //Test Xic construction with mass indexing on fake scans which contain two different masses
            double[] intensityMultipliers = { 1, 2, 3, 2, 1 };
            var fakeScans1 = GetSimpleFakeScans("PEPTIDE", intensityMultipliers, 1e6, 1.0, 1, out IsotopicDistribution preDist1);
            var fakeScans2 = GetSimpleFakeScans("PEPTIDEP", intensityMultipliers, 1e6, 1.0 + 0.5, 1, out IsotopicDistribution preDist2);
            var fakeScans = fakeScans1.Concat(fakeScans2).ToArray();

            var deconParameters = new ClassicDeconvolutionParameters(1, 20, 4, 3);
            var massXicConstructor = new NeutralMassXicConstructor(new PpmTolerance(20), 2, 1, 3, deconParameters);
            var massXics = massXicConstructor.GetAllXics(fakeScans);

            //Two different masses should give 2 XICs in the neutral mass space and each XIC contains 5 peaks
            Assert.That(massXics.Count, Is.EqualTo(2));
            foreach (var xic in massXics)
            {
                Assert.That(xic.Peaks.Count, Is.EqualTo(5));
            }
            Assert.That(massXics.Any(xic => Math.Abs(xic.Peaks.First().M - preDist1.Masses.First()) < 0.001));
            Assert.That(massXics.Any(xic => Math.Abs(xic.Peaks.First().M - preDist2.Masses.First()) < 0.001));

            //Test Xic construction with mz peak indexing with the same fake scans
            var mzPeakXicConstructor = new MzPeakXicConstructor(new PpmTolerance(5), 2, 1, 3);
            var mzXics = mzPeakXicConstructor.GetAllXics(fakeScans);
            //There are six peaks in each scan of fakeScans, which should give 10 XICs and each XIC contains 5 peaks
            Assert.That(mzXics.Count, Is.EqualTo(preDist1.Masses.Count() + preDist2.Masses.Count()));
            foreach(var xic in mzXics)
            {
                Assert.That(xic.Peaks.Count, Is.EqualTo(5));
            }
            var allMzs = fakeScans.SelectMany(s => s.MassSpectrum.XArray).Distinct().ToList();
            //Ensure that the xic constructor can find XICs for all mz peaks in the ms2 scans
            Assert.That(mzXics.Sum(xic => xic.Peaks.Count), Is.EqualTo(fakeScans.Sum(s => s.MassSpectrum.Size)));
            foreach (var mz in allMzs)
            {
                Assert.That(mzXics.Any(xic => xic.Peaks.First().M == (float)mz));
            }

            //Test excetion handling in XicConstructor
            var emptyScans = new MsDataScan[0];
            var ex = Assert.Throws<MetaMorpheusException>(() => massXicConstructor.GetAllXics(emptyScans));
            Assert.That(ex.Message, Is.EqualTo("XIC construction failed."));
        }

        [Test]
        public static void TestPfGroupingEngine()
        {
            //Create MS1 scans with two precursors that do not have any overlap in their retention time
            double[] intensityMultipliers = { 1, 2, 3, 2, 1 };
            var fakeaMs1Scans1 = GetSimpleFakeScans("PEPTIDE", intensityMultipliers, 1e6, 1.0, 1, out IsotopicDistribution preDist1);
            var fakeMs1Scans2 = GetSimpleFakeScans("PEPTIDEP", intensityMultipliers, 1e6, 1.0 + 0.5, 1, out IsotopicDistribution preDist2);
            var fakeMs1Scans = fakeaMs1Scans1.Concat(fakeMs1Scans2).ToArray();

            //construct all ms1 xics using mass indexing
            var deconParameters = new ClassicDeconvolutionParameters(1, 20, 4, 3);
            var massXicConstructor = new NeutralMassXicConstructor(new PpmTolerance(20), 2, 1, 3, deconParameters);
            var ms1Xics = massXicConstructor.GetAllXics(fakeMs1Scans);

            //Create MS2 scans with two sets of fragments each aligning perfectly with one of the precursors in MS1 scans
            var fakeaMs2Scans1 = GetSimpleFakeScans("PEP", intensityMultipliers, 1e6, 1.0, 2, out IsotopicDistribution fragDist1);
            var fakeMs2Scans2 = GetSimpleFakeScans("TIDE", intensityMultipliers, 1e6, 1.0 + 0.5, 2, out IsotopicDistribution fragDist2);
            var fakeMs2Scans = fakeaMs2Scans1.Concat(fakeMs2Scans2).ToArray();

            //construct all ms2 xics using mz peak indexing
            var mzPeakXicConstructor = new MzPeakXicConstructor(new PpmTolerance(5), 2, 1, 3);
            var ms2Xics = mzPeakXicConstructor.GetAllXics(fakeMs2Scans);

            //create a XicGroupingEngine to find all precursor-fragment groups in the MS1 and MS2 scans
            var xicGroupingEngine = new XicGroupingEngine(0.1f, 0.5, 0.5);
            var allGroups = xicGroupingEngine.PrecursorFragmentGrouping(ms1Xics, ms2Xics);
            //There should be two pfGroups, one for each precursor in FakeMs1Scans
            Assert.That(allGroups.Count, Is.EqualTo(2));
            foreach (var pfGroup in allGroups)
            {
                //For each precursor XIC, half of the fragment XICs should be paired with it (which is 5 in this case)
                Assert.That(pfGroup.PFpairs.Count, Is.EqualTo(5));
                //All the fragment XICs should have the same apex as the precursor XIC
                Assert.That(pfGroup.PFpairs.All(pf => pf.FragmentXic.ApexScanIndex == pf.PrecursorXic.ApexScanIndex));
            }
        }

        [Test]
        public static void PseudoScanConversionTest()
        {
            //Make fake MS1 and MS2 scans
            double[] intensityMultipliers = { 1, 2, 3, 2, 1 };
            var fakeaMs1Scans = GetSimpleFakeScans("PEPTIDE", intensityMultipliers, 1e6, 1.0, 1, out IsotopicDistribution preDist);
            var fakeaMs2Scans = GetSimpleFakeScans("PEP", intensityMultipliers, 1e6, 1.0, 2, out IsotopicDistribution fragDist);

            var deconParameters = new ClassicDeconvolutionParameters(1, 20, 4, 3);
            var ms1XicConstructor = new NeutralMassXicConstructor(new PpmTolerance(20), 2, 1, 3, deconParameters);
            var ms1Xics = ms1XicConstructor.GetAllXics(fakeaMs1Scans);

            //Test with PseudoMs2ConstructionType.MzPeak when we use mzPeak indexing on Ms2 scans
            var ms2XicConstructor = new MzPeakXicConstructor(new PpmTolerance(5), 2, 1, 3);
            var ms2Xics = ms2XicConstructor.GetAllXics(fakeaMs2Scans);
            var xicGroupingEngine = new XicGroupingEngine(0.1f, 0.5, 0.5);
            var allGroups = xicGroupingEngine.PrecursorFragmentGrouping(ms1Xics, ms2Xics);
            var pseudoScan = PrecursorFragmentsGroup.GetPseudoMs2ScanFromPfGroup(allGroups[0], PseudoMs2ConstructionType.MzPeak, new CommonParameters(), "test");

            //The precursor information should match with PrecursorDist
            Assert.That(pseudoScan.PrecursorCharge, Is.EqualTo(1));
            Assert.That(pseudoScan.PrecursorMass, Is.EqualTo(preDist.Masses.First()).Within(0.001));
            Assert.That(pseudoScan.PrecursorMonoisotopicPeakMz, Is.EqualTo(preDist.Masses.First().ToMz(1)).Within(0.01));
            //There are 5 fragments so the pseudo spectrum should contain 5 peaks
            Assert.That(pseudoScan.TheScan.MassSpectrum.XArray.Length, Is.EqualTo(5));
            //The 5 peaks that belong to the same isotopic distribution should give one deconvolution result, and the deconvoluted mass should agree with FragDist from which it is generated
            Assert.That(pseudoScan.ExperimentalFragments.Count(), Is.EqualTo(1));
            Assert.That(pseudoScan.ExperimentalFragments.First().MonoisotopicMass, Is.EqualTo(fragDist.Masses.First()).Within(0.001));

            //Test with PseudoMs2ConstructionType.Mass when we use mass indexing on Ms2 scans
            ms2Xics = ms1XicConstructor.GetAllXics(fakeaMs2Scans);
            allGroups = xicGroupingEngine.PrecursorFragmentGrouping(ms1Xics, ms2Xics);
            pseudoScan = PrecursorFragmentsGroup.GetPseudoMs2ScanFromPfGroup(allGroups[0], PseudoMs2ConstructionType.Mass, new CommonParameters(), "test");

            //The precursor information should still match
            Assert.That(pseudoScan.PrecursorCharge, Is.EqualTo(1));
            Assert.That(pseudoScan.PrecursorMass, Is.EqualTo(preDist.Masses.First()).Within(0.001));
            Assert.That(pseudoScan.PrecursorMonoisotopicPeakMz, Is.EqualTo(preDist.Masses.First().ToMz(1)).Within(0.01));
            //Because ms2Xics are now in neutral mass space, there should be only one "peak" in the pseudo spectrum
            Assert.That(pseudoScan.TheScan.MassSpectrum.XArray.Length, Is.EqualTo(1));
            Assert.That(pseudoScan.ExperimentalFragments.Count(), Is.EqualTo(1));
            Assert.That(pseudoScan.ExperimentalFragments.First().MonoisotopicMass, Is.EqualTo(fragDist.Masses.First()).Within(0.001));
        }

        [Test]
        public static void TestDIAEngine()
        {
            //Create fake MS1 scans where there are two precursors with same mass but different charge states
            int numberOfScansPerCycle = 3;
            string precursorSequence = "PEPTIDEPEPTIDE";
            double[] intensityMultipliers = { 1, 2, 3, 2, 1 };
            ChemicalFormula cf = new Proteomics.AminoAcidPolymer.Peptide(precursorSequence).GetChemicalFormula();
            var preDist = IsotopicDistribution.GetDistribution(cf, 0.125, 1e-4);
            var fakeMs1Scans = new MsDataScan[intensityMultipliers.Length];
            for (int s = 0; s < fakeMs1Scans.Length; s++)
            {
                double[] mzs = preDist.Masses.Select(v => v.ToMz(3)).Concat(preDist.Masses.Select(v => v.ToMz(2))).ToArray();
                double[] intensities = preDist.Intensities.Select(v => v * 1e6 * intensityMultipliers[s]).Concat(preDist.Intensities.Select(v => v * 1e6 * intensityMultipliers[s])).ToArray();
                fakeMs1Scans[s] = new MsDataScan(massSpectrum: new MzSpectrum(mzs, intensities, false), 0, 1, true, Polarity.Positive, 1.0 + s / 10.0, new MzRange(400, 1600), "f", MZAnalyzerType.Orbitrap, intensities.Sum(), 1.0, null, "nativeId");
            }

            //Create two sets of fake MS2 scans, each having a different isolation window that contains one of the precursors in MS1 scans
            string fragSequence = "PEPTIDE";
            var ms2ScanGroup1 = GetSimpleFakeScans(fragSequence, intensityMultipliers, 1e5, 1.01, 2, out IsotopicDistribution fragDist1, 2); 
            var ms2ScanGroup2 = GetSimpleFakeScans(fragSequence, intensityMultipliers, 1e5, 1.02, 2, out IsotopicDistribution fragDist2, 1);

            //Put all scans together
            var allScans = new MsDataScan[intensityMultipliers.Length * numberOfScansPerCycle];
            for (int i = 0; i < intensityMultipliers.Length; i++)
            {
                allScans[numberOfScansPerCycle * i] = fakeMs1Scans[i];
                allScans[numberOfScansPerCycle * i + 1] = ms2ScanGroup1[i];
                allScans[numberOfScansPerCycle * i + 2] = ms2ScanGroup2[i];
            }

            //Set oneBasedScanNumber for all scans
            int oneBasedScanNumber = 1;
            foreach (var scan in allScans) scan.SetOneBasedScanNumber(oneBasedScanNumber++);
            //Set isolation range and precursor scan number for MS2 scans; assuming that ms2ScanGroup1 are fragments from precursor with charge 3 and ms2ScanGroup2 are fragments from precursor with charge 2
            foreach (var scan in ms2ScanGroup1)
            {
                scan.SetIsolationRange(preDist.Masses.First().ToMz(3) - 5, preDist.Masses.First().ToMz(3) + 5);
                scan.SetOneBasedPrecursorScanNumber(scan.OneBasedScanNumber - 1);
            }
            foreach (var scan in ms2ScanGroup2)
            {
                scan.SetIsolationRange(preDist.Masses.First().ToMz(2) - 5, preDist.Masses.First().ToMz(2) + 5);
                scan.SetOneBasedPrecursorScanNumber(scan.OneBasedScanNumber - 2);
            }
            var testMsDataFile = new GenericMsDataFile(allScans, new SourceFile("no nativeID format", "mzML format", null, null, null));

            var DIAparams = new DIAparameters(new NeutralMassXicConstructor(new PpmTolerance(20), 2, 1, 3, new ClassicDeconvolutionParameters(1, 20, 4, 3)), new MzPeakXicConstructor(new PpmTolerance(5), 2, 1, 3), new XicGroupingEngine(0.1f, 0.5, 0.5), PseudoMs2ConstructionType.MzPeak);
            var diaEngine = new DIAEngine(DIAparams, testMsDataFile, new CommonParameters(), new List<(string FileName, CommonParameters Parameters)>(), new List<string> { "test" });
            diaEngine.Run();

            //There should be two precursor-fragment groups, resulting in two pseudo scans
            Assert.That(diaEngine.PseudoMs2Scans.Count, Is.EqualTo(2));
            foreach(var pseudoScan in diaEngine.PseudoMs2Scans)
            {
                //The precursor charge should match with the charge of the precursor in MS1 scans
                Assert.That(pseudoScan.PrecursorCharge, Is.EqualTo(2).Or.EqualTo(3));
                //The precursor mass should match with the precursor mass in MS1 scans
                Assert.That(pseudoScan.PrecursorMass, Is.EqualTo(preDist.Masses.First()).Within(0.01));
                //The number of peaks in pseudo scan should match with the number of fake fragments
                Assert.That(pseudoScan.TheScan.MassSpectrum.XArray.Length, Is.EqualTo(fragDist1.Masses.Count()));
                //The fragments in each pseudo scan are all from the same isotopic envelope, which should lead to one neutral mass; so there should be one neutral experimental fragment in each Ms2WithSpecificMass
                Assert.That(pseudoScan.ExperimentalFragments.Count(), Is.EqualTo(1));
                //The neutral mass of the fragment should match with the fragment mass
                Assert.That(pseudoScan.ExperimentalFragments.First().MonoisotopicMass, Is.EqualTo(fragDist1.Masses.First()).Within(0.001));
            }
        }
    }
}