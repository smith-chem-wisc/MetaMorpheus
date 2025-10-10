using EngineLayer.DIA;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using EngineLayer;
using TaskLayer;
using System.IO;
using Chemistry;
using Readers;
using System.Drawing.Imaging;

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
        public static void DIAScanWindowMapTestForInaccurateIsoWindow()
        {
            int numberOfScansPerCycle = 6;
            int numberOfCycles = 5;
            var isoMz = new double[] { 500, 600, 700, 800, 900, 1000 };
            double isoWidth = 50;
            //This mannually sets the isolation windows to be slightly different for each cycle to simulate the fluctuation of iso window selection or instability of doubles
            var isoMzOffset = new double[] { 0.0, 0.01, 0.02, 0.03, -0.01};
            var ms2Scans = new MsDataScan[numberOfScansPerCycle * numberOfCycles];
            for (int i = 0; i < numberOfCycles; i++)
            {
                for (int j = 0; j < numberOfScansPerCycle; j++)
                {
                    ms2Scans[i * numberOfScansPerCycle + j] = new MsDataScan(new MzSpectrum(new double[] { }, new double[] { }, false), i*j + j, 2, true, Polarity.Positive, 1.0 + (i * j + j)/10, new MzRange(400, 1600), "f", MZAnalyzerType.Orbitrap, 0, 1.0, null, null, isolationWidth: isoWidth, isolationMZ: isoMz[j] + isoMzOffset[i]);
                }
            }
            var diaScanWindowMap = DIAEngine.ConstructMs2Groups(ms2Scans);

            Assert.That(diaScanWindowMap.Count, Is.EqualTo(6));
            foreach (var group in diaScanWindowMap)
            {
                //each window has 5 MS2 scans
                Assert.That(group.Value.Count, Is.EqualTo(5));
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

            //Test exception handling in XicConstructor when the input scans are empty
            var emptyScans = new MsDataScan[0];
            var ex = Assert.Throws<MetaMorpheusException>(() => massXicConstructor.GetAllXics(emptyScans));
            Assert.That(ex.Message, Is.EqualTo("XIC construction failed."));
        }

        [Test]
        public static void TestXicSplineForAllXics()
        {
            //Test Xic spline with fake scans
            double[] intensityMultipliers = { 1, 2, 3, 2, 1 };
            var fakeScans1 = GetSimpleFakeScans("PEPTIDE", intensityMultipliers, 1e6, 1.0, 1, out IsotopicDistribution preDist1);
            var xicLinearSpline = new XicLinearSpline(0.05);
            var mzPeakXicConstructor = new MzPeakXicConstructor(new PpmTolerance(5), 2, 1, 3, xicLinearSpline);
            var mzXics = mzPeakXicConstructor.GetAllXicsWithXicSpline(fakeScans1);
            //GetAllXicsAndXicSpline should return all Xics in a given set of scans and set XYData for each XIC if XicSplineEngine is defined in the XicConstructor
            foreach (var xic in mzXics)
            {
                Assert.That(xic.XYData, Is.Not.Null);
                Assert.That(xic.XYData.Length, Is.GreaterThan(intensityMultipliers.Length));
            }
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
            var xicGroupingEngine = new XicGroupingEngine(0.1f, 0.5, 0.5, maxThreadsForGrouping: 1);
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
            var allGroups = xicGroupingEngine.PrecursorFragmentGrouping(ms1Xics, ms2Xics).ToArray();
            var pseudoScan = allGroups[0].GetPseudoMs2ScanFromPfGroup(PseudoMs2ConstructionType.MzPeak, new CommonParameters(), "test");

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
            allGroups = xicGroupingEngine.PrecursorFragmentGrouping(ms1Xics, ms2Xics).ToArray();
            pseudoScan = allGroups[0].GetPseudoMs2ScanFromPfGroup(PseudoMs2ConstructionType.Mass, new CommonParameters(), "test");

            //The precursor information should still match
            Assert.That(pseudoScan.PrecursorCharge, Is.EqualTo(1));
            Assert.That(pseudoScan.PrecursorMass, Is.EqualTo(preDist.Masses.First()).Within(0.001));
            Assert.That(pseudoScan.PrecursorMonoisotopicPeakMz, Is.EqualTo(preDist.Masses.First().ToMz(1)).Within(0.01));
            //Because ms2Xics are now in neutral mass space, there should be only one "peak" in the pseudo spectrum
            Assert.That(pseudoScan.TheScan.MassSpectrum.XArray.Length, Is.EqualTo(1));
            Assert.That(pseudoScan.ExperimentalFragments.Count(), Is.EqualTo(1));
            Assert.That(pseudoScan.ExperimentalFragments.First().MonoisotopicMass, Is.EqualTo(fragDist.Masses.First()).Within(0.001));

            //Test exception handling
            var invalidType = (PseudoMs2ConstructionType)999;
            var ex = Assert.Throws<MetaMorpheusException>(() => allGroups[0].GetPseudoMs2ScanFromPfGroup(invalidType, new CommonParameters(), "test"));
            Assert.That(ex.Message, Is.EqualTo("Invalid pseudo MS2 construction type specified."));
        }

        [Test]
        public static void TestDIAEngine()
        {
            //Create fake MS1 scans where there are two precursors with same mass but different charge states
            //The two precursors have same elution profiles but are fragmented in different DIA MS2 isolation windows
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

            //put the data into the search pipeline to see if it returns pseudo scans
            var DIAparams = new DIAparameters(DIAanalysisType.DIA, new NeutralMassXicConstructor(new PpmTolerance(20), 2, 1, 3, new ClassicDeconvolutionParameters(1, 20, 4, 3)), new MzPeakXicConstructor(new PpmTolerance(5), 2, 1, 3), new XicGroupingEngine(0.15f, 0.5, 0.5, 1), PseudoMs2ConstructionType.MzPeak);
            var commonParams = new CommonParameters { DIAparameters = DIAparams };
            var pseudoScans = MetaMorpheusTask.GetMs2Scans(testMsDataFile, null, commonParams).ToArray();

            //There should be two precursor-fragment groups, resulting in two pseudo scans
            Assert.That(pseudoScans.Length, Is.EqualTo(2));
            foreach(var pseudoScan in pseudoScans)
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

        [Test]
        public static void TestExceptionHandlingForGetMs2Scans()
        {
            var emptyScans = new MsDataScan[0];
            var fakeFile = new GenericMsDataFile(emptyScans, new SourceFile("no nativeID format", "mzML format", null, null, null));
            //put the data into the search pipeline to see if it returns pseudo scans
            var invalidType = (DIAanalysisType)999;
            var DIAparams = new DIAparameters(invalidType, new NeutralMassXicConstructor(new PpmTolerance(20), 2, 1, 3, new ClassicDeconvolutionParameters(1, 20, 4, 3)), new MzPeakXicConstructor(new PpmTolerance(5), 2, 1, 3), new XicGroupingEngine(0.15f, 0.5, 0.5, 1), PseudoMs2ConstructionType.MzPeak);
            var commonParams = new CommonParameters { DIAparameters = DIAparams };

            var ex = Assert.Throws<NotImplementedException>(() => MetaMorpheusTask.GetMs2Scans(fakeFile, null, commonParams));
            Assert.That(ex.Message, Is.EqualTo("DIA analysis type not implemented."));
        }

        [Test]
        public static void IsdVoltageScanTest()
        {
            string dataPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData\\DIA\\08-12-24_PEPPI_FractionD_orbiMS1_ISD60-80-100_RT8.1-9.23.mzML");
            var myFileManager = new MyFileManager(true);
            var dataFile = myFileManager.LoadFile(dataPath, new CommonParameters());
            var allScans = dataFile.GetAllScansList().ToArray();

            //Read in all scans, relabel ISD fragment scans and construct "MS2" groups
            //In this data file, one scan cycle contains four scans: MS1, 60, 80, 100
            var isdVoltageMap = ISDEngine.ConstructIsdGroups(allScans, out MsDataScan[] ms1Scans);
            ISDEngine.ReLabelIsdScans(isdVoltageMap, ms1Scans);

            //There are three different ISD voltages 60, 80, and 100, so there should be three kvps in the isdVoltageMap
            Assert.That(isdVoltageMap.Count, Is.EqualTo(3));
            //The number of scans in each ISD voltage group should be all the same as the number of ms1 scans in this file
            foreach(var kvp in isdVoltageMap)
            {
                Assert.That(kvp.Value.Count, Is.EqualTo(ms1Scans.Length));
            }
            //The ISD fragment scans should be relabeled as ms2 scan
            Assert.That(isdVoltageMap.Values.All(v => v.All(scan => scan.MsnOrder == 2)));
            //The ISD fragment scans should have a valid precursor scan number
            Assert.That(isdVoltageMap.Values.All(v => v.All(scan => scan.OneBasedPrecursorScanNumber.HasValue)));
            //All ms1 scans should stay as ms1 scans
            Assert.That(ms1Scans.All(scan => scan.MsnOrder == 1));
            //The isolation range of each ISD fragment scan should match with the scan window of the corresponding MS1 scan
            foreach(var scan in isdVoltageMap.Values.SelectMany(v => v))
            {
                var ms1Scan = ms1Scans.Where(s => s.OneBasedScanNumber == scan.OneBasedPrecursorScanNumber.Value).FirstOrDefault();
                Assert.That(scan.IsolationRange.Minimum, Is.EqualTo(ms1Scan.ScanWindowRange.Minimum));
                Assert.That(scan.IsolationRange.Maximum, Is.EqualTo(ms1Scan.ScanWindowRange.Maximum));
            }
        }

        [Test]
        public static void ISDEngineTest()
        {
            //Create fake scans with a fake precursor and three fake fragments that appear in isd60, 80, 100 scans respectively
            //There are four masses which should give four distinct XICs
            var sequences = new List<string> { "PEPTIDEPEPTIDE", "PEPTIDEPEP", "PEPTIDE", "PEP" };
            //create five scan cycles, each with a different intensity level; (ms1, isd60, isd80, isd100) x 5
            double[] intensityMultipliers = { 1, 2, 3, 2, 1 };
            var cfs = sequences.Select(p => new Proteomics.AminoAcidPolymer.Peptide(p).GetChemicalFormula()).ToArray();
            var dists = cfs.Select(cf => IsotopicDistribution.GetDistribution(cf, 0.125, 1e-4)).ToArray();
            var fakeScans = new MsDataScan[intensityMultipliers.Length * cfs.Length];
            var fakeScanFilters = new string[] { "sid=15", "sid=60", "sid=80", "sid=100" };
            var rtStart = new double[] { 1.0, 1.1, 1.2, 1.3, 1.4};

            //There should be 20 fake scans in total
            int oneBasedScanNumber = 1;
            for (int i = 0; i < intensityMultipliers.Length; i++)
            {
                for (int j = 0; j < fakeScanFilters.Length; j++)
                {
                    double[] mzs = dists[j].Masses.Select(v => v.ToMz(2)).ToArray();
                    double[] intensities = dists[j].Intensities.Select(v => v * 1e6 * intensityMultipliers[i]).ToArray();
                    fakeScans[oneBasedScanNumber - 1] = new MsDataScan(massSpectrum: new MzSpectrum(mzs, intensities, false), oneBasedScanNumber: oneBasedScanNumber, msnOrder: 1, isCentroid: true, polarity: Polarity.Positive, retentionTime: rtStart[i] + j / 100.0, scanWindowRange: new MzRange(100, 2000), scanFilter: fakeScanFilters[j], MZAnalyzerType.Orbitrap, totalIonCurrent: intensities.Sum(), injectionTime: 1.0, noiseData: null, nativeId: "scan=" + oneBasedScanNumber, dissociationType: DissociationType.Unknown); 
                    oneBasedScanNumber++;
                }
            }
            var testMsDataFile = new GenericMsDataFile(fakeScans, new SourceFile("no nativeID format", "mzML format", null, null, null));
            
            //put the data into the search pipeline to see if it returns pseudo scans
            var DIAparams = new DIAparameters(DIAanalysisType.ISD, new NeutralMassXicConstructor(new PpmTolerance(20), 2, 1, 3, new ClassicDeconvolutionParameters(1, 20, 4, 3)), new NeutralMassXicConstructor(new PpmTolerance(20), 2, 1, 3, new ClassicDeconvolutionParameters(1, 20, 4, 3)), new XicGroupingEngine(0.2f, 0.5, 0.5, 1), PseudoMs2ConstructionType.Mass);
            var commonParams = new CommonParameters { DIAparameters = DIAparams };
            var pseudoScans = MetaMorpheusTask.GetMs2Scans(testMsDataFile, null, commonParams).ToArray();

            //Without combining the fragments, each ISD level should give one pseudo scan and each scan contains a different fake fragment, so there should be three pseudo scans in total
            Assert.That(pseudoScans.Length, Is.EqualTo(3));
            for (int i = 0; i < pseudoScans.Length; i++)
            {
                //The precursor charge should be 2
                Assert.That(pseudoScans[i].PrecursorCharge, Is.EqualTo(2));
                //The precursor mass should match with the fake precursor
                Assert.That(pseudoScans[i].PrecursorMass, Is.EqualTo(cfs[0].MonoisotopicMass).Within(0.01));
                //Each pseudo scan should have one fragment mass
                Assert.That(pseudoScans[i].ExperimentalFragments.Count, Is.EqualTo(1));
                //The fragment mass in each pseudo scan should match with the fake fragments above
                Assert.That(pseudoScans[i].ExperimentalFragments.First().MonoisotopicMass, Is.EqualTo(cfs[i + 1].MonoisotopicMass).Within(0.01));
            }

            //If we combine the fragments, there should be one pseudo scan with three fragments
            DIAparams = new DIAparameters(DIAanalysisType.ISD, new NeutralMassXicConstructor(new PpmTolerance(20), 2, 1, 3, new ClassicDeconvolutionParameters(1, 20, 4, 3)), new NeutralMassXicConstructor(new PpmTolerance(20), 2, 1, 3, new ClassicDeconvolutionParameters(1, 20, 4, 3)), new XicGroupingEngine(0.2f, 0.5, 0.5, 1), PseudoMs2ConstructionType.Mass, combineFragments: true);
            commonParams = new CommonParameters { DIAparameters = DIAparams };
            pseudoScans = MetaMorpheusTask.GetMs2Scans(testMsDataFile, null, commonParams).ToArray();
            Assert.That(pseudoScans.Length, Is.EqualTo(1));
            Assert.That(pseudoScans[0].PrecursorCharge, Is.EqualTo(2));
            Assert.That(pseudoScans[0].PrecursorMass, Is.EqualTo(cfs[0].MonoisotopicMass).Within(0.01));
            Assert.That(pseudoScans[0].ExperimentalFragments.Count, Is.EqualTo(3));
        }
    }
}