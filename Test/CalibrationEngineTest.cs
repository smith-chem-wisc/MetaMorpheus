using EngineLayer;
using EngineLayer.Calibration;
using IO.MzML;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    [TestFixture]
    public class CalibrationEngineTests
    {
        #region Public Methods

        [Test]
        public static void TestCalibrationEngine()
        {
            Protein ParentProtein = new Protein("MQQQQQQQ", null);
            IEnumerable<ModificationWithMass> fixedModifications = new List<ModificationWithMass>();
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);

            PeptideWithPossibleModifications modPep = ParentProtein.Digest(protease, 0, null, null, InitiatorMethionineBehavior.Variable, fixedModifications).First();
            //Dictionary<int, MorpheusModification> twoBasedVariableAndLocalizeableModificationss = new Dictionary<int, MorpheusModification>();
            List<ModificationWithMass> variableModifications = new List<ModificationWithMass>();
            PeptideWithSetModifications pepWithSetMods = modPep.GetPeptidesWithSetModifications(variableModifications, 4096, 3).First();

            IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile = new TestDataFile(pepWithSetMods);

            Tolerance fragmentTolerance = new AbsoluteTolerance(0.01);

            List<Psm> identifications = new List<Psm>();
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(new MzmlScanWithPrecursor(2, new MzmlMzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 1, null, null), new MzPeak(0, 0), 2, null);
            Psm newPsm = new Psm(pepWithSetMods.CompactPeptide(TerminusType.None), 0, 0, 0, scan);

            Dictionary<ModificationWithMass, ushort> modsDictionary = new Dictionary<ModificationWithMass, ushort>();
            Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> matching = new Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>
            {
                {pepWithSetMods.CompactPeptide(TerminusType.None), new HashSet<PeptideWithSetModifications>{ pepWithSetMods } }
            };
            List<ProductType> lp = new List<ProductType> { ProductType.B, ProductType.Y };
            newPsm.MatchToProteinLinkedPeptides(matching);

            newPsm.SetFdrValues(1, 0, 0, 1, 0, 0);
            identifications.Add(newPsm);

            int minMS1isotopicPeaksNeededForConfirmedIdentification = 3;
            int minMS2isotopicPeaksNeededForConfirmedIdentification = 2;
            int numFragmentsNeededForEveryIdentification = 10;

            var calibrationEngine = new CalibrationEngine(myMsDataFile, fragmentTolerance, identifications, minMS1isotopicPeaksNeededForConfirmedIdentification, minMS2isotopicPeaksNeededForConfirmedIdentification, numFragmentsNeededForEveryIdentification, new PpmTolerance(10), FragmentTypes.b | FragmentTypes.y, (List<LabeledMs1DataPoint> theList, string s) => {; }, (List<LabeledMs2DataPoint> theList, string s) => {; }, true, new List<string>());

            var res = calibrationEngine.Run();
            Assert.IsTrue(res is CalibrationResults);
        }

        [Test]
        public static void TestQuadratic()
        {
            PeptideWithSetModifications pepWithSetMods = new Protein("MQQQQQQQ", null).Digest(new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null), 0, null, null, InitiatorMethionineBehavior.Variable, new List<ModificationWithMass>()).First().GetPeptidesWithSetModifications(new List<ModificationWithMass>(), 4096, 3).First();

            IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile = new TestDataFile(pepWithSetMods, "quadratic");

            Tolerance fragmentTolerance = new AbsoluteTolerance(0.1);

            IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>> dfd = new MzmlScanWithPrecursor(2, new MzmlMzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 1, null, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfd, new MzPeak(2, 2), 2, null);
            Psm newPsm = new Psm(pepWithSetMods.CompactPeptide(TerminusType.None), 0, 0, 0, scan);

            Dictionary<ModificationWithMass, ushort> modsDictionary = new Dictionary<ModificationWithMass, ushort>();
            Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> matching = new Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>
            {
                {pepWithSetMods.CompactPeptide(TerminusType.None), new HashSet<PeptideWithSetModifications>{ pepWithSetMods } }
            };
            List<ProductType> lp = new List<ProductType> { ProductType.B, ProductType.Y };
            newPsm.MatchToProteinLinkedPeptides(matching);

            newPsm.SetFdrValues(1, 0, 0, 1, 0, 0);

            var res = new CalibrationEngine(myMsDataFile, fragmentTolerance, new List<Psm> { newPsm }, 3, 2, 10, new PpmTolerance(10), FragmentTypes.b | FragmentTypes.y, (List<LabeledMs1DataPoint> theList, string s) => {; }, (List<LabeledMs2DataPoint> theList, string s) => {; }, true, new List<string>()).Run();
            Assert.IsTrue(res is CalibrationResults);
        }

        #endregion Public Methods
    }
}