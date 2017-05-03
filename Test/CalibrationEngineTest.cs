using EngineLayer;
using EngineLayer.Calibration;
using EngineLayer.ClassicSearch;
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
            var oneBasedPossibleLocalizedModifications = new Dictionary<int, List<Modification>>();
            Protein ParentProtein = new Protein("MQQQQQQQ", null, null, oneBasedPossibleLocalizedModifications, null, null, null, null, null, false, false, null);
            IEnumerable<ModificationWithMass> fixedModifications = new List<ModificationWithMass>();
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);

            PeptideWithPossibleModifications modPep = ParentProtein.Digest(protease, 0, null, null, InitiatorMethionineBehavior.Variable, fixedModifications).First();
            //Dictionary<int, MorpheusModification> twoBasedVariableAndLocalizeableModificationss = new Dictionary<int, MorpheusModification>();
            List<ModificationWithMass> variableModifications = new List<ModificationWithMass>();
            PeptideWithSetModifications pepWithSetMods = modPep.GetPeptidesWithSetModifications(variableModifications, 4096, 3).First();

            IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile = new TestDataFile(pepWithSetMods);

            Tolerance fragmentTolerance = new Tolerance(ToleranceUnit.Absolute, 0.01);

            List<NewPsmWithFdr> identifications = new List<NewPsmWithFdr>();
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(new MzmlScanWithPrecursor(2, new MzmlMzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 1, null, null), new MzPeak(0, 0), 2, null);
            PsmParent newPsm = new PsmClassic(pepWithSetMods, 0, 0, 0, scan);

            Dictionary<ModificationWithMass, ushort> modsDictionary = new Dictionary<ModificationWithMass, ushort>();
            Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> matching = new Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>>
            {
                {newPsm.GetCompactPeptide(modsDictionary), new HashSet<PeptideWithSetModifications>{ pepWithSetMods } }
            };
            List<ProductType> lp = new List<ProductType> { ProductType.B, ProductType.Y };
            newPsm.GetTheActualPeptidesWithSetModificationsAndComputeStuff(matching, fragmentTolerance, scan, lp, modsDictionary);

            NewPsmWithFdr thePsmwithfdr = new NewPsmWithFdr(newPsm);
            thePsmwithfdr.SetValues(1, 0, 0, 1, 0, 0);
            identifications.Add(thePsmwithfdr);

            int minMS1isotopicPeaksNeededForConfirmedIdentification = 3;
            int minMS2isotopicPeaksNeededForConfirmedIdentification = 2;
            int numFragmentsNeededForEveryIdentification = 10;

            var calibrationEngine = new CalibrationEngine(myMsDataFile, fragmentTolerance, identifications, minMS1isotopicPeaksNeededForConfirmedIdentification, minMS2isotopicPeaksNeededForConfirmedIdentification, numFragmentsNeededForEveryIdentification, new Tolerance(ToleranceUnit.PPM, 10), FragmentTypes.b | FragmentTypes.y, (List<LabeledMs1DataPoint> theList, string s) => {; }, (List<LabeledMs2DataPoint> theList, string s) => {; }, true, new List<string>());

            var res = calibrationEngine.Run();
            Assert.IsTrue(res is CalibrationResults);
        }

        [Test]
        public static void TestQuadratic()
        {
            PeptideWithSetModifications pepWithSetMods = new Protein("MQQQQQQQ", null, null, new Dictionary<int, List<Modification>>(), null, null, null, null, null, false, false, null).Digest(new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null), 0, null, null, InitiatorMethionineBehavior.Variable, new List<ModificationWithMass>()).First().GetPeptidesWithSetModifications(new List<ModificationWithMass>(), 4096, 3).First();

            IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile = new TestDataFile(pepWithSetMods, "quadratic");

            Tolerance fragmentTolerance = new Tolerance(ToleranceUnit.Absolute, 0.1);

            IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>> dfd = new MzmlScanWithPrecursor(2, new MzmlMzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 1, null, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfd, new MzPeak(2, 2), 2, null);
            PsmParent newPsm = new PsmClassic(pepWithSetMods, 0, 0, 0, scan);

            Dictionary<ModificationWithMass, ushort> modsDictionary = new Dictionary<ModificationWithMass, ushort>();
            Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> matching = new Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>>
            {
                {newPsm.GetCompactPeptide(modsDictionary), new HashSet<PeptideWithSetModifications>{ pepWithSetMods } }
            };
            List<ProductType> lp = new List<ProductType> { ProductType.B, ProductType.Y };
            newPsm.GetTheActualPeptidesWithSetModificationsAndComputeStuff(matching, fragmentTolerance, scan, lp, modsDictionary);

            NewPsmWithFdr thePsmwithfdr = new NewPsmWithFdr(newPsm);
            thePsmwithfdr.SetValues(1, 0, 0, 1, 0, 0);

            var res = new CalibrationEngine(myMsDataFile, fragmentTolerance, new List<NewPsmWithFdr> { thePsmwithfdr }, 3, 2, 10, new Tolerance(ToleranceUnit.PPM, 10), FragmentTypes.b | FragmentTypes.y, (List<LabeledMs1DataPoint> theList, string s) => {; }, (List<LabeledMs2DataPoint> theList, string s) => {; }, true, new List<string>()).Run();
            Assert.IsTrue(res is CalibrationResults);
        }

        #endregion Public Methods

    }
}