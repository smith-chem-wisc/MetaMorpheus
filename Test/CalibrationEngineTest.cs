using EngineLayer;
using EngineLayer.Calibration;
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
            Protein ParentProtein = new Protein("MQQQQQQQ", null, oneBasedPossibleLocalizedModifications, null, null, null, null, null, false, false, null);
            IEnumerable<ModificationWithMass> fixedModifications = new List<ModificationWithMass>();
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);

            PeptideWithPossibleModifications modPep = ParentProtein.Digest(protease, 0, InitiatorMethionineBehavior.Variable, fixedModifications).First();
            //Dictionary<int, MorpheusModification> twoBasedVariableAndLocalizeableModificationss = new Dictionary<int, MorpheusModification>();
            List<ModificationWithMass> variableModifications = new List<ModificationWithMass>();
            PeptideWithSetModifications pepWithSetMods = modPep.GetPeptidesWithSetModifications(variableModifications, 4096, 3).First();

            IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile = new TestDataFile(pepWithSetMods);

            Tolerance fragmentTolerance = new Tolerance(ToleranceUnit.Absolute, 0.01);

            List<NewPsmWithFdr> identifications = new List<NewPsmWithFdr>();
            PsmParent newPsm = new TestParentSpectrumMatch(2, 2);
            PsmWithMultiplePossiblePeptides thisPSM = new PsmWithMultiplePossiblePeptides(newPsm, new HashSet<PeptideWithSetModifications>() { pepWithSetMods }, fragmentTolerance, myMsDataFile, new List<ProductType> { ProductType.B, ProductType.Y });
            NewPsmWithFdr thePsmwithfdr = new NewPsmWithFdr(thisPSM);
            thePsmwithfdr.SetValues(1, 0, 0, 1, 0, 0);
            identifications.Add(thePsmwithfdr);

            int minMS1isotopicPeaksNeededForConfirmedIdentification = 3;
            int minMS2isotopicPeaksNeededForConfirmedIdentification = 2;
            int numFragmentsNeededForEveryIdentification = 10;

            var calibrationEngine = new CalibrationEngine(myMsDataFile, fragmentTolerance, identifications, minMS1isotopicPeaksNeededForConfirmedIdentification, minMS2isotopicPeaksNeededForConfirmedIdentification, numFragmentsNeededForEveryIdentification, new Tolerance(ToleranceUnit.PPM, 10), FragmentTypes.b | FragmentTypes.y, (List<LabeledMs1DataPoint> theList, string s) => {; }, (List<LabeledMs2DataPoint> theList, string s) => {; }, true, new List<string> ());

            var res = calibrationEngine.Run();
            Assert.IsTrue(res is CalibrationResults);
        }

        [Test]
        public static void TestQuadratic()
        {
            PeptideWithSetModifications pepWithSetMods = new Protein("MQQQQQQQ", null, new Dictionary<int, List<Modification>>(), null, null, null, null, null, false, false, null).Digest(new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null), 0, InitiatorMethionineBehavior.Variable, new List<ModificationWithMass>()).First().GetPeptidesWithSetModifications(new List<ModificationWithMass>(), 4096, 3).First();

            IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile = new TestDataFile(pepWithSetMods, "quadratic");

            Tolerance fragmentTolerance = new Tolerance(ToleranceUnit.Absolute, 0.1);

            NewPsmWithFdr thePsmwithfdr = new NewPsmWithFdr(new PsmWithMultiplePossiblePeptides(new TestParentSpectrumMatch(2, 2), new HashSet<PeptideWithSetModifications>() { pepWithSetMods }, fragmentTolerance, myMsDataFile, new List<ProductType> { ProductType.B, ProductType.Y }));
            thePsmwithfdr.SetValues(1, 0, 0, 1, 0, 0);

            var res = new CalibrationEngine(myMsDataFile, fragmentTolerance, new List<NewPsmWithFdr> { thePsmwithfdr }, 3, 2, 10, new Tolerance(ToleranceUnit.PPM, 10), FragmentTypes.b | FragmentTypes.y, (List<LabeledMs1DataPoint> theList, string s) => {; }, (List<LabeledMs2DataPoint> theList, string s) => {; }, true, new List<string>()).Run();
            Assert.IsTrue(res is CalibrationResults);
        }

        #endregion Public Methods

    }
}