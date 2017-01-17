using InternalLogicCalibration;
using InternalLogicEngineLayer;
using MassSpectrometry;
using NUnit.Framework;
using OldInternalLogic;
using Spectra;
using System.Collections.Generic;

namespace Test
{
    [TestFixture]
    public class CalibrationEngineTests
    {
        #region Public Methods

        [Test]
        public void TestCalibrationEngine()
        {
            Dictionary<int, List<MorpheusModification>> oneBasedPossibleLocalizedModifications = new Dictionary<int, List<MorpheusModification>>();
            Protein ParentProtein = new Protein("MQQQQQQQ", null, null, oneBasedPossibleLocalizedModifications, null, null, null, null, null, 0, false);
            PeptideWithPossibleModifications modPep = new PeptideWithPossibleModifications(1, 8, ParentProtein, 0, "kk");
            Dictionary<int, MorpheusModification> twoBasedVariableAndLocalizeableModificationss = new Dictionary<int, MorpheusModification>();
            PeptideWithSetModifications pepWithSetMods = new PeptideWithSetModifications(modPep, twoBasedVariableAndLocalizeableModificationss);

            IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile = new TestDataFile(pepWithSetMods);

            Tolerance fragmentTolerance = new Tolerance(ToleranceUnit.Absolute, 0.01);
            double toleranceInMZforMS2Search = fragmentTolerance.Value;

            List<NewPsmWithFDR> identifications = new List<NewPsmWithFDR>();
            ParentSpectrumMatch newPsm = new TestParentSpectrumMatch(2, 2);
            PSMwithTargetDecoyKnown thisPSM = new PSMwithTargetDecoyKnown(newPsm, new HashSet<PeptideWithSetModifications>() { pepWithSetMods }, fragmentTolerance, myMsDataFile);
            NewPsmWithFDR thePsmwithfdr = new NewPsmWithFDR(thisPSM, 1, 0, 0);
            identifications.Add(thePsmwithfdr);

            int randomSeed = 0;

            var calibrationEngine = new CalibrationEngine(myMsDataFile, randomSeed, toleranceInMZforMS2Search, identifications);

            var res = calibrationEngine.Run();
            Assert.IsTrue(res is CalibrationResults);
        }

        #endregion Public Methods
    }
}