using EngineLayer.Calibration;
using EngineLayer;
using MassSpectrometry;
using NUnit.Framework;

using Proteomics;
using Spectra;
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
            Dictionary<int, List<MetaMorpheusModification>> oneBasedPossibleLocalizedModifications = new Dictionary<int, List<MetaMorpheusModification>>();
            Protein ParentProtein = new Protein("MQQQQQQQ", null, oneBasedPossibleLocalizedModifications, null, null, null, null, null, 0, false, false);
            PeptideWithPossibleModifications modPep = new PeptideWithPossibleModifications(1, 8, ParentProtein, 0, "kk");
            //Dictionary<int, MorpheusModification> twoBasedVariableAndLocalizeableModificationss = new Dictionary<int, MorpheusModification>();
            List<MetaMorpheusModification> variableModifications = new List<MetaMorpheusModification>();
            PeptideWithSetModifications pepWithSetMods = modPep.GetPeptideWithSetModifications(variableModifications, 4096, 3, new List<MetaMorpheusModification>()).First();

            IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile = new TestDataFile(pepWithSetMods);

            Tolerance fragmentTolerance = new Tolerance(ToleranceUnit.Absolute, 0.01);
            double toleranceInMZforMS2Search = fragmentTolerance.Value;

            List<NewPsmWithFdr> identifications = new List<NewPsmWithFdr>();
            ParentSpectrumMatch newPsm = new TestParentSpectrumMatch(2, 2);
            PsmWithMultiplePossiblePeptides thisPSM = new PsmWithMultiplePossiblePeptides(newPsm, new HashSet<PeptideWithSetModifications>() { pepWithSetMods }, fragmentTolerance, myMsDataFile, new List<ProductType> { ProductType.B, ProductType.Y });
            NewPsmWithFdr thePsmwithfdr = new NewPsmWithFdr(thisPSM, 1, 0, 0);
            identifications.Add(thePsmwithfdr);

            int randomSeed = 0;

            int minMS1isotopicPeaksNeededForConfirmedIdentification = 3;
            int minMS2isotopicPeaksNeededForConfirmedIdentification = 2;
            int numFragmentsNeededForEveryIdentification = 10;
            double toleranceInMZforMS1Search = 0.01;

            var calibrationEngine = new CalibrationEngine(myMsDataFile, randomSeed, toleranceInMZforMS2Search, identifications, minMS1isotopicPeaksNeededForConfirmedIdentification, minMS2isotopicPeaksNeededForConfirmedIdentification, numFragmentsNeededForEveryIdentification, toleranceInMZforMS1Search, FragmentTypes.b | FragmentTypes.y);

            var res = calibrationEngine.Run();
            Assert.IsTrue(res is CalibrationResults);
        }

        #endregion Public Methods

    }
}