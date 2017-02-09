using EngineLayer;
using EngineLayer.ClassicSearch;
using EngineLayer.Gptmd;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using Spectra;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    [TestFixture]
    public class LocalizationTest
    {

        #region Public Methods

        [Test]
        public static void TestLocalization()
        {
            Dictionary<int, HashSet<BaseModification>> oneBasedPossibleLocalizedModifications = new Dictionary<int, HashSet<BaseModification>>();
            string accession = null;
            int?[] beginPositions = null;
            int?[] endPositions = null;
            string[] bigPeptideTypes = null;
            string name = null;
            string fullName = null;
            IEnumerable<ModificationWithMass> allKnownFixedModifications = new List<ModificationWithMass>();

            Protein parentProteinForMatch = new Protein("MEK", accession, oneBasedPossibleLocalizedModifications, beginPositions, endPositions, bigPeptideTypes, name, fullName, 0, false, false);
            PeptideWithPossibleModifications pwpm = new PeptideWithPossibleModifications(1, 3, parentProteinForMatch, 0, null, allKnownFixedModifications);

            List<ModificationWithMass> variableModifications = new List<ModificationWithMass> { new ModificationWithMass("21", null, "E", ModificationSites.Any, 21.981943, 0, null) };

            List<PeptideWithSetModifications> allPeptidesWithSetModifications = pwpm.GetPeptideWithSetModifications(variableModifications, 2, 1).ToList();
            Assert.AreEqual(2, allPeptidesWithSetModifications.Count());
            PeptideWithSetModifications ps = allPeptidesWithSetModifications.First();

            List<ProductType> lp = new List<ProductType> { ProductType.B, ProductType.Y };

            PeptideWithSetModifications pepWithSetModsForSpectrum = allPeptidesWithSetModifications.Last();
            IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile = new TestDataFile(new List<PeptideWithSetModifications> { pepWithSetModsForSpectrum });
            Tolerance fragmentTolerance = new Tolerance(ToleranceUnit.Absolute, 0.01);

            HashSet<PeptideWithSetModifications> peptidesWithSetModifications = new HashSet<PeptideWithSetModifications> { ps };

            string fileName = null;
            double scanRetentionTime = double.NaN;
            double scanPrecursorIntensity = double.NaN;
            double scanPrecursorMass = pepWithSetModsForSpectrum.MonoisotopicMass;
            int scanPrecursorCharge = 0;
            int scanExperimentalPeaks = 0;
            double totalIonCurrent = double.NaN;
            double scanPrecursorMZ = double.NaN;
            double score = double.NaN;
            PsmParent newPsm = new PsmClassic(ps, fileName, scanRetentionTime, scanPrecursorIntensity, scanPrecursorMass, 2, scanPrecursorCharge, scanExperimentalPeaks, totalIonCurrent, scanPrecursorMZ, score, 0);

            Assert.IsNull(newPsm.LocalizedScores);
            Assert.IsNull(newPsm.matchedIonsList);
            new PsmWithMultiplePossiblePeptides(newPsm, peptidesWithSetModifications, fragmentTolerance, myMsDataFile, lp);

            // Was single peak, now three match!!!
            Assert.AreEqual(0, newPsm.matchedIonsList[ProductType.B].Count(b => b > 0));
            Assert.AreEqual(1, newPsm.matchedIonsList[ProductType.Y].Count(b => b > 0));
            Assert.IsTrue(newPsm.LocalizedScores[1] > 3 && newPsm.LocalizedScores[1] < 4);


        }

        #endregion Public Methods

    }
}