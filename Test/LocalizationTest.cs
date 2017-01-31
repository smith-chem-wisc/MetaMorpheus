using EngineLayer;
using EngineLayer.ClassicSearch;
using EngineLayer.Gptmd;
using MassSpectrometry;
using NUnit.Framework;

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
            Dictionary<int, List<MetaMorpheusModification>> oneBasedPossibleLocalizedModifications = new Dictionary<int, List<MetaMorpheusModification>>();
            string accession = null;
            int[] beginPositions = null;
            int[] endPositions = null;
            string[] bigPeptideTypes = null;
            string name = null;
            string fullName = null;
            Protein parentProteinForMatch = new Protein("MEK", accession, oneBasedPossibleLocalizedModifications, beginPositions, endPositions, bigPeptideTypes, name, fullName, 0, false, false);
            PeptideWithPossibleModifications pwpm = new PeptideWithPossibleModifications(1, 3, parentProteinForMatch, 0, null);

            List<MetaMorpheusModification> variableModifications = new List<MetaMorpheusModification> { new MetaMorpheusModification(null, ModificationType.AminoAcidResidue, 'E', null, '\0', 21.981943, 21.981943, 21.981943, new Chemistry.ChemicalFormula("Na-1 H1")) };

            IEnumerable<MetaMorpheusModification> allKnownFixedModifications = new List<MetaMorpheusModification>();

            List<PeptideWithSetModifications> allPeptidesWithSetModifications = pwpm.GetPeptideWithSetModifications(variableModifications, 2, 1, allKnownFixedModifications).ToList();
            Assert.AreEqual(2, allPeptidesWithSetModifications.Count());
            PeptideWithSetModifications ps = allPeptidesWithSetModifications.First();

            List<ProductType> lp = new List<ProductType> { ProductType.B, ProductType.Y };

            PeptideWithSetModifications pepWithSetModsForSpectrum = allPeptidesWithSetModifications.Last();
            IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile = new TestDataFile(new List<PeptideWithSetModifications> { pepWithSetModsForSpectrum });
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
            PsmParent newPsm = new PsmClassic(ps, fileName, scanRetentionTime, scanPrecursorIntensity, scanPrecursorMass, 2, scanPrecursorCharge, scanExperimentalPeaks, totalIonCurrent, scanPrecursorMZ, score);

            Assert.IsNull(newPsm.LocalizedScores);
            Assert.IsNull(newPsm.matchedIonsList);
            var cool = new PsmWithMultiplePossiblePeptides(newPsm, peptidesWithSetModifications, fragmentTolerance, myMsDataFile, lp);

            // Was single peak, now three match!!!
            Assert.AreEqual(0, newPsm.matchedIonsList[ProductType.B].Count(b => b > 0));
            Assert.AreEqual(1, newPsm.matchedIonsList[ProductType.Y].Count(b => b > 0));
            Assert.IsTrue(newPsm.LocalizedScores[1] > 3 && newPsm.LocalizedScores[1] < 4);


        }

        #endregion Public Methods

    }
}