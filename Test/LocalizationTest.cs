using Chemistry;
using EngineLayer;
using EngineLayer.ClassicSearch;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
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
            Dictionary<int, List<Modification>> oneBasedPossibleLocalizedModifications = new Dictionary<int, List<Modification>>();
            string accession = null;
            int?[] beginPositions = null;
            int?[] endPositions = null;
            string[] bigPeptideTypes = null;
            string name = null;
            string fullName = null;
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);

            Protein parentProteinForMatch = new Protein("MEK", accession, null, oneBasedPossibleLocalizedModifications, beginPositions, endPositions, bigPeptideTypes, name, fullName, false, false, null);
            PeptideWithPossibleModifications pwpm = parentProteinForMatch.Digest(protease, 0, null, null, InitiatorMethionineBehavior.Variable, new List<ModificationWithMass>()).First();
            ModificationMotif motif;
            ModificationMotif.TryGetMotif("E", out motif);
            List<ModificationWithMass> variableModifications = new List<ModificationWithMass> { new ModificationWithMass("21", null, motif, ModificationSites.Any, 21.981943, null, new List<double> { 0 }, new List<double> { 21.981943 }, null) };

            List<PeptideWithSetModifications> allPeptidesWithSetModifications = pwpm.GetPeptidesWithSetModifications(variableModifications, 2, 1).ToList();
            Assert.AreEqual(2, allPeptidesWithSetModifications.Count());
            PeptideWithSetModifications ps = allPeptidesWithSetModifications.First();

            List<ProductType> lp = new List<ProductType> { ProductType.B, ProductType.Y };

            PeptideWithSetModifications pepWithSetModsForSpectrum = allPeptidesWithSetModifications.Last();
            IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile = new TestDataFile(new List<PeptideWithSetModifications> { pepWithSetModsForSpectrum });
            Tolerance fragmentTolerance = new Tolerance(ToleranceUnit.Absolute, 0.01);

            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(myMsDataFile.Last() as IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>>, new MzPeak(pepWithSetModsForSpectrum.MonoisotopicMass.ToMz(1), 1), 1, null);
            PsmParent newPsm = new PsmClassic(ps, 0, 0, 2, scan);

            Assert.IsNull(newPsm.Pli);

            Dictionary<ModificationWithMass, ushort> modsDictionary = new Dictionary<ModificationWithMass, ushort>();
            Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> matching = new Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>>
            {
                {newPsm.GetCompactPeptide(modsDictionary), new HashSet<PeptideWithSetModifications>{ ps} }
            };

            newPsm.ComputeProteinLevelInfo(matching, fragmentTolerance, scan, lp, modsDictionary);
            // Was single peak!!!
            Assert.AreEqual(0, newPsm.Pli.MatchedIonsListPositiveIsMatch[ProductType.B].Count(b => b > 0));
            Assert.AreEqual(1, newPsm.Pli.MatchedIonsListPositiveIsMatch[ProductType.Y].Count(b => b > 0));
            // If localizing, three match!!!
            Assert.IsTrue(newPsm.Pli.LocalizedScores[1] > 3 && newPsm.Pli.LocalizedScores[1] < 4);
        }

        #endregion Public Methods

    }
}