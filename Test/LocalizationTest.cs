using Chemistry;
using EngineLayer;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using System.Collections.Generic;
using System.Linq;
using TaskLayer;

namespace Test
{
    [TestFixture]
    public class LocalizationTest
    {

        #region Public Methods

        [Test]
        public static void TestNonSpecific()
        {
            Protease p = GlobalTaskLevelSettings.ProteaseDictionary["non-specific"];
            Protein prot = new Protein("MABCDEFGH", null);

            Assert.AreEqual(1 + 2 + 3 + 4 + 5 + 6 + 7 + 8, prot.Digest(p, 8, 1, 9, InitiatorMethionineBehavior.Retain, new List<ModificationWithMass>(), false).Count());
        }

        [Test]
        public static void TestLocalization()
        {
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);

            Protein parentProteinForMatch = new Protein("MEK", null);
            PeptideWithPossibleModifications pwpm = parentProteinForMatch.Digest(protease, 0, null, null, InitiatorMethionineBehavior.Variable, new List<ModificationWithMass>(), false).First();
            ModificationMotif motif;
            ModificationMotif.TryGetMotif("E", out motif);
            List<ModificationWithMass> variableModifications = new List<ModificationWithMass> { new ModificationWithMass("21", null, motif, TerminusLocalization.Any, 21.981943, null, new List<double> { 0 }, new List<double> { 21.981943 }, null) };

            List<PeptideWithSetModifications> allPeptidesWithSetModifications = pwpm.GetPeptidesWithSetModifications(variableModifications, 2, 1).ToList();
            Assert.AreEqual(2, allPeptidesWithSetModifications.Count());
            PeptideWithSetModifications ps = allPeptidesWithSetModifications.First();

            List<ProductType> lp = new List<ProductType> { ProductType.B, ProductType.Y };

            PeptideWithSetModifications pepWithSetModsForSpectrum = allPeptidesWithSetModifications.Last();
            IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile = new TestDataFile(new List<PeptideWithSetModifications> { pepWithSetModsForSpectrum });
            Tolerance fragmentTolerance = new AbsoluteTolerance(0.01);

            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(myMsDataFile.Last() as IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>>, new MzPeak(pepWithSetModsForSpectrum.MonoisotopicMass.ToMz(1), 1), 1, null);
            Psm newPsm = new Psm(ps.CompactPeptide, 0, 0, 2, scan);

            Assert.IsNull(newPsm.MostProbableProteinInfo);

            Dictionary<ModificationWithMass, ushort> modsDictionary = new Dictionary<ModificationWithMass, ushort>();
            Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> matching = new Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>>
            {
                {ps.CompactPeptide, new HashSet<PeptideWithSetModifications>{ ps} }
            };

            newPsm.MatchToProteinLinkedPeptides(matching);

            LocalizationEngine f = new LocalizationEngine(new List<Psm> { newPsm }, lp, myMsDataFile, fragmentTolerance, null, false);
            f.Run();

            // Was single peak!!!
            Assert.AreEqual(0, newPsm.LocalizationResults.MatchedIonMassesListPositiveIsMatch[ProductType.B].Count(b => b > 0));
            Assert.AreEqual(1, newPsm.LocalizationResults.MatchedIonMassesListPositiveIsMatch[ProductType.Y].Count(b => b > 0));
            // If localizing, three match!!!
            Assert.IsTrue(newPsm.LocalizationResults.LocalizedScores[1] > 3 && newPsm.LocalizationResults.LocalizedScores[1] < 4);
        }

        #endregion Public Methods

    }
}