using Chemistry;
using EngineLayer;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    [TestFixture]
    public static class LocalizationTest
    {
        #region Public Methods

        [Test]
        public static void TestNonSpecific()
        {
            Protease p = GlobalEngineLevelSettings.ProteaseDictionary["non-specific"];
            Protein prot = new Protein("MABCDEFGH", null);

            DigestionParams digestionParams = new DigestionParams
            {
                InitiatorMethionineBehavior = InitiatorMethionineBehavior.Retain,
                MaxMissedCleavages = 8,
                MinPeptideLength = 1,
                MaxPeptideLength = 9,
                Protease = p
            };

            Assert.AreEqual(1 + 2 + 3 + 4 + 5 + 6 + 7 + 8, prot.Digest(digestionParams, new List<ModificationWithMass>()).Count());
        }

        [Test]
        public static void TestLocalization()
        {
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);

            Protein parentProteinForMatch = new Protein("MEK", null);
            DigestionParams digestionParams = new DigestionParams
            {
                MinPeptideLength = 1,
            };
            PeptideWithPossibleModifications pwpm = parentProteinForMatch.Digest(digestionParams, new List<ModificationWithMass>()).First();
            ModificationMotif.TryGetMotif("E", out ModificationMotif motif);
            List<ModificationWithMass> variableModifications = new List<ModificationWithMass> { new ModificationWithMass("21", null, motif, TerminusLocalization.Any, 21.981943) };

            List<PeptideWithSetModifications> allPeptidesWithSetModifications = pwpm.GetPeptidesWithSetModifications(digestionParams, variableModifications).ToList();
            Assert.AreEqual(2, allPeptidesWithSetModifications.Count());
            PeptideWithSetModifications ps = allPeptidesWithSetModifications.First();

            List<ProductType> lp = new List<ProductType> { ProductType.BnoB1ions, ProductType.Y };

            PeptideWithSetModifications pepWithSetModsForSpectrum = allPeptidesWithSetModifications.Last();
            IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile = new TestDataFile(new List<PeptideWithSetModifications> { pepWithSetModsForSpectrum });
            Tolerance fragmentTolerance = new AbsoluteTolerance(0.01);

            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(myMsDataFile.Last() as IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>>, new MzPeak(pepWithSetModsForSpectrum.MonoisotopicMass.ToMz(1), 1), 1, null);
            Psm newPsm = new Psm(ps.CompactPeptide(TerminusType.None), 0, 0, 2, scan);

            Assert.IsNull(newPsm.MostProbableProteinInfo);

            Dictionary<ModificationWithMass, ushort> modsDictionary = new Dictionary<ModificationWithMass, ushort>();
            Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> matching = new Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>
            {
                {ps.CompactPeptide(TerminusType.None), new HashSet<PeptideWithSetModifications>{ ps} }
            };

            newPsm.MatchToProteinLinkedPeptides(matching);

            LocalizationEngine f = new LocalizationEngine(new List<Psm> { newPsm }, lp, myMsDataFile, fragmentTolerance, new List<string>(), false);
            f.Run();

            // Was single peak!!!
            Assert.AreEqual(0, newPsm.MatchedIonDictOnlyMatches[ProductType.BnoB1ions].Count(b => b > 0));
            Assert.AreEqual(1, newPsm.MatchedIonDictOnlyMatches[ProductType.Y].Count(b => b > 0));
            // If localizing, three match!!!
            Assert.IsTrue(newPsm.LocalizedScores[1] > 3 && newPsm.LocalizedScores[1] < 4);
        }

        [Test]
        public static void TestSeparateIonsByTerminus()
        {
            List<ProductType> allIonTypes = new List<ProductType> { ProductType.B, ProductType.C, ProductType.Zdot, ProductType.Y };
            List<List<ProductType>> separated = ProductTypeMethod.SeparateIonsByTerminus(allIonTypes);
            Assert.IsTrue(separated.Count == 2);
            Assert.IsTrue(separated[0].Count == 2);
            Assert.IsTrue(separated[1].Count == 2);
            Assert.IsTrue(separated[0].Contains(ProductType.B));
            Assert.IsTrue(separated[0].Contains(ProductType.C));
            Assert.IsTrue(separated[1].Contains(ProductType.Y));
            Assert.IsTrue(separated[1].Contains(ProductType.Zdot));
            List<List<ProductType>> nOnly = ProductTypeMethod.SeparateIonsByTerminus(separated[0]);
            Assert.IsTrue(nOnly.Count == 1);
            Assert.IsTrue(nOnly[0].Count == 2);
            Assert.IsTrue(nOnly[0].Contains(ProductType.B));
            Assert.IsTrue(nOnly[0].Contains(ProductType.C));
            List<List<ProductType>> cOnly = ProductTypeMethod.SeparateIonsByTerminus(separated[1]);
            Assert.IsTrue(cOnly.Count == 1);
            Assert.IsTrue(cOnly[0].Count == 2);
            Assert.IsTrue(cOnly[0].Contains(ProductType.Y));
            Assert.IsTrue(cOnly[0].Contains(ProductType.Zdot));
        }

        #endregion Public Methods
    }
}