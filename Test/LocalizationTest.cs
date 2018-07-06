using Chemistry;
using EngineLayer;
using EngineLayer.Localization;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    [TestFixture]
    public static class LocalizationTest
    {
        [Test]
        public static void TestNonSpecific()
        {
            Protease p = ProteaseDictionary.Dictionary["non-specific"];
            Protein prot = new Protein("MABCDEFGH", null);

            DigestionParams digestionParams = new DigestionParams(protease: p.Name, maxMissedCleavages: 8, minPeptideLength: 1, maxPeptideLength: 9, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);

            Assert.AreEqual(1 + 2 + 3 + 4 + 5 + 6 + 7 + 8, prot.Digest(digestionParams, new List<ModificationWithMass>(), new List<ModificationWithMass>()).Count());
        }

        [Test]
        public static void TestLocalization()
        {
            Protein parentProteinForMatch = new Protein("MEK", null);
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 1);
            ModificationMotif.TryGetMotif("E", out ModificationMotif motif);
            List<ModificationWithMass> variableModifications = new List<ModificationWithMass> { new ModificationWithMass("21", null, motif, TerminusLocalization.Any, 21.981943) };

            List<PeptideWithSetModifications> allPeptidesWithSetModifications = parentProteinForMatch.Digest(digestionParams, new List<ModificationWithMass>(), variableModifications).ToList();
            Assert.AreEqual(4, allPeptidesWithSetModifications.Count);
            PeptideWithSetModifications ps = allPeptidesWithSetModifications.First();

            List<ProductType> productTypes = new List<ProductType> { ProductType.BnoB1ions, ProductType.Y };

            PeptideWithSetModifications pepWithSetModsForSpectrum = allPeptidesWithSetModifications[1];
            MsDataFile myMsDataFile = new TestDataFile(new List<PeptideWithSetModifications> { pepWithSetModsForSpectrum });
            Tolerance fragmentTolerance = new AbsoluteTolerance(0.01);

            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(myMsDataFile.GetAllScansList().Last(), pepWithSetModsForSpectrum.MonoisotopicMass.ToMz(1), 1, null);
            PeptideSpectralMatch newPsm = new PeptideSpectralMatch(ps.CompactPeptide(TerminusType.None), 0, 0, 2, scan, digestionParams);

            Dictionary<ModificationWithMass, ushort> modsDictionary = new Dictionary<ModificationWithMass, ushort>();
            Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> matching = new Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>
            {
                {ps.CompactPeptide(TerminusType.None), new HashSet<PeptideWithSetModifications>{ ps} }
            };

            newPsm.MatchToProteinLinkedPeptides(matching);

            CommonParameters commonParameters = new CommonParameters(productMassTolerance: fragmentTolerance);
            LocalizationEngine f = new LocalizationEngine(new List<PeptideSpectralMatch> { newPsm }, productTypes, myMsDataFile, commonParameters, new List<string>());
            f.Run();

            // Was single peak!!!
            Assert.AreEqual(0, newPsm.MatchedIonMassToChargeRatioDict[ProductType.BnoB1ions].Count(b => b > 0));
            Assert.AreEqual(1, newPsm.MatchedIonMassToChargeRatioDict[ProductType.Y].Count(b => b > 0));
            Assert.AreEqual(0, newPsm.MatchedIonIntensitiesDict[ProductType.BnoB1ions].Count(b => b > 0));
            Assert.AreEqual(1, newPsm.MatchedIonIntensitiesDict[ProductType.Y].Count(b => b > 0));
            // If localizing, three match!!!
            Assert.IsTrue(newPsm.LocalizedScores[1] > 3 && newPsm.LocalizedScores[1] < 4);
        }

        [Test]
        public static void TestSeparateIonsByTerminus()
        {
            List<ProductType> allIonTypes = new List<ProductType> { ProductType.B, ProductType.C, ProductType.Zdot, ProductType.Y };
            List<List<ProductType>> separated = ProductTypeMethods.SeparateIonsByTerminus(allIonTypes);
            Assert.IsTrue(separated.Count == 2);
            Assert.IsTrue(separated[0].Count == 2);
            Assert.IsTrue(separated[1].Count == 2);
            Assert.IsTrue(separated[0].Contains(ProductType.B));
            Assert.IsTrue(separated[0].Contains(ProductType.C));
            Assert.IsTrue(separated[1].Contains(ProductType.Y));
            Assert.IsTrue(separated[1].Contains(ProductType.Zdot));
            List<List<ProductType>> nOnly = ProductTypeMethods.SeparateIonsByTerminus(separated[0]);
            Assert.IsTrue(nOnly.Count == 1);
            Assert.IsTrue(nOnly[0].Count == 2);
            Assert.IsTrue(nOnly[0].Contains(ProductType.B));
            Assert.IsTrue(nOnly[0].Contains(ProductType.C));
            List<List<ProductType>> cOnly = ProductTypeMethods.SeparateIonsByTerminus(separated[1]);
            Assert.IsTrue(cOnly.Count == 1);
            Assert.IsTrue(cOnly[0].Count == 2);
            Assert.IsTrue(cOnly[0].Contains(ProductType.Y));
            Assert.IsTrue(cOnly[0].Contains(ProductType.Zdot));
        }
    }
}