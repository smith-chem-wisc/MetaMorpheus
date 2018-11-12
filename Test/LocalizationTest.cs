using Chemistry;
using EngineLayer;
using EngineLayer.Localization;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using Proteomics.Fragmentation;
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
            var peps = prot.Digest(digestionParams, new List<Modification>(), new List<Modification>()).ToList();

             Assert.AreEqual(1 + 2 + 3 + 4 + 5 + 6 + 7 + 8 + 9, peps.Count);
        }

        [Test]
        public static void TestLocalization()
        {
            Protein parentProteinForMatch = new Protein("MEK", null);
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 1);
            ModificationMotif.TryGetMotif("E", out ModificationMotif motif);
            List<Modification> variableModifications = new List<Modification> { new Modification(_originalId: "21", _target: motif, _locationRestriction: "Anywhere.", _monoisotopicMass: 21.981943) };

            List<PeptideWithSetModifications> allPeptidesWithSetModifications = parentProteinForMatch.Digest(digestionParams, new List<Modification>(), variableModifications).ToList();
            Assert.AreEqual(4, allPeptidesWithSetModifications.Count());
            PeptideWithSetModifications ps = allPeptidesWithSetModifications.First();

            PeptideWithSetModifications pepWithSetModsForSpectrum = allPeptidesWithSetModifications[1];
            MsDataFile myMsDataFile = new TestDataFile(new List<PeptideWithSetModifications> { pepWithSetModsForSpectrum });
            Tolerance fragmentTolerance = new AbsoluteTolerance(0.01);

            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(myMsDataFile.GetAllScansList().Last(), pepWithSetModsForSpectrum.MonoisotopicMass.ToMz(1), 1, null, new CommonParameters());

            var theoreticalProducts = ps.Fragment(DissociationType.HCD, FragmentationTerminus.Both).ToList();
            var matchedIons = MetaMorpheusEngine.MatchFragmentIons(scan, theoreticalProducts, new CommonParameters());
            PeptideSpectralMatch newPsm = new PeptideSpectralMatch(ps, 0, 0, 2, scan, digestionParams, matchedIons);
            newPsm.ResolveAllAmbiguities();

            CommonParameters commonParameters = new CommonParameters(productMassTolerance: fragmentTolerance);

            LocalizationEngine f = new LocalizationEngine(new List<PeptideSpectralMatch> { newPsm }, myMsDataFile, commonParameters, new List<string>());
            f.Run();

            // single peak matches
            Assert.AreEqual(1, newPsm.MatchedFragmentIons.Where(p => p.NeutralTheoreticalProduct.ProductType == ProductType.b).Count());//including b1 now
            Assert.AreEqual(1, newPsm.MatchedFragmentIons.Where(p => p.NeutralTheoreticalProduct.ProductType == ProductType.y).Count());

            // when localizing, three peaks match
            Assert.IsTrue(newPsm.LocalizedScores[1] > 4 && newPsm.LocalizedScores[1] < 5);//we have another matched ion
        }

        [Test]
        public static void TestSeparateIonsByTerminus()
        {
            List<ProductType> allIonTypes = new List<ProductType> { ProductType.b, ProductType.c, ProductType.zPlusOne, ProductType.y };
            List<List<ProductType>> separated = ProductTypeMethods.SeparateIonsByTerminus(allIonTypes);
            Assert.IsTrue(separated.Count == 2);
            Assert.IsTrue(separated[0].Count == 2);
            Assert.IsTrue(separated[1].Count == 2);
            Assert.IsTrue(separated[0].Contains(ProductType.b));
            Assert.IsTrue(separated[0].Contains(ProductType.c));
            Assert.IsTrue(separated[1].Contains(ProductType.y));
            Assert.IsTrue(separated[1].Contains(ProductType.zPlusOne));
            List<List<ProductType>> nOnly = ProductTypeMethods.SeparateIonsByTerminus(separated[0]);
            Assert.IsTrue(nOnly.Count == 1);
            Assert.IsTrue(nOnly[0].Count == 2);
            Assert.IsTrue(nOnly[0].Contains(ProductType.b));
            Assert.IsTrue(nOnly[0].Contains(ProductType.c));
            List<List<ProductType>> cOnly = ProductTypeMethods.SeparateIonsByTerminus(separated[1]);
            Assert.IsTrue(cOnly.Count == 1);
            Assert.IsTrue(cOnly[0].Count == 2);
            Assert.IsTrue(cOnly[0].Contains(ProductType.y));
            Assert.IsTrue(cOnly[0].Contains(ProductType.zPlusOne));
        }
    }
}