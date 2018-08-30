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

            Assert.AreEqual(1 + 2 + 3 + 4 + 5 + 6 + 7 + 8, prot.Digest(digestionParams, new List<Modification>(), new List<Modification>()).Count());
        }

        [Test]
        public static void TestLocalization()
        {
            Protein parentProteinForMatch = new Protein("MEK", null);
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 1);
            ModificationMotif.TryGetMotif("E", out ModificationMotif motif);
            List<Modification> variableModifications = new List<Modification> { new Modification(_id: "21", _target: motif, _locationRestriction: "Anywhere.", _monoisotopicMass: 21.981943) };

            List<PeptideWithSetModifications> allPeptidesWithSetModifications = parentProteinForMatch.Digest(digestionParams, new List<Modification>(), variableModifications).ToList();
            Assert.AreEqual(4, allPeptidesWithSetModifications.Count());
            PeptideWithSetModifications ps = allPeptidesWithSetModifications.First();

            List<ProductType> lp = new List<ProductType> { ProductType.B, ProductType.Y };

            PeptideWithSetModifications pepWithSetModsForSpectrum = allPeptidesWithSetModifications[1];
            MsDataFile myMsDataFile = new TestDataFile(new List<PeptideWithSetModifications> { pepWithSetModsForSpectrum });
            Tolerance fragmentTolerance = new AbsoluteTolerance(0.01);

            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(myMsDataFile.GetAllScansList().Last(), pepWithSetModsForSpectrum.MonoisotopicMass.ToMz(1), 1, null);

            var theoreticalProducts = ps.Fragment(DissociationType.HCD, FragmentationTerminus.Both).ToList();
            var matchedIons = MetaMorpheusEngine.MatchFragmentIons(scan.TheScan.MassSpectrum, theoreticalProducts, new CommonParameters());
            PeptideSpectralMatch newPsm = new PeptideSpectralMatch(ps, 0, 0, 2, scan, digestionParams, matchedIons);
            newPsm.ResolveAllAmbiguities();

            CommonParameters commonParameters = new CommonParameters(productMassTolerance: fragmentTolerance);

            ////DEBUG CODE
            //List<PeptideSpectralMatch> l = new List<PeptideSpectralMatch>() { newPsm };
            //PeptideSpectralMatch[] a = l.ToArray();

            //List<Ms2ScanWithSpecificMass> m = new List<Ms2ScanWithSpecificMass> { scan };
            //Ms2ScanWithSpecificMass[] marray = m.ToArray();

            //EngineLayer.ClassicSearch.ClassicSearchEngine c = new EngineLayer.ClassicSearch.ClassicSearchEngine(a, marray, null, null, new List<Protein> { parentProteinForMatch }, lp, null, commonParameters, new List<string>());
            //c.Run();

            //int j = newPsm.MatchedFragmentIons.Count();

            ////END DEBUG CODE

            LocalizationEngine f = new LocalizationEngine(new List<PeptideSpectralMatch> { newPsm }, lp, myMsDataFile, commonParameters, new List<string>());
            f.Run();

            // single peak matches
            Assert.AreEqual(0, newPsm.MatchedFragmentIons.Where(p => p.NeutralTheoreticalProduct.ProductType == ProductType.B).Count());
            Assert.AreEqual(1, newPsm.MatchedFragmentIons.Where(p => p.NeutralTheoreticalProduct.ProductType == ProductType.Y).Count());

            // when localizing, three peaks match
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