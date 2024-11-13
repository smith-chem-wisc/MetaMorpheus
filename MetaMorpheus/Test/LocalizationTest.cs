using Chemistry;
using EngineLayer;
using EngineLayer.Localization;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.Linq;
using Omics.Modifications;

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

            Assert.That(peps.Count, Is.EqualTo(1 + 2 + 3 + 4 + 5 + 6 + 7 + 8 + 9));
        }

        [Test]
        public static void TestLocalization()
        {
            Protein parentProteinForMatch = new Protein("MEK", null);
            CommonParameters commonParameters = new CommonParameters(digestionParams: new DigestionParams(minPeptideLength: 1));
            var fsp = new List<(string fileName, CommonParameters fileSpecificParameters)>();
            fsp.Add(("", commonParameters));
            ModificationMotif.TryGetMotif("E", out ModificationMotif motif);
            List<Modification> variableModifications = new List<Modification> { new Modification(_originalId: "21", _target: motif, _locationRestriction: "Anywhere.", _monoisotopicMass: 21.981943) };

            List<PeptideWithSetModifications> allPeptidesWithSetModifications = parentProteinForMatch.Digest(commonParameters.DigestionParams, new List<Modification>(), variableModifications).ToList();
            Assert.That(allPeptidesWithSetModifications.Count(), Is.EqualTo(4));
            PeptideWithSetModifications ps = allPeptidesWithSetModifications.First();

            PeptideWithSetModifications pepWithSetModsForSpectrum = allPeptidesWithSetModifications[1];
            MsDataFile myMsDataFile = new TestDataFile(new List<PeptideWithSetModifications> { pepWithSetModsForSpectrum });
            Tolerance fragmentTolerance = new AbsoluteTolerance(0.01);

            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(myMsDataFile.GetAllScansList().Last(), pepWithSetModsForSpectrum.MonoisotopicMass.ToMz(1), 1, null, new CommonParameters());

            var theoreticalProducts = new List<Product>();
            ps.Fragment(DissociationType.HCD, FragmentationTerminus.Both, theoreticalProducts);
            var matchedIons = MetaMorpheusEngine.MatchFragmentIons(scan, theoreticalProducts, new CommonParameters());
            SpectralMatch newPsm = new PeptideSpectralMatch(ps, 0, 0, 2, scan, commonParameters, matchedIons);
            newPsm.ResolveAllAmbiguities();

            LocalizationEngine f = new LocalizationEngine(new List<SpectralMatch> { newPsm }, myMsDataFile, commonParameters, fsp, new List<string>());
            f.Run();

            // single peak matches
            Assert.That(newPsm.MatchedFragmentIons.Where(p => p.NeutralTheoreticalProduct.ProductType == ProductType.b).Count(), Is.EqualTo(1));//including b1 now
            Assert.That(newPsm.MatchedFragmentIons.Where(p => p.NeutralTheoreticalProduct.ProductType == ProductType.y).Count(), Is.EqualTo(1));

            // when localizing, three peaks match
            Assert.That(newPsm.LocalizedScores[1] > 4 && newPsm.LocalizedScores[1] < 5);//we have another matched ion
        }

        [Test]
        public static void TestSeparateIonsByTerminus()
        {
            List<ProductType> allIonTypes = new List<ProductType> { ProductType.b, ProductType.c, ProductType.zPlusOne, ProductType.y };
            List<List<ProductType>> separated = ProductTypeMethods.SeparateIonsByTerminus(allIonTypes);
            Assert.That(separated.Count == 2);
            Assert.That(separated[0].Count == 2);
            Assert.That(separated[1].Count == 2);
            Assert.That(separated[0].Contains(ProductType.b));
            Assert.That(separated[0].Contains(ProductType.c));
            Assert.That(separated[1].Contains(ProductType.y));
            Assert.That(separated[1].Contains(ProductType.zPlusOne));
            List<List<ProductType>> nOnly = ProductTypeMethods.SeparateIonsByTerminus(separated[0]);
            Assert.That(nOnly.Count == 1);
            Assert.That(nOnly[0].Count == 2);
            Assert.That(nOnly[0].Contains(ProductType.b));
            Assert.That(nOnly[0].Contains(ProductType.c));
            List<List<ProductType>> cOnly = ProductTypeMethods.SeparateIonsByTerminus(separated[1]);
            Assert.That(cOnly.Count == 1);
            Assert.That(cOnly[0].Count == 2);
            Assert.That(cOnly[0].Contains(ProductType.y));
            Assert.That(cOnly[0].Contains(ProductType.zPlusOne));
        }
    }
}