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
using Omics;
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

            List<IBioPolymerWithSetMods> allPeptidesWithSetModifications = parentProteinForMatch.Digest(commonParameters.DigestionParams, new List<Modification>(), variableModifications).ToList();
            Assert.That(allPeptidesWithSetModifications.Count(), Is.EqualTo(4));
            var ps = allPeptidesWithSetModifications.First();

            var pepWithSetModsForSpectrum = allPeptidesWithSetModifications[1];
            MsDataFile myMsDataFile = new TestDataFile(new List<IBioPolymerWithSetMods> { pepWithSetModsForSpectrum });
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

            // the evidence holds one candidate per residue, in order; the termini are not separate candidates
            LocalizationEvidence evidence = newPsm.LocalizationEvidence;
            Assert.That(evidence, Is.Not.Null);
            Assert.That(evidence.CandidateCount, Is.EqualTo(ps.Length));
            for (int r = 0; r < ps.Length; r++)
            {
                Assert.That(evidence.Placements[r].Sites, Is.EqualTo(new[] { new LocalizationSite(r, ps.BaseSequence[r]) }));
                Assert.That(evidence.Scores[r], Is.EqualTo(newPsm.LocalizedScores[r]));
            }
            Assert.That(evidence.MassDifference, Is.EqualTo(21.981943).Within(1e-3));
            Assert.That(evidence.IsProteinNTerminal, Is.True);

            // the E placement matches a peak the M placement does not, and every relative intensity is in (0, 1]
            Assert.That(Enumerable.Range(0, evidence.PeakCount).Any(p => evidence.Matches(1, p) && !evidence.Matches(0, p)));
            Assert.That(Enumerable.Range(0, evidence.PeakCount).All(p => evidence.PeakRelativeIntensity(p) > 0 && evidence.PeakRelativeIntensity(p) <= 1));
        }

        [Test]
        public static void TerminalMassFragmentsAsTheSameMassOnTheTerminalResidue()
        {
            // the localization engine treats a peptide-terminal placement and its terminal residue as one site; that is
            // only right if the two fragment to the same masses
            var peptide = (PeptideWithSetModifications)new Protein("PEPTIDEK", "P1")
                .Digest(new DigestionParams(), new List<Modification>(), new List<Modification>()).First();
            const double mass = 42.010565;

            foreach (var (residueIndex, terminusKey) in new[] { (0, 1), (peptide.Length - 1, peptide.Length + 2) })
            {
                var onResidue = peptide.Localize(residueIndex, mass);
                var onTerminus = new PeptideWithSetModifications(peptide.Protein, peptide.DigestionParams,
                    peptide.OneBasedStartResidueInProtein, peptide.OneBasedEndResidueInProtein, peptide.CleavageSpecificityForFdrCategory,
                    peptide.PeptideDescription, peptide.MissedCleavages,
                    new Dictionary<int, Modification> { { terminusKey, new Modification(_locationRestriction: "Anywhere.", _monoisotopicMass: mass) } },
                    peptide.NumFixedMods);

                foreach (DissociationType dissociation in new[] { DissociationType.HCD, DissociationType.ETD })
                {
                    var residueProducts = new List<Product>();
                    var terminusProducts = new List<Product>();
                    onResidue.Fragment(dissociation, FragmentationTerminus.Both, residueProducts);
                    onTerminus.Fragment(dissociation, FragmentationTerminus.Both, terminusProducts);

                    Assert.That(terminusProducts.Select(p => (p.ProductType, p.FragmentNumber, System.Math.Round(p.NeutralMass, 6))),
                        Is.EquivalentTo(residueProducts.Select(p => (p.ProductType, p.FragmentNumber, System.Math.Round(p.NeutralMass, 6)))));
                }
            }
        }

        [Test]
        public static void EvidenceBuilderCountsAPeakOncePerCandidateAndAcrossCandidates()
        {
            var product = new Product(ProductType.b, FragmentationTerminus.N, 100, 1, 1, 0);
            var otherProduct = new Product(ProductType.y, FragmentationTerminus.C, 100, 1, 1, 0);
            var builder = new LocalizationEvidenceBuilder(basePeakIntensity: 200);

            // 70 distinct peaks, so the match masks run past one 64-bit word
            var manyIons = Enumerable.Range(0, 70).Select(i => new MatchedFragmentIon(product, 100 + i, 10 + i, 1)).ToList();
            builder.AddCandidate(new LocalizationPlacement(new[] { new LocalizationSite(0, 'P') }), 5, manyIons);

            // two products matching one peak are one peak; the second candidate matches peak 0 and peak 69 only
            builder.AddCandidate(new LocalizationPlacement(new[] { new LocalizationSite(1, 'E') }), 2, new List<MatchedFragmentIon>
            {
                new MatchedFragmentIon(product, 100, 10, 1),
                new MatchedFragmentIon(otherProduct, 100, 10, 1),
                new MatchedFragmentIon(product, 169, 79, 1),
            });

            LocalizationEvidence evidence = builder.Build(massDifference: 1, peptideLength: 2, isProteinNTerminal: false, isProteinCTerminal: false);

            Assert.That(evidence.PeakCount, Is.EqualTo(70));
            Assert.That(evidence.CandidateCount, Is.EqualTo(2));
            Assert.That(evidence.PeakRelativeIntensity(0), Is.EqualTo(0.05).Within(1e-6));
            Assert.That(evidence.PeakRelativeIntensity(69), Is.EqualTo(79 / 200.0).Within(1e-6));
            Assert.That(Enumerable.Range(0, 70).All(p => evidence.Matches(0, p)));
            Assert.That(Enumerable.Range(0, 70).Where(p => evidence.Matches(1, p)), Is.EqualTo(new[] { 0, 69 }));
            Assert.That(evidence.Scores, Is.EqualTo(new[] { 5.0, 2.0 }));
            Assert.That(evidence.IsSharedByAllCandidates(0), Is.True);
            Assert.That(evidence.IsSharedByAllCandidates(1), Is.False);
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