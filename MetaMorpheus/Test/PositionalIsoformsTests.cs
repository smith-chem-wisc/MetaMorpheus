using Chemistry;
using EngineLayer;
using EngineLayer.Localization;
using MassSpectrometry;
using NUnit.Framework;
using Omics;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    [TestFixture]
    public static class PositionalIsoformsTests
    {
        private static readonly Modification PhosphoS = Make("Phosphorylation", "S", "H1O3P1");
        private static readonly Modification PhosphoT = Make("Phosphorylation", "T", "H1O3P1");
        private static readonly Modification OxidationM = Make("Oxidation", "M", "O1");
        private static readonly Modification AcetylNTerm = Make("Acetylation", "X", "C2H2O1", "N-terminal.");

        private static Modification Make(string id, string residue, string formula, string location = "Anywhere.")
        {
            ModificationMotif.TryGetMotif(residue, out ModificationMotif motif);
            return new Modification(_originalId: id, _modificationType: "Common Biological", _target: motif,
                _locationRestriction: location, _chemicalFormula: ChemicalFormula.ParseFormula(formula));
        }

        [Test]
        public static void ModificationsAnnotatedInTheDatabaseArePlaceable()
        {
            // a GPTMD database annotates phosphotyrosine on one protein and phosphothreonine on another, and the
            // search varies only phosphoserine and oxidation
            ModificationMotif.TryGetMotif("Y", out ModificationMotif y);
            var phosphoY = new Modification(_originalId: "Phosphorylation", _modificationType: "Common Biological", _target: y,
                _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H1O3P1"));
            var proteins = new List<IBioPolymer>
            {
                new Protein("PEPSYK", "P1", oneBasedModifications: new Dictionary<int, List<Modification>> { { 5, new List<Modification> { phosphoY } } }),
                new Protein("PETYK", "P2", oneBasedModifications: new Dictionary<int, List<Modification>> { { 3, new List<Modification> { PhosphoT } } }),
                new Protein("KYSPEP", "DECOY_P1", isDecoy: true, oneBasedModifications: new Dictionary<int, List<Modification>> { { 2, new List<Modification> { phosphoY } } }),
            };

            var placeable = TaskLayer.PostSearchAnalysisTask.PlaceableModifications(new[] { PhosphoS, OxidationM }, proteins);

            // the variable ones and the annotated ones, each once
            Assert.That(placeable, Is.EquivalentTo(new[] { PhosphoS, OxidationM, phosphoY, PhosphoT }));
        }

        /// <summary>
        /// The peptide of <paramref name="protein"/> with the given modifications, keyed as AllModsOneIsNterminus is.
        /// </summary>
        private static PeptideWithSetModifications Peptide(string protein, Dictionary<int, Modification> mods)
        {
            var unmodified = (PeptideWithSetModifications)new Protein(protein, "P1")
                .Digest(new DigestionParams(minPeptideLength: 1), new List<Modification>(), new List<Modification>())
                .First(p => p.Length == protein.Length);
            return (PeptideWithSetModifications)PositionalIsoforms.WithModifications(unmodified, mods);
        }

        private static List<PositionalIsoform> Enumerate(IBioPolymerWithSetMods peptide, int max = PositionalIsoforms.DefaultMaxPlacements,
            params Modification[] fixedMods)
        {
            Assert.That(PositionalIsoforms.TryEnumerate(peptide, new[] { PhosphoS, PhosphoT, OxidationM }, m => fixedMods.Contains(m), max, out var isoforms), Is.True);
            return isoforms;
        }

        private static string[] Placements(List<PositionalIsoform> isoforms) => isoforms.Select(i => i.Placement.ToString()).ToArray();

        [Test]
        public static void APhosphateMovesToEveryResidueAPhosphorylationFits()
        {
            // reported on S3; T4 takes phosphorylation on T, and P, E and K take none
            var isoforms = Enumerate(Peptide("PESTSK", new() { { 4, PhosphoS } }));

            Assert.That(Placements(isoforms), Is.EquivalentTo(new[] { "S3", "S5", "T4" }));
            Assert.That(isoforms.Single(i => i.Placement.ToString() == "T4").Placement.Sites[0].Modification, Is.SameAs(PhosphoT));
            Assert.That(isoforms.Single(i => i.Placement.ToString() == "S5").AllModsOneIsNterminus, Is.EquivalentTo(new Dictionary<int, Modification> { { 6, PhosphoS } }));
        }

        [Test]
        public static void TwoPhosphatesAreAnArrangementNotAnOrder()
        {
            // four sites, two phosphates: 4 choose 2
            var isoforms = Enumerate(Peptide("SESTSK", new() { { 2, PhosphoS }, { 5, PhosphoT } }));

            Assert.That(Placements(isoforms), Is.EquivalentTo(new[] { "S1,S3", "S1,S5", "S1,T4", "S3,S5", "S3,T4", "T4,S5" }));
        }

        [Test]
        public static void ModificationsOfDifferentChemistryMoveIndependently()
        {
            var isoforms = Enumerate(Peptide("MSMSK", new() { { 2, OxidationM }, { 3, PhosphoS } }));

            Assert.That(Placements(isoforms), Is.EquivalentTo(new[] { "M1,S2", "M1,S4", "S2,M3", "M3,S4" }));
        }

        [Test]
        public static void FixedAndTerminalModificationsStayAndHoldTheirResidues()
        {
            // a fixed modification on S2 is not moved, and S2 is not open to the phosphate
            var fixedOnS = Make("FixedThing", "S", "C2H3N1O1");
            var isoforms = Enumerate(Peptide("ASESK", new() { { 1, AcetylNTerm }, { 3, fixedOnS }, { 5, PhosphoS } }), fixedMods: fixedOnS);

            Assert.That(Placements(isoforms), Is.EquivalentTo(new[] { "S4" }));
            Assert.That(isoforms[0].AllModsOneIsNterminus[1], Is.SameAs(AcetylNTerm));
            Assert.That(isoforms[0].AllModsOneIsNterminus[3], Is.SameAs(fixedOnS));
        }

        [Test]
        public static void AnUnmodifiedPeptideHasOneEmptyPlacement()
        {
            var isoforms = Enumerate(Peptide("PESTSK", new()));

            Assert.That(isoforms.Count, Is.EqualTo(1));
            Assert.That(isoforms[0].Placement.Sites, Is.Empty);
        }

        [Test]
        public static void MoreArrangementsThanTheLimitAreRefused()
        {
            var peptide = Peptide("PESTSK", new() { { 4, PhosphoS } });

            Assert.That(PositionalIsoforms.TryEnumerate(peptide, new[] { PhosphoS, PhosphoT }, _ => false, 2, out var isoforms), Is.False);
            Assert.That(isoforms, Is.Null);
            Assert.That(PositionalIsoforms.TryEnumerate(peptide, new[] { PhosphoS, PhosphoT }, _ => false, 3, out isoforms), Is.True);
        }

        [Test]
        public static void EngineFindsThePhosphateWhereTheSpectrumHasIt()
        {
            // the spectrum is of the phosphate on T4; the match reports it on S3
            var commonParameters = new CommonParameters(digestionParams: new DigestionParams(minPeptideLength: 1));
            var fsp = new List<(string fileName, CommonParameters fileSpecificParameters)> { ("", commonParameters) };
            var onT4 = Peptide("PESTSK", new() { { 5, PhosphoT } });
            var onS3 = Peptide("PESTSK", new() { { 4, PhosphoS } });

            MsDataFile dataFile = new TestDataFile(new List<IBioPolymerWithSetMods> { onT4 });
            var scan = new Ms2ScanWithSpecificMass(dataFile.GetAllScansList().Last(), onT4.MonoisotopicMass.ToMz(1), 1, null, new CommonParameters());
            var products = new List<Product>();
            onS3.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
            var matchedIons = MetaMorpheusEngine.MatchFragmentIons(scan, products, new CommonParameters());
            SpectralMatch psm = new PeptideSpectralMatch(onS3, 0, 0, 2, scan, commonParameters, matchedIons);
            psm.ResolveAllAmbiguities();

            new LocalizationEngine(new List<SpectralMatch> { psm }, dataFile, commonParameters, fsp, new List<string>(),
                localizeMassDifferences: false, placeableModifications: new[] { PhosphoS, PhosphoT }).Run();

            LocalizationEvidence evidence = psm.LocalizationEvidence;
            Assert.That(psm.LocalizedScores, Is.Null);
            Assert.That(evidence, Is.Not.Null);
            Assert.That(evidence.Placements.Select(p => p.ToString()), Is.EquivalentTo(new[] { "S3", "T4", "S5" }));
            Assert.That(evidence.MassDifference, Is.EqualTo(0).Within(1e-3));

            int best = Enumerable.Range(0, evidence.CandidateCount).MaxBy(c => evidence.Scores[c]);
            Assert.That(evidence.Placements[best].ToString(), Is.EqualTo("T4"));

            // T4 matches a peak that S3 does not
            int s3 = evidence.Placements.ToList().FindIndex(p => p.ToString() == "S3");
            Assert.That(Enumerable.Range(0, evidence.PeakCount).Any(p => evidence.Matches(best, p) && !evidence.Matches(s3, p)));
        }

        private static List<PositionalIsoform> EnumerateAnyResidue(IBioPolymerWithSetMods peptide, params Modification[] fixedMods)
        {
            Assert.That(PositionalIsoforms.TryEnumerate(peptide, null, m => fixedMods.Contains(m), PositionalIsoforms.DefaultMaxPlacements,
                out var isoforms, LocalizationSiteRule.AnyResidue), Is.True);
            return isoforms;
        }

        [Test]
        public static void AnyResiduePutsTheModificationOnEveryResidueAsItself()
        {
            var isoforms = EnumerateAnyResidue(Peptide("PESTSK", new() { { 4, PhosphoS } }));

            Assert.That(Placements(isoforms), Is.EquivalentTo(new[] { "P1", "E2", "S3", "T4", "S5", "K6" }));
            Assert.That(isoforms.All(i => i.Placement.Sites.Single().Modification == PhosphoS));
        }

        [Test]
        public static void AnyResidueMovesATerminalModificationAndKeepsTheTerminusForTheTerminalResidue()
        {
            var acetylated = Peptide("AKESK", new() { { 1, AcetylNTerm } });

            // under the motif rule a terminal modification stays where it is
            Assert.That(Placements(Enumerate(acetylated)), Is.EquivalentTo(new[] { "" }));

            // under the any-residue rule it moves; A1 is the N-terminus itself, not a second, identical placement
            var isoforms = EnumerateAnyResidue(acetylated);
            Assert.That(Placements(isoforms), Is.EquivalentTo(new[] { "A1", "K2", "E3", "S4", "K5" }));
            Assert.That(isoforms.Single(i => i.Placement.ToString() == "A1").AllModsOneIsNterminus.Keys, Is.EquivalentTo(new[] { 1 }));
            Assert.That(isoforms.Single(i => i.Placement.ToString() == "K2").AllModsOneIsNterminus.Keys, Is.EquivalentTo(new[] { 3 }));
        }

        [Test]
        public static void AnyResidueLetsATerminalAndAResidueModificationShareTheTerminalResidue()
        {
            // the acetyl has keys 1, 3..6 and the phosphate keys 2..6; they may not share a key, but key 1 (the
            // N-terminus) and key 2 (the side chain of A1) are different positions: 5 x 5 - 4 shared keys
            var isoforms = EnumerateAnyResidue(Peptide("AKESK", new() { { 1, AcetylNTerm }, { 5, PhosphoS } }));

            Assert.That(isoforms.Count, Is.EqualTo(21));
            Assert.That(isoforms.Any(i => i.AllModsOneIsNterminus.Count == 2 && i.AllModsOneIsNterminus[1] == AcetylNTerm && i.AllModsOneIsNterminus[2] == PhosphoS));
            Assert.That(isoforms.All(i => i.AllModsOneIsNterminus.Count == 2));
        }

        [Test]
        public static void AnyResidueStillLeavesFixedModificationsWhereTheyAre()
        {
            var fixedOnC = Make("Carbamidomethyl", "C", "C2H3N1O1");
            var isoforms = EnumerateAnyResidue(Peptide("PECSK", new() { { 4, fixedOnC }, { 5, PhosphoS } }), fixedOnC);

            Assert.That(Placements(isoforms), Is.EquivalentTo(new[] { "P1", "E2", "S4", "K5" }));
            Assert.That(isoforms.All(i => i.AllModsOneIsNterminus[4] == fixedOnC));
        }

        [Test]
        public static void AnyResidueFindsAModificationWhereItsMotifSaysItCannotBe()
        {
            // the spectrum has the phosphate's mass on K6, which no phosphorylation motif allows; the match reports S3
            var commonParameters = new CommonParameters(digestionParams: new DigestionParams(minPeptideLength: 1));
            var fsp = new List<(string fileName, CommonParameters fileSpecificParameters)> { ("", commonParameters) };
            var onK6 = Peptide("PESTSK", new() { { 7, PhosphoS } });
            var onS3 = Peptide("PESTSK", new() { { 4, PhosphoS } });

            MsDataFile dataFile = new TestDataFile(new List<IBioPolymerWithSetMods> { onK6 });
            var scan = new Ms2ScanWithSpecificMass(dataFile.GetAllScansList().Last(), onK6.MonoisotopicMass.ToMz(1), 1, null, new CommonParameters());
            var products = new List<Product>();
            onS3.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);

            string BestPlacement(LocalizationSiteRule rule)
            {
                SpectralMatch psm = new PeptideSpectralMatch(onS3, 0, 0, 2, scan, commonParameters, MetaMorpheusEngine.MatchFragmentIons(scan, products, new CommonParameters()));
                psm.ResolveAllAmbiguities();
                new LocalizationEngine(new List<SpectralMatch> { psm }, dataFile, commonParameters, fsp, new List<string>(),
                    localizeMassDifferences: false, placeableModifications: rule == LocalizationSiteRule.Motif ? new[] { PhosphoS, PhosphoT } : null,
                    siteRule: rule).Run();

                LocalizationEvidence evidence = psm.LocalizationEvidence;
                return evidence.Placements[Enumerable.Range(0, evidence.CandidateCount).MaxBy(c => evidence.Scores[c])].ToString();
            }

            Assert.That(BestPlacement(LocalizationSiteRule.Motif), Is.Not.EqualTo("K6"));
            Assert.That(BestPlacement(LocalizationSiteRule.AnyResidue), Is.EqualTo("K6"));
        }

        [Test]
        public static void EngineLeavesMassDifferenceLocalizationAloneByDefault()
        {
            // the original constructor places no identified modifications, so open-search output is unchanged
            var commonParameters = new CommonParameters(digestionParams: new DigestionParams(minPeptideLength: 1));
            var fsp = new List<(string fileName, CommonParameters fileSpecificParameters)> { ("", commonParameters) };
            var onS3 = Peptide("PESTSK", new() { { 4, PhosphoS } });

            MsDataFile dataFile = new TestDataFile(new List<IBioPolymerWithSetMods> { onS3 });
            var scan = new Ms2ScanWithSpecificMass(dataFile.GetAllScansList().Last(), onS3.MonoisotopicMass.ToMz(1), 1, null, new CommonParameters());
            var products = new List<Product>();
            onS3.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
            SpectralMatch psm = new PeptideSpectralMatch(onS3, 0, 0, 2, scan, commonParameters, MetaMorpheusEngine.MatchFragmentIons(scan, products, new CommonParameters()));
            psm.ResolveAllAmbiguities();

            new LocalizationEngine(new List<SpectralMatch> { psm }, dataFile, commonParameters, fsp, new List<string>()).Run();

            Assert.That(psm.LocalizedScores.Count, Is.EqualTo(6));
            // the peptide explains the precursor, so there is no mass difference to place: one unmodified placement,
            // not six that match the same peaks and would each look 1/6 likely
            Assert.That(psm.LocalizationEvidence.CandidateCount, Is.EqualTo(1));
            Assert.That(psm.LocalizationEvidence.Placements.Single().Sites, Is.Empty);
        }
    }
}
