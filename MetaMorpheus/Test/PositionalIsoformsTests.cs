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
            Assert.That(psm.LocalizationEvidence.CandidateCount, Is.EqualTo(6));
            Assert.That(psm.LocalizationEvidence.Placements.All(p => p.Sites.Single().Modification == null));
        }
    }
}
