using EngineLayer;
using EngineLayer.SpectrumMatch;
using MassSpectrometry;
using NUnit.Framework;
using Omics;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    /// <summary>
    /// Which peptides protein parsimony calls unique to a protein. A peptide is unique when only one protein in the
    /// database contains it (among the proteins with an observed peptide); "shared" means another protein contains it
    /// too. These tests pin that a peptide found in only ONE protein is unique however many times that protein is
    /// listed for it: once per modified form of the peptide; once per full sequence of an ambiguous match (two or more
    /// different full sequences matching the same MS2 scan within the score tolerance); or once per position when the
    /// same sequence occurs twice in the protein.
    /// </summary>
    [TestFixture]
    public static class ParsimonyUniquePeptideTests
    {
        private static readonly CommonParameters Common = new(digestionParams: new DigestionParams(protease: "trypsin", minPeptideLength: 1));

        private static Modification Oxidation
        {
            get
            {
                ModificationMotif.TryGetMotif("M", out var motif);
                return new Modification(_originalId: "Oxidation on M", _modificationType: "Test", _target: motif,
                    _locationRestriction: "Anywhere.", _monoisotopicMass: 15.994915);
            }
        }

        private static Protein Protein(string sequence, string accession) =>
            new(sequence, accession, null, new List<System.Tuple<string, string>>(), new Dictionary<int, List<Modification>>());

        /// <summary>A peptide of <paramref name="protein"/>; <paramref name="oxidizedResidues"/> are 1-based positions in the peptide.</summary>
        private static PeptideWithSetModifications Peptide(Protein protein, int start, int end, params int[] oxidizedResidues) =>
            new(protein: protein, digestionParams: Common.DigestionParams, oneBasedStartResidueInProtein: start,
                oneBasedEndResidueInProtein: end, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "",
                missedCleavages: 0, allModsOneIsNterminus: oxidizedResidues.ToDictionary(r => r + 1, _ => Oxidation), numFixedMods: 0);

        private static int _scanNumber;

        /// <summary>
        /// One confident PSM. Several peptides make it ambiguous: each matches the same MS2 scan with the same score.
        /// </summary>
        private static SpectralMatch Psm(params PeptideWithSetModifications[] peptides)
        {
            int scanNumber = ++_scanNumber;
            var dataScan = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), scanNumber, 2, true,
                Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=" + scanNumber,
                double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            var scan = new Ms2ScanWithSpecificMass(dataScan, 2, 0, "File", new CommonParameters());

            SpectralMatch psm = new PeptideSpectralMatch(peptides[0], 0, 10, 0, scan, Common, new List<MatchedFragmentIon>());
            foreach (var peptide in peptides.Skip(1))
                psm.AddOrReplace(peptide, 10, 0, true, new List<MatchedFragmentIon>());
            psm.ResolveAllAmbiguities();
            psm.SetFdrValues(1, 0, 0, 1, 0, 0, 0, 0);
            return psm;
        }

        /// <summary>Parsimony, then scoring with indistinguishable groups merged, as a search does.</summary>
        private static List<ProteinGroup> InferProteins(List<SpectralMatch> psms, bool modPeptidesAreDifferent = false)
        {
            var filtered = FilteredPsms.Filter(psms, Common);
            var parsimony = (ProteinParsimonyResults)new ProteinParsimonyEngine(filtered.FilteredPsmsList,
                modPeptidesAreDifferent, new CommonParameters(), null, null).Run();
            var scored = (ProteinScoringAndFdrResults)new ProteinScoringAndFdrEngine(parsimony.ProteinGroups,
                filtered.FilteredPsmsList, false, modPeptidesAreDifferent, true, new CommonParameters(), null, new List<string>()).Run();
            return scored.SortedAndScoredProteinGroups;
        }

        /// <summary>The group's "Unique Peptides" and "Shared Peptides" cells, as the protein table writes them.</summary>
        private static (string Unique, string Shared, string NumberUnique) Columns(ProteinGroup group)
        {
            group.GetIdentifiedPeptidesOutput(null);
            var header = group.GetTabSeparatedHeader().Split('\t').ToList();
            var row = group.ToString().Split('\t');
            return (row[header.IndexOf("Unique Peptides")], row[header.IndexOf("Shared Peptides")],
                row[header.IndexOf("Number of Unique Peptides")]);
        }

        private static string[] BaseSequences(IEnumerable<IBioPolymerWithSetMods> peptides) =>
            peptides.Select(p => p.BaseSequence).Distinct().OrderBy(s => s).ToArray();

        [Test]
        public static void APeptideSeenInTwoModifiedFormsInOneProteinIsUnique()
        {
            var protein = Protein("ACDMEKGGGR", "P1");
            var unmodified = Peptide(protein, 1, 6);
            var oxidized = Peptide(protein, 1, 6, 4);

            var group = InferProteins(new List<SpectralMatch> { Psm(unmodified), Psm(oxidized) }).Single();

            Assert.That(group.UniquePeptides, Is.EquivalentTo(new[] { unmodified, oxidized }),
                "both forms are unique: no other protein contains ACDMEK");
            Assert.That(Columns(group), Is.EqualTo(("ACDMEK", "", "1")));
        }

        /// <summary>
        /// An ambiguous match: two full sequences (oxidation on M3, or on M5) match one MS2 scan with the same score.
        /// Both come from the one protein, so the base sequence is unique to it.
        /// </summary>
        [Test]
        public static void AnAmbiguousMatchBetweenTwoFullSequencesOfOneProteinIsUnique()
        {
            var protein = Protein("ACMDMEKGGGR", "P1");
            var oxidizedAt3 = Peptide(protein, 1, 7, 3);
            var oxidizedAt5 = Peptide(protein, 1, 7, 5);
            var ambiguous = Psm(oxidizedAt3, oxidizedAt5);
            Assume.That(ambiguous.FullSequence, Is.Null, "two full sequences match the scan");
            Assume.That(ambiguous.BaseSequence, Is.EqualTo("ACMDMEK"), "one base sequence, so parsimony uses the match");

            var group = InferProteins(new List<SpectralMatch> { ambiguous }).Single();

            Assert.That(group.UniquePeptides, Is.EquivalentTo(new[] { oxidizedAt3, oxidizedAt5 }));
            Assert.That(Columns(group), Is.EqualTo(("ACMDMEK", "", "1")));
        }

        [Test]
        public static void APeptideOccurringTwiceInOneProteinIsUnique()
        {
            var protein = Protein("ACDEKACDEKR", "P1");
            var first = Peptide(protein, 1, 5);
            var second = Peptide(protein, 6, 10);

            var group = InferProteins(new List<SpectralMatch> { Psm(first, second) }).Single();

            Assert.That(BaseSequences(group.UniquePeptides), Is.EqualTo(new[] { "ACDEK" }));
            Assert.That(Columns(group).Shared, Is.Empty);
        }

        [Test]
        public static void WithModifiedPeptidesTreatedAsDifferent_ARepeatedPeptideInOneProteinIsUnique()
        {
            var protein = Protein("ACDEKACDEKR", "P1");
            var first = Peptide(protein, 1, 5);
            var second = Peptide(protein, 6, 10);

            var group = InferProteins(new List<SpectralMatch> { Psm(first, second) }, modPeptidesAreDifferent: true).Single();

            Assert.That(BaseSequences(group.UniquePeptides), Is.EqualTo(new[] { "ACDEK" }));
        }

        /// <summary>Control: a peptide in two proteins stays shared, in every modified form.</summary>
        [Test]
        public static void APeptideInTwoProteinsStaysShared()
        {
            var p1 = Protein("ACDMEKGGGR", "P1");
            var p2 = Protein("ACDMEKWWWR", "P2");
            var psms = new List<SpectralMatch>
            {
                Psm(Peptide(p1, 1, 6), Peptide(p2, 1, 6)),
                Psm(Peptide(p1, 1, 6, 4), Peptide(p2, 1, 6, 4)),
                Psm(Peptide(p1, 7, 10)),
                Psm(Peptide(p2, 7, 10)),
            };

            var groups = InferProteins(psms).ToDictionary(g => g.ProteinGroupName);

            Assert.That(groups.Keys, Is.EquivalentTo(new[] { "P1", "P2" }));
            Assert.That(BaseSequences(groups["P1"].UniquePeptides), Is.EqualTo(new[] { "GGGR" }));
            Assert.That(BaseSequences(groups["P2"].UniquePeptides), Is.EqualTo(new[] { "WWWR" }));
            Assert.That(Columns(groups["P1"]), Is.EqualTo(("GGGR", "ACDMEK", "1")));
            Assert.That(Columns(groups["P2"]), Is.EqualTo(("WWWR", "ACDMEK", "1")));
        }

        /// <summary>Control: a peptide seen once in one protein was already unique.</summary>
        [Test]
        public static void APeptideSeenOnceInOneProteinIsUnique()
        {
            var protein = Protein("ACDMEKGGGR", "P1");
            var peptide = Peptide(protein, 1, 6);

            var group = InferProteins(new List<SpectralMatch> { Psm(peptide) }).Single();

            Assert.That(group.UniquePeptides, Is.EquivalentTo(new[] { peptide }));
            Assert.That(Columns(group), Is.EqualTo(("ACDMEK", "", "1")));
        }

        /// <summary>
        /// Control, and the meaning this fix keeps: two proteins with the same peptides merge into one group, and a
        /// peptide both contain is not unique to either protein, so the merged group lists it as shared.
        /// </summary>
        [Test]
        public static void IndistinguishableProteinsStillShareTheirPeptides()
        {
            var p1 = Protein("ACDMEKGGGR", "P1");
            var p2 = Protein("ACDMEKGGGR", "P2");
            var psms = new List<SpectralMatch>
            {
                Psm(Peptide(p1, 1, 6), Peptide(p2, 1, 6)),
                Psm(Peptide(p1, 1, 6, 4), Peptide(p2, 1, 6, 4)),
            };

            var group = InferProteins(psms).Single();

            Assert.That(group.ProteinGroupName, Is.EqualTo("P1|P2"));
            Assert.That(group.UniquePeptides, Is.Empty);
            Assert.That(Columns(group), Is.EqualTo(("", "ACDMEK", "0")));
        }
    }
}
