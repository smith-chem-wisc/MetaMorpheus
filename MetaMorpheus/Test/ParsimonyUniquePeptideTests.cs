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
    /// Which peptides a protein group reports as unique. A peptide is unique to a protein group when no other protein
    /// group of the run lists its sequence, however many proteins in the database contain it; "shared" means another
    /// group lists it too. These tests pin that a peptide found in only ONE protein is unique however many times that
    /// protein is listed for it: once per modified form of the peptide; once per full sequence of an ambiguous match
    /// (two or more different full sequences matching the same MS2 scan within the score tolerance); or once per
    /// position when the same sequence occurs twice in the protein. They also pin the group-level meaning: a peptide of
    /// several indistinguishable proteins is unique to the one group they form, and a peptide whose other protein
    /// parsimony did not need is unique to the group that explains it. Parsimony itself still decides which proteins
    /// to keep from database-level uniqueness; only the reported labels are judged among the groups.
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

        private static Protein Protein(string sequence, string accession, bool isDecoy = false, bool isContaminant = false) =>
            new(sequence, accession, null, new List<System.Tuple<string, string>>(), new Dictionary<int, List<Modification>>(),
                isDecoy: isDecoy, isContaminant: isContaminant, databaseFilePath: "test.fasta");

        /// <summary>A peptide of <paramref name="protein"/>; <paramref name="oxidizedResidues"/> are 1-based positions in the peptide.</summary>
        private static PeptideWithSetModifications Peptide(Protein protein, int start, int end, params int[] oxidizedResidues) =>
            new(protein: protein, digestionParams: Common.DigestionParams, oneBasedStartResidueInProtein: start,
                oneBasedEndResidueInProtein: end, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "",
                missedCleavages: 0, allModsOneIsNterminus: oxidizedResidues.ToDictionary(r => r + 1, _ => Oxidation), numFixedMods: 0);

        private static int _scanNumber;

        /// <summary>
        /// One confident PSM. Several peptides make it ambiguous: each matches the same MS2 scan with the same score.
        /// </summary>
        private static SpectralMatch Psm(params PeptideWithSetModifications[] peptides) => Psm("File", 0, peptides);

        /// <summary>One PSM from <paramref name="file"/> at q-value <paramref name="qValue"/>.</summary>
        private static SpectralMatch Psm(string file, double qValue, params PeptideWithSetModifications[] peptides)
        {
            int scanNumber = ++_scanNumber;
            var dataScan = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), scanNumber, 2, true,
                Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=" + scanNumber,
                double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            var scan = new Ms2ScanWithSpecificMass(dataScan, 2, 0, file, new CommonParameters());

            SpectralMatch psm = new PeptideSpectralMatch(peptides[0], 0, 10, 0, scan, Common, new List<MatchedFragmentIon>());
            foreach (var peptide in peptides.Skip(1))
                psm.AddOrReplace(peptide, 10, 0, true, new List<MatchedFragmentIon>());
            psm.ResolveAllAmbiguities();
            psm.SetFdrValues(1, 0, qValue, 1, 0, qValue, 0, 0);
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

        /// <summary>
        /// As a search runs it: parsimony and scoring see every PSM, whatever its q-value, and scoring drops the peptides
        /// with no PSM at 1% (and, with <paramref name="noOneHitWonders"/>, the groups then left with one peptide).
        /// </summary>
        private static List<ProteinGroup> InferProteinsFromEveryPsm(List<SpectralMatch> psms, bool noOneHitWonders)
        {
            var parsimony = (ProteinParsimonyResults)new ProteinParsimonyEngine(psms, false, new CommonParameters(), null, null).Run();
            var scored = (ProteinScoringAndFdrResults)new ProteinScoringAndFdrEngine(parsimony.ProteinGroups, psms,
                noOneHitWonders, false, true, new CommonParameters(), null, new List<string>()).Run();
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
        /// Two proteins with the same peptides merge into one group (calmodulin: three genes, one sequence). No other
        /// group lists ACDMEK, so it is unique to that group in every form, though two database proteins contain it.
        /// </summary>
        [TestCase(false, "ACDMEK", TestName = "IndistinguishableProteinsShareOneGroupAndItsPeptidesAreUnique")]
        [TestCase(true, "ACDM[Test:Oxidation on M]EK|ACDMEK", TestName = "IndistinguishableProteinsShareOneGroupAndItsPeptidesAreUnique_ModifiedPeptidesDifferent")]
        public static void IndistinguishableProteinsShareOneGroupAndItsPeptidesAreUnique(bool modPeptidesAreDifferent, string uniqueCell)
        {
            var p1 = Protein("ACDMEKGGGR", "P1");
            var p2 = Protein("ACDMEKGGGR", "P2");
            var psms = new List<SpectralMatch>
            {
                Psm(Peptide(p1, 1, 6), Peptide(p2, 1, 6)),
                Psm(Peptide(p1, 1, 6, 4), Peptide(p2, 1, 6, 4)),
            };

            var group = InferProteins(psms, modPeptidesAreDifferent).Single();

            Assert.That(group.ProteinGroupName, Is.EqualTo("P1|P2"));
            Assert.That(group.UniquePeptides, Is.EquivalentTo(group.AllPeptides), "no other group lists ACDMEK");
            Assert.That(group.AllPeptides, Has.Count.EqualTo(4), "both forms, from both proteins");
            var columns = Columns(group);
            Assert.That(columns.Unique.Split('|').OrderBy(s => s), Is.EqualTo(uniqueCell.Split('|').OrderBy(s => s)));
            Assert.That(columns.Shared, Is.Empty);
            Assert.That(columns.NumberUnique, Is.EqualTo(modPeptidesAreDifferent ? "2" : "1"));
        }

        /// <summary>
        /// ACDMEK is in P1 and P2, but parsimony keeps only P1: P1 has GGGR of its own and explains ACDMEK too, so P2 is
        /// not needed and forms no group. ACDMEK is then listed by one group only, P1's, and is unique to it. Parsimony is
        /// unchanged: it still kept P1 for GGGR and dropped P2.
        /// </summary>
        [Test]
        public static void APeptideWhoseOtherProteinParsimonyDidNotNeedIsUniqueToTheGroupThatKeptIt()
        {
            var p1 = Protein("ACDMEKGGGR", "P1");
            var p2 = Protein("ACDMEKWWWR", "P2");
            var psms = new List<SpectralMatch>
            {
                Psm(Peptide(p1, 1, 6), Peptide(p2, 1, 6)),
                Psm(Peptide(p1, 7, 10)),
            };

            var group = InferProteins(psms).Single();

            Assert.That(group.ProteinGroupName, Is.EqualTo("P1"), "parsimony keeps P1 and drops P2, as before");
            Assert.That(BaseSequences(group.UniquePeptides), Is.EqualTo(new[] { "ACDMEK", "GGGR" }));
            Assert.That(Columns(group), Is.EqualTo(("ACDMEK|GGGR", "", "2")));
        }

        /// <summary>
        /// The per-file protein table is a subset of the run's group (<see cref="ProteinGroup.ConstructSubsetProteinGroup"/>),
        /// so it reports the run's labels: the indistinguishable pair's peptide is unique in each file's table too.
        /// </summary>
        [Test]
        public static void APerFileProteinGroupKeepsTheGroupLevelLabels()
        {
            var p1 = Protein("ACDMEKGGGR", "P1");
            var p2 = Protein("ACDMEKGGGR", "P2");
            var psms = new List<SpectralMatch> { Psm(Peptide(p1, 1, 6), Peptide(p2, 1, 6)) };

            var group = InferProteins(psms).Single();
            var subset = group.ConstructSubsetProteinGroup(psms[0].FullFilePath);

            Assert.That(subset.UniquePeptides, Is.EquivalentTo(group.AllPeptides));
            Assert.That(Columns(subset), Is.EqualTo(("ACDMEK", "", "1")));
        }

        /// <summary>
        /// The sequence is the one parsimony compares. P1 lists ACDMEK unmodified and P2 lists it oxidized. With modified
        /// peptides treated as different they are two peptides, each listed by one group, so each is unique; otherwise
        /// they are one peptide, listed by both groups, so it is shared. The groups are built by hand: with modified
        /// peptides not different, parsimony gives every protein carrying a base sequence every form of it seen, so the
        /// groups it builds never split the forms like this, and only a hand-built pair shows which key is used.
        /// </summary>
        [TestCase(true, TestName = "WithModifiedPeptidesTreatedAsDifferent_TwoFormsOfOneSequenceInTwoGroupsAreEachUnique")]
        [TestCase(false, TestName = "WithModifiedPeptidesNotDifferent_TwoFormsOfOneSequenceInTwoGroupsAreShared")]
        public static void TwoFormsOfOneSequenceInTwoGroups(bool modPeptidesAreDifferent)
        {
            var p1 = Protein("ACDMEKGGGR", "P1");
            var p2 = Protein("ACDMEKWWWR", "P2");
            var unmodifiedOfP1 = Peptide(p1, 1, 6);
            var oxidizedOfP2 = Peptide(p2, 1, 6, 4);
            var gggr = Peptide(p1, 7, 10);
            var wwwr = Peptide(p2, 7, 10);
            var psms = new List<SpectralMatch> { Psm(unmodifiedOfP1), Psm(oxidizedOfP2), Psm(gggr), Psm(wwwr) };
            var handBuilt = new List<ProteinGroup>
            {
                new(new HashSet<IBioPolymer> { p1 }, new HashSet<IBioPolymerWithSetMods> { unmodifiedOfP1, gggr }, new HashSet<IBioPolymerWithSetMods>()),
                new(new HashSet<IBioPolymer> { p2 }, new HashSet<IBioPolymerWithSetMods> { oxidizedOfP2, wwwr }, new HashSet<IBioPolymerWithSetMods>()),
            };

            var groups = ((ProteinScoringAndFdrResults)new ProteinScoringAndFdrEngine(handBuilt, psms, false, modPeptidesAreDifferent, true,
                new CommonParameters(), null, new List<string>()).Run()).SortedAndScoredProteinGroups.ToDictionary(g => g.ProteinGroupName);

            Assert.That(groups.Keys, Is.EquivalentTo(new[] { "P1", "P2" }));
            var expectedP1 = modPeptidesAreDifferent ? new[] { "ACDMEK", "GGGR" } : new[] { "GGGR" };
            var expectedP2 = modPeptidesAreDifferent ? new[] { "ACDMEK", "WWWR" } : new[] { "WWWR" };
            Assert.That(BaseSequences(groups["P1"].UniquePeptides), Is.EqualTo(expectedP1));
            Assert.That(BaseSequences(groups["P2"].UniquePeptides), Is.EqualTo(expectedP2));
        }

        /// <summary>
        /// Every group counts, whatever its q-value or decoy or contaminant status. P1 and a decoy (or contaminant)
        /// group both list ACDMEK, so it is shared in P1, though a decoy group is never reported.
        /// </summary>
        [TestCase(true, false, TestName = "APeptideAlsoListedByADecoyGroupIsShared")]
        [TestCase(false, true, TestName = "APeptideAlsoListedByAContaminantGroupIsShared")]
        public static void APeptideAlsoListedByADecoyOrContaminantGroupIsShared(bool otherIsDecoy, bool otherIsContaminant)
        {
            var p1 = Protein("ACDMEKGGGR", "P1");
            var other = Protein("ACDMEKWWWR", otherIsDecoy ? "DECOY_P2" : "P2", isDecoy: otherIsDecoy, isContaminant: otherIsContaminant);
            var psms = new List<SpectralMatch>
            {
                Psm(Peptide(p1, 1, 6)),
                Psm(Peptide(other, 1, 6)),
                Psm(Peptide(p1, 7, 10)),
                Psm(Peptide(other, 7, 10)),
            };

            var groups = InferProteins(psms).ToDictionary(g => g.ProteinGroupName);

            Assert.That(groups.Keys, Is.EquivalentTo(new[] { "P1", other.Accession }));
            Assert.That(groups[other.Accession].IsDecoy, Is.EqualTo(otherIsDecoy));
            Assert.That(groups[other.Accession].IsContaminant, Is.EqualTo(otherIsContaminant));
            Assert.That(BaseSequences(groups["P1"].UniquePeptides), Is.EqualTo(new[] { "GGGR" }));
            Assert.That(Columns(groups["P1"]), Is.EqualTo(("GGGR", "ACDMEK", "1")));
        }

        /// <summary>
        /// The groups judged are the ones scoring keeps. P2 is handed to scoring with ACDMEK alone, which P1 also lists.
        /// With the two-peptide rule off, P2 stays and ACDMEK is shared in P1; with it on, P2 is dropped and ACDMEK is
        /// unique to P1, the one group that reports it. FlashLFQ and the isobaric engine quantify against the same kept
        /// groups. Parsimony itself never forms such a P2 (P1 already explains ACDMEK), so the groups are built by hand.
        /// </summary>
        [TestCase(false, TestName = "AGroupKeptWithOnePeptideMakesThatPeptideShared")]
        [TestCase(true, TestName = "AGroupTheTwoPeptideRuleDropsDoesNotMakeAPeptideShared")]
        public static void TheTwoPeptideRuleDecidesWhichGroupsAreJudged(bool noOneHitWonders)
        {
            var p1 = Protein("ACDMEKGGGR", "P1");
            var p2 = Protein("ACDMEKWWWR", "P2");
            var acdmekOfP1 = Peptide(p1, 1, 6);
            var acdmekOfP2 = Peptide(p2, 1, 6);
            var gggr = Peptide(p1, 7, 10);
            var psms = new List<SpectralMatch> { Psm(acdmekOfP1, acdmekOfP2), Psm(gggr) };
            var handBuilt = new List<ProteinGroup>
            {
                new(new HashSet<IBioPolymer> { p1 }, new HashSet<IBioPolymerWithSetMods> { acdmekOfP1, gggr }, new HashSet<IBioPolymerWithSetMods> { gggr }),
                new(new HashSet<IBioPolymer> { p2 }, new HashSet<IBioPolymerWithSetMods> { acdmekOfP2 }, new HashSet<IBioPolymerWithSetMods>()),
            };

            var groups = ((ProteinScoringAndFdrResults)new ProteinScoringAndFdrEngine(handBuilt, psms, noOneHitWonders, false, true,
                new CommonParameters(), null, new List<string>()).Run()).SortedAndScoredProteinGroups.ToDictionary(g => g.ProteinGroupName);

            Assert.That(groups.Keys, Is.EquivalentTo(noOneHitWonders ? new[] { "P1" } : new[] { "P1", "P2" }));
            Assert.That(BaseSequences(groups["P1"].UniquePeptides),
                Is.EqualTo(noOneHitWonders ? new[] { "ACDMEK", "GGGR" } : new[] { "GGGR" }));
        }

        /// <summary>
        /// One group is judged once, however often the list holds it. The picked FDR pairs a group with the decoy of each
        /// of its accessions, and a target group that loses to the decoys of two accessions is added back once for each
        /// (it is rescued by the classic FDR). P1|P2 then appears twice, and must not make its own peptides shared.
        /// </summary>
        [Test]
        public static void AGroupListedTwiceAfterThePickedFdrKeepsItsUniquePeptides()
        {
            var p1 = Protein("ACDMEKGGGR", "P1");
            var p2 = Protein("ACDMEKGGGR", "P2");
            var decoy1 = Protein("EEEKWWWR", "DECOY_P1", isDecoy: true);
            var decoy2 = Protein("DDDKYYYR", "DECOY_P2", isDecoy: true);
            var psms = new List<SpectralMatch>
            {
                Psm("File", 0.005, Peptide(p1, 1, 6), Peptide(p2, 1, 6)),
                Psm("File", 0.005, Peptide(p1, 7, 10), Peptide(p2, 7, 10)),
                Psm(Peptide(decoy1, 1, 4)), Psm(Peptide(decoy1, 5, 8)),
                Psm(Peptide(decoy2, 1, 4)), Psm(Peptide(decoy2, 5, 8)),
            };

            var groups = InferProteins(psms);

            var targets = groups.Where(g => g.ProteinGroupName == "P1|P2").ToList();
            Assume.That(targets, Has.Count.EqualTo(2), "the picked FDR adds the rescued target back once per accession");
            Assume.That(targets[0], Is.SameAs(targets[1]));
            Assert.That(BaseSequences(targets[0].UniquePeptides), Is.EqualTo(new[] { "ACDMEK", "GGGR" }));
            var columns = Columns(targets[0]);
            Assert.That(columns.Unique.Split('|').OrderBy(s => s), Is.EqualTo(new[] { "ACDMEK", "GGGR" }));
            Assert.That(columns.Shared, Is.Empty);
            Assert.That(columns.NumberUnique, Is.EqualTo("2"));
        }

        /// <summary>
        /// A per-file scoring call works on subsets of the search's groups and keeps the search's labels, even where its
        /// own two-peptide rule drops a subset. ACDMEK is shared by P1 and P2 in the search; in file A, P2's subset has
        /// ACDMEK only and is dropped, but P1's subset still reports ACDMEK as shared.
        /// </summary>
        [Test]
        public static void APerFileScoringCallKeepsTheSearchsLabels()
        {
            var p1 = Protein("ACDMEKGGGR", "P1");
            var p2 = Protein("ACDMEKWWWR", "P2");
            var psms = new List<SpectralMatch>
            {
                Psm("A", 0, Peptide(p1, 1, 6), Peptide(p2, 1, 6)),
                Psm("A", 0, Peptide(p1, 7, 10)),
                Psm("B", 0, Peptide(p2, 7, 10)),
            };
            var groups = InferProteinsFromEveryPsm(psms, noOneHitWonders: false);
            Assume.That(BaseSequences(groups.Single(g => g.ProteinGroupName == "P1").UniquePeptides), Is.EqualTo(new[] { "GGGR" }));

            var subsets = groups.Select(g => g.ConstructSubsetProteinGroup("A")).ToList();
            var psmsOfA = psms.Where(p => p.FullFilePath == "A").ToList();
            var perFile = ((ProteinScoringAndFdrResults)new ProteinScoringAndFdrEngine(subsets, psmsOfA, true, false, false,
                new CommonParameters(), null, new List<string>(), labelPeptidesUniqueToOneGroup: false).Run()).SortedAndScoredProteinGroups;

            var p1InA = perFile.Single();
            Assert.That(p1InA.ProteinGroupName, Is.EqualTo("P1"), "P2's subset in A has one peptide, so the rule drops it");
            Assert.That(Columns(p1InA), Is.EqualTo(("GGGR", "ACDMEK", "1")));
        }

        /// <summary>
        /// mzIdentML's MS:1001097 for a group is its number of distinct unique peptides, as the protein table counts them.
        /// A group of two indistinguishable proteins holds each peptide twice (once per protein), and must still report 2.
        /// </summary>
        [Test]
        public static void MzIdentMlCountsEachUniquePeptideOnceInAGroupOfIndistinguishableProteins()
        {
            var p1 = Protein("ACDMEKGGGR", "P1");
            var p2 = Protein("ACDMEKGGGR", "P2");
            var psms = new List<SpectralMatch>
            {
                Psm(Peptide(p1, 1, 6), Peptide(p2, 1, 6)),
                Psm(Peptide(p1, 7, 10), Peptide(p2, 7, 10)),
            };
            var group = InferProteins(psms).Single();
            Assume.That(group.ProteinGroupName, Is.EqualTo("P1|P2"));
            Assume.That(group.UniquePeptides, Has.Count.EqualTo(4), "two peptides, once per protein");
            Assume.That(Columns(group).NumberUnique, Is.EqualTo("2"));

            string path = System.IO.Path.Combine(TestContext.CurrentContext.TestDirectory, "unique_count.mzID");
            try
            {
                TaskLayer.MzIdentMLWriter.WriteMzIdentMl(psms, new List<ProteinGroup> { group }, new List<Modification>(),
                    new List<Modification>(), null, new List<DigestionAgent>(), new MzLibUtil.PpmTolerance(20),
                    new MzLibUtil.PpmTolerance(20), 0, path, true);

                var counts = System.Xml.Linq.XDocument.Load(path).Descendants()
                    .Where(e => e.Name.LocalName == "cvParam" && (string)e.Attribute("accession") == "MS:1001097")
                    .Select(e => (string)e.Attribute("value"))
                    .ToList();
                Assert.That(counts, Is.EqualTo(new[] { "2", "2" }), "one hypothesis per protein, each with the group's count");
            }
            finally
            {
                if (System.IO.File.Exists(path)) System.IO.File.Delete(path);
            }
        }
    }
}
