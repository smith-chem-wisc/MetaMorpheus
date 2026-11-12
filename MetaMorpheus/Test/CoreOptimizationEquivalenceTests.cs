using EngineLayer;
using MassSpectrometry;
using NUnit.Framework;
using Omics;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;
using Omics.SpectralMatch;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    /// <summary>
    /// The allocation and speed work on fragment coverage, protein sequence coverage and the indistinguishable-group
    /// merge is meant to change no output. Each test here keeps a copy of the implementation it replaced and compares
    /// the two on many generated inputs. The one intended difference, in the merge, has its own tests.
    /// </summary>
    [TestFixture]
    public static class CoreOptimizationEquivalenceTests
    {
        private static readonly CommonParameters Common = new(digestionParams: new DigestionParams(protease: "trypsin", minPeptideLength: 1));
        private const string Residues = "ACDEFGHIKLMNPQRSTVWY";

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
            new(sequence, accession, null, new List<Tuple<string, string>>(), new Dictionary<int, List<Modification>>(),
                databaseFilePath: "test.fasta");

        /// <summary>A peptide of <paramref name="protein"/>; <paramref name="oxidizedResidues"/> are 1-based positions in the peptide.</summary>
        private static PeptideWithSetModifications Peptide(Protein protein, int start, int end, params int[] oxidizedResidues) =>
            new(protein: protein, digestionParams: Common.DigestionParams, oneBasedStartResidueInProtein: start,
                oneBasedEndResidueInProtein: end, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "",
                missedCleavages: 0, allModsOneIsNterminus: oxidizedResidues.ToDictionary(r => r + 1, _ => Oxidation), numFixedMods: 0);

        private static int _scanNumber;

        /// <summary>One confident PSM of score 10. Several peptides make it ambiguous between them.</summary>
        private static SpectralMatch Psm(List<MatchedFragmentIon> matchedIons, params PeptideWithSetModifications[] peptides)
        {
            int scanNumber = ++_scanNumber;
            var dataScan = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), scanNumber, 2, true,
                Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=" + scanNumber,
                double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            var scan = new Ms2ScanWithSpecificMass(dataScan, 2, 0, "File", new CommonParameters());

            SpectralMatch psm = new PeptideSpectralMatch(peptides[0], 0, 10, 0, scan, Common, matchedIons);
            foreach (var peptide in peptides.Skip(1))
            {
                psm.AddOrReplace(peptide, 10, 0, true, matchedIons);
            }
            psm.ResolveAllAmbiguities();
            psm.SetFdrValues(1, 0, 0, 1, 0, 0, 0, 0);
            return psm;
        }

        private static string RandomSequence(Random random, int length) =>
            new(Enumerable.Range(0, length).Select(_ => Residues[random.Next(Residues.Length)]).ToArray());

        /// <summary>Matched ions at random residue positions 1..<paramref name="length"/>, from either terminus.</summary>
        private static List<MatchedFragmentIon> RandomMatchedIons(Random random, int length)
        {
            var ions = new List<MatchedFragmentIon>();
            double density = random.NextDouble();
            for (int position = 1; position <= length; position++)
            {
                if (random.NextDouble() < density)
                {
                    ions.Add(new MatchedFragmentIon(new Product(ProductType.b, FragmentationTerminus.N, 1, position, position, 0), 1, 1, 1));
                }
                if (random.NextDouble() < density)
                {
                    ions.Add(new MatchedFragmentIon(new Product(ProductType.y, FragmentationTerminus.C, 1, position, position, 0), 1, 1, 1));
                }
            }
            return ions;
        }

        #region Fragment coverage of a PSM

        /// <summary>SpectralMatch.GetAminoAcidCoverage as it was before the pooled-array rewrite.</summary>
        private static List<int> FragmentCoverageBeforeRewrite(SpectralMatch psm)
        {
            var nTermFragmentAAPositions = psm.MatchedFragmentIons.Where(p =>
                    p.NeutralTheoreticalProduct.Terminus is FragmentationTerminus.N or FragmentationTerminus.FivePrime)
                .Select(j => j.NeutralTheoreticalProduct.AminoAcidPosition).Distinct().ToList();

            var cTermFragmentAAPositions = psm.MatchedFragmentIons.Where(p =>
                    p.NeutralTheoreticalProduct.Terminus is FragmentationTerminus.C or FragmentationTerminus.ThreePrime)
                .Select(j => j.NeutralTheoreticalProduct.AminoAcidPosition).Distinct().ToList();

            HashSet<int> fragmentCoveredAminoAcids = new();

            if (nTermFragmentAAPositions.Any())
            {
                nTermFragmentAAPositions.Sort();

                if (nTermFragmentAAPositions.Contains(psm.BaseSequence.Length - 1))
                {
                    fragmentCoveredAminoAcids.Add(psm.BaseSequence.Length);
                }

                if (nTermFragmentAAPositions.Contains(1))
                {
                    fragmentCoveredAminoAcids.Add(1);
                }

                for (int i = 0; i < nTermFragmentAAPositions.Count - 1; i++)
                {
                    if (nTermFragmentAAPositions[i + 1] - nTermFragmentAAPositions[i] == 1)
                    {
                        fragmentCoveredAminoAcids.Add(nTermFragmentAAPositions[i + 1]);
                    }

                    if (cTermFragmentAAPositions.Contains(nTermFragmentAAPositions[i + 1]))
                    {
                        fragmentCoveredAminoAcids.Add(nTermFragmentAAPositions[i + 1]);
                    }

                    if (cTermFragmentAAPositions.Contains(nTermFragmentAAPositions[i + 1] + 2))
                    {
                        fragmentCoveredAminoAcids.Add(nTermFragmentAAPositions[i + 1] + 1);
                    }
                }
            }

            if (cTermFragmentAAPositions.Any())
            {
                cTermFragmentAAPositions.Sort();

                if (cTermFragmentAAPositions.Contains(2))
                {
                    fragmentCoveredAminoAcids.Add(1);
                }

                if (cTermFragmentAAPositions.Contains(psm.BaseSequence.Length))
                {
                    fragmentCoveredAminoAcids.Add(psm.BaseSequence.Length);
                }

                for (int i = 0; i < cTermFragmentAAPositions.Count - 1; i++)
                {
                    if (cTermFragmentAAPositions[i + 1] - cTermFragmentAAPositions[i] == 1)
                    {
                        fragmentCoveredAminoAcids.Add(cTermFragmentAAPositions[i]);
                    }
                }
            }

            var fragmentCoveredAminoAcidsList = fragmentCoveredAminoAcids.ToList();
            fragmentCoveredAminoAcidsList.Sort();
            return fragmentCoveredAminoAcidsList;
        }

        [Test]
        public static void FragmentCoverageOfAPsmIsUnchanged()
        {
            var random = new Random(1);
            for (int trial = 0; trial < 3000; trial++)
            {
                int length = random.Next(1, 40);
                var protein = Protein(RandomSequence(random, length), "P1");
                var matchedIons = RandomMatchedIons(random, length);
                if (matchedIons.Count == 0)
                {
                    continue; // the method returns early and leaves the coverage unset
                }
                var psm = Psm(matchedIons, Peptide(protein, 1, length));

                psm.GetAminoAcidCoverage();

                Assert.That(psm.FragmentCoveragePositionInPeptide, Is.EqualTo(FragmentCoverageBeforeRewrite(psm)),
                    $"trial {trial}: {protein.BaseSequence} with ions at " +
                    string.Join(",", matchedIons.Select(i => i.NeutralTheoreticalProduct.Terminus + i.NeutralTheoreticalProduct.AminoAcidPosition.ToString())));
            }
        }

        #endregion

        #region Sequence coverage of a protein group

        /// <summary>ProteinGroup.CalculateSequenceCoverage as it was before the pooled-array rewrite.</summary>
        private static (List<double> Fraction, List<string> Display, List<string> DisplayWithMods, List<string> FragmentDisplay)
            SequenceCoverageBeforeRewrite(ProteinGroup group)
        {
            var fractions = new List<double>();
            var displays = new List<string>();
            var displaysWithMods = new List<string>();
            var fragmentDisplays = new List<string>();

            var proteinsWithUnambigSeqPsms = new Dictionary<IBioPolymer, List<IBioPolymerWithSetMods>>();
            var proteinsWithPsmsWithLocalizedMods = new Dictionary<IBioPolymer, List<IBioPolymerWithSetMods>>();

            foreach (var protein in group.Proteins)
            {
                proteinsWithUnambigSeqPsms.Add(protein, new List<IBioPolymerWithSetMods>());
                proteinsWithPsmsWithLocalizedMods.Add(protein, new List<IBioPolymerWithSetMods>());
            }

            foreach (var psm in group.AllPsmsBelowOnePercentFDR.Cast<SpectralMatch>())
            {
                if (psm.BaseSequence != null)
                {
                    psm.GetAminoAcidCoverage();

                    foreach (var peptide in psm.BestMatchingBioPolymersWithSetMods
                        .Select(p => p.SpecificBioPolymer).DistinctBy(pep => pep.FullSequence))
                    {
                        if (group.Proteins.Contains(peptide.Parent))
                        {
                            proteinsWithUnambigSeqPsms[peptide.Parent].Add(peptide);

                            if (peptide.FullSequence != null)
                            {
                                proteinsWithPsmsWithLocalizedMods[peptide.Parent].Add(peptide);
                            }
                        }
                    }
                }
            }

            foreach (IBioPolymer protein in group.ListOfProteinsOrderedByAccession)
            {
                HashSet<int> coveredResiduesInProteinOneBased = new();

                foreach (SpectralMatch psm in group.AllPsmsBelowOnePercentFDR.Cast<SpectralMatch>()
                    .Where(psm => psm.BaseSequence != null))
                {
                    psm.GetAminoAcidCoverage();
                    if (psm.FragmentCoveragePositionInPeptide == null) continue;

                    IEnumerable<IBioPolymerWithSetMods> pwsms = psm.BestMatchingBioPolymersWithSetMods
                        .Select(p => p.SpecificBioPolymer)
                        .Where(p => p.Parent.Accession == protein.Accession);

                    foreach (var pwsm in pwsms)
                    {
                        HashSet<int> coveredResiduesInPeptide = new();
                        foreach (var position in psm.FragmentCoveragePositionInPeptide)
                        {
                            coveredResiduesInPeptide.Add(position + pwsm.OneBasedStartResidue - 1);
                        }

                        coveredResiduesInProteinOneBased.UnionWith(coveredResiduesInPeptide);
                    }
                }

                char[] fragmentCoverageArray = protein.BaseSequence.ToLower().ToCharArray();
                foreach (var residue in coveredResiduesInProteinOneBased)
                {
                    fragmentCoverageArray[residue - 1] = char.ToUpper(fragmentCoverageArray[residue - 1]);
                }

                fragmentDisplays.Add(new string(fragmentCoverageArray));
            }

            foreach (var protein in group.ListOfProteinsOrderedByAccession)
            {
                HashSet<int> coveredOneBasedResidues = new HashSet<int>();

                foreach (var peptide in proteinsWithUnambigSeqPsms[protein])
                {
                    for (int i = peptide.OneBasedStartResidue; i <= peptide.OneBasedEndResidue; i++)
                    {
                        coveredOneBasedResidues.Add(i);
                    }
                }

                fractions.Add((double)coveredOneBasedResidues.Count / protein.Length);

                string sequenceCoverageDisplay = protein.BaseSequence.ToLower();
                var coverageArray = sequenceCoverageDisplay.ToCharArray();
                foreach (var obsResidueLocation in coveredOneBasedResidues)
                {
                    coverageArray[obsResidueLocation - 1] = char.ToUpper(coverageArray[obsResidueLocation - 1]);
                }

                sequenceCoverageDisplay = new string(coverageArray);
                displays.Add(sequenceCoverageDisplay);

                var modsOnThisProtein = new HashSet<KeyValuePair<int, Modification>>();
                foreach (var pep in proteinsWithPsmsWithLocalizedMods[protein])
                {
                    foreach (var mod in pep.AllModsOneIsNterminus)
                    {
                        if (!mod.Value.ModificationType.Contains("PeptideTermMod")
                            && !mod.Value.ModificationType.Contains("Common Variable")
                            && !mod.Value.ModificationType.Contains("Common Fixed"))
                        {
                            modsOnThisProtein.Add(
                                new KeyValuePair<int, Modification>(pep.OneBasedStartResidue + mod.Key - 2, mod.Value));
                        }
                    }
                }

                var tempMods = modsOnThisProtein.OrderBy(p => p.Key).ToList();
                foreach (var mod in tempMods)
                {
                    if (mod.Value.LocationRestriction.Equals("N-terminal."))
                    {
                        sequenceCoverageDisplay = sequenceCoverageDisplay.Insert(0, $"[{mod.Value.IdWithMotif}]-");
                    }
                    else if (mod.Value.LocationRestriction.Equals("Anywhere."))
                    {
                        int modStringIndex = sequenceCoverageDisplay.Length - (protein.Length - mod.Key);
                        sequenceCoverageDisplay = sequenceCoverageDisplay.Insert(modStringIndex, $"[{mod.Value.IdWithMotif}]");
                    }
                    else if (mod.Value.LocationRestriction.Equals("C-terminal."))
                    {
                        sequenceCoverageDisplay = sequenceCoverageDisplay.Insert(sequenceCoverageDisplay.Length, $"-[{mod.Value.IdWithMotif}]");
                    }
                }

                displaysWithMods.Add(sequenceCoverageDisplay);
            }

            return (fractions, displays, displaysWithMods, fragmentDisplays);
        }

        /// <summary>
        /// One or two proteins of the same sequence, with random overlapping peptides (some oxidized, some seen twice)
        /// and random matched ions.
        /// </summary>
        [Test]
        public static void SequenceCoverageOfAProteinGroupIsUnchanged()
        {
            var random = new Random(2);
            for (int trial = 0; trial < 500; trial++)
            {
                int proteinLength = random.Next(5, 80);
                string sequence = RandomSequence(random, proteinLength);
                var proteins = Enumerable.Range(1, random.Next(1, 3)).Select(i => Protein(sequence, "P" + i)).ToList();

                var peptides = new HashSet<IBioPolymerWithSetMods>();
                var psms = new HashSet<ISpectralMatch>();
                int peptideCount = random.Next(0, 8);
                for (int p = 0; p < peptideCount; p++)
                {
                    int start = random.Next(1, proteinLength + 1);
                    int end = random.Next(start, Math.Min(proteinLength, start + 25) + 1);
                    int[] oxidized = Enumerable.Range(1, end - start + 1)
                        .Where(i => sequence[start + i - 2] == 'M' && random.Next(2) == 0).ToArray();
                    var forms = proteins.Select(protein => Peptide(protein, start, end, oxidized)).ToArray();
                    foreach (var form in forms)
                    {
                        peptides.Add(form);
                    }

                    int psmCount = random.Next(1, 3);
                    for (int s = 0; s < psmCount; s++)
                    {
                        psms.Add(Psm(RandomMatchedIons(random, end - start + 1), forms));
                    }
                }

                var group = new ProteinGroup(new HashSet<IBioPolymer>(proteins), peptides, new HashSet<IBioPolymerWithSetMods>())
                {
                    AllPsmsBelowOnePercentFDR = psms
                };

                var expected = SequenceCoverageBeforeRewrite(group);
                group.CalculateSequenceCoverage();

                Assert.Multiple(() =>
                {
                    Assert.That(group.SequenceCoverageFraction, Is.EqualTo(expected.Fraction), $"trial {trial}");
                    Assert.That(group.SequenceCoverageDisplayList, Is.EqualTo(expected.Display), $"trial {trial}");
                    Assert.That(group.SequenceCoverageDisplayListWithMods, Is.EqualTo(expected.DisplayWithMods), $"trial {trial}");
                    Assert.That(group.FragmentSequenceCoverageDisplayList, Is.EqualTo(expected.FragmentDisplay), $"trial {trial}");
                });
            }
        }

        #endregion

        #region Merging indistinguishable protein groups

        /// <summary>
        /// Proteins built from tryptic blocks, one hand-built group per protein, one PSM per block sequence shared by
        /// every protein that has the block. Every PSM scores 10, so two groups tie when they have the same number of
        /// blocks, and are indistinguishable when they have the same blocks.
        /// </summary>
        private static (List<ProteinGroup> Groups, List<SpectralMatch> Psms) GroupsFromBlocks(IReadOnlyList<string[]> blocksOfEachProtein)
        {
            var proteins = blocksOfEachProtein.Select((blocks, i) => Protein(string.Concat(blocks), $"P{i + 1:00}")).ToList();
            var peptidesOfBlock = new Dictionary<string, List<PeptideWithSetModifications>>();
            var groups = new List<ProteinGroup>();
            for (int i = 0; i < proteins.Count; i++)
            {
                var peptides = new HashSet<IBioPolymerWithSetMods>();
                int start = 1;
                foreach (string block in blocksOfEachProtein[i])
                {
                    var peptide = Peptide(proteins[i], start, start + block.Length - 1);
                    peptides.Add(peptide);
                    if (!peptidesOfBlock.TryGetValue(block, out var list))
                    {
                        peptidesOfBlock.Add(block, list = new List<PeptideWithSetModifications>());
                    }
                    list.Add(peptide);
                    start += block.Length;
                }
                groups.Add(new ProteinGroup(new HashSet<IBioPolymer> { proteins[i] }, peptides, new HashSet<IBioPolymerWithSetMods>()));
            }

            var psms = peptidesOfBlock.Values.Select(p => Psm(new List<MatchedFragmentIon>(), p.ToArray())).ToList();
            return (groups, psms);
        }

        /// <summary>The accessions of each group the engine returns, one sorted "|"-joined string per group.</summary>
        private static List<string> MergedGroups(IReadOnlyList<string[]> blocksOfEachProtein)
        {
            var (groups, psms) = GroupsFromBlocks(blocksOfEachProtein);
            var scored = (ProteinScoringAndFdrResults)new ProteinScoringAndFdrEngine(groups, psms, false, false, true,
                new CommonParameters(), null, new List<string>()).Run();
            return scored.SortedAndScoredProteinGroups
                .Select(g => string.Join("|", g.Proteins.Select(p => p.Accession).OrderBy(a => a, StringComparer.Ordinal)))
                .OrderBy(n => n, StringComparer.Ordinal).ToList();
        }

        private sealed class GroupModel
        {
            public List<string> Accessions;
            public string Blocks;
            public int Score;
        }

        private static List<GroupModel> Model(IReadOnlyList<string[]> blocksOfEachProtein) =>
            blocksOfEachProtein.Select((blocks, i) => new GroupModel
            {
                Accessions = new List<string> { $"P{i + 1:00}" },
                Blocks = string.Join(",", blocks.OrderBy(b => b, StringComparer.Ordinal)),
                Score = 10 * blocks.Distinct().Count()
            }).OrderByDescending(g => g.Score).ToList();

        private static List<string> Names(IEnumerable<GroupModel> groups) =>
            groups.Where(g => g.Score != 0)
                .Select(g => string.Join("|", g.Accessions.OrderBy(a => a, StringComparer.Ordinal)))
                .OrderBy(n => n, StringComparer.Ordinal).ToList();

        /// <summary>
        /// The merge loop as it was: for each adjacent score tie, the first group of the pair absorbs every group of
        /// that score with the same peptides. A group already merged away has score 0, so a pair that begins or ends
        /// with one is not a tie, and the groups after it in the run are never offered to each other.
        /// </summary>
        private static List<string> MergeBeforeRewrite(IReadOnlyList<string[]> blocksOfEachProtein)
        {
            var pg = Model(blocksOfEachProtein);
            for (int i = 0; i < pg.Count - 1; i++)
            {
                if (pg[i].Score == pg[i + 1].Score && pg[i].Score != 0)
                {
                    var pgsWithThisScore = pg.Where(p => p.Score == pg[i].Score).ToList();
                    foreach (var p in pgsWithThisScore)
                    {
                        if (p != pg[i] && p.Blocks == pg[i].Blocks)
                        {
                            pg[i].Accessions.AddRange(p.Accessions);
                            p.Score = 0;
                        }
                    }
                }
            }
            return Names(pg);
        }

        /// <summary>Every set of groups with the same score and the same peptides becomes one group.</summary>
        private static List<string> EveryIndistinguishableSetMerged(IReadOnlyList<string[]> blocksOfEachProtein) =>
            Model(blocksOfEachProtein).GroupBy(g => g.Blocks)
                .Select(set => string.Join("|", set.SelectMany(g => g.Accessions).OrderBy(a => a, StringComparer.Ordinal)))
                .OrderBy(n => n, StringComparer.Ordinal).ToList();

        /// <summary>
        /// The intended difference. P01..P04 tie; P01 and P03 have the same peptides, and so do P02 and P04. The loop
        /// used to merge P03 into P01 and then stop: the pairs P02-P03 and P03-P04 each contain the merged-away P03,
        /// so neither counted as a tie and P02 and P04 stayed separate. They are now merged as well.
        /// </summary>
        [Test]
        public static void TwoInterleavedIndistinguishablePairsAreBothMerged()
        {
            var blocks = new List<string[]>
            {
                new[] { "AAAK", "CCCK" },
                new[] { "DDDK", "EEEK" },
                new[] { "AAAK", "CCCK" },
                new[] { "DDDK", "EEEK" },
            };

            Assert.That(MergeBeforeRewrite(blocks), Is.EqualTo(new[] { "P01|P03", "P02", "P04" }));
            Assert.That(MergedGroups(blocks), Is.EqualTo(new[] { "P01|P03", "P02|P04" }));
        }

        /// <summary>
        /// A group that was merged away is not merged back into, or again out of, the group that absorbed it: three
        /// indistinguishable groups become one group with each protein once.
        /// </summary>
        [Test]
        public static void ThreeIndistinguishableGroupsBecomeOne()
        {
            var blocks = new List<string[]>
            {
                new[] { "AAAK", "CCCK" },
                new[] { "AAAK", "CCCK" },
                new[] { "AAAK", "CCCK" },
            };

            Assert.That(MergeBeforeRewrite(blocks), Is.EqualTo(new[] { "P01|P02|P03" }));
            Assert.That(MergedGroups(blocks), Is.EqualTo(new[] { "P01|P02|P03" }));
        }

        /// <summary>
        /// On random sets of tied groups the engine now merges every indistinguishable set. That equals what the loop
        /// used to give in every case but the one above: where they differ, the old result leaves two groups with the
        /// same score and the same peptides unmerged.
        /// </summary>
        [Test]
        public static void MergingMatchesTheOldLoopExceptWhereItLeftIndistinguishableGroupsApart()
        {
            string[] blockPool = { "AAAK", "CCCK", "DDDK", "EEEK", "FFFK", "GGGK" };
            var random = new Random(3);
            int same = 0, fixedCases = 0;
            for (int trial = 0; trial < 400; trial++)
            {
                var templates = Enumerable.Range(0, random.Next(1, 4))
                    .Select(_ => blockPool.OrderBy(_ => random.Next()).Take(random.Next(1, 4)).ToArray()).ToList();
                var blocks = Enumerable.Range(0, random.Next(2, 9))
                    .Select(_ => templates[random.Next(templates.Count)]).ToList();

                var merged = MergedGroups(blocks);
                var before = MergeBeforeRewrite(blocks);
                string description = string.Join(" ; ", blocks.Select(b => string.Join(",", b)));

                Assert.That(merged, Is.EqualTo(EveryIndistinguishableSetMerged(blocks)), $"trial {trial}: {description}");
                if (before.SequenceEqual(merged))
                {
                    same++;
                }
                else
                {
                    fixedCases++;
                    Assert.That(before.Count, Is.GreaterThan(merged.Count), $"trial {trial}: {description}");
                }
            }

            Assert.That(same, Is.GreaterThan(0));
            Assert.That(fixedCases, Is.GreaterThan(0), "no trial exercised the case the rewrite changes");
        }

        #endregion
    }
}
