using EngineLayer;
using EngineLayer.DatabaseLoading;
using EngineLayer.GlycoSearch;
using EngineLayer.Indexing;
using Nett;
using NUnit.Framework;
using Omics;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Readers;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace Test
{
    /// <summary>
    /// Glyco search sends only each scan's TopN coarse candidates on to glycan matching. When the database is
    /// split into partitions that cut has to be taken over the whole database, not per partition, or the
    /// identifications depend on the partition count -- which RaisePartitionsToFitMemory can change on its own.
    /// </summary>
    [TestFixture]
    public static class GlycoPartitionInvarianceTest
    {
        /// <summary>
        /// Pooling each partition's own TopN cut and cutting again must select exactly the candidates a single
        /// partition would. Scores are drawn from a narrow range so ties at the cut are common, since keeping every
        /// tie with the TopN-th score is where a pooled cut could most easily diverge.
        /// </summary>
        [Test]
        public static void PooledPartitionCuts_SelectTheSameCandidatesAsOnePartition()
        {
            var random = new Random(20260915);
            int[] topNs = { 0, 1, 3, 5, 50 };
            int trialsWithMoreCandidatesThanTopN = 0;

            for (int trial = 0; trial < 2000; trial++)
            {
                int peptideCount = random.Next(1, 300);
                byte cutoff = (byte)random.Next(0, 4);
                int topN = topNs[random.Next(topNs.Length)];
                int partitions = random.Next(2, 9);

                var scores = new byte[peptideCount];
                for (int i = 0; i < peptideCount; i++)
                {
                    scores[i] = (byte)random.Next(0, 12);
                }

                // observed order is arbitrary in the engine, so shuffle it rather than use id order
                List<int> observed = Enumerable.Range(0, peptideCount).OrderBy(_ => random.Next()).ToList();

                var singlePartition = new List<int>();
                GlycoSearchEngine.SelectTopN(observed, scores, cutoff, topN, singlePartition);
                if (topN > 0 && observed.Count(id => scores[id] >= cutoff) > topN)
                {
                    trialsWithMoreCandidatesThanTopN++;
                }

                // contiguous partitions, as the task slices the protein list
                var pooled = new List<(int Partition, int PeptideId, byte Score)>();
                for (int partition = 0; partition < partitions; partition++)
                {
                    int start = partition * peptideCount / partitions;
                    int end = (partition + 1) * peptideCount / partitions;
                    List<int> observedHere = observed.Where(id => id >= start && id < end).ToList();

                    var cutHere = new List<int>();
                    GlycoSearchEngine.SelectTopN(observedHere, scores, cutoff, topN, cutHere);
                    if (cutHere.Count == 0)
                    {
                        continue;
                    }

                    pooled.AddRange(cutHere.Select(id => (partition, id, scores[id])));
                    pooled = GlycoSearchEngine.KeepGlobalTopN(pooled, topN);
                }

                Assert.That(pooled.Select(c => c.PeptideId).OrderBy(id => id), Is.EqualTo(singlePartition.OrderBy(id => id)),
                    $"trial {trial}: {peptideCount} peptides, cutoff {cutoff}, TopN {topN}, {partitions} partitions");
                Assert.That(pooled.Select(c => c.Score), Is.Ordered.Descending, $"trial {trial}: candidates are not highest score first");
            }

            Assert.That(trialsWithMoreCandidatesThanTopN, Is.GreaterThan(500), "premise: the cut must actually remove candidates in most trials");
        }

        [Test]
        public static void KeepGlobalTopN_KeepsEveryTieWithTheTopNthScore()
        {
            var candidates = new List<(int Partition, int PeptideId, byte Score)>
            {
                (0, 0, 5), (0, 1, 9), (1, 0, 7), (1, 1, 7), (2, 0, 7), (2, 1, 3),
            };

            var kept = GlycoSearchEngine.KeepGlobalTopN(candidates, 2);

            Assert.That(kept, Is.EqualTo(new List<(int, int, byte)> { (0, 1, 9), (1, 0, 7), (1, 1, 7), (2, 0, 7) }),
                "the second-highest score is 7, so all three 7s stay, in the order they were added");
        }

        [Test]
        public static void KeepGlobalTopN_NonPositiveTopNKeepsEverything()
        {
            var candidates = new List<(int Partition, int PeptideId, byte Score)> { (0, 0, 1), (1, 0, 4), (1, 1, 2) };

            Assert.That(GlycoSearchEngine.KeepGlobalTopN(candidates, 0).Count, Is.EqualTo(3));
            Assert.That(GlycoSearchEngine.KeepGlobalTopN(candidates, -1).Count, Is.EqualTo(3));
        }

        /// <summary>
        /// The engine-level check on the two rounds: FirstRoundSearch over each partition filling one shared
        /// Candidates array, then Run() per partition matching only that partition's share (the Partition tag).
        /// Every match the engine keeps per scan (up to 10) is compared, not just the best one the task writes, so
        /// an extra or missing candidate shows up even when it would not change the reported result.
        ///
        /// The third run is the premise: cutting each partition on its own, as glyco search did before the two
        /// rounds, must give a different answer on this data. Without that, agreement between the first two
        /// would prove nothing.
        /// </summary>
        [Test]
        public static void GlycoSearchEngine_TwoRoundsOverTwoPartitions_KeepExactlyTheOnePassMatches()
        {
            const int topN = 1;
            CommonParameters one = Toml.ReadFile<GlycoSearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, "GlycoTestData", "GlycoSnip.toml"),
                MetaMorpheusTask.tomlConfig).CommonParameters.CloneWithNewTotalPartitions(1);
            CommonParameters two = one.CloneWithNewTotalPartitions(2);

            string spectra = Path.Combine(TestContext.CurrentContext.TestDirectory, "GlycoTestData", "GlycoPepMix_snip.mzML");
            var file = new MyFileManager(true).LoadFile(spectra, one);
            Ms2ScanWithSpecificMass[] scans = MetaMorpheusTask.GetMs2Scans(file, spectra, one).OrderBy(s => s.PrecursorMass).ToArray();

            // glycoproteins first, then HeLa competitors, so the two partitions hold different kinds of candidate
            var proteins = new List<Protein>();
            foreach (string fasta in new[] { Path.Combine("GlycoTestData", "GlycoProteinFASTA_7proteins.fasta"), Path.Combine("TestData", "hela_snip_for_unitTest.fasta") })
            {
                proteins.AddRange(ProteinDbLoader.LoadProteinFasta(Path.Combine(TestContext.CurrentContext.TestDirectory, fasta), true, DecoyType.None, false, out _));
            }

            // one partition, one pass
            var onePass = new List<GlycoSpectralMatch>[scans.Length];
            IndexingResults whole = BuildIndex(proteins, 0, one);
            MakeGlycoEngine(onePass, scans, whole.PeptideIndex, whole.FragmentIndex, 0, one, topN, null).Run();

            // two partitions, sliced the way GlycoSearchTask slices them
            IndexingResults[] halves = Enumerable.Range(0, 2)
                .Select(p => BuildIndex(proteins.GetRange(p * proteins.Count / 2, (p + 1) * proteins.Count / 2 - p * proteins.Count / 2), p, two))
                .ToArray();

            var twoRounds = new List<GlycoSpectralMatch>[scans.Length];
            var candidates = new List<(int Partition, int PeptideId, byte Score)>[scans.Length];
            for (int p = 0; p < 2; p++)
            {
                MakeGlycoEngine(twoRounds, scans, halves[p].PeptideIndex, halves[p].FragmentIndex, p, two, topN, candidates).FirstRoundSearch();
            }
            for (int p = 0; p < 2; p++)
            {
                // round 2 gets no fragment index, as in the task
                MakeGlycoEngine(twoRounds, scans, halves[p].PeptideIndex, null, p, two, topN, candidates).Run();
            }

            // the old way: each partition cut to TopN and matched on its own
            var perPartitionCut = new List<GlycoSpectralMatch>[scans.Length];
            for (int p = 0; p < 2; p++)
            {
                MakeGlycoEngine(perPartitionCut, scans, halves[p].PeptideIndex, halves[p].FragmentIndex, p, two, topN, null).Run();
            }

            List<string> expected = DescribeMatches(onePass);
            Assert.That(expected.Count, Is.GreaterThan(0), "premise: the one-pass search must match something");
            Assert.That(DescribeMatches(perPartitionCut), Is.Not.EqualTo(expected),
                "premise: a per-partition cut must change the matches on this data, or the next assertion proves nothing");
            Assert.That(DescribeMatches(twoRounds), Is.EqualTo(expected));
        }

        /// <summary>
        /// A partitioned glyco search scores partition p in round 1 and matches it in round 2, and between the two it
        /// holds only peptide ids. If round 2 gets a differently ordered index for p -- a cache written by a run with
        /// another thread count, then rebuilt -- the ids point at other peptides, and the count alone does not show
        /// it: the peptides are the same, only equal-mass ones change places. The fingerprint has to see exactly
        /// that, and also a shared peptide whose two proteins swap.
        /// </summary>
        [Test]
        public static void PeptideOrderFingerprint_SeesEqualMassPeptidesChangePlaces()
        {
            var noMods = new List<Modification>();
            var digestion = new DigestionParams(minPeptideLength: 7);

            // PEPTIDEK and EDITPEPK are the same residues in a different order, so they have the same mass
            List<PeptideWithSetModifications> isomers = new Protein("PEPTIDEKEDITPEPK", "P1").Digest(digestion, noMods, noMods)
                .Where(p => p.BaseSequence.Length == 8).OrderBy(p => p.BaseSequence).ToList();
            Assert.That(isomers.Select(p => p.BaseSequence), Is.EqualTo(new[] { "EDITPEPK", "PEPTIDEK" }));
            Assert.That(isomers[0].MonoisotopicMass, Is.EqualTo(isomers[1].MonoisotopicMass).Within(1e-9), "premise: equal mass");

            // the same peptide from two proteins
            List<PeptideWithSetModifications> shared = new[] { "A1", "B1" }
                .Select(accession => new Protein("PEPTIDEK", accession).Digest(digestion, noMods, noMods).Single()).ToList();

            List<PeptideWithSetModifications> index = isomers.Concat(shared).ToList();
            int fingerprint = IndexingEngine.PeptideOrderFingerprint(index);

            Assert.That(IndexingEngine.PeptideOrderFingerprint(index.ToList()), Is.EqualTo(fingerprint), "same peptides, same order");
            Assert.That(IndexingEngine.PeptideOrderFingerprint(new List<PeptideWithSetModifications> { isomers[1], isomers[0], shared[0], shared[1] }),
                Is.Not.EqualTo(fingerprint), "equal-mass isomers swapped");
            Assert.That(IndexingEngine.PeptideOrderFingerprint(new List<PeptideWithSetModifications> { isomers[0], isomers[1], shared[1], shared[0] }),
                Is.Not.EqualTo(fingerprint), "one sequence from two proteins, proteins swapped");
            Assert.That(IndexingEngine.PeptideOrderFingerprint(index.Take(3).ToList()), Is.Not.EqualTo(fingerprint), "a peptide missing");
        }

        private static IndexingResults BuildIndex(List<Protein> proteins, int partition, CommonParameters parameters)
        {
            return (IndexingResults)new IndexingEngine(proteins, new List<Modification>(), new List<Modification>(), null, null, null, partition,
                DecoyType.Reverse, parameters, null, 30000.0, false, new List<FileInfo>(), TargetContaminantAmbiguity.RemoveContaminant, new List<string>()).Run();
        }

        private static GlycoSearchEngine MakeGlycoEngine(List<GlycoSpectralMatch>[] gsms, Ms2ScanWithSpecificMass[] scans, List<IBioPolymerWithSetMods> peptideIndex,
            FragmentIndex fragmentIndex, int partition, CommonParameters parameters, int topN, List<(int Partition, int PeptideId, byte Score)>[] candidates)
        {
            return new GlycoSearchEngine(gsms, scans, peptideIndex, fragmentIndex, null, partition, parameters, null,
                "OGlycan.gdb", null, GlycoSearchType.OGlycanSearch, topN, 4, true, new List<string>(), candidates);
        }

        /// <summary>One line per kept match: scan, peptide, score and glycan box, sorted. Rank is left out on purpose.</summary>
        private static List<string> DescribeMatches(List<GlycoSpectralMatch>[] gsms)
        {
            return gsms.SelectMany((matches, scan) => (matches ?? new List<GlycoSpectralMatch>()).Select(g =>
                    $"{scan}\t{g.FullSequence}\t{g.Score:R}\t" +
                    string.Join(",", (g.LocalizationGraphs ?? new List<LocalizationGraph>()).Select(l => l.ModBoxId))))
                .OrderBy(s => s, StringComparer.Ordinal)
                .ToList();
        }

        /// <summary>
        /// The whole partitioned path through the task -- both partition loops, FirstRoundSearch filling Candidates,
        /// SecondRoundSearch picking each partition's share by its tag -- must report exactly what one partition
        /// reports. No real run needs 2 partitions on a database this small, so TotalPartitions is set directly; in
        /// real use RaisePartitionsToFitMemory can make the same change on a machine short of memory.
        ///
        /// What this can and cannot catch: the task writes only each scan's best match, and on this data the extra
        /// candidates a per-partition cut lets through never beat it, so a per-partition cut passes here (the engine
        /// test above is the one that catches that). A wrong partition tag, a skipped second loop or a crash in the
        /// wiring fails it. Everything written is compared except Rank, a candidate's position in the cut, which can
        /// legitimately differ between candidates tied on coarse score.
        /// </summary>
        [Test]
        [TestCase(1)]
        [TestCase(3)]
        public static void GlycoSearchTask_TwoPartitions_ReportExactlyWhatOnePartitionReports(int topN)
        {
            string root = MakeEmptyFolder("TestGlycoPartitionInvariance_" + topN);
            try
            {
                string database = WriteGlycoproteinsPlusHela(root);
                string onePartition = RunGlycoSnip(root, Path.Combine(root, "p1"), database, topN, totalPartitions: 1);
                string twoPartitions = RunGlycoSnip(root, Path.Combine(root, "p2"), database, topN, totalPartitions: 2);

                foreach (string file in new[] { "AllPSMs.psmtsv", "oglyco.psmtsv" })
                {
                    List<string> expected = ReadWithoutRank(Path.Combine(onePartition, file));
                    List<string> actual = ReadWithoutRank(Path.Combine(twoPartitions, file));

                    Assert.That(expected.Count, Is.GreaterThan(1), $"premise: {file} must report something to compare");
                    Assert.That(actual, Is.EquivalentTo(expected), $"{file} differs between 1 and 2 partitions at TopN {topN}");
                }
            }
            finally
            {
                Directory.Delete(root, true);
            }
        }

        /// <summary>
        /// Round 2 of a partitioned glyco search holds only peptide ids from round 1. The index cache key does not
        /// include the thread count, but the order of equal-mass peptides does depend on it. So a run can read, in
        /// round 1, a cache written by a run with another thread count; if that cache is then unreadable in round 2,
        /// the rebuild orders equal-mass peptides differently and the ids point at other peptides. That must stop
        /// the search, not report wrong glycopeptides.
        ///
        /// Reproduced exactly: a 1-thread run writes the cache; a 7-thread run reads it in round 1; the cache is
        /// deleted as round 1 finishes, so round 2 rebuilds with 7 threads.
        /// </summary>
        [Test]
        [NonParallelizable] // subscribes to a process-wide task event
        public static void GlycoSearchTask_IndexReorderedBetweenRounds_StopsInsteadOfReportingWrongPeptides()
        {
            string root = MakeEmptyFolder("TestGlycoPartitionIndexReordered");
            try
            {
                string database = WriteGlycoproteinsPlusHela(root);

                // premise: on this database, 1 and 7 threads order at least one partition's index differently
                List<Protein> proteins = ProteinDbLoader.LoadProteinFasta(database, true, DecoyType.None, false, out _);
                bool orderDependsOnThreads = Enumerable.Range(0, 2).Any(p =>
                {
                    List<Protein> slice = proteins.GetRange(p * proteins.Count / 2, (p + 1) * proteins.Count / 2 - p * proteins.Count / 2);
                    return IndexingEngine.PeptideOrderFingerprint(BuildIndex(slice, p, LoadGlycoSnipTask(root, 1, 2, maxThreads: 1).CommonParameters).PeptideIndex)
                        != IndexingEngine.PeptideOrderFingerprint(BuildIndex(slice, p, LoadGlycoSnipTask(root, 1, 2, maxThreads: 7).CommonParameters).PeptideIndex);
                });
                Assert.That(orderDependsOnThreads, "premise: the thread count must change the peptide order on this database");

                RunGlycoSnip(root, Path.Combine(root, "seed"), database, topN: 1, totalPartitions: 2, maxThreads: 1);
                string cache = Path.Combine(root, MetaMorpheusTask.IndexFolderName);
                Assert.That(Directory.Exists(cache), "premise: the seeding run must leave a cache for the next run to read");

                bool deleted = false;
                EventHandler<ProgressEventArgs> deleteCacheAfterRoundOne = (_, e) =>
                {
                    if (!deleted && e.V == "Done with search 2/2!")
                    {
                        Directory.Delete(cache, true);
                        deleted = true;
                    }
                };

                MetaMorpheusTask.OutProgressHandler += deleteCacheAfterRoundOne;
                try
                {
                    GlycoSearchTask task = LoadGlycoSnipTask(root, topN: 1, totalPartitions: 2, maxThreads: 7);
                    string output = Path.Combine(root, "reordered");
                    Directory.CreateDirectory(output);

                    var e = Assert.Throws<MetaMorpheusException>(() => task.RunTask(output,
                        new List<DbForTask> { new DbForTask(database, false) },
                        new List<string> { Path.Combine(TestContext.CurrentContext.TestDirectory, "GlycoTestData", "GlycoPepMix_snip.mzML") }, "Task"));
                    Assert.That(deleted, "premise: the cache must have been deleted between the rounds");
                    Assert.That(e.Message, Does.Contain("changed between the two rounds"));
                    Assert.That(e.Message, Does.Contain(MetaMorpheusTask.IndexFolderName), "must say how to recover");
                }
                finally
                {
                    MetaMorpheusTask.OutProgressHandler -= deleteCacheAfterRoundOne;
                }
            }
            finally
            {
                Directory.Delete(root, true);
            }
        }

        private static string MakeEmptyFolder(string name)
        {
            string folder = Path.Combine(TestContext.CurrentContext.TestDirectory, name);
            if (Directory.Exists(folder))
            {
                Directory.Delete(folder, true);
            }
            Directory.CreateDirectory(folder);
            return folder;
        }

        /// <summary>
        /// The 7 glycoproteins followed by 512 HeLa proteins. The glycoproteins alone are too few for partitioning to
        /// matter: every scan's strong candidates come from them. The HeLa proteins fill the second partition with
        /// competitors, as a real proteome would. Written into <paramref name="folder"/> because the index cache is
        /// written beside the database, and it must not land in GlycoTestData.
        /// </summary>
        private static string WriteGlycoproteinsPlusHela(string folder)
        {
            string database = Path.Combine(folder, "glycoproteins_plus_hela.fasta");
            File.WriteAllText(database,
                File.ReadAllText(Path.Combine(TestContext.CurrentContext.TestDirectory, "GlycoTestData", "GlycoProteinFASTA_7proteins.fasta")).TrimEnd() + Environment.NewLine +
                File.ReadAllText(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "hela_snip_for_unitTest.fasta")));
            return database;
        }

        /// <summary>GlycoSnip.toml with TopN, partitions and threads replaced. Threads have no public setter, hence the text edit.</summary>
        private static GlycoSearchTask LoadGlycoSnipTask(string folder, int topN, int totalPartitions, int maxThreads = 27)
        {
            string toml = File.ReadAllText(Path.Combine(TestContext.CurrentContext.TestDirectory, "GlycoTestData", "GlycoSnip.toml"));
            Assert.That(toml, Does.Contain("MaxThreadsToUsePerFile = 27"), "GlycoSnip.toml changed; update this helper");
            string path = Path.Combine(folder, $"GlycoSnip_{maxThreads}threads.toml");
            File.WriteAllText(path, toml.Replace("MaxThreadsToUsePerFile = 27", "MaxThreadsToUsePerFile = " + maxThreads));

            var task = Toml.ReadFile<GlycoSearchTask>(path, MetaMorpheusTask.tomlConfig);
            task._glycoSearchParameters.GlycoSearchTopNum = topN;
            task.CommonParameters = task.CommonParameters.CloneWithNewTotalPartitions(totalPartitions);
            return task;
        }

        private static string RunGlycoSnip(string folder, string outputFolder, string database, int topN, int totalPartitions, int maxThreads = 27)
        {
            GlycoSearchTask task = LoadGlycoSnipTask(folder, topN, totalPartitions, maxThreads);
            Directory.CreateDirectory(outputFolder);
            new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("Task", task) },
                new List<string> { Path.Combine(TestContext.CurrentContext.TestDirectory, "GlycoTestData", "GlycoPepMix_snip.mzML") },
                new List<DbForTask> { new DbForTask(database, false) }, outputFolder).Run();

            return Path.Combine(outputFolder, "Task");
        }

        /// <summary>Every line of a results file, with the Rank column removed.</summary>
        private static List<string> ReadWithoutRank(string path)
        {
            string[] lines = File.ReadAllLines(path);
            int rankColumn = Array.IndexOf(lines[0].Split('\t'), SpectrumMatchFromTsvHeader.RankLabel);
            Assert.That(rankColumn, Is.GreaterThanOrEqualTo(0), $"no Rank column in {path}");

            return lines.Select(line => string.Join('\t', line.Split('\t').Where((_, column) => column != rankColumn))).ToList();
        }
    }
}
