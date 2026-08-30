using System.Collections.Generic;
using System.IO;
using System.Linq;
using EngineLayer;
using EngineLayer.DatabaseLoading;
using EngineLayer.Indexing;
using EngineLayer.ModernSearch;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Omics;
using Omics.Fragmentation;
using Omics.Modifications;
using Readers;
using TaskLayer;
using Transcriptomics;
using Transcriptomics.Digestion;
using UsefulProteomicsDatabases;
using UsefulProteomicsDatabases.Transcriptomics;

namespace Test.Transcriptomics
{
    /// <summary>
    /// Modern (index-based) search over a nucleic acid database. The classic-search equivalents live in
    /// TestRnaSearchEngine and SearchTaskWithRna; these assert the indexed path reaches the same answer.
    /// </summary>
    [TestFixture]
    public class TestRnaModernSearch
    {
        private static readonly string SixmerFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics", "TestData", "GUACUG_NegativeMode_Sliced.mzML");

        private static CommonParameters RnaCommonParameters => new CommonParameters(
            dissociationType: DissociationType.CID,
            deconvolutionMaxAssumedChargeState: -20,
            deconvolutionIntensityRatio: 3,
            deconvolutionMassTolerance: new PpmTolerance(20),
            precursorMassTolerance: new PpmTolerance(10),
            productMassTolerance: new PpmTolerance(20),
            scoreCutoff: 5,
            totalPartitions: 1,
            maxThreadsToUsePerFile: 1,
            doPrecursorDeconvolution: true,
            useProvidedPrecursorInfo: false,
            digestionParams: new RnaDigestionParams());

        /// <summary>
        /// The index used to refuse RNA outright: IndexingEngine threw "Not yet implemented for Rna
        /// Digestion" and SearchTask cast the biopolymer list to Protein. It now digests and fragments
        /// oligos like any other biopolymer.
        /// </summary>
        [Test]
        public static void IndexingEngine_IndexesOligos()
        {
            var commonParameters = RnaCommonParameters;
            List<IBioPolymer> targets = [new RNA("GUACUG")];

            var indexEngine = new IndexingEngine(targets, [], [], null, null, null, 0, DecoyType.None,
                commonParameters, null, 30000, false, new List<FileInfo>(),
                TargetContaminantAmbiguity.RemoveContaminant, new List<string>());

            var results = (IndexingResults)indexEngine.Run();

            Assert.That(results.PeptideIndex, Is.Not.Empty, "no oligos were digested into the index");
            Assert.That(results.PeptideIndex.First(), Is.InstanceOf<OligoWithSetMods>());
            Assert.That(results.FragmentIndex.Count(bin => bin != null), Is.GreaterThan(0), "no fragments were binned");
        }

        /// <summary>
        /// Peptide ids must ascend within a fragment bin: the search binary-searches a bin by mass, and
        /// that only holds because the index is sorted by mass. True for oligos as much as for peptides.
        /// </summary>
        [Test]
        public static void OligoIndex_BinsAreSortedByOligoId()
        {
            List<IBioPolymer> targets = [new RNA("GUACUGGUACUGAAUUCC"), new RNA("CAUGCAUGCAUGAAUU")];

            var indexEngine = new IndexingEngine(targets, [], [], null, null, null, 0, DecoyType.Reverse,
                RnaCommonParameters, null, 30000, false, new List<FileInfo>(),
                TargetContaminantAmbiguity.RemoveContaminant, new List<string>());

            var results = (IndexingResults)indexEngine.Run();

            int multiEntryBins = 0;
            foreach (var bin in results.FragmentIndex.Where(b => b != null))
            {
                for (int i = 1; i < bin.Count; i++)
                {
                    Assert.That(bin[i], Is.GreaterThan(bin[i - 1]), "a fragment bin is not ascending by oligo id");
                }
                if (bin.Count > 1)
                {
                    multiEntryBins++;
                }
            }

            Assert.That(multiEntryBins, Is.GreaterThan(0), "no multi-entry bins; the test proves nothing");
            Assert.That(results.PeptideIndex.Select(p => p.MonoisotopicMass), Is.Ordered, "index is not sorted by mass");
        }

        /// <summary>
        /// Non-specific search refuses a nucleic acid database with an explanation, rather than casting
        /// and letting the runtime say it.
        ///
        /// It is built around proteases -- terminal mod placement, the "single" agents, the FDR
        /// categories -- and none of those has an oligo counterpart yet. Without the guard the next line
        /// is `bioPolymerList.Cast&lt;Protein&gt;()`, which throws InvalidCastException naming Protein and
        /// RNA and explaining neither. The distinction the test pins is not that it fails, but that it
        /// fails as a MetaMorpheusException carrying a message a user can act on.
        ///
        /// RunTask writes the failure into results.txt and rethrows, so the exception is observable here.
        /// </summary>
        [Test]
        [NonParallelizable]
        public static void NonSpecificSearch_RefusesANucleicAcidDatabase()
        {
            var original = GlobalVariables.AnalyteType;
            string outputDir = Path.Combine(TestContext.CurrentContext.TestDirectory, "RnaNonSpecificRefused");
            if (Directory.Exists(outputDir)) Directory.Delete(outputDir, true);
            Directory.CreateDirectory(outputDir);

            try
            {
                GlobalVariables.AnalyteType = AnalyteType.Oligo;

                var task = new SearchTask
                {
                    CommonParameters = new CommonParameters(
                        dissociationType: DissociationType.CID,
                        deconvolutionMaxAssumedChargeState: -20,
                        precursorMassTolerance: new PpmTolerance(10),
                        productMassTolerance: new PpmTolerance(20),
                        scoreCutoff: 5,
                        totalPartitions: 1,
                        maxThreadsToUsePerFile: 1,
                        digestionParams: new RnaDigestionParams(
                            rnase: "top-down", fragmentationTerminus: FragmentationTerminus.N)),
                    SearchParameters = new RnaSearchParameters
                    {
                        SearchType = SearchType.NonSpecific,
                        DecoyType = DecoyType.None,
                        MassDiffAcceptorType = MassDiffAcceptorType.Custom,
                        CustomMdac = "Custom interval [-5,5]",
                        DisposeOfFileWhenDone = true
                    }
                };

                var databases = new List<DbForTask>
                {
                    new(Path.Combine(TestContext.CurrentContext.TestDirectory,
                        "Transcriptomics", "TestData", "6mer.fasta"), false)
                };
                var spectra = new List<string> { SixmerFilePath };

                var thrown = Assert.Throws<MetaMorpheusException>(
                    () => task.RunTask(outputDir, databases, spectra, "RefuseNonSpecific"));

                Assert.That(thrown.Message, Does.Contain("only implemented for proteins"));
                Assert.That(thrown.Message, Does.Contain("Classic or Modern"),
                    "the message has to say what to do instead, not just what failed");
            }
            finally
            {
                GlobalVariables.AnalyteType = original;
                if (Directory.Exists(outputDir)) Directory.Delete(outputDir, true);
            }
        }

        /// <summary>
        /// The index cache is skipped for oligos, and this is the branch that decides it.
        ///
        /// The on-disk index only round-trips peptides: OligoWithSetMods is neither [Serializable] nor
        /// able to restore its parent the way SetNonSerializedPeptideInfo does for
        /// PeptideWithSetModifications, so a cached oligo index would come back unusable. GenerateIndexes
        /// therefore builds in memory when the analyte type is Oligo.
        ///
        /// Asserting the index came back is not enough -- it comes back either way. What separates the
        /// two paths is that the cacheable one WRITES, into a "DatabaseIndex" folder beside the database
        /// file. So the database is copied somewhere of its own first and the test asserts nothing
        /// appeared next to it; watching a folder GenerateIndexes was never given would have asserted
        /// nothing at all. Verified against a build with the branch forced the other way, which writes
        /// the folder and fails this.
        /// </summary>
        [Test]
        [NonParallelizable]
        public static void GenerateIndexes_ForOligos_BuildsInMemoryAndCachesNothing()
        {
            var original = GlobalVariables.AnalyteType;
            string databaseFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, "OligoIndexNotCached");
            if (Directory.Exists(databaseFolder)) Directory.Delete(databaseFolder, true);
            Directory.CreateDirectory(databaseFolder);

            try
            {
                GlobalVariables.AnalyteType = AnalyteType.Oligo;

                // A copy of its own, so "was anything cached" is a question about this run only. The
                // shared TestData folder would answer it with somebody else's leftovers.
                string dbPath = Path.Combine(databaseFolder, "20mer1.fasta");
                File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory,
                    "Transcriptomics", "TestData", "20mer1.fasta"), dbPath);

                var commonParameters = RnaCommonParameters;
                var targets = RnaDbLoader.LoadRnaFasta(dbPath, true, DecoyType.None, false, out _)
                    .Cast<IBioPolymer>().ToList();
                Assert.That(targets, Is.Not.Empty, "premise: the fixture database has entries");

                var indexEngine = new IndexingEngine(targets, [], [], null, null, null, 0, DecoyType.None,
                    commonParameters, null, 30000, false, new List<FileInfo>(),
                    TargetContaminantAmbiguity.RemoveContaminant, new List<string>());

                var task = new SearchTask { CommonParameters = commonParameters };
                List<IBioPolymerWithSetMods> oligoIndex = null;
                List<int>[] fragmentIndex = null;
                List<int>[] precursorIndex = null;

                task.GenerateIndexes(indexEngine, new List<DbForTask> { new DbForTask(dbPath, false) },
                    ref oligoIndex, ref fragmentIndex, ref precursorIndex, targets, "TestTaskId");

                Assert.That(oligoIndex, Is.Not.Null.And.Not.Empty, "no oligos came back from the index");
                Assert.That(oligoIndex.First(), Is.InstanceOf<OligoWithSetMods>());
                Assert.That(fragmentIndex.Count(bin => bin != null), Is.GreaterThan(0), "no fragments were binned");

                // The part that distinguishes the two paths.
                Assert.That(Directory.Exists(Path.Combine(databaseFolder, MetaMorpheusTask.IndexFolderName)), Is.False,
                    "an oligo index must not be cached to disk -- it cannot be read back");
            }
            finally
            {
                GlobalVariables.AnalyteType = original;
                if (Directory.Exists(databaseFolder)) Directory.Delete(databaseFolder, true);
            }
        }

        /// <summary>
        /// The payoff: modern search finds the same sixmer classic search finds, and reports it as an
        /// OligoSpectralMatch rather than mislabelling it a PSM.
        /// </summary>
        [Test]
        public static void ModernSearch_FindsSimpleSixmer()
        {
            var commonParameters = RnaCommonParameters;
            var searchParameters = new RnaSearchParameters
            {
                MassDiffAcceptorType = MassDiffAcceptorType.Custom,
                CustomMdac = "Custom interval [-5,5]",
            };

            var dataFile = MsDataFileReader.GetDataFile(SixmerFilePath);
            var ms2Scans = MetaMorpheusTask.GetMs2Scans(dataFile, SixmerFilePath, commonParameters)
                .OrderBy(b => b.PrecursorMass)
                .ToArray();

            MassDiffAcceptor massDiffAcceptor = SearchTask.GetMassDiffAcceptor(commonParameters.PrecursorMassTolerance,
                searchParameters.MassDiffAcceptorType, searchParameters.CustomMdac);

            List<IBioPolymer> targets = [new RNA("GUACUG")];

            var indexEngine = new IndexingEngine(targets, [], [], null, null, null, 0, DecoyType.None,
                commonParameters, null, 30000, false, new List<FileInfo>(),
                TargetContaminantAmbiguity.RemoveContaminant, new List<string>());
            var indexResults = (IndexingResults)indexEngine.Run();

            var osms = new SpectralMatch[ms2Scans.Length];
            new ModernSearchEngine(osms, ms2Scans, indexResults.PeptideIndex, indexResults.FragmentIndex, 0,
                commonParameters, null, massDiffAcceptor, 0, ["modern rna search"]).Run();

            var matches = osms.Where(p => p != null).OrderByDescending(p => p.Score).ToList();

            Assert.That(matches, Is.Not.Empty, "modern search found no matches for the GUACUG sixmer");

            var match = matches.First();
            Assert.That(match, Is.TypeOf<OligoSpectralMatch>(), "match is not an OligoSpectralMatch");
            Assert.That(match.BaseSequence, Is.EqualTo("GUACUG"));
            Assert.That(match.Score, Is.GreaterThan(22), "modern score is below what classic search reaches");
        }
    }
}
