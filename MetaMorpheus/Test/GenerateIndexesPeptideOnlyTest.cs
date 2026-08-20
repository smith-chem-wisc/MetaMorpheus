using EngineLayer;
using EngineLayer.DatabaseLoading;
using EngineLayer.Indexing;
using NUnit.Framework;
using Omics.Digestion;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace Test
{
    /// <summary>
    /// Regression tests for issue #2412.
    ///
    /// XLSearchTask's second round gets its peptide index from GenerateIndexes_PeptideOnly. That method
    /// used to leave the caller's index null whenever it could not read one from disk -- the branch that
    /// should have rebuilt it was empty, on the assumption that the first round had always just written
    /// one. The null then travelled to the CrosslinkSearchEngine constructor, which does
    /// peptideIndex.Select(...), and surfaced as "Value cannot be null. (Parameter 'source')" with a
    /// stack trace pointing at the crosslink engine rather than at the missing index.
    /// </summary>
    [TestFixture]
    public static class GenerateIndexesPeptideOnlyTest
    {
        private static string _testDirectory;

        [SetUp]
        public static void SetUp()
        {
            _testDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory, "GenerateIndexesPeptideOnly");
            if (Directory.Exists(_testDirectory))
            {
                Directory.Delete(_testDirectory, true);
            }
            Directory.CreateDirectory(_testDirectory);
        }

        [TearDown]
        public static void TearDown()
        {
            if (Directory.Exists(_testDirectory))
            {
                Directory.Delete(_testDirectory, true);
            }
        }

        private static (IndexingEngine engine, List<DbForTask> dbList, List<Protein> proteins) SetUpIndexing(
            string fastaContents, List<Protein> proteins)
        {
            string dbPath = Path.Combine(_testDirectory, "test.fasta");
            File.WriteAllText(dbPath, fastaContents);

            var commonParameters = new CommonParameters(scoreCutoff: 1,
                digestionParams: new DigestionParams(minPeptideLength: 1));
            var fsp = new List<(string fileName, CommonParameters fileSpecificParameters)> { ("", commonParameters) };

            var engine = new IndexingEngine(proteins, new List<Modification>(), new List<Modification>(),
                null, null, null, 1, DecoyType.None, commonParameters, fsp, 30000, false,
                new List<FileInfo> { new FileInfo(dbPath) }, TargetContaminantAmbiguity.RemoveContaminant,
                new List<string>());

            return (engine, new List<DbForTask> { new DbForTask(dbPath, false) }, proteins);
        }

        /// <summary>
        /// No index has been written next to the database, so the read path cannot run at all. Before the
        /// fix this returned with peptideIndex still null.
        /// </summary>
        [Test]
        public static void GenerateIndexesPeptideOnly_WithNoIndexOnDisk_StillReturnsAPeptideIndex()
        {
            var proteins = new List<Protein> { new Protein("MNNNKQQQ", "prot1") };
            var (engine, dbList, allKnownProteins) = SetUpIndexing(">prot1\nMNNNKQQQ\n", proteins);

            List<PeptideWithSetModifications> peptideIndex = null;
            List<int>[] precursorIndex = null;

            new XLSearchTask().GenerateIndexes_PeptideOnly(engine, dbList, ref peptideIndex, ref precursorIndex,
                allKnownProteins, "test");

            Assert.That(peptideIndex, Is.Not.Null,
                "a null peptide index reaches CrosslinkSearchEngine as ArgumentNullException(source)");
            Assert.That(peptideIndex, Is.Not.Empty);
            Assert.That(peptideIndex.Select(p => p.BaseSequence), Does.Contain("MNNNKQQQ"));
        }

        /// <summary>
        /// Exercises the other route into the same branch: an index exists and matches the parameters, so
        /// the read runs, but ReadPeptideIndex rejects a database carrying one accession twice with two
        /// different sequences. That exception used to be swallowed with nothing put in the index's place.
        ///
        /// This is the mechanism the issue reporter's overlapping human/contaminant databases would meet.
        /// Note it hands the duplicated list straight to the method rather than going through database
        /// loading, so it demonstrates the mechanism without establishing that loading passes duplicate
        /// accessions through.
        /// </summary>
        [Test]
        public static void GenerateIndexesPeptideOnly_WhenTheIndexCannotBeRead_RebuildsAndWarns()
        {
            var proteins = new List<Protein> { new Protein("MNNNKQQQ", "prot1") };
            var (engine, dbList, _) = SetUpIndexing(">prot1\nMNNNKQQQ\n", proteins);

            // write a real index next to the database so the read path is reached
            List<PeptideWithSetModifications> written = null;
            List<int>[] fragmentIndex = null;
            List<int>[] precursorIndex = null;
            new XLSearchTask().GenerateIndexes(engine, dbList, ref written, ref fragmentIndex, ref precursorIndex,
                proteins, "test");
            Assert.That(written, Is.Not.Null, "precondition: the index was written");

            // same accession, different sequence -- ReadPeptideIndex throws MetaMorpheusException on this
            var duplicated = new List<Protein>
            {
                new Protein("MNNNKQQQ", "prot1"),
                new Protein("MNNNKQQQPPP", "prot1"),
            };

            string warning = null;
            EventHandler<StringEventArgs> handler = (o, e) => warning ??= e.S;
            MetaMorpheusTask.WarnHandler += handler;

            List<PeptideWithSetModifications> peptideIndex = null;
            List<int>[] precursorIndex2 = null;
            try
            {
                new XLSearchTask().GenerateIndexes_PeptideOnly(engine, dbList, ref peptideIndex, ref precursorIndex2,
                    duplicated, "test");
            }
            finally
            {
                MetaMorpheusTask.WarnHandler -= handler;
            }

            Assert.Multiple(() =>
            {
                Assert.That(peptideIndex, Is.Not.Null, "the unreadable index must be rebuilt, not left null");
                Assert.That(peptideIndex, Is.Not.Empty);
                Assert.That(warning, Is.Not.Null, "the swallowed reason should be reported");
                Assert.That(warning, Does.Contain("rebuilt"));
                Assert.That(warning, Does.Contain("accession"),
                    "the warning should carry the reason the read failed");
            });
        }
    }
}
