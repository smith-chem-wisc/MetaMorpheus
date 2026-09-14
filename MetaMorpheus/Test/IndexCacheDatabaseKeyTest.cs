using EngineLayer;
using EngineLayer.DatabaseLoading;
using EngineLayer.Indexing;
using NUnit.Framework;
using Omics.Modifications;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace Test
{
    /// <summary>
    /// The database half of the index cache key, and the lookup that consumes it.
    ///
    /// The key used to identify a database by file name and CreationTime. CreationTime does not move
    /// when a file is edited in place, so regenerating or hand-editing a FASTA at the same path left
    /// every byte of the key unchanged and the search reused an index built from the old proteins.
    /// The key now carries a hash of the file's content.
    /// </summary>
    [TestFixture]
    public static class IndexCacheDatabaseKeyTest
    {
        private static string _scratch;

        [OneTimeSetUp]
        public static void SetUp()
        {
            _scratch = Path.Combine(TestContext.CurrentContext.TestDirectory, "IndexCacheDatabaseKeyTest", Guid.NewGuid().ToString("N"));
            Directory.CreateDirectory(_scratch);
        }

        [OneTimeTearDown]
        public static void TearDown()
        {
            try
            {
                if (Directory.Exists(_scratch)) Directory.Delete(_scratch, recursive: true);
            }
            catch (IOException)
            {
                // a leftover scratch directory is not worth failing a test run over
            }
        }

        private static FileInfo WriteDatabase(string relativePath, string contents)
        {
            string full = Path.Combine(_scratch, relativePath);
            Directory.CreateDirectory(Path.GetDirectoryName(full));
            File.WriteAllText(full, contents);
            return new FileInfo(full);
        }

        private const string OneProtein = ">sp|P00001|ONE_HUMAN Protein one\nPEPTIDEKMNNNK\n";
        private const string OtherProtein = ">sp|P00002|TWO_HUMAN Protein two\nCOMPLETELYDIFFERENTSEQVENCEK\n";

        private static IndexingEngine EngineFor(params FileInfo[] databases)
        {
            CommonParameters commonParameters = new CommonParameters();
            var fileSpecificParameters = new List<(string, CommonParameters)> { ("", commonParameters) };

            return new IndexingEngine(
                new List<Protein> { new Protein("MNNNKQQQSTC", "accession") },
                new List<Modification>(), new List<Modification>(),
                null, null, null, 1, DecoyType.None, commonParameters, fileSpecificParameters,
                30000, false, databases.ToList(), TargetContaminantAmbiguity.RemoveContaminant, new List<string>());
        }

        private static string KeyFor(params FileInfo[] databases) => EngineFor(databases).ToString();

        private static string LineStartingWith(string key, string prefix) =>
            key.Split((char)10).Select(line => line.TrimEnd((char)13)).Single(line => line.StartsWith(prefix));

        // ----- the defect this change exists to fix -----

        [Test]
        public static void EditingADatabaseInPlaceChangesTheKey()
        {
            FileInfo database = WriteDatabase("edited/proteins.fasta", OneProtein);
            string before = KeyFor(database);

            File.WriteAllText(database.FullName, OtherProtein);
            string after = KeyFor(database);

            Assert.That(after, Is.Not.EqualTo(before),
                "Same path, different proteins. A shared key here reuses an index built from the old database.");
        }

        [Test]
        public static void CreationTimeAloneWouldNotHaveCaughtThat()
        {
            // Pins the reason the old key failed, so nobody reintroduces it. Editing content leaves
            // CreationTime untouched; only LastWriteTime and Length move, and neither was in the key.
            FileInfo database = WriteDatabase("ctime/proteins.fasta", OneProtein);
            DateTime creationBefore = database.CreationTimeUtc;

            File.WriteAllText(database.FullName, OtherProtein);
            database.Refresh();

            Assert.That(database.CreationTimeUtc, Is.EqualTo(creationBefore),
                "If this ever fails, CreationTime became a content signal and this test can go.");
        }

        [Test]
        public static void TouchingADatabaseWithoutChangingItKeepsTheKey()
        {
            // Content decides the key, not the timestamp. A naive (length, last write time) key would
            // rebuild the index here for nothing.
            FileInfo database = WriteDatabase("touched/proteins.fasta", OneProtein);
            string before = KeyFor(database);

            File.SetLastWriteTimeUtc(database.FullName, DateTime.UtcNow.AddMinutes(5));

            Assert.That(KeyFor(database), Is.EqualTo(before));
        }

        [Test]
        public static void TwoDatabasesWithTheSameNameInDifferentFoldersGiveDifferentKeys()
        {
            // The key names a bare file name, not a path, so content is what separates these.
            FileInfo first = WriteDatabase("samename/a/proteins.fasta", OneProtein);
            FileInfo second = WriteDatabase("samename/b/proteins.fasta", OtherProtein);

            Assert.That(first.Name, Is.EqualTo(second.Name));
            Assert.That(KeyFor(second), Is.Not.EqualTo(KeyFor(first)));
        }

        // ----- reuse must still happen when nothing changed -----

        [Test]
        public static void TheSameDatabaseGivesTheSameKey()
        {
            FileInfo database = WriteDatabase("stable/proteins.fasta", OneProtein);

            Assert.That(KeyFor(database), Is.EqualTo(KeyFor(database)),
                "Equal inputs must reuse the cached index, or every search re-indexes.");
        }

        [Test]
        public static void IdenticalContentUnderDifferentNamesHashesTheSame()
        {
            FileInfo first = WriteDatabase("copies/one.fasta", OneProtein);
            FileInfo second = WriteDatabase("copies/two.fasta", OneProtein);

            Assert.That(IndexingEngine.ContentHash(second), Is.EqualTo(IndexingEngine.ContentHash(first)),
                "The hash is over content, so two copies agree...");
            Assert.That(KeyFor(second), Is.Not.EqualTo(KeyFor(first)),
                "...while the name still separates them in the key.");
        }

        // ----- multiple databases, including contaminants -----

        [Test]
        public static void AddingAContaminantDatabaseChangesTheKey()
        {
            FileInfo sample = WriteDatabase("multi/sample.fasta", OneProtein);
            FileInfo contaminants = WriteDatabase("multi/contaminants.fasta", OtherProtein);

            Assert.That(KeyFor(sample, contaminants), Is.Not.EqualTo(KeyFor(sample)));
        }

        [Test]
        public static void EditingOnlyTheContaminantDatabaseChangesTheKey()
        {
            FileInfo sample = WriteDatabase("multi2/sample.fasta", OneProtein);
            FileInfo contaminants = WriteDatabase("multi2/contaminants.fasta", OtherProtein);
            string before = KeyFor(sample, contaminants);

            File.WriteAllText(contaminants.FullName, OtherProtein + ">sp|P00003|THREE\nMORECONTAMINANTK\n");

            Assert.That(KeyFor(sample, contaminants), Is.Not.EqualTo(before));
        }

        [Test]
        public static void ReorderingDatabasesChangesTheKey()
        {
            // Not cosmetic, and a change from the old sorted key. LoadBioPolymers appends each database's
            // entries in list order and SearchTask slices that list by index to build partitions, so with
            // more than one partition a reordered list puts different proteins in partition k.
            FileInfo first = WriteDatabase("order/a.fasta", OneProtein);
            FileInfo second = WriteDatabase("order/b.fasta", OtherProtein);

            Assert.That(KeyFor(second, first), Is.Not.EqualTo(KeyFor(first, second)));
        }

        [Test]
        public static void TheKeyNamesEachDatabaseAndQualifiesItWithAContentHash()
        {
            FileInfo first = WriteDatabase("render/a.fasta", OneProtein);
            FileInfo second = WriteDatabase("render/b.fasta", OtherProtein);

            string line = LineStartingWith(KeyFor(first, second), "Databases: ");

            Assert.That(line, Is.EqualTo("Databases: " +
                "a.fasta[" + IndexingEngine.ContentHash(first) + "]," +
                "b.fasta[" + IndexingEngine.ContentHash(second) + "]"));
        }

        // ----- degenerate inputs -----

        [Test]
        public static void AMissingDatabaseIsDescribedNotThrown()
        {
            FileInfo missing = new FileInfo(Path.Combine(_scratch, "gone", "never-written.fasta"));

            Assert.That(IndexingEngine.ContentHash(missing), Is.EqualTo("missing"));
            Assert.That(() => KeyFor(missing), Throws.Nothing);
        }

        [Test]
        public static void ADatabaseThatAppearsLaterIsNoticed()
        {
            // FileInfo caches its metadata from construction, and the same instance is reused across
            // candidate folders, so ContentHash refreshes before reading.
            FileInfo database = new FileInfo(Path.Combine(_scratch, "appears", "proteins.fasta"));
            string whileMissing = IndexingEngine.ContentHash(database);

            Directory.CreateDirectory(Path.GetDirectoryName(database.FullName));
            File.WriteAllText(database.FullName, OneProtein);

            Assert.That(whileMissing, Is.EqualTo("missing"));
            Assert.That(IndexingEngine.ContentHash(database), Is.Not.EqualTo("missing"));
        }

        [Test]
        public static void NullDatabaseListsAndEntriesAreDescribed()
        {
            Assert.That(IndexingEngine.DescribeDatabases(null), Is.EqualTo("none"));
            Assert.That(IndexingEngine.DescribeDatabases(new List<FileInfo>()), Is.Empty);
            Assert.That(IndexingEngine.DescribeDatabase(null), Is.EqualTo("none"));
        }

        [Test]
        public static void TheStampSeparatesPathFromLengthFromWriteTime()
        {
            // The memo is keyed on path|length|write time. Drop the separators and "x1.fasta" of length
            // 23 collides with "x1.fasta2" of length 3, handing the second file the first one's hash.
            FileInfo first = WriteDatabase("stamp/x1.fasta", new string('A', 23));
            FileInfo second = WriteDatabase("stamp/x1.fasta2", new string('B', 3));

            Assert.That(first.Length, Is.EqualTo(23));
            Assert.That(second.Length, Is.EqualTo(3));

            DateTime shared = DateTime.UtcNow.AddMinutes(-1);
            File.SetLastWriteTimeUtc(first.FullName, shared);
            File.SetLastWriteTimeUtc(second.FullName, shared);

            Assert.That(IndexingEngine.ContentHash(second), Is.Not.EqualTo(IndexingEngine.ContentHash(first)));
        }

        [Test]
        public static void AnUnreadableDatabaseFallsBackInsteadOfThrowing()
        {
            FileInfo database = WriteDatabase("locked/proteins.fasta", OneProtein);

            string whileLocked;
            using (File.Open(database.FullName, FileMode.Open, FileAccess.Read, FileShare.None))
            {
                whileLocked = IndexingEngine.ContentHash(database);
            }

            Assert.That(whileLocked, Has.Length.EqualTo(16), "the fallback still has to look like a hash");

            // The fallback must not be memoised, or a momentary lock would pin this database's key for
            // the life of the process.
            Assert.That(IndexingEngine.ContentHash(database), Is.Not.EqualTo(whileLocked));
        }

        [Test]
        public static void ContentHashIsSixteenLowercaseHexCharacters()
        {
            FileInfo database = WriteDatabase("hexshape/proteins.fasta", OneProtein);
            string hash = IndexingEngine.ContentHash(database);

            Assert.That(hash, Has.Length.EqualTo(16));
            Assert.That(hash.All(c => (c >= '0' && c <= '9') || (c >= 'a' && c <= 'f')), Is.True,
                "Expected lowercase hex, got: " + hash);
        }

        // ----- the lookup that consumes the key -----

        private static void PlantIndexFolder(FileInfo besideDatabase, IndexingEngine engine)
        {
            string folder = Path.Combine(besideDatabase.DirectoryName, MetaMorpheusTask.IndexFolderName, "2026-09-10-00-00-00");
            Directory.CreateDirectory(folder);

            File.WriteAllText(Path.Combine(folder, MetaMorpheusTask.IndexEngineParamsFileName), engine.ToString());
            File.WriteAllText(Path.Combine(folder, MetaMorpheusTask.PeptideIndexFileName), string.Empty);
            File.WriteAllText(Path.Combine(folder, MetaMorpheusTask.FragmentIndexFileName), string.Empty);
        }

        [Test]
        public static void AnIndexUnderALaterDatabaseIsFoundWhenTheFirstHasNoIndexFolder()
        {
            // The index of a multi-database search is written beside the FIRST database only. Reorder the
            // databases and the cache now sits under a later one; the lookup used to give up at the first
            // database without an index folder and rebuild from scratch.
            FileInfo first = WriteDatabase("lookup/first/a.fasta", OneProtein);
            FileInfo second = WriteDatabase("lookup/second/b.fasta", OtherProtein);

            IndexingEngine engine = EngineFor(first, second);
            PlantIndexFolder(second, engine);

            Assert.That(Directory.Exists(Path.Combine(first.DirectoryName, MetaMorpheusTask.IndexFolderName)), Is.False,
                "the first database must have no index folder for this to test anything");

            string found = MetaMorpheusTask.GetExistingFolderWithIndices(engine,
                new List<DbForTask> { new DbForTask(first.FullName, false), new DbForTask(second.FullName, false) });

            Assert.That(found, Is.Not.Null, "a usable index under the second database was stranded");
            Assert.That(Directory.GetParent(found).Parent.FullName, Is.EqualTo(second.DirectoryName));
        }

        [Test]
        public static void AnIndexUnderTheFirstDatabaseIsStillFound()
        {
            FileInfo first = WriteDatabase("lookupfirst/first/a.fasta", OneProtein);
            FileInfo second = WriteDatabase("lookupfirst/second/b.fasta", OtherProtein);

            IndexingEngine engine = EngineFor(first, second);
            PlantIndexFolder(first, engine);

            string found = MetaMorpheusTask.GetExistingFolderWithIndices(engine,
                new List<DbForTask> { new DbForTask(first.FullName, false), new DbForTask(second.FullName, false) });

            Assert.That(found, Is.Not.Null);
            Assert.That(Directory.GetParent(found).Parent.FullName, Is.EqualTo(first.DirectoryName));
        }

        [Test]
        public static void AnIndexBuiltForADifferentDatabaseIsNotReused()
        {
            FileInfo first = WriteDatabase("mismatch/first/a.fasta", OneProtein);
            FileInfo second = WriteDatabase("mismatch/second/b.fasta", OtherProtein);

            PlantIndexFolder(second, EngineFor(first, second));

            // now the search asks for the second database only, so the planted key must not match
            string found = MetaMorpheusTask.GetExistingFolderWithIndices(EngineFor(second),
                new List<DbForTask> { new DbForTask(second.FullName, false) });

            Assert.That(found, Is.Null);
        }

        [Test]
        public static void NoIndexAnywhereReturnsNull()
        {
            FileInfo first = WriteDatabase("empty/first/a.fasta", OneProtein);
            FileInfo second = WriteDatabase("empty/second/b.fasta", OtherProtein);

            string found = MetaMorpheusTask.GetExistingFolderWithIndices(EngineFor(first, second),
                new List<DbForTask> { new DbForTask(first.FullName, false), new DbForTask(second.FullName, false) });

            Assert.That(found, Is.Null);
        }
    }
}
