using EngineLayer;
using EngineLayer.DatabaseLoading;
using EngineLayer.Indexing;
using MassSpectrometry;
using NUnit.Framework;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace Test
{
    /// <summary>
    /// A cached index folder is selected by <c>CheckFiles</c> before anything reads it. These cover the
    /// header check that happens there, and the crosslink path that has no recovery if it is skipped.
    /// </summary>
    [TestFixture]
    public static class IndexCacheHeaderTest
    {
        private const BindingFlags PrivateStatic = BindingFlags.NonPublic | BindingFlags.Static;

        private static int Constant(string name) =>
            (int)typeof(MetaMorpheusTask).GetField(name, PrivateStatic).GetRawConstantValue();

        private static readonly int FragmentMagic = Constant("FragmentIndexMagic");
        private static readonly int PrecursorMagic = Constant("PrecursorIndexMagic");
        private static readonly int CurrentVersion = Constant("FragmentIndexFormatVersion");

        /// <summary>Minimal task so the index-cache methods can be exercised without running a search.</summary>
        private class IndexProbeTask : MetaMorpheusTask
        {
            public IndexProbeTask() : base(MyTask.XLSearch) { CommonParameters = new CommonParameters(); }

            protected override MyTaskResults RunSpecific(string outputFolder, List<DbForTask> dbFilenameList,
                List<string> currentRawFileList, string taskId, FileSpecificParameters[] fileSettingsList)
                => throw new NotSupportedException("probe only");
        }

        private static List<Protein> Proteins() => new()
        {
            new Protein("MNNNKQQQMNNNKQQQPEPTIDEKMSSSRTTTKAAAWWWKGGGYYYK", "prot1"),
            new Protein("MKVLINGYGTIGKRVADAVSQQDDMKVIGVSKTRPDFEARMALQKGYDLYVAIPK", "prot2"),
        };

        private static IndexingEngine MakeEngine(CommonParameters parameters, List<DbForTask> databases, bool generatePrecursorIndex) =>
            new IndexingEngine(Proteins(), new List<Modification>(), new List<Modification>(), null, null, null, 0,
                DecoyType.Reverse, parameters, null, 30000.0, generatePrecursorIndex,
                databases.Select(p => new FileInfo(p.FilePath)).ToList(),
                TargetContaminantAmbiguity.RemoveContaminant, new List<string>());

        /// <summary>A scratch database directory; the fasta itself is never parsed, only located.</summary>
        private static string NewDatabaseFolder(out List<DbForTask> databases)
        {
            string folder = Path.Combine(TestContext.CurrentContext.TestDirectory, "IndexCacheHeader-" + Guid.NewGuid().ToString("N"));
            Directory.CreateDirectory(folder);
            string fasta = Path.Combine(folder, "probe.fasta");
            File.WriteAllLines(fasta, Proteins().SelectMany(p => new[] { ">sp|" + p.Accession + "|PROBE", p.BaseSequence }));
            databases = new List<DbForTask> { new DbForTask(fasta, false) };
            return folder;
        }

        private static void WriteHeader(string path, int magic, int version)
        {
            var bytes = new byte[16];
            BitConverter.TryWriteBytes(bytes.AsSpan(0), magic);
            BitConverter.TryWriteBytes(bytes.AsSpan(4), version);
            File.WriteAllBytes(path, bytes);
        }

        /// <summary>Overwrite the version word of an existing index file, leaving the rest of it intact.</summary>
        private static void SetVersion(string path, int version)
        {
            byte[] bytes = File.ReadAllBytes(path);
            BitConverter.TryWriteBytes(bytes.AsSpan(4), version);
            File.WriteAllBytes(path, bytes);
        }

        private static string CheckFiles(IndexingEngine engine, string folder) =>
            (string)typeof(MetaMorpheusTask).GetMethod("CheckFiles", PrivateStatic)
                .Invoke(null, new object[] { engine, new DirectoryInfo(folder) });

        /// <summary>
        /// A folder holding the right file names but the wrong payload format has to be a cache miss here.
        /// If it is selected, the reader throws several layers down: GenerateIndexes swallows that and writes
        /// a new folder, but the stale one still matches first on the next run, so the rebuild repeats forever.
        /// </summary>
        [Test]
        [TestCase(true, 0, 0, ExpectedResult = true, TestName = "CheckFiles accepts a current fragment index")]
        [TestCase(true, -1, 0, ExpectedResult = false, TestName = "CheckFiles rejects a previous format version")]
        [TestCase(true, 0, 1, ExpectedResult = false, TestName = "CheckFiles rejects a foreign magic number")]
        [TestCase(false, 0, 0, ExpectedResult = false, TestName = "CheckFiles rejects a header too short to read")]
        public static bool CheckFiles_ValidatesTheFragmentIndexHeader(bool fullHeader, int versionOffset, int magicOffset)
        {
            string folder = NewDatabaseFolder(out List<DbForTask> databases);
            var engine = MakeEngine(new CommonParameters(), databases, generatePrecursorIndex: false);

            string indexFolder = Path.Combine(folder, MetaMorpheusTask.IndexFolderName, "stamp");
            Directory.CreateDirectory(indexFolder);
            File.WriteAllText(Path.Combine(indexFolder, MetaMorpheusTask.IndexEngineParamsFileName), engine.ToString());
            File.WriteAllBytes(Path.Combine(indexFolder, MetaMorpheusTask.PeptideIndexFileName), new byte[8]);

            string fragmentIndexFile = Path.Combine(indexFolder, MetaMorpheusTask.FragmentIndexFileName);
            if (fullHeader)
            {
                WriteHeader(fragmentIndexFile, FragmentMagic + magicOffset, CurrentVersion + versionOffset);
            }
            else
            {
                File.WriteAllBytes(fragmentIndexFile, new byte[4]);
            }

            return CheckFiles(engine, indexFolder) != null;
        }

        /// <summary>
        /// The precursor index carries its own magic, and is only required when the engine asks for one.
        /// </summary>
        [Test]
        public static void CheckFiles_ValidatesThePrecursorIndexHeaderOnlyWhenItIsNeeded()
        {
            string folder = NewDatabaseFolder(out List<DbForTask> databases);
            var withPrecursor = MakeEngine(new CommonParameters(), databases, generatePrecursorIndex: true);
            var withoutPrecursor = MakeEngine(new CommonParameters(), databases, generatePrecursorIndex: false);

            string indexFolder = Path.Combine(folder, MetaMorpheusTask.IndexFolderName, "stamp");
            Directory.CreateDirectory(indexFolder);
            File.WriteAllBytes(Path.Combine(indexFolder, MetaMorpheusTask.PeptideIndexFileName), new byte[8]);
            WriteHeader(Path.Combine(indexFolder, MetaMorpheusTask.FragmentIndexFileName), FragmentMagic, CurrentVersion);
            string precursorIndexFile = Path.Combine(indexFolder, MetaMorpheusTask.PrecursorIndexFileName);
            WriteHeader(precursorIndexFile, PrecursorMagic, CurrentVersion - 1);

            File.WriteAllText(Path.Combine(indexFolder, MetaMorpheusTask.IndexEngineParamsFileName), withPrecursor.ToString());
            Assert.That(CheckFiles(withPrecursor, indexFolder), Is.Null, "a stale precursor index must be a miss when one is required");

            File.WriteAllText(Path.Combine(indexFolder, MetaMorpheusTask.IndexEngineParamsFileName), withoutPrecursor.ToString());
            Assert.That(CheckFiles(withoutPrecursor, indexFolder), Is.Not.Null, "an unread precursor index must not veto the folder");

            WriteHeader(precursorIndexFile, PrecursorMagic, CurrentVersion);
            File.WriteAllText(Path.Combine(indexFolder, MetaMorpheusTask.IndexEngineParamsFileName), withPrecursor.ToString());
            Assert.That(CheckFiles(withPrecursor, indexFolder), Is.Not.Null, "a current precursor index is a hit");
        }

        /// <summary>
        /// The second fragment index is written on demand, so absent is normal and stale is not. It is the one
        /// that has to be caught here: GenerateSecondIndexes reads it with no try/catch of its own.
        /// </summary>
        [Test]
        public static void CheckFiles_RejectsAStaleSecondFragmentIndexButNotAMissingOne()
        {
            string folder = NewDatabaseFolder(out List<DbForTask> databases);
            var engine = MakeEngine(new CommonParameters(), databases, generatePrecursorIndex: false);

            string indexFolder = Path.Combine(folder, MetaMorpheusTask.IndexFolderName, "stamp");
            Directory.CreateDirectory(indexFolder);
            File.WriteAllText(Path.Combine(indexFolder, MetaMorpheusTask.IndexEngineParamsFileName), engine.ToString());
            File.WriteAllBytes(Path.Combine(indexFolder, MetaMorpheusTask.PeptideIndexFileName), new byte[8]);
            WriteHeader(Path.Combine(indexFolder, MetaMorpheusTask.FragmentIndexFileName), FragmentMagic, CurrentVersion);

            Assert.That(CheckFiles(engine, indexFolder), Is.Not.Null, "no second index yet is the normal case");

            string secondFragmentIndexFile = Path.Combine(indexFolder, MetaMorpheusTask.SecondFragmentIndexFileName);
            WriteHeader(secondFragmentIndexFile, FragmentMagic, CurrentVersion - 1);
            Assert.That(CheckFiles(engine, indexFolder), Is.Null, "a stale second index must retire the folder");

            WriteHeader(secondFragmentIndexFile, FragmentMagic, CurrentVersion);
            Assert.That(CheckFiles(engine, indexFolder), Is.Not.Null, "a current second index is a hit");
        }

        /// <summary>
        /// The crosslink sequence end to end, over a cache left by an older format: GenerateIndexes then
        /// GenerateSecondIndexes, as XLSearchTask runs them. Without the header check the stale folder is
        /// selected by both, and the second read throws where nothing catches it.
        /// </summary>
        [Test]
        public static void GenerateSecondIndexes_StaleCachedIndexIsRebuiltRatherThanRead()
        {
            string folder = NewDatabaseFolder(out List<DbForTask> databases);
            var task = new IndexProbeTask();
            var parameters = new CommonParameters();
            var secondParameters = parameters.CloneWithNewDissociationType(DissociationType.ETD);

            FragmentIndex fragmentIndex = null;
            FragmentIndex secondFragmentIndex = null;
            List<PeptideWithSetModifications> peptideIndex = null;
            List<int>[] precursorIndex = null;

            var engine = MakeEngine(parameters, databases, generatePrecursorIndex: false);
            var secondEngine = MakeEngine(secondParameters, databases, generatePrecursorIndex: false);
            task.GenerateIndexes(engine, databases, ref peptideIndex, ref fragmentIndex, ref precursorIndex, Proteins(), "probe");
            task.GenerateSecondIndexes(engine, secondEngine, databases, ref secondFragmentIndex, Proteins(), "probe");

            string indexRoot = Path.Combine(folder, MetaMorpheusTask.IndexFolderName);
            string staleFolder = Directory.GetDirectories(indexRoot).Single();
            string staleSecondIndex = Path.Combine(staleFolder, MetaMorpheusTask.SecondFragmentIndexFileName);
            Assert.That(File.Exists(staleSecondIndex), Is.True, "the first pass must have written a second index to read back");
            int expectedEntries = secondFragmentIndex.EntryCount;
            Assert.That(expectedEntries, Is.GreaterThan(0), "an empty index would make the reread prove nothing");

            // exactly what a build from before the format bump left behind: right names, previous payload
            // version. Only the second index is aged, so without the header check the folder still reads
            // clean in GenerateIndexes, stays the only folder there is, and GenerateSecondIndexes has to
            // read the stale file -- no dependence on the order GetDirectories happens to return.
            SetVersion(staleSecondIndex, CurrentVersion - 1);

            FragmentIndex rebuiltFragmentIndex = null;
            FragmentIndex rebuiltSecondFragmentIndex = null;
            List<PeptideWithSetModifications> rebuiltPeptideIndex = null;
            List<int>[] rebuiltPrecursorIndex = null;

            Assert.DoesNotThrow(() =>
            {
                task.GenerateIndexes(engine, databases, ref rebuiltPeptideIndex, ref rebuiltFragmentIndex, ref rebuiltPrecursorIndex, Proteins(), "probe");
                task.GenerateSecondIndexes(engine, secondEngine, databases, ref rebuiltSecondFragmentIndex, Proteins(), "probe");
            });

            Assert.That(rebuiltSecondFragmentIndex, Is.Not.Null);
            Assert.That(rebuiltSecondFragmentIndex.EntryCount, Is.EqualTo(expectedEntries), "the rebuilt second index must match the one the stale file replaced");

            string freshFolder = Directory.GetDirectories(indexRoot).Single(d => d != staleFolder);
            Assert.That(File.Exists(Path.Combine(freshFolder, MetaMorpheusTask.SecondFragmentIndexFileName)), Is.True,
                "the second index must be written beside the fresh fragment index, not into the retired folder");
        }
    }
}
