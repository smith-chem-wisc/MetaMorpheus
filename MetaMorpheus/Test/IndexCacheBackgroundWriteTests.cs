using EngineLayer;
using EngineLayer.DatabaseLoading;
using EngineLayer.Indexing;
using NUnit.Framework;
using Omics;
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
    /// GenerateIndexes can leave the index cache writing on another thread while the search goes ahead with the index in memory.
    /// The cache must come out exactly as a write on the calling thread makes it. Nothing may look for the folder before the write
    /// is finished. A write that fails must still fail the task, at the point where it is waited for.
    /// </summary>
    [TestFixture]
    public static class IndexCacheBackgroundWriteTests
    {
        /// <summary>Minimal task so the index-cache methods can be exercised without running a search.</summary>
        private class IndexProbeTask : MetaMorpheusTask
        {
            public IndexProbeTask() : base(MyTask.Search) { CommonParameters = new CommonParameters(); }

            protected override MyTaskResults RunSpecific(string outputFolder, List<DbForTask> dbFilenameList,
                List<string> currentRawFileList, string taskId, FileSpecificParameters[] fileSettingsList)
                => throw new NotSupportedException("probe only");
        }

        /// <summary>A database folder holding a copy of the small yeast fasta, so each call caches into its own folder.</summary>
        private static (string Folder, List<DbForTask> Databases, List<Protein> Proteins) NewDatabaseFolder()
        {
            string folder = Path.Combine(TestContext.CurrentContext.TestDirectory, "IndexCacheBackground-" + Guid.NewGuid().ToString("N"));
            Directory.CreateDirectory(folder);
            string fasta = Path.Combine(folder, "smalldb.fasta");
            File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "smalldb.fasta"), fasta);
            var proteins = ProteinDbLoader.LoadProteinFasta(fasta, true, DecoyType.Reverse, false, out _);
            return (folder, new List<DbForTask> { new DbForTask(fasta, false) }, proteins);
        }

        private static IndexingEngine MakeEngine(List<Protein> proteins, List<DbForTask> databases, bool generatePrecursorIndex, int maxMissedCleavages = 2)
        {
            var parameters = new CommonParameters(maxThreadsToUsePerFile: 4,
                digestionParams: new DigestionParams(protease: "trypsin", maxMissedCleavages: maxMissedCleavages, minPeptideLength: 7));
            var variable = GlobalVariables.AllModsKnown.Where(m => m.ModificationType == "Common Variable" && m.IdWithMotif == "Oxidation on M").ToList();
            var fixedMods = GlobalVariables.AllModsKnown.Where(m => m.ModificationType == "Common Fixed" && m.IdWithMotif == "Carbamidomethyl on C").ToList();
            return new IndexingEngine(proteins, variable, fixedMods, null, null, null, 0, DecoyType.Reverse, parameters,
                new List<(string, CommonParameters)> { ("", parameters) }, 30000, generatePrecursorIndex,
                databases.Select(p => new FileInfo(p.FilePath)).ToList(), TargetContaminantAmbiguity.RemoveContaminant, new List<string>());
        }

        private static string OnlyCacheFolder(string databaseFolder)
        {
            string[] folders = Directory.GetDirectories(Path.Combine(databaseFolder, "DatabaseIndex"));
            Assert.That(folders, Has.Length.EqualTo(1), "exactly one cache folder");
            return folders[0];
        }

        /// <summary>
        /// The same database indexed twice, once writing on the calling thread and once in the background: every cache file is
        /// byte for byte the same, for a fragment-only index and for one with a precursor index.
        /// </summary>
        [Test]
        public static void ABackgroundWriteLeavesTheSameCacheFilesAsAWriteOnTheCallingThread([Values(false, true)] bool generatePrecursorIndex)
        {
            var (syncFolder, syncDatabases, syncProteins) = NewDatabaseFolder();
            var (backgroundFolder, backgroundDatabases, backgroundProteins) = NewDatabaseFolder();
            try
            {
                List<IBioPolymerWithSetMods> syncPeptides = null, backgroundPeptides = null;
                FragmentIndex syncFragments = null, backgroundFragments = null;
                List<int>[] syncPrecursors = null, backgroundPrecursors = null;

                new IndexProbeTask().GenerateIndexes(MakeEngine(syncProteins, syncDatabases, generatePrecursorIndex), syncDatabases,
                    ref syncPeptides, ref syncFragments, ref syncPrecursors, syncProteins, "sync");

                var task = new IndexProbeTask();
                task.GenerateIndexes(MakeEngine(backgroundProteins, backgroundDatabases, generatePrecursorIndex), backgroundDatabases,
                    ref backgroundPeptides, ref backgroundFragments, ref backgroundPrecursors, backgroundProteins, "background", writeCacheInBackground: true);

                // the index is returned before its cache is written, and it is the one being written
                Assert.That(backgroundPeptides.Select(p => p.FullSequence), Is.EqualTo(syncPeptides.Select(p => p.FullSequence)));
                task.FinishWritingIndexCache();

                string syncCache = OnlyCacheFolder(syncFolder);
                string backgroundCache = OnlyCacheFolder(backgroundFolder);
                var names = new List<string> { MetaMorpheusTask.PeptideIndexFileName, MetaMorpheusTask.FragmentIndexFileName, "DigestionParameters.toml" };
                if (generatePrecursorIndex)
                {
                    names.Add(MetaMorpheusTask.PrecursorIndexFileName);
                }

                Assert.That(Directory.GetFiles(backgroundCache).Select(Path.GetFileName).OrderBy(n => n),
                    Is.EqualTo(Directory.GetFiles(syncCache).Select(Path.GetFileName).OrderBy(n => n)), "the same files, and no partial file left");
                foreach (string name in names)
                {
                    Assert.That(File.ReadAllBytes(Path.Combine(backgroundCache, name)), Is.EqualTo(File.ReadAllBytes(Path.Combine(syncCache, name))), name);
                }
            }
            finally
            {
                Directory.Delete(syncFolder, true);
                Directory.Delete(backgroundFolder, true);
            }
        }

        /// <summary>
        /// The next GenerateIndexes reads the cache the write is leaving, and does not build and write the same index a second
        /// time.
        /// </summary>
        [Test]
        public static void TheNextLookupWaitsForTheWriteAndReadsWhatItWrote()
        {
            var (folder, databases, proteins) = NewDatabaseFolder();
            try
            {
                var task = new IndexProbeTask();
                List<IBioPolymerWithSetMods> built = null, read = null;
                FragmentIndex builtFragments = null, readFragments = null;
                List<int>[] builtPrecursors = null, readPrecursors = null;

                task.GenerateIndexes(MakeEngine(proteins, databases, false), databases, ref built, ref builtFragments, ref builtPrecursors, proteins, "first",
                    writeCacheInBackground: true);
                task.GenerateIndexes(MakeEngine(proteins, databases, false), databases, ref read, ref readFragments, ref readPrecursors, proteins, "second",
                    writeCacheInBackground: true);
                task.FinishWritingIndexCache();

                OnlyCacheFolder(folder);
                Assert.That(read, Is.Not.SameAs(built), "the second call must read the cache, not build again");
                Assert.That(read.Select(p => p.FullSequence), Is.EqualTo(built.Select(p => p.FullSequence)));
                Assert.That(readFragments.PeptideIds, Is.EqualTo(builtFragments.PeptideIds));
            }
            finally
            {
                Directory.Delete(folder, true);
            }
        }

        /// <summary>
        /// Everything a peptide brings to a search and to the results written from it, modifications in their dictionary order
        /// and the mass to the bit, so two indexes compare equal only if no output could tell them apart.
        /// </summary>
        private static List<string> Describe(List<IBioPolymerWithSetMods> peptides)
            => peptides.Cast<PeptideWithSetModifications>().Select(p => string.Join("|",
                p.FullSequence, p.BaseSequence, BitConverter.DoubleToInt64Bits(p.MonoisotopicMass), p.Parent.Accession, p.Parent.IsDecoy,
                p.OneBasedStartResidue, p.OneBasedEndResidue, p.MissedCleavages, p.CleavageSpecificityForFdrCategory, p.PreviousResidue,
                p.NextResidue, p.NumFixedMods, p.DigestionParams.DigestionAgent.Name,
                string.Join(";", p.AllModsOneIsNterminus.Select(m => m.Key + "=" + m.Value.IdWithMotif)))).ToList();

        /// <summary>
        /// The file after the first gets, without waiting for the whole cache to be written, exactly what a lookup of the finished
        /// cache returns: the peptides read from the cache's peptide file, not the built ones (whose modifications can be in a
        /// different order), and the fragment index the cache holds, which is the built one.
        /// </summary>
        [Test]
        public static void TheNextFileGetsWhatALookupOfTheFinishedCacheReturns()
        {
            var (folder, databases, proteins) = NewDatabaseFolder();
            try
            {
                var task = new IndexProbeTask();
                List<IBioPolymerWithSetMods> built = null, next = null, lookedUp = null;
                FragmentIndex builtFragments = null, nextFragments = null, lookedUpFragments = null;
                List<int>[] builtPrecursors = null, nextPrecursors = null, lookedUpPrecursors = null;

                task.GenerateIndexes(MakeEngine(proteins, databases, false), databases, ref built, ref builtFragments, ref builtPrecursors, proteins, "first",
                    writeCacheInBackground: true);
                task.GenerateIndexes(MakeEngine(proteins, databases, false), databases, ref next, ref nextFragments, ref nextPrecursors, proteins, "second",
                    writeCacheInBackground: true);
                task.FinishWritingIndexCache();

                // the lookup the second call stands in for, on a task with nothing kept
                new IndexProbeTask().GenerateIndexes(MakeEngine(proteins, databases, false), databases, ref lookedUp, ref lookedUpFragments,
                    ref lookedUpPrecursors, proteins, "lookup", writeCacheInBackground: true);

                OnlyCacheFolder(folder);
                Assert.That(next, Is.Not.SameAs(built), "the peptides come from the cache file");
                Assert.That(next.Zip(built, ReferenceEquals).Any(same => same), Is.False, "none of them is a built peptide");
                Assert.That(Describe(next), Is.EqualTo(Describe(lookedUp)));
                Assert.That(nextFragments, Is.SameAs(builtFragments), "the fragment index is not read back");
                Assert.That(nextFragments.BinStart, Is.EqualTo(lookedUpFragments.BinStart));
                Assert.That(nextFragments.PeptideIds, Is.EqualTo(lookedUpFragments.PeptideIds));
                Assert.That(nextPrecursors, Is.Null);
            }
            finally
            {
                Directory.Delete(folder, true);
            }
        }

        /// <summary>
        /// The built index is kept only for a lookup that would find the folder being written: one with the same settings, without
        /// a precursor index, before the write has been waited for. Otherwise the lookup goes to disk, as before, and never returns
        /// the kept fragment index.
        /// </summary>
        [Test]
        public static void TheBuiltIndexIsKeptOnlyForALookupThatWouldFindIt([Values("other settings", "precursor index", "after the wait")] string change)
        {
            var (folder, databases, proteins) = NewDatabaseFolder();
            try
            {
                var task = new IndexProbeTask();
                bool precursors = change == "precursor index";
                List<IBioPolymerWithSetMods> built = null, next = null;
                FragmentIndex builtFragments = null, nextFragments = null;
                List<int>[] builtPrecursors = null, nextPrecursors = null;

                task.GenerateIndexes(MakeEngine(proteins, databases, precursors), databases, ref built, ref builtFragments, ref builtPrecursors, proteins, "first",
                    writeCacheInBackground: true);
                if (change == "after the wait")
                {
                    task.FinishWritingIndexCache();
                }
                task.GenerateIndexes(MakeEngine(proteins, databases, precursors, change == "other settings" ? 1 : 2), databases, ref next, ref nextFragments,
                    ref nextPrecursors, proteins, "second", writeCacheInBackground: true);
                task.FinishWritingIndexCache();

                Assert.That(nextFragments, Is.Not.SameAs(builtFragments));
                Assert.That(Directory.GetDirectories(Path.Combine(folder, "DatabaseIndex")),
                    Has.Length.EqualTo(change == "other settings" ? 2 : 1), "other settings build their own cache; the same settings read the first one");
            }
            finally
            {
                Directory.Delete(folder, true);
            }
        }

        /// <summary>
        /// A modern search task waits for its index cache write only when it ends, after post-search, and still leaves the cache
        /// complete for the next task: the next task finds and reads it rather than building another, and no partial file is
        /// left. Two spectra files, so the first task's second file reads back the index its first file built.
        /// </summary>
        [Test]
        public static void AModernSearchTaskEndsWithItsCacheCompleteForTheNextTask()
        {
            var (folder, databases, _) = NewDatabaseFolder();
            try
            {
                string spectra = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "SmallCalibratible_Yeast.mzML");
                string secondSpectra = Path.Combine(folder, "SmallCalibratible_Yeast_copy.mzML");
                File.Copy(spectra, secondSpectra);
                SearchTask Modern() => new SearchTask
                {
                    SearchParameters = new SearchParameters { SearchType = SearchType.Modern },
                    CommonParameters = new CommonParameters(digestionParams: new DigestionParams(protease: "trypsin", minPeptideLength: 7)),
                };
                string output = Path.Combine(folder, "out");
                new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("first", Modern()), ("second", Modern()) },
                    new List<string> { spectra, secondSpectra }, databases, output).Run();

                string cache = OnlyCacheFolder(folder);
                Assert.That(Directory.GetFiles(cache, "*.partial"), Is.Empty, "no write left unfinished");
                Assert.That(File.ReadAllLines(Path.Combine(output, "first", "AllPSMs.psmtsv")), Has.Length.GreaterThan(1));
                Assert.That(File.ReadAllLines(Path.Combine(output, "second", "AllPSMs.psmtsv")), Has.Length.GreaterThan(1));
            }
            finally
            {
                Directory.Delete(folder, true);
            }
        }

        /// <summary>
        /// A write that fails in the background fails the task where the write is waited for: the next GenerateIndexes, before
        /// it looks for or builds anything, or the explicit wait. It fails with the exception itself, not one wrapped around it, and
        /// only once. The real writer still cleans up its partial file.
        /// </summary>
        [Test]
        public static void AFailedBackgroundWriteSurfacesAtTheWait()
        {
            var (folder, databases, proteins) = NewDatabaseFolder();
            try
            {
                var task = new IndexProbeTask();
                string target = Path.Combine(folder, MetaMorpheusTask.FragmentIndexFileName);
                var writeThroughTemporaryFile = typeof(MetaMorpheusTask).GetMethod("WriteThroughTemporaryFile", BindingFlags.NonPublic | BindingFlags.Static);

                // the real writer, failing part way as a full disk would
                task.StartWritingIndexCache(() =>
                {
                    try
                    {
                        writeThroughTemporaryFile.Invoke(null, new object[] { target, (Action<FileStream>)(file =>
                        {
                            file.Write(new byte[64]);
                            throw new IOException("disk full");
                        }) });
                    }
                    catch (TargetInvocationException e)
                    {
                        System.Runtime.ExceptionServices.ExceptionDispatchInfo.Capture(e.InnerException).Throw();
                    }
                });

                List<IBioPolymerWithSetMods> peptides = null;
                FragmentIndex fragments = null;
                List<int>[] precursors = null;
                var thrown = Assert.Throws<IOException>(() => task.GenerateIndexes(MakeEngine(proteins, databases, false), databases,
                    ref peptides, ref fragments, ref precursors, proteins, "next"));
                Assert.That(thrown.Message, Is.EqualTo("disk full"));
                Assert.That(Directory.Exists(Path.Combine(folder, "DatabaseIndex")), Is.False, "the failure surfaces before anything is looked for or built");
                Assert.That(File.Exists(target) || File.Exists(target + ".partial"), Is.False, "the failed write leaves no file");

                Assert.DoesNotThrow(() => task.FinishWritingIndexCache(), "a failure is reported once");

                task.StartWritingIndexCache(() => throw new InvalidOperationException("second failure"));
                Assert.That(() => task.FinishWritingIndexCache(), Throws.InstanceOf<InvalidOperationException>().With.Message.EqualTo("second failure"));
            }
            finally
            {
                Directory.Delete(folder, true);
            }
        }
    }
}
