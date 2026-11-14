using EngineLayer;
using EngineLayer.DatabaseLoading;
using EngineLayer.FdrAnalysis;
using EngineLayer.ModernSearch;
using EngineLayer.NonSpecificEnzymeSearch;
using NUnit.Framework;
using Omics.Digestion;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Reflection;
using TaskLayer;

namespace Test
{
    /// <summary>
    /// A search task's indexes must be unreachable once its searches are done, during post-search analysis (FDR, PEP,
    /// quantification). Held, they are walked by every full collection that ML.NET's PEP training forces: a Modern search of
    /// the bottom-up nightly reached PEP with 9.6 GB reachable against Classic's 1.9 GB, and spent 63 of its 92 s there in GC.
    /// </summary>
    [TestFixture]
    public static class SearchIndexLetGoTests
    {
        /// <summary>
        /// Two spectra files through a Modern search, and through a semi-specific search (the non-specific loop, N and C
        /// termini): no index a search engine was given may still be reachable when post-search's first FDR finishes.
        /// GenerateIndexes hands the indexes back through ref locals of SearchTask.RunSpecific, and a local whose address is
        /// taken is reported live for the whole method, post-search included, so the last ones searched stayed reachable.
        /// Needs an optimized build: see <see cref="RequireOptimizedSearchTask"/>.
        /// </summary>
        [Test]
        public static void IndexesAreLetGoBeforePostSearch([Values(SearchType.Modern, SearchType.NonSpecific)] SearchType searchType)
        {
            RequireOptimizedSearchTask();

            string folder = Path.Combine(TestContext.CurrentContext.TestDirectory, "SearchIndexLetGo-" + Guid.NewGuid().ToString("N"));
            Directory.CreateDirectory(folder);
            string testData = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData");
            // the index cache goes next to the database: a copy, so the index is built here, not read from an earlier run's cache
            string fasta = Path.Combine(folder, "smalldb.fasta");
            File.Copy(Path.Combine(testData, "smalldb.fasta"), fasta);
            var spectra = new List<string> { Path.Combine(folder, "fileA.mzML"), Path.Combine(folder, "fileB.mzML") };
            spectra.ForEach(path => File.Copy(Path.Combine(testData, "SmallCalibratible_Yeast.mzML"), path));

            var searched = new List<WeakReference>();
            int reachableDuringPostSearch = -1;
            void Watch(object sender, SingleEngineFinishedEventArgs e)
            {
                if (sender is ModernSearchEngine engine) // NonSpecificEnzymeSearchEngine is one too
                {
                    searched.AddRange(IndexesHandedTo(engine).Select(index => new WeakReference(index)));
                }
                else if (sender is FdrAnalysisEngine && reachableDuringPostSearch < 0)
                {
                    reachableDuringPostSearch = StillReachable(searched);
                }
            }

            MetaMorpheusEngine.FinishedSingleEngineHandler += Watch;
            try
            {
                var search = new SearchTask
                {
                    SearchParameters = new SearchParameters { SearchType = searchType },
                    CommonParameters = new CommonParameters(maxThreadsToUsePerFile: 4, digestionParams: searchType == SearchType.Modern
                        ? new DigestionParams(protease: "trypsin", minPeptideLength: 7)
                        : new DigestionParams(protease: "trypsin", minPeptideLength: 7, searchModeType: CleavageSpecificity.Semi)),
                };
                string output = Path.Combine(folder, "out");
                Directory.CreateDirectory(output);
                search.RunTask(output, new List<DbForTask> { new(fasta, false) }, spectra, "search");

                Assert.That(searched, Has.Count.GreaterThanOrEqualTo(2), "a search engine per file");
                Assert.That(reachableDuringPostSearch, Is.EqualTo(0), "indexes still reachable during post-search analysis");
                Assert.That(File.ReadAllLines(Path.Combine(output, "AllPSMs.psmtsv")), Has.Length.GreaterThan(1));
            }
            finally
            {
                MetaMorpheusEngine.FinishedSingleEngineHandler -= Watch;
                Directory.Delete(folder, true);
            }
        }

        /// <summary>
        /// A Modern search that fails part way must keep no index on the task. The GUI runs the same task object again, and the
        /// index a search keeps for its next file is the whole fragment index; its cache write must not be left running either.
        /// Runs in Debug too: it looks once RunSpecific has thrown, when no frame of the task's is left to hold anything.
        /// </summary>
        [Test]
        public static void IndexesAreLetGoWhenTheSearchFails()
        {
            string folder = Path.Combine(TestContext.CurrentContext.TestDirectory, "SearchIndexLetGo-" + Guid.NewGuid().ToString("N"));
            Directory.CreateDirectory(folder);
            string testData = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData");
            // a copy of the database, so the index is built and kept here, not read from an earlier run's cache
            string fasta = Path.Combine(folder, "smalldb.fasta");
            File.Copy(Path.Combine(testData, "smalldb.fasta"), fasta);
            var spectra = new List<string> { Path.Combine(folder, "fileA.mzML"), Path.Combine(folder, "fileB.mzML") };
            spectra.ForEach(path => File.Copy(Path.Combine(testData, "SmallCalibratible_Yeast.mzML"), path));

            var searched = new List<WeakReference>();
            void Watch(object sender, SingleEngineFinishedEventArgs e)
            {
                if (sender is ModernSearchEngine engine)
                {
                    const BindingFlags Instance = BindingFlags.NonPublic | BindingFlags.Instance | BindingFlags.DeclaredOnly;
                    searched.Add(new WeakReference(typeof(ModernSearchEngine).GetProperty("FragmentIndex", Instance).GetValue(engine)));
                }
            }
            // fail as the second file starts, once the first file's index is built and kept for it
            void FailAtSecondFile(object sender, StringEventArgs e)
            {
                if (e.S == spectra[1])
                {
                    throw new InvalidOperationException("the search fails part way");
                }
            }

            MetaMorpheusEngine.FinishedSingleEngineHandler += Watch;
            MetaMorpheusTask.StartingDataFileHandler += FailAtSecondFile;
            try
            {
                var search = new SearchTask
                {
                    SearchParameters = new SearchParameters { SearchType = SearchType.Modern },
                    CommonParameters = new CommonParameters(maxThreadsToUsePerFile: 4, digestionParams: new DigestionParams(protease: "trypsin", minPeptideLength: 7)),
                };
                string output = Path.Combine(folder, "out");
                Directory.CreateDirectory(output);
                Assert.That(() => search.RunTask(output, new List<DbForTask> { new(fasta, false) }, spectra, "search"),
                    Throws.InvalidOperationException.With.Message.EqualTo("the search fails part way"));

                Assert.That(searched, Has.Count.GreaterThanOrEqualTo(1), "the first file was searched");
                Assert.That(StillReachable(searched), Is.EqualTo(0), "the index is still reachable after the search failed");
                GC.KeepAlive(search);
            }
            finally
            {
                MetaMorpheusEngine.FinishedSingleEngineHandler -= Watch;
                MetaMorpheusTask.StartingDataFileHandler -= FailAtSecondFile;
                Directory.Delete(folder, true);
            }
        }

        /// <summary>
        /// The index objects the task handed the engine, as the engine keeps them. A Modern engine keeps the task's peptide list
        /// itself; its fragment index may still be in the cache write, which goes on through post-search. A non-specific engine
        /// copies the peptides into lists of its own, but keeps the task's fragment and precursor indexes, which it writes before
        /// searching.
        /// </summary>
        private static IEnumerable<object> IndexesHandedTo(ModernSearchEngine engine)
        {
            const BindingFlags Instance = BindingFlags.NonPublic | BindingFlags.Instance | BindingFlags.DeclaredOnly;
            if (engine is NonSpecificEnzymeSearchEngine)
            {
                yield return typeof(ModernSearchEngine).GetProperty("FragmentIndex", Instance).GetValue(engine);
                yield return typeof(NonSpecificEnzymeSearchEngine).GetField("PrecursorIndex", Instance).GetValue(engine);
            }
            else
            {
                yield return typeof(ModernSearchEngine).GetField("PeptideIndex", Instance).GetValue(engine);
            }
        }

        /// <summary>
        /// Skips a test that looks while <see cref="SearchTask"/>'s RunSpecific is still running, unless TaskLayer is JIT-optimized
        /// (a Release build). Unoptimized code reports every local and JIT temporary live until the method returns. In a Debug
        /// build that keeps the last index searched reachable through post-search, which runs inside RunSpecific, even after
        /// RunSpecific clears its index locals: the JIT copies those locals (their address is taken) into temporaries for the
        /// search engine's constructor call, and nothing clears those. RunSpecific asks to be optimized
        /// (MethodImplOptions.AggressiveOptimization), but a Debug build's DebuggableAttribute overrides that. A debugger set to
        /// suppress JIT optimization defeats it the same way.
        /// </summary>
        private static void RequireOptimizedSearchTask()
        {
            DebuggableAttribute debuggable = typeof(SearchTask).Assembly.GetCustomAttribute<DebuggableAttribute>();
            if (debuggable != null && debuggable.IsJITOptimizerDisabled)
            {
                Assert.Ignore("Needs an optimized (Release) build of TaskLayer. In a Debug build SearchTask.RunSpecific's frame keeps " +
                    "the last index searched reachable while post-search runs inside it, so this lifetime check cannot pass there.");
            }
        }

        private static int StillReachable(List<WeakReference> references)
        {
            GC.Collect();
            GC.WaitForPendingFinalizers();
            GC.Collect();
            return references.Count(reference => reference.IsAlive);
        }
    }
}
