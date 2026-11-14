using EngineLayer;
using EngineLayer.ClassicSearch;
using EngineLayer.DatabaseLoading;
using EngineLayer.FdrAnalysis;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Omics;
using Omics.Modifications;
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
    /// Classic searches of one task keep the database's digestion products for the next spectra file, or calibration
    /// round, instead of digesting the database again. These tests hold the searches to fresh digestion: the same PSMs,
    /// written out, for every file, with the cache, without it, and with too little memory to keep everything; on one
    /// thread and several; and through the search, G-PTM-D and calibration tasks, whose output files must be the same.
    /// </summary>
    [TestFixture]
    [NonParallelizable]
    public static class ClassicSearchDigestionCacheTest
    {
        private static string TestData(string name) => Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", name);

        private sealed record Inputs(List<IBioPolymer> Proteins, List<Modification> Fixed, List<Modification> Variable,
            Ms2ScanWithSpecificMass[][] Files, CommonParameters Parameters, MassDiffAcceptor SearchMode);

        private static Inputs Load(int threads)
        {
            var parameters = new CommonParameters(maxThreadsToUsePerFile: threads, scoreCutoff: 3,
                digestionParams: new DigestionParams(maxMissedCleavages: 2, minPeptideLength: 5, maxModsForPeptides: 2));
            List<Modification> fixedMods = GlobalVariables.AllModsKnown.Where(m => m.ModificationType == "Common Fixed" && m.IdWithMotif == "Carbamidomethyl on C").ToList();
            List<Modification> variableMods = GlobalVariables.AllModsKnown.Where(m => m.ModificationType == "Common Variable" && m.IdWithMotif == "Oxidation on M").ToList();
            List<IBioPolymer> proteins = DatabaseLoadingEngine.LoadBioPolymers("test", new List<DbForTask> { new(TestData("smalldb.fasta"), false) },
                true, DecoyType.Reverse, GlobalVariables.AllModTypesKnown.ToList(), parameters, out _);

            // three "files": every third scan of one file
            string spectra = TestData("SmallCalibratible_Yeast.mzML");
            MsDataFile data = MsDataFileReader.GetDataFile(spectra).LoadAllStaticData();
            List<Ms2ScanWithSpecificMass> scans = MetaMorpheusTask.GetMs2Scans(data, spectra, parameters).ToList();
            var files = Enumerable.Range(0, 3)
                .Select(k => scans.Where((_, i) => i % 3 == k).OrderBy(s => s.GetPrecursorMassForSearch(parameters)).ToArray())
                .ToArray();
            var searchMode = new DotMassDiffAcceptor("", new[] { 0.0, 1.0029 }, new PpmTolerance(10));
            return new Inputs(proteins, fixedMods, variableMods, files, parameters, searchMode);
        }

        /// <summary>Searches every file in order and writes each file's PSMs out as the PSM table would.</summary>
        private static List<string[]> SearchAll(Inputs inputs, ClassicSearchDigestionCache cache, out List<KeyValuePair<(string, string), int>> digestionCounts)
        {
            var written = new List<string[]>();
            digestionCounts = null;
            for (int k = 0; k < inputs.Files.Length; k++)
            {
                var psms = new SpectralMatch[inputs.Files[k].Length];
                var engine = new ClassicSearchEngine(psms, inputs.Files[k], inputs.Variable, inputs.Fixed, null, null, null, inputs.Proteins,
                    inputs.SearchMode, inputs.Parameters, null, null, new List<string>(), false, writeDigestionCounts: k == 0, digestionCache: cache);
                engine.Run();
                if (k == 0)
                {
                    digestionCounts = engine.DigestionCountDictionary.OrderBy(kv => kv.Key.Accession, StringComparer.Ordinal)
                        .ThenBy(kv => kv.Key.BaseSequence, StringComparer.Ordinal).ToList();
                }
                written.Add(psms.Select(p => p == null ? "-" : p.ToString()).ToArray());
            }
            return written;
        }

        [Test]
        [TestCase(1)]
        [TestCase(6)]
        public static void SearchesReadingKeptProductsFindTheSamePsms(int threads)
        {
            Inputs inputs = Load(threads);
            List<string[]> expected = SearchAll(inputs, null, out var expectedCounts);
            Assert.That(expected.Sum(f => f.Count(p => p != "-")), Is.GreaterThan(100), "the files must have PSMs to compare");

            try
            {
                // everything kept
                ClassicSearchDigestionCache.BudgetOverrideBytes = long.MaxValue;
                ClassicSearchDigestionCache cache = ClassicSearchDigestionCache.Create(inputs.Proteins, inputs.Parameters.DigestionParams, inputs.Fixed, inputs.Variable, ClassicSearchDigestionCache.MinimumSearches);
                List<string[]> cached = SearchAll(inputs, cache, out var cachedCounts);
                Assert.That(cache.ProteinsKept, Is.EqualTo(inputs.Proteins.Count));
                AssertSamePsms(expected, cached);
                Assert.That(cachedCounts, Is.EqualTo(expectedCounts), "the first file's digestion counts");

                // about a third kept, the rest digested every time
                long third = cache.EstimatedBytes / 3;
                ClassicSearchDigestionCache.BudgetOverrideBytes = third;
                ClassicSearchDigestionCache partial = ClassicSearchDigestionCache.Create(inputs.Proteins, inputs.Parameters.DigestionParams, inputs.Fixed, inputs.Variable, ClassicSearchDigestionCache.MinimumSearches);
                List<string[]> partlyCached = SearchAll(inputs, partial, out var partialCounts);
                Assert.That(partial.ProteinsKept, Is.GreaterThan(0).And.LessThan(inputs.Proteins.Count));
                Assert.That(partial.EstimatedBytes, Is.GreaterThan(third), "the budget ran out");
                AssertSamePsms(expected, partlyCached);
                Assert.That(partialCounts, Is.EqualTo(expectedCounts));
            }
            finally
            {
                ClassicSearchDigestionCache.BudgetOverrideBytes = null;
            }
        }

        [Test]
        public static void SearchesWithOtherInputsDigestAfresh()
        {
            Inputs inputs = Load(4);
            try
            {
                ClassicSearchDigestionCache.BudgetOverrideBytes = long.MaxValue;
                ClassicSearchDigestionCache cache = ClassicSearchDigestionCache.Create(inputs.Proteins, inputs.Parameters.DigestionParams, inputs.Fixed, inputs.Variable, ClassicSearchDigestionCache.MinimumSearches);

                // equal parameters in another object are served
                var equalParameters = new DigestionParams(maxMissedCleavages: 2, minPeptideLength: 5, maxModsForPeptides: 2);
                Assert.That(cache.Serves(inputs.Proteins, equalParameters, inputs.Fixed, inputs.Variable, null, null));

                // anything that changes digestion is not
                Assert.That(cache.Serves(inputs.Proteins, new DigestionParams(maxMissedCleavages: 1, minPeptideLength: 5, maxModsForPeptides: 2), inputs.Fixed, inputs.Variable, null, null), Is.False);
                Assert.That(cache.Serves(inputs.Proteins.ToList(), inputs.Parameters.DigestionParams, inputs.Fixed, inputs.Variable, null, null), Is.False);
                Assert.That(cache.Serves(inputs.Proteins, inputs.Parameters.DigestionParams, inputs.Fixed.ToList(), inputs.Variable, null, null), Is.False);
                Assert.That(cache.Serves(inputs.Proteins, inputs.Parameters.DigestionParams, inputs.Fixed, inputs.Variable.ToList(), null, null), Is.False);
                Assert.That(cache.Serves(inputs.Proteins, inputs.Parameters.DigestionParams, inputs.Fixed, inputs.Variable, new List<SilacLabel>(), null), Is.False);

                // SILAC and turnover labels, which make digestion build new labelled proteins, get no cache at all
                var label = new SilacLabel('K', 'a', "C-6", 6.020129);
                Assert.That(ClassicSearchDigestionCache.Create(inputs.Proteins, inputs.Parameters.DigestionParams, inputs.Fixed, inputs.Variable, ClassicSearchDigestionCache.MinimumSearches, new List<SilacLabel> { label }), Is.Null);
                Assert.That(ClassicSearchDigestionCache.Create(inputs.Proteins, inputs.Parameters.DigestionParams, inputs.Fixed, inputs.Variable, ClassicSearchDigestionCache.MinimumSearches, null, label, null), Is.Null);

                // too few searches to pay for keeping the products, no cache
                Assert.That(ClassicSearchDigestionCache.Create(inputs.Proteins, inputs.Parameters.DigestionParams, inputs.Fixed, inputs.Variable,
                    ClassicSearchDigestionCache.MinimumSearches - 1), Is.Null);

                // no memory free, no cache
                ClassicSearchDigestionCache.BudgetOverrideBytes = 0;
                Assert.That(ClassicSearchDigestionCache.Create(inputs.Proteins, inputs.Parameters.DigestionParams, inputs.Fixed, inputs.Variable, ClassicSearchDigestionCache.MinimumSearches), Is.Null);
            }
            finally
            {
                ClassicSearchDigestionCache.BudgetOverrideBytes = null;
            }
        }

        [Test]
        public static void TheDefaultBudgetIsAQuarterOfFreeMemory()
        {
            Inputs inputs = Load(1);
            ClassicSearchDigestionCache cache = ClassicSearchDigestionCache.Create(inputs.Proteins, inputs.Parameters.DigestionParams, inputs.Fixed, inputs.Variable, ClassicSearchDigestionCache.MinimumSearches);
            GCMemoryInfo info = GC.GetGCMemoryInfo();
            Assert.That(cache, Is.Not.Null);
            Assert.That(cache.BudgetBytes, Is.GreaterThan(0).And.LessThanOrEqualTo(info.TotalAvailableMemoryBytes / 4));
        }

        private static void AssertSamePsms(List<string[]> expected, List<string[]> actual)
        {
            Assert.That(actual.Count, Is.EqualTo(expected.Count));
            for (int k = 0; k < expected.Count; k++)
            {
                Assert.That(actual[k].Length, Is.EqualTo(expected[k].Length));
                for (int i = 0; i < expected[k].Length; i++)
                {
                    Assert.That(actual[k][i], Is.EqualTo(expected[k][i]), $"file {k}, scan slot {i}");
                }
            }
        }

        #region tasks

        /// <summary>
        /// A task's kept products must be unreachable once its searches are done: during the search task's post-search
        /// analysis (FDR, PEP, quantification), and after the G-PTM-D and calibration tasks return. Held, they double the
        /// live heap that ML.NET's PEP training collects again and again.
        /// </summary>
        [Test]
        public static void KeptProductsAreLetGoWhenTheSearchesAreDone()
        {
            string root = Path.Combine(TestContext.CurrentContext.TestDirectory, "ClassicSearchDigestionCacheTest-release");
            if (Directory.Exists(root))
            {
                Directory.Delete(root, true);
            }
            Directory.CreateDirectory(root);
            var spectra = new List<string> { Path.Combine(root, "fileA.mzML"), Path.Combine(root, "fileB.mzML") };
            spectra.ForEach(path => File.Copy(TestData("SmallCalibratible_Yeast.mzML"), path));
            var databases = new List<DbForTask> { new(TestData("smalldb.fasta"), false) };

            int keptDuringPostSearch = -1;
            void AtFirstFdr(object sender, SingleEngineFinishedEventArgs e)
            {
                if (sender is FdrAnalysisEngine && keptDuringPostSearch < 0)
                {
                    keptDuringPostSearch = ProductsStillReachable();
                }
            }

            MetaMorpheusEngine.FinishedSingleEngineHandler += AtFirstFdr;
            ClassicSearchDigestionCache.BudgetOverrideBytes = long.MaxValue;
            ClassicSearchDigestionCache.MinimumSearchesOverride = 2;
            try
            {
                long served = ClassicSearchDigestionCache.ProteinsServed;
                Directory.CreateDirectory(Path.Combine(root, "search"));
                var search = new SearchTask { CommonParameters = new CommonParameters(maxThreadsToUsePerFile: 4) };
                search.SearchParameters.SearchType = SearchType.Classic;
                search.RunTask(Path.Combine(root, "search"), databases, spectra, "search");
                Assert.That(ClassicSearchDigestionCache.ProteinsServed, Is.GreaterThan(served), "the search task must have used a cache");
                Assert.That(keptDuringPostSearch, Is.EqualTo(0), "products still reachable during post-search analysis");

                MetaMorpheusEngine.FinishedSingleEngineHandler -= AtFirstFdr;
                Directory.CreateDirectory(Path.Combine(root, "gptmd"));
                new GptmdTask { CommonParameters = new CommonParameters(maxThreadsToUsePerFile: 4) }
                    .RunTask(Path.Combine(root, "gptmd"), databases, spectra, "gptmd");
                Assert.That(ProductsStillReachable(), Is.EqualTo(0), "products still reachable after G-PTM-D");

                Directory.CreateDirectory(Path.Combine(root, "calibration"));
                var calibration = new CalibrationTask { CommonParameters = new CommonParameters(maxThreadsToUsePerFile: 4) };
                calibration.RunTask(Path.Combine(root, "calibration"), databases, spectra.Take(1).ToList(), "calibration");
                Assert.That(ProductsStillReachable(), Is.EqualTo(0), "products still reachable after calibration");
                GC.KeepAlive(calibration);
            }
            finally
            {
                MetaMorpheusEngine.FinishedSingleEngineHandler -= AtFirstFdr;
                ClassicSearchDigestionCache.BudgetOverrideBytes = null;
                ClassicSearchDigestionCache.MinimumSearchesOverride = null;
            }
            Directory.Delete(root, true);
        }

        /// <summary>
        /// A calibration task that fails part way must let its kept products go too. The GUI runs the same task object again,
        /// so anything left on it stays reachable until the next run.
        /// </summary>
        [Test]
        public static void KeptProductsAreLetGoWhenCalibrationFails()
        {
            string root = Path.Combine(TestContext.CurrentContext.TestDirectory, "ClassicSearchDigestionCacheTest-failure");
            if (Directory.Exists(root))
            {
                Directory.Delete(root, true);
            }
            Directory.CreateDirectory(root);
            var spectra = new List<string> { Path.Combine(root, "fileA.mzML"), Path.Combine(root, "fileB.mzML") };
            spectra.ForEach(path => File.Copy(TestData("SmallCalibratible_Yeast.mzML"), path));
            var databases = new List<DbForTask> { new(TestData("smalldb.fasta"), false) };

            // fail as the second file starts, after the first file's rounds have made the cache and read from it
            void FailAtSecondFile(object sender, StringEventArgs e)
            {
                if (e.S == spectra[1])
                {
                    throw new InvalidOperationException("calibration fails part way");
                }
            }

            MetaMorpheusTask.StartingDataFileHandler += FailAtSecondFile;
            ClassicSearchDigestionCache.BudgetOverrideBytes = long.MaxValue;
            ClassicSearchDigestionCache.MinimumSearchesOverride = 2;
            try
            {
                long served = ClassicSearchDigestionCache.ProteinsServed;
                var calibration = new CalibrationTask { CommonParameters = new CommonParameters(maxThreadsToUsePerFile: 4) };
                Assert.That(() => calibration.RunTask(root, databases, spectra, "calibration"),
                    Throws.InvalidOperationException.With.Message.EqualTo("calibration fails part way"));
                Assert.That(ClassicSearchDigestionCache.ProteinsServed, Is.GreaterThan(served), "the first file's rounds must have used a cache");
                Assert.That(ProductsStillReachable(), Is.EqualTo(0), "products still reachable after calibration failed");
                GC.KeepAlive(calibration);
            }
            finally
            {
                MetaMorpheusTask.StartingDataFileHandler -= FailAtSecondFile;
                ClassicSearchDigestionCache.BudgetOverrideBytes = null;
                ClassicSearchDigestionCache.MinimumSearchesOverride = null;
            }
            Directory.Delete(root, true);
        }

        private static int ProductsStillReachable()
        {
            GC.Collect();
            GC.WaitForPendingFinalizers();
            GC.Collect();
            return ClassicSearchDigestionCache.ReachableInstances().Sum(cache => cache.ProteinsKept);
        }

        /// <summary>
        /// Two files through the search, G-PTM-D and calibration tasks, with the cache and without it: every output file
        /// that does not report a time must be byte-identical.
        /// </summary>
        [Test]
        public static void TasksWriteTheSameFilesWithAndWithoutKeptProducts()
        {
            string root = Path.Combine(TestContext.CurrentContext.TestDirectory, "ClassicSearchDigestionCacheTest");
            if (Directory.Exists(root))
            {
                Directory.Delete(root, true);
            }
            Directory.CreateDirectory(root);
            string spectraA = Path.Combine(root, "fileA.mzML");
            string spectraB = Path.Combine(root, "fileB.mzML");
            File.Copy(TestData("SmallCalibratible_Yeast.mzML"), spectraA);
            File.Copy(TestData("SmallCalibratible_Yeast.mzML"), spectraB);
            var databases = new List<DbForTask> { new(TestData("smalldb.fasta"), false) };
            var spectra = new List<string> { spectraA, spectraB };

            string withCache = Path.Combine(root, "with");
            string without = Path.Combine(root, "without");
            try
            {
                ClassicSearchDigestionCache.BudgetOverrideBytes = long.MaxValue;
                ClassicSearchDigestionCache.MinimumSearchesOverride = 2;
                long served = ClassicSearchDigestionCache.ProteinsServed;
                RunTasks(withCache, databases, spectra);
                // every protein for the second file of the search and G-PTM-D tasks, and for every calibration round after the first
                Assert.That(ClassicSearchDigestionCache.ProteinsServed - served, Is.GreaterThanOrEqualTo(4L * databases[0].BioPolymerCount.Value));

                ClassicSearchDigestionCache.BudgetOverrideBytes = 0;
                served = ClassicSearchDigestionCache.ProteinsServed;
                RunTasks(without, databases, spectra);
                Assert.That(ClassicSearchDigestionCache.ProteinsServed, Is.EqualTo(served));
            }
            finally
            {
                ClassicSearchDigestionCache.BudgetOverrideBytes = null;
                ClassicSearchDigestionCache.MinimumSearchesOverride = null;
            }

            int compared = 0;
            foreach (string file in Directory.GetFiles(withCache, "*", SearchOption.AllDirectories))
            {
                string name = Path.GetFileName(file);
                if (name is "results.txt" or "allResults.txt" or "AutoGeneratedManuscriptProse.txt" || name.EndsWith(".toml") && file.Contains("Task Settings"))
                {
                    continue; // these report times
                }
                string other = file.Replace(withCache, without);
                Assert.That(File.Exists(other), other);
                Assert.That(File.ReadAllBytes(other), Is.EqualTo(File.ReadAllBytes(file)), name);
                compared++;
            }
            Assert.That(compared, Is.GreaterThan(10));
            Directory.Delete(root, true);
        }

        private static void RunTasks(string outputFolder, List<DbForTask> databases, List<string> spectra)
        {
            foreach (string task in new[] { "search", "gptmd", "calibration" })
            {
                Directory.CreateDirectory(Path.Combine(outputFolder, task));
            }
            var search = new SearchTask { CommonParameters = new CommonParameters(maxThreadsToUsePerFile: 4) };
            search.SearchParameters.SearchType = SearchType.Classic;
            search.RunTask(Path.Combine(outputFolder, "search"), databases, spectra, "search");

            var gptmd = new GptmdTask { CommonParameters = new CommonParameters(maxThreadsToUsePerFile: 4) };
            gptmd.RunTask(Path.Combine(outputFolder, "gptmd"), databases, spectra, "gptmd");

            var calibration = new CalibrationTask { CommonParameters = new CommonParameters(maxThreadsToUsePerFile: 4) };
            calibration.RunTask(Path.Combine(outputFolder, "calibration"), databases, spectra.Take(1).ToList(), "calibration");
        }

        #endregion
    }
}
