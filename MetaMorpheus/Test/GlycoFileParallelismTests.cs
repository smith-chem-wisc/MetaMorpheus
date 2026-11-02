using EngineLayer;
using EngineLayer.DatabaseLoading;
using EngineLayer.GlycoSearch;
using EngineLayer.Util;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Threading;
using System.Threading.Tasks;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    public class GlycoFileParallelismTests
    {
        private const long GB = 1_000_000_000;

        [Test]
        public static void DecideSearchesOneFileWithTheWholeBudget()
        {
            var plan = FileParallelism.Decide(fileCount: 1, threadBudget: 32, availableBytes: 500 * GB, bytesPerFile: GB);
            Assert.That(plan.FilesInParallel, Is.EqualTo(1));
            Assert.That(plan.ThreadsPerFile, Is.EqualTo(32));
            Assert.That(plan.LimitedBy, Does.Contain("one spectra file"));
        }

        [Test]
        public static void DecideHonorsTheSequentialSetting()
        {
            var plan = FileParallelism.Decide(fileCount: 10, threadBudget: 64, availableBytes: 500 * GB, bytesPerFile: GB, maximumFilesInParallel: 1);
            Assert.That(plan.FilesInParallel, Is.EqualTo(1));
            Assert.That(plan.ThreadsPerFile, Is.EqualTo(64));
        }

        [Test]
        public static void DecideKeepsAMinimumOfThreadsPerFile()
        {
            // 7 threads cannot give two files the minimum each, so the files are searched one after another.
            var plan = FileParallelism.Decide(fileCount: 2, threadBudget: 7, availableBytes: 500 * GB, bytesPerFile: GB);
            Assert.That(plan.FilesInParallel, Is.EqualTo(1));
            Assert.That(plan.ThreadsPerFile, Is.EqualTo(7));
            Assert.That(plan.LimitedBy, Does.Contain("thread budget"));

            plan = FileParallelism.Decide(fileCount: 10, threadBudget: 60, availableBytes: 500 * GB, bytesPerFile: GB);
            Assert.That(plan.FilesInParallel, Is.EqualTo(10));
            Assert.That(plan.ThreadsPerFile, Is.EqualTo(6));
            Assert.That(plan.LimitedBy, Does.Contain("number of spectra files"));

            plan = FileParallelism.Decide(fileCount: 30, threadBudget: 60, availableBytes: 500 * GB, bytesPerFile: GB);
            Assert.That(plan.FilesInParallel, Is.EqualTo(60 / FileParallelism.MinimumThreadsPerFile));
            Assert.That(plan.ThreadsPerFile, Is.EqualTo(FileParallelism.MinimumThreadsPerFile));
        }

        [Test]
        public static void DecideFitsTheFilesInFreeMemory()
        {
            // The first file is already loaded, so 10 GB free at 80% fits it plus two more 4 GB files.
            var plan = FileParallelism.Decide(fileCount: 8, threadBudget: 64, availableBytes: 10 * GB, bytesPerFile: 4 * GB);
            Assert.That(plan.FilesInParallel, Is.EqualTo(3));
            Assert.That(plan.ThreadsPerFile, Is.EqualTo(21));
            Assert.That(plan.LimitedBy, Is.EqualTo("free memory"));

            // With no memory to spare the first file still runs.
            plan = FileParallelism.Decide(fileCount: 8, threadBudget: 64, availableBytes: 0, bytesPerFile: 4 * GB);
            Assert.That(plan.FilesInParallel, Is.EqualTo(1));
        }

        [Test]
        public static void DecideChargesTheScoringTablesOnceForAllFiles()
        {
            // Two tables of one byte per peptide for each thread; the threads are one budget across every file.
            Assert.That(FileParallelism.ScoringTableBytes(threadBudget: 32, peptideCount: 50_000_000), Is.EqualTo(3_200_000_000L));
            Assert.That(FileParallelism.ScoringTableBytes(threadBudget: 0, peptideCount: 10), Is.EqualTo(20));

            // 10 GB free at 80% is 8 GB; the 3.2 GB of tables come off it once, leaving room for four more 1 GB files beside the
            // first. Charged once per file instead (4.2 GB each), only one more would have fit.
            long tables = FileParallelism.ScoringTableBytes(threadBudget: 32, peptideCount: 50_000_000);
            var plan = FileParallelism.Decide(fileCount: 8, threadBudget: 32, availableBytes: 10 * GB, bytesPerFile: GB, fixedBytes: tables);
            Assert.That(plan.FilesInParallel, Is.EqualTo(5));
            Assert.That(plan.LimitedBy, Is.EqualTo("free memory"));

            // Tables that take the whole budget leave only the first file, even with scans too small to measure.
            plan = FileParallelism.Decide(fileCount: 8, threadBudget: 32, availableBytes: 3 * GB, bytesPerFile: 1, fixedBytes: tables);
            Assert.That(plan.FilesInParallel, Is.EqualTo(1));
            Assert.That(plan.ThreadsPerFile, Is.EqualTo(32));
            Assert.That(plan.LimitedBy, Is.EqualTo("free memory"));
        }

        /// <summary>
        /// Each file searched at once has its own search engine, and each engine holds two copies of the peptide index, so one more
        /// file costs its engine's copies of the index as well as its scans.
        /// </summary>
        [Test]
        public static void BytesPerFileChargesTheEnginesCopiesOfThePeptideIndex()
        {
            // 1 GB of scans is doubled for the results; 50M peptides add 24 bytes each for the two copies.
            Assert.That(FileParallelism.BytesPerFile(scanBytes: GB, peptideCount: 50_000_000), Is.EqualTo(2 * GB + 1_200_000_000L));

            // Scans too small to measure still charge the index copies, and nothing at all still charges a byte, since Decide reads 0 as "unknown".
            Assert.That(FileParallelism.BytesPerFile(scanBytes: 0, peptideCount: 1000), Is.EqualTo(24_000));
            Assert.That(FileParallelism.BytesPerFile(scanBytes: 0, peptideCount: 0), Is.EqualTo(1));

            // 10 GB free at 80% is 8 GB, less 3.2 GB of tables leaves 4.8 GB. At 3.2 GB a file, that is room for one more file beside
            // the first; charging the scans alone (2 GB a file) would have run three at once.
            long tables = FileParallelism.ScoringTableBytes(threadBudget: 32, peptideCount: 50_000_000);
            long perFile = FileParallelism.BytesPerFile(scanBytes: GB, peptideCount: 50_000_000);
            var plan = FileParallelism.Decide(fileCount: 8, threadBudget: 32, availableBytes: 10 * GB, bytesPerFile: perFile, fixedBytes: tables);
            Assert.That(plan.FilesInParallel, Is.EqualTo(2));
            Assert.That(plan.LimitedBy, Is.EqualTo("free memory"));
        }

        /// <summary>
        /// The per-file charge for the index copies must cover what a glyco engine really allocates for them. If the engine stops
        /// copying the index, this fails on the lower bound, and the charge in <see cref="FileParallelism.BytesPerFile"/> should come down.
        /// </summary>
        [Test]
        [NonParallelizable] // measures the managed heap, and the glyco engine writes process-wide glycan state (GlycanBox statics)
        public static void BytesPerFileCoversTheIndexCopiesAGlycoEngineMakes()
        {
            const int peptideCount = 2_000_000;
            var peptide = new PeptideWithSetModifications("PEPTIDE", new Dictionary<string, Omics.Modifications.Modification>());
            var peptideIndex = new List<PeptideWithSetModifications>(peptideCount);
            for (int i = 0; i < peptideCount; i++)
            {
                peptideIndex.Add(peptide);
            }
            var commonParameters = new CommonParameters(dissociationType: DissociationType.HCD, maxThreadsToUsePerFile: 1);
            GlycanSearchSpace searchSpace = GlycanSearchSpace.Build("OGlycan.gdb", null, GlycoSearchType.OGlycanSearch, 3, GlycanBox.DefaultMaximumGlycanBoxMass);

            long before = GC.GetTotalMemory(forceFullCollection: true);
            var engine = new GlycoSearchEngine(new List<GlycoSpectralMatch>[0], new Ms2ScanWithSpecificMass[0], peptideIndex, null, null, 0,
                commonParameters, null, searchSpace, "OGlycan.gdb", null, 30, 3, false, new List<string>());
            long engineBytes = GC.GetTotalMemory(forceFullCollection: true) - before;
            GC.KeepAlive(engine);
            GC.KeepAlive(peptideIndex); // the task holds its index for the whole search; collected here, it would hide one copy
            TestContext.WriteLine("bytes a glyco engine holds for its copies of a " + peptideCount + "-peptide index: " + engineBytes);

            long charged = FileParallelism.BytesPerFile(scanBytes: 0, peptideCount: peptideCount);
            Assert.That(engineBytes, Is.LessThanOrEqualTo(charged), "an engine holds more for its copies of the index than one more file is charged");
            Assert.That(engineBytes, Is.GreaterThanOrEqualTo(16L * peptideCount), "the engine no longer copies the index twice; lower the charge to match");
        }

        [Test]
        public static void DecideHonorsAUserCapAndACappedBudget()
        {
            var plan = FileParallelism.Decide(fileCount: 8, threadBudget: 64, availableBytes: 500 * GB, bytesPerFile: GB, maximumFilesInParallel: 3);
            Assert.That(plan.FilesInParallel, Is.EqualTo(3));
            Assert.That(plan.ThreadsPerFile, Is.EqualTo(21));
            Assert.That(plan.LimitedBy, Does.Contain("MaximumSpectraFilesInParallel = 3"));

            plan = FileParallelism.Decide(fileCount: 4, threadBudget: 0, availableBytes: 500 * GB, bytesPerFile: GB);
            Assert.That(plan.FilesInParallel, Is.EqualTo(1));
            Assert.That(plan.ThreadsPerFile, Is.EqualTo(1));
        }

        /// <summary>
        /// No spectra file may sit unstarted while a file slot is free - that idle tail is what searching files in parallel is
        /// here to remove. Parallel.For claims file indices in chunks that double in size (1, then 2, then 4...), and a worker
        /// takes its whole chunk with it, so a file claimed behind a slow one cannot be picked up by the worker that is free:
        /// with six files at two in parallel a worker claims file 0, then files 1 and 2 together, and file 2 waits for file 1
        /// while the other worker runs out of files of its own and exits. Handing the indices out one at a time, only as a
        /// worker asks for one, is what keeps both slots busy.
        /// </summary>
        [Test]
        public static void ForEachFileLeavesNoFileWaitingWhileAFileSlotIsFree()
        {
            const int fileCount = 6;
            const int filesInParallel = 2;
            var slowFileMayFinish = new ManualResetEventSlim(false);
            var everyOtherFileStarted = new CountdownEvent(fileCount - 1);
            var started = new bool[fileCount];

            var run = Task.Run(() => FileParallelism.ForEachFile(fileCount, filesInParallel, fileIndex =>
            {
                started[fileIndex] = true;
                if (fileIndex == 1)
                {
                    // One slow file, holding one of the two slots until the test lets it go.
                    slowFileMayFinish.Wait(TimeSpan.FromSeconds(60));
                }
                else
                {
                    everyOtherFileStarted.Signal();
                }
            }));

            bool theFreeSlotSearchedEveryOtherFile = everyOtherFileStarted.Wait(TimeSpan.FromSeconds(20));
            string filesLeftWaiting = string.Join(", ", Enumerable.Range(0, fileCount).Where(i => !started[i]));
            slowFileMayFinish.Set();

            Assert.That(run.Wait(TimeSpan.FromSeconds(60)), Is.True, "the run finished");
            Assert.That(theFreeSlotSearchedEveryOtherFile, Is.True,
                $"files not started while the slow file searched: {filesLeftWaiting}");
            Assert.That(started, Is.All.True);
        }

        /// <summary>
        /// Mirrors CloneWithNewTotalPartitions_PreservesEveryOtherSetting: every public setting other than the thread count must
        /// survive the copy, and the fixture must differ from a default instance so a dropped setting would be visible.
        /// </summary>
        [Test]
        public static void CloneWithNewMaxThreadsToUsePerFilePreservesEveryOtherSetting()
        {
            var original = new CommonParameters(
                taskDescriptor: "descriptor",
                dissociationType: DissociationType.ETD,
                ms2childScanDissociationType: DissociationType.EThcD,
                ms3childScanDissociationType: DissociationType.HCD,
                separationType: "CZE",
                doPrecursorDeconvolution: false,
                useProvidedPrecursorInfo: false,
                deconvolutionIntensityRatio: 7,
                deconvolutionMaxAssumedChargeState: 9,
                reportAllAmbiguity: false,
                addCompIons: true,
                totalPartitions: 3,
                qValueThreshold: 0.02,
                pepQValueThreshold: 0.5,
                qValueCutoffForPepCalculation: 0.004,
                scoreCutoff: 4,
                numberOfPeaksToKeepPerWindow: 111,
                minimumAllowedIntensityRatioToBasePeak: 0.02,
                windowWidthThomsons: 12,
                numberOfWindows: 6,
                normalizePeaksAccrossAllWindows: true,
                trimMs1Peaks: true,
                trimMsMsPeaks: false,
                productMassTolerance: new PpmTolerance(11),
                precursorMassTolerance: new PpmTolerance(6),
                productMassTolerance_LowRes: new AbsoluteTolerance(0.44),
                deconvolutionMassTolerance: new PpmTolerance(7),
                maxThreadsToUsePerFile: 2,
                digestionParams: new DigestionParams(protease: "Asp-N", maxMissedCleavages: 3, minPeptideLength: 6),
                listOfModsVariable: new List<(string, string)> { ("Common Biological", "Phosphorylation on S") },
                listOfModsFixed: new List<(string, string)> { ("Common Fixed", "Carbamidomethyl on C") },
                assumeOrphanPeaksAreZ1Fragments: false,
                maxHeterozygousVariants: 2,
                minVariantDepth: 5,
                addTruncations: true,
                precursorDeconParams: new ClassicDeconvolutionParameters(2, 9, 5, 4),
                productDeconParams: new ClassicDeconvolutionParameters(1, 7, 6, 3),
                useMostAbundantPrecursorIntensity: false,
                // a distinct instance without running DIAparameters' constructor, as IndexPartitioningTest does
                diaParameters: (EngineLayer.DIA.DIAparameters)System.Runtime.CompilerServices.RuntimeHelpers.GetUninitializedObject(typeof(EngineLayer.DIA.DIAparameters)),
                fragmentationParams: new FragmentationParams(),
                precursorMassMatchMode: PrecursorMassMatchMode.MostAbundant,
                rtPredictorName: "SSRCalc3");
            typeof(CommonParameters).GetProperty(nameof(CommonParameters.CustomIons))
                .SetValue(original, new List<ProductType> { ProductType.c, ProductType.zDot });

            CommonParameters clone = original.CloneWithNewMaxThreadsToUsePerFile(9);

            Assert.That(clone.MaxThreadsToUsePerFile, Is.EqualTo(9));
            Assert.That(original.MaxThreadsToUsePerFile, Is.EqualTo(2), "the original must not be mutated");

            var mismatches = new List<string>();
            foreach (PropertyInfo property in typeof(CommonParameters).GetProperties(BindingFlags.Public | BindingFlags.Instance))
            {
                if (property.Name == nameof(CommonParameters.MaxThreadsToUsePerFile) || property.GetIndexParameters().Length > 0)
                {
                    continue;
                }
                object a = property.GetValue(original);
                object b = property.GetValue(clone);
                bool equal = a == null || property.PropertyType.IsValueType || property.PropertyType == typeof(string)
                    ? Equals(a, b)
                    : ReferenceEquals(a, b);
                if (!equal)
                {
                    mismatches.Add($"{property.Name}: '{a}' -> '{b}'");
                }
            }
            Assert.That(mismatches, Is.Empty, "clone lost settings: " + string.Join("; ", mismatches));
        }

        /// <summary>
        /// Searching two spectra files at once must write exactly what searching them one after another writes, and the manuscript
        /// prose must say which was done and why.
        /// </summary>
        [Test]
        [NonParallelizable] // the glyco search writes process-wide glycan state (GlycanBox, GlycoSpectralMatch statics)
        public static void SearchingSpectraFilesInParallelMatchesSearchingThemInTurn()
        {
            string testData = Path.Combine(TestContext.CurrentContext.TestDirectory, "GlycoTestData", "N_O_glycoWithFileSpecific");
            string root = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestGlycoFileParallelism");
            if (Directory.Exists(root))
            {
                Directory.Delete(root, true);
            }

            // The spectra alone, without their file-specific tomls, so both files share one index.
            string spectraFolder = Path.Combine(root, "spectra");
            Directory.CreateDirectory(spectraFolder);
            var spectra = Directory.GetFiles(testData, "*.mzML").Select(source =>
            {
                string target = Path.Combine(spectraFolder, Path.GetFileName(source));
                File.Copy(source, target);
                return target;
            }).OrderBy(p => p).ToList();
            Assert.That(spectra, Has.Count.EqualTo(2));

            string nGlycanDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, "GlycoTestData", "NGlycan_ForNoSearch.gdb");
            if (!GlobalVariables.NGlycanDatabasePaths.Any(p => Path.GetFileName(p) == Path.GetFileName(nGlycanDatabase)))
            {
                GlobalVariables.NGlycanDatabasePaths.Add(nGlycanDatabase);
            }
            var database = new List<DbForTask> { new DbForTask(Path.Combine(testData, "FourMucins_NoSigPeps_FASTA.fasta"), false) };

            string Run(string name, int maximumFilesInParallel)
            {
                string output = Path.Combine(root, name);
                Directory.CreateDirectory(output);
                var task = new GlycoSearchTask
                {
                    CommonParameters = new CommonParameters(dissociationType: DissociationType.HCD, ms2childScanDissociationType: DissociationType.EThcD, maxThreadsToUsePerFile: 8),
                    _glycoSearchParameters = new GlycoSearchParameters
                    {
                        OGlycanDatabasefile = "OGlycan.gdb",
                        NGlycanDatabasefile = "NGlycan_ForNoSearch.gdb",
                        GlycoSearchType = GlycoSearchType.N_O_GlycanSearch,
                        MaximumOGlycanAllowed = 2,
                        MaximumSpectraFilesInParallel = maximumFilesInParallel,
                    }
                };
                task.RunTask(output, database, spectra, "task");
                return output;
            }

            string inParallel = Run("parallel", 0);
            string inTurn = Run("sequential", 1);

            string parallelProse = File.ReadAllText(Path.Combine(inParallel, "AutoGeneratedManuscriptProse.txt"));
            string sequentialProse = File.ReadAllText(Path.Combine(inTurn, "AutoGeneratedManuscriptProse.txt"));
            // Threads move between files as they finish, so the prose states the shared budget, not a fixed number per file.
            Assert.That(parallelProse, Does.Contain("spectra files searched in parallel = 2 of 2 (limited by the number of spectra files), each loaded with 4 threads"
                + " and searched with threads shared from a budget of 8 that move from files that finish to files still searching"));
            Assert.That(parallelProse, Does.Contain("memory free after building the index and loading the first spectra file"));
            Assert.That(sequentialProse, Does.Contain("spectra files searched one after another because the task setting MaximumSpectraFilesInParallel = 1"));

            string[] resultFiles = { "AllPSMs.psmtsv", "no_glyco.psmtsv", "_AllProteinGroups.tsv", "seen_no_glyco_localization.tsv", "protein_no_glyco_localization.tsv" };
            foreach (string resultFile in resultFiles)
            {
                string parallelFile = Path.Combine(inParallel, resultFile);
                string sequentialFile = Path.Combine(inTurn, resultFile);
                Assert.That(File.Exists(sequentialFile), Is.True, resultFile + " was not written");
                Assert.That(File.ReadAllBytes(parallelFile), Is.EqualTo(File.ReadAllBytes(sequentialFile)), resultFile + " differs between parallel and sequential searches");
            }

            Directory.Delete(root, true);
        }
    }
}
