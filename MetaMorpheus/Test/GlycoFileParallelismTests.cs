using EngineLayer;
using EngineLayer.DatabaseLoading;
using EngineLayer.GlycoSearch;
using EngineLayer.Util;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
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
