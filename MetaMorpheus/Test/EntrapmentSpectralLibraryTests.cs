using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using EngineLayer;
using EngineLayer.DatabaseLoading;
using Nett;
using NUnit.Framework;
using TaskLayer;

namespace Test
{
    /// <summary>
    /// Entrapment is searched as a target, so a spectral library written from a search with entrapment
    /// would carry entrapment spectra. MetaMorpheus writes no library from such a search and says why;
    /// the search itself still runs. Each case also runs the same search without entrapment, so a missing
    /// library is the refusal and not a search that found nothing.
    /// </summary>
    [TestFixture]
    [NonParallelizable]
    public static class EntrapmentSpectralLibraryTests
    {
        private const string RefusalWarning = "No spectral library was written: the database contains entrapment proteins";

        [TestCase(false)]
        [TestCase(true)]
        public static void ASearchWithEntrapmentWritesNoSpectralLibrary(bool withEntrapment)
        {
            string outputFolder = FreshFolder("EntrapmentWriteLibrary");
            string database = Database(@"TestData\hela_snip_for_unitTest.fasta", outputFolder, withEntrapment);
            SearchTask task = Toml.ReadFile<SearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"TestData\SpectralLibrarySearch\SpectralSearchTask.toml"), MetaMorpheusTask.tomlConfig);
            task.SearchParameters.WriteSpectralLibrary = true;
            string raw = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\TaGe_SA_HeLa_04_subset_longestSeq.mzML");

            List<string> warnings = RunCapturingWarnings(() => new EverythingRunnerEngine(
                new List<(string, MetaMorpheusTask)> { ("Search", task) }, new List<string> { raw },
                new List<DbForTask> { new(database, false) }, outputFolder).Run());

            AssertLibraryRefusedExactlyWhenEntrapment(outputFolder, "SpectralLibrary", withEntrapment, warnings);
            Assert.That(File.Exists(Path.Combine(outputFolder, "Search", "AllPSMs.psmtsv")), "the search itself still runs");
            Directory.Delete(outputFolder, true);
        }

        [TestCase(false)]
        [TestCase(true)]
        public static void ASearchWithEntrapmentDoesNotUpdateASpectralLibrary(bool withEntrapment)
        {
            string outputFolder = FreshFolder("EntrapmentUpdateLibrary");
            string database = Database(@"TestData\hela_snip_for_unitTest.fasta", outputFolder, withEntrapment);
            SearchTask task = Toml.ReadFile<SearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"TestData\SpectralLibrarySearch\SpectralSearchTask.toml"), MetaMorpheusTask.tomlConfig);
            task.SearchParameters.UpdateSpectralLibrary = true;
            string raw = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\TaGe_SA_HeLa_04_subset_longestSeq.mzML");
            string library = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SpectralLibrarySearch\SpectralLibrary.msp");

            List<string> warnings = RunCapturingWarnings(() => new EverythingRunnerEngine(
                new List<(string, MetaMorpheusTask)> { ("Search", task) }, new List<string> { raw },
                new List<DbForTask> { new(library, false), new(database, false) }, outputFolder).Run());

            AssertLibraryRefusedExactlyWhenEntrapment(outputFolder, "updateSpectralLibrary", withEntrapment, warnings);
            Directory.Delete(outputFolder, true);
        }

        [TestCase(false)]
        [TestCase(true)]
        public static void ACrosslinkSearchWithEntrapmentWritesNoSpectralLibrary(bool withEntrapment)
        {
            string outputFolder = FreshFolder("EntrapmentXlLibrary");
            string database = Database(@"XlTestData\RibosomeGO.fasta", outputFolder, withEntrapment);
            var task = new XLSearchTask
            {
                XlSearchParameters = new XlSearchParameters { WriteSpectralLibrary = true, CrosslinkAtCleavageSite = true }
            };
            string raw = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData\2017-11-21_XL_DSSO_Ribosome_RT60min_28800-28898.mzML");

            List<string> warnings = RunCapturingWarnings(() =>
                task.RunTask(outputFolder, new List<DbForTask> { new(database, false) }, new List<string> { raw }, "test"));

            AssertLibraryRefusedExactlyWhenEntrapment(outputFolder, "SpectralLibrary", withEntrapment, warnings);
            Directory.Delete(outputFolder, true);
        }

        private static void AssertLibraryRefusedExactlyWhenEntrapment(string outputFolder, string libraryPrefix,
            bool withEntrapment, List<string> warnings)
        {
            string[] libraries = Directory.GetFiles(outputFolder, libraryPrefix + "_*.msp", SearchOption.AllDirectories);
            bool warned = warnings.Any(w => w.Contains(RefusalWarning));
            if (withEntrapment)
            {
                Assert.That(libraries, Is.Empty);
                Assert.That(warned, "the refusal must say why no library was written");
            }
            else
            {
                Assert.That(libraries, Is.Not.Empty, "without entrapment the same search writes a library");
                Assert.That(warned, Is.False);
            }
        }

        /// <summary>
        /// Writes a copy of a test FASTA in which, when <paramref name="withEntrapment"/> is set, every second
        /// entry's accession gains the "Random_" prefix that mzLib's loader reads as entrapment.
        /// </summary>
        private static string Database(string relativePath, string outputFolder, bool withEntrapment)
        {
            string copy = Path.Combine(outputFolder, Path.GetFileName(relativePath));
            int entry = 0;
            IEnumerable<string> lines = File.ReadLines(Path.Combine(TestContext.CurrentContext.TestDirectory, relativePath))
                .Select(line => line.StartsWith(">sp|") && entry++ % 2 == 1 && withEntrapment
                    ? ">sp|Random_" + line.Substring(">sp|".Length)
                    : line);
            File.WriteAllLines(copy, lines);
            return copy;
        }

        private static string FreshFolder(string name)
        {
            string folder = Path.Combine(TestContext.CurrentContext.TestDirectory, name);
            if (Directory.Exists(folder))
                Directory.Delete(folder, true);
            Directory.CreateDirectory(folder);
            return folder;
        }

        private static List<string> RunCapturingWarnings(Action run)
        {
            var warnings = new List<string>();
            EventHandler<StringEventArgs> handler = (_, e) => { lock (warnings) warnings.Add(e.S); };
            MetaMorpheusTask.WarnHandler += handler;
            try
            {
                run();
            }
            finally
            {
                MetaMorpheusTask.WarnHandler -= handler;
            }
            return warnings;
        }
    }
}
