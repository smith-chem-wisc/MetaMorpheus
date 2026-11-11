using Chemistry;
using EngineLayer;
using EngineLayer.ClassicSearch;
using EngineLayer.DatabaseLoading;
using EngineLayer.FdrAnalysis;
using EngineLayer.Util;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Omics;
using Omics.Modifications;
using Proteomics.ProteolyticDigestion;
using Readers;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace Test
{
    /// <summary>
    /// Result writing builds .psmtsv rows on several threads (<see cref="OrderedLines"/>). It may not change a byte of any
    /// table: these tests hold the rows to one thread's rows, in five cultures.
    /// </summary>
    [TestFixture]
    [NonParallelizable] // sets the current culture
    public static class ResultWritingTests
    {
        private static string TestData(string name) => Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", name);

        #region OrderedLines

        private static string WriteAll<T>(IReadOnlyList<T> items, Func<T, string> line, int threads)
        {
            using var writer = new StringWriter();
            OrderedLines.Write(writer, items, line, threads);
            return writer.ToString();
        }

        [Test]
        public static void LinesAreWrittenInItemOrderOnAnyNumberOfThreads([Values(1, 2, 8)] int threads)
        {
            // two full batches and a part of a third
            var items = Enumerable.Range(0, 2 * OrderedLines.BatchSize + 517).ToList();
            static string Line(int i) => i.ToString(CultureInfo.InvariantCulture) + "\t" + Math.Sqrt(i).ToString("F5");

            string expected = string.Concat(items.Select(i => Line(i) + Environment.NewLine));
            Assert.That(WriteAll(items, Line, threads), Is.EqualTo(expected));
        }

        [Test]
        public static void LinesAreBuiltInTheCallersCulture()
        {
            var items = Enumerable.Range(0, 3 * OrderedLines.MinimumItemsForThreads).Select(i => i + 0.5).ToList();
            CultureInfo original = CultureInfo.CurrentCulture;
            try
            {
                CultureInfo.CurrentCulture = new CultureInfo("de-DE");
                string text = WriteAll(items, x => x.ToString("F1"), 8);
                Assert.That(text, Does.StartWith("0,5" + Environment.NewLine));
                Assert.That(text, Is.EqualTo(WriteAll(items, x => x.ToString("F1"), 1)));
            }
            finally
            {
                CultureInfo.CurrentCulture = original;
            }

            // the pool threads that built the lines have their own culture back
            var cultures = Enumerable.Range(0, 64).Select(_ => Task.Run(() => CultureInfo.CurrentCulture.Name)).Select(t => t.Result);
            Assert.That(cultures, Has.None.EqualTo("de-DE"));
        }

        [Test]
        public static void AFailingLineIsThrownAsItself()
        {
            var items = Enumerable.Range(0, 3 * OrderedLines.BatchSize).ToList();
            using var writer = new StringWriter();
            int failing = OrderedLines.BatchSize + 7;
            Assert.That(() => OrderedLines.Write(writer, items, i => i == failing ? throw new InvalidOperationException("row " + i) : i.ToString(), 8),
                Throws.InvalidOperationException.With.Message.EqualTo("row " + failing));
        }

        #endregion

        #region rows of a real search

        /// <summary>Minimal task, to reach the protected writer.</summary>
        private sealed class WriterProbe : MetaMorpheusTask
        {
            public WriterProbe() : base(MyTask.Search) { CommonParameters = new CommonParameters(); }

            protected override MyTaskResults RunSpecific(string outputFolder, List<DbForTask> dbFilenameList,
                List<string> currentRawFileList, string taskId, FileSpecificParameters[] fileSettingsList)
                => throw new NotSupportedException("probe only");

            public static void Write(List<SpectralMatch> psms, string path, IReadOnlyDictionary<string, int> mods, bool peptideLevel, int threads)
                => WritePsmsToTsv(psms, path, mods, peptideLevel, threads);
        }

        /// <summary>A Classic search of the yeast test spectra against the small database, with FDR.</summary>
        private static List<SpectralMatch> SearchedPsms()
        {
            var parameters = new CommonParameters(maxThreadsToUsePerFile: 4, scoreCutoff: 3,
                digestionParams: new DigestionParams(maxMissedCleavages: 2, minPeptideLength: 5, maxModsForPeptides: 2));
            List<Modification> fixedMods = GlobalVariables.AllModsKnown.Where(m => m.ModificationType == "Common Fixed" && m.IdWithMotif == "Carbamidomethyl on C").ToList();
            List<Modification> variableMods = GlobalVariables.AllModsKnown.Where(m => m.ModificationType == "Common Variable" && m.IdWithMotif == "Oxidation on M").ToList();
            List<IBioPolymer> proteins = DatabaseLoadingEngine.LoadBioPolymers("test", new List<DbForTask> { new(TestData("smalldb.fasta"), false) },
                true, DecoyType.Reverse, GlobalVariables.AllModTypesKnown.ToList(), parameters, out _);

            string spectra = TestData("SmallCalibratible_Yeast.mzML");
            MsDataFile data = MsDataFileReader.GetDataFile(spectra).LoadAllStaticData();
            Ms2ScanWithSpecificMass[] scans = MetaMorpheusTask.GetMs2Scans(data, spectra, parameters).OrderBy(s => s.GetPrecursorMassForSearch(parameters)).ToArray();
            var searchMode = new DotMassDiffAcceptor("", new[] { 0.0, 1.0029 }, new PpmTolerance(10));
            var psms = new SpectralMatch[scans.Length];
            new ClassicSearchEngine(psms, scans, variableMods, fixedMods, null, null, null, proteins, searchMode, parameters, null, null,
                new List<string>(), false).Run();

            List<SpectralMatch> found = psms.Where(p => p != null).ToList();
            new FdrAnalysisEngine(found, searchMode.NumNotches, parameters, new List<(string, CommonParameters)>(), new List<string>(), doPEP: false).Run();
            return found;
        }

        /// <summary>
        /// The PSM and peptide tables of a real search, written with one thread and with eight, three times over, in five
        /// cultures: every byte the same. Each match is listed several times, so a match's row is built on several threads
        /// at once, which no table asks for.
        /// </summary>
        [Test]
        public static void TablesAreTheSameOnAnyNumberOfThreads()
        {
            List<SpectralMatch> found = SearchedPsms();
            Assert.That(found.Count, Is.GreaterThan(100), "the search must find matches to write");
            // repeated until the table takes more than two batches
            int copies = (2 * OrderedLines.BatchSize + 100) / found.Count + 1;
            List<SpectralMatch> rows = Enumerable.Repeat(found, copies).SelectMany(x => x).ToList();
            IReadOnlyDictionary<string, int> mods = new SearchParameters().ModsToWriteSelection;

            string dir = Path.Combine(TestContext.CurrentContext.TestDirectory, "ResultWritingTests-" + Guid.NewGuid().ToString("N"));
            Directory.CreateDirectory(dir);
            string one = Path.Combine(dir, "one.psmtsv");
            string many = Path.Combine(dir, "many.psmtsv");
            CultureInfo original = CultureInfo.CurrentCulture;
            byte[] english = null;
            try
            {
                foreach (string culture in new[] { "en-US", "de-DE", "fr-FR", "tr-TR", "" })
                {
                    CultureInfo.CurrentCulture = new CultureInfo(culture);
                    foreach (bool peptideLevel in new[] { false, true })
                    {
                        WriterProbe.Write(rows, one, mods, peptideLevel, 1);
                        byte[] expected = File.ReadAllBytes(one);
                        for (int repeat = 0; repeat < 3; repeat++)
                        {
                            WriterProbe.Write(rows, many, mods, peptideLevel, 8);
                            Assert.That(File.ReadAllBytes(many), Is.EqualTo(expected), $"culture '{culture}', peptide level {peptideLevel}, repeat {repeat}");
                        }
                        if (culture == "en-US" && !peptideLevel)
                        {
                            english = expected;
                        }
                        else if (culture == "de-DE" && !peptideLevel)
                        {
                            Assert.That(expected, Is.Not.EqualTo(english), "some columns are written in the current culture, so the cultures must differ");
                        }
                    }
                }
            }
            finally
            {
                CultureInfo.CurrentCulture = original;
                Directory.Delete(dir, true);
            }
        }

        #endregion
    }
}
