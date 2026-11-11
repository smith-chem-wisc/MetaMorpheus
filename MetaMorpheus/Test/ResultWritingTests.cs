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
    /// Result writing builds .psmtsv rows on several threads (<see cref="OrderedLines"/>), prepares each individual file's FDR
    /// one file ahead of writing its tables (<see cref="PostSearchAnalysisTask.PrepareOneAhead{TItem,TPrepared}"/>), and
    /// resolves a column's one-hypothesis case without building a list (<see cref="PsmTsvWriter"/>.Resolve). None of it may
    /// change a byte of any table: these tests hold the rows to one thread's rows in five cultures, and Resolve to a verbatim
    /// copy of the code it replaced.
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

        #region PrepareOneAhead

        [Test]
        public static void PreparationsRunOneAheadOfTheWritesAndInOrder()
        {
            var events = new ConcurrentQueue<string>();
            int preparing = 0;
            int mostPreparing = 0;
            PostSearchAnalysisTask.PrepareOneAhead(Enumerable.Range(0, 6).ToList(),
                i =>
                {
                    int now = Interlocked.Increment(ref preparing);
                    InterlockedMax(ref mostPreparing, now);
                    events.Enqueue("prepare start " + i);
                    Thread.Sleep(10);
                    Interlocked.Decrement(ref preparing);
                    return i * 10;
                },
                (i, prepared) =>
                {
                    Assert.That(prepared, Is.EqualTo(i * 10));
                    Thread.Sleep(10);
                    events.Enqueue("write end " + i);
                });

            List<string> order = events.ToList();
            Assert.That(order.Where(e => e.StartsWith("prepare")), Is.EqualTo(Enumerable.Range(0, 6).Select(i => "prepare start " + i)));
            Assert.That(order.Where(e => e.StartsWith("write")), Is.EqualTo(Enumerable.Range(0, 6).Select(i => "write end " + i)));
            Assert.That(mostPreparing, Is.EqualTo(1), "preparations never overlap each other");
            for (int i = 2; i < 6; i++)
            {
                Assert.That(order.IndexOf("prepare start " + i), Is.GreaterThan(order.IndexOf("write end " + (i - 2))),
                    "a preparation starts only once the write two items back has finished: one item ahead, no more");
            }
        }

        private static void InterlockedMax(ref int target, int value)
        {
            int seen;
            while (value > (seen = Volatile.Read(ref target)) && Interlocked.CompareExchange(ref target, value, seen) != seen)
            {
            }
        }

        [Test]
        public static void AFailedPreparationIsThrownAsItselfBeforeItsWrite()
        {
            var written = new List<int>();
            Assert.That(() => PostSearchAnalysisTask.PrepareOneAhead(Enumerable.Range(0, 6).ToList(),
                    i => i == 3 ? throw new InvalidOperationException("prepare 3") : i,
                    (i, _) => written.Add(i)),
                Throws.InvalidOperationException.With.Message.EqualTo("prepare 3"));
            Assert.That(written, Is.EqualTo(new[] { 0, 1, 2 }));
        }

        [Test]
        public static void AFailedWriteWaitsForThePreparationBehindIt()
        {
            int finished = 0;
            Assert.That(() => PostSearchAnalysisTask.PrepareOneAhead(Enumerable.Range(0, 6).ToList(),
                    i =>
                    {
                        Thread.Sleep(i == 3 ? 200 : 0);
                        Interlocked.Increment(ref finished);
                        return i;
                    },
                    (i, _) =>
                    {
                        if (i == 2) throw new IOException("write 2");
                    }),
                Throws.TypeOf<IOException>().With.Message.EqualTo("write 2"));
            Assert.That(Volatile.Read(ref finished), Is.EqualTo(4), "the preparation of item 3 had finished when the failure was thrown");
        }

        #endregion

        #region Resolve

        /// <summary>
        /// Every Resolve overload against a verbatim copy of the code before the one-element shortcut: random lists of 0-4
        /// values, with nulls, NaN, infinities, -0.0, near-integers and real modifications, give the same string and value, or
        /// the same exception.
        /// </summary>
        [Test]
        public static void ResolveAnswersAsBefore()
        {
            var random = new Random(20261111);
            string[] strings = { "A", "B", "PEPTIDE", "", " ", null };
            int[] ints = { 0, 1, -1, 7, int.MaxValue };
            double[] doubles = { 0.0, -0.0, 1.0, 1.0004, 0.9996, 2.5, 1e-7, -3.25, 3.0000004, 12345.678901, double.NaN, double.PositiveInfinity, double.NegativeInfinity };
            List<Modification> known = GlobalVariables.AllModsKnown.Where(m => m.ChemicalFormula != null).Take(40).ToList();
            List<Modification> noFormula = GlobalVariables.AllModsKnown.Where(m => m.ChemicalFormula == null).Take(3).ToList();
            List<Modification> modPool = known.Concat(noFormula).Concat(new Modification[] { null }).ToList();

            List<T> Some<T>(T[] pool, int most) => Enumerable.Range(0, random.Next(0, most + 1)).Select(_ => pool[random.Next(pool.Length)]).ToList();

            for (int trial = 0; trial < 3000; trial++)
            {
                List<string> s = Some(strings, 4);
                Same(() => PsmTsvWriter.Resolve(s), () => OldResolve(s), r => $"{r.ResolvedString}|{r.ResolvedValue}");
                Same(() => PsmTsvWriter.Resolve(s, "ambiguous"), () => OldResolve(s, "ambiguous"), r => $"{r.ResolvedString}|{r.ResolvedValue}");
                Same(() => PsmTsvWriter.Resolve(s, null), () => OldResolve(s, null), r => $"{r.ResolvedString}|{r.ResolvedValue}");

                List<int> n = Some(ints, 4);
                Same(() => PsmTsvWriter.Resolve(n), () => OldResolve(n), r => $"{r.ResolvedString}|{r.ResolvedValue}");

                List<double> d = Some(doubles, 4);
                Same(() => PsmTsvWriter.Resolve(d), () => OldResolve(d),
                    r => $"{r.ResolvedString}|{(r.ResolvedValue.HasValue ? BitConverter.DoubleToInt64Bits(r.ResolvedValue.Value).ToString() : "null")}");

                // dictionaries of modifications with formulas (as on a peptide: never null values)
                List<Dictionary<int, Modification>> dicts = Enumerable.Range(0, random.Next(0, 4)).Select(_ => Enumerable.Range(0, random.Next(0, 4))
                    .Select(k => (position: random.Next(1, 12), mod: known[random.Next(known.Count)])).GroupBy(x => x.position)
                    .ToDictionary(g => g.Key, g => g.First().mod)).ToList();
                if (random.Next(3) == 0 && dicts.Count > 1)
                {
                    dicts[1] = dicts[0]; // equal dictionaries are the usual ambiguity
                }
                Same(() => PsmTsvWriter.Resolve(dicts), () => OldResolve(dicts),
                    r => $"{r.ResolvedString}|{(r.ResolvedValue == null ? "null" : string.Join(",", r.ResolvedValue.Select(kv => kv.Key + "=" + kv.Value)))}");

                List<IEnumerable<Modification>> modLists = Enumerable.Range(0, random.Next(0, 4)).Select(_ => (IEnumerable<Modification>)Some(modPool.ToArray(), 3)).ToList();
                if (random.Next(3) == 0 && modLists.Count > 1)
                {
                    modLists[1] = modLists[0];
                }
                Same(() => PsmTsvWriter.Resolve(modLists), () => OldResolve(modLists),
                    r => $"{r.ResolvedString}|{r.ResolvedValue?.Formula ?? "null"}");
            }

            Same(() => PsmTsvWriter.Resolve((IEnumerable<string>)null), () => OldResolve((IEnumerable<string>)null), r => r.ResolvedString);
        }

        private static void Same<T>(Func<T> now, Func<T> before, Func<T, string> show)
        {
            static string Outcome(Func<T> f, Func<T, string> show)
            {
                try
                {
                    return "returned " + show(f());
                }
                catch (Exception e)
                {
                    return "threw " + e.GetType().Name;
                }
            }
            Assert.That(Outcome(now, show), Is.EqualTo(Outcome(before, show)));
        }

        // ---- verbatim copy of PsmTsvWriter.Resolve before the one-element shortcut (master 11313d07e) ----

        private static (string ResolvedString, ChemicalFormula ResolvedValue) OldResolve(IEnumerable<IEnumerable<Modification>> enumerable)
        {
            var list = enumerable.ToList();
            ChemicalFormula firstChemFormula = new ChemicalFormula();
            foreach (var firstMods in list[0])
            {
                if (firstMods == null || firstMods.ChemicalFormula == null)
                {
                    return ("unknown", null);
                }
                firstChemFormula.Add(firstMods.ChemicalFormula);
            }

            bool equals = true;
            List<ChemicalFormula> formulas = new List<ChemicalFormula>();
            foreach (var anEnum in list)
            {
                ChemicalFormula fhere = new ChemicalFormula();
                foreach (var mod in anEnum)
                {
                    if (mod == null || mod.ChemicalFormula == null)
                    {
                        return ("unknown", null);
                    }
                    fhere.Add(mod.ChemicalFormula);
                }
                if (!firstChemFormula.Equals(fhere))
                {
                    equals = false;
                }
                formulas.Add(fhere);
            }
            if (!equals)
            {
                var returnString = GlobalVariables.CheckLengthOfOutput(string.Join("|", formulas.Select(b => b.Formula)));
                return (returnString, null);
            }
            else
            {
                return (firstChemFormula.Formula, firstChemFormula);
            }
        }

        private static (string ResolvedString, Dictionary<string, int> ResolvedValue) OldResolve(IEnumerable<Dictionary<int, Modification>> enumerable)
        {
            var list = enumerable.ToList();
            Dictionary<string, int> firstDict = list[0].Values.OrderBy(b => b.IdWithMotif).GroupBy(b => b.IdWithMotif).ToDictionary(b => b.Key, b => b.Count());

            bool equals = true;
            foreach (var dict in list)
            {
                Dictionary<string, int> okTest = dict.Values.OrderBy(b => b.IdWithMotif).GroupBy(b => b.IdWithMotif).ToDictionary(b => b.Key, b => b.Count());
                if (!firstDict.SequenceEqual(okTest))
                {
                    equals = false;
                    break;
                }
            }
            if (!equals)
            {
                var returnString = string.Join("|", list.Select(b => string.Join(" ", b.Values.Select(c => c.IdWithMotif).OrderBy(c => c))));
                returnString = GlobalVariables.CheckLengthOfOutput(returnString);
                return (returnString, null);
            }
            else
            {
                return (string.Join(" ", list[0].Values.Select(c => c.IdWithMotif).OrderBy(c => c)), firstDict);
            }
        }

        private static (string ResolvedString, double? ResolvedValue) OldResolve(IEnumerable<double> enumerable)
        {
            var list = enumerable.ToList();
            double min = list.Min();
            double max = list.Max();
            double avg = list.Average();

            // If all values are within tolerance, check if all are integers
            if (max - min < PsmTsvWriter.ToleranceForDoubleResolutionF5)
            {
                bool allInts = list.All(x => Math.Abs(x - Math.Round(x)) < PsmTsvWriter.ToleranceForDoubleResolutionF2);
                if (!allInts)
                    return (avg.ToString("F5", CultureInfo.InvariantCulture), avg);

                int intVal = (int)Math.Round(avg);
                return (intVal.ToString(CultureInfo.InvariantCulture), intVal);
            }

            var returnString = GlobalVariables.CheckLengthOfOutput(string.Join("|", list.Select(b => b.ToString("F5", CultureInfo.InvariantCulture))));
            return (returnString, null);
        }

        private static (string ResolvedString, int? ResolvedValue) OldResolve(IEnumerable<int> enumerable)
        {
            var list = enumerable.ToList();
            var first = list[0];
            if (list.All(b => first.Equals(b)))
            {
                return (first.ToString(CultureInfo.InvariantCulture), first);
            }
            else
            {
                var returnString = GlobalVariables.CheckLengthOfOutput(string.Join("|", list.Select(b => b.ToString(CultureInfo.InvariantCulture))));
                return (returnString, null);
            }
        }

        private static (string ResolvedString, string ResolvedValue) OldResolve(IEnumerable<string> enumerable)
        {
            var list = enumerable.ToList();
            string first = list.FirstOrDefault(b => b != null);
            // Only first if list is either all null or all equal to the first
            if (list.All(b => b == null) || list.All(b => first.Equals(b)))
            {
                return (first, first);
            }
            else
            {
                var returnString = GlobalVariables.CheckLengthOfOutput(string.Join("|", list));
                return (returnString, null);
            }
        }

        private static (string ResolvedString, string ResolvedValue) OldResolve(IEnumerable<string> enumerable, string ambiguousIfNull)
        {
            var list = enumerable.ToList();
            string first = list.FirstOrDefault(b => b != null);
            // Only first if list is either all null or all equal to the first
            if (list.All(b => b == null) || list.All(b => first.Equals(b)))
            {
                return (first, first);
            }
            // use only distinct names if all of the base sequences are the same
            else if (ambiguousIfNull != null)
            {
                var returnString = GlobalVariables.CheckLengthOfOutput(string.Join("|", list.Distinct()));
                return (returnString, null);
            }
            else
            {
                var returnString = GlobalVariables.CheckLengthOfOutput(string.Join("|", list));
                return (returnString, null);
            }
        }

        #endregion
    }
}
