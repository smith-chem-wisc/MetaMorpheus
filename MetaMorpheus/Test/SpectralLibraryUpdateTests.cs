using Chemistry;
using EngineLayer;
using EngineLayer.DatabaseLoading;
using MassSpectrometry;
using NUnit.Framework;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;
using Omics.SpectrumMatch;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Readers.SpectralLibrary;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using TaskLayer;

namespace Test
{
    [TestFixture]
    public static class SpectralLibraryUpdateTests
    {
        /// <summary>
        /// Issue #2291. Asking to update a spectral library without giving one used to run the whole search
        /// and then throw NullReferenceException out of UpdateSpectralLibrary, which surfaced as the task
        /// hanging on "Writing PSM results" with the exception only visible in results.txt afterwards.
        ///
        /// The combination is refused before searching now, so the wasted run is what this pins. The check
        /// also precedes the background protein load, so reaching it costs no I/O at all.
        /// </summary>
        [Test]
        public static void UpdatingASpectralLibraryWithoutOneIsRefusedBeforeSearching()
        {
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory,
                "SpectralLibraryUpdateTests", "NoLibrary");
            Directory.CreateDirectory(outputFolder);

            var task = new SearchTask();
            task.SearchParameters.UpdateSpectralLibrary = true;

            string database = Path.Combine(TestContext.CurrentContext.TestDirectory,
                "TestData", "hela_snip_for_unitTest.fasta");
            string spectra = Path.Combine(TestContext.CurrentContext.TestDirectory,
                "TestData", "TaGe_SA_A549_3_snip.mzML");

            // "Before searching" is half the fix, so observe it rather than trusting the method name.
            // Absence of output files does not show it: every result file is written by
            // PostSearchAnalysisTask, so a check placed anywhere ahead of that leaves the folder empty even
            // though the search ran. Engine start events are the signal that discriminates -- the guard
            // precedes the database load, so no engine may start at all.
            var enginesStarted = new List<string>();
            void OnEngineStarting(object sender, SingleEngineEventArgs e)
                => enginesStarted.Add(e.MyEngine.GetType().Name);

            MetaMorpheusEngine.StartingSingleEngineHander += OnEngineStarting;
            MetaMorpheusException thrown;
            try
            {
                thrown = Assert.Throws<MetaMorpheusException>(() => task.RunTask(
                    outputFolder,
                    new List<DbForTask> { new DbForTask(database, false) },
                    new List<string> { spectra },
                    "TestUpdateWithoutLibrary"));
            }
            finally
            {
                MetaMorpheusEngine.StartingSingleEngineHander -= OnEngineStarting;
            }

            Assert.That(thrown.Message, Does.Contain("no spectral library was given").IgnoreCase,
                "say what is missing, not just that something went wrong");
            Assert.That(thrown.Message, Does.Contain("list of databases").IgnoreCase,
                "and say what to do about it, so the message cannot be trimmed to the diagnosis alone");

            Assert.That(enginesStarted, Is.Empty,
                "the refusal has to come before any work starts, but these ran: " + string.Join(", ", enginesStarted));

            // No delete here: RunTask writes its toml one level above this folder, and a failing assertion
            // would skip it anyway. OneTimeTearDown removes the whole root.
        }

        /// <summary>
        /// The guard SearchTask now applies is only reachable when a spectral library is actually in the
        /// database list, so this pins the precondition it reads: LoadSpectralLibraries returns null when
        /// nothing in the list is a library, and a real library when something is. If that ever changed to
        /// return an empty SpectralLibrary instead of null, the refusal above would silently stop firing.
        /// </summary>
        [Test]
        public static void LoadSpectralLibrariesReturnsNullOnlyWhenNoLibraryIsInTheDatabaseList()
        {
            var task = new SearchTask();
            var load = typeof(MetaMorpheusTask).GetMethod("LoadSpectralLibraries",
                BindingFlags.NonPublic | BindingFlags.Instance);
            Assert.That(load, Is.Not.Null, "LoadSpectralLibraries is the seam SearchTask's guard depends on");

            var fastaOnly = new List<DbForTask> { new DbForTask(HelaSnipFasta, false) };
            Assert.That(load.Invoke(task, new object[] { "id", fastaOnly }), Is.Null);

            Assert.That(load.Invoke(task, new object[] { "id", new List<DbForTask>() }), Is.Null,
                "an empty database list has no library either");

            var withLibrary = new List<DbForTask>
            {
                new DbForTask(HelaSnipFasta, false),
                new DbForTask(SmallLibrary, false)
            };
            var loaded = load.Invoke(task, new object[] { "id", withLibrary });
            Assert.That(loaded, Is.Not.Null);
            ((SpectralLibrary)loaded).CloseConnections();
        }

        /// <summary>
        /// The null guard PR #2721 adds to PostSearchAnalysisTask.UpdateSpectralLibrary. SearchTask refuses
        /// the combination before searching, so this branch is only reachable by a caller that assembles
        /// PostSearchAnalysisParameters itself — which is exactly what this test does.
        ///
        /// Without the guard the method dereferences the missing library, the surrounding catch turns that
        /// into an UpdateSpectralLibrary_crash.txt nobody reads, and the task reports success. Asserting the
        /// crash file is absent is what makes this a regression test rather than a restatement of the code.
        /// </summary>
        [Test]
        public static void UpdateSpectralLibraryWithoutALibraryWarnsAndLeavesNothingBehind()
        {
            string outputFolder = MakeOutputFolder("NoLibraryPostSearch");

            var task = new PostSearchAnalysisTask { CommonParameters = new CommonParameters() };
            task.Parameters = new PostSearchAnalysisParameters
            {
                SearchParameters = new SearchParameters { UpdateSpectralLibrary = true },
                OutputFolder = outputFolder,
                SearchTaskId = "NoLibraryPostSearch",
                AllSpectralMatches = new List<SpectralMatch>(),
                DatabaseFilenameList = new List<DbForTask>(),
                SpectralLibrary = null,
                SearchTaskResults = NewTaskResults(task)
            };

            var warnings = CaptureWarnings(() => InvokePrivate(task, "UpdateSpectralLibrary"));

            Assert.That(warnings.Any(w => w.Contains("spectral library", StringComparison.OrdinalIgnoreCase)),
                Is.True, "the user has to be told why no library came out: " + string.Join(" | ", warnings));
            Assert.That(File.Exists(Path.Combine(outputFolder, "UpdateSpectralLibrary_crash.txt")), Is.False,
                "a missing library is a configuration problem, not an engine crash");
            Assert.That(Directory.GetFiles(outputFolder, "updateSpectralLibrary*.msp"), Is.Empty,
                "nothing to update means nothing written");
            Assert.That(task.Parameters.SearchTaskResults.NewDatabases, Is.Null,
                "no library was produced, so none may be advertised to the next task in the run");

            Directory.Delete(outputFolder, recursive: true);
        }

        /// <summary>
        /// The other side of the same guard: given a library and no search results, every original spectrum
        /// has to survive into the updated library, and the updated library plus the original protein
        /// database have to be handed to the next task. This is the path the guard exists to protect, and
        /// until now it was covered only by a full two-file search that takes about seven seconds.
        /// </summary>
        [Test]
        public static void UpdateSpectralLibraryWithNoSearchResultsCarriesEveryOriginalSpectrumThrough()
        {
            string outputFolder = MakeOutputFolder("LibraryNoPsms");
            var library = new SpectralLibrary(new List<string> { SmallLibrary });

            var task = new PostSearchAnalysisTask { CommonParameters = new CommonParameters() };
            task.Parameters = new PostSearchAnalysisParameters
            {
                SearchParameters = new SearchParameters { UpdateSpectralLibrary = true },
                OutputFolder = outputFolder,
                SearchTaskId = "LibraryNoPsms",
                AllSpectralMatches = new List<SpectralMatch>(),
                DatabaseFilenameList = new List<DbForTask>
                {
                    new DbForTask(SmallLibrary, false),
                    new DbForTask(HelaSnipFasta, false)
                },
                SpectralLibrary = library,
                SearchTaskResults = NewTaskResults(task)
            };

            var originalSequences = library.GetAllLibrarySpectra().Select(s => s.Sequence).OrderBy(s => s).ToList();
            Assert.That(originalSequences.Count, Is.EqualTo(3), "fixture changed; the assertions below assume three spectra");

            var warnings = CaptureWarnings(() => InvokePrivate(task, "UpdateSpectralLibrary"));
            library.CloseConnections();

            Assert.That(warnings, Is.Empty, "nothing went wrong: " + string.Join(" | ", warnings));

            var written = Directory.GetFiles(outputFolder, "updateSpectralLibrary*.msp");
            Assert.That(written.Length, Is.EqualTo(1), "exactly one updated library should be written");

            var updated = new SpectralLibrary(new List<string> { written[0] });
            var updatedSequences = updated.GetAllLibrarySpectra().Select(s => s.Sequence).OrderBy(s => s).ToList();
            updated.CloseConnections();
            Assert.That(updatedSequences, Is.EqualTo(originalSequences),
                "with no PSMs to compare against, an update must not lose or invent a spectrum");

            var newDatabases = task.Parameters.SearchTaskResults.NewDatabases;
            Assert.That(newDatabases.Count, Is.EqualTo(2));
            Assert.That(newDatabases[0].FilePath, Is.EqualTo(written[0]));
            Assert.That(newDatabases[0].IsSpectralLibrary, Is.True);
            Assert.That(newDatabases[1].FilePath, Is.EqualTo(HelaSnipFasta),
                "the protein database has to travel with the library or the next task has nothing to search");

            Directory.Delete(outputFolder, recursive: true);
        }

        /// <summary>
        /// The selection rule at the heart of the update: for a peptide/charge already in the library, the
        /// search result replaces the stored spectrum only when the stored one has fewer matched ions than
        /// the PSM's score. Both outcomes are exercised here from the same fixture, so a change to the
        /// comparison shows up as a failure rather than as a quietly different library.
        /// </summary>
        [Test]
        [TestCase(30d, true, TestName = "UpdateSpectralLibraryReplacesTheStoredSpectrumWhenThePsmScoresHigher")]
        [TestCase(2d, false, TestName = "UpdateSpectralLibraryKeepsTheStoredSpectrumWhenThePsmScoresLower")]
        public static void UpdateSpectralLibraryPrefersTheBetterOfTheStoredSpectrumAndTheSearchResult(
            double psmScore, bool expectReplacement)
        {
            string outputFolder = MakeOutputFolder("LibraryVsPsm_" + psmScore);
            var library = new SpectralLibrary(new List<string> { SmallLibrary });

            var stored = library.GetAllLibrarySpectra().Single(s => s.Sequence == RepeatedSequence);
            Assert.That(stored.ChargeState, Is.EqualTo(2), "fixture changed; the PSM below is built at charge 2");

            var psm = MakePsm(RepeatedSequence, charge: 2, score: psmScore);
            Assert.That(psm.MatchedFragmentIons.Count, Is.Not.EqualTo(stored.MatchedFragmentIons.Count),
                "the two spectra must be distinguishable for this test to mean anything");

            var task = new PostSearchAnalysisTask { CommonParameters = new CommonParameters() };
            task.Parameters = new PostSearchAnalysisParameters
            {
                SearchParameters = new SearchParameters { UpdateSpectralLibrary = true },
                OutputFolder = outputFolder,
                SearchTaskId = "LibraryVsPsm",
                AllSpectralMatches = new List<SpectralMatch> { psm },
                DatabaseFilenameList = new List<DbForTask>
                {
                    new DbForTask(SmallLibrary, false),
                    new DbForTask(HelaSnipFasta, false)
                },
                SpectralLibrary = library,
                SearchTaskResults = NewTaskResults(task)
            };

            InvokePrivate(task, "UpdateSpectralLibrary");
            library.CloseConnections();

            var written = Directory.GetFiles(outputFolder, "updateSpectralLibrary*.msp").Single();
            var updated = new SpectralLibrary(new List<string> { written });
            var all = updated.GetAllLibrarySpectra().ToList();
            var result = all.Single(s => s.Sequence == RepeatedSequence);
            int resultIonCount = result.MatchedFragmentIons.Count;
            updated.CloseConnections();

            Assert.That(all.Count, Is.EqualTo(3),
                "the peptide was already in the library, so it must be replaced rather than appended");
            Assert.That(resultIonCount, Is.EqualTo(expectReplacement
                    ? psm.MatchedFragmentIons.Count
                    : stored.MatchedFragmentIons.Count),
                expectReplacement
                    ? "a PSM scoring above the stored ion count should win"
                    : "a PSM scoring below the stored ion count should not displace it");

            Directory.Delete(outputFolder, recursive: true);
        }

        /// <summary>
        /// A confidently identified peptide that the library has never seen has to be added, not dropped —
        /// growing the library is the point of asking for an update at all.
        /// </summary>
        [Test]
        public static void UpdateSpectralLibraryAppendsAPeptideTheLibraryHasNotSeen()
        {
            string outputFolder = MakeOutputFolder("LibraryNewPeptide");
            var library = new SpectralLibrary(new List<string> { SmallLibrary });

            const string newSequence = "PEPTIDEK";
            Assert.That(library.GetAllLibrarySpectra().Any(s => s.Sequence == newSequence), Is.False,
                "fixture changed; this peptide is supposed to be absent from the library");

            var task = new PostSearchAnalysisTask { CommonParameters = new CommonParameters() };
            task.Parameters = new PostSearchAnalysisParameters
            {
                SearchParameters = new SearchParameters { UpdateSpectralLibrary = true },
                OutputFolder = outputFolder,
                SearchTaskId = "LibraryNewPeptide",
                AllSpectralMatches = new List<SpectralMatch> { MakePsm(newSequence, charge: 2, score: 12) },
                DatabaseFilenameList = new List<DbForTask>
                {
                    new DbForTask(SmallLibrary, false),
                    new DbForTask(HelaSnipFasta, false)
                },
                SpectralLibrary = library,
                SearchTaskResults = NewTaskResults(task)
            };

            InvokePrivate(task, "UpdateSpectralLibrary");
            library.CloseConnections();

            var written = Directory.GetFiles(outputFolder, "updateSpectralLibrary*.msp").Single();
            var updated = new SpectralLibrary(new List<string> { written });
            var sequences = updated.GetAllLibrarySpectra().Select(s => s.Sequence).ToList();
            updated.CloseConnections();

            Assert.That(sequences.Count, Is.EqualTo(4), "three originals plus the newly identified peptide");
            Assert.That(sequences, Does.Contain(newSequence));

            Directory.Delete(outputFolder, recursive: true);
        }

        /// <summary>
        /// A decoy or a low-confidence hit must never reach the library — an updated library is used as
        /// ground truth by the next search, so anything that leaks in is amplified.
        /// </summary>
        [Test]
        public static void UpdateSpectralLibraryIgnoresDecoyAndLowConfidenceMatches()
        {
            string outputFolder = MakeOutputFolder("LibraryFiltering");
            var library = new SpectralLibrary(new List<string> { SmallLibrary });

            var decoy = MakePsm("DECOYPEPTIDEK", charge: 2, score: 40, isDecoy: true);
            var lowConfidence = MakePsm("BADQVALUEPEPTIDEK", charge: 2, score: 40, qValue: 0.5);

            var task = new PostSearchAnalysisTask { CommonParameters = new CommonParameters() };
            task.Parameters = new PostSearchAnalysisParameters
            {
                SearchParameters = new SearchParameters { UpdateSpectralLibrary = true },
                OutputFolder = outputFolder,
                SearchTaskId = "LibraryFiltering",
                AllSpectralMatches = new List<SpectralMatch> { decoy, lowConfidence },
                DatabaseFilenameList = new List<DbForTask>
                {
                    new DbForTask(SmallLibrary, false),
                    new DbForTask(HelaSnipFasta, false)
                },
                SpectralLibrary = library,
                SearchTaskResults = NewTaskResults(task)
            };

            InvokePrivate(task, "UpdateSpectralLibrary");
            library.CloseConnections();

            var written = Directory.GetFiles(outputFolder, "updateSpectralLibrary*.msp").Single();
            var updated = new SpectralLibrary(new List<string> { written });
            var sequences = updated.GetAllLibrarySpectra().Select(s => s.Sequence).ToList();
            updated.CloseConnections();

            Assert.That(sequences.Count, Is.EqualTo(3), "neither match qualifies, so the library is unchanged");
            Assert.That(sequences, Does.Not.Contain("DECOYPEPTIDEK"));
            Assert.That(sequences, Does.Not.Contain("BADQVALUEPEPTIDEK"));

            Directory.Delete(outputFolder, recursive: true);
        }

        #region helpers

        /// <summary>KAPAGGAADAAAK is one of the three spectra in spectralLibraryForTestingLibraryUpdate.msp.</summary>
        private const string RepeatedSequence = "KAPAGGAADAAAK";

        private static string SmallLibrary => Path.Combine(TestContext.CurrentContext.TestDirectory,
            "TestData", "SpectralLibrarySearch", "spectralLibraryForTestingLibraryUpdate.msp");

        private static string HelaSnipFasta => Path.Combine(TestContext.CurrentContext.TestDirectory,
            "TestData", "hela_snip_for_unitTest.fasta");

        private static string MakeOutputFolder(string name)
        {
            string folder = Path.Combine(TestContext.CurrentContext.TestDirectory,
                "SpectralLibraryUpdateTests", name);
            if (Directory.Exists(folder))
            {
                Directory.Delete(folder, true);
            }
            Directory.CreateDirectory(folder);
            return folder;
        }

        /// <summary>
        /// A scored, FDR-resolved PSM for the given sequence, built from a synthetic scan of that peptide's
        /// own fragments so the matched ion list is real rather than empty. Score is set independently of
        /// the ion count because the update rule compares the two against each other.
        /// </summary>
        private static SpectralMatch MakePsm(string sequence, int charge, double score,
            bool isDecoy = false, double qValue = 0)
        {
            var commonParameters = new CommonParameters();
            var protein = new Protein(sequence, isDecoy ? "DECOY_acc" : "acc", isDecoy: isDecoy);
            var peptide = new PeptideWithSetModifications(protein, commonParameters.DigestionParams,
                1, sequence.Length, CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0);

            MsDataFile dataFile = new TestDataFile(peptide, "quadratic");
            var scan = new Ms2ScanWithSpecificMass(dataFile.GetOneBasedScan(2),
                peptide.MonoisotopicMass.ToMz(charge), charge, null, commonParameters);

            var theoretical = new List<Product>();
            peptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, theoretical);
            var matched = MetaMorpheusEngine.MatchFragmentIons(scan, theoretical, commonParameters);

            SpectralMatch psm = new PeptideSpectralMatch(peptide, 0, score, 1, scan, commonParameters, matched);
            psm.ResolveAllAmbiguities();
            psm.SetFdrValues(0, 0, qValue, 0, 0, 0, 0, 0);
            return psm;
        }

        /// <summary>MyTaskResults has an internal constructor and TaskLayer does not expose internals to Test.</summary>
        private static MyTaskResults NewTaskResults(MetaMorpheusTask task) =>
            (MyTaskResults)Activator.CreateInstance(typeof(MyTaskResults),
                BindingFlags.Instance | BindingFlags.NonPublic, null, new object[] { task }, null);

        private static void InvokePrivate(PostSearchAnalysisTask task, string methodName)
        {
            // The parameterless overload has to be named explicitly: MetaMorpheusTask carries a
            // protected UpdateSpectralLibrary(List<LibrarySpectrum>, string) that would otherwise tie.
            var method = typeof(PostSearchAnalysisTask).GetMethod(methodName,
                BindingFlags.NonPublic | BindingFlags.Instance, null, Type.EmptyTypes, null);
            Assert.That(method, Is.Not.Null, methodName + " was renamed or removed");
            try
            {
                method.Invoke(task, null);
            }
            catch (TargetInvocationException e)
            {
                throw e.InnerException ?? e;
            }
        }

        private static List<string> CaptureWarnings(Action action)
        {
            var warnings = new List<string>();
            void Handler(object sender, StringEventArgs e) => warnings.Add(e.S);

            MetaMorpheusTask.WarnHandler += Handler;
            try
            {
                action();
            }
            finally
            {
                MetaMorpheusTask.WarnHandler -= Handler;
            }
            return warnings;
        }

        /// <summary>
        /// RunTask writes its settings toml to a "Task Settings" folder beside the output folder, so the
        /// fixture's own directory outlives the tests that created it.
        /// </summary>
        [OneTimeTearDown]
        public static void TearDown()
        {
            string root = Path.Combine(TestContext.CurrentContext.TestDirectory, "SpectralLibraryUpdateTests");
            if (Directory.Exists(root))
            {
                Directory.Delete(root, true);
            }
        }

        #endregion
    }
}
