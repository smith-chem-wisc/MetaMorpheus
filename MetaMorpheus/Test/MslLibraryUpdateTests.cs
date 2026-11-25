using Chemistry;
using EngineLayer;
using EngineLayer.DatabaseLoading;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.SpectralMatch.MslSpectralLibrary;
using Omics.SpectrumMatch;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Readers.SpectralLibrary;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Net.Http;
using System.Text;
using TaskLayer;

namespace Test
{
    /// <summary>
    /// The spectral-library update when the library is an .msl (PostSearchAnalysisTaskMsl.cs): spectra, then times
    /// calibrated onto the iRT scale per raw file, then predicted times for entries never measured.
    ///
    /// Nothing here touches the network: the Koina predictor is replaced through PostSearchAnalysisTask.IrtPredictor
    /// with a table of iRT values, and each raw file's scan times are a straight line of those values. The live
    /// Koina test is in MslLibraryUpdateLiveTests. NonParallelizable: some tests set the process-wide StopLoops and
    /// AnalyteType, and every test listens to the static WarnHandler.
    /// </summary>
    [TestFixture]
    [NonParallelizable]
    [ExcludeFromCodeCoverage]
    public class MslLibraryUpdateTests
    {
        private const int Charge = 2;
        private const string FileA = @"C:\data\runA.raw";
        private const string FileB = @"C:\data\runB.raw";
        private const string NeverMeasured = "NEVERMEASVREDK";
        private const string NoPrediction = "NQPREDICTIQNK";
        private const string KeptElsewhere = "KEPTELSEWHEREK";
        private const string Decoy = "YQEDITPEPK";
        /// <summary>Identified in both files in the void at the start of the run, so no time is taken for it.</summary>
        private const string InTheVoid = "VQIDFLEMAGHK";
        /// <summary>Never measured, unmodified and 35 residues long: past Prosit, so Chronologer predicts it.</summary>
        private const string LongNeverMeasured = "ADEFGHILMNQSTVWYADEFGHILMNQSTVWYADE";

        private static readonly CommonParameters SearchParameters = new(
            digestionParams: new DigestionParams(protease: "trypsin"),
            productMassTolerance: new PpmTolerance(20),
            precursorMassTolerance: new PpmTolerance(5));

        private string _folder;
        private readonly List<string> _warnings = new();
        private void CaptureWarning(object sender, StringEventArgs e) => _warnings.Add(e.S);

        /// <summary>The peptides the search identifies, with the iRT the substitute predictor gives each.</summary>
        private static readonly List<(string Sequence, double Irt)> Peptides = SearchPeptides(150);

        [SetUp]
        public void SetUp()
        {
            _folder = Path.Combine(TestContext.CurrentContext.TestDirectory, "MslLibraryUpdate");
            if (Directory.Exists(_folder))
                Directory.Delete(_folder, true);
            Directory.CreateDirectory(_folder);
            _warnings.Clear();
            MetaMorpheusTask.WarnHandler += CaptureWarning;
        }

        [TearDown]
        public void TearDown()
        {
            MetaMorpheusTask.WarnHandler -= CaptureWarning;
            GlobalVariables.StopLoops = false;
            GlobalVariables.AnalyteType = AnalyteType.Peptide;
            Directory.Delete(_folder, true);
        }

        [Test]
        public void AnMslLibrary_IsUpdatedAsMsl_WithEachFilesTimesOnTheIrtScale()
        {
            var asked = new List<string>();
            var task = TaskUpdating(WriteLibrary(), sequences =>
            {
                asked.AddRange(sequences);
                return sequences.Select(PredictedIrt).ToArray();
            });

            task.UpdateSpectralLibrary();

            string written = Directory.GetFiles(_folder, "updateSpectralLibrary_*").Single();
            Assert.That(Path.GetExtension(written), Is.EqualTo(".msl"));
            Assert.That(task.Parameters.SearchTaskResults.NewDatabases.Select(d => d.FilePath),
                Is.EqualTo(new[] { written, "proteins.fasta" }));
            Assert.That(File.Exists(Path.Combine(_folder, "UpdateSpectralLibrary_crash.txt")), Is.False);
            Assert.That(_warnings, Is.Empty);
            AssertTheDiaLibrarySearchOpens(written);

            using var library = MslLibrary.Load(written);
            var entries = library.GetAllEntries().ToDictionary(e => (e.FullSequence, e.IsDecoy));

            // Each identified peptide: the search's spectrum where it has more ions, the file's dissociation type
            // where the spectrum is the search's, and a time from both files, mapped onto the iRT scale.
            foreach (var (sequence, irt) in Peptides)
            {
                var entry = entries[(sequence, false)];
                Assert.That(entry.RtStatus, Is.Not.EqualTo(MslFormat.RetentionTimeStatus.None), sequence);
                Assert.That(entry.RtIsCalibrated, Is.False, sequence);
                Assert.That(entry.RetentionTime, Is.EqualTo(irt).Within(2), sequence);
            }
            var replaced = entries[(Peptides[0].Sequence, false)];
            Assert.That(replaced.MatchedFragmentIons, Has.Count.EqualTo(5));
            Assert.That(replaced.DissociationType, Is.EqualTo(DissociationType.CID)); // file A's, file-specific
            var keptForItsIons = entries[(Peptides[1].Sequence, false)];
            Assert.That(keptForItsIons.MatchedFragmentIons, Has.Count.EqualTo(6));
            Assert.That(keptForItsIons.DissociationType, Is.EqualTo(DissociationType.ETD)); // the library's own
            Assert.That(entries[(Peptides[2].Sequence, false)].DissociationType, Is.EqualTo(DissociationType.CID));

            Assert.That(replaced.RtStatus, Is.EqualTo(MslFormat.RetentionTimeStatus.Confirmed), "a better spectrum keeps the confirmed time");
            Assert.That(entries[(InTheVoid, false)].RetentionTime, Is.NaN);
            Assert.That(entries[(InTheVoid, false)].RtStatus, Is.EqualTo(MslFormat.RetentionTimeStatus.None));

            // Entries the search did not see: a measured time is kept, a never-measured one gets a prediction
            Assert.That(entries[(KeptElsewhere, false)].RetentionTime, Is.EqualTo(55.5).Within(1e-4));
            Assert.That(entries[(KeptElsewhere, false)].RtStatus, Is.EqualTo(MslFormat.RetentionTimeStatus.Confirmed));
            Assert.That(entries[(KeptElsewhere, false)].DissociationType, Is.EqualTo(DissociationType.ETD));
            Assert.That(entries[(NeverMeasured, false)].RetentionTime, Is.EqualTo(42.0).Within(1e-4));
            Assert.That(entries[(NeverMeasured, false)].RtStatus, Is.EqualTo(MslFormat.RetentionTimeStatus.None));
            Assert.That(entries[(NoPrediction, false)].RetentionTime, Is.NaN);
            Assert.That(entries[(LongNeverMeasured, false)].RetentionTime, Is.Not.NaN, () => "Chronologer, mapped onto the library's iRT: " + task.Parameters.SearchTaskResults);
            Assert.That(asked, Does.Not.Contain(LongNeverMeasured), "Prosit is not asked past 30 residues");
            Assert.That(entries.ContainsKey((Decoy, true)), "the library's decoy is kept as a decoy");

            // Each sequence is asked once for the whole update: the peptide in the void is asked as a calibration
            // candidate, refused a time, and not asked again for a predicted one
            Assert.That(asked, Is.Unique);
            Assert.That(asked, Does.Contain(InTheVoid));
            Assert.That(asked, Does.Contain(NeverMeasured));

            string summary = task.Parameters.SearchTaskResults.ToString();
            Assert.That(summary, Does.Contain("Spectral library updated as .msl"));
            Assert.That(summary, Does.Contain($"{Peptides.Count - 1} added, 1 replaced, 6 kept"));
            Assert.That(summary, Does.Contain("2 of 2 files calibrated"));
            Assert.That(summary, Does.Contain("Predicted times given to 2 entries never measured"));
        }

        [Test]
        public void AKoinaOutage_KeepsTheSpectra_LeavesTheTimesForAnotherUpdate_AndSaysWhy()
        {
            int calls = 0;
            var task = TaskUpdating(WriteLibrary(), _ =>
            {
                calls++;
                throw new HttpRequestException("Koina is unreachable");
            });

            task.UpdateSpectralLibrary();

            Assert.That(calls, Is.EqualTo(1), "a predictor that failed is not asked again");
            string written = Directory.GetFiles(_folder, "updateSpectralLibrary_*").Single();
            Assert.That(task.Parameters.SearchTaskResults.NewDatabases[0].FilePath, Is.EqualTo(written));
            Assert.That(File.Exists(Path.Combine(_folder, "UpdateSpectralLibrary_crash.txt")), Is.False);
            AssertTheDiaLibrarySearchOpens(written);

            using var library = MslLibrary.Load(written);
            var entries = library.GetAllEntries().ToDictionary(e => (e.FullSequence, e.IsDecoy));
            foreach (var (sequence, _) in Peptides.Skip(2))
            {
                Assert.That(entries[(sequence, false)].RetentionTime, Is.NaN, sequence);
                Assert.That(entries[(sequence, false)].RtStatus, Is.EqualTo(MslFormat.RetentionTimeStatus.None), sequence);
            }
            Assert.That(entries[(Peptides[0].Sequence, false)].MatchedFragmentIons, Has.Count.EqualTo(5));
            Assert.That(entries[(KeptElsewhere, false)].RetentionTime, Is.EqualTo(55.5).Within(1e-4));
            Assert.That(entries[(NeverMeasured, false)].RetentionTime, Is.NaN);

            Assert.That(_warnings.Single(), Does.Contain("retention times were not updated: HttpRequestException: Koina is unreachable"));
            string summary = task.Parameters.SearchTaskResults.ToString();
            Assert.That(summary, Does.Contain("Retention times not updated: HttpRequestException: Koina is unreachable"));
            Assert.That(summary, Does.Contain("No predicted times given: the predictor was unavailable"));
        }

        [Test]
        public void APredictorFailing_AfterTheTimes_KeepsTheTimes_AndSaysWhy()
        {
            var task = TaskUpdating(WriteLibrary(), sequences => sequences.Contains(NeverMeasured)
                ? throw new HttpRequestException("Koina went down")
                : sequences.Select(PredictedIrt).ToArray());

            task.UpdateSpectralLibrary();

            using var library = MslLibrary.Load(Directory.GetFiles(_folder, "updateSpectralLibrary_*").Single());
            var entries = library.GetAllEntries().ToDictionary(e => (e.FullSequence, e.IsDecoy));
            Assert.That(entries[(Peptides[5].Sequence, false)].RetentionTime, Is.EqualTo(Peptides[5].Irt).Within(2));
            Assert.That(entries[(NeverMeasured, false)].RetentionTime, Is.NaN);
            Assert.That(_warnings.Single(), Does.Contain("were not given a predicted one: HttpRequestException: Koina went down"));
            Assert.That(task.Parameters.SearchTaskResults.ToString(), Does.Contain("2 of 2 files calibrated"));
        }

        [Test]
        public void Stop_DuringThePredictions_WritesNoLibrary()
        {
            var task = TaskUpdating(WriteLibrary(), sequences =>
            {
                GlobalVariables.StopLoops = true;
                return sequences.Select(PredictedIrt).ToArray();
            });

            task.UpdateSpectralLibrary();

            Assert.That(Directory.GetFiles(_folder, "updateSpectralLibrary_*"), Is.Empty);
            Assert.That(task.Parameters.SearchTaskResults.NewDatabases, Is.Null);
        }

        [Test]
        public void Stop_BeforeTheTimes_WritesNoLibrary()
        {
            // No substitute predictor, so this is the Koina model; Stop is checked before each call, so it is never asked
            var task = TaskUpdating(WriteLibrary(), predictor: null);
            GlobalVariables.StopLoops = true;

            task.UpdateSpectralLibrary();

            Assert.That(Directory.GetFiles(_folder, "updateSpectralLibrary_*"), Is.Empty);
            Assert.That(task.Parameters.SearchTaskResults.NewDatabases, Is.Null);
        }

        [Test]
        public void Stop_DuringThePredictedTimes_WritesNoLibrary()
        {
            // Only step 3 asks for a peptide the search never measured
            var task = TaskUpdating(WriteLibrary(), sequences =>
            {
                if (sequences.Contains(NeverMeasured))
                    GlobalVariables.StopLoops = true;
                return sequences.Select(PredictedIrt).ToArray();
            });

            task.UpdateSpectralLibrary();

            Assert.That(Directory.GetFiles(_folder, "updateSpectralLibrary_*"), Is.Empty);
            Assert.That(task.Parameters.SearchTaskResults.NewDatabases, Is.Null);
        }

        [Test]
        public void ANonPeptideSearch_UpdatesTheSpectra_WithoutAskingProsit()
        {
            var task = TaskUpdating(WriteLibrary(), _ => throw new AssertionException("Prosit must not be asked"));
            GlobalVariables.AnalyteType = AnalyteType.Oligo;

            task.UpdateSpectralLibrary();

            string written = Directory.GetFiles(_folder, "updateSpectralLibrary_*").Single();
            using var library = MslLibrary.Load(written);
            Assert.That(library.GetAllEntries().Count(), Is.EqualTo(Peptides.Count + 6));
            Assert.That(task.Parameters.SearchTaskResults.ToString(),
                Does.Contain("Retention times not updated: Prosit predicts peptides, and this is a Oligo search"));
        }

        [Test]
        public void OnlyALoneMslLibrary_IsUpdatedAsMsl()
        {
            string msl = WriteLibrary();
            string msp = Path.Combine(_folder, "other.msp");

            Assert.That(TaskWithLibraries(msl).SingleMslLibraryPath(), Is.EqualTo(msl));
            Assert.That(TaskWithLibraries(msp).SingleMslLibraryPath(), Is.Null);
            Assert.That(_warnings, Is.Empty);

            Assert.That(TaskWithLibraries(msl, msp).SingleMslLibraryPath(), Is.Null);
            Assert.That(_warnings.Single(), Does.Contain("written together as .msp"));
        }

        [Test]
        public void NoFileCalibrating_WritesNoMinutes_AndSaysWhichFilesAndWhy()
        {
            // The predictor knows none of the search's peptides: no anchors, so neither file can be put on the iRT scale
            var task = TaskUpdating(WriteLibrary(), sequences => sequences.Select(_ => double.NaN).ToArray());

            task.UpdateSpectralLibrary();

            string written = Directory.GetFiles(_folder, "updateSpectralLibrary_*").Single();
            AssertTheDiaLibrarySearchOpens(written);
            using var library = MslLibrary.Load(written);
            var entries = library.GetAllEntries().ToDictionary(e => (e.FullSequence, e.IsDecoy));
            foreach (var (sequence, _) in Peptides.Skip(2))
                Assert.That(entries[(sequence, false)].RetentionTime, Is.NaN, sequence);
            Assert.That(entries[(Peptides[0].Sequence, false)].RetentionTime, Is.EqualTo(Peptides[0].Irt).Within(1e-4));
            Assert.That(entries[(LongNeverMeasured, false)].RetentionTime, Is.NaN, "two known times are too few to map Chronologer");

            string summary = task.Parameters.SearchTaskResults.ToString();
            Assert.That(summary, Does.Contain("0 of 2 files calibrated (not: runA, too few anchors, runB, too few anchors)"));
            Assert.That(summary, Does.Contain("Chronologer (peptides of 31-50 residues) not used: ArgumentException"));
        }

        [Test]
        public void APredictorGivingTheWrongNumberOfTimes_UpdatesNoTimes_AndSaysWhy()
        {
            var task = TaskUpdating(WriteLibrary(), _ => Array.Empty<double>());

            task.UpdateSpectralLibrary();

            using var library = MslLibrary.Load(Directory.GetFiles(_folder, "updateSpectralLibrary_*").Single());
            var entries = library.GetAllEntries().ToDictionary(e => (e.FullSequence, e.IsDecoy));
            Assert.That(entries[(Peptides[0].Sequence, false)].MatchedFragmentIons, Has.Count.EqualTo(5), "the spectra still update");
            Assert.That(entries[(Peptides[5].Sequence, false)].RetentionTime, Is.NaN);
            Assert.That(_warnings[0], Does.Contain("retention times were not updated: InvalidOperationException: The iRT predictor returned 0 values for"));
        }

        [Test]
        public void AFileSearchedWithAutodetect_GivesItsSpectraNoDissociationType()
        {
            var task = TaskUpdating(WriteLibrary(), sequences => sequences.Select(PredictedIrt).ToArray());
            task.FileSpecificParameters = new List<(string, CommonParameters)>
            {
                (FileA, new CommonParameters(dissociationType: DissociationType.Autodetect))
            };

            task.UpdateSpectralLibrary();

            using var library = MslLibrary.Load(Directory.GetFiles(_folder, "updateSpectralLibrary_*").Single());
            var replaced = library.GetAllEntries().Single(e => e.FullSequence == Peptides[0].Sequence && !e.IsDecoy);
            Assert.That(replaced.DissociationType, Is.EqualTo(DissociationType.Unknown), "Autodetect decides per scan, and a PSM does not record it");
        }

        /// <summary>
        /// The DIA library search refuses a library holding run minutes (DiaLibrarySearchTask.OpenLibrary); what the update
        /// writes must open there, times or not.
        /// </summary>
        private static void AssertTheDiaLibrarySearchOpens(string path)
        {
            MslLibrary opened = null;
            Assert.DoesNotThrow(() => opened = DiaLibrarySearchTask.OpenLibrary(path), "the DIA library search refused it");
            opened?.Dispose();
        }

        /// <summary>
        /// The library before the search: the first identified peptide with fewer ions than the search finds (replaced),
        /// with a confirmed time, and the second with more (kept), both from an ETD library; a peptide the search does not see, with a
        /// confirmed time; three never measured, one with no prediction and one past Prosit's 30 residues; and a decoy.
        /// </summary>
        private string WriteLibrary()
        {
            var entries = new List<MslLibraryEntry>
            {
                Entry(Peptides[0].Sequence, ions: 3, time: Peptides[0].Irt, MslFormat.RetentionTimeStatus.Confirmed),
                Entry(Peptides[1].Sequence, ions: 6),
                Entry(KeptElsewhere, ions: 4, time: 55.5, MslFormat.RetentionTimeStatus.Confirmed),
                Entry(NeverMeasured, ions: 4),
                Entry(NoPrediction, ions: 4),
                Entry(LongNeverMeasured, ions: 4),
                Entry(Decoy, ions: 4, decoy: true),
            };
            string path = Path.Combine(_folder, "library.msl");
            MslLibrary.Save(path, entries);
            return path;
        }

        internal static MslLibraryEntry Entry(string sequence, int ions, double time = double.NaN,
            MslFormat.RetentionTimeStatus status = MslFormat.RetentionTimeStatus.None, bool decoy = false)
        {
            var peptide = new PeptideWithSetModifications(sequence, new Dictionary<string, Omics.Modifications.Modification>());
            var entry = MslLibraryEntry.FromLibrarySpectrum(new LibrarySpectrum(sequence, peptide.MonoisotopicMass.ToMz(Charge), Charge,
                Fragments(peptide, ions), double.IsNaN(time) ? null : time, decoy));
            entry.RtStatus = status;
            entry.RtIsCalibrated = false;
            entry.DissociationType = DissociationType.ETD;
            return entry;
        }

        /// <summary>
        /// A task updating <paramref name="libraryPath"/> from two raw files that identify every peptide in
        /// <see cref="Peptides"/> and <see cref="InTheVoid"/>, file A with the higher score (so its spectrum is the library's) and searched as CID.
        /// </summary>
        private PostSearchAnalysisTask TaskUpdating(string libraryPath, Func<IReadOnlyList<string>, IReadOnlyList<double>> predictor)
        {
            var psms = new List<SpectralMatch>();
            for (int i = 0; i < Peptides.Count; i++)
            {
                double irt = Peptides[i].Irt;
                double wobble = 0.02 * ((i * 7919 % 11) - 5);
                psms.Add(Psm(Peptides[i].Sequence, FileA, 10 + 0.5 * (irt + 10) + wobble, score: 60));
                psms.Add(Psm(Peptides[i].Sequence, FileB, 12 + 0.55 * (irt + 10) - wobble, score: 50));
            }
            psms.Add(Psm(InTheVoid, FileA, 1.0, score: 60));
            psms.Add(Psm(InTheVoid, FileB, 1.0, score: 50));
            var task = TaskWithLibraries(libraryPath);
            task.IrtPredictor = predictor;
            task.Parameters.AllSpectralMatches = psms;
            task.FileSpecificParameters = new List<(string, CommonParameters)>
            {
                (FileA, new CommonParameters(dissociationType: DissociationType.CID))
            };
            return task;
        }

        private PostSearchAnalysisTask TaskWithLibraries(params string[] libraries) => LibraryUpdateTask(_folder, libraries);

        internal static PostSearchAnalysisTask LibraryUpdateTask(string folder, params string[] libraries)
        {
            var task = new PostSearchAnalysisTask { CommonParameters = SearchParameters };
            task.Parameters = new PostSearchAnalysisParameters
            {
                SearchParameters = new SearchParameters { UpdateSpectralLibrary = true },
                AllSpectralMatches = new List<SpectralMatch>(),
                DatabaseFilenameList = new[] { "proteins.fasta" }.Concat(libraries).Select(p => new DbForTask(p, false)).ToList(),
                OutputFolder = folder,
                SearchTaskId = "test",
                SearchTaskResults = new MyTaskResults(task)
            };
            return task;
        }

        /// <summary>A confident, resolved PSM of <paramref name="sequence"/> at <paramref name="minutes"/> in <paramref name="file"/>, with five fragments.</summary>
        private static SpectralMatch Psm(string sequence, string file, double minutes, double score)
        {
            var peptide = new Protein(sequence, "accession")
                .Digest(SearchParameters.DigestionParams, new List<Omics.Modifications.Modification>(), new List<Omics.Modifications.Modification>())
                .First(p => p.BaseSequence == sequence);
            var matched = Fragments(peptide, 5);
            var scan = new MsDataScan(
                massSpectrum: new MzSpectrum(matched.Select(m => m.Mz).ToArray(), matched.Select(_ => 100.0).ToArray(), false),
                oneBasedScanNumber: 1, msnOrder: 2, isCentroid: true, polarity: Polarity.Positive,
                retentionTime: minutes, scanWindowRange: new MzRange(50, 2000), scanFilter: "FTMS + p NSI",
                mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: 1000, injectionTime: null,
                noiseData: null, nativeId: "scan=1");
            var psm = new PeptideSpectralMatch(peptide, 0, score, 0,
                new Ms2ScanWithSpecificMass(scan, peptide.MonoisotopicMass.ToMz(Charge), Charge, file, SearchParameters),
                SearchParameters, matched);
            psm.ResolveAllAmbiguities();
            psm.SetFdrValues(1, 0, 0, 1, 0, 0, 0, 0);
            return psm;
        }

        private static List<MatchedFragmentIon> Fragments(Omics.IBioPolymerWithSetMods peptide, int count)
        {
            var products = new List<Product>();
            peptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
            return products.Take(count).Select(p => new MatchedFragmentIon(p, p.ToMz(1), 100, 1)).ToList();
        }

        private static double PredictedIrt(string sequence)
        {
            if (sequence == NeverMeasured)
                return 42.0;
            int i = Peptides.FindIndex(p => p.Sequence == sequence);
            return i < 0 ? double.NaN : Peptides[i].Irt;
        }

        /// <summary>
        /// <paramref name="count"/> distinct tryptic peptides of 8-18 residues (no K, R or P inside, no C), with iRT
        /// spread evenly from -10 to 120: enough to calibrate a file (100 anchors) in every third of it.
        /// </summary>
        private static List<(string, double)> SearchPeptides(int count)
        {
            const string residues = "ADEFGHILMNQSTVWY";
            var peptides = new List<(string, double)>();
            var seen = new HashSet<string>();
            uint x = 12345;
            while (peptides.Count < count)
            {
                var sequence = new StringBuilder();
                for (int i = 0; i < 7 + peptides.Count % 11; i++)
                {
                    x = x * 1664525 + 1013904223;
                    sequence.Append(residues[(int)(x >> 24) % residues.Length]);
                }
                sequence.Append('K');
                if (seen.Add(sequence.ToString()))
                    peptides.Add((sequence.ToString(), -10 + 130.0 * peptides.Count / (count - 1)));
            }
            return peptides;
        }
    }
    /// <summary>
    /// The .msl update against the live Koina service: Prosit 2025 puts the Biognosys iRT peptides in their
    /// published order. ExternalService, so it runs in its own CI job; skipped when Koina does not answer. Its own
    /// fixture, so the category does not pull the offline tests out of the required job.
    /// </summary>
    [TestFixture]
    [Category("ExternalService")]
    [NonParallelizable]
    [ExcludeFromCodeCoverage]
    public class MslLibraryUpdateLiveTests
    {
        private static readonly (string Sequence, double Irt)[] BiognosysIrt =
        {
            ("LGGNEQVTR", -24.92), ("GAGSSEPVTGLDAK", 0.00), ("VEATFGVDESNAK", 12.39), ("YILAGVENSK", 19.79),
            ("TPVISGGPYEYR", 28.71), ("DGLDAASYYAPVR", 42.26), ("ADVTPADFSEWSK", 54.62), ("GTFIIDPGGVIR", 70.52),
            ("LFLQFGAQGSPFLK", 100.00),
        };

        [Test]
        public void Prosit2025_PutsTheIrtPeptidesInTheirOrder()
        {
            string folder = Path.Combine(TestContext.CurrentContext.TestDirectory, "MslLibraryUpdateLive");
            if (Directory.Exists(folder))
                Directory.Delete(folder, true);
            Directory.CreateDirectory(folder);
            try
            {
                string path = Path.Combine(folder, "library.msl");
                MslLibrary.Save(path, BiognosysIrt.Select(p => MslLibraryUpdateTests.Entry(p.Sequence, ions: 4)).ToList());
                var task = MslLibraryUpdateTests.LibraryUpdateTask(folder, path);

                task.UpdateSpectralLibrary();

                string summary = task.Parameters.SearchTaskResults.ToString();
                if (summary.Contains("No predicted times given"))
                    Assert.Ignore($"Koina unavailable: {summary}");
                Assert.That(summary, Does.Contain("Prosit_2025_irt_40PTM"));
                using var library = MslLibrary.Load(Directory.GetFiles(folder, "updateSpectralLibrary_*").Single());
                var times = library.GetAllEntries().ToDictionary(e => e.FullSequence, e => e.RetentionTime);
                // Prosit 2025's iRT is its own scale, not Biognosys's (on 2026-11-24 it gave LGGNEQVTR -11.1,
                // GAGSSEPVTGLDAK 24.2, LFLQFGAQGSPFLK 137.2), so only the peptides' order is checked.
                var predicted = BiognosysIrt.Select(p => times[p.Sequence]).ToList();
                Assert.That(predicted, Is.All.Not.NaN);
                Assert.That(predicted, Is.Ordered.Ascending, string.Join("; ", predicted));
            }
            finally
            {
                Directory.Delete(folder, true);
            }
        }
    }
}
