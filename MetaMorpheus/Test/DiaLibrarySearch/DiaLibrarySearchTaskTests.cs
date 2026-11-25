#nullable enable
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Net.Http;
using System.Text.RegularExpressions;
using System.Threading;
using System.Threading.Tasks;
using EngineLayer;
using EngineLayer.DatabaseLoading;
using MassSpectrometry;
using Nett;
using Newtonsoft.Json;
using Newtonsoft.Json.Linq;
using NUnit.Framework;
using Omics.Digestion;
using Omics.SpectralMatch.MslSpectralLibrary;
using PredictionClients.Koina.Client;
using PredictionClients.Koina.SupportedModels.FragmentIntensityModels;
using PredictionClients.Koina.SupportedModels.RetentionTimeModels;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Readers;
using Readers.SpectralLibrary;
using TaskLayer;
using Transcriptomics.Digestion;

namespace Test.DiaLibrarySearch;

/// <summary>
/// M7 slice 1 (design/M7.md): the DIA library search as a MetaMorpheus task. It reads and writes TOML like every task,
/// takes the library as a database (`-d library.msl`) with no FASTA, loads DIA scans with every peak, and reports the
/// headline count PXReprise reads from results.txt.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class DiaLibrarySearchTaskTests
{
    private string _directory = "";

    [OneTimeSetUp]
    public void SetUp()
    {
        _directory = Path.Combine(TestContext.CurrentContext.TestDirectory, "DiaLibrarySearchTaskTests");
        Directory.CreateDirectory(_directory);
    }

    [OneTimeTearDown]
    public void TearDown()
    {
        Directory.Delete(_directory, true);
    }

    /// <summary>Writes a synthetic run as mzML and its library as .msl, as a user would hand them to the task.</summary>
    /// <param name="accessions">
    /// Gives each target the protein accession P00000, P00001, ... and its decoy DECOY_ and the same, as mzLib's library
    /// builder names a decoy protein.
    /// </param>
    /// <param name="edit">Changes the library's entries before it is written (the run's scans stay as built).</param>
    private (string Mzml, string Library, SyntheticDiaRun Run) WriteInputs(string name, double noisePeaksPerScan = 60,
        Func<Omics.SpectralMatch.MslSpectralLibrary.MslLibraryEntry, double>? abundance = null, bool accessions = false,
        Action<List<Omics.SpectralMatch.MslSpectralLibrary.MslLibraryEntry>>? edit = null)
    {
        var run = SyntheticDiaRun.Build(300, entry => !entry.IsDecoy && SyntheticDiaRun.Bucket(entry, 4) != 0, withMs1: true, noisePeaksPerScan: noisePeaksPerScan,
            abundance: abundance);
        if (accessions)
        {
            // Build adds each target, then its decoy
            for (int i = 0; i < run.Library.Count; i += 2)
            {
                run.Library[i].ProteinAccession = $"P{i / 2:D5}";
                run.Library[i + 1].ProteinAccession = $"DECOY_P{i / 2:D5}";
            }
        }
        edit?.Invoke(run.Library);
        string folder = Path.Combine(_directory, name);
        Directory.CreateDirectory(folder);
        string mzml = Path.Combine(folder, name + ".mzML");
        MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(new GenericMsDataFile(AsWritten(run.Scans), new SourceFile("no nativeID format", "mzML format", null, null, null)), mzml, false);
        return (mzml, run.WriteLibrary(folder), run);
    }

    /// <summary>
    /// The scans as a converter writes a DIA run: each MS2 names the window centre as its selected ion and the cycle's MS1
    /// as its precursor scan (the mzML writer requires a selected ion).
    /// </summary>
    private static MsDataScan[] AsWritten(MsDataScan[] scans)
    {
        int? lastMs1 = null;
        return scans.Select(scan =>
        {
            if (scan.MsnOrder == 1)
            {
                lastMs1 = scan.OneBasedScanNumber;
                return scan;
            }
            return new MsDataScan(scan.MassSpectrum, scan.OneBasedScanNumber, scan.MsnOrder, scan.IsCentroid, scan.Polarity, scan.RetentionTime,
                scan.ScanWindowRange, scan.ScanFilter, scan.MzAnalyzer, scan.TotalIonCurrent, scan.InjectionTime, scan.NoiseData, scan.NativeId,
                selectedIonMz: scan.IsolationMz, isolationMZ: scan.IsolationMz, isolationWidth: scan.IsolationWidth,
                dissociationType: scan.DissociationType, oneBasedPrecursorScanNumber: lastMs1);
        }).ToArray();
    }

    /// <summary>The task's settings survive a TOML write and read, and the default task is a DIA library search.</summary>
    [Test]
    public void TheTaskRoundTripsThroughToml()
    {
        var task = new DiaLibrarySearchTask(new DiaLibrarySearchTaskParameters
        {
            FragmentTolerancePpm = 15, Ms1TolerancePpm = 4, ClassifierSeed = 3,
            LibraryPrecursorCharges = [2, 3, 4], LibraryNce = 30, LibraryDecoys = DiaLibraryDecoys.Reversed,
        })
        {
            CommonParameters = new CommonParameters(filterThreshold: 0.05),
        };
        string path = Path.Combine(_directory, "DiaLibrarySearchTask.toml");
        Toml.WriteFile(task, path, MetaMorpheusTask.tomlConfig);

        var read = Toml.ReadFile<DiaLibrarySearchTask>(path, MetaMorpheusTask.tomlConfig);

        Assert.That(read.TaskType, Is.EqualTo(MyTask.DiaLibrarySearch));
        Assert.That(read.Parameters.FragmentTolerancePpm, Is.EqualTo(15));
        Assert.That(read.Parameters.Ms1TolerancePpm, Is.EqualTo(4));
        Assert.That(read.Parameters.ClassifierSeed, Is.EqualTo(3));
        Assert.That(read.Parameters.LibraryPrecursorCharges, Is.EqualTo(new[] { 2, 3, 4 }));
        Assert.That(read.Parameters.LibraryNce, Is.EqualTo(30));
        Assert.That(read.Parameters.LibraryDecoys, Is.EqualTo(DiaLibraryDecoys.Reversed));
        var defaults = new DiaLibrarySearchTaskParameters();
        Assert.That(defaults.LibraryPrecursorCharges, Is.EqualTo(new[] { 2, 3 }), "the benchmark libraries' charges");
        Assert.That(defaults.LibraryNce, Is.EqualTo(27));
        Assert.That(defaults.LibraryDecoys, Is.EqualTo(DiaLibraryDecoys.Permuted));
        Assert.That(File.ReadAllText(path), Does.Contain("LibraryPrecursorCharges = [2, 3, 4]"), "an int list in TOML, as LocalFdrCategories is a list");
        Assert.That(read.CommonParameters.FilterThreshold, Is.EqualTo(0.05), "the search's one threshold");
        Assert.That(new DiaLibrarySearchTask().CommonParameters.FilterThreshold, Is.EqualTo(0.01));
        Assert.That(typeof(DiaLibrarySearchTaskParameters).GetProperty("QValueThreshold"), Is.Null,
            "no second threshold: precursors, peptides, protein evidence and protein groups all use FilterThreshold");
        Assert.That(new DiaLibrarySearchTask().CommonParameters, Is.Not.Null);
    }

    /// <summary>
    /// The engine's parameters come from the task's, with everything a user cannot set left at the engine's (benchmarked)
    /// defaults.
    /// </summary>
    [Test]
    public void TheEngineParametersFollowTheTasks()
    {
        var engine = new DiaLibrarySearchTaskParameters { FragmentTolerancePpm = 15, Ms1TolerancePpm = 4, ClassifierSeed = 3 }.ToEngineParameters();
        var defaults = new EngineLayer.DiaLibrarySearch.DiaLibrarySearchParameters();

        Assert.That(engine.FragmentTolerancePpm, Is.EqualTo(15));
        Assert.That(engine.Ms1TolerancePpm, Is.EqualTo(4));
        Assert.That(engine.ClassifierSeed, Is.EqualTo(3));
        Assert.That(engine with { FragmentTolerancePpm = defaults.FragmentTolerancePpm, Ms1TolerancePpm = defaults.Ms1TolerancePpm, ClassifierSeed = 0 },
            Is.EqualTo(defaults));
        Assert.That(new DiaLibrarySearchTaskParameters().ToEngineParameters(), Is.EqualTo(defaults));
    }

    /// <summary>
    /// No peak trimming or filtering on DIA input (user, 2026-09-30): the task keeps every MS1 and MS2 peak, although
    /// CommonParameters' default trims MS2 to the top 200 peaks per window, which would bite on this run.
    /// </summary>
    [Test]
    public void DiaScansAreLoadedWithEveryPeak()
    {
        var (mzml, _, _) = WriteInputs("untrimmed", noisePeaksPerScan: 3000);
        var all = MsDataFileReader.GetDataFile(mzml).LoadAllStaticData().GetAllScansList();

        var loaded = DiaLibrarySearchTask.LoadScans(mzml);

        Assert.That(loaded.Select(s => s.MassSpectrum.Size), Is.EqualTo(all.Select(s => s.MassSpectrum.Size)));
        var trimmed = new MyFileManager(true).LoadFile(mzml, new CommonParameters()).GetAllScansList();
        Assert.That(trimmed.Sum(s => s.MassSpectrum.Size), Is.LessThan(all.Sum(s => s.MassSpectrum.Size)), "the default filter would have trimmed");
    }

    /// <summary>
    /// The task searches a run against a library with no FASTA, end to end, and results.txt carries the headline PXReprise
    /// reads: the target precursors at q &lt;= 0.01, per run and in total. Most planted targets are found and few others.
    /// </summary>
    [Test]
    public void TheTaskSearchesARunAgainstALibraryEndToEnd()
    {
        var (count, results, run, matches) = RunTheTask("endtoend");

        Assert.That(count, Is.EqualTo(DiaLibrarySearchTask.CountPassingTargets(matches, 0.01)), "the headline is the engine's count");
        int planted = run.Library.Count(e => !e.IsDecoy && run.PlantedSequences.Contains(e.FullSequence));
        Assert.That(count, Is.InRange(planted * 8 / 10, planted * 11 / 10));
        Assert.That(results, Does.Contain("endtoend.mzML: " + count + " target precursors with q-value <= 0.01"));
    }

    /// <summary>
    /// Several runs (M7 slice 3, PLAN D4): each row carries its run's q-value and the global one (best score across the runs,
    /// then target-decoy competition); a row is reported when both pass; the headline counts distinct precursors, not rows.
    /// </summary>
    [Test]
    public void SeveralRunsShareAGlobalQValue()
    {
        string folder = Path.Combine(_directory, "tworuns");
        Directory.CreateDirectory(folder);
        string WriteRun(string name, Func<Omics.SpectralMatch.MslSpectralLibrary.MslLibraryEntry, bool> plant, out SyntheticDiaRun run)
        {
            run = SyntheticDiaRun.Build(300, plant, withMs1: true, noisePeaksPerScan: 1500);
            string mzml = Path.Combine(folder, name + ".mzML");
            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(new GenericMsDataFile(AsWritten(run.Scans), new SourceFile("no nativeID format", "mzML format", null, null, null)), mzml, false);
            return mzml;
        }
        string a = WriteRun("runA", e => !e.IsDecoy && SyntheticDiaRun.Bucket(e, 4) != 0, out var runA);
        string b = WriteRun("runB", e => !e.IsDecoy && SyntheticDiaRun.Bucket(e, 3) != 0, out _);
        string library = runA.WriteLibrary(folder);
        string output = Path.Combine(folder, "output", "Task1");
        Directory.CreateDirectory(output);

        new DiaLibrarySearchTask().RunTask(output, [new DbForTask(library, false)], [a, b], "Task1");

        // The same per-run searches by hand, and the global q-values over them
        using var msl = Readers.SpectralLibrary.MslLibrary.Load(library);
        var parameters = new DiaLibrarySearchTaskParameters().ToEngineParameters();
        var byRun = new[] { a, b }.Select(path =>
        {
            var scans = DiaLibrarySearchTask.LoadScans(path);
            var calibration = EngineLayer.DiaLibrarySearch.DiaIrtSelfCalibration.Calibrate(scans, msl, parameters, new CommonParameters());
            return ((EngineLayer.DiaLibrarySearch.DiaLibrarySearchResults)new EngineLayer.DiaLibrarySearch.DiaLibrarySearchEngine(scans, msl,
                calibration.Model, calibration.ApplyTo(parameters), new CommonParameters(), [], []).Run()).Matches;
        }).ToArray();
        var global = EngineLayer.DiaLibrarySearch.DiaPrecursorFdr.GlobalQValues(byRun);

        var table = new DiaPrecursorFile(Path.Combine(output, "AllDiaPrecursors.tsv"));
        table.LoadResults();
        Assert.That(table.Results.Select(r => r.FileName).Distinct(), Is.EquivalentTo(new[] { "runA.mzML", "runB.mzML" }));
        var expected = new[] { ("runA.mzML", byRun[0]), ("runB.mzML", byRun[1]) }
            .SelectMany(r => r.Item2.Where(m => !m.IsDecoy && m.QValue <= 0.01 && global[m.PrecursorIndex] <= 0.01)
                .Select(m => (File: r.Item1, m.FullSequence, m.Charge, Run: m.QValue, Global: global[m.PrecursorIndex])))
            .ToList();
        Assert.That(table.Results.Select(r => (File: r.FileName, r.FullSequence, Charge: r.PrecursorCharge, Run: r.QValuePrecursorRun, Global: r.QValuePrecursorGlobal)),
            Is.EquivalentTo(expected));
        Assert.That(table.Results.Any(r => r.QValuePrecursorGlobal != r.QValuePrecursorRun), "the global level is not the run's");

        string results = File.ReadAllText(Path.Combine(output, "results.txt"));
        int distinct = expected.Select(e => (e.FullSequence, e.Charge)).Distinct().Count();
        Assert.That(distinct, Is.LessThan(expected.Count), "some precursors are reported in both runs");
        Assert.That(results, Does.Contain("All target precursors with q-value <= 0.01: " + distinct));
        Assert.That(results, Does.Contain("runA.mzML: " + expected.Count(e => e.File == "runA.mzML") + " target precursors"));
    }

    /// <summary>
    /// DiaIrtCalibration.tsv (D12): each run's calibration anchors, so a bad fit can be seen. One row per anchor: the run,
    /// its apex in minutes, its library iRT, the calibrated iRT at that apex and the residual; as many rows as anchors.
    /// </summary>
    [Test]
    public void TheTaskWritesEachRunsCalibration()
    {
        var (mzml, library, _) = WriteInputs("calibration");
        string output = Path.Combine(_directory, "calibration", "output", "Task1");
        Directory.CreateDirectory(output);
        new DiaLibrarySearchTask().RunTask(output, [new DbForTask(library, false)], [mzml], "Task1");

        var lines = File.ReadAllLines(Path.Combine(output, "DiaIrtCalibration.tsv"));
        Assert.That(lines[0].Split('\t'), Is.EqualTo(new[] { "File Name", "AnchorRtMin", "LibraryIrt", "CalibratedIrt", "ResidualIrt" }));

        using var msl = Readers.SpectralLibrary.MslLibrary.Load(library);
        var calibration = EngineLayer.DiaLibrarySearch.DiaIrtSelfCalibration.Calibrate(DiaLibrarySearchTask.LoadScans(mzml), msl,
            new DiaLibrarySearchTaskParameters().ToEngineParameters(), new CommonParameters());
        Assert.That(lines.Length - 1, Is.EqualTo(calibration.Anchors.Count).And.GreaterThan(0));
        for (int i = 0; i < calibration.Anchors.Count; i++)
        {
            var cells = lines[i + 1].Split('\t');
            var (rt, irt) = calibration.Anchors[i];
            double calibrated = calibration.Model.ToIrt(rt).Value;
            Assert.That(cells[0], Is.EqualTo("calibration.mzML"));
            Assert.That(double.Parse(cells[1], System.Globalization.CultureInfo.InvariantCulture), Is.EqualTo(rt.Value));
            Assert.That(double.Parse(cells[2], System.Globalization.CultureInfo.InvariantCulture), Is.EqualTo(irt.Value));
            Assert.That(double.Parse(cells[3], System.Globalization.CultureInfo.InvariantCulture), Is.EqualTo(calibrated));
            Assert.That(double.Parse(cells[4], System.Globalization.CultureInfo.InvariantCulture), Is.EqualTo(irt.Value - calibrated).Within(1e-9));
        }
    }

    /// <summary>
    /// DiaCalibrationSummary.tsv: what each run's calibration settled on, one row per run, so a record of the search can say
    /// which window and tolerances it searched with (PXReprise DIA-Q6). The anchors and rounds; the residual SD in iRT and
    /// the iRT half-window it sets; the MS1 offset applied (its median over the retention time) and its residual spread; and
    /// the MS1 and fragment tolerances the search used.
    /// </summary>
    [Test]
    public void TheTaskWritesEachRunsCalibrationSummary()
    {
        var (mzml, library, _) = WriteInputs("summary");
        string output = Path.Combine(_directory, "summary", "output", "Task1");
        Directory.CreateDirectory(output);
        new DiaLibrarySearchTask().RunTask(output, [new DbForTask(library, false)], [mzml], "Task1");

        var lines = File.ReadAllLines(Path.Combine(output, "DiaCalibrationSummary.tsv"));
        Assert.That(lines[0].Split('\t'), Is.EqualTo(new[] { "File Name", "Anchors", "Rounds", "ResidualSdIrt", "IrtHalfWindow",
            "Ms1OffsetPpm", "Ms1SpreadPpm", "Ms1TolerancePpm", "FragmentTolerancePpm" }));
        Assert.That(lines, Has.Length.EqualTo(2), "one row per run");

        using var msl = Readers.SpectralLibrary.MslLibrary.Load(library);
        var parameters = new DiaLibrarySearchTaskParameters().ToEngineParameters();
        var calibration = EngineLayer.DiaLibrarySearch.DiaIrtSelfCalibration.Calibrate(DiaLibrarySearchTask.LoadScans(mzml), msl, parameters, new CommonParameters());
        Assert.That(calibration.Ms1Offset, Is.Not.Null, "the synthetic run calibrates over several rounds, so it has an MS1 offset to report");
        var cells = lines[1].Split('\t');
        double Cell(int i) => double.Parse(cells[i], System.Globalization.CultureInfo.InvariantCulture);
        Assert.That(cells[0], Is.EqualTo("summary.mzML"));
        Assert.That(int.Parse(cells[1]), Is.EqualTo(calibration.AnchorCount).And.GreaterThan(0));
        Assert.That(int.Parse(cells[2]), Is.EqualTo(calibration.Rounds));
        Assert.That(Cell(3), Is.EqualTo(calibration.Model.ResidualSd));
        Assert.That(Cell(4), Is.EqualTo(calibration.IrtHalfWindow));
        Assert.That(Cell(5), Is.EqualTo(calibration.Ms1Offset!.MedianOffsetPpm));
        Assert.That(Cell(6), Is.EqualTo(calibration.Ms1Offset!.ResidualSpreadPpm));
        Assert.That(Cell(7), Is.EqualTo(calibration.ApplyTo(parameters).EffectiveMs1TolerancePpm));
        Assert.That(Cell(8), Is.EqualTo(parameters.FragmentTolerancePpm));
    }

    /// <summary>
    /// A summary cell with nothing to report is left empty rather than written as NaN: no MS1 offset (a single-round
    /// calibration), or one with no measured spread.
    /// </summary>
    [Test]
    public void ACalibrationSummaryLeavesWhatWasNotMeasuredEmpty()
    {
        var model = Chromatography.RetentionTimeCalibration.IrtCalibration.Line((new MzLibUtil.RtMinutes(1), new MzLibUtil.Irt(0)), (new MzLibUtil.RtMinutes(10), new MzLibUtil.Irt(100)));
        var single = new EngineLayer.DiaLibrarySearch.DiaIrtCalibration(model, 12.5, 600);
        var header = DiaLibrarySearchTask.CalibrationSummaryHeader.Split('\t');
        int offset = Array.IndexOf(header, "Ms1OffsetPpm"), spread = Array.IndexOf(header, "Ms1SpreadPpm");
        var row = DiaLibrarySearchTask.CalibrationSummaryRow("a.mzML", single, new DiaLibrarySearchTaskParameters().ToEngineParameters()).Split('\t');
        Assert.That(row, Has.Length.EqualTo(header.Length));
        Assert.That(row[offset], Is.Empty, "no MS1 offset");
        Assert.That(row[spread], Is.Empty, "no MS1 offset, so no spread");

        var constant = single with { Ms1Offset = EngineLayer.DiaLibrarySearch.Ms1OffsetModel.Constant(1.5) };
        row = DiaLibrarySearchTask.CalibrationSummaryRow("a.mzML", constant, new DiaLibrarySearchTaskParameters().ToEngineParameters()).Split('\t');
        Assert.That(double.Parse(row[offset], System.Globalization.CultureInfo.InvariantCulture), Is.EqualTo(1.5), "a constant offset");
        Assert.That(row[spread], Is.Empty, "a constant offset has no measured spread");
    }

    /// <summary>
    /// AllDiaPrecursors.tsv (M7 slice 2) is written through mzLib's DiaPrecursorFile, in the schema agreed with dataRepo: one
    /// row per target precursor at the threshold, so its rows are the headline; the run key verbatim; each row's apex scan,
    /// retention in minutes and iRT, and quantity; and, with no contaminant database, "contaminants: not assessed".
    /// </summary>
    [Test]
    public void TheTaskWritesThePrecursorTable()
    {
        var (count, _, run, matches) = RunTheTask("table");
        string path = Path.Combine(_directory, "table", "output", "Task1", "AllDiaPrecursors.tsv");

        var table = new DiaPrecursorFile(path);
        table.LoadResults();

        Assert.That(table.ContaminantsAssessed, Is.False);
        Assert.That(table.Results.Count, Is.EqualTo(count));
        var bySequence = matches.Where(m => !m.IsDecoy).ToDictionary(m => (m.FullSequence, m.Charge));
        foreach (var row in table.Results)
        {
            var match = bySequence[(row.FullSequence, row.PrecursorCharge)];
            Assert.That(row.FileName, Is.EqualTo("table.mzML"));
            Assert.That(row.Label, Is.EqualTo("T"));
            Assert.That(row.QValuePrecursorRun, Is.EqualTo(match.QValue).And.LessThanOrEqualTo(0.01));
            Assert.That(row.QValuePrecursorGlobal, Is.EqualTo(row.QValuePrecursorRun), "one run: global is the run's");
            Assert.That(row.ApexScanNumber, Is.EqualTo(match.ApexScanNumber).And.GreaterThan(0));
            Assert.That(row.ApexRtMin, Is.EqualTo(match.ApexRt.Value));
            Assert.That(row.ApexIrt, Is.EqualTo(match.ApexIrt.Value));
            Assert.That(row.LibraryIrt, Is.EqualTo(match.LibraryIrt.Value));
            Assert.That(row.PrecursorMz, Is.EqualTo(match.PrecursorMz));
            Assert.That(row.Score, Is.EqualTo(match.Score));
            Assert.That(row.PrecursorQuantity, Is.EqualTo(double.IsNaN(match.Quantity) ? null : match.Quantity));
            var entry = run.Library.Single(e => e.FullSequence == row.FullSequence && e.ChargeState == row.PrecursorCharge);
            Assert.That(row.BaseSequence, Is.EqualTo(entry.BaseSequence));
            Assert.That(row.ProteinAccession, Is.EqualTo(entry.ProteinAccession ?? ""), "no accession is an empty cell");
        }
    }

    /// <summary>
    /// AllDiaPeptides.tsv (PLAN D4, "peptide, global"; mzLib's DiaPeptideFile): the target peptides at the threshold from a
    /// fresh competition among peptides, each with its best precursor's run and charge, and their count in results.txt.
    /// With no protein database there is no protein-group table and no protein line.
    /// </summary>
    [Test]
    public void TheTaskWritesThePeptideTable()
    {
        var (_, results, run, matches) = RunTheTask("peptides");
        string output = Path.Combine(_directory, "peptides", "output", "Task1");

        var table = new DiaPeptideFile(Path.Combine(output, "AllDiaPeptides.tsv"));
        table.LoadResults();

        var expected = EngineLayer.DiaLibrarySearch.DiaPeptideFdr.Assign(matches).Where(p => !p.IsDecoy && p.QValue <= 0.01).ToList();
        Assert.That(expected, Is.Not.Empty);
        Assert.That(table.ContaminantsAssessed, Is.False);
        Assert.That(table.Results.Select(r => (r.FullSequence, r.Score, r.QValuePeptideGlobal)),
            Is.EquivalentTo(expected.Select(p => (p.FullSequence, p.Score, p.QValue))));
        var precursors = matches.ToDictionary(m => m.PrecursorIndex);
        foreach (var row in table.Results)
        {
            var best = precursors[expected.Single(p => p.FullSequence == row.FullSequence).BestPrecursorIndex];
            Assert.That(row.FileName, Is.EqualTo("peptides.mzML"));
            Assert.That(row.PrecursorCharge, Is.EqualTo(best.Charge));
            Assert.That(row.Label, Is.EqualTo("T"));
            Assert.That(row.BaseSequence, Is.EqualTo(run.Library.First(e => e.FullSequence == row.FullSequence).BaseSequence));
        }
        Assert.That(results, Does.Contain("All target peptides with q-value <= 0.01: " + expected.Count));
        Assert.That(results, Does.Not.Contain("protein groups"));
        Assert.That(File.Exists(Path.Combine(output, "AllDiaProteinGroups.tsv")), Is.False);
    }

    /// <summary>
    /// Given a protein database beside the library, the task writes AllDiaProteinGroups.tsv in MetaMorpheus's protein-group
    /// format (read back by mzLib's protein-group reader): the target groups at the threshold that MetaMorpheus's parsimony
    /// and picked protein FDR infer from the peptides, decoys placed on the decoy proteins the library names. Precursors,
    /// peptides and groups whose proteins are all in the contaminant database are labelled C, and the tables say
    /// contaminants were assessed.
    /// </summary>
    [Test]
    public void GivenAProteinDatabaseTheTaskWritesProteinGroups()
    {
        // Enough noise that decoys find peaks and are scored, as in a real run
        var (mzml, library, run) = WriteInputs("proteins", noisePeaksPerScan: 1500, accessions: true);
        string folder = Path.GetDirectoryName(mzml)!;
        var targets = run.Library.Where(e => !e.IsDecoy).ToList();
        bool IsContaminant(Omics.SpectralMatch.MslSpectralLibrary.MslLibraryEntry e) => SyntheticDiaRun.Bucket(e, 5) == 0;
        string fasta = WriteFasta(Path.Combine(folder, "proteins.fasta"), targets.Where(e => !IsContaminant(e)));
        string contaminants = WriteFasta(Path.Combine(folder, "contaminants.fasta"), targets.Where(IsContaminant));
        var contaminantAccessions = targets.Where(IsContaminant).Select(e => e.ProteinAccession).ToHashSet();
        string output = Path.Combine(folder, "output", "Task1");
        Directory.CreateDirectory(output);
        var commonParameters = new CommonParameters(listOfModsFixed: new List<(string, string)>());

        new DiaLibrarySearchTask { CommonParameters = commonParameters }.RunTask(output,
            [new DbForTask(library, false), new DbForTask(fasta, false), new DbForTask(contaminants, true)], [mzml], "Task1");

        // The same inference by hand, on the hand-run engine's peptides
        var scans = DiaLibrarySearchTask.LoadScans(mzml);
        using var msl = Readers.SpectralLibrary.MslLibrary.Load(library);
        var parameters = new DiaLibrarySearchTaskParameters().ToEngineParameters();
        var calibration = EngineLayer.DiaLibrarySearch.DiaIrtSelfCalibration.Calibrate(scans, msl, parameters, commonParameters);
        var matches = ((EngineLayer.DiaLibrarySearch.DiaLibrarySearchResults)new EngineLayer.DiaLibrarySearch.DiaLibrarySearchEngine(scans, msl,
            calibration.Model, calibration.ApplyTo(parameters), commonParameters, [], []).Run()).Matches;
        var peptides = EngineLayer.DiaLibrarySearch.DiaPeptideFdr.Assign(matches);
        var proteins = UsefulProteomicsDatabases.ProteinDbLoader.LoadProteinFasta(fasta, true, UsefulProteomicsDatabases.DecoyType.None, false, out _)
            .Concat(UsefulProteomicsDatabases.ProteinDbLoader.LoadProteinFasta(contaminants, true, UsefulProteomicsDatabases.DecoyType.None, true, out _)).ToList();
        var variable = GlobalVariables.AllModsKnown.Where(m => commonParameters.ListOfModsVariable.Contains((m.ModificationType, m.IdWithMotif))).ToList();
        var decoyAccessions = run.Library.Where(e => e.IsDecoy).ToDictionary(e => e.FullSequence, e => (IReadOnlyList<string>)[e.ProteinAccession]);
        var inference = EngineLayer.DiaLibrarySearch.DiaProteinInference.Infer(peptides, proteins, commonParameters, mzml, [], variable,
            decoyAccessionsOf: s => decoyAccessions.GetValueOrDefault(s));
        var expected = inference.ProteinGroups.Where(g => !g.IsDecoy && g.QValue <= 0.01).ToList();
        Assert.That(expected.Count, Is.GreaterThan(10));
        Assert.That(inference.SpectralMatches.Count(m => m.IsDecoy), Is.GreaterThan(0),
            $"decoys were placed: {peptides.Count(p => p.IsDecoy)} decoy peptides, {inference.UnmappedPeptides} unmapped, {inference.SpectralMatches.Count} matches, " +
            $"{peptides.Where(p => p.IsDecoy).Count(p => decoyAccessions.ContainsKey(p.FullSequence))} decoys with accessions");

        string groupsPath = Path.Combine(output, "AllDiaProteinGroups.tsv");
        var groups = FileReader.ReadFile<ProteinGroupFromTsvFile>(groupsPath).Results;
        Assert.That(groups.Select(g => (g.ProteinGroupName, g.QValue, g.DecoyContaminantTarget)),
            Is.EquivalentTo(expected.Select(g => (g.ProteinGroupName, g.QValue, g.IsContaminant ? "C" : "T"))));
        Assert.That(groups.Where(g => g.DecoyContaminantTarget == "C").Select(g => g.ProteinGroupName), Is.Not.Empty.And.All.Matches<string>(contaminantAccessions.Contains));
        string results = File.ReadAllText(Path.Combine(output, "results.txt"));
        Assert.That(results, Does.Contain("All target protein groups with q-value <= 0.01: " + groups.Count));

        var precursorTable = new DiaPrecursorFile(Path.Combine(output, "AllDiaPrecursors.tsv"));
        precursorTable.LoadResults();
        var peptideTable = new DiaPeptideFile(Path.Combine(output, "AllDiaPeptides.tsv"));
        peptideTable.LoadResults();
        Assert.That(precursorTable.ContaminantsAssessed && peptideTable.ContaminantsAssessed);
        foreach (var (accession, label) in precursorTable.Results.Select(r => (r.ProteinAccession, r.Label)).Concat(peptideTable.Results.Select(r => (r.ProteinAccession, r.Label))))
            Assert.That(label, Is.EqualTo(contaminantAccessions.Contains(accession) ? "C" : "T"), accession);
        Assert.That(peptideTable.Results.Count(r => r.Label == "C"), Is.GreaterThan(0));
    }

    /// <summary>
    /// The task opens every library index-only, whatever its size. Reading entries on demand costs a search no time it
    /// would notice, and a full load keeps several times the file in memory: 7.1 GB against 1.1 GB for a 2 GB library
    /// (results/msl_in_dia).
    /// </summary>
    [Test]
    public void TheTaskOpensEveryLibraryIndexOnly()
    {
        var (_, library, run) = WriteInputs("indexOnly");

        using var opened = DiaLibrarySearchTask.OpenLibrary(library);

        Assert.That(opened.IsIndexOnly, Is.True);
        Assert.That(opened.PrecursorCount, Is.EqualTo(run.Library.Count));
    }

    /// <summary>
    /// A library whose retention times are run minutes (.msl rt_is_calibrated) is refused before any run is read. The
    /// search calibrates each run onto library iRT, so it would read minutes as iRT without a word.
    /// </summary>
    [Test]
    public void ALibraryInRunMinutesIsRefused()
    {
        var (mzml, library, run) = WriteInputs("runMinutes", edit: entries => entries[0].RtIsCalibrated = true);
        string output = Path.Combine(_directory, "runMinutes", "output", "Task1");
        Directory.CreateDirectory(output);

        var refusal = Assert.Throws<MetaMorpheusException>(() => DiaLibrarySearchTask.OpenLibrary(library));
        Assert.That(refusal!.Message, Does.Contain($"1 of {run.Library.Count}").And.Contain("iRT"));
        Assert.That(() => new DiaLibrarySearchTask().RunTask(output, [new DbForTask(library, false)], [mzml], "Task1"),
            Throws.InstanceOf<MetaMorpheusException>().With.Message.Contains("iRT"));
    }

    /// <summary>
    /// The Decoy/Contaminant/Target label is mzLib's (DecoyContaminantTargetLabel). An entrapment precursor (.msl flag bit 3)
    /// is ET, or ED for a decoy, ahead of D and C. A target whose proteins are all contaminants is C. Any other is T.
    /// </summary>
    [TestCase(false, false, "P1", "T")]
    [TestCase(true, false, "DECOY_P1", "D")]
    [TestCase(false, true, "Random_P1", "ET")]
    [TestCase(true, true, "DECOY_Random_P1", "ED")]
    [TestCase(false, false, "C1", "C")]
    [TestCase(false, false, "C1|P1", "T")]
    [TestCase(false, true, "C1", "ET")]
    public void TheLabelIsMzLibsWithEntrapment(bool isDecoy, bool isEntrapment, string accession, string label) =>
        Assert.That(DiaLibrarySearchTask.Label(isDecoy, isEntrapment, accession, new HashSet<string> { "C1" }), Is.EqualTo(label));

    /// <summary>
    /// Entrapment precursors in the library (.msl flag bit 3, set by mzLib's library builder) are searched as targets and
    /// labelled ET in the precursor and peptide tables, so an entrapment check can read the task's output as it is.
    /// </summary>
    [Test]
    public void EntrapmentPrecursorsAreLabelledEtInTheTables()
    {
        static bool IsEntrapment(Omics.SpectralMatch.MslSpectralLibrary.MslLibraryEntry e) => SyntheticDiaRun.Bucket(e, 5) == 1;
        var (mzml, library, run) = WriteInputs("entrapment", edit: entries =>
        {
            foreach (var e in entries.Where(IsEntrapment))
                e.IsEntrapment = true;
        });
        string output = Path.Combine(_directory, "entrapment", "output", "Task1");
        Directory.CreateDirectory(output);

        new DiaLibrarySearchTask().RunTask(output, [new DbForTask(library, false)], [mzml], "Task1");

        var precursors = new DiaPrecursorFile(Path.Combine(output, "AllDiaPrecursors.tsv"));
        precursors.LoadResults();
        var peptides = new DiaPeptideFile(Path.Combine(output, "AllDiaPeptides.tsv"));
        peptides.LoadResults();
        var entryOf = run.Library.Where(e => !e.IsDecoy).ToDictionary(e => e.FullSequence + "/" + e.ChargeState);
        foreach (var (sequence, charge, label) in precursors.Results.Select(r => (r.FullSequence, r.PrecursorCharge, r.Label))
            .Concat(peptides.Results.Select(r => (r.FullSequence, r.PrecursorCharge, r.Label))))
            Assert.That(label, Is.EqualTo(IsEntrapment(entryOf[sequence + "/" + charge]) ? "ET" : "T"), sequence);
        Assert.That(precursors.Results.Count(r => r.Label == "ET"), Is.GreaterThan(0));
        Assert.That(peptides.Results.Count(r => r.Label == "ET"), Is.GreaterThan(0));
    }

    /// <summary>
    /// A library whose decoys name no protein gives protein FDR no decoys, so every group would pass. The task warns and
    /// writes no protein-group table rather than q-values with no decoy evidence; the peptide table is still written.
    /// </summary>
    [Test]
    public void ALibraryWhoseDecoysNameNoProteinWritesNoProteinGroups()
    {
        var (mzml, library, run) = WriteInputs("noDecoyProteins", noisePeaksPerScan: 1500);
        string folder = Path.GetDirectoryName(mzml)!;
        string fasta = WriteFasta(Path.Combine(folder, "proteins.fasta"), run.Library.Where(e => !e.IsDecoy));
        string output = Path.Combine(folder, "output", "Task1");
        Directory.CreateDirectory(output);
        var warnings = new List<string>();
        EventHandler<StringEventArgs> onWarn = (_, e) => warnings.Add(e.S);
        MetaMorpheusTask.WarnHandler += onWarn;
        try
        {
            new DiaLibrarySearchTask { CommonParameters = new CommonParameters(listOfModsFixed: new List<(string, string)>()) }
                .RunTask(output, [new DbForTask(library, false), new DbForTask(fasta, false)], [mzml], "Task1");
        }
        finally
        {
            MetaMorpheusTask.WarnHandler -= onWarn;
        }

        // Decoys were scored, so it is their missing accessions, not their absence, that stops the protein groups
        var warning = warnings.Single(w => w.Contains("No decoy peptide names a decoy protein"));
        Assert.That(int.Parse(Regex.Match(warning, @"of (\d+) decoy peptides scored").Groups[1].Value), Is.GreaterThan(0), warning);
        Assert.That(File.Exists(Path.Combine(output, "AllDiaProteinGroups.tsv")), Is.False);
        Assert.That(File.ReadAllText(Path.Combine(output, "results.txt")), Does.Not.Contain("protein groups"));
        var peptides = new DiaPeptideFile(Path.Combine(output, "AllDiaPeptides.tsv"));
        peptides.LoadResults();
        Assert.That(peptides.Results, Is.Not.Empty);
    }

    /// <summary>
    /// The task's one threshold is the search's FilterThreshold (MetaMorpheus #2982: one identification filter, used at
    /// every step), so protein evidence and the groups written are cut alike (acesnik's #2984 review), and results.txt
    /// says which threshold it used.
    /// </summary>
    [Test]
    public void TheSearchsFilterThresholdIsTheTasksOneThreshold()
    {
        // Enough noise that decoys are scored and some targets fall between 1% and 5%
        var (mzml, library, _) = WriteInputs("threshold5", noisePeaksPerScan: 1500);
        List<double> Run(double threshold, string folder)
        {
            string output = Path.Combine(_directory, "threshold5", folder, "Task1");
            Directory.CreateDirectory(output);
            new DiaLibrarySearchTask { CommonParameters = new CommonParameters(filterThreshold: threshold) }
                .RunTask(output, [new DbForTask(library, false)], [mzml], "Task1");
            Assert.That(File.ReadAllText(Path.Combine(output, "results.txt")),
                Does.Contain($"All target precursors with q-value <= {threshold.ToString(System.Globalization.CultureInfo.InvariantCulture)}: "));
            var table = new DiaPrecursorFile(Path.Combine(output, "AllDiaPrecursors.tsv"));
            table.LoadResults();
            return table.Results.Select(r => Math.Max(r.QValuePrecursorRun, r.QValuePrecursorGlobal)).ToList();
        }

        var atOne = Run(0.01, "one");
        var atFive = Run(0.05, "five");

        Assert.That(atOne, Is.Not.Empty.And.All.LessThanOrEqualTo(0.01));
        Assert.That(atFive, Is.All.LessThanOrEqualTo(0.05));
        Assert.That(atFive.Count(q => q > 0.01), Is.GreaterThan(0), "the fixture has precursors between 1% and 5%");
        Assert.That(atFive.Count, Is.EqualTo(atOne.Count + atFive.Count(q => q > 0.01)));
    }

    /// <summary>
    /// Targets the database does not yield (a modification or cleavage the settings do not allow, or a library built from
    /// another database) are no protein evidence, while their decoys still map by the library's accessions, so the target
    /// side alone shrinks. The task says how many (acesnik's #2984 review), and says so too when none maps at all, rather
    /// than write an empty protein table in silence.
    /// </summary>
    [TestCase(true)]
    [TestCase(false)]
    public void TargetsTheDatabaseDoesNotYieldAreWarnedAbout(bool half)
    {
        string name = half ? "unmappedHalf" : "unmappedAll";
        var (mzml, library, run) = WriteInputs(name, noisePeaksPerScan: 1500, accessions: true);
        string folder = Path.GetDirectoryName(mzml)!;
        var targets = run.Library.Where(e => !e.IsDecoy).ToList();
        string fasta = half
            ? WriteFasta(Path.Combine(folder, "half.fasta"), targets.Where((_, i) => i % 2 == 0))
            : WriteUnrelatedFasta(Path.Combine(folder, "unrelated.fasta"));
        string output = Path.Combine(folder, "output", "Task1");
        Directory.CreateDirectory(output);
        var warnings = new List<string>();
        EventHandler<StringEventArgs> onWarn = (_, e) => warnings.Add(e.S);
        MetaMorpheusTask.WarnHandler += onWarn;
        try
        {
            new DiaLibrarySearchTask { CommonParameters = new CommonParameters(listOfModsFixed: new List<(string, string)>()) }
                .RunTask(output, [new DbForTask(library, false), new DbForTask(fasta, false)], [mzml], "Task1");
        }
        finally
        {
            MetaMorpheusTask.WarnHandler -= onWarn;
        }

        var warning = warnings.Single(w => w.Contains("target peptides map to no protein"));
        var counts = Regex.Match(warning, @"(\d+) of (\d+) target peptides");
        int unmapped = int.Parse(counts.Groups[1].Value), scored = int.Parse(counts.Groups[2].Value);
        Assert.That(unmapped, Is.GreaterThan(0), warning);
        if (half)
            Assert.That(unmapped, Is.LessThan(scored), warning);
        else
            Assert.That(unmapped, Is.EqualTo(scored), warning);
    }

    /// <summary>
    /// One run calibration cannot use (here, one with no DIA MS2 scans) is skipped with a warning naming it and why, and
    /// the other runs are searched and written; results.txt lists the run skipped (acesnik's #2984 review). Only when no
    /// run can be searched does the task fail, naming each run and its reason.
    /// </summary>
    [Test]
    public void ARunCalibrationCannotUseIsSkippedAndTheOthersAreWritten()
    {
        var (mzml, library, run) = WriteInputs("badRun");
        string folder = Path.GetDirectoryName(mzml)!;
        string ms1Only = Path.Combine(folder, "ms1only.mzML");
        MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(new GenericMsDataFile(run.Scans.Where(s => s.MsnOrder == 1).ToArray(),
            new SourceFile("no nativeID format", "mzML format", null, null, null)), ms1Only, false);
        string output = Path.Combine(folder, "output", "Task1");
        Directory.CreateDirectory(output);
        var warnings = new List<string>();
        EventHandler<StringEventArgs> onWarn = (_, e) => warnings.Add(e.S);
        MetaMorpheusTask.WarnHandler += onWarn;
        try
        {
            new DiaLibrarySearchTask().RunTask(output, [new DbForTask(library, false)], [mzml, ms1Only], "Task1");
        }
        finally
        {
            MetaMorpheusTask.WarnHandler -= onWarn;
        }

        Assert.That(warnings.Single(w => w.Contains("ms1only.mzML")), Does.Contain("no DIA MS2 scans"));
        string results = File.ReadAllText(Path.Combine(output, "results.txt"));
        Assert.That(results, Does.Contain("ms1only.mzML: skipped").And.Contain("no DIA MS2 scans"));
        Assert.That(int.Parse(Regex.Match(results, @"All target precursors with q-value <= 0\.01: (\d+)").Groups[1].Value), Is.GreaterThan(0));
        var table = new DiaPrecursorFile(Path.Combine(output, "AllDiaPrecursors.tsv"));
        table.LoadResults();
        Assert.That(table.Results.Select(r => r.FileName).Distinct(), Is.EqualTo(new[] { "badRun.mzML" }));
        Assert.That(File.ReadAllLines(Path.Combine(output, "DiaCalibrationSummary.tsv")).Skip(1).Select(l => l.Split('\t')[0]),
            Is.EqualTo(new[] { "badRun.mzML" }));

        string nothingUsable = Path.Combine(folder, "output", "Task2");
        Directory.CreateDirectory(nothingUsable);
        Assert.That(() => new DiaLibrarySearchTask().RunTask(nothingUsable, [new DbForTask(library, false)], [ms1Only], "Task2"),
            Throws.InstanceOf<MetaMorpheusException>().With.Message.Contains("ms1only.mzML").And.Message.Contains("no DIA MS2 scans"));
    }

    /// <summary>A FASTA of one protein no library target comes from.</summary>
    private static string WriteUnrelatedFasta(string path)
    {
        File.WriteAllLines(path, new[] { ">sp|U00001|UNREL_SYNTH Unrelated protein OS=Synthetic OX=1 GN=UNREL", "MWWHHCCYYWWHHCCYYK" });
        return path;
    }

    /// <summary>
    /// A FASTA with one protein per library target, its accession the entry's (or one made from its position), so digesting
    /// it yields exactly the library's target peptides.
    /// </summary>
    private static string WriteFasta(string path, IEnumerable<Omics.SpectralMatch.MslSpectralLibrary.MslLibraryEntry> targets)
    {
        File.WriteAllLines(path, targets.SelectMany((e, i) =>
        {
            string accession = string.IsNullOrEmpty(e.ProteinAccession) ? $"Q{i:D5}" : e.ProteinAccession;
            return new[] { $">sp|{accession}|{accession}_SYNTH Synthetic protein OS=Synthetic OX=1 GN=G{accession}", e.BaseSequence };
        }));
        return path;
    }

    /// <summary>
    /// The headline counts targets at the q-value threshold, not every target scored, and never a decoy (whose q is NaN).
    /// </summary>
    [Test]
    public void TheHeadlineCountsOnlyTargetsPassingTheThreshold()
    {
        static EngineLayer.DiaLibrarySearch.DiaPrecursorMatch Match(bool decoy, double q) =>
            new(0, "PEPTIDE", 2, 400, decoy, new MzLibUtil.Irt(0), new MzLibUtil.RtMinutes(1), new MzLibUtil.Irt(0), 1, [], q);

        var matches = new[] { Match(false, 0.001), Match(false, 0.01), Match(false, 0.02), Match(true, double.NaN), Match(true, 0.0) };

        Assert.That(DiaLibrarySearchTask.CountPassingTargets(matches, 0.01), Is.EqualTo(2));
        Assert.That(DiaLibrarySearchTask.CountPassingTargets(matches, 0.05), Is.EqualTo(3));
        Assert.That(DiaLibrarySearchTask.PassingTargets(matches, 0.01), Is.EqualTo(new[] { matches[0], matches[1] }), "the table's rows");
    }

    /// <summary>A run's passing target is reported only when its global q-value passes too (PLAN D4).</summary>
    [Test]
    public void ARunsPassingTargetIsReportedOnlyWhenItPassesGlobally()
    {
        static EngineLayer.DiaLibrarySearch.DiaPrecursorMatch Match(int index, double q) =>
            new(index, "PEPTIDE" + index, 2, 400, false, new MzLibUtil.Irt(0), new MzLibUtil.RtMinutes(1), new MzLibUtil.Irt(0), 1, [], q);
        var runPassing = new[] { Match(1, 0.001), Match(2, 0.005), Match(3, 0.009) };
        var global = new Dictionary<int, double> { [1] = 0.002, [2] = 0.03, [3] = 0.01 };

        Assert.That(DiaLibrarySearchTask.Reported(runPassing, global, 0.01), Is.EqualTo(new[] { runPassing[0], runPassing[2] }));
    }

    /// <summary>
    /// Runs the task on a synthetic run, and the same pipeline by hand on the same files. Returns the headline count,
    /// results.txt, the run, and the hand-run engine's matches.
    /// </summary>
    private (int Count, string Results, SyntheticDiaRun Run, List<EngineLayer.DiaLibrarySearch.DiaPrecursorMatch> Matches) RunTheTask(string name,
        Func<Omics.SpectralMatch.MslSpectralLibrary.MslLibraryEntry, double>? abundance = null)
    {
        var (mzml, library, run) = WriteInputs(name, abundance: abundance);
        string output = Path.Combine(_directory, name, "output", "Task1");
        Directory.CreateDirectory(output);

        new DiaLibrarySearchTask().RunTask(output, [new DbForTask(library, false)], [mzml], "Task1");

        string results = File.ReadAllText(Path.Combine(output, "results.txt"));
        var headline = Regex.Match(results, @"All target precursors with q-value <= 0\.01: (\d+)");
        Assert.That(headline.Success, results);

        var scans = DiaLibrarySearchTask.LoadScans(mzml);
        using var msl = Readers.SpectralLibrary.MslLibrary.Load(library);
        var parameters = new DiaLibrarySearchTaskParameters().ToEngineParameters();
        var calibration = EngineLayer.DiaLibrarySearch.DiaIrtSelfCalibration.Calibrate(scans, msl, parameters, new CommonParameters());
        var matches = ((EngineLayer.DiaLibrarySearch.DiaLibrarySearchResults)new EngineLayer.DiaLibrarySearch.DiaLibrarySearchEngine(scans, msl,
            calibration.Model, calibration.ApplyTo(parameters), new CommonParameters(), [], []).Run()).Matches;
        return (int.Parse(headline.Groups[1].Value), results, run, matches);
    }

    /// <summary>
    /// The runner refuses a search with no protein database, but a DIA library search needs only its library: a .msl alone
    /// must reach the task.
    /// </summary>
    [Test]
    public void TheRunnerAcceptsALibraryAloneForTheDiaTask()
    {
        var (mzml, library, _) = WriteInputs("runner");
        string output = Path.Combine(_directory, "runner", "output");
        var warnings = new List<string>();
        EventHandler<StringEventArgs> onWarn = (_, e) => warnings.Add(e.S);
        EverythingRunnerEngine.WarnHandler += onWarn;
        try
        {
            new EverythingRunnerEngine([("Task1", new DiaLibrarySearchTask())], [mzml], [new DbForTask(library, false)], output).Run();
        }
        finally
        {
            EverythingRunnerEngine.WarnHandler -= onWarn;
        }

        Assert.That(warnings, Has.None.Contains("No protein database"));
        Assert.That(File.Exists(Path.Combine(output, "Task1", "results.txt")));
    }

    /// <summary>
    /// Koina answered from memory, the same in every build: each peptide's fragments and iRT come from a seed its sequence
    /// gives, so a library predicted here and one the task predicts are the same library.
    /// </summary>
    internal sealed class CannedProsit
    {
        public bool Unreachable { get; set; }
        public int Calls;
        public ConcurrentBag<string> Requested { get; } = [];
        public FakeIntensity Intensity { get; }
        public FakeIrt Irt { get; }

        public CannedProsit()
        {
            Intensity = new FakeIntensity(this);
            Irt = new FakeIrt(this);
        }

        /// <summary>The residues of a peptide as Koina is sent it, modifications and terminal marks left out.</summary>
        public static string Residues(string koinaSequence) => Regex.Replace(koinaSequence, @"\[[^\]]*\]|-", "");

        private static Random SeededBy(string text) => new(text.Aggregate(17, (h, c) => unchecked(h * 31 + c)));

        private static string Numbers(IEnumerable<double> values) => string.Join(",", values.Select(v => v.ToString("F4", CultureInfo.InvariantCulture)));

        private Task<string> Answer(Dictionary<string, object> request, Func<string[], string> outputs)
        {
            Interlocked.Increment(ref Calls);
            if (Unreachable)
                throw new HttpRequestException("Koina is unreachable");
            var peptides = JObject.Parse(JsonConvert.SerializeObject(request))["inputs"]![0]!["data"]!.Select(t => t.ToString()).ToArray();
            foreach (string peptide in peptides)
                Requested.Add(peptide);
            return Task.FromResult("{\"outputs\":[" + outputs(peptides) + "]}");
        }

        public sealed class FakeIntensity(CannedProsit owner) : Prosit2025Intensity40PTM
        {
            private static readonly string[] Ions = ["y3+1", "y4+1", "y5+1", "y6+1", "y7+1", "y8+1", "b3+1", "b4+1"];

            protected override Task<string> SendInferenceRequestAsync(string modelName, Dictionary<string, object> request, CancellationToken cancellationToken) =>
                owner.Answer(request, peptides =>
                {
                    int n = peptides.Length, k = Ions.Length;
                    string annotations = string.Join(",", Enumerable.Repeat(string.Join(",", Ions.Select(i => $"\"{i}\"")), n));
                    string mz = Numbers(peptides.SelectMany(p => { var r = SeededBy(p); return Ions.Select(_ => 200 + 1200 * r.NextDouble()); }));
                    string intensities = Numbers(peptides.SelectMany(p => { var r = SeededBy(p + "/intensity"); return Ions.Select(_ => 0.1 + 0.9 * r.NextDouble()); }));
                    return $"{{\"name\":\"annotation\",\"datatype\":\"BYTES\",\"shape\":[{n},{k}],\"data\":[{annotations}]}}," +
                        $"{{\"name\":\"mz\",\"datatype\":\"FP32\",\"shape\":[{n},{k}],\"data\":[{mz}]}}," +
                        $"{{\"name\":\"intensities\",\"datatype\":\"FP32\",\"shape\":[{n},{k}],\"data\":[{intensities}]}}";
                });
        }

        public sealed class FakeIrt(CannedProsit owner) : Prosit2025iRT40PTM
        {
            protected override Task<string> SendInferenceRequestAsync(string modelName, Dictionary<string, object> request, CancellationToken cancellationToken) =>
                owner.Answer(request, peptides =>
                    $"{{\"name\":\"irt\",\"datatype\":\"FP32\",\"shape\":[{peptides.Length}],\"data\":[{Numbers(peptides.Select(p => -15 + 130 * SeededBy(p + "/irt").NextDouble()))}]}}");
        }
    }

    /// <summary>
    /// Proteins of two tryptic pieces each, P00000, P00001, ...: random residues with no K, R, C or M inside (so the task's
    /// default carbamidomethyl and oxidation never apply), each piece ending in K.
    /// </summary>
    private static List<Protein> Proteins(int count, int seed = 7, int pieceLength = 12)
    {
        var random = new Random(seed);
        const string residues = "ADEFGHILNPQSTVWY";
        string Piece(int length) => new string(Enumerable.Range(0, length - 1).Select(_ => residues[random.Next(residues.Length)]).ToArray()) + "K";
        return Enumerable.Range(0, count).Select(i => new Protein(Piece(pieceLength) + Piece(pieceLength), $"P{i:D5}")).ToList();
    }

    private static string WriteProteins(string folder, IEnumerable<Protein> proteins)
    {
        string path = Path.Combine(folder, "proteins.fasta");
        File.WriteAllLines(path, proteins.SelectMany(p => new[] { $">sp|{p.Accession}|{p.Accession}_HUMAN Protein {p.Accession} OS=Homo sapiens OX=9606 GN=G{p.Accession} PE=1 SV=1", p.BaseSequence }));
        return path;
    }

    private static string WriteRunMzml(string folder, SyntheticDiaRun run)
    {
        string mzml = Path.Combine(folder, "run.mzML");
        MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(new GenericMsDataFile(AsWritten(run.Scans), new SourceFile("no nativeID format", "mzML format", null, null, null)), mzml, false);
        return mzml;
    }

    private static List<MslLibraryEntry> PredictFrom(IReadOnlyList<Protein> proteins, CannedProsit models, out PredictionClients.SpectralLibraryGeneration.MslLibraryBuildResult result,
        CommonParameters? commonParameters = null, DiaLibrarySearchTaskParameters? parameters = null) =>
        DiaLibrarySearchTask.PredictLibraryEntries(proteins, commonParameters ?? new CommonParameters(), parameters ?? new DiaLibrarySearchTaskParameters(), [], [],
            models.Intensity, models.Irt, out result);

    private static List<string> WarningsOf(Action action)
    {
        var warnings = new List<string>();
        EventHandler<StringEventArgs> warned = (_, e) => warnings.Add(e.S);
        MetaMorpheusTask.WarnHandler += warned;
        try
        {
            action();
        }
        finally
        {
            MetaMorpheusTask.WarnHandler -= warned;
        }
        return warnings;
    }

    /// <summary>
    /// End to end: given a FASTA and no library, the task predicts one from the proteins (targets and permuted decoys), says
    /// so before it starts, writes it as DiaPredictedLibrary.msl beside its tables, searches the run against it, infers
    /// protein groups from the decoys' accessions, and says in results.txt what it built and with which models.
    /// </summary>
    [Test]
    public void WithProteinsAndNoLibraryTheTaskPredictsOneAndSearchesIt()
    {
        string folder = Path.Combine(_directory, "predicted");
        Directory.CreateDirectory(folder);
        var proteins = Proteins(150);
        string fasta = WriteProteins(folder, proteins);
        // The library the task will predict, predicted here too, so the run can be planted from it
        var expected = PredictFrom(proteins, new CannedProsit(), out var expectedResult);
        var run = SyntheticDiaRun.Plant(expected, e => !e.IsDecoy && SyntheticDiaRun.Bucket(e, 4) != 0, withMs1: true);
        string mzml = WriteRunMzml(folder, run);
        string output = Path.Combine(folder, "output", "Task1");
        Directory.CreateDirectory(output);
        var models = new CannedProsit();

        var warnings = WarningsOf(() => new DiaLibrarySearchTask { PredictionModels = () => (models.Intensity, models.Irt) }
            .RunTask(output, [new DbForTask(fasta, false)], [mzml], "Task1"));

        Assert.That(warnings, Has.Some.Contains("No .msl library was given, so one is predicted from the protein databases, complete with targets and decoys"),
            "the user is told, before the long network step, that a library is being made and what it holds");
        string libraryPath = Path.Combine(output, DiaLibrarySearchTask.PredictedLibraryFileName);
        Assert.That(File.Exists(libraryPath));
        using var library = MslLibrary.Load(libraryPath);
        var entries = Enumerable.Range(0, library.PrecursorCount).Select(i => library.GetEntry(i, includeFragments: false)!).ToList();
        Assert.That(entries.Select(e => (e.FullSequence, e.ChargeState, e.IsDecoy)), Is.EquivalentTo(expected.Select(e => (e.FullSequence, e.ChargeState, e.IsDecoy))),
            "the task predicts the library predicted here");
        var accessions = proteins.Select(p => p.Accession).ToHashSet();
        Assert.That(entries.Where(e => !e.IsDecoy).SelectMany(e => e.ProteinAccession.Split('|')).Distinct(), Is.SubsetOf(accessions));
        var decoys = entries.Where(e => e.IsDecoy).ToList();
        Assert.That(decoys, Is.Not.Empty);
        Assert.That(decoys, Has.None.Matches<MslLibraryEntry>(e => e.IsEntrapment), "a permuted decoy is a decoy, not entrapment");
        Assert.That(decoys.SelectMany(e => e.ProteinAccession.Split('|')).Select(a => Regex.Match(a, @"^DECOY_(P\d{5})_f1$")),
            Has.All.Matches<Match>(m => m.Success && accessions.Contains(m.Groups[1].Value)), "each decoy names its owner's permuted decoy protein");

        string results = File.ReadAllText(Path.Combine(output, "results.txt"));
        Assert.That(results, Does.Contain($"Library: {DiaLibrarySearchTask.PredictedLibraryFileName}, predicted from the protein databases (150 proteins), " +
            $"complete with targets and decoys: {expectedResult.TargetPrecursors} target and {expectedResult.DecoyPrecursors} permuted decoy precursors"));
        Assert.That(results, Does.Contain(models.Intensity.ModelName).And.Contain(models.Irt.ModelName).And.Contain(HTTP.ModelsURL));
        int found = int.Parse(Regex.Match(results, @"All target precursors with q-value <= 0\.01: (\d+)").Groups[1].Value);
        int planted = run.Library.Count(e => !e.IsDecoy && run.PlantedSequences.Contains(e.FullSequence)
            && e.PrecursorMz >= SyntheticDiaRun.FirstWindowLowMz && e.PrecursorMz < SyntheticDiaRun.FirstWindowLowMz + SyntheticDiaRun.WindowCount * SyntheticDiaRun.WindowWidth);
        Assert.That(found, Is.InRange(planted * 7 / 10, planted * 11 / 10), "most planted targets are found, few others");
        Assert.That(File.Exists(Path.Combine(output, "AllDiaProteinGroups.tsv")), "the decoys name decoy proteins, so protein groups are inferred");
    }

    /// <summary>
    /// The library's digestion stops at the longest peptide the models take (30 for Prosit 2025), though MetaMorpheus's
    /// default maximum length is unlimited: nothing longer is asked of Koina or counted as not predicted.
    /// </summary>
    [Test]
    public void NoPeptideLongerThanTheModelsTakeIsAskedFor()
    {
        var models = new CannedProsit();

        var entries = PredictFrom(Proteins(5, pieceLength: 20), models, out var result);

        Assert.That(models.Requested.Select(CannedProsit.Residues).Max(p => p.Length), Is.EqualTo(20), "the 40-residue joins of two pieces are not digested");
        Assert.That(result.NotPredicted, Is.Zero);
        Assert.That(entries.Where(e => !e.IsDecoy).Select(e => e.BaseSequence.Length).Distinct(), Is.EqualTo(new[] { 20 }));
    }

    /// <summary>Reversed decoys, for a digestion the permuter refuses: mzLib's DECOY_ and the target's accession.</summary>
    [Test]
    public void ReversedDecoysCanBeAskedFor()
    {
        var entries = PredictFrom(Proteins(5), new CannedProsit(), out var result, parameters: new DiaLibrarySearchTaskParameters { LibraryDecoys = DiaLibraryDecoys.Reversed });

        var decoys = entries.Where(e => e.IsDecoy).ToList();
        Assert.That(decoys, Is.Not.Empty);
        Assert.That(decoys.Select(e => e.ProteinAccession), Has.All.Matches<string>(a => Regex.IsMatch(a, @"^DECOY_P\d{5}$")));
        Assert.That(result.DecoyPrecursors, Is.EqualTo(decoys.Count));
    }

    /// <summary>
    /// A contaminant protein is a target of the library too, and gets a permuted decoy like any other (the permuter itself
    /// refuses a contaminant, so the task hands it a plain copy).
    /// </summary>
    [Test]
    public void AContaminantGetsAPermutedDecoyToo()
    {
        var proteins = Proteins(3);
        proteins.Add(new Protein(Proteins(1, seed: 11)[0].BaseSequence, "Q99999", isContaminant: true));

        var entries = PredictFrom(proteins, new CannedProsit(), out _);

        Assert.That(entries.Where(e => !e.IsDecoy).Select(e => e.ProteinAccession), Does.Contain("Q99999"));
        Assert.That(entries.Where(e => e.IsDecoy).Select(e => e.ProteinAccession), Does.Contain("DECOY_Q99999_f1"));
    }

    /// <summary>A digestion the permuter cannot pair (here a semi-specific one) is refused before Koina is asked, saying what to set instead.</summary>
    [Test]
    public void ADigestionThePermuterRefusesSaysToUseReversedDecoys()
    {
        var models = new CannedProsit();
        var semi = new CommonParameters(digestionParams: new DigestionParams(searchModeType: CleavageSpecificity.Semi));

        var e = Assert.Throws<MetaMorpheusException>(() => PredictFrom(Proteins(3), models, out _, commonParameters: semi))!;

        Assert.That(e.Message, Does.Contain("LibraryDecoys = \"Reversed\""));
        Assert.That(models.Calls, Is.Zero);
    }

    /// <summary>
    /// Koina unreachable: nothing is left to search, so the task fails, saying the library needs network access to Koina and
    /// that an .msl library can be given instead, and writes no library.
    /// </summary>
    [Test]
    public void UnreachableKoinaFailsTheTaskSayingWhatItNeeds()
    {
        string folder = Path.Combine(_directory, "unreachable");
        Directory.CreateDirectory(folder);
        string fasta = WriteProteins(folder, Proteins(5));
        string mzml = WriteRunMzml(folder, SyntheticDiaRun.Build(20, _ => true));
        string output = Path.Combine(folder, "output", "Task1");
        Directory.CreateDirectory(output);
        var models = new CannedProsit { Unreachable = true };

        var e = Assert.Throws<MetaMorpheusException>(() =>
            new DiaLibrarySearchTask { PredictionModels = () => (models.Intensity, models.Irt) }.RunTask(output, [new DbForTask(fasta, false)], [mzml], "Task1"))!;

        Assert.That(e.Message, Does.Contain(HTTP.ModelsURL).And.Contain("network").And.Contain(".msl"));
        Assert.That(e.InnerException, Is.InstanceOf<HttpRequestException>());
        Assert.That(File.Exists(Path.Combine(output, DiaLibrarySearchTask.PredictedLibraryFileName)), Is.False);
    }

    /// <summary>
    /// What the model cannot predict for is refused before Koina is asked, read from the model's own limits: fragmentation
    /// other than its own (HCD or CID for Prosit 2025), charges outside its range, a digestion that is not a protease's, and
    /// an energy that is not one.
    /// </summary>
    [TestCase("ETD", "ETD")]
    [TestCase("Autodetect", "Autodetect")]
    [TestCase("rna", "protease")]
    [TestCase("nocharges", "charge")]
    [TestCase("charge7", "charge")]
    [TestCase("nce0", "NCE")]
    public void WhatTheModelCannotPredictForIsRefused(string setting, string said)
    {
        var model = new CannedProsit().Intensity;
        var commonParameters = setting switch
        {
            "ETD" => new CommonParameters(dissociationType: DissociationType.ETD),
            "Autodetect" => new CommonParameters(dissociationType: DissociationType.Autodetect),
            "rna" => new CommonParameters(digestionParams: new RnaDigestionParams()),
            _ => new CommonParameters(),
        };
        var parameters = setting switch
        {
            "nocharges" => new DiaLibrarySearchTaskParameters { LibraryPrecursorCharges = [] },
            "charge7" => new DiaLibrarySearchTaskParameters { LibraryPrecursorCharges = [2, 7] },
            "nce0" => new DiaLibrarySearchTaskParameters { LibraryNce = 0 },
            _ => new DiaLibrarySearchTaskParameters(),
        };

        Assert.That(DiaLibrarySearchTask.PredictionUnsupported(commonParameters, [], parameters, model), Does.Contain(said));
        Assert.That(DiaLibrarySearchTask.PredictionUnsupported(new CommonParameters(), [], new DiaLibrarySearchTaskParameters(), model), Is.Null, "the defaults can be predicted");
        Assert.That(DiaLibrarySearchTask.PredictionUnsupported(new CommonParameters(dissociationType: DissociationType.CID), [], new DiaLibrarySearchTaskParameters(), model), Is.Null);
        Assert.That(DiaLibrarySearchTask.PredictionUnsupported(new CommonParameters(), [DissociationType.ETD], new DiaLibrarySearchTaskParameters(), model), Does.Contain("ETD"),
            "a file's own fragmentation counts too");
    }

    /// <summary>The task refuses an unpredictable search before asking Koina, and says to give it an .msl library instead.</summary>
    [Test]
    public void TheTaskRefusesAnUnpredictableSearchBeforeAskingKoina()
    {
        string folder = Path.Combine(_directory, "etd");
        Directory.CreateDirectory(folder);
        string fasta = WriteProteins(folder, Proteins(5));
        string mzml = WriteRunMzml(folder, SyntheticDiaRun.Build(20, _ => true));
        string output = Path.Combine(folder, "output", "Task1");
        Directory.CreateDirectory(output);
        var models = new CannedProsit();
        var task = new DiaLibrarySearchTask { CommonParameters = new CommonParameters(dissociationType: DissociationType.ETD), PredictionModels = () => (models.Intensity, models.Irt) };

        var e = Assert.Throws<MetaMorpheusException>(() => task.RunTask(output, [new DbForTask(fasta, false)], [mzml], "Task1"))!;

        Assert.That(e.Message, Does.Contain("ETD").And.Contain(".msl"));
        Assert.That(models.Calls, Is.Zero);
    }

    /// <summary>With neither a library nor proteins to predict one from, the task says it needs one or the other.</summary>
    [Test]
    public void WithNeitherALibraryNorProteinsTheTaskSaysWhatItNeeds()
    {
        string folder = Path.Combine(_directory, "nothing");
        Directory.CreateDirectory(folder);
        string mzml = WriteRunMzml(folder, SyntheticDiaRun.Build(20, _ => true));
        string output = Path.Combine(folder, "output", "Task1");
        Directory.CreateDirectory(output);

        var e = Assert.Throws<MetaMorpheusException>(() => new DiaLibrarySearchTask().RunTask(output, [], [mzml], "Task1"))!;

        Assert.That(e.Message, Does.Contain(".msl").And.Contain("protein database"));
    }

    /// <summary>
    /// Given a library, the task predicts nothing, whatever protein database comes with it, and results.txt says so: only
    /// the library's spectra are searched, and the protein database adds none.
    /// </summary>
    [Test]
    public void WithALibraryNothingIsPredicted()
    {
        string folder = Path.Combine(_directory, "given");
        Directory.CreateDirectory(folder);
        var run = SyntheticDiaRun.Build(100, e => !e.IsDecoy, withMs1: true);
        string mzml = WriteRunMzml(folder, run);
        string library = run.WriteLibrary(folder);
        string fasta = WriteProteins(folder, Proteins(5));
        string output = Path.Combine(folder, "output", "Task1");
        Directory.CreateDirectory(output);
        var models = new CannedProsit { Unreachable = true };

        new DiaLibrarySearchTask { PredictionModels = () => (models.Intensity, models.Irt) }
            .RunTask(output, [new DbForTask(library, false), new DbForTask(fasta, false)], [mzml], "Task1");

        Assert.That(models.Calls, Is.Zero);
        Assert.That(File.Exists(Path.Combine(output, DiaLibrarySearchTask.PredictedLibraryFileName)), Is.False);
        string results = File.ReadAllText(Path.Combine(output, "results.txt"));
        Assert.That(results, Does.Not.Contain(DiaLibrarySearchTask.PredictedLibraryFileName));
        Assert.That(results, Does.Contain($"Library: {Path.GetFileName(library)}, as given: {run.Library.Count(e => !e.IsDecoy)} target and {run.Library.Count(e => e.IsDecoy)} decoy precursors. " +
            "Only its spectra are searched; the protein databases add none."));
    }

    /// <summary>
    /// A library with no decoys gives the search nothing to estimate an FDR from. It is refused once, before any run is
    /// calibrated (not run by run, as each calibration found it), saying what to give instead.
    /// </summary>
    [Test]
    public void ALibraryWithoutDecoysIsRefusedOnceSayingWhatToGiveInstead()
    {
        string folder = Path.Combine(_directory, "nodecoys");
        Directory.CreateDirectory(folder);
        var run = SyntheticDiaRun.Build(50, _ => true, withDecoys: false);
        string mzml = WriteRunMzml(folder, run);
        string library = run.WriteLibrary(folder);
        string output = Path.Combine(folder, "output", "Task1");
        Directory.CreateDirectory(output);
        MetaMorpheusException? e = null;

        var warnings = WarningsOf(() => e = Assert.Throws<MetaMorpheusException>(() => new DiaLibrarySearchTask().RunTask(output, [new DbForTask(library, false)], [mzml], "Task1")));

        Assert.That(e!.Message, Does.Contain($"{Path.GetFileName(library)} holds no decoys").And.Contain("protein database"));
        Assert.That(warnings, Has.None.Contains("is skipped"), "refused once, not run by run");
    }
}

/// <summary>
/// Live: a small library through the real Prosit 2025 models on Koina, with the task's default modifications
/// (carbamidomethyl on C fixed, oxidation on M variable), so MetaMorpheus's own modifications are shown to reach Prosit.
/// Skipped when Koina cannot be reached.
/// </summary>
[TestFixture]
[Category("ExternalService")]
[ExcludeFromCodeCoverage]
public class DiaPredictedLibraryLiveTests
{
    [Test]
    public void ASmallLibraryIsPredictedThroughKoina()
    {
        var (intensity, irt) = DiaLibrarySearchTask.KoinaModels();
        var fixedModifications = GlobalVariables.AllModsKnown.Where(m => m.IdWithMotif == "Carbamidomethyl on C" && m.ModificationType == "Common Fixed").Take(1).ToList();
        var variableModifications = GlobalVariables.AllModsKnown.Where(m => m.IdWithMotif == "Oxidation on M" && m.ModificationType == "Common Variable").Take(1).ToList();
        List<MslLibraryEntry> entries;
        PredictionClients.SpectralLibraryGeneration.MslLibraryBuildResult result;
        try
        {
            entries = DiaLibrarySearchTask.PredictLibraryEntries([new Protein("MSLPDKAGTCVLVEGDRPYVLNSR", "P1")], new CommonParameters(), new DiaLibrarySearchTaskParameters(),
                fixedModifications, variableModifications, intensity, irt, out result);
        }
        catch (MetaMorpheusException e) when (e.InnerException is HttpRequestException or OperationCanceledException)
        {
            Assert.Ignore($"Koina unavailable: {e.Message}");
            return;
        }

        Assert.That(fixedModifications, Has.Count.EqualTo(1));
        Assert.That(result.TargetPrecursors, Is.GreaterThan(0));
        Assert.That(entries.Where(e => !e.IsDecoy), Has.All.Matches<MslLibraryEntry>(e => e.MatchedFragmentIons.Count > 0 && double.IsFinite(e.RetentionTime)));
        Assert.That(entries.Where(e => !e.IsDecoy && e.BaseSequence == "AGTCVLVEGDR").Select(e => e.FullSequence), Has.All.Contain("Carbamidomethyl").And.Not.Empty,
            "MetaMorpheus's carbamidomethyl reaches Prosit");
        Assert.That(entries.Where(e => !e.IsDecoy && e.BaseSequence.StartsWith("MSLPDK")).Select(e => e.FullSequence), Has.Some.Contain("Oxidation"),
            "and its variable oxidation");
        Assert.That(entries.Where(e => e.IsDecoy).Select(e => e.ProteinAccession), Has.All.StartWith("DECOY_P1_f1"));
    }
}
