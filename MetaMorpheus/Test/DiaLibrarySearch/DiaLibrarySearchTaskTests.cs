#nullable enable
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using EngineLayer;
using EngineLayer.DatabaseLoading;
using MassSpectrometry;
using Nett;
using NUnit.Framework;
using Readers;
using TaskLayer;

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
    private (string Mzml, string Library, SyntheticDiaRun Run) WriteInputs(string name, double noisePeaksPerScan = 60,
        Func<Omics.SpectralMatch.MslSpectralLibrary.MslLibraryEntry, double>? abundance = null, bool accessions = false)
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
        var task = new DiaLibrarySearchTask(new DiaLibrarySearchTaskParameters { FragmentTolerancePpm = 15, Ms1TolerancePpm = 4, ClassifierSeed = 3, QValueThreshold = 0.05 });
        string path = Path.Combine(_directory, "DiaLibrarySearchTask.toml");
        Toml.WriteFile(task, path, MetaMorpheusTask.tomlConfig);

        var read = Toml.ReadFile<DiaLibrarySearchTask>(path, MetaMorpheusTask.tomlConfig);

        Assert.That(read.TaskType, Is.EqualTo(MyTask.DiaLibrarySearch));
        Assert.That(read.Parameters.FragmentTolerancePpm, Is.EqualTo(15));
        Assert.That(read.Parameters.Ms1TolerancePpm, Is.EqualTo(4));
        Assert.That(read.Parameters.ClassifierSeed, Is.EqualTo(3));
        Assert.That(read.Parameters.QValueThreshold, Is.EqualTo(0.05));
        Assert.That(new DiaLibrarySearchTask().Parameters.QValueThreshold, Is.EqualTo(0.01));
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
}
