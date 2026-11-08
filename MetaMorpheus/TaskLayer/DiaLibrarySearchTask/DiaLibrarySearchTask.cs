#nullable enable
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using EngineLayer;
using EngineLayer.DatabaseLoading;
using EngineLayer.DiaLibrarySearch;
using MassSpectrometry;
using Omics.Modifications;
using Proteomics;
using Readers;
using Readers.SpectralLibrary;
using UsefulProteomicsDatabases;

namespace TaskLayer;

/// <summary>
/// Searches DIA runs against one spectral library (.msl), given as a database. Each run calibrates its own retention time
/// onto the library's iRT and its MS1 offset, then is searched and scored; results.txt reports the target precursors at
/// the q-value threshold, per run and in total, and the target peptides. A protein database is not needed; given one
/// (FASTA or XML, beside the library), the task also infers protein groups, and labels C the precursors and peptides
/// whose proteins are all in a contaminant database.
/// </summary>
public class DiaLibrarySearchTask : MetaMorpheusTask
{
    public DiaLibrarySearchTaskParameters Parameters { get; set; }

    public DiaLibrarySearchTask(DiaLibrarySearchTaskParameters parameters) : base(MyTask.DiaLibrarySearch)
    {
        CommonParameters = new CommonParameters();
        Parameters = parameters;
    }

    /// <summary>Used when reading TOML and when `CMD -g` writes the default task.</summary>
    public DiaLibrarySearchTask() : this(new DiaLibrarySearchTaskParameters())
    {
    }

    /// <summary>
    /// A run's scans with every peak: never through <see cref="MyFileManager"/>'s filter, so CommonParameters' MS2
    /// trimming (on by default) can never reach DIA scans (user, 2026-09-30).
    /// </summary>
    public static MsDataScan[] LoadScans(string path) =>
        MsDataFileReader.GetDataFile(path).LoadAllStaticData().GetAllScansList().ToArray();

    /// <summary>
    /// A match as a row of AllDiaPrecursors.tsv (mzLib's <see cref="DiaPrecursorFromTsv"/>, the schema agreed with dataRepo).
    /// The run key is the file name as given; the quantity is empty when not measured.
    /// </summary>
    internal static DiaPrecursorFromTsv PrecursorRow(DiaPrecursorMatch match, string fileName, MslLibrary library, double globalQValue,
        IReadOnlySet<string>? contaminantAccessions = null)
    {
        var entry = library.GetEntry(match.PrecursorIndex);
        return new DiaPrecursorFromTsv
        {
            FileName = fileName,
            FullSequence = match.FullSequence,
            BaseSequence = entry?.BaseSequence ?? "",
            PrecursorCharge = match.Charge,
            PrecursorMz = match.PrecursorMz,
            ProteinAccession = entry?.ProteinAccession ?? "",
            Label = Label(match.IsDecoy, entry?.ProteinAccession, contaminantAccessions),
            Score = match.Score,
            QValuePrecursorRun = match.QValue,
            QValuePrecursorGlobal = globalQValue,
            LibraryIrt = match.LibraryIrt.Value,
            ApexIrt = match.ApexIrt.Value,
            ApexRtMin = match.ApexRt.Value,
            ApexScanNumber = match.ApexScanNumber,
            PrecursorQuantity = double.IsNaN(match.Quantity) ? null : match.Quantity,
        };
    }

    /// <summary>
    /// A peptide as a row of AllDiaPeptides.tsv (mzLib's <see cref="DiaPeptideFromTsv"/>): its peptide q-value, and the run
    /// and charge of its best precursor, so the row joins the precursor table.
    /// </summary>
    internal static DiaPeptideFromTsv PeptideRow(DiaPeptideMatch peptide, DiaPrecursorMatch precursor, string fileName, MslLibrary library,
        IReadOnlySet<string>? contaminantAccessions = null)
    {
        var entry = library.GetEntry(precursor.PrecursorIndex);
        return new DiaPeptideFromTsv
        {
            FullSequence = peptide.FullSequence,
            BaseSequence = entry?.BaseSequence ?? "",
            ProteinAccession = entry?.ProteinAccession ?? "",
            Label = Label(peptide.IsDecoy, entry?.ProteinAccession, contaminantAccessions),
            Score = peptide.Score,
            QValuePeptideGlobal = peptide.QValue,
            FileName = fileName,
            PrecursorCharge = precursor.Charge,
        };
    }

    /// <summary>
    /// The T/D/C label of a library entry: D for a decoy; C for a target every one of whose proteins is in a contaminant
    /// database the search was given; otherwise T. With no contaminant database, nothing is C.
    /// </summary>
    internal static string Label(bool isDecoy, string? proteinAccession, IReadOnlySet<string>? contaminantAccessions)
    {
        if (isDecoy)
            return "D";
        var accessions = (proteinAccession ?? "").Split('|', StringSplitOptions.RemoveEmptyEntries);
        return contaminantAccessions is { Count: > 0 } && accessions.Length > 0 && accessions.All(contaminantAccessions.Contains) ? "C" : "T";
    }

    /// <summary>The target precursors at or below the q-value threshold, which the headline counts and the table lists; never a decoy.</summary>
    public static IEnumerable<DiaPrecursorMatch> PassingTargets(IEnumerable<DiaPrecursorMatch> matches, double qValueThreshold) =>
        matches.Where(m => !m.IsDecoy && m.QValue <= qValueThreshold);

    /// <summary>
    /// A run's passing targets that are also reported (PLAN D4): those whose global q-value passes too.
    /// </summary>
    public static List<DiaPrecursorMatch> Reported(IEnumerable<DiaPrecursorMatch> runPassing, IReadOnlyDictionary<int, double> globalQValues,
        double qValueThreshold) =>
        runPassing.Where(m => globalQValues[m.PrecursorIndex] <= qValueThreshold).ToList();

    /// <summary>How many <see cref="PassingTargets"/> there are.</summary>
    public static int CountPassingTargets(IEnumerable<DiaPrecursorMatch> matches, double qValueThreshold) =>
        PassingTargets(matches, qValueThreshold).Count();

    protected override MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId,
        FileSpecificParameters[] fileSettingsList)
    {
        MyTaskResults = new MyTaskResults(this);
        var libraries = dbFilenameList.Where(db => db.IsSpectralLibrary).ToList();
        if (libraries.Count != 1 || !string.Equals(Path.GetExtension(libraries[0].FilePath), ".msl", StringComparison.OrdinalIgnoreCase))
            throw new MetaMorpheusException($"A DIA library search needs exactly one .msl library; {libraries.Count} spectral libraries were given.");

        string libraryPath = libraries[0].FilePath;
        Status("Loading library...", taskId);
        // Above 2 GB the whole-file load cannot hold the library in one array (MSL-Q6); index-only reads entries on demand
        using var library = new FileInfo(libraryPath).Length > int.MaxValue ? MslLibrary.LoadIndexOnly(libraryPath) : MslLibrary.Load(libraryPath);
        var parameters = Parameters.ToEngineParameters();
        string threshold = Parameters.QValueThreshold.ToString(CultureInfo.InvariantCulture);

        // A protein database, if given, is loaded before any run is searched, so a bad one fails at once. Its targets are
        // for protein inference; the library's decoys name their own decoy proteins, so none are generated.
        var proteinDatabases = dbFilenameList.Where(db => !db.IsSpectralLibrary).ToList();
        List<Protein>? proteins = null;
        List<Modification> fixedModifications = [], variableModifications = [];
        HashSet<string>? contaminantAccessions = null;
        if (proteinDatabases.Count > 0)
        {
            Status("Loading proteins...", taskId);
            LoadModifications(taskId, out variableModifications, out fixedModifications, out var localizableModificationTypes);
            var loaded = (DatabaseLoadingEngineResults)new DatabaseLoadingEngine(CommonParameters, FileSpecificParameters, [taskId], proteinDatabases,
                taskId, DecoyType.None, true, localizableModificationTypes).Run();
            proteins = loaded.BioPolymers.OfType<Protein>().ToList();
            contaminantAccessions = proteins.Where(p => p.IsContaminant).Select(p => p.Accession).ToHashSet();
        }

        // Each run's targets passing its own q-value, kept until every run is in and the global q-values are known
        var passing = new List<(string FileName, List<DiaPrecursorMatch> Matches)>();
        var bestScores = new DiaPrecursorFdr.BestScores();
        var bestPeptides = new DiaPeptideFdr.BestPrecursors();
        string? firstRun = null;
        var calibrationLines = new List<string> { "File Name\tAnchorRtMin\tLibraryIrt\tCalibratedIrt\tResidualIrt" };
        foreach (string rawFile in currentRawFileList)
        {
            if (GlobalVariables.StopLoops)
                break;
            string fileName = Path.GetFileName(rawFile);
            var ids = new List<string> { taskId, "Individual Spectra Files", rawFile };
            StartingDataFile(rawFile, ids);

            Status("Loading spectra file...", ids);
            var scans = LoadScans(rawFile);
            Status("Calibrating retention time...", ids);
            var calibration = DiaIrtSelfCalibration.Calibrate(scans, library, parameters, CommonParameters);
            foreach (var (rt, irt) in calibration.Anchors)
            {
                double calibrated = calibration.Model.ToIrt(rt).Value;
                calibrationLines.Add(string.Join('\t', fileName, rt.Value.ToString("R", CultureInfo.InvariantCulture), irt.Value.ToString("R", CultureInfo.InvariantCulture),
                    calibrated.ToString("R", CultureInfo.InvariantCulture), (irt.Value - calibrated).ToString("R", CultureInfo.InvariantCulture)));
            }
            Status("Searching...", ids);
            var results = (DiaLibrarySearchResults)new DiaLibrarySearchEngine(scans, library, calibration.Model, calibration.ApplyTo(parameters),
                CommonParameters, FileSpecificParameters, ids).Run();

            bestScores.Add(results.Matches);
            bestPeptides.Add(fileName, results.Matches);
            firstRun ??= rawFile;
            passing.Add((fileName, PassingTargets(results.Matches, Parameters.QValueThreshold).ToList()));
            FinishedDataFile(rawFile, ids);
        }

        // PLAN D4: a precursor is reported in a run when both its run q-value and its global q-value (best score across the
        // runs, then target-decoy competition) pass. The headline is the distinct precursors reported; with one run the
        // global q-value is the run's, so every run-passing target is reported.
        var global = bestScores.GlobalQValues();
        var rows = new List<DiaPrecursorFromTsv>();
        var reported = new HashSet<int>();
        var perRun = new List<string>();
        foreach (var (fileName, matches) in passing)
        {
            var kept = Reported(matches, global, Parameters.QValueThreshold);
            rows.AddRange(kept.Select(m => PrecursorRow(m, fileName, library, global[m.PrecursorIndex], contaminantAccessions)));
            reported.UnionWith(kept.Select(m => m.PrecursorIndex));
            perRun.Add($"{fileName}: {kept.Count} target precursors with q-value <= {threshold}");
        }
        int total = reported.Count;

        bool contaminantsAssessed = dbFilenameList.Any(db => db.IsContaminant);
        string tablePath = Path.Combine(OutputFolder, "AllDiaPrecursors.tsv");
        new DiaPrecursorFile(tablePath, rows, contaminantsAssessed).WriteResults(tablePath);
        FinishedWritingFile(tablePath, new List<string> { taskId });

        // PLAN D4, "peptide, global": each peptide's best precursor over the runs, then a fresh competition among peptides
        var peptides = bestPeptides.Assign();
        var peptideRows = peptides.Where(p => !p.Peptide.IsDecoy && p.Peptide.QValue <= Parameters.QValueThreshold)
            .Select(p => PeptideRow(p.Peptide, p.Precursor, p.Run, library, contaminantAccessions)).ToList();
        string peptidePath = Path.Combine(OutputFolder, "AllDiaPeptides.tsv");
        new DiaPeptideFile(peptidePath, peptideRows, contaminantsAssessed).WriteResults(peptidePath);
        FinishedWritingFile(peptidePath, new List<string> { taskId });

        int? proteinGroupCount = proteins is null || firstRun is null ? null
            : WriteProteinGroups(OutputFolder, taskId, peptides, proteins, library, firstRun, fixedModifications, variableModifications);

        // Each run's calibration anchors (D12), so a bad retention-time fit can be seen
        string calibrationPath = Path.Combine(OutputFolder, "DiaIrtCalibration.tsv");
        File.WriteAllLines(calibrationPath, calibrationLines);
        FinishedWritingFile(calibrationPath, new List<string> { taskId });

        MyTaskResults.AddTaskSummaryText($"All target precursors with q-value <= {threshold}: {total}");
        MyTaskResults.AddTaskSummaryText($"All target peptides with q-value <= {threshold}: {peptideRows.Count}");
        if (proteinGroupCount is int groups)
            MyTaskResults.AddTaskSummaryText($"All target protein groups with q-value <= {threshold}: {groups}");
        foreach (string line in perRun)
            MyTaskResults.AddTaskSummaryText(line);
        return MyTaskResults;
    }

    /// <summary>
    /// AllDiaProteinGroups.tsv (PLAN D4, "protein group, global"): protein groups inferred from the peptides of every run by
    /// MetaMorpheus's parsimony and picked protein FDR (<see cref="DiaProteinInference"/>), in MetaMorpheus's protein-group
    /// format, so mzLib's protein-group reader reads it. The target groups at the threshold are written. Decoy peptides are
    /// placed on the decoy proteins their library entries name. Returns the groups written, or null, writing nothing, when
    /// no decoy peptide names a decoy protein: the groups' q-values would rest on no decoy evidence.
    /// </summary>
    /// <param name="runPath">
    /// The path given to the stand-in spectral matches. The groups pool every run and nothing is read per run, so any run
    /// searched will do; the first is used.
    /// </param>
    private int? WriteProteinGroups(string outputFolder, string taskId, List<(DiaPeptideMatch Peptide, DiaPrecursorMatch Precursor, string Run)> peptides,
        List<Protein> proteins, MslLibrary library, string runPath, List<Modification> fixedModifications, List<Modification> variableModifications)
    {
        Status("Inferring protein groups...", taskId);
        var decoyAccessions = peptides.Where(p => p.Peptide.IsDecoy).ToDictionary(p => p.Peptide.FullSequence,
            p => (IReadOnlyList<string>)(library.GetEntry(p.Precursor.PrecursorIndex)?.ProteinAccession ?? "").Split('|', StringSplitOptions.RemoveEmptyEntries));
        var inference = DiaProteinInference.Infer(peptides.Select(p => p.Peptide).ToList(), proteins, CommonParameters, runPath,
            fixedModifications, variableModifications, decoyAccessionsOf: sequence => decoyAccessions.GetValueOrDefault(sequence) is { Count: > 0 } a ? a : null);
        if (!inference.SpectralMatches.Any(m => m.IsDecoy))
        {
            Warn($"No decoy peptide names a decoy protein in the library (of {peptides.Count(p => p.Peptide.IsDecoy)} decoy peptides scored), so protein q-values would rest on no decoy evidence; " +
                "AllDiaProteinGroups.tsv is not written. Build the library with decoy protein accessions (mzLib's library builder records them).");
            return null;
        }

        var groups = inference.ProteinGroups.Where(g => !g.IsDecoy && g.QValue <= Parameters.QValueThreshold).ToList();
        string path = Path.Combine(outputFolder, "AllDiaProteinGroups.tsv");
        using (var writer = new StreamWriter(path))
        {
            // With no quantification the header does not depend on the group, so an empty table still has one
            var headerGroup = inference.ProteinGroups.FirstOrDefault() ?? new ProteinGroup([], [], []);
            writer.WriteLine(headerGroup.GetTabSeparatedHeader());
            foreach (var group in groups)
            {
                group.GetIdentifiedPeptidesOutput(null);
                writer.WriteLine(group);
            }
        }
        FinishedWritingFile(path, new List<string> { taskId });
        return groups.Count;
    }
}
