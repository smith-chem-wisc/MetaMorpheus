using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using EngineLayer;
using MassSpectrometry;
using Omics.SpectralMatch.MslSpectralLibrary;
using PredictionClients.Koina.Client;
using PredictionClients.Koina.SupportedModels.RetentionTimeModels;
using PredictionClients.SpectralLibraryGeneration;
using Readers.SpectralLibrary;

namespace TaskLayer
{
    /// <summary>
    /// Updates an <c>.msl</c> spectral library as <c>.msl</c>, with its retention times on the iRT scale.
    ///
    /// The <c>.msp</c> update stores each spectrum's scan time in minutes, which means nothing outside the
    /// run it came from. Here every step is mzLib's (1.0.600):
    /// <list type="number">
    ///   <item>Spectra: <see cref="MslLibrary.UpdateAndSave"/> with the search's best spectrum per precursor and
    ///   no time, so it decides only which spectrum to keep and keeps each time and its status.</item>
    ///   <item>Times: <see cref="MslLibrary.UpdateRetentionTimesAndSave"/> with every confident PSM, grouped by
    ///   raw file. Each file is calibrated onto the iRT scale against Prosit 2025's predictions, and its times
    ///   then fill, confirm or challenge the stored ones.</item>
    ///   <item>Predicted times: <see cref="MslPredictedRetentionTimes.Assign"/> in
    ///   <see cref="MslPredictedIrtMode.FillMissing"/> mode, for entries that have never been measured.</item>
    /// </list>
    /// Steps 2 and 3 ask Koina, a third-party web service. Each either finishes or writes nothing, so when
    /// Koina fails the library written is the one the last finished step left, the search keeps all its other
    /// results, and results.txt says what happened. The search's times did not go in, so the same files can
    /// update the library again. The library keeps no record of which files it has seen, so a file that did
    /// update it must not be fed to it twice.
    /// </summary>
    public partial class PostSearchAnalysisTask
    {
        /// <summary>Sequences sent to the iRT predictor in one call; Stop is checked between calls.</summary>
        internal const int MslPredictionChunkSize = 5_000;

        /// <summary>
        /// Stands in for the Koina iRT predictor when set, so tests can run the <c>.msl</c> update without the
        /// network: iRT for each full sequence given, in order, NaN where there is none. Null in production.
        /// </summary>
        internal Func<IReadOnlyList<string>, IReadOnlyList<double>> IrtPredictor { get; set; }

        /// <summary>
        /// The spectral library to update as <c>.msl</c>: the one library given, when it is an <c>.msl</c>.
        /// Null otherwise, and the update writes <c>.msp</c> as before.
        /// </summary>
        internal string SingleMslLibraryPath()
        {
            var libraries = (Parameters.DatabaseFilenameList ?? new List<EngineLayer.DatabaseLoading.DbForTask>())
                .Where(p => p.IsSpectralLibrary).Select(p => p.FilePath).ToList();
            if (libraries.Count == 1 && MslFileTypeHandler.IsMslFile(libraries[0]))
                return libraries[0];
            if (libraries.Any(MslFileTypeHandler.IsMslFile))
                Warn($"An .msl library is updated as .msl only when it is the one spectral library given; these {libraries.Count} libraries are written together as .msp.");
            return null;
        }

        /// <summary>
        /// Updates the <c>.msl</c> library at <paramref name="libraryPath"/> from <paramref name="psms"/> (the
        /// search's confident target PSMs) and writes it to the output folder. Returns its path, or null when
        /// Stop was pressed.
        /// </summary>
        internal string UpdateMslLibrary(string libraryPath, IReadOnlyList<SpectralMatch> psms)
        {
            Status("Updating the .msl spectral library...", Parameters.SearchTaskId);
            var report = new MslUpdateReport();
            string stamp = DateTime.Now.ToString("yyyy-MM-dd-HH-mm-ss", CultureInfo.InvariantCulture);
            string outputPath = Path.Combine(Parameters.OutputFolder, $"updateSpectralLibrary_{stamp}.msl");
            string spectraPath = Path.Combine(Parameters.OutputFolder, $"updateSpectralLibrary_{stamp}_spectra.msl.partial");
            string timesPath = Path.Combine(Parameters.OutputFolder, $"updateSpectralLibrary_{stamp}_times.msl.partial");
            try
            {
                // 1. Spectra. The best match per precursor, as the .msp update chooses it, without its time.
                Dictionary<(string, int), SpectralMatch> best = psms
                    .Where(p => p.BioPolymerWithSetModsMonoisotopicMass.HasValue)
                    .GroupBy(p => (p.FullSequence, p.ScanPrecursorCharge))
                    .ToDictionary(g => g.Key, g => g.MaxBy(p => p.Score));
                var spectra = best.Values.Select(p =>
                {
                    var spectrum = p.ToLibrarySpectrum();
                    spectrum.RetentionTime = null;
                    return spectrum;
                }).ToList();

                Dictionary<(string, int), MslLibraryEntry> before;
                using (var library = MslLibrary.Load(libraryPath))
                {
                    before = EntriesFor(library.GetAllEntries(includeDecoys: false), best.Keys);
                    report.Spectra = library.UpdateAndSave(spectra, spectraPath);
                }
                string current = spectraPath;

                // 2. Times, and 3. predicted times for entries never measured: both from Koina's predictions.
                Func<IReadOnlyList<string>, IReadOnlyList<double>> predictIrt = null;
                if (GlobalVariables.AnalyteType != AnalyteType.Peptide)
                    report.TimesNotUpdated = $"Prosit predicts peptides, and this is a {GlobalVariables.AnalyteType} search";
                else
                {
                    predictIrt = MslIrtPredictor(report);
                    var runs = psms
                        .GroupBy(p => p.FullFilePath)
                        .Select(g => new MslRunObservations(Path.GetFileNameWithoutExtension(g.Key),
                            g.Select(p => new MslRunPsm(p.FullSequence, p.BaseSequence, p.ScanPrecursorCharge, p.ScanRetentionTime)).ToList()))
                        .ToList();
                    try
                    {
                        using var library = MslLibrary.Load(current);
                        var predictLongPeptideIrt = LongPeptideIrtPredictor(KnownIrt(library.GetAllEntries(includeDecoys: false)), report);
                        Status("Updating library retention times...", Parameters.SearchTaskId);
                        report.Times = library.UpdateRetentionTimesAndSave(runs, predictIrt, timesPath, predictLongPeptideIrt);
                        current = timesPath;
                    }
                    catch (OperationCanceledException) when (GlobalVariables.StopLoops)
                    {
                        return null;
                    }
                    catch (Exception e)
                    {
                        report.TimesNotUpdated = $"{e.GetType().Name}: {e.Message}";
                        Warn($"The library's retention times were not updated: {report.TimesNotUpdated}. The spectra were; " +
                            "this search's times did not go in, so these files can update the library again.");
                    }
                }
                if (GlobalVariables.StopLoops)
                    return null;

                List<MslLibraryEntry> entries;
                using (var library = MslLibrary.Load(current))
                    entries = library.GetAllEntries().ToList();

                // A predictor that just failed is not asked again: an outage would cost its timeout twice.
                if (predictIrt is not null && report.PredictorFailed)
                    report.PredictionsNotGiven = "the predictor was unavailable";
                else if (predictIrt is not null)
                {
                    try
                    {
                        // Chronologer is mapped with the times step 2 just gave, not the ones it started from
                        entries = MslPredictedRetentionTimes.Assign(entries, predictIrt, MslPredictedIrtMode.FillMissing,
                            out var predicted, LongPeptideIrtPredictor(KnownIrt(entries.Where(e => !e.IsDecoy)), report));
                        report.Predicted = predicted;
                    }
                    catch (OperationCanceledException) when (GlobalVariables.StopLoops)
                    {
                        return null;
                    }
                    catch (Exception e)
                    {
                        report.PredictionsNotGiven = $"{e.GetType().Name}: {e.Message}";
                        Warn($"Entries without a retention time were not given a predicted one: {report.PredictionsNotGiven}. " +
                            "The next update gives them one.");
                    }
                    if (GlobalVariables.StopLoops)
                        return null;
                }

                report.DissociationTypesSet = SetDissociationTypes(entries, before, best);
                MslLibrary.Save(outputPath, entries);
                FinishedWritingFile(outputPath, new List<string> { Parameters.SearchTaskId });
                Parameters.SearchTaskResults.AddTaskSummaryText(report.Describe(outputPath));
                return outputPath;
            }
            finally
            {
                File.Delete(spectraPath);
                File.Delete(timesPath);
            }
        }

        /// <summary>
        /// The iRT predictor for both steps: Koina's Prosit 2025 (or <see cref="IrtPredictor"/>). Each full sequence
        /// is asked once per update, and at most <see cref="MslPredictionChunkSize"/> are sent in one call, with Stop
        /// checked between calls. An outage propagates as the predictor throws it.
        /// </summary>
        private Func<IReadOnlyList<string>, IReadOnlyList<double>> MslIrtPredictor(MslUpdateReport report)
        {
            Func<IReadOnlyList<string>, IReadOnlyList<double>> predict = IrtPredictor;
            if (predict is null)
            {
                // As many requests in flight as the spectral-angle step allows (MaxConcurrentPredictionRequests).
                var model = new Prosit2025iRT40PTM(maxNumberOfBatchesPerRequest: MaxConcurrentPredictionRequests);
                report.Model = $"{model.ModelName} at {HTTP.ModelsURL}";
                predict = MslPredictedRetentionTimes.FromKoina(model);
            }

            var known = new Dictionary<string, double>(StringComparer.Ordinal);
            return sequences =>
            {
                var needed = sequences.Where(s => !known.ContainsKey(s)).Distinct(StringComparer.Ordinal).ToList();
                foreach (string[] chunk in needed.Chunk(MslPredictionChunkSize))
                {
                    if (GlobalVariables.StopLoops)
                        throw new OperationCanceledException("Stop was pressed.");
                    Status($"Predicting retention times for {needed.Count} peptides...", Parameters.SearchTaskId);
                    IReadOnlyList<double> values;
                    try
                    {
                        values = predict(chunk);
                    }
                    catch
                    {
                        report.PredictorFailed = true;
                        throw;
                    }
                    if (values is null || values.Count != chunk.Length)
                        throw new InvalidOperationException(
                            $"The iRT predictor returned {values?.Count.ToString() ?? "null"} values for {chunk.Length} sequences.");
                    for (int i = 0; i < chunk.Length; i++)
                    {
                        known[chunk[i]] = double.IsFinite(values[i]) ? values[i] : double.NaN;
                        report.Requested++;
                        if (!double.IsFinite(values[i]))
                            report.WithoutPrediction++;
                    }
                }
                return sequences.Select(s => known[s]).ToArray();
            };
        }

        /// <summary>
        /// Chronologer, run locally and mapped onto the iRT of <paramref name="knownIrt"/>, for unmodified peptides of
        /// 31-50 residues (beyond Prosit's 30). It is fitted only when such a peptide needs it. When it cannot be
        /// (fewer than 20 known peptides) or fails, those peptides have no prediction, and the report says why.
        /// </summary>
        private static Func<IReadOnlyList<string>, IReadOnlyList<double>> LongPeptideIrtPredictor(
            Dictionary<string, double> knownIrt, MslUpdateReport report)
        {
            Func<IReadOnlyList<string>, IReadOnlyList<double>> fitted = null;
            bool failed = false;
            return sequences =>
            {
                try
                {
                    if (fitted is null && !failed)
                    {
                        fitted = MslPredictedRetentionTimes.ChronologerOnIrtScale(knownIrt);
                        report.LongPeptidePredictorUsed = true;
                    }
                    if (fitted is not null)
                        return fitted(sequences);
                }
                catch (Exception e)
                {
                    failed = true;
                    report.LongPeptidePredictorNotUsed = $"{e.GetType().Name}: {e.Message}";
                }
                return Enumerable.Repeat(double.NaN, sequences.Count).ToArray();
            };
        }

        /// <summary>Each target's iRT, by full sequence: what Chronologer is mapped onto the library's scale with.</summary>
        private static Dictionary<string, double> KnownIrt(IEnumerable<MslLibraryEntry> entries)
        {
            var known = new Dictionary<string, double>(StringComparer.Ordinal);
            foreach (var e in entries)
                if (!e.RtIsCalibrated && double.IsFinite(e.RetentionTime))
                    known.TryAdd(e.FullSequence, e.RetentionTime);
            return known;
        }

        /// <summary>The library's target entries for <paramref name="keys"/>, before the update.</summary>
        private static Dictionary<(string, int), MslLibraryEntry> EntriesFor(IEnumerable<MslLibraryEntry> entries,
            IEnumerable<(string, int)> keys)
        {
            var wanted = keys.ToHashSet();
            var found = new Dictionary<(string, int), MslLibraryEntry>();
            foreach (var e in entries)
                if (wanted.Contains((e.FullSequence, e.ChargeState)))
                    found.TryAdd((e.FullSequence, e.ChargeState), e);
            return found;
        }

        /// <summary>
        /// Gives each spectrum this search added or replaced the dissociation type of its raw file, which
        /// <see cref="MslLibrary.UpdateAndSave"/> leaves unknown; a spectrum the library kept keeps its own. A spectrum
        /// is this search's when the library had none for its precursor, or had a different one. Returns how many.
        /// </summary>
        private int SetDissociationTypes(List<MslLibraryEntry> entries, Dictionary<(string, int), MslLibraryEntry> before,
            Dictionary<(string, int), SpectralMatch> best)
        {
            int set = 0;
            foreach (var entry in entries)
            {
                var key = (entry.FullSequence, entry.ChargeState);
                if (entry.IsDecoy || !best.TryGetValue(key, out var psm))
                    continue;
                if (before.TryGetValue(key, out var original) && SameFragments(entry, original))
                    continue;
                entry.DissociationType = DissociationTypeOf(psm.FullFilePath);
                set++;
            }
            return set;
        }

        private static bool SameFragments(MslLibraryEntry a, MslLibraryEntry b)
        {
            if (a.MatchedFragmentIons.Count != b.MatchedFragmentIons.Count)
                return false;
            for (int i = 0; i < a.MatchedFragmentIons.Count; i++)
                if (a.MatchedFragmentIons[i].Mz != b.MatchedFragmentIons[i].Mz
                    || a.MatchedFragmentIons[i].Intensity != b.MatchedFragmentIons[i].Intensity)
                    return false;
            return true;
        }

        /// <summary>
        /// The dissociation type a raw file was searched with. Unknown for Autodetect, which is decided per scan and
        /// not recorded on a PSM, and for Custom, which names no type.
        /// </summary>
        private DissociationType DissociationTypeOf(string rawFilePath)
        {
            var parameters = FileSpecificParameters?.FirstOrDefault(f => f.FileName == rawFilePath).Parameters ?? CommonParameters;
            return parameters.DissociationType is DissociationType.Autodetect or DissociationType.Custom
                ? DissociationType.Unknown
                : parameters.DissociationType;
        }

        /// <summary>What the <c>.msl</c> update did, for results.txt.</summary>
        private sealed class MslUpdateReport
        {
            public string Model = "a substitute predictor";
            public MslUpdateResult Spectra;
            public MslRetentionTimeUpdateResult Times;
            public string TimesNotUpdated;
            public MslPredictedIrtResult Predicted;
            public string PredictionsNotGiven;
            public string LongPeptidePredictorNotUsed;
            public int Requested, WithoutPrediction, DissociationTypesSet;
            public bool PredictorFailed, LongPeptidePredictorUsed;

            public string Describe(string outputPath)
            {
                var text = new StringBuilder();
                text.Append($"Spectral library updated as .msl: {Path.GetFileName(outputPath)}. Spectra: {Spectra.AddedCount} added, " +
                    $"{Spectra.ReplacedCount} replaced, {Spectra.RetainedCount} kept; {DissociationTypesSet} given their file's dissociation type.");

                if (Times is not null)
                {
                    var notCalibrated = Times.Runs.Where(r => !r.Calibrated).ToList();
                    text.Append($" Retention times (iRT from {Model}): {Times.Runs.Count - notCalibrated.Count} of {Times.Runs.Count} files calibrated");
                    if (notCalibrated.Count > 0)
                        text.Append(" (not: " + string.Join(", ", notCalibrated.Select(r => $"{r.RunName}, {r.Gate}")) + ")");
                    text.Append($"; {Times.Filled + Times.FilledChallenged} filled, {Times.Confirmed} confirmed, {Times.Challenged} challenged, " +
                        $"{Times.ReplacedSecondStrike + Times.ReplacedAndStruckAgain + Times.ReplacedTiebreak} replaced" +
                        (Times.OnlyConfirmations ? "; this search disagreed with the library too often, so it only confirmed times" : "") + ".");
                }
                else
                    text.Append($" Retention times not updated: {TimesNotUpdated}.");

                if (Predicted is not null)
                    text.Append($" Predicted times given to {Predicted.FromIrtPredictor + Predicted.FromLongPeptidePredictor} entries never measured; " +
                        $"{Predicted.NoPrediction} have none.");
                else if (PredictionsNotGiven is not null)
                    text.Append($" No predicted times given: {PredictionsNotGiven}.");

                if (Requested > 0)
                    text.Append($" {Requested} peptides sent to the predictor, {WithoutPrediction} without a prediction.");
                if (LongPeptidePredictorNotUsed is not null && !LongPeptidePredictorUsed)
                    text.Append($" Chronologer (peptides of 31-50 residues) not used: {LongPeptidePredictorNotUsed}.");
                return text.ToString();
            }
        }
    }
}
