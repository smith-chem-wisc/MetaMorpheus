using EngineLayer;
using EngineLayer.CircularSearch;
using EngineLayer.DatabaseLoading;
using EngineLayer.FdrAnalysis;
using MassSpectrometry;
using MzLibUtil;
using Omics;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using UsefulProteomicsDatabases;

namespace TaskLayer
{
    /// <summary>
    /// A MetaMorpheus task that searches MS/MS spectra against a database of
    /// head-to-tail cyclic proteins represented as <see cref="CircularProtein"/>
    /// entries.
    ///
    /// WORKFLOW
    /// --------
    /// 1. Load spectra files.
    /// 2. Load protein databases — any <see cref="Protein"/> entry is wrapped in
    ///    a <see cref="CircularProtein"/> so that its sequence is canonicalized
    ///    and the circular digestion logic is applied.
    /// 3. For each spectra file, run <see cref="CircularSearchEngine"/> which
    ///    digests each <see cref="CircularProtein"/>, scores precursor candidates,
    ///    and populates a PSM array.
    /// 4. Perform FDR calculation and write output files.
    ///
    /// OUTPUT FILES
    /// ------------
    /// • CircularPeptideSpectralMatches.psmtsv
    /// • CircularPeptides.psmtsv  (unique peptide-level results)
    /// • (Optional) .mzID
    /// </summary>
    public class CircularSearchTask : MetaMorpheusTask
    {
        public CircularSearchTask() : base(MyTask.CircularSearch)
        {
            CommonParameters = new CommonParameters();
            CircularSearchParameters = new CircularSearchParameters();
        }

        public CircularSearchParameters CircularSearchParameters { get; set; }

        // ── RunSpecific ───────────────────────────────────────────────────────

        protected override MyTaskResults RunSpecific(
            string outputFolder,
            List<DbForTask> dbFilenameList,
            List<string> currentRawFileList,
            string taskId,
            FileSpecificParameters[] fileSettingsList)
        {
            MyTaskResults = new MyTaskResults(this);

            // ── Load modifications ────────────────────────────────────────────
            LoadModifications(taskId,
                out var variableModifications,
                out var fixedModifications,
                out var localizeableModificationTypes);

            // ── Load proteins as CircularProteins ─────────────────────────────
            Status("Loading circular protein database...", taskId);

            var dbLoader = new DatabaseLoadingEngine(
                CommonParameters,
                FileSpecificParameters,
                [taskId],
                dbFilenameList,
                taskId,
                CircularSearchParameters.DecoyType,
                CircularSearchParameters.SearchTarget,
                localizeableModificationTypes,
                TargetContaminantAmbiguity.RemoveContaminant);

            var dbResults = (DatabaseLoadingEngineResults)dbLoader.Run();
            var bioPolymers = dbResults.BioPolymers;

            // Wrap every loaded protein as a CircularProtein.
            // Sequences are canonicalized by the CircularProtein constructor.
            var circularProteins = bioPolymers
                .OfType<Protein>()
                .Select(p => CircularProtein.FromProtein(p))
                .ToList();

            if (!circularProteins.Any())
            {
                Warn("No proteins found in the database — circular search will produce no results.");
                return MyTaskResults;
            }

            // ── Mass-difference acceptor ──────────────────────────────────────
            var massDiffAcceptor = GetMassDiffAcceptor(
                CommonParameters.PrecursorMassTolerance,
                CircularSearchParameters.MassDiffAcceptorType,
                CircularSearchParameters.CustomMdac);

            // ── Search each spectra file ──────────────────────────────────────
            var myFileManager = new MyFileManager(
                CircularSearchParameters.DisposeOfFileWhenDone);

            var allPsms = new List<SpectralMatch>();
            int completedFiles = 0;

            // Pre-load first file
            string firstFile = currentRawFileList[0];
            var nextFileTask = new Task<MsDataFile>(
                () => myFileManager.LoadFile(firstFile,
                    SetAllFileSpecificCommonParams(CommonParameters, fileSettingsList[0])));
            nextFileTask.Start();

            for (int spectraFileIndex = 0;
                 spectraFileIndex < currentRawFileList.Count;
                 spectraFileIndex++)
            {
                if (GlobalVariables.StopLoops) break;

                string origDataFile = currentRawFileList[spectraFileIndex];
                StartingDataFile(origDataFile,
                    [taskId, "Individual Spectra Files", origDataFile]);

                var combinedParams = SetAllFileSpecificCommonParams(
                    CommonParameters, fileSettingsList[spectraFileIndex]);

                var thisId = new List<string> { taskId, "Individual Spectra Files", origDataFile };
                NewCollection(Path.GetFileName(origDataFile), thisId);
                Status("Loading spectra file...", thisId);

                nextFileTask.Wait();
                var myMsDataFile = nextFileTask.Result;

                // Begin loading the next file asynchronously
                if (origDataFile != currentRawFileList.Last())
                {
                    int next = spectraFileIndex + 1;
                    nextFileTask = new Task<MsDataFile>(
                        () => myFileManager.LoadFile(
                            currentRawFileList[next],
                            SetAllFileSpecificCommonParams(
                                CommonParameters, fileSettingsList[next])));
                    nextFileTask.Start();
                }

                Status("Getting MS2 scans...", thisId);
                var ms2Scans = GetMs2Scans(myMsDataFile, origDataFile, combinedParams)
                    .OrderBy(b => b.PrecursorMass)
                    .ToArray();

                myFileManager.DoneWithFile(origDataFile);

                var fileSpecificPsms = new SpectralMatch[ms2Scans.Length];

                Status("Searching...", thisId);

                var engine = new CircularSearchEngine(
                    fileSpecificPsms,
                    ms2Scans,
                    variableModifications,
                    fixedModifications,
                    circularProteins,
                    massDiffAcceptor,
                    combinedParams,
                    FileSpecificParameters,
                    thisId,
                    Math.Max(2, CircularSearchParameters.MinInternalFragmentLength));

                engine.Run();

                ReportProgress(new ProgressEventArgs(100, "Done with search!", thisId));

                allPsms.AddRange(fileSpecificPsms.Where(p => p != null));

                completedFiles++;
                FinishedDataFile(origDataFile,
                    [taskId, "Individual Spectra Files", origDataFile]);
                ReportProgress(new ProgressEventArgs(
                    completedFiles / currentRawFileList.Count,
                    "Searching...",
                    [taskId, "Individual Spectra Files"]));
            }

            // ── FDR ──────────────────────────────────────────────────────────
            Status("Calculating FDR...", taskId);

            new FdrAnalysisEngine(
                allPsms, 1 /* notch count */,
                CommonParameters, FileSpecificParameters,
                [taskId]).Run();

            // ── Write results ─────────────────────────────────────────────────
            Status("Writing results...", taskId);

            IEnumerable<SpectralMatch> writablePsms = CircularSearchParameters.WriteHighQValuePsms
                ? allPsms
                : allPsms.Where(p => p.FdrInfo?.QValue <= 0.01);

            if (!CircularSearchParameters.WriteDecoys)
                writablePsms = writablePsms.Where(p => !p.IsDecoy);

            if (!CircularSearchParameters.WriteContaminants)
                writablePsms = writablePsms.Where(p => !p.IsContaminant);

            // PSM-level TSV — inherited protected static from MetaMorpheusTask
            string psmPath = Path.Combine(outputFolder, "CircularPeptideSpectralMatches.psmtsv");
            WritePsmsToTsv(
                writablePsms,
                psmPath,
                new Dictionary<string, int>());
            FinishedWritingFile(psmPath, [taskId]);

            MyTaskResults.AddResultText(
                $"Circular search complete. {allPsms.Count} PSMs found " +
                $"({allPsms.Count(p => p.FdrInfo?.QValue <= 0.01)} at q≤0.01).");

            return MyTaskResults;
        }

        // ── Mass-diff acceptor factory ────────────────────────────────────────

        public static MassDiffAcceptor GetMassDiffAcceptor(
            Tolerance precursorMassTolerance,
            MassDiffAcceptorType type,
            string customMdac)
        {
            return type switch
            {
                MassDiffAcceptorType.Exact =>
                    precursorMassTolerance is PpmTolerance
                        ? new SinglePpmAroundZeroSearchMode(precursorMassTolerance.Value)
                        : new SingleAbsoluteAroundZeroSearchMode(precursorMassTolerance.Value),

                MassDiffAcceptorType.OneMM =>
                    new DotMassDiffAcceptor("1mm",
                        [0, 1.0029], precursorMassTolerance),

                MassDiffAcceptorType.TwoMM =>
                    new DotMassDiffAcceptor("2mm",
                        [0, 1.0029, 2.0052], precursorMassTolerance),

                MassDiffAcceptorType.ThreeMM =>
                    new DotMassDiffAcceptor("3mm",
                        [0, 1.0029, 2.0052, 3.0077], precursorMassTolerance),

                MassDiffAcceptorType.Open =>
                    new OpenSearchMode(),

                MassDiffAcceptorType.ModOpen =>
                    new IntervalMassDiffAcceptor("-187andUp",
                        [new DoubleRange(-187, double.PositiveInfinity)]),

                MassDiffAcceptorType.Custom =>
                    ParseSearchMode(customMdac),

                _ => throw new MetaMorpheusException(
                    $"Unknown MassDiffAcceptorType: {type}")
            };
        }

        private static MassDiffAcceptor ParseSearchMode(string text)
        {
            try
            {
                var split = text.Split(' ');
                return split[1] switch
                {
                    "dot" => ParseDotMode(split),
                    "interval" => ParseIntervalMode(split),
                    "OpenSearch" => new OpenSearchMode(),
                    "daltonsAroundZero" =>
                        new SingleAbsoluteAroundZeroSearchMode(
                            double.Parse(split[2], CultureInfo.InvariantCulture)),
                    "ppmAroundZero" =>
                        new SinglePpmAroundZeroSearchMode(
                            double.Parse(split[2], CultureInfo.InvariantCulture)),
                    _ => throw new MetaMorpheusException(
                        $"Unrecognized search mode type: {split[1]}")
                };
            }
            catch (Exception e)
            {
                throw new MetaMorpheusException(
                    $"Could not parse search mode string: {e.Message}");
            }
        }

        private static MassDiffAcceptor ParseDotMode(string[] split)
        {
            double[] massShifts = split[4].Split(',')
                .Select(p => double.Parse(p, CultureInfo.InvariantCulture))
                .ToArray();
            double toleranceValue = double.Parse(
                split[2].Replace("±", ""), CultureInfo.InvariantCulture);

            Tolerance tol = split[3].ToUpperInvariant() == "PPM"
                ? new PpmTolerance(toleranceValue)
                : new AbsoluteTolerance(toleranceValue);

            return new DotMassDiffAcceptor(split[0], massShifts, tol);
        }

        private static MassDiffAcceptor ParseIntervalMode(string[] split)
        {
            var ranges = Array.ConvertAll(split[2].Split(';'),
                b => new DoubleRange(
                    double.Parse(b.Trim('[', ']').Split(',')[0],
                        CultureInfo.InvariantCulture),
                    double.Parse(b.Trim('[', ']').Split(',')[1],
                        CultureInfo.InvariantCulture)));
            return new IntervalMassDiffAcceptor(split[0], ranges);
        }
    }
}