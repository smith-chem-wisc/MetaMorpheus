using Chemistry;
using EngineLayer;
using EngineLayer.Indexing;
using MassSpectrometry;
using MzLibUtil;
using Nett;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using UsefulProteomicsDatabases;

namespace TaskLayer
{
    public enum MyTask
    {
        Search,
        Gptmd,
        Calibrate,
        XLSearch,
        GlycoSearch
    }

    public abstract class MetaMorpheusTask
    {
        public static readonly TomlSettings tomlConfig = TomlSettings.Create(cfg => cfg
            .ConfigureType<Tolerance>(type => type
                .WithConversionFor<TomlString>(convert => convert
                    .FromToml(tmlString => Tolerance.ParseToleranceString(tmlString.Value))))
            .ConfigureType<PpmTolerance>(type => type
                .WithConversionFor<TomlString>(convert => convert
                    .ToToml(custom => custom.ToString())))
            .ConfigureType<AbsoluteTolerance>(type => type
                .WithConversionFor<TomlString>(convert => convert
                    .ToToml(custom => custom.ToString())))
            .ConfigureType<Protease>(type => type
                .WithConversionFor<TomlString>(convert => convert
                    .ToToml(custom => custom.ToString())
                    .FromToml(tmlString => ProteaseDictionary.Dictionary[tmlString.Value])))
            .ConfigureType<List<string>>(type => type
                    .WithConversionFor<TomlString>(convert => convert
                        .ToToml(custom => string.Join("\t", custom))
                        .FromToml(tmlString => GetModsTypesFromString(tmlString.Value))))
            .ConfigureType<List<(string, string)>>(type => type
                    .WithConversionFor<TomlString>(convert => convert
                        .ToToml(custom => string.Join("\t\t", custom.Select(b => b.Item1 + "\t" + b.Item2)))
                        .FromToml(tmlString => GetModsFromString(tmlString.Value)))));

        protected readonly StringBuilder ProseCreatedWhileRunning = new StringBuilder();

        protected MyTaskResults MyTaskResults;

        protected MetaMorpheusTask(MyTask taskType)
        {
            this.TaskType = taskType;
        }

        public static event EventHandler<SingleTaskEventArgs> FinishedSingleTaskHandler;

        public static event EventHandler<SingleFileEventArgs> FinishedWritingFileHandler;

        public static event EventHandler<SingleTaskEventArgs> StartingSingleTaskHander;

        public static event EventHandler<StringEventArgs> StartingDataFileHandler;

        public static event EventHandler<StringEventArgs> FinishedDataFileHandler;

        public static event EventHandler<StringEventArgs> OutLabelStatusHandler;

        public static event EventHandler<StringEventArgs> WarnHandler;

        public static event EventHandler<StringEventArgs> LogHandler;

        public static event EventHandler<StringEventArgs> NewCollectionHandler;

        public static event EventHandler<ProgressEventArgs> OutProgressHandler;

        public MyTask TaskType { get; set; }

        public CommonParameters CommonParameters { get; set; }

        public const string IndexFolderName = "DatabaseIndex";
        public const string IndexEngineParamsFileName = "indexEngine.params";
        public const string PeptideIndexFileName = "peptideIndex.ind";
        public const string FragmentIndexFileName = "fragmentIndex.ind";
        public const string SecondIndexEngineParamsFileName = "secondIndexEngine.params";
        public const string SecondFragmentIndexFileName = "secondFragmentIndex.ind";
        public const string PrecursorIndexFileName = "precursorIndex.ind";

        public static IEnumerable<Ms2ScanWithSpecificMass> GetMs2Scans(MsDataFile myMSDataFile, string fullFilePath, CommonParameters commonParameters)
        {
            var msNScans = myMSDataFile.GetAllScansList().Where(x => x.MsnOrder > 1).ToArray();
            var ms2Scans = msNScans.Where(p => p.MsnOrder == 2).ToArray();
            var ms3Scans = msNScans.Where(p => p.MsnOrder == 3).ToArray();
            List<Ms2ScanWithSpecificMass>[] scansWithPrecursors = new List<Ms2ScanWithSpecificMass>[ms2Scans.Length];

            if (!ms2Scans.Any())
            {
                return new List<Ms2ScanWithSpecificMass>();
            }

            Parallel.ForEach(Partitioner.Create(0, ms2Scans.Length), new ParallelOptions { MaxDegreeOfParallelism = commonParameters.MaxThreadsToUsePerFile },
                (partitionRange, loopState) =>
                {
                    for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                    {
                        if (GlobalVariables.StopLoops) { break; }

                        MsDataScan ms2scan = ms2Scans[i];

                        List<(double, int)> precursors = new List<(double, int)>();
                        if (ms2scan.OneBasedPrecursorScanNumber.HasValue)
                        {
                            var precursorSpectrum = myMSDataFile.GetOneBasedScan(ms2scan.OneBasedPrecursorScanNumber.Value);

                            try
                            {
                                ms2scan.RefineSelectedMzAndIntensity(precursorSpectrum.MassSpectrum);
                            }
                            catch (MzLibException ex)
                            {
                                Warn("Could not get precursor ion for MS2 scan #" + ms2scan.OneBasedScanNumber + "; " + ex.Message);
                                continue;
                            }

                            if (ms2scan.SelectedIonMonoisotopicGuessMz.HasValue)
                            {
                                ms2scan.ComputeMonoisotopicPeakIntensity(precursorSpectrum.MassSpectrum);
                            }

                            if (commonParameters.DoPrecursorDeconvolution)
                            {
                                foreach (var envelope in ms2scan.GetIsolatedMassesAndCharges(
                                    precursorSpectrum.MassSpectrum, 1,
                                    commonParameters.DeconvolutionMaxAssumedChargeState,
                                    commonParameters.DeconvolutionMassTolerance.Value,
                                    commonParameters.DeconvolutionIntensityRatio))
                                {
                                    var monoPeakMz = envelope.monoisotopicMass.ToMz(envelope.charge);
                                    precursors.Add((monoPeakMz, envelope.charge));
                                }
                            }
                        }

                        if (commonParameters.UseProvidedPrecursorInfo && ms2scan.SelectedIonChargeStateGuess.HasValue)
                        {
                            var precursorCharge = ms2scan.SelectedIonChargeStateGuess.Value;
                            if (ms2scan.SelectedIonMonoisotopicGuessMz.HasValue)
                            {
                                double precursorMZ = ms2scan.SelectedIonMonoisotopicGuessMz.Value;
                                if (!precursors.Any(b =>
                                    commonParameters.DeconvolutionMassTolerance.Within(
                                        precursorMZ.ToMass(precursorCharge), b.Item1.ToMass(b.Item2))))
                                {
                                    precursors.Add((precursorMZ, precursorCharge));
                                }
                            }
                            else
                            {
                                double precursorMZ = ms2scan.SelectedIonMZ.Value;
                                if (!precursors.Any(b =>
                                    commonParameters.DeconvolutionMassTolerance.Within(
                                        precursorMZ.ToMass(precursorCharge), b.Item1.ToMass(b.Item2))))
                                {
                                    precursors.Add((precursorMZ, precursorCharge));
                                }
                            }
                        }

                        scansWithPrecursors[i] = new List<Ms2ScanWithSpecificMass>();
                        IsotopicEnvelope[] neutralExperimentalFragments = null;

                        if (commonParameters.DissociationType != DissociationType.LowCID)
                        {
                            neutralExperimentalFragments = Ms2ScanWithSpecificMass.GetNeutralExperimentalFragments(ms2scan, commonParameters);
                        }

                        // get child scans
                        List<MsDataScan> ms2ChildScans = new List<MsDataScan>();
                        List<MsDataScan> ms3ChildScans = new List<MsDataScan>();
                        if (commonParameters.ChildScanDissociationType != DissociationType.Unknown)
                        {
                            ms3ChildScans = ms3Scans.Where(p => p.OneBasedPrecursorScanNumber == ms2scan.OneBasedScanNumber).ToList();

                            ms2ChildScans = ms2Scans.Where(p => (p.OneBasedPrecursorScanNumber == ms2scan.OneBasedPrecursorScanNumber || p.OneBasedPrecursorScanNumber== ms2scan.OneBasedScanNumber)
                                && p.OneBasedScanNumber > ms2scan.OneBasedScanNumber
                                && Math.Abs(p.IsolationMz.Value - ms2scan.IsolationMz.Value) < 0.01).ToList();
                        }

                        foreach (var precursor in precursors)
                        {
                            var scan = new Ms2ScanWithSpecificMass(ms2scan, precursor.Item1,
                                precursor.Item2, fullFilePath, commonParameters, neutralExperimentalFragments);

                            foreach (var ms2ChildScan in ms2ChildScans)
                            {
                                IsotopicEnvelope[] childNeutralExperimentalFragments = null;

                                if (commonParameters.ChildScanDissociationType != DissociationType.LowCID)
                                {
                                    childNeutralExperimentalFragments = Ms2ScanWithSpecificMass.GetNeutralExperimentalFragments(ms2ChildScan, commonParameters);
                                }

                                scan.ChildScans.Add(new Ms2ScanWithSpecificMass(ms2ChildScan, precursor.Item1,
                                    precursor.Item2, fullFilePath, commonParameters, childNeutralExperimentalFragments));
                            }

                            foreach (var ms3ChildScan in ms3ChildScans)
                            {
                                int precursorCharge = 1;
                                var precursorSpectrum = ms2scan;

                                try
                                {
                                    ms3ChildScan.RefineSelectedMzAndIntensity(precursorSpectrum.MassSpectrum);
                                }
                                catch (MzLibException ex)
                                {
                                    Warn("Could not get precursor ion for MS3 scan #" + ms3ChildScan.OneBasedScanNumber + "; " + ex.Message);
                                    continue;
                                }

                                if (ms3ChildScan.SelectedIonMonoisotopicGuessMz.HasValue)
                                {
                                    ms3ChildScan.ComputeMonoisotopicPeakIntensity(precursorSpectrum.MassSpectrum);
                                }

                                if (ms3ChildScan.SelectedIonChargeStateGuess.HasValue)
                                {
                                    precursorCharge = ms3ChildScan.SelectedIonChargeStateGuess.Value;
                                }
                                if (!ms3ChildScan.SelectedIonMonoisotopicGuessMz.HasValue)
                                {
                                    Warn("Could not get precursor ion m/z for MS3 scan #" + ms3ChildScan.OneBasedScanNumber);
                                    continue;
                                }

                                IsotopicEnvelope[] childNeutralExperimentalFragments = null;

                                if (commonParameters.ChildScanDissociationType != DissociationType.LowCID)
                                {
                                    childNeutralExperimentalFragments = Ms2ScanWithSpecificMass.GetNeutralExperimentalFragments(ms3ChildScan, commonParameters);
                                }

                                scan.ChildScans.Add(new Ms2ScanWithSpecificMass(ms3ChildScan, ms3ChildScan.SelectedIonMonoisotopicGuessMz.Value,
                                    ms3ChildScan.SelectedIonChargeStateGuess.Value, fullFilePath, commonParameters, childNeutralExperimentalFragments));
                            }

                            scansWithPrecursors[i].Add(scan);
                        }
                    }
                });

            var childScanNumbers = new HashSet<int>(scansWithPrecursors.SelectMany(p => p.SelectMany(v => v.ChildScans.Select(x => x.OneBasedScanNumber))));
            var parentScans = scansWithPrecursors.Where(p => p.Any() && !childScanNumbers.Contains(p.First().OneBasedScanNumber)).SelectMany(v => v);

            // XCorr pre-processing for low-res data. this is here because the parent/child scans may have different
            // resolutions, so this pre-processing must take place after the parent/child scans have been determined
            foreach (var parentScan in parentScans)
            {
                if (commonParameters.DissociationType == DissociationType.LowCID && !parentScan.TheScan.MassSpectrum.XcorrProcessed)
                {
                    parentScan.TheScan.MassSpectrum.XCorrPrePreprocessing(0, 1969, parentScan.TheScan.IsolationMz.Value);
                }

                foreach (var childScan in parentScan.ChildScans)
                {
                    if (commonParameters.ChildScanDissociationType == DissociationType.LowCID && !childScan.TheScan.MassSpectrum.XcorrProcessed)
                    {
                        childScan.TheScan.MassSpectrum.XCorrPrePreprocessing(0, 1969, childScan.TheScan.IsolationMz.Value);
                    }
                }
            }

            return parentScans;
        }

        public static CommonParameters SetAllFileSpecificCommonParams(CommonParameters commonParams, FileSpecificParameters fileSpecificParams)
        {
            if (fileSpecificParams == null)
            {
                return commonParams;
            }

            // set file-specific digestion parameters
            Protease protease = fileSpecificParams.Protease ?? commonParams.DigestionParams.SpecificProtease; //set to specific for nonspecific searches to update
            int minPeptideLength = fileSpecificParams.MinPeptideLength ?? commonParams.DigestionParams.MinPeptideLength;
            int maxPeptideLength = fileSpecificParams.MaxPeptideLength ?? commonParams.DigestionParams.MaxPeptideLength;
            int maxMissedCleavages = fileSpecificParams.MaxMissedCleavages ?? commonParams.DigestionParams.MaxMissedCleavages;
            int maxModsForPeptide = fileSpecificParams.MaxModsForPeptide ?? commonParams.DigestionParams.MaxModsForPeptide;
            DigestionParams fileSpecificDigestionParams = new DigestionParams(
                protease: protease.Name,
                maxMissedCleavages: maxMissedCleavages,
                minPeptideLength: minPeptideLength,
                maxPeptideLength: maxPeptideLength,
                maxModsForPeptides: maxModsForPeptide,

                //NEED THESE OR THEY'LL BE OVERWRITTEN
                maxModificationIsoforms: commonParams.DigestionParams.MaxModificationIsoforms,
                initiatorMethionineBehavior: commonParams.DigestionParams.InitiatorMethionineBehavior,
                fragmentationTerminus: commonParams.DigestionParams.FragmentationTerminus,
                searchModeType: commonParams.DigestionParams.SearchModeType
                );

            // set the rest of the file-specific parameters
            Tolerance precursorMassTolerance = fileSpecificParams.PrecursorMassTolerance ?? commonParams.PrecursorMassTolerance;
            Tolerance productMassTolerance = fileSpecificParams.ProductMassTolerance ?? commonParams.ProductMassTolerance;
            DissociationType dissociationType = fileSpecificParams.DissociationType ?? commonParams.DissociationType;

            CommonParameters returnParams = new CommonParameters(
                dissociationType: dissociationType,
                precursorMassTolerance: precursorMassTolerance,
                productMassTolerance: productMassTolerance,
                digestionParams: fileSpecificDigestionParams,

                //NEED THESE OR THEY'LL BE OVERWRITTEN
                childScanDissociationType: commonParams.ChildScanDissociationType,
                doPrecursorDeconvolution: commonParams.DoPrecursorDeconvolution,
                useProvidedPrecursorInfo: commonParams.UseProvidedPrecursorInfo,
                deconvolutionIntensityRatio: commonParams.DeconvolutionIntensityRatio,
                deconvolutionMaxAssumedChargeState: commonParams.DeconvolutionMaxAssumedChargeState,
                reportAllAmbiguity: commonParams.ReportAllAmbiguity,
                addCompIons: commonParams.AddCompIons,
                totalPartitions: commonParams.TotalPartitions,
                scoreCutoff: commonParams.ScoreCutoff,
                numberOfPeaksToKeepPerWindow: commonParams.NumberOfPeaksToKeepPerWindow,
                minimumAllowedIntensityRatioToBasePeak: commonParams.MinimumAllowedIntensityRatioToBasePeak,
                windowWidthThomsons: commonParams.WindowWidthThomsons,
                trimMs1Peaks: commonParams.TrimMs1Peaks,
                trimMsMsPeaks: commonParams.TrimMsMsPeaks,
                normalizePeaksAccrossAllWindows: commonParams.NormalizePeaksAccrossAllWindows,
                useDeltaScore: commonParams.UseDeltaScore,
                deconvolutionMassTolerance: commonParams.DeconvolutionMassTolerance,
                maxThreadsToUsePerFile: commonParams.MaxThreadsToUsePerFile,
                listOfModsVariable: commonParams.ListOfModsVariable,
                listOfModsFixed: commonParams.ListOfModsFixed,
                qValueOutputFilter: commonParams.QValueOutputFilter,
                taskDescriptor: commonParams.TaskDescriptor,
                assumeOrphanPeaksAreZ1Fragments: commonParams.AssumeOrphanPeaksAreZ1Fragments,
                maxHeterozygousVariants: commonParams.MaxHeterozygousVariants,
                minVariantDepth: commonParams.MinVariantDepth);

            return returnParams;
        }

        public MyTaskResults RunTask(string output_folder, List<DbForTask> currentProteinDbFilenameList, List<string> currentRawDataFilepathList, string displayName)
        {
            StartingSingleTask(displayName);

            var tomlFileName = Path.Combine(Directory.GetParent(output_folder).ToString(), "Task Settings", displayName + "config.toml");
            Toml.WriteFile(this, tomlFileName, tomlConfig);
            FinishedWritingFile(tomlFileName, new List<string> { displayName });

            MetaMorpheusEngine.FinishedSingleEngineHandler += SingleEngineHandlerInTask;
            try
            {
                var stopWatch = new Stopwatch();
                stopWatch.Start();

                FileSpecificParameters[] fileSettingsList = new FileSpecificParameters[currentRawDataFilepathList.Count];
                for (int i = 0; i < currentRawDataFilepathList.Count; i++)
                {
                    if (GlobalVariables.StopLoops) { break; }
                    string rawFilePath = currentRawDataFilepathList[i];
                    string directory = Directory.GetParent(rawFilePath).ToString();
                    string fileSpecificTomlPath = Path.Combine(directory, Path.GetFileNameWithoutExtension(rawFilePath)) + ".toml";
                    if (File.Exists(fileSpecificTomlPath))
                    {
                        TomlTable fileSpecificSettings = Toml.ReadFile(fileSpecificTomlPath, tomlConfig);
                        try
                        {
                            fileSettingsList[i] = new FileSpecificParameters(fileSpecificSettings);
                        }
                        catch (MetaMorpheusException e)
                        {
                            //file - specific toml has already been validated in the GUI when the spectra files were added, so...
                            // probably the only time you can get here is if the user modifies the file-specific parameter file in the middle of a run...
                            Warn("Problem parsing the file-specific toml " + Path.GetFileName(fileSpecificTomlPath) + "; " + e.Message + "; is the toml from an older version of MetaMorpheus?");
                        }
                    }
                }

                RunSpecific(output_folder, currentProteinDbFilenameList, currentRawDataFilepathList, displayName, fileSettingsList);
                stopWatch.Stop();
                MyTaskResults.Time = stopWatch.Elapsed;
                var resultsFileName = Path.Combine(output_folder, "results.txt");
                using (StreamWriter file = new StreamWriter(resultsFileName))
                {
                    file.WriteLine("MetaMorpheus: version " + GlobalVariables.MetaMorpheusVersion);
                    file.Write(MyTaskResults.ToString());
                }
                FinishedWritingFile(resultsFileName, new List<string> { displayName });
                FinishedSingleTask(displayName);
            }
            catch (Exception e)
            {
                MetaMorpheusEngine.FinishedSingleEngineHandler -= SingleEngineHandlerInTask;
                var resultsFileName = Path.Combine(output_folder, "results.txt");
                e.Data.Add("folder", output_folder);
                using (StreamWriter file = new StreamWriter(resultsFileName))
                {
                    file.WriteLine(GlobalVariables.MetaMorpheusVersion.Equals("1.0.0.0") ? "MetaMorpheus: Not a release version" : "MetaMorpheus: version " + GlobalVariables.MetaMorpheusVersion);
                    file.WriteLine(SystemInfo.CompleteSystemInfo()); //OS, OS Version, .Net Version, RAM, processor count, MSFileReader .dll versions X3
                    file.Write("e: " + e);
                    file.Write("e.Message: " + e.Message);
                    file.Write("e.InnerException: " + e.InnerException);
                    file.Write("e.Source: " + e.Source);
                    file.Write("e.StackTrace: " + e.StackTrace);
                    file.Write("e.TargetSite: " + e.TargetSite);
                }
                throw;
            }

            {
                var proseFilePath = Path.Combine(output_folder, "prose.txt");
                using (StreamWriter file = new StreamWriter(proseFilePath))
                {
                    file.Write("The data analysis was performed using MetaMorpheus version " + GlobalVariables.MetaMorpheusVersion + ", available at " + "https://github.com/smith-chem-wisc/MetaMorpheus.");
                    file.Write(ProseCreatedWhileRunning.ToString());
                    file.Write(SystemInfo.SystemProse().Replace(Environment.NewLine, "") + " ");
                    file.WriteLine("The total time to perform the " + TaskType + " task on " + currentRawDataFilepathList.Count + " spectra file(s) was " + String.Format("{0:0.00}", MyTaskResults.Time.TotalMinutes) + " minutes.");
                    file.WriteLine();
                    file.WriteLine("Published works using MetaMorpheus software are encouraged to cite: Solntsev, S. K.; Shortreed, M. R.; Frey, B. L.; Smith, L. M. Enhanced Global Post-translational Modification Discovery with MetaMorpheus. Journal of Proteome Research. 2018, 17 (5), 1844-1851.");

                    file.WriteLine();
                    file.WriteLine("Spectra files: ");
                    file.WriteLine(string.Join(Environment.NewLine, currentRawDataFilepathList.Select(b => '\t' + b)));
                    file.WriteLine("Databases:");
                    file.Write(string.Join(Environment.NewLine, currentProteinDbFilenameList.Select(b => '\t' + (b.IsContaminant ? "Contaminant " : "") + b.FilePath)));
                }
                FinishedWritingFile(proseFilePath, new List<string> { displayName });
            }

            MetaMorpheusEngine.FinishedSingleEngineHandler -= SingleEngineHandlerInTask;
            return MyTaskResults;
        }

        protected List<Protein> LoadProteins(string taskId, List<DbForTask> dbFilenameList, bool searchTarget, DecoyType decoyType, List<string> localizeableModificationTypes, CommonParameters commonParameters)
        {
            Status("Loading proteins...", new List<string> { taskId });
            int emptyProteinEntries = 0;
            List<Protein> proteinList = new List<Protein>();
            foreach (var db in dbFilenameList)
            {
                var dbProteinList = LoadProteinDb(db.FilePath, searchTarget, decoyType, localizeableModificationTypes, db.IsContaminant, out Dictionary<string, Modification> unknownModifications, out int emptyProteinEntriesForThisDb, commonParameters);
                proteinList = proteinList.Concat(dbProteinList).ToList();
                emptyProteinEntries += emptyProteinEntriesForThisDb;
            }
            if (!proteinList.Any())
            {
                Warn("Warning: No protein entries were found in the database");
            }
            else if (emptyProteinEntries > 0)
            {
                Warn("Warning: " + emptyProteinEntries + " empty protein entries ignored");
            }
            return proteinList;
        }

        protected static List<Protein> LoadProteinDb(string fileName, bool generateTargets, DecoyType decoyType, List<string> localizeableModificationTypes, bool isContaminant, out Dictionary<string, Modification> um,
            out int emptyEntriesCount, CommonParameters commonParameters)
        {
            List<string> dbErrors = new List<string>();
            List<Protein> proteinList = new List<Protein>();

            string theExtension = Path.GetExtension(fileName).ToLowerInvariant();
            bool compressed = theExtension.EndsWith("gz"); // allows for .bgz and .tgz, too which are used on occasion
            theExtension = compressed ? Path.GetExtension(Path.GetFileNameWithoutExtension(fileName)).ToLowerInvariant() : theExtension;

            if (theExtension.Equals(".fasta") || theExtension.Equals(".fa"))
            {
                um = null;
                proteinList = ProteinDbLoader.LoadProteinFasta(fileName, generateTargets, decoyType, isContaminant, ProteinDbLoader.UniprotAccessionRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotGeneNameRegex,
                    ProteinDbLoader.UniprotOrganismRegex, out dbErrors, commonParameters.MaxThreadsToUsePerFile);
            }
            else
            {
                List<string> modTypesToExclude = GlobalVariables.AllModTypesKnown.Where(b => !localizeableModificationTypes.Contains(b)).ToList();
                proteinList = ProteinDbLoader.LoadProteinXML(fileName, generateTargets, decoyType, GlobalVariables.AllModsKnown, isContaminant, modTypesToExclude, out um, commonParameters.MaxThreadsToUsePerFile, commonParameters.MaxHeterozygousVariants, commonParameters.MinVariantDepth);
            }

            emptyEntriesCount = proteinList.Count(p => p.BaseSequence.Length == 0);
            return proteinList.Where(p => p.BaseSequence.Length > 0).ToList();
        }

        protected void LoadModifications(string taskId, out List<Modification> variableModifications, out List<Modification> fixedModifications, out List<string> localizableModificationTypes)
        {
            // load modifications
            Status("Loading modifications...", taskId);
            variableModifications = GlobalVariables.AllModsKnown.OfType<Modification>().Where(b => CommonParameters.ListOfModsVariable.Contains((b.ModificationType, b.IdWithMotif))).ToList();
            fixedModifications = GlobalVariables.AllModsKnown.OfType<Modification>().Where(b => CommonParameters.ListOfModsFixed.Contains((b.ModificationType, b.IdWithMotif))).ToList();
            localizableModificationTypes = GlobalVariables.AllModTypesKnown.ToList();

            var recognizedVariable = variableModifications.Select(p => p.IdWithMotif);
            var recognizedFixed = fixedModifications.Select(p => p.IdWithMotif);
            var unknownMods = CommonParameters.ListOfModsVariable.Select(p => p.Item2).Except(recognizedVariable).ToList();
            unknownMods.AddRange(CommonParameters.ListOfModsFixed.Select(p => p.Item2).Except(recognizedFixed));
            foreach (var unrecognizedMod in unknownMods)
            {
                Warn("Unrecognized mod " + unrecognizedMod + "; are you using an old .toml?");
            }
        }

        protected static void WritePsmsToTsv(IEnumerable<PeptideSpectralMatch> psms, string filePath, IReadOnlyDictionary<string, int> modstoWritePruned)
        {
            using (StreamWriter output = new StreamWriter(filePath))
            {
                output.WriteLine(PeptideSpectralMatch.GetTabSeparatedHeader());
                foreach (var psm in psms)
                {
                    output.WriteLine(psm.ToString(modstoWritePruned));
                }
            }
        }

        protected void ReportProgress(ProgressEventArgs v)
        {
            OutProgressHandler?.Invoke(this, v);
        }

        protected abstract MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId, FileSpecificParameters[] fileSettingsList);

        protected void FinishedWritingFile(string path, List<string> nestedIDs)
        {
            FinishedWritingFileHandler?.Invoke(this, new SingleFileEventArgs(path, nestedIDs));
        }

        protected void StartingDataFile(string v, List<string> nestedIDs)
        {
            StartingDataFileHandler?.Invoke(this, new StringEventArgs(v, nestedIDs));
        }

        protected void FinishedDataFile(string v, List<string> nestedIDs)
        {
            FinishedDataFileHandler?.Invoke(this, new StringEventArgs(v, nestedIDs));
        }

        protected void Status(string v, string id)
        {
            OutLabelStatusHandler?.Invoke(this, new StringEventArgs(v, new List<string> { id }));
        }

        protected void Status(string v, List<string> nestedIds)
        {
            OutLabelStatusHandler?.Invoke(this, new StringEventArgs(v, nestedIds));
        }

        protected static void Warn(string v)
        {
            WarnHandler?.Invoke(null, new StringEventArgs(v, null));
        }

        protected void Log(string v, List<string> nestedIds)
        {
            LogHandler?.Invoke(this, new StringEventArgs(v, nestedIds));
        }

        protected void NewCollection(string displayName, List<string> nestedIds)
        {
            NewCollectionHandler?.Invoke(this, new StringEventArgs(displayName, nestedIds));
        }

        private static List<string> GetModsTypesFromString(string value)
        {
            return value.Split(new string[] { "\t" }, StringSplitOptions.RemoveEmptyEntries).ToList();
        }

        private static List<(string, string)> GetModsFromString(string value)
        {
            return value.Split(new string[] { "\t\t" }, StringSplitOptions.RemoveEmptyEntries).Select(b => (b.Split('\t').First(), b.Split('\t').Last())).ToList();
        }

        private void SingleEngineHandlerInTask(object sender, SingleEngineFinishedEventArgs e)
        {
            MyTaskResults.AddResultText(e.ToString());
        }

        private void FinishedSingleTask(string displayName)
        {
            FinishedSingleTaskHandler?.Invoke(this, new SingleTaskEventArgs(displayName));
        }

        private void StartingSingleTask(string displayName)
        {
            StartingSingleTaskHander?.Invoke(this, new SingleTaskEventArgs(displayName));
        }

        private static IEnumerable<Type> GetSubclassesAndItself(Type type)
        {
            yield return type;
        }

        private static bool SameSettings(string pathToOldParamsFile, IndexingEngine indexEngine)
        {
            using (StreamReader reader = new StreamReader(pathToOldParamsFile))
            {
                if (reader.ReadToEnd().Equals(indexEngine.ToString()))
                {
                    return true;
                }
            }
            return false;
        }

        private static void WritePeptideIndex(List<PeptideWithSetModifications> peptideIndex, string peptideIndexFileName)
        {
            var messageTypes = GetSubclassesAndItself(typeof(List<PeptideWithSetModifications>));
            var ser = new NetSerializer.Serializer(messageTypes);

            using (var file = File.Create(peptideIndexFileName))
            {
                ser.Serialize(file, peptideIndex);
            }
        }

        private static List<PeptideWithSetModifications> ReadPeptideIndex(string peptideIndexFileName, List<Protein> allKnownProteins)
        {
            var messageTypes = GetSubclassesAndItself(typeof(List<PeptideWithSetModifications>));
            var ser = new NetSerializer.Serializer(messageTypes);
            List<PeptideWithSetModifications> peptideIndex;
            using (var file = File.OpenRead(peptideIndexFileName))
            {
                peptideIndex = (List<PeptideWithSetModifications>)ser.Deserialize(file);
            }

            // populate dictionaries of known proteins for deserialization
            Dictionary<string, Protein> proteinDictionary = new Dictionary<string, Protein>();
            foreach (Protein protein in allKnownProteins)
            {
                if (!proteinDictionary.ContainsKey(protein.Accession))
                {
                    proteinDictionary.Add(protein.Accession, protein);
                }
                else if (proteinDictionary[protein.Accession].BaseSequence != protein.BaseSequence)
                {
                    throw new MetaMorpheusException($"The protein database contained multiple proteins with accession {protein.Accession} ! This is not allowed for index-based searches (modern, non-specific, crosslink searches)");
                }
            }

            // get non-serialized information for the peptides (proteins, mod info)
            foreach (var peptide in peptideIndex)
            {
                peptide.SetNonSerializedPeptideInfo(GlobalVariables.AllModsKnownDictionary, proteinDictionary);
            }

            return peptideIndex;
        }

        private static void WriteFragmentIndex(List<int>[] fragmentIndex, string fragmentIndexFileName)
        {
            var messageTypes = GetSubclassesAndItself(typeof(List<int>[]));
            var ser = new NetSerializer.Serializer(messageTypes);

            using (var file = File.Create(fragmentIndexFileName))
            {
                ser.Serialize(file, fragmentIndex);
            }
        }

        private static List<int>[] ReadFragmentIndex(string fragmentIndexFileName)
        {
            var messageTypes = GetSubclassesAndItself(typeof(List<int>[]));
            var ser = new NetSerializer.Serializer(messageTypes);

            using (var file = File.OpenRead(fragmentIndexFileName))
            {
                return (List<int>[])ser.Deserialize(file);
            }
        }

        private static string GetExistingFolderWithIndices(IndexingEngine indexEngine, List<DbForTask> dbFilenameList)
        {
            foreach (var database in dbFilenameList)
            {
                string baseDir = Path.GetDirectoryName(database.FilePath);
                DirectoryInfo indexDirectory = new DirectoryInfo(Path.Combine(baseDir, IndexFolderName));

                if (!Directory.Exists(indexDirectory.FullName))
                {
                    return null;
                }

                // all directories in the same directory as the protein database
                DirectoryInfo[] directories = indexDirectory.GetDirectories();

                // look in each subdirectory to find indexes folder
                foreach (DirectoryInfo possibleFolder in directories)
                {
                    string result = CheckFiles(indexEngine, possibleFolder);

                    if (result != null)
                    {
                        return result;
                    }
                }
            }

            return null;
        }

        private static string CheckFiles(IndexingEngine indexEngine, DirectoryInfo folder)
        {
            if (File.Exists(Path.Combine(folder.FullName, IndexEngineParamsFileName)) &&
                File.Exists(Path.Combine(folder.FullName, PeptideIndexFileName)) &&
                File.Exists(Path.Combine(folder.FullName, FragmentIndexFileName)) &&
                (File.Exists(Path.Combine(folder.FullName, PrecursorIndexFileName)) || !indexEngine.GeneratePrecursorIndex) &&
                SameSettings(Path.Combine(folder.FullName, IndexEngineParamsFileName), indexEngine))
            {
                return folder.FullName;
            }
            return null;
        }

        private static void WriteIndexEngineParams(IndexingEngine indexEngine, string fileName)
        {
            using (StreamWriter output = new StreamWriter(fileName))
            {
                output.Write(indexEngine);
            }
        }

        private static string GenerateOutputFolderForIndices(List<DbForTask> dbFilenameList)
        {
            var pathToIndexes = Path.Combine(Path.GetDirectoryName(dbFilenameList.First().FilePath), IndexFolderName);
            if (!File.Exists(pathToIndexes))
            {
                Directory.CreateDirectory(pathToIndexes);
            }
            var folder = Path.Combine(pathToIndexes, DateTime.Now.ToString("yyyy-MM-dd-HH-mm-ss", CultureInfo.InvariantCulture));
            Directory.CreateDirectory(folder);
            return folder;
        }

        public void GenerateIndexes(IndexingEngine indexEngine, List<DbForTask> dbFilenameList, ref List<PeptideWithSetModifications> peptideIndex, ref List<int>[] fragmentIndex, ref List<int>[] precursorIndex, List<Protein> allKnownProteins, string taskId)
        {
            string pathToFolderWithIndices = GetExistingFolderWithIndices(indexEngine, dbFilenameList);
            if (pathToFolderWithIndices == null) //if no indexes exist
            {
                var output_folderForIndices = GenerateOutputFolderForIndices(dbFilenameList);
                Status("Writing params...", new List<string> { taskId });
                var paramsFile = Path.Combine(output_folderForIndices, IndexEngineParamsFileName);
                WriteIndexEngineParams(indexEngine, paramsFile);
                FinishedWritingFile(paramsFile, new List<string> { taskId });

                Status("Running Index Engine...", new List<string> { taskId });
                var indexResults = (IndexingResults)indexEngine.Run();
                peptideIndex = indexResults.PeptideIndex;
                fragmentIndex = indexResults.FragmentIndex;
                precursorIndex = indexResults.PrecursorIndex;

                Status("Writing peptide index...", new List<string> { taskId });
                var peptideIndexFile = Path.Combine(output_folderForIndices, PeptideIndexFileName);
                WritePeptideIndex(peptideIndex, peptideIndexFile);
                FinishedWritingFile(peptideIndexFile, new List<string> { taskId });

                Status("Writing fragment index...", new List<string> { taskId });
                var fragmentIndexFile = Path.Combine(output_folderForIndices, FragmentIndexFileName);
                WriteFragmentIndex(fragmentIndex, fragmentIndexFile);
                FinishedWritingFile(fragmentIndexFile, new List<string> { taskId });

                if (indexEngine.GeneratePrecursorIndex) //If a precursor index is specified (used for speedy semi and non-specific searches)
                {
                    Status("Writing precursor index...", new List<string> { taskId });
                    var precursorIndexFile = Path.Combine(output_folderForIndices, PrecursorIndexFileName);
                    WriteFragmentIndex(precursorIndex, precursorIndexFile);
                    FinishedWritingFile(precursorIndexFile, new List<string> { taskId });
                }
            }
            else //if we found indexes with the same params
            {
                Status("Reading peptide index...", new List<string> { taskId });
                peptideIndex = ReadPeptideIndex(Path.Combine(pathToFolderWithIndices, PeptideIndexFileName), allKnownProteins);

                Status("Reading fragment index...", new List<string> { taskId });
                fragmentIndex = ReadFragmentIndex(Path.Combine(pathToFolderWithIndices, FragmentIndexFileName));

                if (indexEngine.GeneratePrecursorIndex)
                {
                    Status("Reading precursor index...", new List<string> { taskId });
                    precursorIndex = ReadFragmentIndex(Path.Combine(pathToFolderWithIndices, PrecursorIndexFileName));
                }
            }
        }

        public void GenerateSecondIndexes(IndexingEngine indexEngine, IndexingEngine secondIndexEngine, List<DbForTask> dbFilenameList, ref List<int>[] secondFragmentIndex, List<Protein> allKnownProteins, string taskId)
        {
            string pathToFolderWithIndices = GetExistingFolderWithIndices(indexEngine, dbFilenameList);
            if (!File.Exists(Path.Combine(pathToFolderWithIndices, SecondFragmentIndexFileName))) //if no indexes exist
            {
                Status("Writing params...", new List<string> { taskId });
                var paramsFile = Path.Combine(pathToFolderWithIndices, SecondIndexEngineParamsFileName);
                WriteIndexEngineParams(secondIndexEngine, paramsFile);
                FinishedWritingFile(paramsFile, new List<string> { taskId });

                Status("Running Index Engine...", new List<string> { taskId });
                var indexResults = (IndexingResults)secondIndexEngine.Run();
                secondFragmentIndex = indexResults.FragmentIndex;

                Status("Writing fragment index...", new List<string> { taskId });
                var fragmentIndexFile = Path.Combine(pathToFolderWithIndices, SecondFragmentIndexFileName);
                WriteFragmentIndex(secondFragmentIndex, fragmentIndexFile);
                FinishedWritingFile(fragmentIndexFile, new List<string> { taskId });
            }
            else //if we found indexes with the same params
            {
                Status("Reading fragment index...", new List<string> { taskId });
                secondFragmentIndex = ReadFragmentIndex(Path.Combine(pathToFolderWithIndices, SecondFragmentIndexFileName));
            }
        }
    }
}