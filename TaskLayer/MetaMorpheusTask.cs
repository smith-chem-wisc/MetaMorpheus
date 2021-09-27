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

        protected string OutputFolder { get; private set; }

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
        public List<(string FileName, CommonParameters Parameters)> FileSpecificParameters { get; set; }

        public const string IndexFolderName = "DatabaseIndex";
        public const string IndexEngineParamsFileName = "indexEngine.params";
        public const string PeptideIndexFileName = "peptideIndex.ind";
        public const string FragmentIndexFileName = "fragmentIndex.ind";
        public const string SecondIndexEngineParamsFileName = "secondIndexEngine.params";
        public const string SecondFragmentIndexFileName = "secondFragmentIndex.ind";
        public const string PrecursorIndexFileName = "precursorIndex.ind";

        public static List<Ms2ScanWithSpecificMass>[] _GetMs2Scans(MsDataFile myMSDataFile, string fullFilePath, CommonParameters commonParameters)
        {
            var msNScans = myMSDataFile.GetAllScansList().Where(x => x.MsnOrder > 1).ToArray();
            var ms2Scans = msNScans.Where(p => p.MsnOrder == 2).ToArray();
            var ms3Scans = msNScans.Where(p => p.MsnOrder == 3).ToArray();
            List<Ms2ScanWithSpecificMass>[] scansWithPrecursors = new List<Ms2ScanWithSpecificMass>[ms2Scans.Length];

            if (!ms2Scans.Any())
            {
                return scansWithPrecursors;
            }

            Parallel.ForEach(Partitioner.Create(0, ms2Scans.Length), new ParallelOptions { MaxDegreeOfParallelism = commonParameters.MaxThreadsToUsePerFile },
                (partitionRange, loopState) =>
                {
                    List<(double, int)> precursors = new List<(double, int)>();

                    for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                    {
                        if (GlobalVariables.StopLoops) { break; }

                        precursors.Clear();
                        MsDataScan ms2scan = ms2Scans[i];

                        if (ms2scan.OneBasedPrecursorScanNumber.HasValue)
                        {
                            MsDataScan precursorSpectrum = myMSDataFile.GetOneBasedScan(ms2scan.OneBasedPrecursorScanNumber.Value);

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
                                foreach (IsotopicEnvelope envelope in ms2scan.GetIsolatedMassesAndCharges(
                                    precursorSpectrum.MassSpectrum, 1,
                                    commonParameters.DeconvolutionMaxAssumedChargeState,
                                    commonParameters.DeconvolutionMassTolerance.Value,
                                    commonParameters.DeconvolutionIntensityRatio))
                                {
                                    double monoPeakMz = envelope.MonoisotopicMass.ToMz(envelope.Charge);
                                    precursors.Add((monoPeakMz, envelope.Charge));
                                }
                            }
                        }

                        if (commonParameters.UseProvidedPrecursorInfo && ms2scan.SelectedIonChargeStateGuess.HasValue)
                        {
                            int precursorCharge = ms2scan.SelectedIonChargeStateGuess.Value;

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
                        List<MsDataScan> ms2ChildScans = null;
                        List<MsDataScan> ms3ChildScans = null;
                        if (commonParameters.MS2ChildScanDissociationType != DissociationType.Unknown || commonParameters.MS3ChildScanDissociationType != DissociationType.Unknown)
                        {
                            ms3ChildScans = ms3Scans.Where(p => p.OneBasedPrecursorScanNumber == ms2scan.OneBasedScanNumber).ToList();

                            ms2ChildScans = ms2Scans.Where(p => p.OneBasedPrecursorScanNumber == ms2scan.OneBasedScanNumber ||
                            (p.OneBasedPrecursorScanNumber == ms2scan.OneBasedPrecursorScanNumber
                                && p.OneBasedScanNumber > ms2scan.OneBasedScanNumber
                                && Math.Abs(p.IsolationMz.Value - ms2scan.IsolationMz.Value) < 0.01)).ToList();
                        }

                        foreach (var precursor in precursors)
                        {
                            // assign precursor for this MS2 scan
                            var scan = new Ms2ScanWithSpecificMass(ms2scan, precursor.Item1,
                                precursor.Item2, fullFilePath, commonParameters, neutralExperimentalFragments);

                            // assign precursors for MS2 child scans
                            if (ms2ChildScans != null)
                            {
                                foreach (var ms2ChildScan in ms2ChildScans)
                                {
                                    IsotopicEnvelope[] childNeutralExperimentalFragments = null;

                                    if (commonParameters.MS2ChildScanDissociationType != DissociationType.LowCID)
                                    {
                                        childNeutralExperimentalFragments = Ms2ScanWithSpecificMass.GetNeutralExperimentalFragments(ms2ChildScan, commonParameters);
                                    }
                                    var theChildScan = new Ms2ScanWithSpecificMass(ms2ChildScan, precursor.Item1,
                                        precursor.Item2, fullFilePath, commonParameters, childNeutralExperimentalFragments);
                                    scan.ChildScans.Add(theChildScan);
                                }
                            }

                            // assign precursors for MS3 child scans
                            if (ms3ChildScans != null)
                            {
                                foreach (var ms3ChildScan in ms3ChildScans)
                                {
                                    int precursorCharge = 1;
                                    double precursorMz = 0;
                                    var precursorSpectrum = ms2scan;

                                    //In current situation, do we need to perform the following function. 
                                    //In some weird data, the MS3 scan has mis-leading precursor mass. 
                                    //MS3 scan is low res in most of the situation, and the matched ions are not scored in a good way.
                                    //{
                                    //    ms3ChildScan.RefineSelectedMzAndIntensity(precursorSpectrum.MassSpectrum);
                                    //    ms3ChildScan.ComputeMonoisotopicPeakIntensity(precursorSpectrum.MassSpectrum);
                                    //}

                                    if (ms3ChildScan.SelectedIonMonoisotopicGuessMz.HasValue)
                                    {
                                        precursorMz = ms3ChildScan.SelectedIonMonoisotopicGuessMz.Value;
                                    }
                                    else if (ms3ChildScan.SelectedIonMZ.HasValue)
                                    {
                                        precursorMz = ms3ChildScan.SelectedIonMZ.Value;
                                    }

                                    if (ms3ChildScan.SelectedIonChargeStateGuess.HasValue)
                                    {
                                        precursorCharge = ms3ChildScan.SelectedIonChargeStateGuess.Value;
                                    }

                                    IsotopicEnvelope[] childNeutralExperimentalFragments = null;

                                    if (commonParameters.MS3ChildScanDissociationType != DissociationType.LowCID)
                                    {
                                        childNeutralExperimentalFragments = Ms2ScanWithSpecificMass.GetNeutralExperimentalFragments(ms3ChildScan, commonParameters);
                                    }
                                    var theChildScan = new Ms2ScanWithSpecificMass(ms3ChildScan, precursorMz,
                                        precursorCharge, fullFilePath, commonParameters, childNeutralExperimentalFragments);

                                    scan.ChildScans.Add(theChildScan);
                                }
                            }

                            scansWithPrecursors[i].Add(scan);
                        }
                    }
                });

            return scansWithPrecursors;
        }

        public static IEnumerable<Ms2ScanWithSpecificMass> GetMs2Scans(MsDataFile myMSDataFile, string fullFilePath, CommonParameters commonParameters)
        {
            var scansWithPrecursors = _GetMs2Scans(myMSDataFile, fullFilePath, commonParameters);

            if (scansWithPrecursors.Length == 0)
            {
                return new List<Ms2ScanWithSpecificMass>();
            }

            var childScanNumbers = new HashSet<int>(scansWithPrecursors.SelectMany(p => p.SelectMany(v => v.ChildScans.Select(x => x.OneBasedScanNumber))));
            var parentScans = scansWithPrecursors.Where(p => p.Any() && !childScanNumbers.Contains(p.First().OneBasedScanNumber))
                .SelectMany(v => v)
                .OrderBy(p => p.OneBasedScanNumber)
                .ToArray();

            // XCorr pre-processing for low-res data. this is here because the parent/child scans may have different
            // resolutions, so this pre-processing must take place after the parent/child scans have been determined
            if (commonParameters.DissociationType == DissociationType.LowCID || commonParameters.MS2ChildScanDissociationType == DissociationType.LowCID || commonParameters.MS3ChildScanDissociationType == DissociationType.LowCID)
            {
                Parallel.ForEach(Partitioner.Create(0, parentScans.Length), new ParallelOptions { MaxDegreeOfParallelism = commonParameters.MaxThreadsToUsePerFile },
                    (partitionRange, loopState) =>
                    {
                        for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                        {
                            if (GlobalVariables.StopLoops) { break; }

                            var parentScan = parentScans[i];

                            if (commonParameters.DissociationType == DissociationType.LowCID && !parentScan.TheScan.MassSpectrum.XcorrProcessed)
                            {
                                lock (parentScan.TheScan)
                                {
                                    if (!parentScan.TheScan.MassSpectrum.XcorrProcessed)
                                    {
                                        parentScan.TheScan.MassSpectrum.XCorrPrePreprocessing(0, 1969, parentScan.TheScan.IsolationMz.Value);
                                    }
                                }
                            }

                            foreach (var childScan in parentScan.ChildScans)
                            {
                                if (((childScan.TheScan.MsnOrder == 2 && commonParameters.MS2ChildScanDissociationType == DissociationType.LowCID) ||
                                (childScan.TheScan.MsnOrder == 3 && commonParameters.MS3ChildScanDissociationType == DissociationType.LowCID))
                                && !childScan.TheScan.MassSpectrum.XcorrProcessed)
                                {
                                    lock (childScan.TheScan)
                                    {
                                        if (!childScan.TheScan.MassSpectrum.XcorrProcessed)
                                        {
                                            childScan.TheScan.MassSpectrum.XCorrPrePreprocessing(0, 1969, childScan.TheScan.IsolationMz.Value);
                                        }
                                    }
                                }
                            }
                        }
                    });
            }

            return parentScans;
        }

        public static List<Ms2ScanWithSpecificMass> GetMs2ScansWrapByScanNum(MsDataFile myMSDataFile, string fullFilePath, CommonParameters commonParameters, out List<List<(double, int, double)>> precursors)
        {
            var scansWithPrecursors = _GetMs2Scans(myMSDataFile, fullFilePath, commonParameters);

            var parentScans = new List<Ms2ScanWithSpecificMass>();
            precursors = new List<List<(double, int, double)>>();

            if (scansWithPrecursors.Length == 0)
            {
                return parentScans;
            }

            var childScanNumbers = new HashSet<int>(scansWithPrecursors.SelectMany(p => p.SelectMany(v => v.ChildScans.Select(x => x.OneBasedScanNumber))));
            //var parentScans = scansWithPrecursors.Where(p => p.Any() && !childScanNumbers.Contains(p.First().OneBasedScanNumber)).Select(p=>p.First()).ToArray();

            
            for (int i = 0; i < scansWithPrecursors.Length; i++)
            {
                if (scansWithPrecursors[i].Any() && !childScanNumbers.Contains(scansWithPrecursors[i].First().OneBasedScanNumber))
                {
                    parentScans.Add(scansWithPrecursors[i].First());
                    precursors.Add(scansWithPrecursors[i].Select(p => (p.PrecursorMass, p.PrecursorCharge, p.PrecursorMonoisotopicPeakMz)).ToList());
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
            string separationType = fileSpecificParams.SeparationType ?? commonParams.SeparationType;

            CommonParameters returnParams = new CommonParameters(
                dissociationType: dissociationType,
                precursorMassTolerance: precursorMassTolerance,
                productMassTolerance: productMassTolerance,
                digestionParams: fileSpecificDigestionParams,
                separationType: separationType,

                //NEED THESE OR THEY'LL BE OVERWRITTEN
                ms2childScanDissociationType: commonParams.MS2ChildScanDissociationType,
                ms3childScanDissociationType: commonParams.MS3ChildScanDissociationType,
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
            this.OutputFolder = output_folder;
            DetermineAnalyteType(CommonParameters);
            StartingSingleTask(displayName);

            var tomlFileName = Path.Combine(Directory.GetParent(output_folder).ToString(), "Task Settings", displayName + "config.toml");
            Toml.WriteFile(this, tomlFileName, tomlConfig);
            FinishedWritingFile(tomlFileName, new List<string> { displayName });

            FileSpecificParameters = new List<(string FileName, CommonParameters Parameters)>();

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
                        try
                        {
                            TomlTable fileSpecificSettings = Toml.ReadFile(fileSpecificTomlPath, tomlConfig);
                            fileSettingsList[i] = new FileSpecificParameters(fileSpecificSettings);
                            FileSpecificParameters.Add((currentRawDataFilepathList[i], SetAllFileSpecificCommonParams(CommonParameters, fileSettingsList[i])));
                        }
                        catch (MetaMorpheusException e)
                        {
                            //file - specific toml has already been validated in the GUI when the spectra files were added, so...
                            // probably the only time you can get here is if the user modifies the file-specific parameter file in the middle of a run...
                            Warn("Problem parsing the file-specific toml " + Path.GetFileName(fileSpecificTomlPath) + "; " + e.Message + "; is the toml from an older version of MetaMorpheus?");
                        }
                        catch (KeyNotFoundException e)
                        {
                            Warn("Problem parsing the file-specific toml " + Path.GetFileName(fileSpecificTomlPath) + "; " + e.Message + "; please update the proteases.tsv file and restart MetaMorpheus to use this file-specific toml.");
                        }
                    }
                    else // just used common parameters for file specific.
                    {
                        FileSpecificParameters.Add((currentRawDataFilepathList[i], CommonParameters));
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
                var proseFilePath = Path.Combine(output_folder, "AutoGeneratedManuscriptProse.txt");
                using (StreamWriter file = new StreamWriter(proseFilePath))
                {
                    file.WriteLine("The data analysis was performed using MetaMorpheus version " + GlobalVariables.MetaMorpheusVersion + ", available at " + "https://github.com/smith-chem-wisc/MetaMorpheus.");
                    file.Write(ProseCreatedWhileRunning.ToString());
                    file.WriteLine(SystemInfo.SystemProse().Replace(Environment.NewLine, "") + " ");
                    file.WriteLine("The total time to perform the " + TaskType + " task on " + currentRawDataFilepathList.Count + " spectra file(s) was " + String.Format("{0:0.00}", MyTaskResults.Time.TotalMinutes) + " minutes.");
                    file.WriteLine();
                    file.WriteLine("Published works using MetaMorpheus software are encouraged to cite: Solntsev, S. K.; Shortreed, M. R.; Frey, B. L.; Smith, L. M. Enhanced Global Post-translational Modification Discovery with MetaMorpheus. Journal of Proteome Research. 2018, 17 (5), 1844-1851.");

                    file.WriteLine();
                    file.WriteLine("Spectra files: ");
                    file.WriteLine(string.Join(Environment.NewLine, currentRawDataFilepathList.Select(b => '\t' + b)));
                    file.WriteLine("Databases:");

                    foreach (var proteinDb in currentProteinDbFilenameList)
                    {
                        if (proteinDb.IsContaminant)
                        {
                            file.Write(string.Join(Environment.NewLine, '\t' + "Contaminants " + proteinDb.FilePath + " Downloaded on: " + File.GetCreationTime(proteinDb.FilePath).ToString()));
                        }
                        else
                        {
                            file.Write(string.Join(Environment.NewLine, '\t' + proteinDb.FilePath + " Downloaded on: " + File.GetCreationTime(proteinDb.FilePath).ToString()));
                        }
                    }
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
            foreach (var db in dbFilenameList.Where(p => !p.IsSpectralLibrary))
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

        protected SpectralLibrary LoadSpectralLibraries(string taskId, List<DbForTask> dbFilenameList)
        {
            Status("Loading spectral libraries...", new List<string> { taskId });

            var paths = dbFilenameList.Where(p => p.IsSpectralLibrary).Select(p => p.FilePath).ToList();

            if (!paths.Any())
            {
                return null;
            }

            var lib = new SpectralLibrary(paths);

            Status("Done loading spectral libraries", new List<string> { taskId });
            return lib;
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
                proteinList = ProteinDbLoader.LoadProteinFasta(fileName, generateTargets, decoyType, isContaminant, out dbErrors,
                    ProteinDbLoader.UniprotAccessionRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotGeneNameRegex,
                    ProteinDbLoader.UniprotOrganismRegex, commonParameters.MaxThreadsToUsePerFile);
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

            // get digestion info from file
            var storedDigestParams = GetDigestionParamsFromFile(Path.Combine(Path.GetDirectoryName(peptideIndexFileName), "DigestionParameters.toml"));

            // get non-serialized information for the peptides (proteins, mod info)
            foreach (var peptide in peptideIndex)
            {
                peptide.SetNonSerializedPeptideInfo(GlobalVariables.AllModsKnownDictionary, proteinDictionary, storedDigestParams);
            }

            return peptideIndex;
        }

        private static DigestionParams GetDigestionParamsFromFile(string path)
        {
            var digestionParams = Toml.ReadFile<DigestionParams>(path, MetaMorpheusTask.tomlConfig);

            return digestionParams;
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

            Toml.WriteFile(indexEngine.CommonParameters.DigestionParams, Path.Combine(Path.GetDirectoryName(fileName), "DigestionParameters.toml"), tomlConfig);
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
            bool successfullyReadIndices = false;
            string pathToFolderWithIndices = GetExistingFolderWithIndices(indexEngine, dbFilenameList);

            if (pathToFolderWithIndices != null) //if indexes exist
            {
                try
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

                    successfullyReadIndices = true;
                }
                catch
                {
                    // could put something here... this basically is just to prevent a crash if the index was unable to be read.

                    // if the old index couldn't be read, a new one will be generated.

                    // an old index may not be able to be read because of information required by new versions of MetaMorpheus
                    // that wasn't written by old versions.
                }
            }

            if (!successfullyReadIndices) //if we didn't find indexes with the same params
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
                //ShrinkPeptideIndex(peptideIndex);
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
        }

        public void GenerateIndexes_PeptideOnly(IndexingEngine indexEngine, List<DbForTask> dbFilenameList, ref List<PeptideWithSetModifications> peptideIndex, ref List<int>[] precursorIndex, List<Protein> allKnownProteins, string taskId)
        {
            bool successfullyReadIndices = false;
            string pathToFolderWithIndices = GetExistingFolderWithIndices(indexEngine, dbFilenameList);

            if (pathToFolderWithIndices != null) //if indexes exist
            {
                try
                {
                    Status("Reading peptide index...", new List<string> { taskId });
                    peptideIndex = ReadPeptideIndex(Path.Combine(pathToFolderWithIndices, PeptideIndexFileName), allKnownProteins);

                    if (indexEngine.GeneratePrecursorIndex)
                    {
                        Status("Reading precursor index...", new List<string> { taskId });
                        precursorIndex = ReadFragmentIndex(Path.Combine(pathToFolderWithIndices, PrecursorIndexFileName));
                    }

                    successfullyReadIndices = true;
                }
                catch
                {
                    // could put something here... this basically is just to prevent a crash if the index was unable to be read.

                    // if the old index couldn't be read, a new one will be generated.

                    // an old index may not be able to be read because of information required by new versions of MetaMorpheus
                    // that wasn't written by old versions.
                }
            }

            if (!successfullyReadIndices) //if we didn't find indexes with the same params
            {
                // This is for the second round search, so successfullyReadIndices must be true, otherwise there are problems.
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

        public static void DetermineAnalyteType(CommonParameters commonParameters)
        {
            // changes the name of the analytes from "peptide" to "proteoform" if the protease is set to top-down

            // TODO: note that this will not function well if the user is using file-specific settings, but it's assumed
            // that bottom-up and top-down data is not being searched in the same task

            if (commonParameters != null
                && commonParameters.DigestionParams != null
                && commonParameters.DigestionParams.Protease != null
                && commonParameters.DigestionParams.Protease.Name == "top-down")
            {
                GlobalVariables.AnalyteType = "Proteoform";
            }
            else
            {
                GlobalVariables.AnalyteType = "Peptide";
            }
        }

        /// <summary>
        /// Handle ambiguity when two theoretical proteins in the
        /// search space have the same accession number.
        /// The accession must be unique for indexed searches
        ///
        /// RemoveAll is important because it references the location in memory, not the Equals
        /// </summary>
        /// <param name="proteins"></param>
        /// <param name="tcAmbiguity"></param>
        protected static void SanitizeProteinDatabase(List<Protein> proteins, TargetContaminantAmbiguity tcAmbiguity)
        {
            foreach (var accessionGroup in proteins.GroupBy(p => p.Accession))
            {
                if (accessionGroup.Count() != 1) //if multiple proteins with the same accession
                {
                    List<Protein> proteinsWithThisAccession = accessionGroup.OrderBy(p => p.OneBasedPossibleLocalizedModifications.Count).ThenBy(p => p.ProteolysisProducts.Count()).ToList();
                    List<Protein> proteinsToRemove = new List<Protein>();
                    if (tcAmbiguity == TargetContaminantAmbiguity.RenameProtein)
                    {
                        int proteinNumber = 1;
                        Warn("The protein '" + accessionGroup.Key + "' has multiple entries. Protein accessions must be unique. Protein " + accessionGroup.Key + " was renamed.");
                        foreach (var originalProtein in proteinsWithThisAccession)
                        {
                            //accession is private and there's no clone method, so we need to make a whole new protein... TODO: put this in mzlib
                            //use PROTEIN_D1 instead of PROTEIN_1 so it doesn't look like an isoform (D for Duplicate)
                            var renamedProtein = new Protein(originalProtein.BaseSequence, originalProtein + "_D" + proteinNumber.ToString(), originalProtein.Organism,
                                originalProtein.GeneNames.ToList(), originalProtein.OneBasedPossibleLocalizedModifications, originalProtein.ProteolysisProducts.ToList(), originalProtein.Name, originalProtein.FullName,
                                originalProtein.IsDecoy, originalProtein.IsContaminant, originalProtein.DatabaseReferences.ToList(), originalProtein.SequenceVariations.ToList(), originalProtein.AppliedSequenceVariations,
                                originalProtein.SampleNameForVariants, originalProtein.DisulfideBonds.ToList(), originalProtein.SpliceSites.ToList(), originalProtein.DatabaseFilePath);
                            proteins.Add(renamedProtein);
                            proteins.RemoveAll(p => p == originalProtein);
                            proteinNumber++;
                        }
                    }
                    else //if (tcAmbiguity == TargetContaminantAmbiguity.RemoveContaminant || tcAmbiguity == TargetContaminantAmbiguity.RemoveTarget)
                    {
                        if (tcAmbiguity == TargetContaminantAmbiguity.RemoveContaminant)
                        {
                            // remove contaminants
                            proteinsToRemove = proteinsWithThisAccession.Where(p => p.IsContaminant).ToList();
                            if (proteinsToRemove.Any())
                            {
                                Warn("The protein '" + accessionGroup.Key + "' has multiple entries. Protein accessions must be unique. Contaminant protein " + accessionGroup.Key + " was ignored.");
                            }
                        }
                        else //if (tcAmbiguity == TargetContaminantAmbiguity.RemoveTarget)
                        {
                            // remove targets
                            proteinsToRemove = proteinsWithThisAccession.Where(p => !p.IsDecoy && !p.IsContaminant).ToList();
                            if (proteinsToRemove.Any())
                            {
                                Warn("The protein '" + accessionGroup.Key + "' has multiple entries. Protein accessions must be unique. Target protein " + accessionGroup.Key + " was ignored.");
                            }
                        }

                        //remove the proteins specified above
                        foreach (var protein in proteinsToRemove)
                        {
                            if (proteinsWithThisAccession.Count > 1)
                            {
                                proteins.RemoveAll(p => p == protein);
                                proteinsWithThisAccession.RemoveAll(p => p == protein);
                            }
                        }

                        // most ambiguity should be handled by now, but for edge cases and decoys:
                        // remove proteins so that only 1 protein with this accession remains
                        for (int i = 0; i < proteinsWithThisAccession.Count - 1; i++) //-1 to keep the last one (most mods)
                        {
                            proteins.RemoveAll(p => p == proteinsWithThisAccession[i]);
                        }
                    }
                }
            }
        }
    }
}