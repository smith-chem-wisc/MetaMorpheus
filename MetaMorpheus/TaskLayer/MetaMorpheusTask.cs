using Chemistry;
using EngineLayer;
using EngineLayer.Indexing;
using MassSpectrometry;
using MzLibUtil;
using Nett;
using Omics;
using Omics.BioPolymer;
using Omics.Digestion;
using Omics.Modifications;
using Omics.SpectrumMatch;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Readers.SpectralLibrary;
using SpectralAveraging;
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
using UsefulProteomicsDatabases.Transcriptomics;
using Transcriptomics;
using Transcriptomics.Digestion;
using EngineLayer.Util;

namespace TaskLayer
{
    public enum MyTask
    {
        Search,
        Gptmd,
        Calibrate,
        XLSearch,
        GlycoSearch,
        Average
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
                    .FromToml(tmlString => GetModsFromString(tmlString.Value))))
            .ConfigureType<SpectraFileAveragingType>(type => type
                .WithConversionFor<TomlString>(convert => convert
                    .ToToml(custom => custom.ToString())
                    .FromToml(tmlString =>
                        tmlString.Value == "AverageDdaScansWithOverlap"
                            ? SpectraFileAveragingType.AverageDdaScans
                            : Enum.Parse<SpectraFileAveragingType>(tmlString.Value))))
            .ConfigureType<IDigestionParams>(type => type
                .WithConversionFor<TomlTable>(c => c
                    .FromToml(tmlTable =>
                        tmlTable.ContainsKey("Protease")
                            ? tmlTable.Get<DigestionParams>()
                            : tmlTable.Get<RnaDigestionParams>())))
            .ConfigureType<DigestionParams>(type => type
                .IgnoreProperty(p => p.DigestionAgent)
                .IgnoreProperty(p => p.MaxMods)
                .IgnoreProperty(p => p.MaxLength)
                .IgnoreProperty(p => p.MinLength))
            .ConfigureType<RnaDigestionParams>(type => type
                .IgnoreProperty(p => p.DigestionAgent))
            .ConfigureType<Rnase>(type => type
                .WithConversionFor<TomlString>(convert => convert
                    .ToToml(custom => custom.Name)
                    .FromToml(tmlString => RnaseDictionary.Dictionary[tmlString.Value])))
            // Switch on DeconvolutionParameters
            .ConfigureType<DeconvolutionParameters>(type => type
                .WithConversionFor<TomlTable>(c => c
                    .FromToml(tmlTable => tmlTable.Get<string>("DeconvolutionType") switch
                    {
                        "ClassicDeconvolution" => tmlTable.Get<ClassicDeconvolutionParameters>(),
                        "IsoDecDeconvolution" => tmlTable.Get<IsoDecDeconvolutionParameters>(),
                        _ => throw new MetaMorpheusException($"Toml Parsing Failure - Unknown Deconvolution Type: {tmlTable.Get<string>("DeconvolutionType")}")
                    })))
            // Ignore all properties that are not user settable, instantiate with defaults. If the toml differs, defaults will be overridden. 
            .ConfigureType<ClassicDeconvolutionParameters>(type => type
                .CreateInstance(() => new ClassicDeconvolutionParameters(1, 20, 4, 3))
                .IgnoreProperty(p => p.IntensityRatioLimit)
                .IgnoreProperty(p => p.DeconvolutionTolerancePpm))
            .ConfigureType<IsoDecDeconvolutionParameters>(type => type
                .CreateInstance(() => new IsoDecDeconvolutionParameters())
                .IgnoreProperty(p => p.Verbose)
                .IgnoreProperty(p => p.PeakWindow)
                .IgnoreProperty(p => p.PeakThreshold)
                .IgnoreProperty(p => p.MinPeaks)
                .IgnoreProperty(p => p.PlusOneIntWindow)
                .IgnoreProperty(p => p.MinScoreDiff)
                .IgnoreProperty(p => p.IsoLength)
                .IgnoreProperty(p => p.MassDiffC)
                .IgnoreProperty(p => p.MinusOneAreasZero)
                .IgnoreProperty(p => p.IsotopeThreshold)
                .IgnoreProperty(p => p.ZScoreThreshold))

            // Convert average residue models to simple strings instead of tables, Nett makes all objects tables by default
            // The base class AverageResidue is used for Toml Reading. The derived classes are used for toml writing. 
            .ConfigureType<AverageResidue>(type => type
                .WithConversionFor<TomlString>(convert => convert
                    .FromToml(tmlString =>
                        tmlString.Value switch
                        {
                            "Averagine" => new Averagine(),
                            "OxyriboAveragine" => new OxyriboAveragine(),
                            _ => throw new MetaMorpheusException($"Toml Parsing Failure - Unknown AverageResidueModel: {tmlString.Value}")
                        }
                    )))
            .ConfigureType<Averagine>(type => type
                .WithConversionFor<TomlString>(convert => convert
                    .ToToml(custom => custom.GetType().Name)))
            .ConfigureType<OxyriboAveragine>(type => type
                .WithConversionFor<TomlString>(convert => convert
                    .ToToml(custom => custom.GetType().Name)))
            .ConfigureType<List<IGptmdFilter>>(type => type
                .WithConversionFor<TomlString>(convert => convert
                    .ToToml(filters => string.Join("\t", filters.Select(f => f.GetType().Name)))
                    .FromToml(tmlString => tmlString.Value
                        .Split('\t', StringSplitOptions.RemoveEmptyEntries)
                        .Select(typeName =>
                            // Find the type in the current AppDomain by name
                            AppDomain.CurrentDomain.GetAssemblies()
                                .SelectMany(a => a.GetTypes())
                                .FirstOrDefault(t => t.Name == typeName && typeof(IGptmdFilter).IsAssignableFrom(t))
                        )
                        .Where(t => t != null)
                        .Select(t => Activator.CreateInstance(t) as IGptmdFilter)
                        .Where(f => f != null)
                        .ToList()
                    )
                )
            )
        );
       

        protected readonly StringBuilder ProseCreatedWhileRunning = new StringBuilder();

        [TomlIgnore]
        public string OutputFolder { get; private set; }

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
                    var precursorSet = new PrecursorSet(commonParameters.DeconvolutionMassTolerance);

                    for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                    {
                        if (GlobalVariables.StopLoops) { break; }

                        precursorSet.Clear();
                        MsDataScan ms2scan = ms2Scans[i];

                        if (ms2scan.OneBasedPrecursorScanNumber.HasValue)
                        {
                            MsDataScan precursorSpectrum = myMSDataFile.GetOneBasedScan(ms2scan.OneBasedPrecursorScanNumber.Value);

                            if (precursorSpectrum is null)
                                goto PrecursorFromScanHeader;

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
                                    precursorSpectrum.MassSpectrum, commonParameters.PrecursorDeconvolutionParameters))
                                {
                                    double? intensity = null;
                                    if (commonParameters.UseMostAbundantPrecursorIntensity) 
                                        intensity = envelope.Peaks.Max(p => p.intensity); 

                                    var fractionalIntensity = envelope.TotalIntensity /
                                          precursorSpectrum.MassSpectrum.YArray
                                          [
                                              precursorSpectrum.MassSpectrum.GetClosestPeakIndex(ms2scan.IsolationRange.Minimum)
                                              ..
                                              precursorSpectrum.MassSpectrum.GetClosestPeakIndex(ms2scan.IsolationRange.Maximum)
                                          ].Sum();

                                    precursorSet.Add(new(envelope, intensity, fractionalIntensity));
                                }
                            }
                        }

                         // If using precursor info from scan header and scan header has charge state.
                         // MsAlign uses this conditional to construct its precursors. 
                        PrecursorFromScanHeader:
                        if (commonParameters.UseProvidedPrecursorInfo && ms2scan.SelectedIonChargeStateGuess.HasValue && ms2scan.SelectedIonChargeStateGuess != 0) 
                        {
                            int precursorCharge = ms2scan.SelectedIonChargeStateGuess.Value;
                            double precursorIntensity = ms2scan.SelectedIonMonoisotopicGuessIntensity ?? ms2scan.SelectedIonIntensity ?? 1.0;
                            double precursorMz = ms2scan.SelectedIonMonoisotopicGuessMz ?? ms2scan.SelectedIonMZ.Value;

                            precursorSet.Add(new Precursor(precursorMz, precursorCharge, precursorIntensity, 1, null));
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

                        foreach (var precursor in precursorSet)
                        {
                            // assign precursor for this MS2 scan
                            var scan = new Ms2ScanWithSpecificMass(ms2scan, precursor.MonoisotopicPeakMz,
                                precursor.Charge, fullFilePath, commonParameters, neutralExperimentalFragments,
                                precursor.Intensity, precursor.EnvelopePeakCount, precursor.FractionalIntensity);

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
                                    var theChildScan = new Ms2ScanWithSpecificMass(ms2ChildScan, precursor.MonoisotopicPeakMz,
                                        precursor.Charge, fullFilePath, commonParameters, childNeutralExperimentalFragments,
                                        precursor.Intensity, precursor.EnvelopePeakCount, precursor.FractionalIntensity);
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
            int minPeptideLength = fileSpecificParams.MinPeptideLength ?? commonParams.DigestionParams.MinLength;
            int maxPeptideLength = fileSpecificParams.MaxPeptideLength ?? commonParams.DigestionParams.MaxLength;
            int maxMissedCleavages = fileSpecificParams.MaxMissedCleavages ?? commonParams.DigestionParams.MaxMissedCleavages;
            int maxModsForPeptide = fileSpecificParams.MaxModsForPeptide ?? commonParams.DigestionParams.MaxMods;

            // set file-specific digestion params based upon the type of digestion params
            IDigestionParams fileSpecificDigestionParams;
            switch (commonParams.DigestionParams)
            {
                case DigestionParams digestionParams:
                    fileSpecificDigestionParams = new DigestionParams(
                        protease: (fileSpecificParams.DigestionAgent ?? digestionParams.SpecificProtease).Name,
                        maxMissedCleavages: maxMissedCleavages, minPeptideLength: minPeptideLength,
                        maxPeptideLength: maxPeptideLength, maxModsForPeptides: maxModsForPeptide,
                        maxModificationIsoforms: digestionParams.MaxModificationIsoforms,
                        initiatorMethionineBehavior: digestionParams.InitiatorMethionineBehavior,
                        fragmentationTerminus: digestionParams.FragmentationTerminus,
                        searchModeType: digestionParams.SearchModeType);
                    break;
                case RnaDigestionParams:
                    fileSpecificDigestionParams = new RnaDigestionParams(
                        rnase: (fileSpecificParams.DigestionAgent ?? commonParams.DigestionParams.DigestionAgent).Name,
                        maxMissedCleavages: maxMissedCleavages, minLength: minPeptideLength,
                        maxLength: maxPeptideLength, maxMods: maxModsForPeptide,
                        maxModificationIsoforms: commonParams.DigestionParams.MaxModificationIsoforms,
                        fragmentationTerminus: commonParams.DigestionParams.FragmentationTerminus);
                    break;
                default:
                    throw new MetaMorpheusException(
                        $"Unrecognized digestion parameters of type {commonParams.DigestionParams.GetType().FullName} in MetaMorpheusTask.SetAllFileSpecificCommonParams");
            }

            // must be set in this manner as CommonParameters constructor will pull from this dictionary, then clear dictionary
            fileSpecificDigestionParams.ProductsFromDissociationType()[DissociationType.Custom] =
                fileSpecificParams.CustomIons ?? commonParams.CustomIons;

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
                qValueThreshold: commonParams.QValueThreshold,
                pepQValueThreshold: commonParams.PepQValueThreshold,
                scoreCutoff: commonParams.ScoreCutoff,
                numberOfPeaksToKeepPerWindow: commonParams.NumberOfPeaksToKeepPerWindow,
                minimumAllowedIntensityRatioToBasePeak: commonParams.MinimumAllowedIntensityRatioToBasePeak,
                windowWidthThomsons: commonParams.WindowWidthThomsons,
                trimMs1Peaks: commonParams.TrimMs1Peaks,
                trimMsMsPeaks: commonParams.TrimMsMsPeaks,
                normalizePeaksAccrossAllWindows: commonParams.NormalizePeaksAccrossAllWindows,
                deconvolutionMassTolerance: commonParams.DeconvolutionMassTolerance,
                maxThreadsToUsePerFile: commonParams.MaxThreadsToUsePerFile,
                listOfModsVariable: commonParams.ListOfModsVariable,
                listOfModsFixed: commonParams.ListOfModsFixed,
                taskDescriptor: commonParams.TaskDescriptor,
                assumeOrphanPeaksAreZ1Fragments: commonParams.AssumeOrphanPeaksAreZ1Fragments,
                maxHeterozygousVariants: commonParams.MaxHeterozygousVariants,
                minVariantDepth: commonParams.MinVariantDepth,
                addTruncations: commonParams.AddTruncations,
                precursorDeconParams: commonParams.PrecursorDeconvolutionParameters,
                productDeconParams: commonParams.ProductDeconvolutionParameters,
                useMostAbundantPrecursorIntensity: commonParams.UseMostAbundantPrecursorIntensity);

            return returnParams;
        }

        public MyTaskResults RunTask(string output_folder, List<DbForTask> currentProteinDbFilenameList, List<string> currentRawDataFilepathList, string displayName)
        {
            this.OutputFolder = output_folder;
            MetaMorpheusEngine.DetermineAnalyteType(CommonParameters);
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

            var proseFilePath = Path.Combine(output_folder, "AutoGeneratedManuscriptProse.txt");
            using (StreamWriter file = new StreamWriter(proseFilePath))
            {
                file.WriteLine("The data analysis was performed using MetaMorpheus version " + GlobalVariables.MetaMorpheusVersion + ", available at " + "https://github.com/smith-chem-wisc/MetaMorpheus.");
		        file.WriteLine();
                file.Write(ProseCreatedWhileRunning.ToString());
                file.WriteLine(SystemInfo.SystemProse().Replace(Environment.NewLine, "") + " ");
		        file.WriteLine();
                file.WriteLine("The total time to perform the " + TaskType + " task on " + currentRawDataFilepathList.Count + " spectra file(s) was " + String.Format("{0:0.00}", MyTaskResults.Time.TotalMinutes) + " minutes.");
                file.WriteLine();
                file.WriteLine("Published works using MetaMorpheus software are encouraged to cite the appropriate publications listed in the reference guide, found here: https://github.com/smith-chem-wisc/MetaMorpheus/blob/master/README.md.");
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
            MetaMorpheusEngine.FinishedSingleEngineHandler -= SingleEngineHandlerInTask;
            return MyTaskResults;
        }

        #region Database Loading

        protected List<IBioPolymer> LoadBioPolymers(string taskId, List<DbForTask> dbFilenameList, bool searchTarget, DecoyType decoyType, List<string> localizeableModificationTypes, CommonParameters commonParameters)
        {
            Status($"Loading {GlobalVariables.AnalyteType.GetBioPolymerLabel()}s...", new List<string> { taskId });
            int emptyEntries = 0;
            List<IBioPolymer> bioPolymerList = new();
            foreach (var db in dbFilenameList.Where(p => !p.IsSpectralLibrary))
            {
                if (GlobalVariables.AnalyteType == AnalyteType.Oligo)
                {
                    var dbOligoList = LoadOligoDb(db.FilePath, searchTarget, decoyType, localizeableModificationTypes, db.IsContaminant, out Dictionary<string, Modification> unknownModifications, out int emptyOligoEntriesForThisDb, commonParameters);
                    bioPolymerList = bioPolymerList.Concat(dbOligoList).ToList();
                    emptyEntries += emptyOligoEntriesForThisDb;
                }
                else
                {
                    var dbProteinList = LoadProteinDb(db.FilePath, searchTarget, decoyType, localizeableModificationTypes, db.IsContaminant, out Dictionary<string, Modification> unknownModifications, out int emptyProteinEntriesForThisDb, commonParameters);
                    bioPolymerList = bioPolymerList.Concat(dbProteinList).ToList();
                    emptyEntries += emptyProteinEntriesForThisDb;
                }
            }
            if (!bioPolymerList.Any())
            {
                Warn($"Warning: No {GlobalVariables.AnalyteType.GetBioPolymerLabel()} entries were found in the database");
            }
            else if (emptyEntries > 0)
            {
                Warn("Warning: " + emptyEntries + $" empty {GlobalVariables.AnalyteType.GetBioPolymerLabel()} entries ignored");
            }

            // We are not generating decoys, so just return the read in database
            if (!bioPolymerList.Any(p => p.IsDecoy))
            {
                Status($"Done loading {GlobalVariables.AnalyteType.GetBioPolymerLabel()}s", new List<string> { taskId });
                return bioPolymerList;
            }

            // Sanitize the decoys
            // TODO: Fix this so that it accounts for multi-protease searches. Currently, we only consider the first protease
            // when looking for target/decoy collisions
            HashSet<string> targetPeptideSequences = new();
            foreach(var bioPolymer in bioPolymerList.Where(p => !p.IsDecoy))
            {
                // When thinking about decoy collisions, we can ignore modifications
                foreach(var peptide in bioPolymer.Digest(commonParameters.DigestionParams, new List<Modification>(), new List<Modification>()))
                {
                    targetPeptideSequences.Add(peptide.BaseSequence);
                }
            }
            // Now, we iterate through the decoys and scramble the sequences that correspond to target peptides
            for(int i = 0; i < bioPolymerList.Count; i++)
            {
                if(bioPolymerList[i].IsDecoy)
                {
                    var peptidesToReplace = bioPolymerList[i]
                        .Digest(commonParameters.DigestionParams, new List<Modification>(), new List<Modification>())
                        .Select(p => p.BaseSequence)
                        .Where(targetPeptideSequences.Contains)
                        .ToList();
                    if(peptidesToReplace.Any())
                    {
                        bioPolymerList[i] = DecoySequenceValidator.ScrambleDecoyBioPolymer(bioPolymerList[i], commonParameters.DigestionParams, forbiddenSequences: targetPeptideSequences, peptidesToReplace);
                    }
                }
            }

            Status($"Done loading {GlobalVariables.AnalyteType.GetBioPolymerLabel()}s", new List<string> { taskId });
            return bioPolymerList;
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
                    ProteinDbLoader.UniprotOrganismRegex, commonParameters.MaxThreadsToUsePerFile, addTruncations: commonParameters.AddTruncations);
            }
            else
            {
                List<string> modTypesToExclude = GlobalVariables.AllModTypesKnown.Where(b => !localizeableModificationTypes.Contains(b)).ToList();
                proteinList = ProteinDbLoader.LoadProteinXML(fileName, generateTargets, decoyType, GlobalVariables.AllModsKnown, isContaminant, modTypesToExclude, out um, commonParameters.MaxThreadsToUsePerFile, commonParameters.MaxHeterozygousVariants, commonParameters.MinVariantDepth, addTruncations: commonParameters.AddTruncations);
            }

            emptyEntriesCount = proteinList.Count(p => p.BaseSequence.Length == 0);
            return proteinList.Where(p => p.BaseSequence.Length > 0).ToList();
        }

        protected void LoadModifications(string taskId, out List<Modification> variableModifications, out List<Modification> fixedModifications, out List<string> localizableModificationTypes)
        {
            // load modifications
            Status("Loading modifications...", taskId);
            switch (GlobalVariables.AnalyteType)
            {
                case AnalyteType.Oligo:
                    variableModifications = GlobalVariables.AllRnaModsKnown
                        .Where(b => CommonParameters.ListOfModsVariable.Contains((b.ModificationType, b.IdWithMotif)))
                        .ToList();
                    fixedModifications = GlobalVariables.AllRnaModsKnown
                        .Where(b => CommonParameters.ListOfModsFixed.Contains((b.ModificationType, b.IdWithMotif)))
                        .ToList();
                    localizableModificationTypes = GlobalVariables.AllRnaModTypesKnown.ToList();
                    break;

                case AnalyteType.Peptide:
                case AnalyteType.Proteoform:
                default:
                    variableModifications = GlobalVariables.AllModsKnown
                        .Where(b => CommonParameters.ListOfModsVariable.Contains((b.ModificationType, b.IdWithMotif)))
                        .ToList();
                    fixedModifications = GlobalVariables.AllModsKnown
                        .Where(b => CommonParameters.ListOfModsFixed.Contains((b.ModificationType, b.IdWithMotif)))
                        .ToList();
                    localizableModificationTypes = GlobalVariables.AllModTypesKnown.ToList();
                    break;
            }

            var recognizedVariable = variableModifications.Select(p => p.IdWithMotif);
            var recognizedFixed = fixedModifications.Select(p => p.IdWithMotif);
            var unknownMods = CommonParameters.ListOfModsVariable.Select(p => p.Item2).Except(recognizedVariable).ToList();
            unknownMods.AddRange(CommonParameters.ListOfModsFixed.Select(p => p.Item2).Except(recognizedFixed));
            foreach (var unrecognizedMod in unknownMods)
            {
                Warn("Unrecognized mod " + unrecognizedMod + "; are you using an old .toml?");
            }
        }

        protected List<RNA> LoadOligoDb(string fileName, bool generateTargets, DecoyType decoyType,
            List<string> localizeableModificationTypes, bool isContaminant,
            out Dictionary<string, Modification> unknownMods, out int emptyEntriesCount,
            CommonParameters commonParameters)
        {
            List<string> dbErrors = new List<string>();
            List<RNA> rnaList = new List<RNA>();

            string theExtension = Path.GetExtension(fileName).ToLowerInvariant();
            bool compressed = theExtension.EndsWith("gz"); // allows for .bgz and .tgz, too which are used on occasion
            theExtension = compressed ? Path.GetExtension(Path.GetFileNameWithoutExtension(fileName)).ToLowerInvariant() : theExtension;

            if (theExtension.Equals(".fasta") || theExtension.Equals(".fa"))
            {
                unknownMods = null;
                rnaList = RnaDbLoader.LoadRnaFasta(fileName, generateTargets, decoyType, isContaminant, out dbErrors);
            }
            else
            {
                List<string> modTypesToExclude = GlobalVariables.AllRnaModTypesKnown.Where(b => !localizeableModificationTypes.Contains(b)).ToList();
                rnaList = RnaDbLoader.LoadRnaXML(fileName, generateTargets, decoyType, isContaminant, GlobalVariables.AllRnaModsKnown, modTypesToExclude, out unknownMods, commonParameters.MaxThreadsToUsePerFile);
            }

            emptyEntriesCount = rnaList.Count(p => p.BaseSequence.Length == 0);
            return rnaList.Where(p => p.BaseSequence.Length > 0).ToList();
        }

        protected void WritePrunedDatabase(List<SpectralMatch> allSpectralMatches, List<IBioPolymer> bioPolymersToWrite, Dictionary<string, int> modificationsToWrite, List<DbForTask> inputDatabases, string outputDirectory, string taskId)
        {
            Status("Writing Pruned Database...", new List<string> { taskId });
            HashSet<Modification> modificationsToWriteIfBoth = new HashSet<Modification>();
            HashSet<Modification> modificationsToWriteIfInDatabase = new HashSet<Modification>();
            HashSet<Modification> modificationsToWriteIfObserved = new HashSet<Modification>();

            var filteredPsms = FilteredPsms.Filter(allSpectralMatches,
                CommonParameters,
                includeDecoys: false,
                includeContaminants: true,
                includeAmbiguous: false,
                includeHighQValuePsms: false);

            var proteinToConfidentBaseSequences = new Dictionary<IBioPolymer, List<IBioPolymerWithSetMods>>();

            // associate all confident PSMs with all possible proteins they could be digest products of (before or after parsimony)
            foreach (SpectralMatch psm in filteredPsms)
            {
                var myPepsWithSetMods = psm.BestMatchingBioPolymersWithSetMods.Select(p => p.SpecificBioPolymer);

                foreach (IBioPolymerWithSetMods peptide in myPepsWithSetMods)
                {
                    if (proteinToConfidentBaseSequences.TryGetValue(peptide.Parent.ConsensusVariant, out var myPepList))
                    {
                        myPepList.Add(peptide);
                    }
                    else
                    {
                        proteinToConfidentBaseSequences.Add(peptide.Parent.ConsensusVariant, new List<IBioPolymerWithSetMods> { peptide });
                    }
                }
            }

            // Add user mod selection behavours to Pruned DB
            foreach (var modType in modificationsToWrite)
            {
                foreach (Modification mod in GlobalVariables.AllModsKnown.Where(b => b.ModificationType.Equals(modType.Key)))
                {
                    if (modType.Value == 1) // Write if observed and in database
                    {
                        modificationsToWriteIfBoth.Add(mod);
                    }
                    if (modType.Value == 2) // Write if in database
                    {
                        modificationsToWriteIfInDatabase.Add(mod);
                    }
                    if (modType.Value == 3) // Write if observed
                    {
                        modificationsToWriteIfObserved.Add(mod);
                    }
                }
            }

            //generates dictionary of proteins with only localized modifications
            var originalModPsms = FilteredPsms.Filter(filteredPsms,
                CommonParameters,
                includeDecoys: false,
                includeContaminants: true,
                includeAmbiguous: false,
                includeAmbiguousMods: false,
                includeHighQValuePsms: false);


            var proteinToConfidentModifiedSequences = new Dictionary<IBioPolymer, List<IBioPolymerWithSetMods>>();

            HashSet<string> modPsmsFullSeq = originalModPsms.Select(p => p.FullSequence).ToHashSet();
            HashSet<string> originalModPsmsFullSeq = originalModPsms.Select(p => p.FullSequence).ToHashSet();
            modPsmsFullSeq.ExceptWith(originalModPsmsFullSeq);

            foreach (SpectralMatch psm in originalModPsms)
            {
                var myPepsWithSetMods = psm.BestMatchingBioPolymersWithSetMods.Select(p => p.SpecificBioPolymer);

                foreach (IBioPolymerWithSetMods peptide in myPepsWithSetMods)
                {
                    if (proteinToConfidentModifiedSequences.TryGetValue(peptide.Parent.ConsensusVariant, out var myPepList))
                    {
                        myPepList.Add(peptide);
                    }
                    else
                    {
                        proteinToConfidentModifiedSequences.Add(peptide.Parent.ConsensusVariant, new List<IBioPolymerWithSetMods> { peptide });
                    }
                }
            }

            Dictionary<IBioPolymer, Dictionary<int, List<Modification>>> proteinsOriginalModifications = new Dictionary<IBioPolymer, Dictionary<int, List<Modification>>>();
            Dictionary<SequenceVariation, Dictionary<int, List<Modification>>> originalSequenceVariantModifications = new Dictionary<SequenceVariation, Dictionary<int, List<Modification>>>();

            // mods included in pruned database will only be confidently localized mods (peptide's FullSequence != null)
            foreach (var nonVariantProtein in bioPolymersToWrite.Select(p => p.ConsensusVariant).Distinct())
            {
                if (!nonVariantProtein.IsDecoy)
                {
                    proteinToConfidentModifiedSequences.TryGetValue(nonVariantProtein, out var psms);
                    HashSet<(int, Modification, SequenceVariation)> modsObservedOnThisProtein = new HashSet<(int, Modification, SequenceVariation)>(); // sequence variant is null if mod is not on a variant
                    foreach (IBioPolymerWithSetMods psm in psms ?? new List<IBioPolymerWithSetMods>())
                    {
                        foreach (var idxModKV in psm.AllModsOneIsNterminus)
                        {
                            int proteinIdx = GetOneBasedIndexInProtein(idxModKV.Key, psm);
                            SequenceVariation relevantVariant = psm.Parent.AppliedSequenceVariations.FirstOrDefault(sv => VariantApplication.IsSequenceVariantModification(sv, proteinIdx));
                            SequenceVariation unappliedVariant =
                                relevantVariant == null ? null : // it's not a sequence variant mod
                                    psm.Parent.SequenceVariations.FirstOrDefault(sv => sv.Description != null && sv.Description.Equals(relevantVariant.Description));
                            modsObservedOnThisProtein.Add((VariantApplication.RestoreModificationIndex(psm.Parent, proteinIdx), idxModKV.Value, unappliedVariant));
                        }
                    }

                    IDictionary<(SequenceVariation, int), List<Modification>> modsToWrite = new Dictionary<(SequenceVariation, int), List<Modification>>();

                    //Add if observed (regardless if in database)
                    foreach (var observedMod in modsObservedOnThisProtein)
                    {
                        var tempMod = observedMod.Item2;

                        if (modificationsToWriteIfObserved.Contains(tempMod))
                        {
                            var svIdxKey = (observedMod.Item3, observedMod.Item1);
                            if (!modsToWrite.ContainsKey(svIdxKey))
                            {
                                modsToWrite.Add(svIdxKey, new List<Modification> { observedMod.Item2 });
                            }
                            else
                            {
                                modsToWrite[svIdxKey].Add(observedMod.Item2);
                            }
                        }
                    }

                    // Add modification if in database (two cases: always or if observed)
                    foreach (var modkv in nonVariantProtein.OneBasedPossibleLocalizedModifications)
                    {
                        foreach (var mod in modkv.Value)
                        {
                            //Add if always In Database or if was observed and in database and not set to not include
                            if (modificationsToWriteIfInDatabase.Contains(mod) ||
                                (modificationsToWriteIfBoth.Contains(mod) && modsObservedOnThisProtein.Contains((modkv.Key, mod, null))))
                            {
                                if (!modsToWrite.ContainsKey((null, modkv.Key)))
                                {
                                    modsToWrite.Add((null, modkv.Key), new List<Modification> { mod });
                                }
                                else
                                {
                                    modsToWrite[(null, modkv.Key)].Add(mod);
                                }
                            }
                        }
                    }

                    //TODO add unit test here
                    // Add variant modification if in database (two cases: always or if observed)
                    foreach (SequenceVariation sv in nonVariantProtein.SequenceVariations)
                    {
                        foreach (var modkv in sv.OneBasedModifications)
                        {
                            foreach (var mod in modkv.Value)
                            {
                                //Add if always In Database or if was observed and in database and not set to not include
                                if (modificationsToWriteIfInDatabase.Contains(mod) ||
                                    (modificationsToWriteIfBoth.Contains(mod) && modsObservedOnThisProtein.Contains((modkv.Key, mod, sv))))
                                {
                                    if (!modsToWrite.ContainsKey((sv, modkv.Key)))
                                    {
                                        modsToWrite.Add((sv, modkv.Key), new List<Modification> { mod });
                                    }
                                    else
                                    {
                                        modsToWrite[(sv, modkv.Key)].Add(mod);
                                    }
                                }
                            }
                        }
                    }

                    var oldMods = nonVariantProtein.OneBasedPossibleLocalizedModifications.ToDictionary(p => p.Key, v => v.Value);
                    if (proteinsOriginalModifications.ContainsKey(nonVariantProtein.ConsensusVariant))
                    {
                        foreach (var entry in oldMods)
                        {
                            if (proteinsOriginalModifications[nonVariantProtein.ConsensusVariant].ContainsKey(entry.Key))
                            {
                                proteinsOriginalModifications[nonVariantProtein.ConsensusVariant][entry.Key].AddRange(entry.Value);
                            }
                            else
                            {
                                proteinsOriginalModifications[nonVariantProtein.ConsensusVariant].Add(entry.Key, entry.Value);
                            }
                        }
                    }
                    else
                    {
                        proteinsOriginalModifications.Add(nonVariantProtein.ConsensusVariant, oldMods);
                    }

                    // adds confidently localized and identified mods
                    nonVariantProtein.OneBasedPossibleLocalizedModifications.Clear();
                    foreach (var kvp in modsToWrite.Where(kv => kv.Key.Item1 == null))
                    {
                        nonVariantProtein.OneBasedPossibleLocalizedModifications.Add(kvp.Key.Item2, kvp.Value);
                    }
                    foreach (var sv in nonVariantProtein.SequenceVariations)
                    {
                        var oldVariantModifications = sv.OneBasedModifications.ToDictionary(p => p.Key, v => v.Value);
                        if (originalSequenceVariantModifications.ContainsKey(sv))
                        {
                            foreach (var entry in oldVariantModifications)
                            {
                                if (originalSequenceVariantModifications[sv].ContainsKey(entry.Key))
                                {
                                    originalSequenceVariantModifications[sv][entry.Key].AddRange(entry.Value);
                                }
                                else
                                {
                                    originalSequenceVariantModifications[sv].Add(entry.Key, entry.Value);
                                }
                            }
                        }
                        else
                        {
                            originalSequenceVariantModifications.Add(sv, oldVariantModifications);
                        }

                        sv.OneBasedModifications.Clear();
                        foreach (var kvp in modsToWrite.Where(kv => kv.Key.Item1 != null && kv.Key.Item1.Equals(sv)))
                        {
                            sv.OneBasedModifications.Add(kvp.Key.Item2, kvp.Value);
                        }
                    }
                }
            }

            //writes all proteins
            if (inputDatabases.Any(b => !b.IsContaminant))
            {
                string outputXMLdbFullName = Path.Combine(outputDirectory, string.Join("-", inputDatabases.Where(b => !b.IsContaminant).Select(b => Path.GetFileNameWithoutExtension(b.FilePath))) + "pruned.xml");
                ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), bioPolymersToWrite.Select(p => p.ConsensusVariant).Where(b => !b.IsDecoy && !b.IsContaminant).ToList(), outputXMLdbFullName);
                FinishedWritingFile(outputXMLdbFullName, new List<string> { taskId });
            }
            if (inputDatabases.Any(b => b.IsContaminant))
            {
                string outputXMLdbFullNameContaminants = Path.Combine(outputDirectory, string.Join("-", inputDatabases.Where(b => b.IsContaminant).Select(b => Path.GetFileNameWithoutExtension(b.FilePath))) + "pruned.xml");
                ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), bioPolymersToWrite.Select(p => p.ConsensusVariant).Where(b => !b.IsDecoy && b.IsContaminant).ToList(), outputXMLdbFullNameContaminants);
                FinishedWritingFile(outputXMLdbFullNameContaminants, new List<string> { taskId });
            }

            //writes only detected proteins
            if (inputDatabases.Any(b => !b.IsContaminant))
            {
                string outputXMLdbFullName = Path.Combine(outputDirectory, string.Join("-", inputDatabases.Where(b => !b.IsContaminant).Select(b => Path.GetFileNameWithoutExtension(b.FilePath))) + "proteinPruned.xml");
                ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), proteinToConfidentBaseSequences.Keys.Where(b => !b.IsDecoy && !b.IsContaminant).ToList(), outputXMLdbFullName);
                FinishedWritingFile(outputXMLdbFullName, new List<string> { taskId });
            }
            if (inputDatabases.Any(b => b.IsContaminant))
            {
                string outputXMLdbFullNameContaminants = Path.Combine(outputDirectory, string.Join("-", inputDatabases.Where(b => b.IsContaminant).Select(b => Path.GetFileNameWithoutExtension(b.FilePath))) + "proteinPruned.xml");
                ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), proteinToConfidentBaseSequences.Keys.Where(b => !b.IsDecoy && b.IsContaminant).ToList(), outputXMLdbFullNameContaminants);
                FinishedWritingFile(outputXMLdbFullNameContaminants, new List<string> { taskId });
            }

            foreach (var nonVariantProtein in bioPolymersToWrite.Select(p => p.ConsensusVariant).Distinct())
            {
                if (!nonVariantProtein.IsDecoy)
                {
                    nonVariantProtein.OneBasedPossibleLocalizedModifications.Clear();
                    foreach (var originalMod in proteinsOriginalModifications[nonVariantProtein.ConsensusVariant])
                    {
                        nonVariantProtein.OneBasedPossibleLocalizedModifications.Add(originalMod.Key, originalMod.Value);
                    }
                    foreach (var sv in nonVariantProtein.SequenceVariations)
                    {
                        sv.OneBasedModifications.Clear();
                        foreach (var originalVariantMods in originalSequenceVariantModifications[sv])
                        {
                            sv.OneBasedModifications.Add(originalVariantMods.Key, originalVariantMods.Value);
                        }
                    }
                }
            }
        }

        // Pruned Database Helper
        private static int GetOneBasedIndexInProtein(int oneIsNterminus, IBioPolymerWithSetMods peptideWithSetModifications)
        {
            if (oneIsNterminus == 1)
            {
                return peptideWithSetModifications.OneBasedStartResidue;
            }
            if (oneIsNterminus == peptideWithSetModifications.Length + 2)
            {
                return peptideWithSetModifications.OneBasedEndResidue;
            }
            return peptideWithSetModifications.OneBasedStartResidue + oneIsNterminus - 2;
        }

        #endregion

        protected static void WritePsmsToTsv(IEnumerable<SpectralMatch> psms, string filePath, IReadOnlyDictionary<string, int> modstoWritePruned, bool writePeptideLevelResults = false)
        {
            
            using (StreamWriter output = new StreamWriter(filePath))
            {
                bool includeOneOverK0Column = psms.Any(p => p.ScanOneOverK0.HasValue);
                output.WriteLine(SpectralMatch.GetTabSeparatedHeader(includeOneOverK0Column));
                foreach (var psm in psms)
                {
                    output.WriteLine(psm.ToString(modstoWritePruned, writePeptideLevelResults, includeOneOverK0Column));
                }
            }
        }

        protected static void WriteSpectrumLibrary(List<LibrarySpectrum> spectrumLibrary, string outputFolder)
        {
            var startTimeForAllFilenames = DateTime.Now.ToString("yyyy-MM-dd-HH-mm-ss", CultureInfo.InvariantCulture);
            string spectrumFilePath = outputFolder + "\\SpectralLibrary" + "_" + startTimeForAllFilenames + ".msp";
            using (StreamWriter output = new StreamWriter(spectrumFilePath))
            {
                foreach (var x in spectrumLibrary)
                {
                    output.WriteLine(x.ToString());
                }

            }
        }

        protected string UpdateSpectralLibrary(List<LibrarySpectrum> spectrumLibrary, string outputFolder)
        {
            var startTimeForAllFilenames = DateTime.Now.ToString("yyyy-MM-dd-HH-mm-ss", CultureInfo.InvariantCulture);
            string spectrumFilePath = outputFolder + "\\updateSpectralLibrary" +"_" + startTimeForAllFilenames + ".msp";
            using (StreamWriter output = new StreamWriter(spectrumFilePath))
            {

                foreach (var x in spectrumLibrary)
                {
                    output.WriteLine(x.ToString());
                }

            }
            return spectrumFilePath;
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

                // all directories in the same directory as the bioPolymer database
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



        /// <summary>
        /// Handle ambiguity when two theoretical bioPolymers in the
        /// search space have the same accession number.
        /// The accession must be unique for indexed searches
        ///
        /// RemoveAll is important because it references the location in memory, not the Equals
        /// </summary>
        protected static void SanitizeBioPolymerDatabase<TBioPolymer>(List<TBioPolymer> bioPolymers, TargetContaminantAmbiguity tcAmbiguity)
            where TBioPolymer : IBioPolymer
        {
            List<TBioPolymer> toRemove = new();
            foreach (var accessionGroup in bioPolymers.GroupBy(p => p.Accession)
                         .Where(group => group.Count() > 1) // only keep the ones with multiple entries sharing an accession
                         .Select(group => group.OrderBy(p => p.OneBasedPossibleLocalizedModifications.Count) // order by mods then truncation products (this is what was here before)
                             .ThenBy(p => p.TruncationProducts.Count)
                             .ToList()) // Individual ordered accession group to list from IEnumerable
                         .ToList()) // Collapse entire group and sort enumerable to a list so we can modify the bioPolymers collection. 
            {
                toRemove.Clear();
                string accession = accessionGroup.First().Accession;
                
                if (tcAmbiguity == TargetContaminantAmbiguity.RenameProtein)
                {
                    int bioPolymerNumber = 1;
                    Warn("The protein '" + accession + "' has multiple entries. Protein accessions must be unique. Protein " + accession + " was renamed.");
                    foreach (var originalBioPolymer in accessionGroup)
                    {
                        //accession is private and there's no clone method, so we need to make a whole new bioPolymer... TODO: put this in mzlib
                        //use PROTEIN_D1 instead of PROTEIN_1 so it doesn't look like an isoform (D for Duplicate)
                        IBioPolymer renamed;
                        if (originalBioPolymer is RNA r)
                        {
                            renamed = new RNA(originalBioPolymer.BaseSequence, originalBioPolymer.Accession + "_D" + bioPolymerNumber,
                                r.OneBasedPossibleLocalizedModifications, r.FivePrimeTerminus, r.ThreePrimeTerminus, r.Name, r.Organism,
                                r.DatabaseFilePath, r.IsContaminant, r.IsDecoy, r.GeneNames, r.AdditionalDatabaseFields, r.TruncationProducts,
                                r.SequenceVariations, r.AppliedSequenceVariations, r.SampleNameForVariants, r.FullName);
                        }
                        else
                        {
                            Protein p = originalBioPolymer as Protein ?? throw new ArgumentException($"Database sanitization assumed BioPolymer was a protein when it was {originalBioPolymer.GetType()}");
                            renamed = new Protein(originalBioPolymer.BaseSequence, originalBioPolymer.Accession + "_D" + bioPolymerNumber, originalBioPolymer.Organism,
                                originalBioPolymer.GeneNames, originalBioPolymer.OneBasedPossibleLocalizedModifications, p.TruncationProducts, originalBioPolymer.Name, originalBioPolymer.FullName,
                                originalBioPolymer.IsDecoy, originalBioPolymer.IsContaminant, p.DatabaseReferences, p.SequenceVariations, p.AppliedSequenceVariations,
                                p.SampleNameForVariants, p.DisulfideBonds, p.SpliceSites, originalBioPolymer.DatabaseFilePath);
                        }

                        bioPolymers.Add((TBioPolymer)renamed);
                        bioPolymers.RemoveAll(m => ReferenceEquals(m, originalBioPolymer));
                        bioPolymerNumber++;
                    }

                    continue;
                }

                // if we are not renaming, we need to remove the duplicates
                if (tcAmbiguity == TargetContaminantAmbiguity.RemoveContaminant)
                {
                    // remove contaminants
                    toRemove.AddRange(accessionGroup.Where(p => p.IsContaminant));
                    if (toRemove.Any())
                    {
                        Warn($"The {GlobalVariables.AnalyteType.GetBioPolymerLabel()}s '" + accession + $"' has multiple entries. {GlobalVariables.AnalyteType.GetBioPolymerLabel()}s accessions must be unique. Contaminant {GlobalVariables.AnalyteType.GetBioPolymerLabel()}s " + accession + " was ignored.");
                    }
                }
                else if (tcAmbiguity == TargetContaminantAmbiguity.RemoveTarget)
                {
                    // remove targets
                    toRemove.AddRange(accessionGroup.Where(p => !p.IsDecoy && !p.IsContaminant));
                    if (toRemove.Any())
                    {
                        Warn($"The {GlobalVariables.AnalyteType.GetBioPolymerLabel()}s '" + accession + $"' has multiple entries. {GlobalVariables.AnalyteType.GetBioPolymerLabel()}s accessions must be unique. Target {GlobalVariables.AnalyteType.GetBioPolymerLabel()}s " + accession + " was ignored.");
                    }
                }

                //remove the bioPolymers specified above
                foreach (var bioPolymer in toRemove.Where(_ => accessionGroup.Count > 1))
                {
                    bioPolymers.RemoveAll(p => ReferenceEquals(p, bioPolymer));
                    accessionGroup.RemoveAll(p => ReferenceEquals(p, bioPolymer));
                }

                // most ambiguity should be handled by now, but for edge cases and decoys:
                // remove bioPolymers so that only 1 bioPolymer with this accession remains
                for (int i = 0; i < accessionGroup.Count - 1; i++) //-1 to keep the last one (most mods)
                {
                    bioPolymers.RemoveAll(p => ReferenceEquals(p, accessionGroup[i]));
                }
            }
        }
    }
}
