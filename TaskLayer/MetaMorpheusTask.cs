using Chemistry;
using EngineLayer;
using MassSpectrometry;
using MzLibUtil;
using Nett;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using UsefulProteomicsDatabases;

namespace TaskLayer
{
    public enum MyTask
    {
        Search,
        Gptmd,
        Calibrate,
        XLSearch
    }

    public abstract class MetaMorpheusTask
    {
        #region Public Fields

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
                        .ConfigureType<MassDiffAcceptor>(type => type
                            .WithConversionFor<TomlString>(convert => convert
                                .FromToml(tmlString => ParseSearchMode(tmlString.Value))))
                        .ConfigureType<DotMassDiffAcceptor>(type => type
                            .WithConversionFor<TomlString>(convert => convert
                                .ToToml(custom => custom.ToString())))
                        .ConfigureType<IntervalMassDiffAcceptor>(type => type
                            .WithConversionFor<TomlString>(convert => convert
                                .ToToml(custom => custom.ToString())))
                        .ConfigureType<OpenSearchMode>(type => type
                            .WithConversionFor<TomlString>(convert => convert
                                .ToToml(custom => custom.ToString())))
                        .ConfigureType<SingleAbsoluteAroundZeroSearchMode>(type => type
                            .WithConversionFor<TomlString>(convert => convert
                                .ToToml(custom => custom.ToString())))
                        .ConfigureType<SinglePpmAroundZeroSearchMode>(type => type
                            .WithConversionFor<TomlString>(convert => convert
                                .ToToml(custom => custom.ToString())))
                        .ConfigureType<Protease>(type => type
                            .WithConversionFor<TomlString>(convert => convert
                                .ToToml(custom => custom.ToString())
                                .FromToml(tmlString => GlobalEngineLevelSettings.ProteaseDictionary[tmlString.Value])))
                        .ConfigureType<List<Tuple<string, string>>>(type => type
                             .WithConversionFor<TomlString>(convert => convert
                                 .ToToml(custom => string.Join("\t\t", custom.Select(b => b.Item1 + "\t" + b.Item2)))
                                 .FromToml(tmlString => GetModsFromString(tmlString.Value)))));

        #endregion Public Fields

        #region Protected Fields

        protected readonly StringBuilder proseCreatedWhileRunning = new StringBuilder();
        protected MyTaskResults myTaskResults;

        #endregion Protected Fields

        #region Protected Constructors

        protected MetaMorpheusTask(MyTask taskType)
        {
            this.TaskType = taskType;
        }

        #endregion Protected Constructors

        #region Public Events

        public static event EventHandler<SingleTaskEventArgs> FinishedSingleTaskHandler;

        public static event EventHandler<SingleFileEventArgs> FinishedWritingFileHandler;

        public static event EventHandler<SingleTaskEventArgs> StartingSingleTaskHander;

        public static event EventHandler<StringEventArgs> StartingDataFileHandler;

        public static event EventHandler<StringEventArgs> FinishedDataFileHandler;

        public static event EventHandler<StringEventArgs> OutLabelStatusHandler;

        public static event EventHandler<StringEventArgs> WarnHandler;

        public static event EventHandler<StringEventArgs> NewCollectionHandler;

        public static event EventHandler<ProgressEventArgs> OutProgressHandler;

        #endregion Public Events

        #region Public Properties

        public MyTask TaskType { get; set; }

        public CommonParameters CommonParameters { get; set; }

        #endregion Public Properties

        #region Public Methods

        public static IEnumerable<Ms2ScanWithSpecificMass> GetMs2Scans(
         IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMSDataFile, string fullFilePath,
         bool doPrecursorDeconvolution,
         bool useProvidedPrecursorInfo,
         double deconvolutionIntensityRatio,
         int deconvolutionMaxAssumedChargeState,
         Tolerance deconvolutionMassTolerance)
        {
            foreach (var ms2scan in myMSDataFile.OfType<IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>>>())
            {
                List<Tuple<List<IMzPeak>, int>> isolatedStuff = new List<Tuple<List<IMzPeak>, int>>();
                if (ms2scan.OneBasedPrecursorScanNumber.HasValue)
                {
                    var precursorSpectrum = myMSDataFile.GetOneBasedScan(ms2scan.OneBasedPrecursorScanNumber.Value);
                    ms2scan.RefineSelectedMzAndIntensity(precursorSpectrum.MassSpectrum);
                    ms2scan.ComputeMonoisotopicPeakIntensity(precursorSpectrum.MassSpectrum);
                    if (doPrecursorDeconvolution)
                        isolatedStuff.AddRange(ms2scan.GetIsolatedMassesAndCharges(precursorSpectrum.MassSpectrum, deconvolutionMaxAssumedChargeState, deconvolutionMassTolerance, deconvolutionIntensityRatio));
                }

                if (useProvidedPrecursorInfo && ms2scan.SelectedIonChargeStateGuess.HasValue)
                {
                    var precursorCharge = ms2scan.SelectedIonChargeStateGuess.Value;
                    if (ms2scan.SelectedIonMonoisotopicGuessMz.HasValue)
                    {
                        var precursorMZ = ms2scan.SelectedIonMonoisotopicGuessMz.Value;
                        if (!isolatedStuff.Any(b => deconvolutionMassTolerance.Within(precursorMZ.ToMass(precursorCharge), b.Item1.First().Mz.ToMass(b.Item2))))
                            isolatedStuff.Add(new Tuple<List<IMzPeak>, int>(new List<IMzPeak> { new MzPeak(precursorMZ, ms2scan.SelectedIonMonoisotopicGuessIntensity.Value) }, precursorCharge));
                    }
                    else
                    {
                        var PrecursorMZ = ms2scan.SelectedIonMZ;
                        if (!isolatedStuff.Any(b => deconvolutionMassTolerance.Within(PrecursorMZ.ToMass(precursorCharge), b.Item1.First().Mz.ToMass(b.Item2))))
                            isolatedStuff.Add(new Tuple<List<IMzPeak>, int>(new List<IMzPeak> { new MzPeak(PrecursorMZ, ms2scan.SelectedIonIntensity ?? double.NaN) }, precursorCharge));
                    }
                }

                foreach (var heh in isolatedStuff)
                    yield return new Ms2ScanWithSpecificMass(ms2scan, heh.Item1.First(), heh.Item2, fullFilePath);
            }
        }

        public static MassDiffAcceptor ParseSearchMode(string text)
        {
            MassDiffAcceptor ye = null;

            var split = text.Split(' ');

            switch (split[1])
            {
                case "dot":

                    var massShifts = Array.ConvertAll(split[4].Split(','), Double.Parse);
                    var newString = split[2].Replace("±", "");
                    var toleranceValue = double.Parse(newString, CultureInfo.InvariantCulture);
                    if (split[3].ToUpperInvariant().Equals("PPM"))
                        ye = new DotMassDiffAcceptor(split[0], massShifts, new PpmTolerance(toleranceValue));
                    else if (split[3].ToUpperInvariant().Equals("DA"))
                        ye = new DotMassDiffAcceptor(split[0], massShifts, new AbsoluteTolerance(toleranceValue));
                    break;

                case "interval":
                    IEnumerable<DoubleRange> doubleRanges = Array.ConvertAll(split[2].Split(','), b => new DoubleRange(double.Parse(b.Trim(new char[] { '[', ']' }).Split(';')[0], CultureInfo.InvariantCulture), double.Parse(b.Trim(new char[] { '[', ']' }).Split(';')[1], CultureInfo.InvariantCulture)));
                    ye = new IntervalMassDiffAcceptor(split[0], doubleRanges);
                    break;

                case "OpenSearch":
                    ye = new OpenSearchMode();
                    break;

                case "daltonsAroundZero":
                    ye = new SingleAbsoluteAroundZeroSearchMode(double.Parse(split[2], CultureInfo.InvariantCulture));
                    break;

                case "ppmAroundZero":
                    ye = new SinglePpmAroundZeroSearchMode(double.Parse(split[2], CultureInfo.InvariantCulture));
                    break;

                default:
                    throw new MetaMorpheusException("Could not parse search mode string");
            }
            return ye;
        }

        public static CommonParameters SetAllFileSpecificCommonParams(CommonParameters commonParams, FileSpecificSettings currentFileSpecificSettings)
        {
            if (currentFileSpecificSettings == null)
                return commonParams;

            CommonParameters returnParams = new CommonParameters
            {
                DoPrecursorDeconvolution = currentFileSpecificSettings.DoPrecursorDeconvolution ?? commonParams.DoPrecursorDeconvolution,

                UseProvidedPrecursorInfo = currentFileSpecificSettings.UseProvidedPrecursorInfo ?? commonParams.UseProvidedPrecursorInfo,

                DeconvolutionIntensityRatio = currentFileSpecificSettings.DeconvolutionIntensityRatio ?? commonParams.DeconvolutionIntensityRatio,

                DeconvolutionMaxAssumedChargeState = currentFileSpecificSettings.DeconvolutionMaxAssumedChargeState ?? commonParams.DeconvolutionMaxAssumedChargeState,

                DeconvolutionMassTolerance = currentFileSpecificSettings.DeconvolutionMassTolerance ?? commonParams.DeconvolutionMassTolerance,

                TotalPartitions = currentFileSpecificSettings.TotalPartitions ?? commonParams.TotalPartitions,

                ProductMassTolerance = currentFileSpecificSettings.ProductMassTolerance ?? commonParams.ProductMassTolerance,

                ConserveMemory = currentFileSpecificSettings.ConserveMemory ?? commonParams.ConserveMemory,

                ScoreCutoff = currentFileSpecificSettings.ScoreCutoff ?? commonParams.ScoreCutoff,

                TopNpeaks = currentFileSpecificSettings.TopNpeaks ?? commonParams.TopNpeaks,
                MinRatio = currentFileSpecificSettings.MinRatio ?? commonParams.MinRatio,
                TrimMs1Peaks = currentFileSpecificSettings.TrimMs1Peaks ?? commonParams.TrimMs1Peaks,
                TrimMsMsPeaks = currentFileSpecificSettings.TrimMsMsPeaks ?? commonParams.TrimMsMsPeaks,
            };

            returnParams.DigestionParams.InitiatorMethionineBehavior = currentFileSpecificSettings.InitiatorMethionineBehavior.Equals(InitiatorMethionineBehavior.Undefined) ? commonParams.DigestionParams.InitiatorMethionineBehavior : currentFileSpecificSettings.InitiatorMethionineBehavior;
            returnParams.DigestionParams.MaxMissedCleavages = currentFileSpecificSettings.MaxMissedCleavages ?? commonParams.DigestionParams.MaxMissedCleavages;
            returnParams.DigestionParams.MinPeptideLength = currentFileSpecificSettings.MinPeptideLength ?? commonParams.DigestionParams.MinPeptideLength;
            returnParams.DigestionParams.MaxPeptideLength = currentFileSpecificSettings.MaxPeptideLength ?? commonParams.DigestionParams.MaxPeptideLength;
            returnParams.DigestionParams.MaxModificationIsoforms = currentFileSpecificSettings.MaxModificationIsoforms ?? commonParams.DigestionParams.MaxModificationIsoforms;
            returnParams.DigestionParams.Protease = currentFileSpecificSettings.Protease ?? commonParams.DigestionParams.Protease;
            returnParams.DigestionParams.MaxModsForPeptide = currentFileSpecificSettings.Max_mods_for_peptide ?? commonParams.DigestionParams.MaxModsForPeptide;

            return returnParams;
        }

        public MyTaskResults RunTask(string output_folder, List<DbForTask> currentProteinDbFilenameList, List<string> currentRawDataFilepathList, string taskId)
        {
            StartingSingleTask(taskId);

            #region write TOML

            {
                var tomlFileName = Path.Combine(output_folder, GetType().Name + "config.toml");
                Toml.WriteFile(this, tomlFileName, tomlConfig);
                SucessfullyFinishedWritingFile(tomlFileName, new List<string> { taskId });
            }

            #endregion write TOML

            MetaMorpheusEngine.FinishedSingleEngineHandler += SingleEngineHandlerInTask;
#if !DEBUG
            try
            {
#endif
            var stopWatch = new Stopwatch();
            stopWatch.Start();

            FileSpecificSettings[] fileSettingsList = new FileSpecificSettings[currentRawDataFilepathList.Count];
            for (int i = 0; i < currentRawDataFilepathList.Count; i++)
            {
                string rawFilePath = currentRawDataFilepathList[i];
                var fileSpecificToml = Directory.GetFiles(Directory.GetParent(rawFilePath).ToString(), Path.GetFileNameWithoutExtension(rawFilePath) + ".toml");
                //Will only enter if Toml file exists with same name
                if (fileSpecificToml.Length == 1)
                {
                    TomlTable fileSpecificSettings = Toml.ReadFile(fileSpecificToml[0], tomlConfig);
                    var tomlSettingsList = fileSpecificSettings.ToDictionary(p => p.Key);
                    fileSettingsList[i] = new FileSpecificSettings(tomlSettingsList);
                }
            }

            RunSpecific(output_folder, currentProteinDbFilenameList, currentRawDataFilepathList, taskId, fileSettingsList);
            stopWatch.Stop();
            myTaskResults.Time = stopWatch.Elapsed;
            var resultsFileName = Path.Combine(output_folder, "results.txt");
            using (StreamWriter file = new StreamWriter(resultsFileName))
            {
                file.WriteLine("MetaMorpheus: version " + GlobalEngineLevelSettings.MetaMorpheusVersion);
                file.Write(myTaskResults.ToString());
            }
            SucessfullyFinishedWritingFile(resultsFileName, new List<string> { taskId });
            FinishedSingleTask(taskId);
#if !DEBUG
            }
            catch (Exception e)
            {
                MetaMorpheusEngine.FinishedSingleEngineHandler -= SingleEngineHandlerInTask;
                var resultsFileName = Path.Combine(output_folder, "results.txt");
                using (StreamWriter file = new StreamWriter(resultsFileName))
                {
                    file.WriteLine(GlobalEngineLevelSettings.MetaMorpheusVersion.Equals("1.0.0.0") ? "MetaMorpheus: Not a release version" : "MetaMorpheus: version " + GlobalEngineLevelSettings.MetaMorpheusVersion);
                    file.WriteLine(MzLibUtil.SystemInfo.CompleteSystemInfo()); //OS, OS Version, .Net Version, RAM, processor count, MSFileReader .dll versions X3
                    file.Write("e: " + e);
                    file.Write("e.Message: " + e.Message);
                    file.Write("e.InnerException: " + e.InnerException);
                    file.Write("e.Source: " + e.Source);
                    file.Write("e.StackTrace: " + e.StackTrace);
                    file.Write("e.TargetSite: " + e.TargetSite);
                }
                throw;
            }
#endif

            #region Write prose

            {
                var proseFilePath = Path.Combine(output_folder, "prose.txt");
                using (StreamWriter file = new StreamWriter(proseFilePath))
                {
                    file.Write("The data analysis was performed using MetaMorpheus Version: " + GlobalEngineLevelSettings.MetaMorpheusVersion + ", available at " + "https://github.com/smith-chem-wisc/MetaMorpheus." + " [INSERT CITATION] ");
                    file.Write(proseCreatedWhileRunning.ToString());
                    file.Write(SystemInfo.SystemProse().Replace(Environment.NewLine, "") + " ");
                    file.WriteLine("The total time to perform " + this.TaskType + " task on " + currentRawDataFilepathList.Count + " spectra file(s) was " + String.Format("{0:0.00}", myTaskResults.Time.TotalMinutes) + " minutes.");
                    file.WriteLine();
                    file.WriteLine("Published works using MetaMorpheus software are encouraged to cite: STEFAN'S VERY IMPORTANT PAPER");

                    file.WriteLine();
                    file.WriteLine("Spectra files: ");
                    file.WriteLine(string.Join(Environment.NewLine, currentRawDataFilepathList.Select(b => '\t' + b)));
                    file.WriteLine("Databases:");
                    file.Write(string.Join(Environment.NewLine, currentProteinDbFilenameList.Select(b => '\t' + (b.IsContaminant ? "Contaminant " : "") + b.FilePath)));
                }
                SucessfullyFinishedWritingFile(proseFilePath, new List<string> { taskId });
            }

            #endregion Write prose

            MetaMorpheusEngine.FinishedSingleEngineHandler -= SingleEngineHandlerInTask;
            return myTaskResults;
        }

        #endregion Public Methods

        #region Protected Methods

        protected static void WritePsmsToTsv(IEnumerable<Psm> items, string filePath)
        {
            using (StreamWriter output = new StreamWriter(filePath))
            {
                output.WriteLine(Psm.GetTabSeparatedHeader());
                foreach (var heh in items)
                    output.WriteLine(heh);
            }
        }

        protected static List<Protein> LoadProteinDb(string fileName, bool generateDecoys, List<ModificationWithMass> localizeableModifications, bool isContaminant, out Dictionary<string, Modification> um)
        {
            if (Path.GetExtension(fileName).Equals(".fasta"))
            {
                um = null;
                return ProteinDbLoader.LoadProteinFasta(fileName, true, generateDecoys, isContaminant, ProteinDbLoader.uniprot_accession_expression, ProteinDbLoader.uniprot_fullName_expression, ProteinDbLoader.uniprot_fullName_expression, ProteinDbLoader.uniprot_gene_expression);
            }
            else
                return ProteinDbLoader.LoadProteinXML(fileName, true, generateDecoys, localizeableModifications, isContaminant, new List<string>(), out um);
        }

        protected static HashSet<DigestionParams> GetListOfDistinctDigestionParams(CommonParameters commonParameters, IEnumerable<CommonParameters> enumerable)
        {
            HashSet<DigestionParams> okay = new HashSet<DigestionParams>
            {
                commonParameters.DigestionParams
            };

            foreach (var hah in enumerable)
                okay.Add(hah.DigestionParams);

            return okay;
        }

        protected void ReportProgress(ProgressEventArgs v)
        {
            OutProgressHandler?.Invoke(this, v);
        }

        protected abstract MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId, FileSpecificSettings[] fileSettingsList);

        protected void SucessfullyFinishedWritingFile(string path, List<string> nestedIDs)
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

        protected void Warn(string v, List<string> nestedIds)
        {
            WarnHandler?.Invoke(this, new StringEventArgs(v, nestedIds));
        }

        protected void NewCollection(string v, List<string> nestedIds)
        {
            NewCollectionHandler?.Invoke(this, new StringEventArgs(v, nestedIds));
        }

        #endregion Protected Methods

        #region Private Methods

        private static List<Tuple<string, string>> GetModsFromString(string value)
        {
            return value.Split(new string[] { "\t\t" }, StringSplitOptions.RemoveEmptyEntries).Select(b => new Tuple<string, string>(b.Split('\t').First(), b.Split('\t').Last())).ToList();
        }

        private void SingleEngineHandlerInTask(object sender, SingleEngineFinishedEventArgs e)
        {
            myTaskResults.AddResultText(e.ToString());
        }

        private void FinishedSingleTask(string taskId)
        {
            FinishedSingleTaskHandler?.Invoke(this, new SingleTaskEventArgs(taskId));
        }

        private void StartingSingleTask(string taskId)
        {
            StartingSingleTaskHander?.Invoke(this, new SingleTaskEventArgs(taskId));
        }

        #endregion Private Methods
    }
}