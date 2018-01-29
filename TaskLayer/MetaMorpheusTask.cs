using Chemistry;
using EngineLayer;
using MassSpectrometry;
using MzLibUtil;
using Nett;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Diagnostics;
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
        XLSearch,
        Neo
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
                        .ConfigureType<Protease>(type => type
                            .WithConversionFor<TomlString>(convert => convert
                                .ToToml(custom => custom.ToString())
                                .FromToml(tmlString => GlobalVariables.ProteaseDictionary[tmlString.Value])))
                        .ConfigureType<ICommonParameters>(ct => ct
                            .CreateInstance(() => new CommonParameters()))
                        .ConfigureType<IDigestionParams>(ct => ct
                            .CreateInstance(() => new DigestionParams()))
                        .ConfigureType<List<(string, string)>>(type => type
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

        public static event EventHandler<StringEventArgs> LogHandler;

        public static event EventHandler<StringEventArgs> NewCollectionHandler;

        public static event EventHandler<ProgressEventArgs> OutProgressHandler;

        #endregion Public Events

        #region Public Properties

        public MyTask TaskType { get; set; }

        public ICommonParameters CommonParameters { get; set; }

        #endregion Public Properties

        #region Public Methods

        public static IEnumerable<Ms2ScanWithSpecificMass> GetMs2Scans(
         IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMSDataFile,
         string fullFilePath,
         bool doPrecursorDeconvolution,
         bool useProvidedPrecursorInfo,
         double deconvolutionIntensityRatio,
         int deconvolutionMaxAssumedChargeState,
         Tolerance deconvolutionMassTolerance)
        {
            foreach (var ms2scan in myMSDataFile.OfType<IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>>>())
            {
                List<(double, int)> isolatedStuff = new List<(double, int)>();
                if (ms2scan.OneBasedPrecursorScanNumber.HasValue)
                {
                    var precursorSpectrum = myMSDataFile.GetOneBasedScan(ms2scan.OneBasedPrecursorScanNumber.Value);
                    ms2scan.RefineSelectedMzAndIntensity(precursorSpectrum.MassSpectrum);
                    if (ms2scan.SelectedIonMonoisotopicGuessMz.HasValue)
                        ms2scan.ComputeMonoisotopicPeakIntensity(precursorSpectrum.MassSpectrum);
                    if (doPrecursorDeconvolution)
                        foreach (var envelope in ms2scan.GetIsolatedMassesAndCharges(precursorSpectrum.MassSpectrum, 1, deconvolutionMaxAssumedChargeState, deconvolutionMassTolerance.Value, deconvolutionIntensityRatio))
                        {
                            var monoPeakMz = envelope.monoisotopicMass.ToMz(envelope.charge);
                            isolatedStuff.Add((monoPeakMz, envelope.charge));
                        }
                }

                if (useProvidedPrecursorInfo && ms2scan.SelectedIonChargeStateGuess.HasValue)
                {
                    var precursorCharge = ms2scan.SelectedIonChargeStateGuess.Value;
                    if (ms2scan.SelectedIonMonoisotopicGuessMz.HasValue)
                    {
                        var precursorMZ = ms2scan.SelectedIonMonoisotopicGuessMz.Value;
                        if (!isolatedStuff.Any(b => deconvolutionMassTolerance.Within(precursorMZ.ToMass(precursorCharge), b.Item1.ToMass(b.Item2))))
                            isolatedStuff.Add((precursorMZ, precursorCharge));
                    }
                    else
                    {
                        var precursorMZ = ms2scan.SelectedIonMZ;
                        if (!isolatedStuff.Any(b => deconvolutionMassTolerance.Within(precursorMZ.ToMass(precursorCharge), b.Item1.ToMass(b.Item2))))
                            isolatedStuff.Add((precursorMZ, precursorCharge));
                    }
                }

                foreach (var heh in isolatedStuff)
                    yield return new Ms2ScanWithSpecificMass(ms2scan, heh.Item1, heh.Item2, fullFilePath);
            }
        }

        public static ICommonParameters SetAllFileSpecificCommonParams(ICommonParameters commonParams, FileSpecificSettings currentFileSpecificSettings)
        {
            if (currentFileSpecificSettings == null)
                return commonParams;

            IDigestionParams computedDigestionParams = new DigestionParams
            {
                InitiatorMethionineBehavior = currentFileSpecificSettings.InitiatorMethionineBehavior.Equals(InitiatorMethionineBehavior.Undefined) ? commonParams.DigestionParams.InitiatorMethionineBehavior : currentFileSpecificSettings.InitiatorMethionineBehavior,
                MaxMissedCleavages = currentFileSpecificSettings.MaxMissedCleavages ?? commonParams.DigestionParams.MaxMissedCleavages,
                MinPeptideLength = currentFileSpecificSettings.MinPeptideLength ?? commonParams.DigestionParams.MinPeptideLength,
                MaxPeptideLength = currentFileSpecificSettings.MaxPeptideLength ?? commonParams.DigestionParams.MaxPeptideLength,
                MaxModificationIsoforms = currentFileSpecificSettings.MaxModificationIsoforms ?? commonParams.DigestionParams.MaxModificationIsoforms,
                Protease = currentFileSpecificSettings.Protease ?? commonParams.DigestionParams.Protease,
                MaxModsForPeptide = currentFileSpecificSettings.Max_mods_for_peptide ?? commonParams.DigestionParams.MaxModsForPeptide,
                SemiProteaseDigestion = currentFileSpecificSettings.SemiProteaseDigestion ?? commonParams.DigestionParams.SemiProteaseDigestion,
                TerminusTypeSemiProtease = currentFileSpecificSettings.TerminusTypeSemiProtease ?? commonParams.DigestionParams.TerminusTypeSemiProtease,
            };

            CommonParameters returnParams = new CommonParameters
            {
                BIons = commonParams.BIons,
                YIons = commonParams.YIons,
                CIons = commonParams.CIons,
                ZdotIons = commonParams.ZdotIons,
                ReportAllAmbiguity = commonParams.ReportAllAmbiguity,
                TotalPartitions = commonParams.TotalPartitions,
                LocalizeAll = commonParams.LocalizeAll,
                ListOfModsVariable = commonParams.ListOfModsVariable,
                ListOfModsFixed = commonParams.ListOfModsFixed,
                ListOfModTypesLocalize = commonParams.ListOfModTypesLocalize,
                ConserveMemory = commonParams.ConserveMemory,
                MaxParallelFilesToAnalyze = commonParams.MaxParallelFilesToAnalyze,
                MaxThreadsToUsePerFile = commonParams.MaxThreadsToUsePerFile,
                ScoreCutoff = commonParams.ScoreCutoff,
                CalculateEValue = commonParams.CalculateEValue,

                DoPrecursorDeconvolution = currentFileSpecificSettings.DoPrecursorDeconvolution ?? commonParams.DoPrecursorDeconvolution,
                UseProvidedPrecursorInfo = currentFileSpecificSettings.UseProvidedPrecursorInfo ?? commonParams.UseProvidedPrecursorInfo,
                DeconvolutionIntensityRatio = currentFileSpecificSettings.DeconvolutionIntensityRatio ?? commonParams.DeconvolutionIntensityRatio,            
                DeconvolutionMaxAssumedChargeState = currentFileSpecificSettings.DeconvolutionMaxAssumedChargeState ?? commonParams.DeconvolutionMaxAssumedChargeState,
                DeconvolutionMassTolerance = currentFileSpecificSettings.DeconvolutionMassTolerance ?? commonParams.DeconvolutionMassTolerance,
                ProductMassTolerance = currentFileSpecificSettings.ProductMassTolerance ?? commonParams.ProductMassTolerance,
                PrecursorMassTolerance = currentFileSpecificSettings.PrecursorMassTolerance ?? commonParams.PrecursorMassTolerance,
                TopNpeaks = currentFileSpecificSettings.TopNpeaks ?? commonParams.TopNpeaks,
                MinRatio = currentFileSpecificSettings.MinRatio ?? commonParams.MinRatio,
                TrimMs1Peaks = currentFileSpecificSettings.TrimMs1Peaks ?? commonParams.TrimMs1Peaks,
                TrimMsMsPeaks = currentFileSpecificSettings.TrimMsMsPeaks ?? commonParams.TrimMsMsPeaks,
                DigestionParams = computedDigestionParams,
            };

            return returnParams;
        }

        public MyTaskResults RunTask(string output_folder, List<DbForTask> currentProteinDbFilenameList, List<string> currentRawDataFilepathList, string displayName)
        {
            StartingSingleTask(displayName);

            #region write TOML

            var tomlFileName = Path.Combine(output_folder, GetType().Name + "config.toml");
            Toml.WriteFile(this, tomlFileName, tomlConfig);
            SucessfullyFinishedWritingFile(tomlFileName, new List<string> { displayName });

            #endregion write TOML

            MetaMorpheusEngine.FinishedSingleEngineHandler += SingleEngineHandlerInTask;
            try
            {
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

                RunSpecific(output_folder, currentProteinDbFilenameList, currentRawDataFilepathList, displayName, fileSettingsList);
                stopWatch.Stop();
                myTaskResults.Time = stopWatch.Elapsed;
                var resultsFileName = Path.Combine(output_folder, "results.txt");
                using (StreamWriter file = new StreamWriter(resultsFileName))
                {
                    file.WriteLine("MetaMorpheus: version " + GlobalVariables.MetaMorpheusVersion);
                    file.Write(myTaskResults.ToString());
                }
                SucessfullyFinishedWritingFile(resultsFileName, new List<string> { displayName });
                FinishedSingleTask(displayName);
            }
            catch (Exception e)
            {
                MetaMorpheusEngine.FinishedSingleEngineHandler -= SingleEngineHandlerInTask;
                var resultsFileName = Path.Combine(output_folder, "results.txt");
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

            #region Write prose

            {
                var proseFilePath = Path.Combine(output_folder, "prose.txt");
                using (StreamWriter file = new StreamWriter(proseFilePath))
                {
                    file.Write("The data analysis was performed using MetaMorpheus Version: " + GlobalVariables.MetaMorpheusVersion + ", available at " + "https://github.com/smith-chem-wisc/MetaMorpheus." + " [INSERT CITATION] ");
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
                SucessfullyFinishedWritingFile(proseFilePath, new List<string> { displayName });
            }

            #endregion Write prose

            MetaMorpheusEngine.FinishedSingleEngineHandler -= SingleEngineHandlerInTask;
            return myTaskResults;
        }

        #endregion Public Methods

        #region Protected Methods

        protected static void WritePsmsToTsv(IEnumerable<Psm> items, string filePath, IReadOnlyDictionary<string, int> ModstoWritePruned)
        {
            using (StreamWriter output = new StreamWriter(filePath))
            {
                output.WriteLine(Psm.GetTabSeparatedHeader());
                foreach (var heh in items)
                    output.WriteLine(heh.ToString(ModstoWritePruned));
            }
        }

        protected static List<Protein> LoadProteinDb(string fileName, bool generateTargets, DecoyType decoyType, List<string> localizeableModificationTypes, bool isContaminant, out Dictionary<string, Modification> um)
        {
            if (Path.GetExtension(fileName).Equals(".fasta"))
            {
                um = null;
                return ProteinDbLoader.LoadProteinFasta(fileName, generateTargets, decoyType, isContaminant, ProteinDbLoader.uniprot_accession_expression, ProteinDbLoader.uniprot_fullName_expression, ProteinDbLoader.uniprot_fullName_expression, ProteinDbLoader.uniprot_gene_expression);
            }
            else
            {
                List<string> modTypesToExclude = GlobalVariables.AllModTypesKnown.Where(b => !localizeableModificationTypes.Contains(b)).ToList();
                return ProteinDbLoader.LoadProteinXML(fileName, generateTargets, decoyType, GlobalVariables.AllModsKnown, isContaminant, modTypesToExclude, out um);
            }
        }

        protected static HashSet<IDigestionParams> GetListOfDistinctDigestionParams(ICommonParameters commonParameters, IEnumerable<ICommonParameters> enumerable)
        {
            HashSet<IDigestionParams> okay = new HashSet<IDigestionParams>
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

        protected void Warn(string v)
        {
            WarnHandler?.Invoke(this, new StringEventArgs(v, null));
        }

        protected void Log(string v, List<string> nestedIds)
        {
            LogHandler?.Invoke(this, new StringEventArgs(v, nestedIds));
        }

        protected void NewCollection(string displayName, List<string> nestedIds)
        {
            NewCollectionHandler?.Invoke(this, new StringEventArgs(displayName, nestedIds));
        }

        #endregion Protected Methods

        #region Private Methods

        private static List<(string, string)> GetModsFromString(string value)
        {
            return value.Split(new string[] { "\t\t" }, StringSplitOptions.RemoveEmptyEntries).Select(b => (b.Split('\t').First(), b.Split('\t').Last())).ToList();
        }

        private void SingleEngineHandlerInTask(object sender, SingleEngineFinishedEventArgs e)
        {
            myTaskResults.AddResultText(e.ToString());
        }

        private void FinishedSingleTask(string displayName)
        {
            FinishedSingleTaskHandler?.Invoke(this, new SingleTaskEventArgs(displayName));
        }

        private void StartingSingleTask(string displayName)
        {
            StartingSingleTaskHander?.Invoke(this, new SingleTaskEventArgs(displayName));
        }

        #endregion Private Methods
    }
}