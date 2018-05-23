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
                        .ConfigureType<CommonParameters> (type => type
                            .WithConversionFor<TomlString>(convert => convert
                                .ToToml(custom => custom.ToString())
                                .FromToml(tmlString => GetCommonParametersFromString(tmlString.Value))))
                        .ConfigureType<DigestionParams>(type => type
                            .WithConversionFor<TomlString>(convert => convert
                                .ToToml(custom => custom.ToString())
                                .FromToml(tmlString => GetDigestionParamsFromString(tmlString.Value))))
                        .ConfigureType<List<string>>(type => type
                             .WithConversionFor<TomlString>(convert => convert
                                 .ToToml(custom => string.Join("\t", custom))
                                 .FromToml(tmlString => GetModsTypesFromString(tmlString.Value))))
                        .ConfigureType<List<(string, string)>>(type => type
                             .WithConversionFor<TomlString>(convert => convert
                                 .ToToml(custom => string.Join("\t\t", custom.Select(b => b.Item1 + "\t" + b.Item2)))
                                 .FromToml(tmlString => GetModsFromString(tmlString.Value)))));

        #endregion Public Fields

        #region Protected Fields

        protected readonly StringBuilder ProseCreatedWhileRunning = new StringBuilder();

        protected MyTaskResults MyTaskResults;

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

        public CommonParameters CommonParams { get; set; }

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

        public static CommonParameters SetAllFileSpecificCommonParams(CommonParameters commonParams, FileSpecificParameters fileSpecificParams)
        {
            if (fileSpecificParams == null)
                return commonParams;

            // set file-specific digestion parameters
            else
            {
                Protease protease = fileSpecificParams.Protease ?? commonParams.DigestionParams.Protease;
                int MinPeptideLength = fileSpecificParams.MinPeptideLength ?? commonParams.DigestionParams.MinPeptideLength;
                int MaxPeptideLength = fileSpecificParams.MaxPeptideLength ?? commonParams.DigestionParams.MaxPeptideLength;
                int MaxMissedCleavages = fileSpecificParams.MaxMissedCleavages ?? commonParams.DigestionParams.MaxMissedCleavages;
                int MaxModsForPeptide = fileSpecificParams.MaxModsForPeptide ?? commonParams.DigestionParams.MaxModsForPeptide;
                DigestionParams fileSpecificDigestionParams = new DigestionParams(protease: protease.Name, MaxMissedCleavages: MaxMissedCleavages, MinPeptideLength: MinPeptideLength, MaxPeptideLength: MaxPeptideLength, MaxModsForPeptides: MaxModsForPeptide);


                // set the rest of the file-specific parameters
                Tolerance PrecursorMassTolerance = fileSpecificParams.PrecursorMassTolerance ?? commonParams.PrecursorMassTolerance;
                double preMassTol = (double)PrecursorMassTolerance.Value;
                Tolerance ProductMassTolerance = fileSpecificParams.ProductMassTolerance ?? commonParams.ProductMassTolerance;
                double prodMassTol = (double)ProductMassTolerance.Value;
                bool BIons = fileSpecificParams.BIons ?? commonParams.BIons;
                bool YIons = fileSpecificParams.BIons ?? commonParams.YIons;
                bool CIons = fileSpecificParams.CIons ?? commonParams.CIons;
                bool ZdotIons = fileSpecificParams.ZdotIons ?? commonParams.ZdotIons;
                CommonParameters returnParams = new CommonParameters(DigestionParams: fileSpecificDigestionParams, preMassTol: preMassTol, prodMassTol: prodMassTol, BIons: BIons, YIons: YIons, CIons: CIons, ZdotIons: ZdotIons);
                return returnParams;
            }
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

                FileSpecificParameters[] fileSettingsList = new FileSpecificParameters[currentRawDataFilepathList.Count];
                for (int i = 0; i < currentRawDataFilepathList.Count; i++)
                {
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
                            // file-specific toml has already been validated in the GUI when the spectra files were added, so...
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
                SucessfullyFinishedWritingFile(resultsFileName, new List<string> { displayName });
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

            #region Write prose

            {
                var proseFilePath = Path.Combine(output_folder, "prose.txt");
                using (StreamWriter file = new StreamWriter(proseFilePath))
                {
                    file.Write("The data analysis was performed using MetaMorpheus version " + GlobalVariables.MetaMorpheusVersion + ", available at " + "https://github.com/smith-chem-wisc/MetaMorpheus." + " [INSERT CITATION] ");
                    file.Write(ProseCreatedWhileRunning.ToString());
                    file.Write(SystemInfo.SystemProse().Replace(Environment.NewLine, "") + " ");
                    file.WriteLine("The total time to perform the " + this.TaskType.ToString().ToLowerInvariant() + " task on " + currentRawDataFilepathList.Count + " spectra file(s) was " + String.Format("{0:0.00}", MyTaskResults.Time.TotalMinutes) + " minutes.");
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
            return MyTaskResults;
        }

        #endregion Public Methods

        #region Protected Methods

        protected List<Protein> LoadProteins(string taskId, List<DbForTask> dbFilenameList, bool searchTarget, DecoyType decoyType, List<string> localizeableModificationTypes)
        {
            Status("Loading proteins...", new List<string> { taskId });
            int emptyProteinEntries = 0;
            List<Protein> proteinList = new List<Protein>();
            foreach (var db in dbFilenameList)
            {
                int emptyProteinEntriesForThisDb = 0;
                var dbProteinList = LoadProteinDb(db.FilePath, searchTarget, decoyType, localizeableModificationTypes, db.IsContaminant, out Dictionary<string, Modification> unknownModifications, out emptyProteinEntriesForThisDb);
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

        protected static List<Protein> LoadProteinDb(string fileName, bool generateTargets, DecoyType decoyType, List<string> localizeableModificationTypes, bool isContaminant, out Dictionary<string, Modification> um, out int emptyEntriesCount)
        {
            List<string> dbErrors = new List<string>();
            List<Protein> proteinList = new List<Protein>();

            string theExtension = Path.GetExtension(fileName).ToLowerInvariant();
            bool compressed = theExtension.EndsWith("gz"); // allows for .bgz and .tgz, too which are used on occasion
            theExtension = compressed ? Path.GetExtension(Path.GetFileNameWithoutExtension(fileName)).ToLowerInvariant() : theExtension;
            
            if (theExtension.Equals(".fasta") || theExtension.Equals(".fa"))
            {
                um = null;
                proteinList = ProteinDbLoader.LoadProteinFasta(fileName, generateTargets, decoyType, isContaminant, ProteinDbLoader.UniprotAccessionRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotGeneNameRegex, ProteinDbLoader.UniprotOrganismRegex, out dbErrors);
            }
            else
            {
                List<string> modTypesToExclude = GlobalVariables.AllModTypesKnown.Where(b => !localizeableModificationTypes.Contains(b)).ToList();
                proteinList = ProteinDbLoader.LoadProteinXML(fileName, generateTargets, decoyType, GlobalVariables.AllModsKnown, isContaminant, modTypesToExclude, out um);
            }

            emptyEntriesCount = proteinList.Count(p => p.BaseSequence.Length == 0);
            return proteinList.Where(p => p.BaseSequence.Length > 0).ToList();
        }
        protected static void WritePsmsToTsv(IEnumerable<PeptideSpectralMatch> items, string filePath, IReadOnlyDictionary<string, int> ModstoWritePruned)
        {
            using (StreamWriter output = new StreamWriter(filePath))
            {
                output.WriteLine(PeptideSpectralMatch.GetTabSeparatedHeader());
                foreach (var heh in items)
                {
                    output.WriteLine(heh.ToString(ModstoWritePruned));
                }
           }
        }
        
        protected void ReportProgress(ProgressEventArgs v)
        {
            OutProgressHandler?.Invoke(this, v);
        }

        protected abstract MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId, FileSpecificParameters[] fileSettingsList);

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

        private static CommonParameters GetCommonParametersFromString(string value)
        {
            CommonParameters temp = new CommonParameters();

            Dictionary<string, string> parameters = new Dictionary<string, string>();
            string[] input = value.Split(';');
            
            for (int i = 0; i < input.Length-1; i++)
            {
                
                string a = input[i];
                string[] temps = a.Split('=');
                string key = temps[0].Trim();
                string val = temps[1];
                parameters.Add(key, val);
               
            }


            //bool BIons = Convert.ToBoolean(output[0, 1]);
            bool BIons = Convert.ToBoolean(parameters[nameof(temp.BIons)]);
            bool YIons = Convert.ToBoolean(parameters[nameof(temp.YIons)]);
            bool ZdotIons = Convert.ToBoolean(parameters[nameof(temp.ZdotIons)]);
            bool CIons = Convert.ToBoolean(parameters[nameof(temp.CIons)]);
            bool DoPrecursorDeconvolution = Convert.ToBoolean(parameters[nameof(temp.DoPrecursorDeconvolution)]);
            bool UseProvidedPrecursorInfo = Convert.ToBoolean(parameters[nameof(temp.UseProvidedPrecursorInfo)]);
            double DeconvolutionIntensityRatio = Convert.ToDouble(parameters[nameof(temp.DeconvolutionIntensityRatio)]);
            int DeconvolutionMaxAssumedChargeState = Convert.ToInt32(parameters[nameof(temp.DeconvolutionMaxAssumedChargeState)]);
            bool ReportAllAmbiguity = Convert.ToBoolean(parameters[nameof(temp.ReportAllAmbiguity)]);
            bool CompIons = Convert.ToBoolean(parameters[nameof(temp.CompIons)]);
            int TotalPartitions = Convert.ToInt32(parameters[nameof(temp.TotalPartitions)]);
            double ScoreCutoff = Convert.ToDouble(parameters[nameof(temp.ScoreCutoff)]);
            int TopNpeaks = Convert.ToInt32(parameters[nameof(temp.TopNpeaks)]);
            double MinRatio = Convert.ToDouble(parameters[nameof(temp.MinRatio)]);
            bool TrimMs1Peaks = Convert.ToBoolean(parameters[nameof(temp.TrimMs1Peaks)]);
            bool TrimMsMsPeaks = Convert.ToBoolean(parameters[nameof(temp.TrimMsMsPeaks)]);
            bool UseDeltaScore = Convert.ToBoolean(parameters[nameof(temp.UseDeltaScore)]);
            bool CalculateEValue = Convert.ToBoolean(parameters[nameof(temp.CalculateEValue)]);
            double prodMassTol = Convert.ToDouble(parameters[nameof(temp.ProductMassTolerance)]);
            double preMassTol = Convert.ToDouble(parameters[nameof(temp.PrecursorMassTolerance)]);
            double deconMassTol = Convert.ToDouble(parameters[nameof(temp.DeconvolutionMassTolerance)]);
            int MaxThreadsToUsePerFile = Convert.ToInt32(parameters[nameof(temp.MaxThreadsToUsePerFile)]);

            string protease = parameters[nameof(temp.DigestionParams.Protease)];
            int MaxMissedCleavages = Convert.ToInt32(parameters[nameof(temp.DigestionParams.MaxMissedCleavages)]);
            int MinPeptideLength = Convert.ToInt32(parameters[nameof(temp.DigestionParams.MinPeptideLength)]);
            int MaxPeptideLength = Convert.ToInt32(parameters[nameof(temp.DigestionParams.MaxPeptideLength)]);
            int MaxModificationIsoforms = Convert.ToInt32(parameters[nameof(temp.DigestionParams.MaxModificationIsoforms)]);
            InitiatorMethionineBehavior InitiatorMethionineBehavior;
            if (parameters[nameof(temp.DigestionParams.InitiatorMethionineBehavior)] == "InitiatorMethionineBehavior.Cleave")
            {
                InitiatorMethionineBehavior = InitiatorMethionineBehavior.Cleave;
            }
            else if (parameters[nameof(temp.DigestionParams.InitiatorMethionineBehavior)] == "InitiatorMethionineBehavior.Retain")
            {
                InitiatorMethionineBehavior = InitiatorMethionineBehavior.Retain;
            }
            else if (parameters[nameof(temp.DigestionParams.InitiatorMethionineBehavior)] == "InitiatorMethionineBehavior.Variable")
            {
                InitiatorMethionineBehavior = InitiatorMethionineBehavior.Variable;
            }
            else
            {
                InitiatorMethionineBehavior = InitiatorMethionineBehavior.Undefined;
            }
            int MaxModsForPeptides = Convert.ToInt32(parameters[nameof(temp.DigestionParams.MaxModsForPeptide)]);
            bool SemiProteaseDigestion = Convert.ToBoolean(parameters[nameof(temp.DigestionParams.SemiProteaseDigestion)]);
            //string term = ;
            TerminusType TerminusTypeSemiProtease;
            if (parameters[nameof(temp.DigestionParams.TerminusTypeSemiProtease)] == "TerminusType.N")
            {
                TerminusTypeSemiProtease = TerminusType.N;
            }
            else if (parameters[nameof(temp.DigestionParams.TerminusTypeSemiProtease)] == "TerminusType.C")
            {
                TerminusTypeSemiProtease = TerminusType.C;
            }
            else
            {
                TerminusTypeSemiProtease = TerminusType.None;
            }
            
            DigestionParams DigestionParams = new DigestionParams(protease: protease, MaxMissedCleavages: MaxMissedCleavages, MinPeptideLength: MinPeptideLength, MaxPeptideLength: MaxPeptideLength, MaxModificationIsoforms: MaxModificationIsoforms, InitiatorMethionineBehavior: InitiatorMethionineBehavior, MaxModsForPeptides: MaxModsForPeptides, SemiProteaseDigestion: SemiProteaseDigestion, TerminusTypeSemiProtease:TerminusTypeSemiProtease);
            string listofModsVariable = parameters[nameof(temp.ListOfModsVariable)];
            string[] listContentVariable = listofModsVariable.Split(',');
            string[,] listContentStringsVariable = new string[listContentVariable.Length, 2];
            List<(string, string)> ListOfModsVariable = new List<(string, string)>();
            for (int i = 0; i < listContentVariable.Length-1; i++)
            {

                string a = listContentVariable[i];
                string[] temps = a.Split('/');
                string key = temps[0];
                string val = temps[1];

                listContentStringsVariable[i, 0] = key;
                listContentStringsVariable[i, 1] = val;
            }

            for (int i = 0; i < listContentVariable.Length-1; i++)
            {
                string string1 = listContentStringsVariable[i, 0];
                string string2 = listContentStringsVariable[i, 1];
                (string, string) item = (string1, string2);
                ListOfModsVariable.Add(item);

            }


            string listofModsFixed = parameters[nameof(temp.ListOfModsFixed)];
            string[] listContentFixed = listofModsFixed.Split(',');
            string[,] listContentStringsFixed = new string[listContentFixed.Length, 2];
            List<(string, string)> ListOfModsFixed = new List<(string, string)>();

            for (int i = 0; i < listContentFixed.Length-1; i++)
            {

                string a = listContentFixed[i];
                string[] temps = a.Split('/');
                string key = temps[0];
                string val = temps[1];

                listContentStringsFixed[i, 0] = key;
                listContentStringsFixed[i, 1] = val;
            }

            
            for (int i = 0; i < listContentFixed.Length-1; i++)
            {
                string string1 = listContentStringsFixed[i, 0];
                string string2 = listContentStringsFixed[i, 1];
                (string, string) item = (string1, string2);
                ListOfModsFixed.Add(item);

            }

            return new CommonParameters(BIons: BIons , YIons: YIons, ZdotIons: ZdotIons, CIons: CIons, DoPrecursorDeconvolution: DoPrecursorDeconvolution, UseProvidedPrecursorInfo: UseProvidedPrecursorInfo, DeconvolutionIntensityRatio: DeconvolutionIntensityRatio, DeconvolutionMaxAssumedChargeState: DeconvolutionMaxAssumedChargeState, ReportAllAmbiguity: ReportAllAmbiguity, CompIons: CompIons , TotalPartitions: TotalPartitions , ScoreCutoff: ScoreCutoff, TopNpeaks: TopNpeaks, MinRatio: MinRatio, TrimMs1Peaks: TrimMs1Peaks, TrimMsMsPeaks: TrimMsMsPeaks, UseDeltaScore: UseDeltaScore, CalculateEValue: CalculateEValue, prodMassTol: prodMassTol, preMassTol: preMassTol , deconMassTol: deconMassTol, MaxThreadsToUsePerFile: MaxThreadsToUsePerFile , DigestionParams: DigestionParams, ListOfModsVariable: ListOfModsVariable , ListOfModsFixed: ListOfModsFixed);
        }

        private static DigestionParams GetDigestionParamsFromString(string value)
        {
            DigestionParams temp = new DigestionParams();
            Dictionary<string, string> parameters = new Dictionary<string, string>();
            string[] input = value.Split(';');
            
            for (int i = 0; i < input.Length - 1; i++)
            {

                string a = input[i];
                string[] temps = a.Split('=');
                string key = temps[0];
                string val = temps[1];
                parameters.Add(key, val);
            }
            string protease = parameters[nameof(temp.Protease)];
            int MaxMissedCleavages = Convert.ToInt32(parameters[nameof(temp.MaxMissedCleavages)]);
            int MinPeptideLength = Convert.ToInt32(parameters[nameof(temp.MinPeptideLength)]);
            int MaxPeptideLength = Convert.ToInt32(parameters[nameof(temp.MaxPeptideLength)]);
            int MaxModificationIsoforms = Convert.ToInt32(parameters[nameof(temp.MaxModificationIsoforms)]);
            InitiatorMethionineBehavior InitiatorMethionineBehavior;
            if (parameters[nameof(temp.InitiatorMethionineBehavior)] == "InitiatorMethionineBehavior.Cleave")
            {
                InitiatorMethionineBehavior = InitiatorMethionineBehavior.Cleave;
            }
            else if (parameters[nameof(temp.InitiatorMethionineBehavior)] == "InitiatorMethionineBehavior.Retain")
            {
                InitiatorMethionineBehavior = InitiatorMethionineBehavior.Retain;
            }
            else if (parameters[nameof(temp.InitiatorMethionineBehavior)] == "InitiatorMethionineBehavior.Variable")
            {
                InitiatorMethionineBehavior = InitiatorMethionineBehavior.Variable;
            }
            else
            {
                InitiatorMethionineBehavior = InitiatorMethionineBehavior.Undefined;
            }
            int MaxModsForPeptides = Convert.ToInt32(parameters[nameof(temp.MaxModsForPeptide)]);
            bool SemiProteaseDigestion = Convert.ToBoolean(parameters[nameof(temp.SemiProteaseDigestion)]);
            
            TerminusType TerminusTypeSemiProtease;
            if (parameters[nameof(temp.TerminusTypeSemiProtease)] == "TerminusType.N")
            {
                TerminusTypeSemiProtease = TerminusType.N;
            }
            else if (parameters[nameof(temp.TerminusTypeSemiProtease)] == "TerminusType.C")
            {
                TerminusTypeSemiProtease = TerminusType.C;
            }
            else
            {
                TerminusTypeSemiProtease = TerminusType.None;
            }

           return new DigestionParams(protease: protease, MaxMissedCleavages: MaxMissedCleavages, MinPeptideLength: MinPeptideLength, MaxPeptideLength: MaxPeptideLength, MaxModificationIsoforms: MaxModificationIsoforms, InitiatorMethionineBehavior: InitiatorMethionineBehavior, MaxModsForPeptides: MaxModsForPeptides, SemiProteaseDigestion: SemiProteaseDigestion, TerminusTypeSemiProtease: TerminusTypeSemiProtease);
        }

        private static List<string> GetModsTypesFromString(string value)
        {
            Console.WriteLine(value);
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

        #endregion Private Methods
    }
}