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
                                .FromToml(tmlString => GetDigestionParamsFromString(tmlString.Value)/*ProteaseDictionary[tmlString.Value]*/)))
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
            CommonParameters returnParams = new CommonParameters(DigestionParams: fileSpecificDigestionParams, preMassTol: preMassTol, prodMassTol: prodMassTol, BIons: BIons, YIons:YIons, CIons:CIons, ZdotIons:ZdotIons);
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
            string[] input = value.Split(';');
            string[,] output = new string[input.Length, 2];
            for (int i = 0; i < input.Length-1; i++)
            {
                
                string a = input[i];
                string[] temp = a.Split('=');
                string key = temp[0];
                string val = temp[1];
               
                output[i, 0] = key;
                output[i, 1] = val;
            }

                            
            bool BIons = Convert.ToBoolean(output[0, 1]);
            bool YIons = Convert.ToBoolean(output[1, 1]);
            bool ZdotIons = Convert.ToBoolean(output[2, 1]);
            bool CIons = Convert.ToBoolean(output[3, 1]);
            bool DoPrecursorDeconvolution = Convert.ToBoolean(output[4, 1]);
            bool UseProvidedPrecursorInfo = Convert.ToBoolean(output[5, 1]);
            double DeconvolutionIntensityRatio = Convert.ToDouble(output[6, 1]);
            int DeconvolutionMaxAssumedChargeState = Convert.ToInt32(output[7, 1]);
            bool ReportAllAmbiguity = Convert.ToBoolean(output[8, 1]);
            bool CompIons = Convert.ToBoolean(output[9, 1]);
            int TotalPartitions = Convert.ToInt32(output[10, 1]);
            double ScoreCutoff = Convert.ToDouble(output[11, 1]);
            int TopNpeaks = Convert.ToInt32(output[12, 1]);
            double MinRatio = Convert.ToDouble(output[13, 1]);
            bool TrimMs1Peaks = Convert.ToBoolean(output[14, 1]);
            bool TrimMsMsPeaks = Convert.ToBoolean(output[15, 1]);
            bool UseDeltaScore = Convert.ToBoolean(output[16, 1]);
            bool CalculateEValue = Convert.ToBoolean(output[17, 1]);
            double prodMassTol = Convert.ToDouble(output[18, 1]);
            double preMassTol = Convert.ToDouble(output[19, 1]);
            double deconMassTol = Convert.ToDouble(output[20, 1]);
            int MaxThreadsToUsePerFile = Convert.ToInt32(output[21, 1]);

            string protease = output[22, 1];
            int MaxMissedCleavages = Convert.ToInt32(output[23, 1]);
            int MinPeptideLength = Convert.ToInt32(output[24, 1]);
            int MaxPeptideLength = Convert.ToInt32(output[25, 1]);
            int MaxModificationIsoforms = Convert.ToInt32(output[26, 1]);
            InitiatorMethionineBehavior InitiatorMethionineBehavior;
            if (output[27, 1] == "InitiatorMethionineBehavior.Cleave")
            {
                InitiatorMethionineBehavior = InitiatorMethionineBehavior.Cleave;
            }
            else if (output[27, 1] == "InitiatorMethionineBehavior.Retain")
            {
                InitiatorMethionineBehavior = InitiatorMethionineBehavior.Retain;
            }
            else if (output[27, 1] == "InitiatorMethionineBehavior.Variable")
            {
                InitiatorMethionineBehavior = InitiatorMethionineBehavior.Variable;
            }
            else
            {
                InitiatorMethionineBehavior = InitiatorMethionineBehavior.Undefined;
            }
            int MaxModsForPeptides = Convert.ToInt32(output[28, 1]);
            bool SemiProteaseDigestion = Convert.ToBoolean(output[29, 1]);
            //string term = ;
            TerminusType TerminusTypeSemiProtease;
            if (output[30, 1] == "TerminusType.N")
            {
                TerminusTypeSemiProtease = TerminusType.N;
            }
            else if (output[30, 1] == "TerminusType.C")
            {
                TerminusTypeSemiProtease = TerminusType.C;
            }
            else
            {
                TerminusTypeSemiProtease = TerminusType.None;
            }
            
            DigestionParams DigestionParams = new DigestionParams(protease: protease, MaxMissedCleavages: MaxMissedCleavages, MinPeptideLength: MinPeptideLength, MaxPeptideLength: MaxPeptideLength, MaxModificationIsoforms: MaxModificationIsoforms, InitiatorMethionineBehavior: InitiatorMethionineBehavior, MaxModsForPeptides: MaxModsForPeptides, SemiProteaseDigestion: SemiProteaseDigestion, TerminusTypeSemiProtease:TerminusTypeSemiProtease);
            string listofModsVariable = output[31, 1];
            string[] listContentVariable = listofModsVariable.Split(',');
            string[,] listContentStringsVariable = new string[listContentVariable.Length, 2];
            List<(string, string)> ListOfModsVariable = new List<(string, string)>();
            for (int i = 0; i < listContentVariable.Length-1; i++)
            {

                string a = listContentVariable[i];
                string[] temp = a.Split('/');
                string key = temp[0];
                string val = temp[1];

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


            string listofModsFixed = output[32, 1];
            string[] listContentFixed = listofModsFixed.Split(',');
            string[,] listContentStringsFixed = new string[listContentFixed.Length, 2];
            List<(string, string)> ListOfModsFixed = new List<(string, string)>();

            for (int i = 0; i < listContentFixed.Length-1; i++)
            {

                string a = listContentFixed[i];
                string[] temp = a.Split('/');
                string key = temp[0];
                string val = temp[1];

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
            string[] input = value.Split(';');
            string[,] output = new string[input.Length, 2];
            for (int i = 0; i < input.Length - 1; i++)
            {

                string a = input[i];
                string[] temp = a.Split('=');
                string key = temp[0];
                string val = temp[1];

                output[i, 0] = key;
                output[i, 1] = val;
            }
            string protease = output[21, 1];
            int MaxMissedCleavages = Convert.ToInt32(output[22, 1]);
            int MinPeptideLength = Convert.ToInt32(output[23, 1]);
            int MaxPeptideLength = Convert.ToInt32(output[24, 1]);
            int MaxModificationIsoforms = Convert.ToInt32(output[25, 1]);
            InitiatorMethionineBehavior InitiatorMethionineBehavior;
            if (output[26, 1] == "InitiatorMethionineBehavior.Cleave")
            {
                InitiatorMethionineBehavior = InitiatorMethionineBehavior.Cleave;
            }
            else if (output[26, 1] == "InitiatorMethionineBehavior.Retain")
            {
                InitiatorMethionineBehavior = InitiatorMethionineBehavior.Retain;
            }
            else if (output[26, 1] == "InitiatorMethionineBehavior.Variable")
            {
                InitiatorMethionineBehavior = InitiatorMethionineBehavior.Variable;
            }
            else
            {
                InitiatorMethionineBehavior = InitiatorMethionineBehavior.Undefined;
            }
            int MaxModsForPeptides = Convert.ToInt32(output[27, 1]);
            bool SemiProteaseDigestion = Convert.ToBoolean(output[28, 1]);
            
            TerminusType TerminusTypeSemiProtease;
            if (output[29, 1] == "TerminusType.N")
            {
                TerminusTypeSemiProtease = TerminusType.N;
            }
            else if (output[29, 1] == "TerminusType.C")
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