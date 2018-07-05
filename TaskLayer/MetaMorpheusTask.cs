using Chemistry;
using EngineLayer;
using MassSpectrometry;
using MzLibUtil;
using Nett;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;
using UsefulProteomicsDatabases;

namespace TaskLayer
{
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

        protected MetaMorpheusTask(TaskType taskType)
        {
            TaskType = taskType;
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

        public TaskType TaskType { get; set; }
        public CommonParameters CommonParameters { get; set; }

        public static IEnumerable<Ms2ScanWithSpecificMass> GetMs2Scans(
            MsDataFile myMSDataFile,
            string fullFilePath,
            bool doPrecursorDeconvolution,
            bool useProvidedPrecursorInfo,
            double deconvolutionIntensityRatio,
            int deconvolutionMaxAssumedChargeState,
            Tolerance deconvolutionMassTolerance)
        {
            foreach (var ms2scan in myMSDataFile.GetAllScansList().Where(x => x.MsnOrder != 1))
            {
                if (GlobalVariables.StopLoops) { break; }

                List<(double, int)> isolatedStuff = new List<(double, int)>();
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

                    if (doPrecursorDeconvolution)
                    {
                        foreach (var envelope in ms2scan.GetIsolatedMassesAndCharges(precursorSpectrum.MassSpectrum, 1, deconvolutionMaxAssumedChargeState, deconvolutionMassTolerance.Value, deconvolutionIntensityRatio))
                        {
                            var monoPeakMz = envelope.monoisotopicMass.ToMz(envelope.charge);
                            isolatedStuff.Add((monoPeakMz, envelope.charge));
                        }
                    }
                }

                if (useProvidedPrecursorInfo && ms2scan.SelectedIonChargeStateGuess.HasValue)
                {
                    var precursorCharge = ms2scan.SelectedIonChargeStateGuess.Value;
                    if (ms2scan.SelectedIonMonoisotopicGuessMz.HasValue)
                    {
                        var precursorMZ = ms2scan.SelectedIonMonoisotopicGuessMz.Value;
                        if (!isolatedStuff.Any(b => deconvolutionMassTolerance.Within(precursorMZ.ToMass(precursorCharge), b.Item1.ToMass(b.Item2))))
                        {
                            isolatedStuff.Add((precursorMZ, precursorCharge));
                        }
                    }
                    else
                    {
                        var precursorMZ = ms2scan.SelectedIonMZ;
                        if (!isolatedStuff.Any(b => deconvolutionMassTolerance.Within(precursorMZ.Value.ToMass(precursorCharge), b.Item1.ToMass(b.Item2))))
                        {
                            isolatedStuff.Add((precursorMZ.Value, precursorCharge));
                        }
                    }
                }

                foreach (var heh in isolatedStuff)
                {
                    yield return new Ms2ScanWithSpecificMass(ms2scan, heh.Item1, heh.Item2, fullFilePath);
                }
            }
        }

        public static CommonParameters SetAllFileSpecificCommonParams(CommonParameters commonParams, FileSpecificParameters fileSpecificParams)
        {
            if (fileSpecificParams == null)
            {
                return commonParams;
            }

            // set file-specific digestion parameters
            Protease protease = fileSpecificParams.Protease ?? commonParams.DigestionParams.Protease;
            int MinPeptideLength = fileSpecificParams.MinPeptideLength ?? commonParams.DigestionParams.MinPeptideLength;
            int MaxPeptideLength = fileSpecificParams.MaxPeptideLength ?? commonParams.DigestionParams.MaxPeptideLength;
            int MaxMissedCleavages = fileSpecificParams.MaxMissedCleavages ?? commonParams.DigestionParams.MaxMissedCleavages;
            int MaxModsForPeptide = fileSpecificParams.MaxModsForPeptide ?? commonParams.DigestionParams.MaxModsForPeptide;
            DigestionParams fileSpecificDigestionParams = new DigestionParams(
                protease: protease.Name,
                maxMissedCleavages: MaxMissedCleavages,
                minPeptideLength: MinPeptideLength,
                maxPeptideLength: MaxPeptideLength,
                maxModsForPeptides: MaxModsForPeptide);

            // set the rest of the file-specific parameters
            Tolerance precursorMassTolerance = fileSpecificParams.PrecursorMassTolerance ?? commonParams.PrecursorMassTolerance;
            Tolerance productMassTolerance = fileSpecificParams.ProductMassTolerance ?? commonParams.ProductMassTolerance;
            bool bIons = fileSpecificParams.BIons ?? commonParams.BIons;
            bool yIons = fileSpecificParams.YIons ?? commonParams.YIons;
            bool cIons = fileSpecificParams.CIons ?? commonParams.CIons;
            bool zdotIons = fileSpecificParams.ZdotIons ?? commonParams.ZdotIons;

            CommonParameters returnParams = new CommonParameters(
                bIons: bIons,
                yIons: yIons,
                cIons: cIons,
                zDotIons: zdotIons,
                precursorMassTolerance: precursorMassTolerance,
                productMassTolerance: productMassTolerance,
                digestionParams: fileSpecificDigestionParams);

            return returnParams;
        }

        public MyTaskResults RunTask(string outputFolder, List<DbForTask> currentProteinDbFilenameList, List<string> currentRawDataFilepathList, string displayName)
        {
            StartingSingleTask(displayName);

            // write TOML
            var tomlFileName = Path.Combine(outputFolder, GetType().Name + "config.toml");
            Toml.WriteFile(this, tomlFileName, tomlConfig);
            SucessfullyFinishedWritingFile(tomlFileName, new List<string> { displayName });

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
                            // file-specific toml has already been validated in the GUI when the spectra files were added, so...
                            // probably the only time you can get here is if the user modifies the file-specific parameter file in the middle of a run...
                            Warn("Problem parsing the file-specific toml " + Path.GetFileName(fileSpecificTomlPath) + "; " + e.Message + "; is the toml from an older version of MetaMorpheus?");
                        }
                    }
                }

                RunSpecific(outputFolder, currentProteinDbFilenameList, currentRawDataFilepathList, displayName, fileSettingsList);
                stopWatch.Stop();
                MyTaskResults.Time = stopWatch.Elapsed;
                var resultsFileName = Path.Combine(outputFolder, "results.txt");
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
                var resultsFileName = Path.Combine(outputFolder, "results.txt");
                e.Data.Add("folder", outputFolder);
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

            // Write prose
            var proseFilePath = Path.Combine(outputFolder, "prose.txt");
            using (StreamWriter file = new StreamWriter(proseFilePath))
            {
                file.Write("The data analysis was performed using MetaMorpheus version " + GlobalVariables.MetaMorpheusVersion + ", available at " + "https://github.com/smith-chem-wisc/MetaMorpheus." + " [INSERT CITATION] ");
                file.Write(ProseCreatedWhileRunning.ToString());
                file.Write(SystemInfo.SystemProse().Replace(Environment.NewLine, "") + " ");
                file.WriteLine("The total time to perform the " + TaskType.ToString().ToLowerInvariant() + " task on " + currentRawDataFilepathList.Count + " spectra file(s) was " + String.Format("{0:0.00}", MyTaskResults.Time.TotalMinutes) + " minutes.");
                file.WriteLine();
                file.WriteLine("Published works using MetaMorpheus software are encouraged to cite: STEFAN'S VERY IMPORTANT PAPER");

                file.WriteLine();
                file.WriteLine("Spectra files: ");
                file.WriteLine(string.Join(Environment.NewLine, currentRawDataFilepathList.Select(b => '\t' + b)));
                file.WriteLine("Databases:");
                file.Write(string.Join(Environment.NewLine, currentProteinDbFilenameList.Select(b => '\t' + (b.IsContaminant ? "Contaminant " : "") + b.FilePath)));
            }
            SucessfullyFinishedWritingFile(proseFilePath, new List<string> { displayName });

            MetaMorpheusEngine.FinishedSingleEngineHandler -= SingleEngineHandlerInTask;
            return MyTaskResults;
        }

        /// <summary>
        /// Loads target proteins (does not generate decoy proteins)
        /// </summary>
        /// <param name="taskId"></param>
        /// <param name="dbFilenameList"></param>
        /// <param name="localizeableModificationTypes"></param>
        /// <returns></returns>
        protected List<Protein> LoadProteins(string taskId, List<DbForTask> dbFilenameList, bool searchTarget, DecoyType decoyType, List<string> localizeableModificationTypes)
        {
            Status("Loading proteins...", new List<string> { taskId });
            int emptyProteinEntries = 0;
            List<Protein> proteinList = new List<Protein>();
            foreach (var db in dbFilenameList)
            {
                int emptyProteinEntriesForThisDb = 0;
                var dbProteinList = LoadProteinDb(db.FilePath, searchTarget, decoyType, localizeableModificationTypes, db.IsContaminant, out Dictionary<string, Modification> unknownModifications, out emptyProteinEntriesForThisDb);
                dbProteinList = dbProteinList.SelectMany(p => p.GetVariantProteins().OfType<Protein>()).ToList();
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

        protected static void WritePsmsToTsv(IEnumerable<PeptideSpectralMatch> items, string filePath,
            IReadOnlyDictionary<string, int> ModstoWritePruned)
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

        private static List<Protein> LoadProteinDb(string fileName, bool generateTargets, DecoyType decoyType,
            List<string> localizeableModificationTypes, bool isContaminant,
            out Dictionary<string, Modification> um, out int emptyEntriesCount)
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
    }
}