using EngineLayer;
using EngineLayer.Analysis;
using MathNet.Numerics.Distributions;
using MzLibUtil;
using Nett;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using UsefulProteomicsDatabases;

namespace TaskLayer
{
    public enum MyTask
    {
        Search,
        Gptmd,
        Calibrate
    }

    public abstract class MetaMorpheusTask
    {

        #region Public Fields

        public static readonly TomlConfig tomlConfig = TomlConfig.Create(cfg => cfg
                        .ConfigureType<Tolerance>(type => type
                            .WithConversionFor<TomlString>(convert => convert
                                .ToToml(custom => custom.ToString())
                                .FromToml(tmlString => new Tolerance(tmlString.Value))))
                        .ConfigureType<SearchMode>(type => type
                            .WithConversionFor<TomlString>(convert => convert
                                .ToToml(custom => custom.ToString())
                                .FromToml(tmlString => MetaMorpheusTask.ParseSearchMode(tmlString.Value))))
                        .ConfigureType<Protease>(type => type
                            .WithConversionFor<TomlString>(convert => convert
                                .ToToml(custom => custom.ToString())
                                .FromToml(tmlString => GlobalTaskLevelSettings.ProteaseDictionary[tmlString.Value])))
                        .ConfigureType<List<Tuple<string, string>>>(type => type
                             .WithConversionFor<TomlString>(convert => convert
                                 .ToToml(custom => string.Join("\t\t", custom.Select(b => b.Item1 + "\t" + b.Item2)))
                                 .FromToml(tmlString => GetModsFromString(tmlString.Value)))));

        private static List<Tuple<string, string>> GetModsFromString(string value)
        {
            return value.Split(new string[] { "\t\t" }, StringSplitOptions.None).Select(b => new Tuple<string, string>(b.Split('\t').First(), b.Split('\t').Last())).ToList();
        }

        #endregion Public Fields

        #region Protected Fields

        protected MyTaskResults myTaskResults;

        #endregion Protected Fields

        #region Public Constructors

        public MetaMorpheusTask(MyTask taskType)
        {
            this.TaskType = taskType;
        }

        #endregion Public Constructors

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

        #endregion Public Properties

        #region Public Methods

        public static SearchMode ParseSearchMode(string text)
        {
            SearchMode ye = null;

            var split = text.Split(' ');

            switch (split[1])
            {
                case "dot":
                    ToleranceUnit tu = ToleranceUnit.PPM;
                    if (split[3].ToUpperInvariant().Equals("PPM"))
                        tu = ToleranceUnit.PPM;
                    else if (split[3].ToUpperInvariant().Equals("DA"))
                        tu = ToleranceUnit.Absolute;
                    else
                        break;

                    var massShifts = Array.ConvertAll(split[4].Split(','), Double.Parse);
                    var newString = split[2].Replace("±", "");
                    var toleranceValue = double.Parse(newString, CultureInfo.InvariantCulture);
                    ye = new DotSearchMode(split[0], massShifts, new Tolerance(tu, toleranceValue));
                    break;

                case "interval":
                    IEnumerable<DoubleRange> doubleRanges = Array.ConvertAll(split[2].Split(','), b => new DoubleRange(double.Parse(b.Trim(new char[] { '[', ']' }).Split(';')[0], CultureInfo.InvariantCulture), double.Parse(b.Trim(new char[] { '[', ']' }).Split(';')[1], CultureInfo.InvariantCulture)));
                    ye = new IntervalSearchMode(split[0], doubleRanges);
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
                    throw new Exception("Could not parse search mode string");
            }
            return ye;
        }

        public MyTaskResults RunTask(string output_folder, List<DbForTask> currentXmlDbFilenameList, List<string> currentRawDataFilenameList, string taskId)
        {
            StartingSingleTask(taskId);
            var paramsFileName = Path.Combine(output_folder, "prose.txt");
            using (StreamWriter file = new StreamWriter(paramsFileName))
            {
                file.WriteLine("MetaMorpheus version "
                    + (GlobalEngineLevelSettings.MetaMorpheusVersion.Equals("1.0.0.0") ? "NOT A RELEASE" : GlobalEngineLevelSettings.MetaMorpheusVersion)
                    + " is used to run a "
                    + this.TaskType
                    + " task on "
                    + currentRawDataFilenameList.Count
                    + " spectra files.");

                file.WriteLine(ToString());

                file.WriteLine();
                file.WriteLine("taskId: " + taskId);
                file.WriteLine("Spectra files:");
                file.WriteLine(string.Join(Environment.NewLine, currentRawDataFilenameList.Select(b => '\t' + b)));
                file.WriteLine("XML files:");
                file.Write(string.Join(Environment.NewLine, currentXmlDbFilenameList.Select(b => '\t' + (b.IsContaminant ? "Contaminant " : "") + b.FileName)));
            }
            SucessfullyFinishedWritingFile(paramsFileName, new List<string> { taskId });

            // TOML
            var tomlFileName = Path.Combine(output_folder, GetType().Name + "config.toml");

            Toml.WriteFile(this, tomlFileName, tomlConfig);
            SucessfullyFinishedWritingFile(tomlFileName, new List<string> { taskId });

            MetaMorpheusEngine.FinishedSingleEngineHandler += SingleEngineHandlerInTask;

#if !DEBUG
            try
            {
#endif
            var stopWatch = new Stopwatch();
            stopWatch.Start();
            RunSpecific(output_folder, currentXmlDbFilenameList, currentRawDataFilenameList, taskId);
            stopWatch.Stop();
            myTaskResults.Time = stopWatch.Elapsed;
            var resultsFileName = Path.Combine(output_folder, "results.txt");
            using (StreamWriter file = new StreamWriter(resultsFileName))
            {
                file.WriteLine(GlobalEngineLevelSettings.MetaMorpheusVersion.Equals("1.0.0.0") ? "MetaMorpheus: Not a release version" : "MetaMorpheus: version " + GlobalEngineLevelSettings.MetaMorpheusVersion);
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

            MetaMorpheusEngine.FinishedSingleEngineHandler -= SingleEngineHandlerInTask;
            return myTaskResults;
        }

        #endregion Public Methods

        #region Protected Internal Methods

        protected internal void WritePsmsToTsv(List<NewPsmWithFdr> items, string outputFolder, string fileName, List<string> nestedIds)
        {
            var writtenFile = Path.Combine(outputFolder, fileName + ".psmtsv");
            using (StreamWriter output = new StreamWriter(writtenFile))
            {
                output.WriteLine(NewPsmWithFdr.TabSeparatedHeader);
                for (int i = 0; i < items.Count; i++)
                    output.WriteLine(items[i]);
            }
            SucessfullyFinishedWritingFile(writtenFile, nestedIds);
        }

        #endregion Protected Internal Methods

        #region Protected Methods

        protected static List<Protein> LoadProteinDb(string fileName, bool generateDecoys, List<ModificationWithMass> localizeableModifications, bool isContaminant, out Dictionary<string, Modification> um)
        {
            if (Path.GetExtension(fileName).Equals(".fasta"))
            {
                um = null;
                return ProteinDbLoader.LoadProteinFasta(fileName, generateDecoys, isContaminant, ProteinDbLoader.uniprot_accession_expression, ProteinDbLoader.uniprot_fullName_expression, ProteinDbLoader.uniprot_fullName_expression, ProteinDbLoader.uniprot_gene_expression);
            }
            else
                return ProteinDbLoader.LoadProteinXML(fileName, generateDecoys, localizeableModifications, isContaminant, null, out um);
        }

        protected void ReportProgress(ProgressEventArgs v)
        {
            OutProgressHandler?.Invoke(this, v);
        }

        protected abstract MyTaskResults RunSpecific(string output_folder, List<DbForTask> currentXmlDbFilenameList, List<string> currentRawDataFilenameList, string taskId);

        protected void WriteProteinGroupsToTsv(List<ProteinGroup> items, string outputFolder, string fileName, List<string> nestedIds)
        {
            if (items != null)
            {
                var writtenFile = Path.Combine(outputFolder, fileName + ".tsv");

                using (StreamWriter output = new StreamWriter(writtenFile))
                {
                    output.WriteLine(items.First().TabSeparatedHeader);
                    for (int i = 0; i < items.Count; i++)
                        output.WriteLine(items[i]);
                }

                SucessfullyFinishedWritingFile(writtenFile, nestedIds);
            }
        }

        protected void WriteTree(BinTreeStructure myTreeStructure, string output_folder, string fileName, List<string> nestedIds)
        {
            var writtenFile = Path.Combine(output_folder, fileName + ".mytsv");
            using (StreamWriter output = new StreamWriter(writtenFile))
            {
                output.WriteLine("MassShift\tCount\tCountDecoy\tCountTarget\tCountLocalizeableTarget\tCountNonLocalizeableTarget\tFDR\tArea 0.01t\tArea 0.255\tFracLocalizeableTarget\tMine\tUnimodID\tUnimodFormulas\tUnimodDiffs\tAA\tCombos\tModsInCommon\tAAsInCommon\tResidues\tprotNtermLocFrac\tpepNtermLocFrac\tpepCtermLocFrac\tprotCtermLocFrac\tFracWithSingle\tOverlappingFrac\tMedianLength\tUniprot");
                foreach (Bin bin in myTreeStructure.FinalBins.OrderByDescending(b => b.Count))
                {
                    output.WriteLine(bin.MassShift.ToString("F4", CultureInfo.InvariantCulture)
                        + "\t" + bin.Count.ToString(CultureInfo.InvariantCulture)
                        + "\t" + bin.CountDecoy.ToString(CultureInfo.InvariantCulture)
                        + "\t" + bin.CountTarget.ToString(CultureInfo.InvariantCulture)
                        + "\t" + bin.LocalizeableTarget.ToString(CultureInfo.InvariantCulture)
                        + "\t" + (bin.CountTarget - bin.LocalizeableTarget).ToString(CultureInfo.InvariantCulture)
                        + "\t" + (bin.Count == 0 ? double.NaN : (double)bin.CountDecoy / bin.Count).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (Normal.CDF(0, 1, bin.ComputeZ(0.01))).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (Normal.CDF(0, 1, bin.ComputeZ(0.255))).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (bin.CountTarget == 0 ? double.NaN : (double)bin.LocalizeableTarget / bin.CountTarget).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + bin.Mine
                        + "\t" + bin.UnimodId
                        + "\t" + bin.UnimodFormulas
                        + "\t" + bin.UnimodDiffs
                        + "\t" + bin.AA
                        + "\t" + bin.combos
                        + "\t" + string.Join(",", bin.modsInCommon.OrderByDescending(b => b.Value).Where(b => b.Value > bin.CountTarget / 10.0).Select(b => b.Key + ":" + ((double)b.Value / bin.CountTarget).ToString("F3", CultureInfo.InvariantCulture)))
                        + "\t" + string.Join(",", bin.AAsInCommon.OrderByDescending(b => b.Value).Where(b => b.Value > bin.CountTarget / 10.0).Select(b => b.Key + ":" + ((double)b.Value / bin.CountTarget).ToString("F3", CultureInfo.InvariantCulture)))
                        + "\t" + string.Join(",", bin.residueCount.OrderByDescending(b => b.Value).Select(b => b.Key + ":" + b.Value))
                        + "\t" + (bin.LocalizeableTarget == 0 ? double.NaN : (double)bin.protNlocCount / bin.LocalizeableTarget).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (bin.LocalizeableTarget == 0 ? double.NaN : (double)bin.pepNlocCount / bin.LocalizeableTarget).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (bin.LocalizeableTarget == 0 ? double.NaN : (double)bin.pepClocCount / bin.LocalizeableTarget).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (bin.LocalizeableTarget == 0 ? double.NaN : (double)bin.protClocCount / bin.LocalizeableTarget).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (bin.FracWithSingle).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + ((double)bin.Overlapping / bin.CountTarget).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (bin.MedianLength).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + bin.uniprotID);
                }
            }
            SucessfullyFinishedWritingFile(writtenFile, nestedIds);
        }

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