using EngineLayer;
using EngineLayer.Analysis;
using MathNet.Numerics.Distributions;
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
        Calibrate
    }

    public abstract class MetaMorpheusTask
    {

        #region Public Events

        public static event EventHandler<SingleTaskEventArgs> FinishedSingleTaskHandler;

        public static event EventHandler<SingleFileEventArgs> FinishedWritingFileHandler;

        public static event EventHandler<SingleTaskEventArgs> StartingSingleTaskHander;

        public static event EventHandler<StringEventArgs> StartingDataFileHandler;

        public static event EventHandler<StringEventArgs> FinishedDataFileHandler;

        public static event EventHandler<StringEventArgs> OutLabelStatusHandler;

        public static event EventHandler<StringEventArgs> WarnHandler;

        public static event EventHandler<StringEventArgs> NewCollectionHandler;

        #endregion Public Events

        #region Public Properties

        public static List<ModList> AllModLists { get; private set; }

        public MyTask TaskType { get; internal set; }

        public InitiatorMethionineBehavior InitiatorMethionineBehavior { get; set; }

        public int MaxMissedCleavages { get; set; }

        public int MaxModificationIsoforms { get; set; }

        public Protease Protease { get; set; }

        public bool BIons { get; set; }

        public bool YIons { get; set; }

        public bool ZdotIons { get; set; }

        public bool CIons { get; set; }

        public int MaxNumPeaksPerScan { get; set; }

        #endregion Public Properties

        #region Protected Properties

        protected abstract string SpecificTaskInfo { get; }

        #endregion Protected Properties

        #region Public Methods

        public static void AddModList(ModList modList)
        {
            if (AllModLists == null)
                AllModLists = new List<ModList>();
            AllModLists.Add(modList);
        }

        public MyTaskResults RunTask(string output_folder, List<DbForTask> currentXmlDbFilenameList, List<string> currentRawDataFilenameList, string taskId)
        {
            startingSingleTask(taskId);
            var paramsFileName = Path.Combine(output_folder, "params.txt");
            using (StreamWriter file = new StreamWriter(paramsFileName))
            {
                file.WriteLine(MyEngine.MetaMorpheusVersion.Equals("1.0.0.0") ? "MetaMorpheus: Not a release version" : "MetaMorpheus: version " + MyEngine.MetaMorpheusVersion);
                file.Write(ToString());
            }
            SucessfullyFinishedWritingFile(paramsFileName, new List<string> { taskId });
#if !DEBUG
            try
            {
#endif
            var stopWatch = new Stopwatch();
            stopWatch.Start();
            var myResults = RunSpecific(output_folder, currentXmlDbFilenameList, currentRawDataFilenameList, taskId);
            stopWatch.Stop();
            myResults.Time = stopWatch.Elapsed;
            var resultsFileName = Path.Combine(output_folder, "results.txt");
            using (StreamWriter file = new StreamWriter(resultsFileName))
            {
                file.WriteLine(MyEngine.MetaMorpheusVersion.Equals("1.0.0.0") ? "MetaMorpheus: Not a release version" : "MetaMorpheus: version " + MyEngine.MetaMorpheusVersion);
                file.Write(myResults.ToString());
            }
            SucessfullyFinishedWritingFile(resultsFileName, new List<string> { taskId });
            finishedSingleTask(taskId);
            return myResults;
#if !DEBUG
            }
            catch (Exception e)
            {
                var resultsFileName = Path.Combine(output_folder, "results.txt");
                using (StreamWriter file = new StreamWriter(resultsFileName))
                {
                    file.WriteLine(MyEngine.MetaMorpheusVersion.Equals("1.0.0.0") ? "MetaMorpheus: Not a release version" : "MetaMorpheus: version " + MyEngine.MetaMorpheusVersion);
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
        }

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine(SpecificTaskInfo);
            sb.AppendLine(TaskType.ToString());
            sb.AppendLine("Spectra files:");
            //sb.AppendLine(string.Join(Environment.NewLine, rawDataFilenameList.Select(b => '\t' + b)));
            sb.AppendLine("XML files:");
            //sb.AppendLine(string.Join(Environment.NewLine, dbFilenameList.Select(b => '\t' + (b.IsContaminant ? "Contaminant " : "") + b.FileName)));
            sb.AppendLine("initiatorMethionineBehavior: " + InitiatorMethionineBehavior);
            sb.AppendLine("maxMissedCleavages: " + MaxMissedCleavages);
            sb.AppendLine("maxModificationIsoforms: " + MaxModificationIsoforms);
            //sb.AppendLine("output_folder: " + OutputFolder);
            sb.AppendLine("protease: " + Protease);
            sb.AppendLine("bIons: " + BIons);
            sb.AppendLine("yIons: " + YIons);
            sb.AppendLine("cIons: " + CIons);
            sb.AppendLine("zdotIons: " + ZdotIons);
            return sb.ToString();
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
                return ProteinDbLoader.LoadProteinFasta(fileName, generateDecoys, isContaminant, ProteinDbLoader.uniprot_accession_expression, ProteinDbLoader.uniprot_fullName_expression, ProteinDbLoader.uniprot_fullName_expression);
            }
            else
                return ProteinDbLoader.LoadProteinXML(fileName, generateDecoys, localizeableModifications, isContaminant, null, out um);
        }

        protected abstract MyTaskResults RunSpecific(string output_folder, List<DbForTask> currentXmlDbFilenameList, List<string> currentRawDataFilenameList, string taskId);

        protected void WriteProteinGroupsToTsv(List<ProteinGroup> items, string outputFolder, string fileName, List<string> nestedIds)
        {
            if (items != null)
            {
                var writtenFile = Path.Combine(outputFolder, fileName + ".tsv");

                using (StreamWriter output = new StreamWriter(writtenFile))
                {
                    output.WriteLine(ProteinGroup.TabSeparatedHeader);
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
                output.WriteLine("MassShift\tCount\tCountDecoy\tCountTarget\tCountLocalizeableTarget\tCountNonLocalizeableTarget\tFDR\tArea 0.01t\tArea 0.255\tFracLocalizeableTarget\tMine\tUnimodID\tUnimodFormulas\tAA\tCombos\tModsInCommon\tAAsInCommon\tResidues\tNtermLocFrac\tCtermLocFrac\tFracWithSingle\tOverlappingFrac\tMedianLength\tUniprot");
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
                        + "\t" + bin.AA
                        + "\t" + bin.combos
                        + "\t" + string.Join(",", bin.modsInCommon.OrderByDescending(b => b.Value).Where(b => b.Value > bin.CountTarget / 10.0).Select(b => b.Key + ":" + ((double)b.Value / bin.CountTarget).ToString("F3", CultureInfo.InvariantCulture)))
                        + "\t" + string.Join(",", bin.AAsInCommon.OrderByDescending(b => b.Value).Where(b => b.Value > bin.CountTarget / 10.0).Select(b => b.Key + ":" + ((double)b.Value / bin.CountTarget).ToString("F3", CultureInfo.InvariantCulture)))
                        + "\t" + string.Join(",", bin.residueCount.OrderByDescending(b => b.Value).Select(b => b.Key + ":" + b.Value))
                        + "\t" + (bin.LocalizeableTarget == 0 ? double.NaN : (double)bin.NlocCount / bin.LocalizeableTarget).ToString("F3", CultureInfo.InvariantCulture)
                        + "\t" + (bin.LocalizeableTarget == 0 ? double.NaN : (double)bin.ClocCount / bin.LocalizeableTarget).ToString("F3", CultureInfo.InvariantCulture)
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

        private void finishedSingleTask(string taskId)
        {
            FinishedSingleTaskHandler?.Invoke(this, new SingleTaskEventArgs(taskId));
        }

        private void startingSingleTask(string taskId)
        {
            StartingSingleTaskHander?.Invoke(this, new SingleTaskEventArgs(taskId));
        }

        #endregion Private Methods

    }
}