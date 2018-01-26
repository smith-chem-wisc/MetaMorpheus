using EngineLayer;
using EngineLayer.Neo;
using MassSpectrometry;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using UsefulProteomicsDatabases;

namespace TaskLayer
{
    public class NeoSearchTask : MetaMorpheusTask
    {

        public bool AggregateTargetDecoyFiles { get; set; }
        public bool GenerateSplicedPeptides { get; set; }
        public bool AggregateNormalSplicedFiles { get; set; }

        public NeoSearchTask() : base(MyTask.Neo)
        {
            NeoParameters = new NeoParameters();

            CommonParameters = new CommonParameters
            {
                DoPrecursorDeconvolution = false,
                PrecursorMassTolerance = null,
                ProductMassTolerance = null
            };
            CommonParameters.DigestionParams.MinPeptideLength = 8;
            CommonParameters.DigestionParams.MaxPeptideLength = 13;
            CommonParameters.DigestionParams.Protease = GlobalEngineLevelSettings.ProteaseDictionary["non-specific"];
            CommonParameters.DigestionParams.MaxMissedCleavages = 12;
        }

        #region Public Properties

        public NeoParameters NeoParameters { get; set; }

        protected override MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId, FileSpecificSettings[] fileSettingsList)
        {
            myTaskResults = new MyTaskResults(this);

            if (AggregateTargetDecoyFiles)
            {
                //getfolders
                if (NeoParameters.DecoyFilePath == null)
                {
                    NeoParameters.DecoyFilePath = new DirectoryInfo(OutputFolder).Name;
                    string taskString = NeoParameters.DecoyFilePath.Split('-')[0];
                    int taskNum = Convert.ToInt32(taskString.Substring(4, taskString.Length - 4));
                    taskNum--;
                    NeoParameters.DecoyFilePath = OutputFolder.Substring(0, OutputFolder.Length - NeoParameters.DecoyFilePath.Length) + "Task" + taskNum + "-SearchTask\\" + Path.GetFileNameWithoutExtension(currentRawFileList[0]) + "_PSMs.psmtsv";
                    if (NeoParameters.TargetFilePath == null)
                    {
                        NeoParameters.TargetFilePath = new DirectoryInfo(OutputFolder).Name;
                        taskNum--;
                        NeoParameters.TargetFilePath = OutputFolder.Substring(0, OutputFolder.Length - NeoParameters.TargetFilePath.Length) + "Task" + taskNum + "-SearchTask\\" + Path.GetFileNameWithoutExtension(currentRawFileList[0]) + "_PSMs.psmtsv";
                    }
                }
                if (NeoParameters.TargetFilePath == null)
                {
                    NeoParameters.TargetFilePath = new DirectoryInfo(OutputFolder).Name;
                    string taskString = NeoParameters.TargetFilePath.Split('-')[0];
                    int taskNum = Convert.ToInt32(taskString.Substring(4, taskString.Length - 4));
                    taskNum--;
                    NeoParameters.TargetFilePath = OutputFolder.Substring(0, OutputFolder.Length - NeoParameters.TargetFilePath.Length) + "Task" + taskNum + "-SearchTask\\" + Path.GetFileNameWithoutExtension(currentRawFileList[0]) + "_PSMs.psmtsv";
                }
                AggregateSearchFiles.Combine(NeoParameters.TargetFilePath, NeoParameters.DecoyFilePath, OutputFolder+"\\"+ Path.GetFileNameWithoutExtension(currentRawFileList[0]));
            }
            else if (AggregateNormalSplicedFiles)
            {

                string cisPath = new DirectoryInfo(OutputFolder).Name;
                string taskString = cisPath.Split('-')[0];
                int taskNum = Convert.ToInt32(taskString.Substring(4, taskString.Length - 4));
                taskNum-=2;
                string transPath = OutputFolder.Substring(0, OutputFolder.Length - cisPath.Length) + "Task" + (taskNum+1) + "-SearchTask\\" + Path.GetFileNameWithoutExtension(currentRawFileList[0]) + "_PSMs.psmtsv";
                cisPath = OutputFolder.Substring(0, OutputFolder.Length - cisPath.Length) + "Task" + taskNum + "-SearchTask\\" + Path.GetFileNameWithoutExtension(currentRawFileList[0]) + "_PSMs.psmtsv";

            //    AggregateSearchFiles.Combine(NeoParameters.TargetFilePath, cisPath);
              //  AggregateSearchFiles.Combine(NeoParameters.TargetFilePath, transPath);
            }
            else
            {
                NeoMassCalculator.ImportMasses();

                ParallelOptions parallelOptions = new ParallelOptions();
                if (CommonParameters.MaxParallelFilesToAnalyze.HasValue)
                    parallelOptions.MaxDegreeOfParallelism = CommonParameters.MaxParallelFilesToAnalyze.Value;
                MyFileManager myFileManager = new MyFileManager(true);

                //Import Spectra
                Parallel.For(0, currentRawFileList.Count, parallelOptions, spectraFileIndex =>
                {
                    var origDataFile = currentRawFileList[spectraFileIndex];
                    ICommonParameters combinedParams = SetAllFileSpecificCommonParams(CommonParameters, fileSettingsList[spectraFileIndex]);

                    var thisId = new List<string> { taskId, "Individual Spectra Files", origDataFile };
                    NewCollection(Path.GetFileName(origDataFile), thisId);
                    Status("Loading spectra file...", thisId);
                    IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile = myFileManager.LoadFile(origDataFile, combinedParams.TopNpeaks, combinedParams.MinRatio, combinedParams.TrimMs1Peaks, combinedParams.TrimMsMsPeaks);
                    Status("Getting ms2 scans...", thisId);
                    Ms2ScanWithSpecificMass[] arrayOfMs2ScansSortedByMass = GetMs2Scans(myMsDataFile, origDataFile, combinedParams.DoPrecursorDeconvolution, combinedParams.UseProvidedPrecursorInfo, combinedParams.DeconvolutionIntensityRatio, combinedParams.DeconvolutionMaxAssumedChargeState, combinedParams.DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();


                    //Import Database
                    Status("Loading modifications...", taskId);

                    #region Load modifications

                    List<ModificationWithMass> variableModifications = GlobalEngineLevelSettings.AllModsKnown.OfType<ModificationWithMass>().Where(b => CommonParameters.ListOfModsVariable.Contains(new Tuple<string, string>(b.modificationType, b.id))).ToList();
                    List<ModificationWithMass> fixedModifications = GlobalEngineLevelSettings.AllModsKnown.OfType<ModificationWithMass>().Where(b => CommonParameters.ListOfModsFixed.Contains(new Tuple<string, string>(b.modificationType, b.id))).ToList();
                    List<ModificationWithMass> localizeableModifications;
                    if (CommonParameters.LocalizeAll)
                        localizeableModifications = GlobalEngineLevelSettings.AllModsKnown.OfType<ModificationWithMass>().ToList();
                    else
                        localizeableModifications = GlobalEngineLevelSettings.AllModsKnown.OfType<ModificationWithMass>().Where(b => CommonParameters.ListOfModsLocalize.Contains(new Tuple<string, string>(b.modificationType, b.id))).ToList();

                    #endregion Load modifications 

                    var proteinList = dbFilenameList.SelectMany(b => LoadProteinDb(b.FilePath, true, DecoyType.None, localizeableModifications, b.IsContaminant, out Dictionary<string, Modification> unknownModifications)).ToList();


                    //Read N and C files
                    string nPath = NeoParameters.NFilePath;
                    string cPath = NeoParameters.CFilePath;
                    //if termini input

                    if (nPath == null || cPath == null)
                    {
                        //if no termini input
                        string taskHeader = "Task";
                        string[] pathArray = OutputFolder.Split('\\');
                        string basePath = "";
                        for (int i = 0; i < pathArray.Length - 1; i++)
                            basePath += pathArray[i] + '\\';
                        string currentTaskNumber = pathArray[pathArray.Length - 1].Split('-')[0];
                        currentTaskNumber = currentTaskNumber.Substring(taskHeader.Length, currentTaskNumber.Length - taskHeader.Length);
                        string NHeader = "";
                        string CHeader = "";
                        if (cPath == null)
                        {
                            CHeader = taskHeader + (Convert.ToInt16(currentTaskNumber) - 1);
                            if (nPath == null)
                                NHeader = taskHeader + (Convert.ToInt16(currentTaskNumber) - 2);
                        }
                        else
                            NHeader = taskHeader + (Convert.ToInt16(currentTaskNumber) - 1);
                        foreach (string s in Directory.GetDirectories(basePath))
                        {
                            if (s.Contains(NHeader))
                                nPath = s;
                            else if (s.Contains(CHeader))
                                cPath = s;
                        }
                    }

                    if (nPath == null || cPath == null)
                    {
                        //if no termini input
                        string taskHeader = "Task";
                        string[] pathArray = OutputFolder.Split('\\');
                        string basePath = "";
                        for (int i = 0; i < pathArray.Length - 1; i++)
                            basePath += pathArray[i] + '\\';
                        string currentTaskNumber = pathArray[pathArray.Length - 1].Split('-')[0];
                        currentTaskNumber = currentTaskNumber.Substring(taskHeader.Length, currentTaskNumber.Length - taskHeader.Length);
                        string NHeader = "";
                        string CHeader = "";
                        if (cPath == null)
                        {
                            CHeader = taskHeader + (Convert.ToInt16(currentTaskNumber) - 1);
                            if (nPath == null)
                                NHeader = taskHeader + (Convert.ToInt16(currentTaskNumber) - 2);
                        }
                        else
                            NHeader = taskHeader + (Convert.ToInt16(currentTaskNumber) - 1);
                        foreach (string s in Directory.GetDirectories(basePath))
                        {
                            if (s.Contains(NHeader))
                                nPath = s;
                            else if (s.Contains(CHeader))
                                cPath = s;
                        }
                        string fileName= Path.GetFileNameWithoutExtension(currentRawFileList[0]) + "_PSMs.psmtsv";
                        nPath += "\\" + fileName;
                        cPath += "\\" + fileName;
                    }
                    List<NeoPsm> psms = ImportPsmtsv.ImportNeoPsms(nPath, cPath);

                    //Splice
                    List<NeoPsm> candidates = NeoSplicePeptides.SplicePeptides(psms);

                    //Find Ambiguity
                    NeoFindAmbiguity.FindAmbiguity(candidates, proteinList, arrayOfMs2ScansSortedByMass);

                    //Export Results
                    NeoExport.ExportAll(candidates, arrayOfMs2ScansSortedByMass, OutputFolder);
                });
            }

            return myTaskResults;
        }

        public NeoSearchTask Clone()
        {
            return (NeoSearchTask)this.MemberwiseClone();
        }

        #endregion Public Properties
    }
}
