using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using EngineLayer;
using EngineLayer.DatabaseLoading;
using MassSpectrometry;
using SpectralAveraging;
using Nett;
using Readers;
using EngineLayer.Util;

namespace TaskLayer
{
    public class SpectralAveragingTask : MetaMorpheusTask
    {
        public SpectralAveragingParameters Parameters { get; set; }
        public const string AveragingSuffix = "-averaged";

        public SpectralAveragingTask(SpectralAveragingParameters parameters) : base(MyTask.Average)
        {
            Parameters = parameters;
        }

        /// <summary>
        /// Used when reading in toml files, and when writing the default toml that
        /// `CMD -g` emits. Initialises its parameter objects the way the other task
        /// types do in their parameterless constructors, so a default-constructed
        /// task serialises to a complete, runnable config rather than null sections.
        /// </summary>
        public SpectralAveragingTask() : base(MyTask.Average)
        {
            CommonParameters = new CommonParameters();
            Parameters = new SpectralAveragingParameters();
        }

        protected override MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId,
            FileSpecificParameters[] fileSettingsList)
        {
            // comment to force new checks on github

            //Start Averaging Task
            Status("Averaging...", new List<string>() { taskId });
            var myFileManager = new MyFileManager(true);
            List<string> unsuccessfulyAveragedFilePaths = new();
            MyTaskResults = new MyTaskResults(this)
            {
                NewSpectra = new List<string>(),
                NewFileSpecificTomls = new List<string>(),
            };

            for (int spectraFileIndex = 0; spectraFileIndex < currentRawFileList.Count; spectraFileIndex++)
            {
                // files not reached before a cancel keep -averaged rows in the new design; left so, as the user cancelled
                if (GlobalVariables.StopLoops) { break; }


                ReportProgress(new ProgressEventArgs((int)((spectraFileIndex / (double)currentRawFileList.Count) * 100), $"Averaging File {spectraFileIndex + 1}/{currentRawFileList.Count} ", new List<string> { taskId, "Individual Spectra Files" }));

                // get filename stuff
                var originalUnaveragedFilepath = currentRawFileList[spectraFileIndex];
                var originalUnaveragedFilepathWithoutExtenstion = Path.GetFileNameWithoutExtension(originalUnaveragedFilepath);
                var averagedFilepath = Path.Combine(OutputFolder, originalUnaveragedFilepathWithoutExtenstion + AveragingSuffix + ".mzML")
                    .ToSafeOutputPath( AveragingSuffix + ".mzML");
                var originalFileExtension = GlobalVariables.GetFileExtension(originalUnaveragedFilepath);
                if (originalFileExtension.Equals(".mgf", StringComparison.OrdinalIgnoreCase) || originalFileExtension.Equals(".d", StringComparison.OrdinalIgnoreCase) || originalFileExtension.Equals(".msalign", StringComparison.OrdinalIgnoreCase) || BrukerDataDirectory.IsInnerFileExtension(originalFileExtension))
                {
                    // record it, or the new experimental design names a -averaged.mzML that was never written
                    unsuccessfulyAveragedFilePaths.Add(originalUnaveragedFilepath);
                    Warn("Averaging for " + originalFileExtension + " files is not supported.");
                    FinishedDataFile(originalUnaveragedFilepath, new List<string> { taskId, "Individual Spectra Files", originalUnaveragedFilepath });
                    ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { taskId, "Individual Spectra Files", originalUnaveragedFilepathWithoutExtenstion }));
                    continue;
                }

                // mark file as in progress
                StartingDataFile(originalUnaveragedFilepath, new List<string>() {taskId, "Individual Spectra Files", originalUnaveragedFilepathWithoutExtenstion });

                // load the file
                Status("Loading spectra file...", new List<string> { taskId, "Individual Spectra Files", originalUnaveragedFilepathWithoutExtenstion });
                MsDataFile myMsdataFile = myFileManager.LoadFile(originalUnaveragedFilepath, CommonParameters);
                List<MsDataScan> scanList = myMsdataFile.GetAllScansList();

                // Average the spectra
                Status("Averaging spectra file...", new List<string> { taskId, "Individual Spectra Files", originalUnaveragedFilepathWithoutExtenstion });
                try
                {
                    var averagedScans = SpectraFileAveraging.AverageSpectraFile(scanList, Parameters);

                    Status("Writing spectra file...", new List<string> { taskId, "Individual Spectra Files", originalUnaveragedFilepathWithoutExtenstion });
                    SourceFile sourceFile = myMsdataFile.GetSourceFile();
                    MsDataFile dataFile = new GenericMsDataFile(averagedScans, sourceFile);
                    dataFile.ExportAsMzML(averagedFilepath, true);
                }
                catch (Exception e)
                {
                    // no -averaged.mzML was written: record the skip, and do not report or pass on the file
                    unsuccessfulyAveragedFilePaths.Add(originalUnaveragedFilepath);
                    Warn($"Averaging Failure! Could not average spectra for file {originalUnaveragedFilepathWithoutExtenstion} with exception {e.Message}");
                    myFileManager.DoneWithFile(originalUnaveragedFilepath);
                    FinishedDataFile(originalUnaveragedFilepath, new List<string> { taskId, "Individual Spectra Files", originalUnaveragedFilepath });
                    ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { taskId, "Individual Spectra Files", originalUnaveragedFilepathWithoutExtenstion }));
                    continue;
                }
                myFileManager.DoneWithFile(originalUnaveragedFilepath);

                // carry over file-specific parameters from the unaveraged file to the averaged one
                if (fileSettingsList[spectraFileIndex] != null)
                {
                    var fileSpecificParams = fileSettingsList[spectraFileIndex].Clone();// write toml settings for the averaged file if there are file specific parameters
                    var newTomlFileName = Path.Combine(OutputFolder, originalUnaveragedFilepathWithoutExtenstion + AveragingSuffix + ".toml");
                    Toml.WriteFile(fileSpecificParams, newTomlFileName, tomlConfig);
                    FinishedWritingFile(newTomlFileName, new List<string> { taskId, "Individual Spectra Files", originalUnaveragedFilepathWithoutExtenstion });
                    MyTaskResults.NewFileSpecificTomls.Add(newTomlFileName);
                }


                // finished averaging this file
                FinishedWritingFile(averagedFilepath, new List<string> { taskId, "Individual Spectra Files", originalUnaveragedFilepathWithoutExtenstion });
                MyTaskResults.NewSpectra.Add(averagedFilepath);
                FinishedDataFile(originalUnaveragedFilepath, new List<string> { taskId, "Individual Spectra Files", originalUnaveragedFilepath });
                ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { taskId, "Individual Spectra Files", originalUnaveragedFilepathWithoutExtenstion }));
            }

            // re-write experimental design (if it has been defined) with new calibrated file names
            string assumedPathToExperDesign = Directory.GetParent(currentRawFileList.First())?.FullName;
            if (assumedPathToExperDesign != null)
            {
                assumedPathToExperDesign =
                    Path.Combine(assumedPathToExperDesign, GlobalVariables.ExperimentalDesignFileName);
                if (File.Exists(assumedPathToExperDesign))
                {
                    WriteNewExperimentalDesignFile(assumedPathToExperDesign, OutputFolder, currentRawFileList,
                        unsuccessfulyAveragedFilePaths);
                }
            }

            // finished calibrating all files for the task
            ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { taskId, "Individual Spectra Files" }));

            return MyTaskResults;
        }

        private static void WriteNewExperimentalDesignFile(string pathToOldExperDesign, string outputFolder, List<string> originalUnaveragedFileNamesWithExtension,
            List<string> unsuccessfullyAveragedFilePaths)
        {
            var oldExperDesign = ExperimentalDesign.ReadExperimentalDesign(pathToOldExperDesign, originalUnaveragedFileNamesWithExtension, out var errors);

            if (errors.Any())
            {
                foreach (var error in errors)
                {
                    Warn(error);
                }

                return;
            }

            var newExperDesign = new List<SpectraFileInfo>();

            foreach (var unaveragedSpectraFile in oldExperDesign)
            {
                var originalUnaveragedFilePath = unaveragedSpectraFile.FullFilePathWithExtension;
                // as the run loop names it; GlobalVariables.GetFilenameWithoutExtension removes every copy of the extension
                var originalUnaveragedFilenameWithoutExtension = Path.GetFileNameWithoutExtension(originalUnaveragedFilePath);
                // ToSafeOutputPath here too, or a name the run loop trimmed is named untrimmed.
                string averagedFilePath = Path.Combine(outputFolder, originalUnaveragedFilenameWithoutExtension + AveragingSuffix + ".mzML")
                    .ToSafeOutputPath(AveragingSuffix + ".mzML");

                var averagedSpectraFile = new SpectraFileInfo(averagedFilePath,
                    unaveragedSpectraFile.Condition, unaveragedSpectraFile.BiologicalReplicate, unaveragedSpectraFile.TechnicalReplicate, unaveragedSpectraFile.Fraction);

                if (unsuccessfullyAveragedFilePaths.Contains(unaveragedSpectraFile.FullFilePathWithExtension))
                {
                    newExperDesign.Add(unaveragedSpectraFile);
                }
                else
                {
                    newExperDesign.Add(averagedSpectraFile);
                }
            }

            // into the output folder by name: a skipped file's row keeps its input path, so the first row's
            // folder can be the input folder, and the input design must not be overwritten
            ExperimentalDesign.WriteExperimentalDesignToFile(newExperDesign, outputFolder);
        }


    }
}
