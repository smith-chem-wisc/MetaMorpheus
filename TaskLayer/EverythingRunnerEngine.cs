using EngineLayer;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;

namespace TaskLayer
{
    public class EverythingRunnerEngine
    {
        private readonly List<(string, MetaMorpheusTask)> TaskList;
        private string OutputFolder;
        private List<string> CurrentRawDataFilenameList;
        private List<DbForTask> CurrentXmlDbFilenameList;

        public EverythingRunnerEngine(List<(string, MetaMorpheusTask)> taskList, List<string> startingRawFilenameList, List<DbForTask> startingXmlDbFilenameList, string outputFolder)
        {
            TaskList = taskList;
            OutputFolder = outputFolder.Trim('"');

            CurrentRawDataFilenameList = startingRawFilenameList;
            CurrentXmlDbFilenameList = startingXmlDbFilenameList;
        }

        public static event EventHandler<StringEventArgs> FinishedWritingAllResultsFileHandler;

        public static event EventHandler StartingAllTasksEngineHandler;

        public static event EventHandler<StringEventArgs> FinishedAllTasksEngineHandler;

        public static event EventHandler<XmlForTaskListEventArgs> NewDbsHandler;

        public static event EventHandler<StringListEventArgs> NewSpectrasHandler;

        public static event EventHandler<StringListEventArgs> NewFileSpecificTomlHandler;

        public static event EventHandler<StringEventArgs> WarnHandler;

        public void Run()
        {
            StartingAllTasks();
            var stopWatch = new Stopwatch();
            stopWatch.Start();

            if (!CurrentRawDataFilenameList.Any())
            {
                Warn("No spectra files selected");
                FinishedAllTasks(null);
                return;
            }

            var startTimeForAllFilenames = DateTime.Now.ToString("yyyy-MM-dd-HH-mm-ss", CultureInfo.InvariantCulture);

            OutputFolder = OutputFolder.Replace("$DATETIME", startTimeForAllFilenames);

            StringBuilder allResultsText = new StringBuilder();

            for (int i = 0; i < TaskList.Count; i++)
            {
                if (!CurrentRawDataFilenameList.Any())
                {
                    Warn("Cannot proceed. No spectra files selected.");
                    FinishedAllTasks(OutputFolder);
                    return;
                }
                if (!CurrentXmlDbFilenameList.Any())
                {
                    Warn("Cannot proceed. No protein database files selected.");
                    FinishedAllTasks(OutputFolder);
                    return;
                }
                var ok = TaskList[i];

                // reset product types for custom fragmentation
                ok.Item2.CommonParameters.SetCustomProductTypes();

                var outputFolderForThisTask = Path.Combine(OutputFolder, ok.Item1);

                if (!Directory.Exists(outputFolderForThisTask))
                    Directory.CreateDirectory(outputFolderForThisTask);

                // Actual task running code
                var myTaskResults = ok.Item2.RunTask(outputFolderForThisTask, CurrentXmlDbFilenameList, CurrentRawDataFilenameList, ok.Item1);

                if (myTaskResults.NewDatabases != null)
                {
                    CurrentXmlDbFilenameList = myTaskResults.NewDatabases;
                    NewDBs(myTaskResults.NewDatabases);
                }
                if (myTaskResults.NewSpectra != null)
                {
                    if (CurrentRawDataFilenameList.Count == myTaskResults.NewSpectra.Count)
                    {
                        CurrentRawDataFilenameList = myTaskResults.NewSpectra;
                    }
                    else
                    {
                        // at least one file was not successfully calibrated
                        var successfullyCalibFiles = myTaskResults.NewSpectra.Select(p => Path.GetFileNameWithoutExtension(p).Replace(CalibrationTask.CalibSuffix, "")).ToList();
                        var origFiles = CurrentRawDataFilenameList.Select(p => Path.GetFileNameWithoutExtension(p)).ToList();
                        var unsuccessfullyCalibFiles = origFiles.Except(successfullyCalibFiles).ToList();
                        var unsuccessfullyCalibFilePaths = CurrentRawDataFilenameList.Where(p => unsuccessfullyCalibFiles.Contains(Path.GetFileNameWithoutExtension(p))).ToList();
                        CurrentRawDataFilenameList = myTaskResults.NewSpectra;
                        CurrentRawDataFilenameList.AddRange(unsuccessfullyCalibFilePaths);
                    }

                    NewSpectras(myTaskResults.NewSpectra);
                }
                if (myTaskResults.NewFileSpecificTomls != null)
                {
                    NewFileSpecificToml(myTaskResults.NewFileSpecificTomls);
                }
                allResultsText.AppendLine(Environment.NewLine + Environment.NewLine + Environment.NewLine + Environment.NewLine + myTaskResults.ToString());
            }
            stopWatch.Stop();
            var resultsFileName = Path.Combine(OutputFolder, "allResults.txt");
            using (StreamWriter file = new StreamWriter(resultsFileName))
            {
                file.WriteLine("MetaMorpheus: version " + GlobalVariables.MetaMorpheusVersion);
                file.WriteLine("Total time: " + stopWatch.Elapsed);
                file.Write(allResultsText.ToString());
            }
            FinishedWritingAllResultsFileHandler?.Invoke(this, new StringEventArgs(resultsFileName, null));
            FinishedAllTasks(OutputFolder);
        }

        private void Warn(string v)
        {
            WarnHandler?.Invoke(this, new StringEventArgs(v, null));
        }

        private void StartingAllTasks()
        {
            StartingAllTasksEngineHandler?.Invoke(this, EventArgs.Empty);
        }

        private void FinishedAllTasks(string rootOutputDir)
        {
            FinishedAllTasksEngineHandler?.Invoke(this, new StringEventArgs(rootOutputDir, null));
        }

        private void NewSpectras(List<string> newSpectra)
        {
            NewSpectrasHandler?.Invoke(this, new StringListEventArgs(newSpectra));
        }

        private void NewFileSpecificToml(List<string> newFileSpecificTomls)
        {
            NewFileSpecificTomlHandler?.Invoke(this, new StringListEventArgs(newFileSpecificTomls));
        }

        private void NewDBs(List<DbForTask> newDatabases)
        {
            NewDbsHandler?.Invoke(this, new XmlForTaskListEventArgs(newDatabases));
        }
    }
}