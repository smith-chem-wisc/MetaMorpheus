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
        #region Private Fields

        private readonly List<(string, MetaMorpheusTask)> taskList;
        private string outputFolder;
        private List<string> currentRawDataFilenameList;
        private List<DbForTask> currentXmlDbFilenameList;

        #endregion Private Fields

        #region Public Constructors

        public EverythingRunnerEngine(List<(string, MetaMorpheusTask)> taskList, List<string> startingRawFilenameList, List<DbForTask> startingXmlDbFilenameList, string outputFolder)
        {
            this.taskList = taskList;
            this.outputFolder = outputFolder;

            currentRawDataFilenameList = startingRawFilenameList;
            currentXmlDbFilenameList = startingXmlDbFilenameList;
        }

        #endregion Public Constructors

        #region Public Events

        public static event EventHandler<StringEventArgs> FinishedWritingAllResultsFileHandler;

        public static event EventHandler StartingAllTasksEngineHandler;

        public static event EventHandler<StringEventArgs> FinishedAllTasksEngineHandler;

        public static event EventHandler<XmlForTaskListEventArgs> NewDbsHandler;

        public static event EventHandler<StringListEventArgs> NewSpectrasHandler;

        public static event EventHandler<StringEventArgs> WarnHandler;

        #endregion Public Events

        #region Public Methods

        public void Run()
        {
            StartingAllTasks();
            var stopWatch = new Stopwatch();
            stopWatch.Start();

            if (!currentRawDataFilenameList.Any())
            {
                Warn("No data files selected");
                FinishedAllTasks(null);
                return;
            }

            var startTimeForAllFilenames = DateTime.Now.ToString("yyyy-MM-dd-HH-mm-ss", CultureInfo.InvariantCulture);

            outputFolder = outputFolder.Replace("$DATETIME", startTimeForAllFilenames);

            StringBuilder allResultsText = new StringBuilder();

            for (int i = 0; i < taskList.Count; i++)
            {
                if (!currentRawDataFilenameList.Any())
                {
                    Warn("Cannot proceed. No data files selected.");
                    FinishedAllTasks(outputFolder);
                    return;
                }
                if (!currentXmlDbFilenameList.Any())
                {
                    Warn("Cannot proceed. No xml files selected.");
                    FinishedAllTasks(outputFolder);
                    return;
                }
                var ok = taskList[i];

                var outputFolderForThisTask = Path.Combine(outputFolder, ok.Item1);

                if (!Directory.Exists(outputFolderForThisTask))
                    Directory.CreateDirectory(outputFolderForThisTask);

                // Actual task running code
                var myTaskResults = ok.Item2.RunTask(outputFolderForThisTask, currentXmlDbFilenameList, currentRawDataFilenameList, ok.Item1);

                if (myTaskResults.newDatabases != null)
                {
                    currentXmlDbFilenameList = myTaskResults.newDatabases;
                    NewDBs(myTaskResults.newDatabases);
                }
                if (myTaskResults.newSpectra != null)
                {
                    currentRawDataFilenameList = myTaskResults.newSpectra;
                    NewSpectras(myTaskResults.newSpectra);
                }
                allResultsText.AppendLine(Environment.NewLine + Environment.NewLine + Environment.NewLine + Environment.NewLine + myTaskResults.ToString());
            }
            stopWatch.Stop();
            var resultsFileName = Path.Combine(outputFolder, "allResults.txt");
            using (StreamWriter file = new StreamWriter(resultsFileName))
            {
                file.WriteLine("MetaMorpheus: version " + GlobalEngineLevelSettings.MetaMorpheusVersion);
                file.WriteLine("Total time: " + stopWatch.Elapsed);
                file.Write(allResultsText.ToString());
            }
            FinishedWritingAllResultsFileHandler?.Invoke(this, new StringEventArgs(resultsFileName, null));
            FinishedAllTasks(outputFolder);
        }

        #endregion Public Methods

        #region Private Methods

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

        private void NewDBs(List<DbForTask> newDatabases)
        {
            NewDbsHandler?.Invoke(this, new XmlForTaskListEventArgs(newDatabases));
        }

        #endregion Private Methods
    }
}