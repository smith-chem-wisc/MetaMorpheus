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

        private readonly List<Tuple<string, MetaMorpheusTask>> taskList;
        private List<string> currentRawDataFilenameList;
        private List<DbForTask> currentXmlDbFilenameList;

        #endregion Private Fields

        #region Public Constructors

        public EverythingRunnerEngine(List<Tuple<string, MetaMorpheusTask>> taskList, List<string> startingRawFilenameList, List<DbForTask> startingXmlDbFilenameList)
        {
            this.taskList = taskList;
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

            var MatchingChars =
                from len in Enumerable.Range(0, currentRawDataFilenameList.Min(s => s.Length)).Reverse()
                let possibleMatch = currentRawDataFilenameList.First().Substring(0, len)
                where currentRawDataFilenameList.All(f => f.StartsWith(possibleMatch, StringComparison.Ordinal))
                select possibleMatch;

            var longestDir = Path.GetDirectoryName(MatchingChars.First());

            string rootOutputDir = Path.Combine(longestDir, startTimeForAllFilenames);

            StringBuilder allResultsText = new StringBuilder();

            for (int i = 0; i < taskList.Count; i++)
            {
                if (!currentRawDataFilenameList.Any())
                {
                    Warn("Cannot proceed. No data files selected.");
                    FinishedAllTasks(rootOutputDir);
                    return;
                }
                if (!currentXmlDbFilenameList.Any())
                {
                    Warn("Cannot proceed. No xml files selected.");
                    FinishedAllTasks(rootOutputDir);
                    return;
                }
                var ok = taskList[i];

                var outputFolderForThisTask = Path.Combine(rootOutputDir, ok.Item1);

                if (!Directory.Exists(outputFolderForThisTask))
                    Directory.CreateDirectory(outputFolderForThisTask);

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

            var resultsFileName = Path.Combine(rootOutputDir, "allResults.txt");
            using (StreamWriter file = new StreamWriter(resultsFileName))
            {
                file.WriteLine(GlobalEngineLevelSettings.MetaMorpheusVersion.Equals("1.0.0.0") ? "MetaMorpheus: Not a release version" : "MetaMorpheus: version " + GlobalEngineLevelSettings.MetaMorpheusVersion);
                file.WriteLine("Total time: " + stopWatch.Elapsed);
                file.Write(allResultsText.ToString());
            }
            FinishedWritingAllResultsFileHandler?.Invoke(this, new StringEventArgs(resultsFileName, null));

            FinishedAllTasks(rootOutputDir);
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