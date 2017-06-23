using EngineLayer;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;

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

        public static event EventHandler startingAllTasksEngineHandler;

        public static event EventHandler<string> finishedAllTasksEngineHandler;

        public static event EventHandler<XmlForTaskListEventArgs> newDbsHandler;

        public static event EventHandler<StringListEventArgs> newSpectrasHandler;

        public static event EventHandler<StringEventArgs> warnHandler;

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

            for (int i = 0; i < taskList.Count; i++)
            {
                if (!currentRawDataFilenameList.Any())
                {
                    Warn("Cannot proceed. No data files selected.");
                    FinishedAllTasks(rootOutputDir);
                }
                if (!currentXmlDbFilenameList.Any())
                {
                    Warn("Cannot proceed. No xml files selected.");
                    FinishedAllTasks(rootOutputDir);
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
            }
            stopWatch.Stop();
            FinishedAllTasks(rootOutputDir);
        }

        #endregion Public Methods

        #region Private Methods

        private void Warn(string v)
        {
            warnHandler?.Invoke(this, new StringEventArgs(v, null));
        }

        private void StartingAllTasks()
        {
            startingAllTasksEngineHandler?.Invoke(this, EventArgs.Empty);
        }

        private void FinishedAllTasks(string rootOutputDir)
        {
            finishedAllTasksEngineHandler?.Invoke(this, rootOutputDir);
        }

        private void NewSpectras(List<string> newSpectra)
        {
            newSpectrasHandler?.Invoke(this, new StringListEventArgs(newSpectra));
        }

        private void NewDBs(List<DbForTask> newDatabases)
        {
            newDbsHandler?.Invoke(this, new XmlForTaskListEventArgs(newDatabases));
        }

        #endregion Private Methods

    }
}
