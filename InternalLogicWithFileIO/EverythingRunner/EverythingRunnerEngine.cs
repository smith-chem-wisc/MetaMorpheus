using InternalLogicEngineLayer;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;

namespace InternalLogicTaskLayer
{
    public class EverythingRunnerEngine : MyEngine
    {
        private List<MyTaskEngine> taskList;
        private List<string> currentRawDataFilenameList;
        private List<string> currentXmlDbFilenameList;

        public static event EventHandler startingAllTasksEngineHandler;
        public static event EventHandler finishedAllTasksEngineHandler;

        public static event EventHandler<List<string>> newDbsHandler;

        public static event EventHandler<List<string>> newSpectrasHandler;

        public void StartingAllTasks()
        {
            startingAllTasksEngineHandler?.Invoke(this, EventArgs.Empty);
        }

        public void FinishedAllTasks()
        {
            finishedAllTasksEngineHandler?.Invoke(this, EventArgs.Empty);
        }

        private void NewSpectras(List<string> newSpectra)
        {
            newSpectrasHandler?.Invoke(this, newSpectra);
        }

        private void NewDBs(List<string> newDatabases)
        {
            newDbsHandler?.Invoke(this, newDatabases);

        }

        public EverythingRunnerEngine(List<MyTaskEngine> taskList, List<string> startingRawFilenameList, List<string> startingXmlDbFilenameList)
        {
            this.taskList = taskList;
            this.currentRawDataFilenameList = startingRawFilenameList;
            this.currentXmlDbFilenameList = startingXmlDbFilenameList;
        }

        public override void ValidateParams()
        {
            if (taskList == null)
                throw new EngineValidationException("taskList cannot be null");
            if (taskList.Count == 0)
                throw new EngineValidationException("taskList has to contain at least one element");
            if (currentRawDataFilenameList == null)
                throw new EngineValidationException("rawDataFilenameList cannot be null");
            if (currentRawDataFilenameList.Count == 0)
                throw new EngineValidationException("rawDataFilenameList has to contain at least one element");
            if (currentXmlDbFilenameList == null)
                throw new EngineValidationException("xmlDbFilenameList cannot be null");
            if (currentXmlDbFilenameList.Count == 0)
                throw new EngineValidationException("xmlDbFilenameList has to contain at least one element");
        }

        protected override MyResults RunSpecific()
        {
            StartingAllTasks();
            var startTimeForFilename = DateTime.Now.ToString("yyyy-MM-dd-HH-mm-ss", CultureInfo.InvariantCulture);

            var MatchingChars =
                from len in Enumerable.Range(0, currentRawDataFilenameList.Min(s => s.Length)).Reverse()
                let possibleMatch = currentRawDataFilenameList.First().Substring(0, len)
                where currentRawDataFilenameList.All(f => f.StartsWith(possibleMatch))
                select possibleMatch;

            var longestDir = Path.GetDirectoryName(MatchingChars.First());

            for (int i = 0; i < taskList.Count; i++)
            {
                var ok = taskList[i];
                string output_folder = null;
                if (taskList.Count == 1)
                {
                    output_folder = Path.Combine(longestDir, startTimeForFilename);
                }
                else
                {
                    output_folder = Path.Combine(longestDir, startTimeForFilename);
                    output_folder = Path.Combine(output_folder, "Task" + (i + 1) + ok.taskType);
                }

                if (!Directory.Exists(output_folder))
                    Directory.CreateDirectory(output_folder);
                ok.setOutputFolder(output_folder);
                ok.xmlDbFilenameList = currentXmlDbFilenameList;
                ok.rawDataFilenameList = currentRawDataFilenameList;

                MyTaskResults myTaskResults = (MyTaskResults)ok.Run();

                if (myTaskResults.newDatabases != null)
                {
                    currentRawDataFilenameList = myTaskResults.newDatabases;
                    NewDBs(myTaskResults.newDatabases);
                }
                if (myTaskResults.newSpectra != null)
                {
                    currentRawDataFilenameList = myTaskResults.newSpectra;
                    NewSpectras(myTaskResults.newSpectra);
                }
            }
            FinishedAllTasks();
            return null;
        }
    }
}