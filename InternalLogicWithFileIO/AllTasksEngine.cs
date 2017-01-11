using InternalLogic;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;

namespace InternalLogicWithFileIO
{
    public class AllTasksEngine : MyTaskEngine
    {
        List<MyTaskEngine> taskList;
        public AllTasksEngine(List<MyTaskEngine> taskList)
        {
            this.taskList = taskList;
        }

        public override void ValidateParams()
        {
            throw new NotImplementedException();
        }

        protected override MyResults RunSpecific()
        {
            StartingAllTasks();
            var startTimeForFilename = DateTime.Now.ToString("yyyy-MM-dd-HH-mm-ss", CultureInfo.InvariantCulture);

            var MatchingChars =
                from len in Enumerable.Range(0, rawDataAndResultslist.Min(s => s.Length)).Reverse()
                let possibleMatch = rawDataAndResultslist.First().Substring(0, len)
                where rawDataAndResultslist.All(f => f.StartsWith(possibleMatch))
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

                //allTasksParams.startingSingleTask(new SingleTaskEventArgs(ok));
                MyTaskResults myTaskResults = (MyTaskResults)ok.Run();
                if (myTaskResults.newDatabases != null)
                    NewDBs(myTaskResults.newDatabases);
                if (myTaskResults.newSpectra != null)
                    NewSpectras(myTaskResults.newSpectra);
                //allTasksParams.finishedSingleTask(new SingleTaskEventArgs(ok));
            }
            FinishedAllTasks();
            return new AllTasksResults(this);
        }

        private void NewSpectras(List<string> newSpectra)
        {
            throw new NotImplementedException();
        }

        private void NewDBs(List<string> newDatabases)
        {
            throw new NotImplementedException();
        }

        private void FinishedAllTasks()
        {
            throw new NotImplementedException();
        }

        private void StartingAllTasks()
        {
            throw new NotImplementedException();
        }
    }
}