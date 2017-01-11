using InternalLogic;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;

namespace InternalLogicWithFileIO
{
    public static class AllTasksClass
    {
        public static void DoAllTasks(List<MyTask> taskList, AllTasksParams allTasksParams)
        {
            allTasksParams.startingAllTasks();
            var startTimeForFilename = DateTime.Now.ToString("yyyy-MM-dd-HH-mm-ss", CultureInfo.InvariantCulture);

            var MatchingChars =
                from len in Enumerable.Range(0, allTasksParams.rawDataAndResultslist.Min(s => s.Length)).Reverse()
                let possibleMatch = allTasksParams.rawDataAndResultslist.First().Substring(0, len)
                where allTasksParams.rawDataAndResultslist.All(f => f.StartsWith(possibleMatch))
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
                MyTaskResults myTaskResults = ok.DoTask(allTasksParams);
                if (myTaskResults.newDatabases != null)
                    allTasksParams.NewDBs(myTaskResults.newDatabases);
                if (myTaskResults.newSpectra != null)
                    allTasksParams.NewSpectras(myTaskResults.newSpectra);
                //allTasksParams.finishedSingleTask(new SingleTaskEventArgs(ok));
            }
            allTasksParams.FinishedAllTasks();
        }
    }
}