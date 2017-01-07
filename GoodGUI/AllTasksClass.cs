using MetaMorpheus;
using System.Collections.Generic;
using System.Collections.ObjectModel;

namespace GoodGUI
{
    internal static class AllTasksClass
    {
        internal static void DoAllTasks(IEnumerable<MyTask> taskList, ObservableCollection<RawData> rawDataAndResultslist, ObservableCollection<XMLdb> xMLdblist, AllTasksParams allTasksParams)
        {
            allTasksParams.startingAllTasks();
            foreach (var ok in taskList)
            {
                allTasksParams.startingSingleTask(new SingleTaskEventArgs(ok));
                ok.DoTask(rawDataAndResultslist, xMLdblist, allTasksParams);
                allTasksParams.finishedSingleTask(new SingleTaskEventArgs(ok));
            }
            allTasksParams.FinishedAllTasks();
        }
    }
}