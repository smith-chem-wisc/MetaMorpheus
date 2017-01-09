using IndexSearchAndAnalyze;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Linq;

namespace GoodGUI
{
    internal class TaskListWrapper
    {
        private ObservableCollection<MyTask> taskList;
        public int Count { get { return taskList.Count; } }

        public TaskListWrapper(ObservableCollection<MyTask> taskList)
        {
            this.taskList = taskList;
        }

        internal List<MyTask> EnumerateTasks()
        {
            return taskList.ToList();
        }

        internal void Clear()
        {
            taskList.Clear();
        }

        internal void RemoveLast()
        {
            taskList.RemoveAt(taskList.Count - 1);
        }

        internal void Add(MyTask theTask)
        {
            taskList.Insert(taskList.Count, theTask);
        }
    }
}