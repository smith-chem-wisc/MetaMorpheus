using System;

namespace GoodGUI
{
    internal class AllTasksParams
    {
        public event EventHandler startingAllTasksHander;

        public event EventHandler finishedAllTasksHandler;

        public event EventHandler<SingleTaskEventArgs> startingSingleTaskHander;

        public event EventHandler<SingleTaskEventArgs> finishedSingleTaskHandler;

        public event EventHandler<string> outLabelStatusHandler;

        public event EventHandler<string> outTextBoxHandler;

        public event EventHandler<int> outProgressBarHandler;

        public event EventHandler<string> outRichTextBoxHandler;

        public event EventHandler<string> SuccessfullyFinishedFileHandler;

        internal void statusLabel(string v)
        {
            outLabelStatusHandler?.Invoke(this, v);
        }

        internal void textBoxOutput(string v)
        {
            outTextBoxHandler?.Invoke(this, v);
        }

        internal void startingAllTasks()
        {
            startingAllTasksHander?.Invoke(this, null);
        }

        internal void startingSingleTask(SingleTaskEventArgs e)
        {
            startingSingleTaskHander?.Invoke(this, e);
        }

        internal void finishedSingleTask(SingleTaskEventArgs e)
        {
            finishedSingleTaskHandler?.Invoke(this, e);
        }

        internal void ReportProgress(int v)
        {
            outProgressBarHandler?.Invoke(this, v);
        }

        internal void RTBoutput(string v)
        {
            outRichTextBoxHandler?.Invoke(this, v);
        }

        internal void FinishedFile(string v)
        {
            SuccessfullyFinishedFileHandler?.Invoke(this, v);
        }

        internal void FinishedAllTasks()
        {
            finishedAllTasksHandler?.Invoke(this, null);
        }
    }
}