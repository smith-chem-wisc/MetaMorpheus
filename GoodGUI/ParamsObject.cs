using System;

namespace GoodGUI
{
    internal class ParamsObject
    {
        public void AttachoutLabelStatusHandler(Action<EventHandler<string>> attach)
        {
            attach(outLabelStatusHandler);
        }

        public void AttachoSuccessfullyFinishedFileHandler(Action<EventHandler<string>> attach)
        {
            attach(SuccessfullyFinishedFileHandler);
        }

        public void AttachoutRichTextBoxHandlerr(Action<EventHandler<string>> attach)
        {
            attach(outRichTextBoxHandler);
        }

        public void AttachoutProgressBarHandler(Action<EventHandler<int>> attach)
        {
            attach(outProgressBarHandler);
        }

        public event EventHandler<string> outLabelStatusHandler;

        public event EventHandler<string> outTextBoxHandler;

        public event EventHandler<int> outProgressBarHandler;

        public event EventHandler<string> outRichTextBoxHandler;

        public event EventHandler outSuccessfullyStartingTaskHandler;

        public event EventHandler outSuccessfullyFinishedTaskHandler;

        public event EventHandler<string> SuccessfullyFinishedFileHandler;

        internal void statusLabel(string v)
        {
            outLabelStatusHandler?.Invoke(this, v);
        }

        internal void textBoxOutput(string v)
        {
            outTextBoxHandler?.Invoke(this, v);
        }

        internal void startingTask()
        {
            outSuccessfullyStartingTaskHandler?.Invoke(this, null);
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

        internal void FinishedTask()
        {
            outSuccessfullyFinishedTaskHandler?.Invoke(this, null);
        }
    }
}