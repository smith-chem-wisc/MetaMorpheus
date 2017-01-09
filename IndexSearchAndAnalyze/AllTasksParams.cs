using Proteomics;
using System;
using System.Collections.Generic;
using UsefulProteomicsDatabases.Generated;

namespace IndexSearchAndAnalyze
{
    public class AllTasksParams
    {
        public static unimod unimodDeserialized;
        public static Dictionary<int, ChemicalFormulaModification> uniprotDeseralized;

        public event EventHandler startingAllTasksHander;

        public event EventHandler finishedAllTasksHandler;

        public event EventHandler<SingleTaskEventArgs> startingSingleTaskHander;

        public event EventHandler<SingleTaskEventArgs> finishedSingleTaskHandler;

        public event EventHandler<string> outLabelStatusHandler;

        public event EventHandler<ProgressEventArgs> outProgressHandler;

        public event EventHandler<string> outRichTextBoxHandler;

        public event EventHandler<SingleFileEventArgs> SuccessfullyFinishedWritingFileHandler;

        public void status(string v)
        {
            outLabelStatusHandler?.Invoke(this, v);
        }

        public void startingAllTasks()
        {
            startingAllTasksHander?.Invoke(this, null);
        }

        public void startingSingleTask(SingleTaskEventArgs e)
        {
            startingSingleTaskHander?.Invoke(this, e);
        }

        public void finishedSingleTask(SingleTaskEventArgs e)
        {
            finishedSingleTaskHandler?.Invoke(this, e);
        }

        internal void ReportProgress(ProgressEventArgs v)
        {
            outProgressHandler?.Invoke(this, v);
        }

        public void RTBoutput(string v)
        {
            outRichTextBoxHandler?.Invoke(this, v);
        }

        internal void SucessfullyFinishedWritingFile(SingleFileEventArgs v)
        {
            SuccessfullyFinishedWritingFileHandler?.Invoke(this, v);
        }

        public void FinishedAllTasks()
        {
            finishedAllTasksHandler?.Invoke(this, null);
        }
    }
}