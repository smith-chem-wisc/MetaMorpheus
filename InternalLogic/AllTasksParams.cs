using Proteomics;
using System;
using System.Collections.Generic;
using UsefulProteomicsDatabases.Generated;

namespace InternalLogic
{
    public class AllTasksParams
    {

        public List<string> rawDataAndResultslist;

        public List<string> xMLdblist;

        //public void status(string v)
        //{
        //    outLabelStatusHandler?.Invoke(this, v);
        //}

        //public void startingAllTasks()
        //{
        //    startingAllTasksHander?.Invoke(this, null);
        //}

        //public void startingSingleTask(SingleTaskEventArgs e)
        //{
        //    startingSingleTaskHander?.Invoke(this, e);
        //}

        //public void finishedSingleTask(SingleTaskEventArgs e)
        //{
        //    finishedSingleTaskHandler?.Invoke(this, e);
        //}

        //internal void ReportProgress(ProgressEventArgs v)
        //{
        //    outProgressHandler?.Invoke(this, v);
        //}

        //public void output(string v)
        //{
        //    outRichTextBoxHandler?.Invoke(this, v);
        //}

        //internal void SucessfullyFinishedWritingFile(SingleFileEventArgs v)
        //{
        //    SuccessfullyFinishedWritingFileHandler?.Invoke(this, v);
        //}

        //public void NewDBs(List<string> newDatabases)
        //{
        //    xMLdblist = newDatabases;
        //    newDbsHandler?.Invoke(this, newDatabases);
        //}

        //public void FinishedAllTasks()
        //{
        //    finishedAllTasksHandler?.Invoke(this, null);
        //}

        //public void NewSpectras(List<string> newSpectra)
        //{
        //    rawDataAndResultslist = newSpectra;
        //    newSpectrasHandler?.Invoke(this, newSpectra);
        //}
    }
}