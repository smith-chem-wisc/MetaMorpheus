using EngineLayer;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using EngineLayer.DatabaseLoading;

namespace TaskLayer
{
    public class EverythingRunnerEngine
    {
        private readonly List<(string, MetaMorpheusTask)> TaskList;
        private string OutputFolder;
        private List<string> CurrentRawDataFilenameList;
        private List<DbForTask> CurrentXmlDbFilenameList;
        private List<string> _warnings;

        public EverythingRunnerEngine(List<(string, MetaMorpheusTask)> taskList, List<string> startingRawFilenameList, List<DbForTask> startingXmlDbFilenameList, string outputFolder)
        {
            TaskList = taskList;
            OutputFolder = outputFolder.Trim('"');

            CurrentRawDataFilenameList = startingRawFilenameList;
            CurrentXmlDbFilenameList = startingXmlDbFilenameList;
            _warnings = new();
        }

        public static event EventHandler<StringEventArgs> FinishedWritingAllResultsFileHandler;

        public static event EventHandler StartingAllTasksEngineHandler;

        public static event EventHandler<StringEventArgs> FinishedAllTasksEngineHandler;

        public static event EventHandler<XmlForTaskListEventArgs> NewDbsHandler;

        public static event EventHandler<StringListEventArgs> NewSpectrasHandler;

        public static event EventHandler<StringListEventArgs> NewFileSpecificTomlHandler;

        public static event EventHandler<StringEventArgs> WarnHandler;

        public List<string> Warnings { get { return _warnings; } private set { } }

        /// <summary>
        /// True when Run stopped at one of the "Cannot proceed" gates instead of running the task list.
        /// WarnHandler is the only other signal, and the CMD does not subscribe to it, so without this a
        /// refusal is indistinguishable from a completed run to anything but the GUI.
        /// </summary>
        public bool RefusedToProceed { get; private set; }

        public void Run()
        {
            StartingAllTasks();
            var stopWatch = new Stopwatch();
            stopWatch.Start();

            if (!CurrentRawDataFilenameList.Any())
            {
                RefuseToProceed("No spectra files selected", null);
                return;
            }

            var startTimeForAllFilenames = DateTime.Now.ToString("yyyy-MM-dd-HH-mm-ss", CultureInfo.InvariantCulture);

            OutputFolder = OutputFolder.Replace("$DATETIME", startTimeForAllFilenames);

            StringBuilder allResultsText = new StringBuilder();

            // Refused before any task runs, not when the offending one comes up: whether a library is
            // present cannot change during a run. The only writers of NewDatabases are GPTMD passing an
            // existing library along and the update path itself, both of which need one already in the
            // list -- SpectralLibraryGeneration writes the .msp but never registers it, so writing a
            // library and updating it in a later task does not work. A per-task check would refuse the
            // same runs, having first run every task ahead of the search in full.
            if (TaskList.Any(t => t.Item2 is SearchTask searchTask && searchTask.SearchParameters.UpdateSpectralLibrary)
                && !CurrentXmlDbFilenameList.AnySpectralLibrary())
            {
                RefuseToProceed("Cannot proceed. Updating a spectral library was requested, but no spectral "
                    + "library was given. Add one to the list of databases, or select writing a new spectral "
                    + "library instead of updating one.", OutputFolder);
                return;
            }

            for (int i = 0; i < TaskList.Count; i++)
            {
                if (!CurrentRawDataFilenameList.Any())
                {
                    RefuseToProceed("Cannot proceed. No spectra files selected.", OutputFolder);
                    return;
                }
                if (!CurrentXmlDbFilenameList.Any() && !(TaskList.Count == 1 && TaskList.First().Item2 is SpectralAveragingTask))
                {
                    RefuseToProceed("Cannot proceed. No protein database files selected.", OutputFolder);
                    return;
                }
                else if (CurrentXmlDbFilenameList.SpectralLibraries().Count() == CurrentXmlDbFilenameList.Count
                             && !(TaskList.Count == 1 && TaskList.First().Item2 is SpectralAveragingTask))
                {
                    RefuseToProceed("Cannot proceed. No protein database files selected.", OutputFolder);
                    return;
                }

                var ok = TaskList[i];

                // Non-specific search is built around proteases -- terminal mod placement, the "single"
                // agents, the FDR categories -- none of which have a nucleic acid counterpart yet.
                // Refused here rather than only in SearchTask because a MetaMorpheusException out of
                // RunSpecific is dumped into results.txt with a stack trace and rethrown, and the GUI
                // routes the faulted task to EverythingRunnerExceptionHandler -- so the user is told
                // MetaMorpheus crashed and invited to file a bug, and never sees the message that says
                // what to do instead. The throw in SearchTask stays as a backstop for a caller invoking
                // RunTask directly.
                if (ok.Item2 is SearchTask nonSpecificCandidate
                    && nonSpecificCandidate.SearchParameters.SearchType == SearchType.NonSpecific
                    && GlobalVariables.AnalyteType == AnalyteType.Oligo)
                {
                    RefuseToProceed("Cannot proceed. Non-specific search is only implemented for proteins. " +
                        "Use Classic or Modern search for nucleic acid databases.", OutputFolder);
                    return;
                }

                // reset product types for custom fragmentation
                ok.Item2.CommonParameters.SetCustomProductTypes();

                var outputFolderForThisTask = Path.Combine(OutputFolder, ok.Item1);

                if (!Directory.Exists(outputFolderForThisTask))
                    Directory.CreateDirectory(outputFolderForThisTask);

                // Actual task running code
                var myTaskResults = ok.Item2.RunTask(outputFolderForThisTask, CurrentXmlDbFilenameList, CurrentRawDataFilenameList, ok.Item1);

                if (myTaskResults.NewDatabases != null)
                {
                    CurrentXmlDbFilenameList = myTaskResults.NewDatabases;
                    NewDBs(myTaskResults.NewDatabases);
                }
                if (myTaskResults.NewSpectra != null)
                {
                    if (CurrentRawDataFilenameList.Count == myTaskResults.NewSpectra.Count)
                    {
                        CurrentRawDataFilenameList = myTaskResults.NewSpectra;
                    }
                    else
                    {
                        // at least one file was not successfully calibrated
                        var successfulFiles = myTaskResults.NewSpectra.Select(p => Path.GetFileNameWithoutExtension(p)
                            .Replace(CalibrationTask.CalibSuffix, "")
                            .Replace(SpectralAveragingTask.AveragingSuffix, "")).ToList();
                        var origFiles = CurrentRawDataFilenameList.Select(Path.GetFileNameWithoutExtension).ToList();
                        var unsuccessfulFiles = origFiles.Except(successfulFiles).ToList();
                        var unsuccessfulFilePaths = CurrentRawDataFilenameList.Where(p => unsuccessfulFiles.Contains(Path.GetFileNameWithoutExtension(p))).ToList();
                        CurrentRawDataFilenameList = myTaskResults.NewSpectra;
                        CurrentRawDataFilenameList.AddRange(unsuccessfulFilePaths);
                    }

                    NewSpectras(myTaskResults.NewSpectra);
                }
                if (myTaskResults.NewFileSpecificTomls != null)
                {
                    NewFileSpecificToml(myTaskResults.NewFileSpecificTomls);
                }
                allResultsText.AppendLine(Environment.NewLine + Environment.NewLine + Environment.NewLine + Environment.NewLine + myTaskResults.ToString());
            }
            stopWatch.Stop();
            var resultsFileName = Path.Combine(OutputFolder, "allResults.txt");
            using (StreamWriter file = new StreamWriter(resultsFileName))
            {
                file.WriteLine("MetaMorpheus: version " + GlobalVariables.MetaMorpheusVersion);
                file.WriteLine("Total time: " + stopWatch.Elapsed);
                file.Write(allResultsText.ToString());
            }
            FinishedWritingAllResultsFileHandler?.Invoke(this, new StringEventArgs(resultsFileName, null));
            FinishedAllTasks(OutputFolder);
        }

        /// <summary>
        /// Warn and stop cleanly, rather than throwing: a MetaMorpheusException out of here reaches the GUI
        /// as a crash report. <see cref="RefusedToProceed"/> is what lets a caller tell this apart from a
        /// run that finished.
        /// </summary>
        private void RefuseToProceed(string reason, string outputFolder)
        {
            RefusedToProceed = true;
            Warn(reason);
            FinishedAllTasks(outputFolder);
        }

        private void Warn(string v)
        {
            WarnHandler?.Invoke(this, new StringEventArgs(v, null));
            _warnings.Add(v);
        }

        private void StartingAllTasks()
        {
            StartingAllTasksEngineHandler?.Invoke(this, EventArgs.Empty);
        }

        private void FinishedAllTasks(string rootOutputDir)
        {
            FinishedAllTasksEngineHandler?.Invoke(this, new StringEventArgs(rootOutputDir, null));
        }

        private void NewSpectras(List<string> newSpectra)
        {
            NewSpectrasHandler?.Invoke(this, new StringListEventArgs(newSpectra));
        }

        private void NewFileSpecificToml(List<string> newFileSpecificTomls)
        {
            NewFileSpecificTomlHandler?.Invoke(this, new StringListEventArgs(newFileSpecificTomls));
        }

        private void NewDBs(List<DbForTask> newDatabases)
        {
            NewDbsHandler?.Invoke(this, new XmlForTaskListEventArgs(newDatabases));
        }
    }
}