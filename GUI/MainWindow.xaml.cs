using EngineLayer;
using IO.ThermoRawFileReader;
using Microsoft.Win32;
using MzLibUtil;
using Nett;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Navigation;
using TaskLayer;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        private readonly ObservableCollection<RawDataForDataGrid> SpectraFiles = new ObservableCollection<RawDataForDataGrid>();
        private ObservableCollection<ProteinDbForDataGrid> ProteinDatabases = new ObservableCollection<ProteinDbForDataGrid>();
        private readonly ObservableCollection<PreRunTask> PreRunTasks = new ObservableCollection<PreRunTask>();
        private readonly ObservableCollection<RawDataForDataGrid> SelectedSpectraFiles = new ObservableCollection<RawDataForDataGrid>();
        private readonly ObservableCollection<ProteinDbForDataGrid> SelectedProteinDatabaseFiles = new ObservableCollection<ProteinDbForDataGrid>();
        private ObservableCollection<InRunTask> InProgressTasks;
        public static string NewestKnownMetaMorpheusVersion { get; private set; }

        public MainWindow()
        {
            InitializeComponent();
            GlobalVariables.SetUpGlobalVariables();

            Title = "MetaMorpheus: version " + GlobalVariables.MetaMorpheusVersion;

            dataGridProteinDatabases.DataContext = ProteinDatabases;
            proteinDbSummaryDataGrid.DataContext = ProteinDatabases;

            dataGridSpectraFiles.DataContext = SpectraFiles;
            spectraFileSummaryDataGrid.DataContext = SpectraFiles;

            tasksTreeView.DataContext = PreRunTasks;
            taskSummary.DataContext = PreRunTasks;

            EverythingRunnerEngine.NewDbsHandler += AddNewProteinDatabaseFromGptmd;
            EverythingRunnerEngine.NewSpectrasHandler += AddNewSpectraFileFromCalibration;
            EverythingRunnerEngine.NewFileSpecificTomlHandler += AddNewFileSpecificTomlFromCalibration;
            EverythingRunnerEngine.StartingAllTasksEngineHandler += SuccessfullyStartingAllTasks;
            EverythingRunnerEngine.FinishedAllTasksEngineHandler += SuccessfullyFinishedAllTasks;
            EverythingRunnerEngine.WarnHandler += NotificationHandler;
            EverythingRunnerEngine.FinishedWritingAllResultsFileHandler += FinishedWritingAllResultsFileHandler;

            MetaMorpheusTask.StartingSingleTaskHander += StartingTaskHander;
            MetaMorpheusTask.FinishedSingleTaskHandler += FinishedTaskHandler;
            MetaMorpheusTask.FinishedWritingFileHandler += FinishedWritingFile;
            MetaMorpheusTask.StartingDataFileHandler += StartingSpectraFileHandler;
            MetaMorpheusTask.FinishedDataFileHandler += FinishedSpectraFileHandler;
            MetaMorpheusTask.OutLabelStatusHandler += NewoutLabelStatus;
            MetaMorpheusTask.NewCollectionHandler += AddBranchToTreeViewHandler;
            MetaMorpheusTask.OutProgressHandler += NewoutProgressBar;
            MetaMorpheusTask.WarnHandler += NotificationHandler;

            MetaMorpheusEngine.OutProgressHandler += NewoutProgressBar;
            MetaMorpheusEngine.OutLabelStatusHandler += NewoutLabelStatus;
            MetaMorpheusEngine.WarnHandler += NotificationHandler;

            MyFileManager.WarnHandler += NotificationHandler;
            Application.Current.MainWindow.Closing += new CancelEventHandler(MainWindow_Closing);
        }

        private void MyWindow_Loaded(object sender, RoutedEventArgs e)
        {
            UpdateGuiOnPreRunChange();
            UpdateOutputFolderTextbox();
            FileSpecificParameters.ValidateFileSpecificVariableNames();
            SearchModifications.SetUpModSearchBoxes();
            PrintErrorsReadingMods();

            if (!UpdateGUISettings.LoadGUISettings())
            {
                notificationsTextBox.Document = GetWelcomeText();
            }

            if (UpdateGUISettings.Params.AskAboutUpdating)
            {
                UpdateMetaMorpheus();
            }
        }

        #region Events triggered by MetaMorpheus

        private void FinishedWritingAllResultsFileHandler(object sender, StringEventArgs e)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => FinishedWritingAllResultsFileHandler(sender, e)));
            }
            else
            {
                InProgressTasks.Add(new InRunTask("All Task Results", null));
                InProgressTasks.Last().Progress = 100;
                InProgressTasks.Last().Children.Add(new OutputFileForTreeView(e.S, Path.GetFileNameWithoutExtension(e.S)));
            }
        }

        private void NotificationHandler(object sender, StringEventArgs e)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => NotificationHandler(sender, e)));
            }
            else
            {
                notificationsTextBox.AppendText(e.S);
                notificationsTextBox.AppendText(Environment.NewLine);

                NotificationExpander.IsExpanded = true;
            }
        }

        private void FinishedSpectraFileHandler(object sender, StringEventArgs s)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => FinishedSpectraFileHandler(sender, s)));
            }
            else
            {
                RawDataForDataGrid spectraFile = SpectraFiles.First(b => b.FilePath.Equals(s.S));
                spectraFile.SetInProgress(false);

                dataGridSpectraFiles.Items.Refresh();
                spectraFileSummaryDataGrid.Items.Refresh();
            }
        }

        private void StartingSpectraFileHandler(object sender, StringEventArgs s)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => StartingSpectraFileHandler(sender, s)));
            }
            else
            {
                RawDataForDataGrid spectraFile = SpectraFiles.First(b => b.FilePath.Equals(s.S));
                spectraFile.SetInProgress(true);
                dataGridSpectraFiles.Items.Refresh();
                spectraFileSummaryDataGrid.Items.Refresh();
            }
        }

        private void AddNewProteinDatabaseFromGptmd(object sender, XmlForTaskListEventArgs e)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => AddNewProteinDatabaseFromGptmd(sender, e)));
            }
            else
            {
                foreach (var db in ProteinDatabases)
                {
                    db.Use = false;
                }

                foreach (var db in e.NewDatabases)
                {
                    ProteinDatabases.Add(new ProteinDbForDataGrid(db));
                }

                dataGridProteinDatabases.Items.Refresh();
                proteinDbSummaryDataGrid.Items.Refresh();
            }
        }

        private void AddNewSpectraFileFromCalibration(object sender, StringListEventArgs e)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => AddNewSpectraFileFromCalibration(sender, e)));
            }
            else
            {
                var newFiles = e.StringList.ToList();
                foreach (var oldFile in SpectraFiles)
                {
                    if (!newFiles.Contains(oldFile.FilePath))
                    {
                        oldFile.Use = false;
                    }
                }

                var files = SpectraFiles.Select(p => p.FilePath).ToList();
                foreach (var newRawData in newFiles.Where(p => !files.Contains(p)))
                {
                    SpectraFiles.Add(new RawDataForDataGrid(newRawData));
                }

                UpdateOutputFolderTextbox();
            }
        }

        private void AddNewFileSpecificTomlFromCalibration(object sender, StringListEventArgs e)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => AddNewFileSpecificTomlFromCalibration(sender, e)));
            }
            else
            {
                foreach (var path in e.StringList)
                {
                    UpdateFileSpecificParamsDisplay(path);
                }
            }
        }

        private void StartingTaskHander(object sender, SingleTaskEventArgs s)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => StartingTaskHander(sender, s)));
            }
            else
            {
                var theTask = InProgressTasks.First(b => b.DisplayName.Equals(s.DisplayName));
                theTask.Status = "Starting...";

                dataGridSpectraFiles.Items.Refresh();
                spectraFileSummaryDataGrid.Items.Refresh();

                dataGridProteinDatabases.Items.Refresh();
                proteinDbSummaryDataGrid.Items.Refresh();
            }
        }

        private void FinishedTaskHandler(object sender, SingleTaskEventArgs s)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => FinishedTaskHandler(sender, s)));
            }
            else
            {
                var theTask = InProgressTasks.First(b => b.DisplayName.Equals(s.DisplayName));
                theTask.IsIndeterminate = false;
                theTask.Progress = 100;
                theTask.Status = "Done!";

                dataGridSpectraFiles.Items.Refresh();
                spectraFileSummaryDataGrid.Items.Refresh();

                dataGridProteinDatabases.Items.Refresh();
                proteinDbSummaryDataGrid.Items.Refresh();
            }
        }

        private void AddBranchToTreeViewHandler(object sender, StringEventArgs s)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => AddBranchToTreeViewHandler(sender, s)));
            }
            else
            {
                // Find the task or the collection!!!

                ForTreeView theEntityOnWhichToUpdateLabel = InProgressTasks.First(b => b.Id.Equals(s.NestedIDs[0]));

                for (int i = 1; i < s.NestedIDs.Count - 1; i++)
                {
                    var hm = s.NestedIDs[i];
                    try
                    {
                        theEntityOnWhichToUpdateLabel = theEntityOnWhichToUpdateLabel.Children.First(b => b.Id.Equals(hm));
                    }
                    catch
                    {
                        theEntityOnWhichToUpdateLabel.Children.Add(new CollectionForTreeView(hm, hm));
                        theEntityOnWhichToUpdateLabel = theEntityOnWhichToUpdateLabel.Children.First(b => b.Id.Equals(hm));
                    }
                }

                theEntityOnWhichToUpdateLabel.Children.Add(new CollectionForTreeView(s.S, s.NestedIDs.Last()));
            }
        }

        private void NewoutLabelStatus(object sender, StringEventArgs s)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => NewoutLabelStatus(sender, s)));
            }
            else
            {
                // Find the task or the collection!!!

                ForTreeView theEntityOnWhichToUpdateLabel = InProgressTasks.First(b => b.Id.Equals(s.NestedIDs[0]));

                foreach (var hm in s.NestedIDs.Skip(1))
                {
                    try
                    {
                        theEntityOnWhichToUpdateLabel = theEntityOnWhichToUpdateLabel.Children.First(b => b.Id.Equals(hm));
                    }
                    catch
                    {
                        theEntityOnWhichToUpdateLabel.Children.Add(new CollectionForTreeView(hm, hm));
                        theEntityOnWhichToUpdateLabel = theEntityOnWhichToUpdateLabel.Children.First(b => b.Id.Equals(hm));
                    }
                }

                theEntityOnWhichToUpdateLabel.Status = s.S;
                theEntityOnWhichToUpdateLabel.IsIndeterminate = true;
            }
        }

        /// <summary>
        /// Updates the progress bar for a task/file.
        /// </summary>
        private void NewoutProgressBar(object sender, ProgressEventArgs s)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => NewoutProgressBar(sender, s)));
            }
            else
            {
                ForTreeView theEntityOnWhichToUpdateLabel = InProgressTasks.First(b => b.Id.Equals(s.NestedIDs[0]));

                foreach (var hm in s.NestedIDs.Skip(1))
                {
                    try
                    {
                        theEntityOnWhichToUpdateLabel = theEntityOnWhichToUpdateLabel.Children.First(b => b.Id.Equals(hm));
                    }
                    catch
                    {
                        theEntityOnWhichToUpdateLabel.Children.Add(new CollectionForTreeView(hm, hm));
                        theEntityOnWhichToUpdateLabel = theEntityOnWhichToUpdateLabel.Children.First(b => b.Id.Equals(hm));
                    }
                }

                theEntityOnWhichToUpdateLabel.Status = s.V;
                theEntityOnWhichToUpdateLabel.IsIndeterminate = false;
                theEntityOnWhichToUpdateLabel.Progress = s.NewProgress;
            }
        }

        private void RefreshBetweenTasksHandler(object sender, EventArgs e)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => RefreshBetweenTasksHandler(sender, e)));
            }
            else
            {
                dataGridSpectraFiles.Items.Refresh();
                spectraFileSummaryDataGrid.Items.Refresh();

                dataGridProteinDatabases.Items.Refresh();
                proteinDbSummaryDataGrid.Items.Refresh();
            }
        }

        private void SuccessfullyStartingAllTasks(object sender, EventArgs e)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => SuccessfullyStartingAllTasks(sender, e)));
            }
            else
            {
                dataGridSpectraFiles.Items.Refresh();
                spectraFileSummaryDataGrid.Items.Refresh();

                ToggleEnabledButtonsOnStartOrFinishRun(false);

                RunTasksButton.IsEnabled = false;
                RunTasksButton.Visibility = Visibility.Hidden;
                ResetTasksButton.IsEnabled = false;
                CancelTasksButton.IsEnabled = true;
                ResetTasksButton.Visibility = Visibility.Visible;
                CancelTasksButton.Visibility = Visibility.Visible;
            }
        }

        private void SuccessfullyFinishedAllTasks(object sender, StringEventArgs e)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => SuccessfullyFinishedAllTasks(sender, e)));
            }
            else
            {
                ResetTasksButton.IsEnabled = true;
                CancelTasksButton.IsEnabled = false;

                dataGridSpectraFiles.Items.Refresh();
                spectraFileSummaryDataGrid.Items.Refresh();
            }
        }

        private void FinishedWritingFile(object sender, SingleFileEventArgs v)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => FinishedWritingFile(sender, v)));
            }
            else
            {
                ForTreeView AddWrittenFileToThisOne = InProgressTasks.First(b => b.Id.Equals(v.NestedIDs[0]));

                foreach (var hm in v.NestedIDs.Skip(1))
                {
                    try
                    {
                        AddWrittenFileToThisOne = AddWrittenFileToThisOne.Children.First(b => b.Id.Equals(hm));
                    }
                    catch
                    {
                    }
                }
                AddWrittenFileToThisOne.Children.Add(new OutputFileForTreeView(v.WrittenFile, Path.GetFileName(v.WrittenFile)));
            }
        }

        private void EverythingRunnerExceptionHandler(Task obj)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => EverythingRunnerExceptionHandler(obj)));
            }
            else
            {
                Exception e = obj.Exception;
                while (e.InnerException != null)
                {
                    e = e.InnerException;
                }

                var message = "Run failed, Exception: " + e.Message;
                var messageBoxResult = System.Windows.MessageBox.Show(message + "\n\nWould you like to report this crash?", "Runtime Error", MessageBoxButton.YesNo);
                notificationsTextBox.AppendText(message + Environment.NewLine);
                Exception exception = e;
                //Find Output Folder
                string outputFolder = e.Data["folder"].ToString();
                string tomlText = "";
                if (Directory.Exists(outputFolder))
                {
                    var tomls = Directory.GetFiles(outputFolder, "*.toml");
                    //will only be 1 toml per task
                    foreach (var tomlFile in tomls)
                    {
                        tomlText += "\n" + File.ReadAllText(tomlFile);
                    }

                    if (!tomls.Any())
                    {
                        tomlText = "TOML not found";
                    }
                }
                else
                {
                    tomlText = "Directory not found";
                }

                if (messageBoxResult == MessageBoxResult.Yes)
                {
                    string body = exception.Message + "%0D%0A" + exception.Data +
                       "%0D%0A" + exception.StackTrace +
                       "%0D%0A" + exception.Source +
                       "%0D%0A %0D%0A %0D%0A %0D%0A SYSTEM INFO: %0D%0A " +
                        SystemInfo.CompleteSystemInfo() +
                       "%0D%0A%0D%0A MetaMorpheus: version " + GlobalVariables.MetaMorpheusVersion
                       + "%0D%0A %0D%0A %0D%0A %0D%0A TOML: %0D%0A " +
                       tomlText;
                    body = body.Replace('&', ' ');
                    body = body.Replace("\n", "%0D%0A");
                    body = body.Replace("\r", "%0D%0A");
                    string mailto = string.Format("mailto:{0}?Subject=MetaMorpheus. Issue:&Body={1}", "mm_support@chem.wisc.edu", body);
                    GlobalVariables.StartProcess(mailto);
                    Console.WriteLine(body);
                }

                SuccessfullyFinishedAllTasks(null, null);
            }
        }

        #endregion

        #region Events triggered by user interaction

        /// <summary>
        /// Event fires when a file is dragged-and-dropped into MetaMorpheus.
        /// </summary>
        private void Window_Drop(object sender, DragEventArgs e)
        {
            string[] files = (string[])e.Data.GetData(DataFormats.FileDrop);

            if (files != null)
            {
                AddPreRunFiles(files);
            }
        }

        /// <summary>
        /// Event fires when the "Add Spectra" button is clicked.
        /// </summary>
        private void AddSpectraFile_Click(object sender, RoutedEventArgs e)
        {
            var openPicker = StartOpenFileDialog("Spectra Files(*.raw;*.mzML;*.mgf)|*.raw;*.mzML;*.mgf");

            if (openPicker.ShowDialog() == true)
            {
                AddPreRunFiles(openPicker.FileNames);
            }
        }

        /// <summary>
        /// Event fires when a spectra file is selected.
        /// </summary>
        private void AddSelectedSpectra(object sender, RoutedEventArgs e)
        {
            DataGridRow obj = (DataGridRow)sender;

            RawDataForDataGrid selectedSpectraFile = (RawDataForDataGrid)obj.DataContext;
            SelectedSpectraFiles.Add(selectedSpectraFile);
        }

        /// <summary>
        /// Event fires when a spectra file is deselected.
        /// </summary>
        private void RemoveSelectedSpectra(object sender, RoutedEventArgs e)
        {
            DataGridRow obj = (DataGridRow)sender;
            RawDataForDataGrid deselectedSpectraFile = (RawDataForDataGrid)obj.DataContext;
            SelectedSpectraFiles.Remove(deselectedSpectraFile);
        }

        private void SetFileSpecificParameters_Click(object sender, RoutedEventArgs e)
        {
            if (!SelectedSpectraFiles.Any())
            {
                MessageBox.Show("Please select at least one spectra file.");
                return;
            }

            try
            {
                var dialog = new FileSpecificParametersWindow(SelectedSpectraFiles);
                if (dialog.ShowDialog() == true)
                {
                    var tomlPathsForSelectedFiles = SelectedSpectraFiles.Select(p => Path.Combine(Directory.GetParent(p.FilePath).ToString(), Path.GetFileNameWithoutExtension(p.FileName)) + ".toml").ToList();

                    foreach (var toml in tomlPathsForSelectedFiles)
                    {
                        UpdateFileSpecificParamsDisplay(toml);
                    }
                }
            }
            catch (MetaMorpheusException ex)
            {
                NotificationHandler(null, new StringEventArgs("Problem parsing the file-specific toml; " + ex.Message + "; is the toml from an older version of MetaMorpheus?", null));
            }
            catch (KeyNotFoundException ex)
            {
                NotificationHandler(null, new StringEventArgs("Problem parsing the file-specific toml; " + ex.Message + "; please update the proteases.tsv file and restart MetaMorpheus to use this file-specific toml.", null));
            }
        }

        private void SetExperimentalDesign_Click(object sender, RoutedEventArgs e)
        {
            var dialog = new ExperimentalDesignWindow(SpectraFiles);
            dialog.ShowDialog();
        }

        /// <summary>
        /// Event fires when the "Add Protein Database" button is clicked.
        /// </summary>
        private void AddProteinDatabase_Click(object sender, RoutedEventArgs e)
        {
            var openPicker = StartOpenFileDialog("Database Files|*.xml;*.xml.gz;*.fasta;*.fa;*.msp");

            if (openPicker.ShowDialog() == true)
            {
                AddPreRunFiles(openPicker.FileNames);
            }
        }

        /// <summary>
        /// Event fires when a database file is selected.
        /// </summary>
        private void AddSelectedDatabase(object sender, RoutedEventArgs e)
        {
            DataGridRow obj = (DataGridRow)sender;

            ProteinDbForDataGrid selectedProteinDb = (ProteinDbForDataGrid)obj.DataContext;
            SelectedProteinDatabaseFiles.Add(selectedProteinDb);
        }

        /// <summary>
        /// Event fires when a database file is deselected.
        /// </summary>
        private void RemoveSelectedDatabase(object sender, RoutedEventArgs e)
        {
            DataGridRow obj = (DataGridRow)sender;
            ProteinDbForDataGrid deselectedProteinDbFile = (ProteinDbForDataGrid)obj.DataContext;
            SelectedProteinDatabaseFiles.Remove(deselectedProteinDbFile);
        }

        /// <summary>
        /// Event fires when "Add Contaminants" button is clicked.
        /// </summary>
        private void AddDefaultContaminantDatabase_Click(object sender, RoutedEventArgs e)
        {
            string[] contaminantFiles = Directory.GetFiles(Path.Combine(GlobalVariables.DataDir, "Contaminants"));
            AddPreRunFiles(contaminantFiles);
        }

        /// <summary>
        /// Event fires when the "Set as contaminant" button is clicked.
        /// </summary>
        public void SetSelectedDatabaseAsContaminant_Click(object sender, RoutedEventArgs e)
        {
            foreach (ProteinDbForDataGrid db in SelectedProteinDatabaseFiles)
            {
                db.Contaminant = true;
            }

            dataGridProteinDatabases.Items.Refresh();
            proteinDbSummaryDataGrid.Items.Refresh();
        }

        /// <summary>
        /// Event fires when the "Set as non-contaminant" button is clicked.
        /// </summary>
        public void SetSelectedDatabaseAsNonContaminant_Click(object sender, RoutedEventArgs e)
        {
            foreach (ProteinDbForDataGrid db in SelectedProteinDatabaseFiles)
            {
                db.Contaminant = false;
            }

            dataGridProteinDatabases.Items.Refresh();
            proteinDbSummaryDataGrid.Items.Refresh();
        }

        /// <summary>
        /// Event fires when a data grid row (protein DB or spectra file) is double-clicked.
        /// </summary>
        private void DatabaseOrSpectraFile_DoubleClick(object sender, MouseButtonEventArgs e)
        {
            // prevent opening if a run is in progress
            if (!RunTasksButton.IsEnabled)
            {
                return;
            }

            // user is probably just checking or unchecking a checkbox, don't open the file
            if (sender is DataGridCell cell && cell.Column is DataGridCheckBoxColumn)
            {
                return;
            }

            var dataContext = GetItemDataContext(sender, e).FirstOrDefault();

            // open the file with the default process for this file format
            if (dataContext != null)
            {
                if (dataContext is RawDataForDataGrid spectraFile)
                {
                    OpenFile(spectraFile.FilePath);
                }
                else if (dataContext is ProteinDbForDataGrid proteinDb)
                {
                    OpenFile(proteinDb.FilePath);
                }
            }
        }

        private void AddSearchTaskButton_Click(object sender, RoutedEventArgs e)
        {
            OpenNewTaskWindow(MyTask.Search);
        }

        private void AddCalibrateTaskButton_Click(object sender, RoutedEventArgs e)
        {
            OpenNewTaskWindow(MyTask.Calibrate);
        }

        private void AddGPTMDTaskButton_Click(object sender, RoutedEventArgs e)
        {
            OpenNewTaskWindow(MyTask.Gptmd);
        }

        private void AddCrosslinkTask_Click(object sender, RoutedEventArgs e)
        {
            OpenNewTaskWindow(MyTask.XLSearch);
        }

        private void AddGlycoSearchTask_Click(object sender, RoutedEventArgs e)
        {
            OpenNewTaskWindow(MyTask.GlycoSearch);
        }

        /// <summary>
        /// Event fires when the "Add Task" button is clicked.
        /// </summary>
        private void AddTask_Click(object sender, RoutedEventArgs e)
        {
            ContextMenu contextMenu = this.FindResource("AddTaskMenu") as ContextMenu;
            contextMenu.PlacementTarget = sender as Button;
            contextMenu.IsOpen = true;
        }

        private void LoadTask_Click(object sender, RoutedEventArgs e)
        {
            var openPicker = StartOpenFileDialog("TOML files(*.toml)|*.toml");

            if (openPicker.ShowDialog() == true)
            {
                AddPreRunFiles(openPicker.FileNames);
            }
        }

        private void EditTask_Click(object sender, RoutedEventArgs e)
        {
            var item = GetItemDataContext(sender, e).FirstOrDefault();

            if (item is PreRunTask preRunTask)
            {
                OpenPreRunTaskForEditing(preRunTask);
            }
        }

        /// <summary>
        /// Event fires when the "Save as .toml" context menu item is clicked.
        /// Can occur in the task tree view.
        /// </summary>
        private void SaveTaskAsToml_Click(object sender, RoutedEventArgs e)
        {
            var item = GetItemDataContext(sender, e).FirstOrDefault();
            MetaMorpheusTask task;

            if (item is PreRunTask)
            {
                task = ((PreRunTask)item).metaMorpheusTask;
            }
            else if (item is InRunTask)
            {
                task = ((InRunTask)item).Task;
            }
            else
            {
                // if this message ever appears, it's a bug...
                MessageBox.Show("Unable to save this item as .toml.");
                return;
            }

            string filename = task.CommonParameters.TaskDescriptor + ".toml";

            SaveFileDialog save = new SaveFileDialog { FileName = filename, AddExtension = true, DefaultExt = ".toml" };

            if (save.ShowDialog() == true)
            {
                Toml.WriteFile(task, save.FileName, MetaMorpheusTask.tomlConfig);
            }
        }

        /// <summary>
        /// Deletes the selected item(s).
        /// </summary>
        private void Delete_Click(object sender, RoutedEventArgs e)
        {
            var item = GetItemDataContext(sender, e).FirstOrDefault();

            if (item == null)
            {
                return;
            }

            if (item is PreRunTask)
            {
                PreRunTasks.Remove(item as PreRunTask);
                UpdateGuiOnPreRunChange();
            }
            else if (item is RawDataForDataGrid)
            {
                foreach (var selectedFile in SelectedSpectraFiles.ToList())
                {
                    SpectraFiles.Remove(selectedFile);
                }

                if (!SpectraFiles.Any())
                {
                    OutputFolderTextBox.Text = string.Empty;
                }
            }
            else if (item is ProteinDbForDataGrid)
            {
                foreach (var selectedFile in SelectedProteinDatabaseFiles.ToList())
                {
                    ProteinDatabases.Remove(selectedFile);
                }
            }
        }

        /// <summary>
        /// Event fires when the "reset" button is clicked (to reset the tasks after the run completes).
        /// </summary>
        private void ResetTasks_Click(object sender, RoutedEventArgs e)
        {
            ToggleEnabledButtonsOnStartOrFinishRun(true);

            RunTasksButton.IsEnabled = true;
            RunTasksButton.Visibility = Visibility.Visible;

            tasksTreeView.DataContext = PreRunTasks;
            taskSummary.DataContext = PreRunTasks;
            UpdateGuiOnPreRunChange();

            var pathOfFirstSpectraFile = Path.GetDirectoryName(SpectraFiles.First().FilePath);
            OutputFolderTextBox.Text = Path.Combine(pathOfFirstSpectraFile, @"$DATETIME");
        }

        /// <summary>
        /// Event fires when the "cancel" button is clicked (i.e., cancel the run).
        /// </summary>
        private void CancelTasks_Click(object sender, RoutedEventArgs e)
        {
            string grammar = PreRunTasks.Count <= 1 ? "this task" : "these tasks";
            if (MessageBox.Show("Are you sure you want to cancel " + grammar + "?", "Cancel Tasks", MessageBoxButton.OKCancel) == MessageBoxResult.OK)
            {
                GlobalVariables.StopLoops = true;
                //CancelButton.IsEnabled = false;
                notificationsTextBox.AppendText("Canceling...\n");
            }
        }

        /// <summary>
        /// Moves the task up or down in the GUI.
        /// </summary>
        private void MoveSelectedTask_Click(object sender, RoutedEventArgs e, bool moveTaskUp)
        {
            PreRunTask taskToMove = (PreRunTask)GetItemDataContext(sender, e).FirstOrDefault();

            if (taskToMove == null)
            {
                return;
            }

            int oldPosition = PreRunTasks.IndexOf(taskToMove);
            int newPosition = oldPosition - 1;
            if (!moveTaskUp)
            {
                newPosition = oldPosition + 1;
            }

            if (newPosition >= 0 && newPosition < PreRunTasks.Count)
            {
                var taskThatUsedToBeInThatPosition = PreRunTasks[newPosition];
                PreRunTasks[newPosition] = taskToMove;
                PreRunTasks[oldPosition] = taskThatUsedToBeInThatPosition;

                // renames the tasks in the new correct order (Task1, Task2, etc.)
                UpdateGuiOnPreRunChange();

                var item = tasksTreeView.ItemContainerGenerator.ContainerFromItem(taskToMove);
                if (item != null)
                {
                    ((TreeViewItem)item).IsSelected = true;
                }

                item = taskSummary.ItemContainerGenerator.ContainerFromItem(taskToMove);
                if (item != null)
                {
                    ((TreeViewItem)item).IsSelected = true;
                }
            }
        }

        /// <summary>
        /// Event fires when the "Run all tasks" button is clicked.
        /// </summary>
        private void RunAllTasks_Click(object sender, RoutedEventArgs e)
        {
            GlobalVariables.StopLoops = false;

            // check for valid tasks/spectra files/protein databases
            if (!PreRunTasks.Any())
            {
                NotificationHandler(null, new StringEventArgs("You need to add at least one task!", null));
                return;
            }
            if (!SpectraFiles.Any())
            {
                NotificationHandler(null, new StringEventArgs("You need to add at least one spectra file!", null));
                return;
            }
            if (!ProteinDatabases.Any())
            {
                NotificationHandler(null, new StringEventArgs("You need to add at least one protein database!", null));
                return;
            }

            // check that experimental design is defined if normalization is enabled
            var searchTasks = PreRunTasks
                .Where(p => p.metaMorpheusTask.TaskType == MyTask.Search)
                .Select(p => (SearchTask)p.metaMorpheusTask);

            string pathToExperDesign = Directory.GetParent(SpectraFiles.First().FilePath).FullName;
            pathToExperDesign = Path.Combine(pathToExperDesign, GlobalVariables.ExperimentalDesignFileName);

            if (!File.Exists(pathToExperDesign))
            {
                if (searchTasks.Any(p => p.SearchParameters.Normalize))
                {
                    MessageBox.Show("Experimental design must be defined for normalization!\n" +
                    "Click the \"Set Experimental Design\" button in the the spectra files tab");
                    return;
                }
            }
            else
            {
                ExperimentalDesign.ReadExperimentalDesign(pathToExperDesign, SpectraFiles.Select(p => p.FilePath).ToList(), out var errors);

                if (errors.Any())
                {
                    if (searchTasks.Any(p => p.SearchParameters.Normalize))
                    {
                        MessageBox.Show(errors.First());
                        return;
                    }
                    else
                    {
                        var result = MessageBox.Show("An experimental design file was found, but an error " +
                            "occurred reading it. Do you wish to continue with an empty experimental design?" +
                            "\nThe error was: " + errors.First(), "Error", MessageBoxButton.YesNo);

                        if (result == MessageBoxResult.Yes)
                        {
                            File.Delete(pathToExperDesign);
                        }
                        else
                        {
                            return;
                        }
                    }
                }
            }

            // settings are OK to run
            InProgressTasks = new ObservableCollection<InRunTask>();

            for (int i = 0; i < PreRunTasks.Count; i++)
            {
                InProgressTasks.Add(new InRunTask("Task" + (i + 1) + "-" + PreRunTasks[i].metaMorpheusTask.CommonParameters.TaskDescriptor, PreRunTasks[i].metaMorpheusTask));
            }
            tasksTreeView.DataContext = InProgressTasks;
            taskSummary.DataContext = InProgressTasks;

            notificationsTextBox.Document.Blocks.Clear();

            // set up output folder
            if (string.IsNullOrEmpty(OutputFolderTextBox.Text))
            {
                var pathOfFirstSpectraFile = Path.GetDirectoryName(SpectraFiles.First().FilePath);
                OutputFolderTextBox.Text = Path.Combine(pathOfFirstSpectraFile, @"$DATETIME");
            }

            var startTimeForAllFilenames = DateTime.Now.ToString("yyyy-MM-dd-HH-mm-ss", CultureInfo.InvariantCulture);
            string outputFolder = OutputFolderTextBox.Text.Replace("$DATETIME", startTimeForAllFilenames);
            OutputFolderTextBox.Text = outputFolder;

            // everything is ready to run
            EverythingRunnerEngine a = new EverythingRunnerEngine(InProgressTasks.Select(b => (b.DisplayName, b.Task)).ToList(),
                SpectraFiles.Where(b => b.Use).Select(b => b.FilePath).ToList(),
                ProteinDatabases.Where(b => b.Use).Select(b => new DbForTask(b.FilePath, b.Contaminant)).ToList(),
                outputFolder);

            var t = new Task(a.Run);
            t.ContinueWith(EverythingRunnerExceptionHandler, TaskContinuationOptions.OnlyOnFaulted);
            t.Start();
        }

        /// <summary>
        /// Event fires when an item in the task treeview is double-clicked.
        /// </summary>
        private void TasksTreeView_MouseDoubleClick(object sender, EventArgs e)
        {
            var a = sender as TreeView;
            if (a.SelectedItem is PreRunTask preRunTask)
            {
                OpenPreRunTaskForEditing(preRunTask);
            }
            else if (a.SelectedItem is OutputFileForTreeView writtenFile)
            {
                OpenFile(writtenFile.FullPath);
            }
        }

        /// <summary>
        /// Event fires when the "Open containing item" context menu item is clicked.
        /// Can occur on a protein DB, spectra file, or written file.
        /// </summary>
        private void OpenContainingFolder_Click(object sender, RoutedEventArgs e)
        {
            var path = GetPathOfItem(sender, e);
            OpenFolder(Path.GetDirectoryName(path));
        }

        /// <summary>
        /// Event fires when the "Open file" context menu item is clicked.
        /// Can occur on a protein DB, spectra file, or written file.
        /// </summary>
        private void OpenFile_Click(object sender, RoutedEventArgs e)
        {
            var path = GetPathOfItem(sender, e);
            OpenFile(path);
        }

        /// <summary>
        /// Event fires when the "Open" button is clicked (referring to the output folder).
        /// </summary>
        private void OpenOutputFolder_Click(object sender, RoutedEventArgs e)
        {
            string outputFolder = OutputFolderTextBox.Text;
            if (outputFolder.Contains("$DATETIME"))
            {
                // the exact file path isn't known yet, so just open the parent directory
                outputFolder = Directory.GetParent(outputFolder).FullName;
            }

            // create the directory if it doesn't exist yet
            if (!Directory.Exists(outputFolder) && !string.IsNullOrEmpty(outputFolder))
            {
                try
                {
                    Directory.CreateDirectory(outputFolder);
                }
                catch (Exception ex)
                {
                    NotificationHandler(null, new StringEventArgs("Error opening directory: " + ex.Message, null));
                }
            }

            OpenFolder(outputFolder);
        }

        /// <summary>
        /// Handles keyboard input.
        /// Note: Window_KeyDown is NOT used because it does not seem to capture keyboard input from a DataGrid.
        /// </summary>
        private void Window_PreviewKeyDown(object sender, KeyEventArgs e)
        {
            if (RunTasksButton.IsEnabled)
            {
                // delete selected task/db/spectra
                if (e.Key == Key.Delete || e.Key == Key.Back)
                {
                    Delete_Click(sender, e);
                    e.Handled = true;
                }

                // move task down
                if (e.Key == Key.Add || e.Key == Key.OemPlus)
                {
                    MoveSelectedTask_Click(sender, e, false);
                    e.Handled = true;
                }

                // move task up
                if (e.Key == Key.Subtract || e.Key == Key.OemMinus)
                {
                    MoveSelectedTask_Click(sender, e, true);
                    e.Handled = true;
                }
            }
        }

        /// <summary>
        /// Event fires when MetaMorpheus is closed.
        /// </summary>
        private void MainWindow_Closing(object sender, CancelEventArgs e)
        {
            if (UpdateGUISettings.Params.AskBeforeExitingMetaMorpheus && !GlobalVariables.MetaMorpheusVersion.Contains("DEBUG"))
            {
                var exit = ExitMsgBox.Show("Exit MetaMorpheus", "Are you sure you want to exit MetaMorpheus?", "Yes", "No", "Yes and don't ask me again");

                if (exit == MessageBoxResult.Yes) // yes, exit MetaMorpheus
                {
                    e.Cancel = false;
                }
                else if (exit == MessageBoxResult.OK) // yes and don't ask me again
                {
                    UpdateGUISettings.Params.AskBeforeExitingMetaMorpheus = false;
                    Toml.WriteFile(UpdateGUISettings.Params, Path.Combine(GlobalVariables.DataDir, @"GUIsettings.toml"), MetaMorpheusTask.tomlConfig);
                    e.Cancel = false;
                }
                else // no, do not exit MetaMorpheus
                {
                    e.Cancel = true;
                }
            }
        }

        /// <summary>
        /// Event fires when the tab item (on the left side) is changed.
        /// </summary>
        private void TabControl_SelectionChanged(object sender, SelectionChangedEventArgs e)
        {
            // this "if" statement that checks the type of the sender is here because this event can somehow 
            // be mistakenly triggered by other events (see FlashLFQ's GUI)
            var senderType = e.OriginalSource.GetType().Name;

            if (senderType == "TabControl")
            {
                var selectedItem = (TabItem)MainWindowTabControl.SelectedItem;
                var selectedItemHeader = selectedItem.Header.ToString();

                if (selectedItemHeader == "Visualize")
                {
                    MenuItem_MetaDraw_Click(sender, e);
                }
            }
        }

        private void MenuItem_YouTube_Click(object sender, RoutedEventArgs e)
        {
            GlobalVariables.StartProcess(@"https://www.youtube.com/playlist?list=PLVk5tTSZ1aWlYiTvJbRj6hjVDq4qH3w__");
        }

        private void MenuItem_ProteomicsNewsBlog_Click(object sender, RoutedEventArgs e)
        {
            GlobalVariables.StartProcess(@"https://proteomicsnews.blogspot.com/");
        }

        private void MenuItem_UpdateMetaMorpheus_Click(object sender, RoutedEventArgs e)
        {
            UpdateMetaMorpheus(printMessageIfThisIsLatestVersion: true);
        }

        /// <summary>
        /// Opens the requested URL with the user's web browser.
        /// </summary>
        private void Hyperlink_RequestNavigate(object sender, RequestNavigateEventArgs e)
        {
            GlobalVariables.StartProcess(e.Uri.ToString());
        }

        private void MenuItem_EmailHelp_Click(object sender, RequestNavigateEventArgs e)
        {
            string mailto = string.Format("mailto:{0}?Subject=MetaMorpheus. Issue:", "mm_support@chem.wisc.edu");
            GlobalVariables.StartProcess(mailto);
        }

        private void MenuItem_Twitter_Click(object sender, RoutedEventArgs e)
        {
            GlobalVariables.StartProcess(@"https://twitter.com/Smith_Chem_Wisc");
        }

        private void MenuItem_Slack_Click(object sender, RoutedEventArgs e)
        {
            GlobalVariables.StartProcess(@"https://join.slack.com/t/smith-chem-public/shared_invite/enQtNDYzNTM5Mzg5NzY0LTRiYWQ5MzVmYmExZWIyMTcyZmNlODJjMWI0YjVhNGM2MmQ2NjE4ZDAzNmM4NWYxMDFhNTQyNDBiM2E0MWE0NGU");
        }

        private void MenuItem_Proxl_Click(object sender, RoutedEventArgs e)
        {
            GlobalVariables.StartProcess(@"http://proxl-ms.org/");
        }

        private void MenuItem_Spritz_Click(object sender, RoutedEventArgs e)
        {
            GlobalVariables.StartProcess(@"https://smith-chem-wisc.github.io/Spritz/");
        }

        private void MenuItem_OpenDataDir_Click(object sender, RoutedEventArgs e)
        {
            OpenFolder(GlobalVariables.DataDir);
        }

        private void MenuItem_OpenSettings_Click(object sender, RoutedEventArgs e)
        {
            GlobalVariables.StartProcess(Path.Combine(GlobalVariables.DataDir, @"settings.toml"), useNotepadToOpenToml: true);
            Application.Current.Shutdown();
        }

        private void MenuItem_GuiSettings_Click(object sender, RoutedEventArgs e)
        {
            GlobalVariables.StartProcess(Path.Combine(GlobalVariables.DataDir, @"GUIsettings.toml"), useNotepadToOpenToml: true);
            Application.Current.Shutdown();
        }

        private void MenuItem_MetaDraw_Click(object sender, RoutedEventArgs e)
        {
            MetaDraw metaDrawGui = new MetaDraw();
            metaDrawGui.Show();
        }

        private void AddCustomMod_Click(object sender, RoutedEventArgs e)
        {
            var dialog = new CustomModWindow();
            dialog.ShowDialog();
        }

        private void AddCustomAminoAcid_Click(object sender, RoutedEventArgs e)
        {
            var dialog = new CustomAminoAcidWindow();
            dialog.ShowDialog();
        }

        private void AddCustomCrosslinker_Click(object sender, RoutedEventArgs e)
        {
            var dialog = new CustomCrosslinkerWindow();
            dialog.ShowDialog();
        }

        private void MoveTaskUp_Click(object sender, RoutedEventArgs e)
        {
            var item = GetItemDataContext(sender, e).FirstOrDefault();
            MoveSelectedTask_Click(item, null, moveTaskUp: true);
        }

        private void MoveTaskDown_Click(object sender, RoutedEventArgs e)
        {
            var item = GetItemDataContext(sender, e).FirstOrDefault();
            MoveSelectedTask_Click(item, null, moveTaskUp: false);
        }

        private void DeleteAll_Click(object sender, RoutedEventArgs e)
        {
            var item = GetItemDataContext(sender, e).FirstOrDefault();

            if (item is PreRunTask task)
            {
                PreRunTasks.Clear();
            }
            else if (item is ProteinDbForDataGrid db)
            {
                ProteinDatabases.Clear();
            }
            else if (item is RawDataForDataGrid spectra)
            {
                SpectraFiles.Clear();
            }
        }

        #endregion

        #region Helper methods called by events

        /// <summary>
        /// Opens a file with the specified path.
        /// </summary>
        private void OpenFile(string filePath)
        {
            if (File.Exists(filePath))
            {
                GlobalVariables.StartProcess(filePath);
            }
            else
            {
                MessageBox.Show("File does not exist: " + filePath);
            }
        }

        /// <summary>
        /// Opens a folder with the specified path.
        /// </summary>
        private void OpenFolder(string folderPath)
        {
            if (Directory.Exists(folderPath))
            {
                // open the directory
                System.Diagnostics.Process.Start(new System.Diagnostics.ProcessStartInfo()
                {
                    FileName = folderPath,
                    UseShellExecute = true,
                    Verb = "open"
                });
            }
            else
            {
                MessageBox.Show("Folder does not exist: " + folderPath);
            }
        }

        /// <summary>
        /// Checks for a MetaMorpheus update via the Internet.
        /// </summary>
        private void UpdateMetaMorpheus(bool printMessageIfThisIsLatestVersion = false)
        {
            try
            {
                NewestKnownMetaMorpheusVersion = MetaUpdater.GetVersionNumbersFromWeb();

                if (NewestKnownMetaMorpheusVersion == null)
                {
                    throw new MetaMorpheusException("Web connection appears to be functional but something else went wrong");
                }
            }
            catch (Exception e)
            {
                NotificationHandler(null, new StringEventArgs("Could not get newest version from web: " + e.Message, null));
                return;
            }

            if (!GlobalVariables.MetaMorpheusVersion.Equals(NewestKnownMetaMorpheusVersion))
            {
                try
                {
                    MetaUpdater newwind = new MetaUpdater();
                    newwind.ShowDialog();
                }
                catch (Exception ex)
                {
                    MessageBox.Show(ex.ToString());
                }
            }
            else if (printMessageIfThisIsLatestVersion)
            {
                MessageBox.Show("You have the most updated version!");
            }
        }

        /// <summary>
        /// This text is displayed the first time a user opens MetaMorpheus.
        /// </summary>
        private FlowDocument GetWelcomeText()
        {
            FlowDocument doc = notificationsTextBox.Document;
            Paragraph p = new Paragraph();
            Run run1 = new Run("Visit our ");
            Run run2 = new Run("Wiki");
            Run run3 = new Run(" or ");
            Run run4 = new Run("Youtube channel");
            Run run5 = new Run(" to check out what MetaMorpheus can do!" + System.Environment.NewLine);

            Hyperlink wikiLink = new Hyperlink(run2);
            wikiLink.NavigateUri = new Uri(@"https://github.com/smith-chem-wisc/MetaMorpheus/wiki");

            Hyperlink youtubeLink = new Hyperlink(run4);
            youtubeLink.NavigateUri = new Uri(@"https://www.youtube.com/playlist?list=PLVk5tTSZ1aWlhNPh7jxPQ8pc0ElyzSUQb");

            var links = new List<Hyperlink> { wikiLink, youtubeLink };

            p.Inlines.Add(run1);
            p.Inlines.Add(wikiLink);
            p.Inlines.Add(run3);
            p.Inlines.Add(youtubeLink);
            p.Inlines.Add(run5);

            foreach (Hyperlink link in links)
            {
                link.RequestNavigate += (sender, e) =>
                {
                    GlobalVariables.StartProcess(e.Uri.ToString());
                };
            }

            doc.Blocks.Add(p);
            return doc;
        }

        private void UpdateGuiOnPreRunChange()
        {
            if (PreRunTasks.Count > 0)
            {
                for (int i = 0; i < PreRunTasks.Count; i++)
                {
                    string newName = "Task" + (i + 1) + "-" + PreRunTasks[i].metaMorpheusTask.CommonParameters.TaskDescriptor;
                    PreRunTasks[i].DisplayName = newName;
                }
                tasksTreeView.Items.Refresh();
                taskSummary.Items.Refresh();
            }

            dataGridSpectraFiles.CommitEdit(DataGridEditingUnit.Row, true);
            spectraFileSummaryDataGrid.CommitEdit(DataGridEditingUnit.Row, true);

            dataGridProteinDatabases.CommitEdit(DataGridEditingUnit.Row, true);
            proteinDbSummaryDataGrid.CommitEdit(DataGridEditingUnit.Row, true);

            dataGridSpectraFiles.Items.Refresh();
            spectraFileSummaryDataGrid.Items.Refresh();

            dataGridProteinDatabases.Items.Refresh();
            proteinDbSummaryDataGrid.Items.Refresh();

            if (RunTasksButton.IsEnabled)
            {
                ResetTasksButton.Visibility = Visibility.Hidden;
                CancelTasksButton.Visibility = Visibility.Hidden;
            }
        }

        /// <summary>
        /// Called when a new spectra file is added.
        /// Checks for existing file-specific .toml file and displays its contents in the GUI.
        /// </summary>
        private void UpdateFileSpecificParamsDisplay(string possibleTomlLocation)
        {
            string correspondingSpectraFileForToml = Path.Combine(Directory.GetParent(possibleTomlLocation).ToString(), Path.GetFileNameWithoutExtension(possibleTomlLocation));

            foreach (var spectraFile in SpectraFiles)
            {
                string spectraFilePathWithoutExtension = Path.Combine(Directory.GetParent(spectraFile.FilePath).ToString(), Path.GetFileNameWithoutExtension(spectraFile.FilePath));

                if (correspondingSpectraFileForToml.Equals(spectraFilePathWithoutExtension))
                {
                    if (File.Exists(possibleTomlLocation))
                    {
                        try
                        {
                            // parse to make sure toml is readable
                            TomlTable fileSpecificSettings = Toml.ReadFile(possibleTomlLocation, MetaMorpheusTask.tomlConfig);
                            var temp = new FileSpecificParameters(fileSpecificSettings);

                            // toml is ok; display the file-specific settings in the gui
                            spectraFile.SetParametersText(File.ReadAllText(possibleTomlLocation));
                        }
                        catch (MetaMorpheusException e)
                        {
                            NotificationHandler(null, new StringEventArgs("Problem parsing the file-specific toml " + Path.GetFileName(possibleTomlLocation) + "; " + e.Message + "; is the toml from an older version of MetaMorpheus?", null));
                        }
                        catch (KeyNotFoundException e)
                        {
                            NotificationHandler(null, new StringEventArgs("Problem parsing the file-specific toml " + Path.GetFileName(possibleTomlLocation) + "; " + e.Message + "; please update the proteases.tsv file and restart MetaMorpheus to use this file-specific toml.", null));
                        }
                    }
                    else
                    {
                        // file does not have a file-specific toml; set its displayed file-specific settings to null
                        spectraFile.SetParametersText(null);
                    }
                }
            }

            UpdateGuiOnPreRunChange();
        }

        private void PrintErrorsReadingMods()
        {
            // print any error messages reading the mods to the notifications area
            foreach (var error in GlobalVariables.ErrorsReadingMods)
            {
                NotificationHandler(null, new StringEventArgs(error, null));
            }
            GlobalVariables.ErrorsReadingMods.Clear();
        }

        private void UpdateOutputFolderTextbox()
        {
            if (SpectraFiles.Any())
            {
                // if current output folder is blank and there is a spectra file, use the spectra file's path as the output path
                if (string.IsNullOrWhiteSpace(OutputFolderTextBox.Text))
                {
                    var pathOfFirstSpectraFile = Path.GetDirectoryName(SpectraFiles.First().FilePath);
                    OutputFolderTextBox.Text = Path.Combine(pathOfFirstSpectraFile, @"$DATETIME");
                }
                // else do nothing (do not override if there is a path already there; might clear user-defined path)
            }
            else
            {
                // no spectra files; clear the output folder from the GUI
                OutputFolderTextBox.Clear();
            }
        }

        private bool DatabaseExists(ObservableCollection<ProteinDbForDataGrid> pDOC, ProteinDbForDataGrid uuu)
        {
            foreach (ProteinDbForDataGrid pdoc in pDOC)
            {
                if (pdoc.FilePath == uuu.FilePath) { return true; }
            }

            return false;
        }

        private bool SpectraFileExists(ObservableCollection<RawDataForDataGrid> rDOC, RawDataForDataGrid zzz)
        {
            foreach (RawDataForDataGrid rdoc in rDOC)
            {
                if (rdoc.FileName == zzz.FileName) { return true; }
            }

            return false;
        }

        private void AddPreRunFiles(IEnumerable<string> paths)
        {
            foreach (string path in paths.OrderBy(p => Path.GetFileName(p)))
            {
                if (Directory.Exists(path))
                {
                    foreach (string file in Directory.EnumerateFiles(path, "*.*", SearchOption.AllDirectories))
                    {
                        AddPreRunFile(file);
                    }
                }
                else if (File.Exists(path))
                {
                    AddPreRunFile(path);
                }
            }

            UpdateGuiOnPreRunChange();
        }

        private void AddPreRunFile(string filePath)
        {
            if (!RunTasksButton.IsEnabled)
            {
                return;
            }

            // this line is NOT used because .xml.gz (extensions with two dots) mess up with Path.GetExtension
            //var theExtension = Path.GetExtension(draggedFilePath).ToLowerInvariant();

            // we need to get the filename before parsing out the extension because if we assume that everything after the dot
            // is the extension and there are dots in the file path (i.e. in a folder name), this will mess up
            var filename = Path.GetFileName(filePath);
            var theExtension = Path.GetExtension(filename).ToLowerInvariant();
            bool compressed = theExtension.EndsWith("gz"); // allows for .bgz and .tgz, too which are used on occasion
            theExtension = compressed ? Path.GetExtension(Path.GetFileNameWithoutExtension(filename)).ToLowerInvariant() : theExtension;

            switch (theExtension)
            {
                case ".raw":
                    if (!GlobalVariables.GlobalSettings.UserHasAgreedToThermoRawFileReaderLicence)
                    {
                        // open the Thermo RawFileReader licence agreement
                        var thermoLicenceWindow = new ThermoLicenceAgreementWindow();
                        thermoLicenceWindow.LicenceText.AppendText(ThermoRawFileReaderLicence.ThermoLicenceText);
                        var dialogResult = thermoLicenceWindow.ShowDialog();

                        var newGlobalSettings = new GlobalSettings
                        {
                            UserHasAgreedToThermoRawFileReaderLicence = dialogResult.Value,
                            WriteExcelCompatibleTSVs = GlobalVariables.GlobalSettings.WriteExcelCompatibleTSVs
                        };

                        Toml.WriteFile<GlobalSettings>(newGlobalSettings, Path.Combine(GlobalVariables.DataDir, @"settings.toml"));
                        GlobalVariables.GlobalSettings = newGlobalSettings;

                        // user declined agreement
                        if (!GlobalVariables.GlobalSettings.UserHasAgreedToThermoRawFileReaderLicence)
                        {
                            return;
                        }
                    }

                    goto case ".mzml";

                case ".mgf":
                    NotificationHandler(null, new StringEventArgs(".mgf files lack MS1 spectra, which are needed for quantification and searching for coisolated peptides. All other features of MetaMorpheus will function.", null));
                    goto case ".mzml";

                case ".mzml":
                    if (compressed) // not implemented yet
                    {
                        NotificationHandler(null, new StringEventArgs("Cannot read, try uncompressing: " + filePath, null));
                        break;
                    }
                    RawDataForDataGrid zz = new RawDataForDataGrid(filePath);
                    if (!SpectraFileExists(SpectraFiles, zz))
                    {
                        SpectraFiles.Add(zz);
                    }
                    UpdateFileSpecificParamsDisplay(Path.ChangeExtension(filePath, ".toml"));
                    UpdateOutputFolderTextbox();
                    break;

                case ".xml":
                case ".fasta":
                case ".fa":
                case ".msp":
                    ProteinDbForDataGrid uu = new ProteinDbForDataGrid(filePath);
                    if (!DatabaseExists(ProteinDatabases, uu))
                    {
                        ProteinDatabases.Add(uu);
                        if (theExtension.Equals(".xml"))
                        {
                            try
                            {
                                GlobalVariables.AddMods(UsefulProteomicsDatabases.ProteinDbLoader.GetPtmListFromProteinXml(filePath).OfType<Modification>(), true);

                                PrintErrorsReadingMods();
                            }
                            catch (Exception ee)
                            {
                                MessageBox.Show(ee.ToString());
                                NotificationHandler(null, new StringEventArgs("Cannot parse modification info from: " + filePath, null));
                                ProteinDatabases.Remove(uu);
                            }
                        }
                    }
                    break;

                case ".toml":
                    TomlTable tomlFile = null;
                    try
                    {
                        tomlFile = Toml.ReadFile(filePath, MetaMorpheusTask.tomlConfig);
                    }
                    catch (Exception)
                    {
                        NotificationHandler(null, new StringEventArgs("Cannot read toml: " + filePath, null));
                        break;
                    }

                    if (tomlFile.ContainsKey("TaskType"))
                    {
                        try
                        {
                            switch (tomlFile.Get<string>("TaskType"))
                            {
                                case "Search":
                                    var search = Toml.ReadFile<SearchTask>(filePath, MetaMorpheusTask.tomlConfig);
                                    AddTaskToCollection(search);
                                    break;

                                case "Calibrate":
                                    var calib = Toml.ReadFile<CalibrationTask>(filePath, MetaMorpheusTask.tomlConfig);
                                    AddTaskToCollection(calib);
                                    break;

                                case "Gptmd":
                                    var gptmd = Toml.ReadFile<GptmdTask>(filePath, MetaMorpheusTask.tomlConfig);
                                    AddTaskToCollection(gptmd);
                                    break;

                                case "XLSearch":
                                    var xl = Toml.ReadFile<XLSearchTask>(filePath, MetaMorpheusTask.tomlConfig);
                                    AddTaskToCollection(xl);
                                    break;

                                case "GlycoSearch":
                                    var glyco = Toml.ReadFile<GlycoSearchTask>(filePath, MetaMorpheusTask.tomlConfig);
                                    AddTaskToCollection(glyco);
                                    break;
                            }
                        }
                        catch (Exception e)
                        {
                            NotificationHandler(null, new StringEventArgs("Cannot read task toml: " + e.Message, null));
                        }
                    }
                    break;

                default:
                    NotificationHandler(null, new StringEventArgs("Unrecognized file type: " + theExtension, null));
                    break;
            }
        }

        private void AddTaskToCollection(MetaMorpheusTask taskToAdd)
        {
            PreRunTasks.Add(new PreRunTask(taskToAdd));
            UpdateGuiOnPreRunChange();
        }

        private OpenFileDialog StartOpenFileDialog(string filter)
        {
            OpenFileDialog openFileDialog = new OpenFileDialog
            {
                Filter = filter,
                FilterIndex = 1,
                RestoreDirectory = true,
                Multiselect = true
            };

            return openFileDialog;
        }

        private void OpenNewTaskWindow(MyTask taskType)
        {
            Window dialog = null;
            MetaMorpheusTask task = null;
            string defaultTomlName = null;

            // determine if there is a default .toml for this task
            switch (taskType)
            {
                case MyTask.Search: defaultTomlName = "SearchTaskDefault.toml"; break;
                case MyTask.Calibrate: defaultTomlName = "CalibrationTaskDefault.toml"; break;
                case MyTask.Gptmd: defaultTomlName = "GptmdTaskDefault.toml"; break;
                case MyTask.XLSearch: defaultTomlName = "XLSearchTaskDefault.toml"; break;
                case MyTask.GlycoSearch: defaultTomlName = "GlycoSearchTaskDefault.toml"; break;
            }

            string defaultTomlFilePath = Path.Combine(GlobalVariables.DataDir, "DefaultParameters", defaultTomlName);

            if (File.Exists(defaultTomlFilePath))
            {
                try
                {
                    switch (taskType)
                    {
                        case MyTask.Search: task = Toml.ReadFile<SearchTask>(defaultTomlFilePath, MetaMorpheusTask.tomlConfig); break;
                        case MyTask.Calibrate: task = Toml.ReadFile<CalibrationTask>(defaultTomlFilePath, MetaMorpheusTask.tomlConfig); break;
                        case MyTask.Gptmd: task = Toml.ReadFile<GptmdTask>(defaultTomlFilePath, MetaMorpheusTask.tomlConfig); break;
                        case MyTask.XLSearch: task = Toml.ReadFile<XLSearchTask>(defaultTomlFilePath, MetaMorpheusTask.tomlConfig); break;
                        case MyTask.GlycoSearch: task = Toml.ReadFile<GlycoSearchTask>(defaultTomlFilePath, MetaMorpheusTask.tomlConfig); break;
                    }
                }
                catch (Exception)
                {
                    NotificationHandler(null, new StringEventArgs("Cannot read toml: " + defaultTomlFilePath, null));
                }
            }

            // open the new task window
            switch (taskType)
            {
                case MyTask.Search: dialog = new SearchTaskWindow((SearchTask)task); break;
                case MyTask.Calibrate: dialog = new CalibrateTaskWindow((CalibrationTask)task); break;
                case MyTask.Gptmd: dialog = new GptmdTaskWindow((GptmdTask)task); break;
                case MyTask.XLSearch: dialog = new XLSearchTaskWindow((XLSearchTask)task); break;
                case MyTask.GlycoSearch: dialog = new GlycoSearchTaskWindow((GlycoSearchTask)task); break;
            }

            // save the task to the task collection
            if (dialog.ShowDialog() == true)
            {
                switch (taskType)
                {
                    case MyTask.Search: AddTaskToCollection(((SearchTaskWindow)dialog).TheTask); break;
                    case MyTask.Calibrate: AddTaskToCollection(((CalibrateTaskWindow)dialog).TheTask); break;
                    case MyTask.Gptmd: AddTaskToCollection(((GptmdTaskWindow)dialog).TheTask); break;
                    case MyTask.XLSearch: AddTaskToCollection(((XLSearchTaskWindow)dialog).TheTask); break;
                    case MyTask.GlycoSearch: AddTaskToCollection(((GlycoSearchTaskWindow)dialog).TheTask); break;
                }

                UpdateGuiOnPreRunChange();
            }
        }

        private IEnumerable<object> GetItemDataContext(object sender, RoutedEventArgs e)
        {
            if (sender is PreRunTask || sender is FinishedFileForDataGrid || sender is InRunTask || sender is RawDataForDataGrid || sender is ProteinDbForDataGrid)
            {
                yield return sender;
            }
            else if (sender is MenuItem)
            {
                var menuItem = (MenuItem)sender;
                var dataContext = (ContextMenu)menuItem.Parent;

                if (dataContext.PlacementTarget is DataGridRow dataGridRow)
                {
                    yield return dataGridRow.Item;
                }
                else if (dataContext.PlacementTarget is TreeViewItem treeViewItem)
                {
                    yield return treeViewItem.Header;
                }
                else
                {
                    yield return menuItem.DataContext;
                }
            }
            else if (sender is TreeViewItem treeViewItem)
            {
                yield return treeViewItem.Header;
            }
            else if (e.Source is DataGrid grid)
            {
                foreach (var item in grid.Items)
                {
                    yield return GetItemDataContext(item, null).FirstOrDefault();
                }
            }
            else if (sender is DataGridRow dataGridRow)
            {
                yield return dataGridRow.Item;
            }
            else if (e.Source is DataGridCell cell)
            {
                yield return cell.DataContext;
            }
            else if (e.Source is TreeViewItem treeViewItem2)
            {
                yield return treeViewItem2.Header;
            }
            else if (e.OriginalSource is DataGridCell cell2)
            {
                yield return cell2.DataContext;
            }
            else if (e.OriginalSource is TreeViewItem treeViewItem3)
            {
                yield return treeViewItem3.Header;
            }

            yield return null;
        }

        private string GetPathOfItem(object sender, RoutedEventArgs e)
        {
            // right now this will only get one filepath... could change multiple if multi-selected
            var item = GetItemDataContext(sender, e).FirstOrDefault();

            string path = null;

            if (item is ProteinDbForDataGrid db)
            {
                path = db.FilePath;
            }
            else if (item is RawDataForDataGrid spectra)
            {
                path = spectra.FilePath;
            }
            else if (item is OutputFileForTreeView writtenFile)
            {
                path = writtenFile.FullPath;
            }

            return path;
        }

        private void OpenPreRunTaskForEditing(PreRunTask preRunTask)
        {
            switch (preRunTask.metaMorpheusTask.TaskType)
            {
                case MyTask.Search:

                    var searchDialog = new SearchTaskWindow(preRunTask.metaMorpheusTask as SearchTask);
                    searchDialog.ShowDialog();
                    break;

                case MyTask.Gptmd:
                    var gptmddialog = new GptmdTaskWindow(preRunTask.metaMorpheusTask as GptmdTask);
                    gptmddialog.ShowDialog();
                    break;

                case MyTask.Calibrate:
                    var calibratedialog = new CalibrateTaskWindow(preRunTask.metaMorpheusTask as CalibrationTask);
                    calibratedialog.ShowDialog();
                    break;

                case MyTask.XLSearch:
                    var XLSearchdialog = new XLSearchTaskWindow(preRunTask.metaMorpheusTask as XLSearchTask);
                    XLSearchdialog.ShowDialog();
                    break;

                case MyTask.GlycoSearch:
                    var GlycoSearchdialog = new GlycoSearchTaskWindow(preRunTask.metaMorpheusTask as GlycoSearchTask);
                    GlycoSearchdialog.ShowDialog();
                    break;
            }

            UpdateGuiOnPreRunChange();
        }

        private void ToggleEnabledButtonsOnStartOrFinishRun(bool enable)
        {
            foreach (MenuItem item in ((ContextMenu)this.Resources["ProteinDatabaseContextMenu"]).Items)
            {
                switch (item.Header.ToString())
                {
                    case "Set as contaminant database": item.IsEnabled = enable; break;
                    case "Set as non-contaminant database": item.IsEnabled = enable; break;
                    case "Open file": item.IsEnabled = enable; break;
                    case "Delete": item.IsEnabled = enable; break;
                    case "Delete all": item.IsEnabled = enable; break;
                }
            }

            foreach (MenuItem item in ((ContextMenu)this.Resources["SpectraFileContextMenu"]).Items)
            {
                switch (item.Header.ToString())
                {
                    case "Set file-specific parameters": item.IsEnabled = enable; break;
                    case "Open file": item.IsEnabled = enable; break;
                    case "Delete": item.IsEnabled = enable; break;
                    case "Delete all": item.IsEnabled = enable; break;
                }
            }

            foreach (MenuItem item in ((ContextMenu)this.Resources["TaskContextMenu"]).Items)
            {
                switch (item.Header.ToString())
                {
                    case "Move task up": item.IsEnabled = enable; break;
                    case "Move task down": item.IsEnabled = enable; break;
                    case "Edit task": item.IsEnabled = enable; break;
                    case "Delete": item.IsEnabled = enable; break;
                    case "Delete all": item.IsEnabled = enable; break;
                }
            }

            dataGridProteinDatabases.IsReadOnly = !enable;
            dataGridSpectraFiles.IsReadOnly = !enable;
            spectraFileSummaryDataGrid.IsReadOnly = !enable;
            proteinDbSummaryDataGrid.IsReadOnly = !enable;

            AddDatabaseButton.IsEnabled = enable;
            AddDefaultContaminantsButton.IsEnabled = enable;
            AddSpectraButton.IsEnabled = enable;
            SetFileSpecificSettingsButton.IsEnabled = enable;
            SetExperimentalDesignButton.IsEnabled = enable;
            AddSearchTaskButton.IsEnabled = enable;
            AddCalibTaskButton.IsEnabled = enable;
            AddGptmdTaskButton.IsEnabled = enable;
            AddXlTaskButton.IsEnabled = enable;
            AddGlycoTaskButton.IsEnabled = enable;
            MiniAddProteinDbButton.IsEnabled = enable;
            MiniAddSpectraButton.IsEnabled = enable;
            MiniAddTaskButton.IsEnabled = enable;
            OutputFolderTextBox.IsEnabled = enable;
            RunTasksButton.IsEnabled = enable;
        }

        #endregion

        private void DownloadUniProtDatabase_Click(object sender, RoutedEventArgs e)
        {
            DownloadUniProtDatabaseWindow uniProtDatabaseWindow = new DownloadUniProtDatabaseWindow(ProteinDatabases);
            uniProtDatabaseWindow.WindowStartupLocation = WindowStartupLocation.CenterScreen;            
            uniProtDatabaseWindow.Show();
            uniProtDatabaseWindow.Closed += UniProtDatabaseWindow_Closed;         
        }

        private void UniProtDatabaseWindow_Closed(object sender, EventArgs e)
        {
            dataGridProteinDatabases.Items.Refresh();
            proteinDbSummaryDataGrid.Items.Refresh();
        }

        private void OpenProteomesFolder_Click(object sender, RoutedEventArgs e)
        {
            OpenFolder(Path.Combine(GlobalVariables.DataDir, @"Proteomes"));
        }
    }
}