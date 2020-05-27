using EngineLayer;
using Microsoft.Win32;
using MzLibUtil;
using Nett;
using Newtonsoft.Json.Linq;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Net.Http;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using TaskLayer;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        private readonly ObservableCollection<RawDataForDataGrid> SpectraFilesObservableCollection = new ObservableCollection<RawDataForDataGrid>();
        private readonly ObservableCollection<ProteinDbForDataGrid> ProteinDbObservableCollection = new ObservableCollection<ProteinDbForDataGrid>();
        private readonly ObservableCollection<PreRunTask> StaticTasksObservableCollection = new ObservableCollection<PreRunTask>();
        private readonly ObservableCollection<RawDataForDataGrid> SelectedRawFiles = new ObservableCollection<RawDataForDataGrid>();
        private ObservableCollection<InRunTask> DynamicTasksObservableCollection;

        public static string NewestKnownMetaMorpheusVersion { get; private set; }

        public MainWindow()
        {
            InitializeComponent();

            Title = "MetaMorpheus: version " + GlobalVariables.MetaMorpheusVersion;

            dataGridProteinDatabases.DataContext = ProteinDbObservableCollection;
            dataGridSpectraFiles.DataContext = SpectraFilesObservableCollection;
            tasksTreeView.DataContext = StaticTasksObservableCollection;
            proteinDbSummaryDataGrid.DataContext = ProteinDbObservableCollection;
            spectraFileSummaryDataGrid.DataContext = SpectraFilesObservableCollection;
            taskSummaryDataGrid.DataContext = StaticTasksObservableCollection;

            EverythingRunnerEngine.NewDbsHandler += AddNewProteinDatabaseFromGptmd;
            EverythingRunnerEngine.NewSpectrasHandler += AddNewSpectraFileFromCalibration;
            EverythingRunnerEngine.NewFileSpecificTomlHandler += AddNewFileSpecificTomlFromCalibration;
            EverythingRunnerEngine.StartingAllTasksEngineHandler += NewSuccessfullyStartingAllTasks;
            EverythingRunnerEngine.FinishedAllTasksEngineHandler += NewSuccessfullyFinishedAllTasks;
            EverythingRunnerEngine.WarnHandler += NotificationHandler;
            EverythingRunnerEngine.FinishedWritingAllResultsFileHandler += EverythingRunnerEngine_FinishedWritingAllResultsFileHandler;

            MetaMorpheusTask.StartingSingleTaskHander += Po_startingSingleTaskHander;
            MetaMorpheusTask.FinishedSingleTaskHandler += Po_finishedSingleTaskHandler;
            MetaMorpheusTask.FinishedWritingFileHandler += NewSuccessfullyFinishedWritingFile;
            MetaMorpheusTask.StartingDataFileHandler += MyTaskEngine_StartingSpectraFileHandler;
            MetaMorpheusTask.FinishedDataFileHandler += MyTaskEngine_FinishedSpectraFileHandler;
            MetaMorpheusTask.OutLabelStatusHandler += NewoutLabelStatus;
            MetaMorpheusTask.NewCollectionHandler += NewCollectionHandler;
            MetaMorpheusTask.OutProgressHandler += NewoutProgressBar;
            MetaMorpheusTask.WarnHandler += NotificationHandler;

            MetaMorpheusEngine.OutProgressHandler += NewoutProgressBar;
            MetaMorpheusEngine.OutLabelStatusHandler += NewoutLabelStatus;
            MetaMorpheusEngine.WarnHandler += NotificationHandler;

            MyFileManager.WarnHandler += NotificationHandler;
            Application.Current.MainWindow.Closing += new CancelEventHandler(MainWindow_Closing);

            KeyDown += new KeyEventHandler(Window_KeyDown);
        }

        private void MyWindow_Loaded(object sender, RoutedEventArgs e)
        {
            UpdateSpectraFileGuiStuff();
            UpdateTaskGuiStuff();
            UpdateOutputFolderTextbox();
            FileSpecificParameters.ValidateFileSpecificVariableNames();
            SearchModifications.SetUpModSearchBoxes();
            PrintErrorsReadingMods();

            if (!UpdateGUISettings.LoadGUISettings())
            {
                notificationsTextBox.Document = WelcomeText();
            }

            if (UpdateGUISettings.Params.AskAboutUpdating)
            {
                UpdateMetaMorpheus();
            }

            // hide the "InProgress" column
            //dataGridProteinDatabases.Columns.Where(p => p.Header.Equals(nameof(ProteinDbForDataGrid.InProgress))).First().Visibility = Visibility.Hidden;
            //dataGridSpectraFiles.Columns.Where(p => p.Header.Equals(nameof(RawDataForDataGrid.InProgress))).First().Visibility = Visibility.Hidden;
        }

        #region Events triggered by MetaMorpheus

        private void EverythingRunnerEngine_FinishedWritingAllResultsFileHandler(object sender, StringEventArgs e)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => EverythingRunnerEngine_FinishedWritingAllResultsFileHandler(sender, e)));
            }
            else
            {
                DynamicTasksObservableCollection.Add(new InRunTask("All Task Results", null));
                DynamicTasksObservableCollection.Last().Progress = 100;
                DynamicTasksObservableCollection.Last().Children.Add(new OutputFileForTreeView(e.S, Path.GetFileNameWithoutExtension(e.S)));
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

        private void MyTaskEngine_FinishedSpectraFileHandler(object sender, StringEventArgs s)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => MyTaskEngine_FinishedSpectraFileHandler(sender, s)));
            }
            else
            {
                var huh = SpectraFilesObservableCollection.First(b => b.FilePath.Equals(s.S));
                huh.SetInProgress(false);

                dataGridSpectraFiles.Items.Refresh();
            }
        }

        private void MyTaskEngine_StartingSpectraFileHandler(object sender, StringEventArgs s)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => MyTaskEngine_StartingSpectraFileHandler(sender, s)));
            }
            else
            {
                var huh = SpectraFilesObservableCollection.First(b => b.FilePath.Equals(s.S));
                huh.SetInProgress(true);
                dataGridSpectraFiles.Items.Refresh();
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
                foreach (var uu in ProteinDbObservableCollection)
                {
                    uu.Use = false;
                }

                foreach (var uu in e.NewDatabases)
                {
                    ProteinDbObservableCollection.Add(new ProteinDbForDataGrid(uu));
                }

                dataGridProteinDatabases.Items.Refresh();
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
                foreach (var oldFile in SpectraFilesObservableCollection)
                {
                    if (!newFiles.Contains(oldFile.FilePath))
                    {
                        oldFile.Use = false;
                    }
                }

                var files = SpectraFilesObservableCollection.Select(p => p.FilePath).ToList();
                foreach (var newRawData in newFiles.Where(p => !files.Contains(p)))
                {
                    SpectraFilesObservableCollection.Add(new RawDataForDataGrid(newRawData));
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
                    UpdateFileSpecificParamsDisplayJustAdded(path);
                }
            }
        }

        private void Po_startingSingleTaskHander(object sender, SingleTaskEventArgs s)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => Po_startingSingleTaskHander(sender, s)));
            }
            else
            {
                var theTask = DynamicTasksObservableCollection.First(b => b.DisplayName.Equals(s.DisplayName));
                theTask.Status = "Starting...";

                dataGridSpectraFiles.Items.Refresh();
                dataGridProteinDatabases.Items.Refresh();
            }
        }

        private void Po_finishedSingleTaskHandler(object sender, SingleTaskEventArgs s)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => Po_finishedSingleTaskHandler(sender, s)));
            }
            else
            {
                var theTask = DynamicTasksObservableCollection.First(b => b.DisplayName.Equals(s.DisplayName));
                theTask.IsIndeterminate = false;
                theTask.Progress = 100;
                theTask.Status = "Done!";

                dataGridSpectraFiles.Items.Refresh();
                dataGridProteinDatabases.Items.Refresh();
            }
        }

        private void NewCollectionHandler(object sender, StringEventArgs s)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => NewCollectionHandler(sender, s)));
            }
            else
            {
                // Find the task or the collection!!!

                ForTreeView theEntityOnWhichToUpdateLabel = DynamicTasksObservableCollection.First(b => b.Id.Equals(s.NestedIDs[0]));

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

                ForTreeView theEntityOnWhichToUpdateLabel = DynamicTasksObservableCollection.First(b => b.Id.Equals(s.NestedIDs[0]));

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

        // update progress bar for task/file
        private void NewoutProgressBar(object sender, ProgressEventArgs s)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => NewoutProgressBar(sender, s)));
            }
            else
            {
                ForTreeView theEntityOnWhichToUpdateLabel = DynamicTasksObservableCollection.First(b => b.Id.Equals(s.NestedIDs[0]));

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

        private void NewRefreshBetweenTasks(object sender, EventArgs e)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => NewRefreshBetweenTasks(sender, e)));
            }
            else
            {
                dataGridSpectraFiles.Items.Refresh();
                dataGridProteinDatabases.Items.Refresh();
            }
        }

        private void NewSuccessfullyStartingAllTasks(object sender, EventArgs e)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => NewSuccessfullyStartingAllTasks(sender, e)));
            }
            else
            {
                dataGridSpectraFiles.Items.Refresh();

                RunTasksButton.IsEnabled = false;
            }
        }

        private void NewSuccessfullyFinishedAllTasks(object sender, StringEventArgs e)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => NewSuccessfullyFinishedAllTasks(sender, e)));
            }
            else
            {
                //ResetTasksButton.IsEnabled = true;

                dataGridSpectraFiles.Items.Refresh();
            }
        }

        private void NewSuccessfullyFinishedWritingFile(object sender, SingleFileEventArgs v)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => NewSuccessfullyFinishedWritingFile(sender, v)));
            }
            else
            {
                ForTreeView AddWrittenFileToThisOne = DynamicTasksObservableCollection.First(b => b.Id.Equals(v.NestedIDs[0]));

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

        #endregion

        #region Events triggered by user interaction

        /// <summary>
        /// Event fires when a file is dragged-and-dropped into MetaMorpheus.
        /// </summary>
        private void Window_Drop(object sender, DragEventArgs e)
        {
            if (RunTasksButton.IsEnabled)
            {
                string[] files = ((string[])e.Data.GetData(DataFormats.FileDrop)).OrderBy(p => p).ToArray();

                if (files != null)
                {
                    foreach (var draggedFilePath in files)
                    {
                        if (Directory.Exists(draggedFilePath))
                        {
                            foreach (string file in Directory.EnumerateFiles(draggedFilePath, "*.*", SearchOption.AllDirectories))
                            {
                                AddAFile(file);
                            }
                        }
                        else
                        {
                            AddAFile(draggedFilePath);
                        }
                        dataGridSpectraFiles.CommitEdit(DataGridEditingUnit.Row, true);
                        dataGridProteinDatabases.CommitEdit(DataGridEditingUnit.Row, true);
                        dataGridSpectraFiles.Items.Refresh();
                        dataGridProteinDatabases.Items.Refresh();
                    }
                }
                UpdateTaskGuiStuff();
            }
        }

        /// <summary>
        /// Event fires when the "Add Spectra" button is clicked.
        /// </summary>
        private void AddSpectraFile_Click(object sender, RoutedEventArgs e)
        {
            var openPicker = OpenFileDialog("Spectra Files(*.raw;*.mzML;*.mgf)|*.raw;*.mzML;*.mgf");

            if (openPicker.ShowDialog() == true)
            {
                foreach (var rawDataFromSelected in openPicker.FileNames.OrderBy(p => p))
                {
                    AddAFile(rawDataFromSelected);
                }
            }

            dataGridSpectraFiles.Items.Refresh();
        }

        private void ClearSpectraFiles_Click(object sender, RoutedEventArgs e)
        {
            SpectraFilesObservableCollection.Clear();
            UpdateOutputFolderTextbox();
        }

        private void ChangeFileParameters_Click(object sender, RoutedEventArgs e)
        {
            try
            {
                var dialog = new FileSpecificParametersWindow(SelectedRawFiles);
                if (dialog.ShowDialog() == true)
                {
                    var tomlPathsForSelectedFiles = SelectedRawFiles.Select(p => Path.Combine(Directory.GetParent(p.FilePath).ToString(), Path.GetFileNameWithoutExtension(p.FileName)) + ".toml").ToList();
                    UpdateFileSpecificParamsDisplay(tomlPathsForSelectedFiles.ToArray());
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
            var dialog = new ExperimentalDesignWindow(SpectraFilesObservableCollection);
            dialog.ShowDialog();
        }

        /// <summary>
        /// Event fires when the "Add Protein Database" button is clicked.
        /// </summary>
        private void AddProteinDatabase_Click(object sender, RoutedEventArgs e)
        {
            var openPicker = OpenFileDialog("Database Files|*.xml;*.xml.gz;*.fasta;*.fa");

            if (openPicker.ShowDialog() == true)
            {
                foreach (var filepath in openPicker.FileNames.OrderBy(p => p))
                {
                    AddAFile(filepath);
                }
            }

            dataGridProteinDatabases.Items.Refresh();
        }

        private void AddDefaultContaminantDatabase_Click(object sender, RoutedEventArgs e)
        {
            string[] contaminantFiles = Directory.GetFiles(Path.Combine(GlobalVariables.DataDir, "Contaminants"));
            foreach (string contaminantFile in contaminantFiles)
            {
                AddAFile(contaminantFile);
            }
            dataGridProteinDatabases.Items.Refresh();
        }

        private void ClearProteinDatabases_Click(object sender, RoutedEventArgs e)
        {
            ProteinDbObservableCollection.Clear();
        }

        /// <summary>
        /// Event fires when the "delete" button is clicked on a protein DB or spectra file.
        /// </summary>
        private void DeleteDatabaseOrSpectra_Click(object sender, RoutedEventArgs e)
        {
            RawDataForDataGrid spectraFile = (sender as Button).DataContext as RawDataForDataGrid;
            if (spectraFile != null)
            {
                SpectraFilesObservableCollection.Remove(spectraFile);
                return;
            }

            ProteinDbForDataGrid proteinDbFile = (sender as Button).DataContext as ProteinDbForDataGrid;
            if (proteinDbFile != null)
            {
                ProteinDbObservableCollection.Remove(proteinDbFile);
                return;
            }
        }

        /// <summary>
        /// Event fires when a data grid row (protein DB or spectra file) is double-clicked.
        /// </summary>
        private void DatabaseOrSpectraFile_DoubleClick(object sender, MouseButtonEventArgs e)
        {
            var dataGridCell = sender as DataGridCell;

            // prevent opening protein DB or spectra files if a run is in progress
            if ((dataGridCell.DataContext is ProteinDbForDataGrid || dataGridCell.DataContext is RawDataForDataGrid) && !RunTasksButton.IsEnabled)
            {
                return;
            }

            // open the file with the default process for this file format
            if (dataGridCell.Content is TextBlock filePath && filePath != null && !string.IsNullOrEmpty(filePath.Text))
            {
                OpenFile(filePath.Text);
            }
        }

        private void AddSearchTaskButton_Click(object sender, RoutedEventArgs e)
        {
            //check if the default toml has been overwritten
            SearchTask task = null;
            string defaultFilePath = Path.Combine(GlobalVariables.DataDir, "DefaultParameters", @"SearchTaskDefault.toml");
            if (File.Exists(defaultFilePath))
            {
                try
                {
                    task = Toml.ReadFile<SearchTask>(defaultFilePath, MetaMorpheusTask.tomlConfig);
                }
                catch (Exception)
                {
                    NotificationHandler(null, new StringEventArgs("Cannot read toml: " + defaultFilePath, null));
                }
            }

            var dialog = new SearchTaskWindow(task);
            if (dialog.ShowDialog() == true)
            {
                AddTaskToCollection(dialog.TheTask);
                UpdateTaskGuiStuff();
            }
        }

        private void AddCalibrateTaskButton_Click(object sender, RoutedEventArgs e)
        {
            //check if the default toml has been overwritten
            CalibrationTask task = null;
            string defaultFilePath = Path.Combine(GlobalVariables.DataDir, "DefaultParameters", @"CalibrationTaskDefault.toml");
            if (File.Exists(defaultFilePath))
            {
                try
                {
                    task = Toml.ReadFile<CalibrationTask>(defaultFilePath, MetaMorpheusTask.tomlConfig);
                }
                catch (Exception)
                {
                    NotificationHandler(null, new StringEventArgs("Cannot read toml: " + defaultFilePath, null));
                }
            }

            var dialog = new CalibrateTaskWindow(task);
            if (dialog.ShowDialog() == true)
            {
                AddTaskToCollection(dialog.TheTask);
                UpdateTaskGuiStuff();
            }
        }

        private void AddGPTMDTaskButton_Click(object sender, RoutedEventArgs e)
        {
            //check if the default toml has been overwritten
            GptmdTask task = null;
            string defaultFilePath = Path.Combine(GlobalVariables.DataDir, "DefaultParameters", @"GptmdTaskDefault.toml");
            if (File.Exists(defaultFilePath))
            {
                try
                {
                    task = Toml.ReadFile<GptmdTask>(defaultFilePath, MetaMorpheusTask.tomlConfig);
                }
                catch (Exception)
                {
                    NotificationHandler(null, new StringEventArgs("Cannot read toml: " + defaultFilePath, null));
                }
            }

            var dialog = new GptmdTaskWindow(task);
            if (dialog.ShowDialog() == true)
            {
                AddTaskToCollection(dialog.TheTask);
                UpdateTaskGuiStuff();
            }
        }

        private void AddCrosslinkTask_Click(object sender, RoutedEventArgs e)
        {
            //check if the default toml has been overwritten
            XLSearchTask task = null;
            string defaultFilePath = Path.Combine(GlobalVariables.DataDir, "DefaultParameters", @"XLSearchTaskDefault.toml");
            if (File.Exists(defaultFilePath))
            {
                try
                {
                    task = Toml.ReadFile<XLSearchTask>(defaultFilePath, MetaMorpheusTask.tomlConfig);
                }
                catch (Exception)
                {
                    NotificationHandler(null, new StringEventArgs("Cannot read toml: " + defaultFilePath, null));
                }
            }

            var dialog = new XLSearchTaskWindow(task);
            if (dialog.ShowDialog() == true)
            {
                AddTaskToCollection(dialog.TheTask);
                UpdateTaskGuiStuff();
            }
        }

        private void AddGlycoSearchTask_Click(object sender, RoutedEventArgs e)
        {
            var dialog = new GlycoSearchTaskWindow();
            if (dialog.ShowDialog() == true)
            {
                AddTaskToCollection(dialog.TheTask);
                UpdateTaskGuiStuff();
            }
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
            var openPicker = OpenFileDialog("TOML files(*.toml)|*.toml");

            if (openPicker.ShowDialog() == true)
            {
                foreach (var tomlFromSelected in openPicker.FileNames.OrderBy(p => p))
                {
                    AddAFile(tomlFromSelected);
                }
            }

            UpdateTaskGuiStuff();
        }

        /// <summary>
        /// Event fires when the "Save as .toml" context menu item is clicked.
        /// Can occur in the task tree view.
        /// </summary>
        private void SaveTask_Click(object sender, RoutedEventArgs e)
        {
            var menuItem = (MenuItem)sender;
            var dataContext = (ContextMenu)menuItem.Parent;
            var treeViewItem = (TreeViewItem)dataContext.PlacementTarget;

            MetaMorpheusTask task;

            if (treeViewItem.Header.GetType() == typeof(PreRunTask))
            {
                task = ((PreRunTask)treeViewItem.Header).metaMorpheusTask;
            }
            else if (treeViewItem.Header.GetType() == typeof(InRunTask))
            {
                task = ((InRunTask)treeViewItem.Header).Task;
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
        /// Deletes the selected task.
        /// </summary>
        private void DeleteTask_Click(object sender, RoutedEventArgs e)
        {
            var selectedTask = (PreRunTask)tasksTreeView.SelectedItem;
            if (selectedTask != null)
            {
                StaticTasksObservableCollection.Remove(selectedTask);
                UpdateTaskGuiStuff();
            }
        }

        private void ClearTasks_Click(object sender, RoutedEventArgs e)
        {
            StaticTasksObservableCollection.Clear();
            UpdateTaskGuiStuff();
        }

        private void ResetTasks_Click(object sender, RoutedEventArgs e)
        {
            RunTasksButton.IsEnabled = true;

            tasksTreeView.DataContext = StaticTasksObservableCollection;
            UpdateSpectraFileGuiStuff();

            var pathOfFirstSpectraFile = Path.GetDirectoryName(SpectraFilesObservableCollection.First().FilePath);
            OutputFolderTextBox.Text = Path.Combine(pathOfFirstSpectraFile, @"$DATETIME");
        }

        private void CancelTasks_Click(object sender, RoutedEventArgs e)
        {
            string grammar = StaticTasksObservableCollection.Count <= 1 ? "this task" : "these tasks";
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
            var selectedTask = (PreRunTask)tasksTreeView.SelectedItem;
            if (selectedTask == null)
            {
                return;
            }

            int indexOfSelected = StaticTasksObservableCollection.IndexOf(selectedTask);
            int indexToMoveTo = indexOfSelected - 1;
            if (moveTaskUp)
            {
                indexToMoveTo = indexOfSelected + 1;
            }

            if (indexToMoveTo >= 0 && indexToMoveTo < StaticTasksObservableCollection.Count)
            {
                var temp = StaticTasksObservableCollection[indexToMoveTo];
                StaticTasksObservableCollection[indexToMoveTo] = selectedTask;
                StaticTasksObservableCollection[indexOfSelected] = temp;

                UpdateTaskGuiStuff();

                var item = tasksTreeView.ItemContainerGenerator.ContainerFromItem(selectedTask);
                ((TreeViewItem)item).IsSelected = true;
            }
        }

        private void RunAllTasks_Click(object sender, RoutedEventArgs e)
        {
            GlobalVariables.StopLoops = false;
            //CancelButton.IsEnabled = true;

            // check for valid tasks/spectra files/protein databases
            if (!StaticTasksObservableCollection.Any())
            {
                NotificationHandler(null, new StringEventArgs("You need to add at least one task!", null));
                return;
            }
            if (!SpectraFilesObservableCollection.Any())
            {
                NotificationHandler(null, new StringEventArgs("You need to add at least one spectra file!", null));
                return;
            }
            if (!ProteinDbObservableCollection.Any())
            {
                NotificationHandler(null, new StringEventArgs("You need to add at least one protein database!", null));
                return;
            }

            DynamicTasksObservableCollection = new ObservableCollection<InRunTask>();

            for (int i = 0; i < StaticTasksObservableCollection.Count; i++)
            {
                DynamicTasksObservableCollection.Add(new InRunTask("Task" + (i + 1) + "-" + StaticTasksObservableCollection[i].metaMorpheusTask.CommonParameters.TaskDescriptor, StaticTasksObservableCollection[i].metaMorpheusTask));
            }
            tasksTreeView.DataContext = DynamicTasksObservableCollection;

            notificationsTextBox.Document.Blocks.Clear();

            // output folder
            if (string.IsNullOrEmpty(OutputFolderTextBox.Text))
            {
                var pathOfFirstSpectraFile = Path.GetDirectoryName(SpectraFilesObservableCollection.First().FilePath);
                OutputFolderTextBox.Text = Path.Combine(pathOfFirstSpectraFile, @"$DATETIME");
            }

            var startTimeForAllFilenames = DateTime.Now.ToString("yyyy-MM-dd-HH-mm-ss", CultureInfo.InvariantCulture);
            string outputFolder = OutputFolderTextBox.Text.Replace("$DATETIME", startTimeForAllFilenames);
            OutputFolderTextBox.Text = outputFolder;

            // check that experimental design is defined if normalization is enabled
            // TODO: move all of this over to EverythingRunnerEngine
            var searchTasks = StaticTasksObservableCollection
                .Where(p => p.metaMorpheusTask.TaskType == MyTask.Search)
                .Select(p => (SearchTask)p.metaMorpheusTask);

            string pathToExperDesign = Directory.GetParent(SpectraFilesObservableCollection.First().FilePath).FullName;
            pathToExperDesign = Path.Combine(pathToExperDesign, GlobalVariables.ExperimentalDesignFileName);

            foreach (var searchTask in searchTasks.Where(p => p.SearchParameters.Normalize))
            {
                if (!File.Exists(pathToExperDesign))
                {
                    MessageBox.Show("Experimental design must be defined for normalization!\n" +
                        "Click the \"Experimental Design\" button in the bottom left by the spectra files");
                    return;
                }

                // check that experimental design is OK (spectra files may have been added after exper design was defined)
                // TODO: experimental design might still have flaws if user edited the file manually, need to check for this
                var experDesign = File.ReadAllLines(pathToExperDesign).ToDictionary(p => p.Split('\t')[0], p => p);
                var filesToUse = new HashSet<string>(SpectraFilesObservableCollection.Select(p => Path.GetFileNameWithoutExtension(p.FileName)));
                var experDesignFilesDefined = new HashSet<string>(experDesign.Keys);

                var undefined = filesToUse.Except(experDesignFilesDefined);

                if (undefined.Any())
                {
                    MessageBox.Show("Need to define experimental design parameters for file: " + undefined.First());
                    return;
                }
            }
            //BtnQuantSet.IsEnabled = false;

            // everything is OK to run
            EverythingRunnerEngine a = new EverythingRunnerEngine(DynamicTasksObservableCollection.Select(b => (b.DisplayName, b.Task)).ToList(),
                SpectraFilesObservableCollection.Where(b => b.Use).Select(b => b.FilePath).ToList(),
                ProteinDbObservableCollection.Where(b => b.Use).Select(b => new DbForTask(b.FilePath, b.Contaminant)).ToList(),
                outputFolder);

            var t = new Task(a.Run);
            t.ContinueWith(EverythingRunnerExceptionHandler, TaskContinuationOptions.OnlyOnFaulted);
            t.Start();
        }

        /// <summary>
        /// Event fires when an item in the task treeview is right-clicked.
        /// Can occur on a task or written file.
        /// </summary>
        private void TreeViewItem_RightClick(object sender, MouseButtonEventArgs e)
        {
            var treeViewItem = (TreeViewItem)sender;
            var header = treeViewItem.Header.GetType();
            string contextMenuName;

            if (header == typeof(PreRunTask) || header == typeof(InRunTask))
            {
                contextMenuName = "TaskContextMenu";
            }
            else if (header == typeof(OutputFileForTreeView))
            {
                contextMenuName = "WrittenFileContextMenu";
            }
            else
            {
                return;
            }

            ContextMenu contextMenu = FindResource(contextMenuName) as ContextMenu;
            contextMenu.PlacementTarget = sender as TreeViewItem;
            contextMenu.IsOpen = true;
        }

        /// <summary>
        /// Event fires when the "Open containing item" context menu item is clicked.
        /// Can occur on a protein DB, spectra file, or written file.
        /// </summary>
        private void OpenContainingFolder_Click(object sender, RoutedEventArgs e)
        {
            var menuItem = (MenuItem)sender;
            var dataContext = (ContextMenu)menuItem.Parent;
            var treeViewItem = (TreeViewItem)dataContext.PlacementTarget;

            if (treeViewItem.Header is OutputFileForTreeView writtenFile)
            {
                OpenFolder(Path.GetDirectoryName(writtenFile.FullPath));
            }
        }

        /// <summary>
        /// Event fires when the "Open file" context menu item is clicked.
        /// Can occur on a protein DB, spectra file, or written file.
        /// </summary>
        private void OpenFile_Click(object sender, RoutedEventArgs e)
        {
            var menuItem = (MenuItem)sender;
            var dataContext = (ContextMenu)menuItem.Parent;
            var treeViewItem = (TreeViewItem)dataContext.PlacementTarget;

            if (treeViewItem.Header is OutputFileForTreeView writtenFile)
            {
                OpenFile(writtenFile.FullPath);
            }
        }

        private void TasksTreeView_MouseDoubleClick(object sender, EventArgs e)
        {
            var a = sender as TreeView;
            if (a.SelectedItem is PreRunTask preRunTask)
            {
                switch (preRunTask.metaMorpheusTask.TaskType)
                {
                    case MyTask.Search:

                        var searchDialog = new SearchTaskWindow(preRunTask.metaMorpheusTask as SearchTask);
                        searchDialog.ShowDialog();
                        preRunTask.DisplayName = "Task" + (StaticTasksObservableCollection.IndexOf(preRunTask) + 1) + "-" + searchDialog.TheTask.CommonParameters.TaskDescriptor;
                        tasksTreeView.Items.Refresh();

                        return;

                    case MyTask.Gptmd:
                        var gptmddialog = new GptmdTaskWindow(preRunTask.metaMorpheusTask as GptmdTask);
                        gptmddialog.ShowDialog();
                        preRunTask.DisplayName = "Task" + (StaticTasksObservableCollection.IndexOf(preRunTask) + 1) + "-" + gptmddialog.TheTask.CommonParameters.TaskDescriptor;
                        tasksTreeView.Items.Refresh();

                        return;

                    case MyTask.Calibrate:
                        var calibratedialog = new CalibrateTaskWindow(preRunTask.metaMorpheusTask as CalibrationTask);
                        calibratedialog.ShowDialog();
                        preRunTask.DisplayName = "Task" + (StaticTasksObservableCollection.IndexOf(preRunTask) + 1) + "-" + calibratedialog.TheTask.CommonParameters.TaskDescriptor;
                        tasksTreeView.Items.Refresh();
                        return;

                    case MyTask.XLSearch:
                        var XLSearchdialog = new XLSearchTaskWindow(preRunTask.metaMorpheusTask as XLSearchTask);
                        XLSearchdialog.ShowDialog();
                        preRunTask.DisplayName = "Task" + (StaticTasksObservableCollection.IndexOf(preRunTask) + 1) + "-" + XLSearchdialog.TheTask.CommonParameters.TaskDescriptor;
                        tasksTreeView.Items.Refresh();
                        return;

                    case MyTask.GlycoSearch:
                        var GlycoSearchdialog = new GlycoSearchTaskWindow(preRunTask.metaMorpheusTask as GlycoSearchTask);
                        GlycoSearchdialog.ShowDialog();
                        preRunTask.DisplayName = "Task" + (StaticTasksObservableCollection.IndexOf(preRunTask) + 1) + "-" + GlycoSearchdialog.TheTask.CommonParameters.TaskDescriptor;
                        tasksTreeView.Items.Refresh();
                        return;
                }
            }

            if (a.SelectedItem is OutputFileForTreeView writtenFile)
            {
                OpenFile(writtenFile.FullPath);
            }
        }

        /// <summary>
        /// Event fires when the "Open" button is clicked (referring to the output folder).
        /// </summary>
        private void OpenOutputFolder_Click(object sender, RoutedEventArgs e)
        {
            string outputFolder = OutputFolderTextBox.Text;
            if (outputFolder.Contains("$DATETIME"))
            {
                // the exact file path isn't known, so just open the parent directory
                outputFolder = Directory.GetParent(outputFolder).FullName;
            }

            if (!Directory.Exists(outputFolder) && !string.IsNullOrEmpty(outputFolder))
            {
                // create the directory if it doesn't exist yet
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
        /// </summary>
        private void Window_KeyDown(object sender, KeyEventArgs e)
        {
            if (RunTasksButton.IsEnabled)
            {
                // delete selected task
                if (e.Key == Key.Delete || e.Key == Key.Back)
                {
                    DeleteTask_Click(sender, e);
                    e.Handled = true;
                }

                // move task up
                if (e.Key == Key.Add || e.Key == Key.OemPlus)
                {
                    MoveSelectedTask_Click(sender, e, true);
                    e.Handled = true;
                }

                // move task down
                if (e.Key == Key.Subtract || e.Key == Key.OemMinus)
                {
                    MoveSelectedTask_Click(sender, e, false);
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
                //e.cancel is if the event should be canceled (where the event is closing MetaMorpheus). Example: if(e.cancel){keep MetaMorpheus open}
                var exit = ExitMsgBox.Show("Exit MetaMorpheus", "Are you sure you want to exit MetaMorpheus?", "Yes", "No", "Yes and don't ask me again");

                if (exit == MessageBoxResult.Yes)
                {
                    e.Cancel = false;
                }
                else if (exit == MessageBoxResult.OK) //don't ask me again
                {
                    UpdateGUISettings.Params.AskBeforeExitingMetaMorpheus = false;
                    Toml.WriteFile(UpdateGUISettings.Params, Path.Combine(GlobalVariables.DataDir, @"GUIsettings.toml"), MetaMorpheusTask.tomlConfig);
                    e.Cancel = false;
                }
                else //assume the event should be canceled and MetaMorpheus should stay open
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

        private void MenuItem_Wiki_Click(object sender, RoutedEventArgs e)
        {
            GlobalVariables.StartProcess(@"https://github.com/smith-chem-wisc/MetaMorpheus/wiki");
        }

        private void MenuItem_YouTube_Click(object sender, RoutedEventArgs e)
        {
            GlobalVariables.StartProcess(@"https://www.youtube.com/playlist?list=PLVk5tTSZ1aWlhNPh7jxPQ8pc0ElyzSUQb");
        }

        private void MenuItem_ProteomicsNewsBlog_Click(object sender, RoutedEventArgs e)
        {
            GlobalVariables.StartProcess(@"https://proteomicsnews.blogspot.com/");
        }

        private void MenuItem_UpdateMetaMorpheus_Click(object sender, RoutedEventArgs e)
        {
            UpdateMetaMorpheus(printMessageIfThisIsLatestVersion: true);
        }

        private void MenuItem_EmailHelp_Click(object sender, RoutedEventArgs e)
        {
            string mailto = string.Format("mailto:{0}?Subject=MetaMorpheus. Issue:", "mm_support@chem.wisc.edu");
            GlobalVariables.StartProcess(mailto);
        }

        private void MenuItem_GitHubIssues_Click(object sender, RoutedEventArgs e)
        {
            GlobalVariables.StartProcess(@"https://github.com/smith-chem-wisc/MetaMorpheus/issues/new");
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

        private void MenuItem_OpenDataDir_Click(object sender, RoutedEventArgs e)
        {
            GlobalVariables.StartProcess(GlobalVariables.DataDir);
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

        private FlowDocument WelcomeText()
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

        private void UpdateTaskGuiStuff()
        {
            if (StaticTasksObservableCollection.Count == 0)
            {
                RunTasksButton.IsEnabled = false;
            }
            else
            {
                RunTasksButton.IsEnabled = true;

                // this exists so that when a task is deleted, the remaining tasks are renamed to keep the task numbers correct
                for (int i = 0; i < StaticTasksObservableCollection.Count; i++)
                {
                    string newName = "Task" + (i + 1) + "-" + StaticTasksObservableCollection[i].metaMorpheusTask.CommonParameters.TaskDescriptor;
                    StaticTasksObservableCollection[i].DisplayName = newName;
                }
                tasksTreeView.Items.Refresh();
            }
        }

        private void UpdateSpectraFileGuiStuff()
        {
            //ChangeFileParameters.IsEnabled = SelectedRawFiles.Count > 0 && LoadTaskButton.IsEnabled;
        }

        //run if fileSpecificParams are changed from GUI
        private void UpdateFileSpecificParamsDisplay(string[] tomlLocations)
        {
            string[] fullPathofTomls = tomlLocations;

            foreach (var file in SelectedRawFiles)
            {
                for (int j = 0; j < fullPathofTomls.Count(); j++)
                {
                    if (Path.GetFileNameWithoutExtension(file.FileName) == Path.GetFileNameWithoutExtension(fullPathofTomls[j]))
                    {
                        if (File.Exists(fullPathofTomls[j]))
                        {
                            TomlTable fileSpecificSettings = Toml.ReadFile(fullPathofTomls[j], MetaMorpheusTask.tomlConfig);
                            try
                            {
                                // parse to make sure toml is readable
                                var temp = new FileSpecificParameters(fileSpecificSettings);

                                // toml is ok; display the file-specific settings in the gui
                                file.SetParametersText(File.ReadAllText(fullPathofTomls[j]));
                            }
                            catch (MetaMorpheusException e)
                            {
                                NotificationHandler(null, new StringEventArgs("Problem parsing the file-specific toml " + Path.GetFileName(fullPathofTomls[j]) + "; " + e.Message + "; is the toml from an older version of MetaMorpheus?", null));
                            }
                        }
                        else
                        {
                            file.SetParametersText(null);
                        }
                    }
                }
            }
            UpdateSpectraFileGuiStuff();
            dataGridSpectraFiles.Items.Refresh();
        }

        //run if data file has just been added with and checks for Existing fileSpecficParams
        private void UpdateFileSpecificParamsDisplayJustAdded(string tomlLocation)
        {
            string assumedPathOfSpectraFileWithoutExtension = Path.Combine(Directory.GetParent(tomlLocation).ToString(), Path.GetFileNameWithoutExtension(tomlLocation));

            for (int i = 0; i < SpectraFilesObservableCollection.Count; i++)
            {
                string thisFilesPathWihoutExtension = Path.Combine(Directory.GetParent(SpectraFilesObservableCollection[i].FilePath).ToString(), Path.GetFileNameWithoutExtension(SpectraFilesObservableCollection[i].FilePath));
                if (File.Exists(tomlLocation) && assumedPathOfSpectraFileWithoutExtension.Equals(thisFilesPathWihoutExtension))
                {
                    TomlTable fileSpecificSettings = Toml.ReadFile(tomlLocation, MetaMorpheusTask.tomlConfig);
                    try
                    {
                        // parse to make sure toml is readable
                        var temp = new FileSpecificParameters(fileSpecificSettings);

                        // toml is ok; display the file-specific settings in the gui
                        SpectraFilesObservableCollection[i].SetParametersText(File.ReadAllText(tomlLocation));
                    }
                    catch (MetaMorpheusException e)
                    {
                        NotificationHandler(null, new StringEventArgs("Problem parsing the file-specific toml " + Path.GetFileName(tomlLocation) + "; " + e.Message + "; is the toml from an older version of MetaMorpheus?", null));
                    }
                    catch (KeyNotFoundException e)
                    {
                        NotificationHandler(null, new StringEventArgs("Problem parsing the file-specific toml " + Path.GetFileName(tomlLocation) + "; " + e.Message + "; please update the proteases.tsv file and restart MetaMorpheus to use this file-specific toml.", null));
                    }
                }
            }
            UpdateSpectraFileGuiStuff();
            dataGridSpectraFiles.Items.Refresh();
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
                //ResetTasksButton.IsEnabled = true;
            }
        }

        private void UpdateOutputFolderTextbox()
        {
            if (SpectraFilesObservableCollection.Any())
            {
                // if current output folder is blank and there is a spectra file, use the spectra file's path as the output path
                if (string.IsNullOrWhiteSpace(OutputFolderTextBox.Text))
                {
                    var pathOfFirstSpectraFile = Path.GetDirectoryName(SpectraFilesObservableCollection.First().FilePath);
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

        private void AddAFile(string filePath)
        {
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
                        thermoLicenceWindow.LicenceText.AppendText(ThermoRawFileReader.ThermoRawFileReaderLicence.ThermoLicenceText);
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
                    if (!SpectraFileExists(SpectraFilesObservableCollection, zz))
                    {
                        SpectraFilesObservableCollection.Add(zz);
                    }
                    UpdateFileSpecificParamsDisplayJustAdded(Path.ChangeExtension(filePath, ".toml"));
                    UpdateOutputFolderTextbox();
                    break;

                case ".xml":
                case ".fasta":
                case ".fa":
                    ProteinDbForDataGrid uu = new ProteinDbForDataGrid(filePath);
                    if (!DatabaseExists(ProteinDbObservableCollection, uu))
                    {
                        ProteinDbObservableCollection.Add(uu);
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
                                ProteinDbObservableCollection.Remove(uu);
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
                                    var ye1 = Toml.ReadFile<SearchTask>(filePath, MetaMorpheusTask.tomlConfig);
                                    AddTaskToCollection(ye1);
                                    break;

                                case "Calibrate":
                                    var ye2 = Toml.ReadFile<CalibrationTask>(filePath, MetaMorpheusTask.tomlConfig);
                                    AddTaskToCollection(ye2);
                                    break;

                                case "Gptmd":
                                    var ye3 = Toml.ReadFile<GptmdTask>(filePath, MetaMorpheusTask.tomlConfig);
                                    AddTaskToCollection(ye3);
                                    break;

                                case "XLSearch":
                                    var ye4 = Toml.ReadFile<XLSearchTask>(filePath, MetaMorpheusTask.tomlConfig);
                                    AddTaskToCollection(ye4);
                                    break;

                                case "GlycoSearch":
                                    var ye5 = Toml.ReadFile<GlycoSearchTask>(filePath, MetaMorpheusTask.tomlConfig);
                                    AddTaskToCollection(ye5);
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

        private void AddTaskToCollection(MetaMorpheusTask ye)
        {
            PreRunTask te = new PreRunTask(ye);
            StaticTasksObservableCollection.Add(te);
            StaticTasksObservableCollection.Last().DisplayName = "Task" + (StaticTasksObservableCollection.IndexOf(te) + 1) + "-" + ye.CommonParameters.TaskDescriptor;
        }

        private OpenFileDialog OpenFileDialog(string filter)
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

        #endregion
    }
}