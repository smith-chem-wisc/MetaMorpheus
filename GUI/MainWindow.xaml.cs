using EngineLayer;
using Nett;
using System;
using System.Collections.ObjectModel;
using System.IO;
using System.Linq;
using System.Threading;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;
using TaskLayer;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {

        #region Private Fields

        private readonly ObservableCollection<RawDataForDataGrid> rawDataObservableCollection = new ObservableCollection<RawDataForDataGrid>();
        private readonly ObservableCollection<ProteinDbForDataGrid> proteinDbObservableCollection = new ObservableCollection<ProteinDbForDataGrid>();
        private readonly ObservableCollection<FinishedFileForDataGrid> finishedFileObservableCollection = new ObservableCollection<FinishedFileForDataGrid>();
        private readonly ObservableCollection<PreRunTask> staticTasksObservableCollection = new ObservableCollection<PreRunTask>();
        private ObservableCollection<InRunTask> dynamicTasksObservableCollection;

        #endregion Private Fields

        #region Public Constructors

        public MainWindow()
        {
            InitializeComponent();

            Title = MyEngine.MetaMorpheusVersion.Equals("1.0.0.0") ?
                "MetaMorpheus: Not a release version" :
                "MetaMorpheus: version " + MyEngine.MetaMorpheusVersion;

            dataGridXMLs.DataContext = proteinDbObservableCollection;

            dataGridDatafiles.DataContext = rawDataObservableCollection;
            tasksTreeView.DataContext = staticTasksObservableCollection;

            EverythingRunnerEngine.newDbsHandler += AddNewDB;
            EverythingRunnerEngine.newSpectrasHandler += AddNewSpectra;

            EverythingRunnerEngine.startingAllTasksEngineHandler += NewSuccessfullyStartingAllTasks;
            EverythingRunnerEngine.finishedAllTasksEngineHandler += NewSuccessfullyFinishedAllTasks;

            foreach (var modFile in Directory.GetFiles(@"Mods"))
            {
                MetaMorpheusTask.AddModList(modFile);
            }

            MetaMorpheusTask.StartingSingleTaskHander += Po_startingSingleTaskHander;
            MetaMorpheusTask.FinishedSingleTaskHandler += Po_finishedSingleTaskHandler;
            MetaMorpheusTask.FinishedWritingFileHandler += NewSuccessfullyFinishedWritingFile;
            MetaMorpheusTask.StartingDataFileHandler += MyTaskEngine_StartingDataFileHandler;
            MetaMorpheusTask.FinishedDataFileHandler += MyTaskEngine_FinishedDataFileHandler;
            MetaMorpheusTask.OutLabelStatusHandler += NewoutLabelStatus;
            MetaMorpheusTask.NewCollectionHandler += NewCollectionHandler;

            MyEngine.OutProgressHandler += NewoutProgressBar;
            MyEngine.OutLabelStatusHandler += NewoutLabelStatus;

            UpdateTaskGuiStuff();
        }

        #endregion Public Constructors

        #region Private Methods

        private void MyTaskEngine_FinishedDataFileHandler(object sender, StringEventArgs s)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => MyTaskEngine_FinishedDataFileHandler(sender, s)));
            }
            else
            {
                var huh = rawDataObservableCollection.First(b => b.FileName.Equals(s.s));
                huh.SetInProgress(false);
                dataGridDatafiles.Items.Refresh();
            }
        }

        private void MyTaskEngine_StartingDataFileHandler(object sender, StringEventArgs s)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => MyTaskEngine_StartingDataFileHandler(sender, s)));
            }
            else
            {
                var huh = rawDataObservableCollection.First(b => b.FileName.Equals(s.s));
                huh.SetInProgress(true);
                dataGridDatafiles.Items.Refresh();
            }
        }

        private void AddNewDB(object sender, XmlForTaskListEventArgs e)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => AddNewDB(sender, e)));
            }
            else
            {
                foreach (var uu in proteinDbObservableCollection)
                    uu.Use = false;
                foreach (var uu in e.newDatabases)
                    proteinDbObservableCollection.Add(new ProteinDbForDataGrid(uu));
            }
        }

        private void AddNewSpectra(object sender, StringListEventArgs e)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => AddNewSpectra(sender, e)));
            }
            else
            {
                foreach (var uu in rawDataObservableCollection)
                    uu.Use = false;
                foreach (var newRawData in e.StringList)
                    rawDataObservableCollection.Add(new RawDataForDataGrid(newRawData));
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
                var theTask = dynamicTasksObservableCollection.First(b => b.Id.Equals(s.TaskId));
                theTask.InProgress = true;
                theTask.IsIndeterminate = true;
                theTask.Status = "Starting...";

                dataGridDatafiles.Items.Refresh();
                dataGridXMLs.Items.Refresh();
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
                var theTask = dynamicTasksObservableCollection.First(b => b.Id.Equals(s.TaskId));
                theTask.InProgress = false;
                theTask.IsIndeterminate = false;
                theTask.Progress = 100;
                theTask.Status = "Done!";

                dataGridDatafiles.Items.Refresh();
                dataGridXMLs.Items.Refresh();
            }
        }

        private void AddFinishedFile(string filepath)
        {
            finishedFileObservableCollection.Add(new FinishedFileForDataGrid(filepath));
        }

        private void ClearRaw_Click(object sender, RoutedEventArgs e)
        {
            rawDataObservableCollection.Clear();
        }

        private void AddXML_Click(object sender, RoutedEventArgs e)
        {
            Microsoft.Win32.OpenFileDialog openPicker = new Microsoft.Win32.OpenFileDialog()
            {
                Filter = "Database Files|*.xml;*.xml.gz;*.fasta",
                FilterIndex = 1,
                RestoreDirectory = true,
                Multiselect = true
            };
            if (openPicker.ShowDialog() == true)
                foreach (var filepath in openPicker.FileNames)
                    proteinDbObservableCollection.Add(new ProteinDbForDataGrid(filepath));
            dataGridXMLs.Items.Refresh();
        }

        private void AddRaw_Click(object sender, RoutedEventArgs e)
        {
            Microsoft.Win32.OpenFileDialog openFileDialog1 = new Microsoft.Win32.OpenFileDialog()
            {
                Filter = "Spectra Files(*.raw;*.mzML)|*.raw;*.mzML",
                FilterIndex = 1,
                RestoreDirectory = true,
                Multiselect = true
            };
            if (openFileDialog1.ShowDialog() == true)
                foreach (var rawDataFromSelected in openFileDialog1.FileNames)
                    rawDataObservableCollection.Add(new RawDataForDataGrid(rawDataFromSelected));
            dataGridDatafiles.Items.Refresh();
        }

        private void Window_Drop(object sender, DragEventArgs e)
        {
            string[] files = (string[])e.Data.GetData(DataFormats.FileDrop);
            if (files != null)
                foreach (var draggedFilePath in files)
                {
                    var theExtension = Path.GetExtension(draggedFilePath).ToLowerInvariant();
                    switch (theExtension)
                    {
                        case ".raw":
                        case ".mzml":
                            rawDataObservableCollection.Add(new RawDataForDataGrid(draggedFilePath));
                            break;

                        case ".xml":
                        case ".fasta":
                        case ".gz":
                            proteinDbObservableCollection.Add(new ProteinDbForDataGrid(draggedFilePath));
                            break;

                        case ".toml":
                            var uhum = Toml.ReadFile(draggedFilePath, MetaMorpheusTask.tomlConfig);
                            switch (uhum.Get<string>("TaskType"))
                            {
                                case "Search":
                                    var ye1 = Toml.ReadFile<SearchTask>(draggedFilePath, MetaMorpheusTask.tomlConfig);
                                    staticTasksObservableCollection.Add(new PreRunTask(ye1));
                                    break;

                                case "Calibrate":
                                    var ye2 = Toml.ReadFile<CalibrationTask>(draggedFilePath, MetaMorpheusTask.tomlConfig);
                                    staticTasksObservableCollection.Add(new PreRunTask(ye2));
                                    break;

                                case "Gptmd":
                                    var ye3 = Toml.ReadFile<GptmdTask>(draggedFilePath, MetaMorpheusTask.tomlConfig);
                                    staticTasksObservableCollection.Add(new PreRunTask(ye3));
                                    break;
                            }
                            break;
                    }
                    dataGridDatafiles.Items.Refresh();
                    dataGridXMLs.Items.Refresh();
                }
            UpdateTaskGuiStuff();
        }

        private void Row_DoubleClick(object sender, MouseButtonEventArgs e)
        {
            var ye = sender as DataGridCell;
            if (ye.Content is TextBlock hm && !string.IsNullOrEmpty(hm.Text))
            {
                System.Diagnostics.Process.Start(hm.Text);
            }
        }

        private void RunAllTasks_Click(object sender, RoutedEventArgs e)
        {
            dynamicTasksObservableCollection = new ObservableCollection<InRunTask>();

            for (int i = 0; i < staticTasksObservableCollection.Count; i++)
                dynamicTasksObservableCollection.Add(new InRunTask("Task" + (i + 1) + staticTasksObservableCollection[i].metaMorpheusTask.TaskType, staticTasksObservableCollection[i].metaMorpheusTask));
            tasksTreeView.DataContext = dynamicTasksObservableCollection;

            EverythingRunnerEngine a = new EverythingRunnerEngine(dynamicTasksObservableCollection.Select(b => new Tuple<string, MetaMorpheusTask>(b.Id, b.task)).ToList(), rawDataObservableCollection.Where(b => b.Use).Select(b => b.FileName).ToList(), proteinDbObservableCollection.Where(b => b.Use).Select(b => new DbForTask(b.FileName, b.Contaminant)).ToList());
            var t = new Thread(() => a.Run())
            {
                IsBackground = true
            };
            t.Start();
        }

        private void ClearTasks_Click(object sender, RoutedEventArgs e)
        {
            staticTasksObservableCollection.Clear();
            UpdateTaskGuiStuff();
        }

        private void UpdateTaskGuiStuff()
        {
            if (staticTasksObservableCollection.Count == 0)
            {
                RunTasksButton.IsEnabled = false;
                RemoveLastTaskButton.IsEnabled = false;
                ClearTasksButton.IsEnabled = false;
                ResetTasksButton.IsEnabled = false;
            }
            else
            {
                RunTasksButton.IsEnabled = true;
                RemoveLastTaskButton.IsEnabled = true;
                ClearTasksButton.IsEnabled = true;
            }
        }

        private void AddSearchTaskButton_Click(object sender, RoutedEventArgs e)
        {
            var dialog = new SearchTaskWindow();
            if (dialog.ShowDialog() == true)
            {
                staticTasksObservableCollection.Add(new PreRunTask(dialog.TheTask));
                UpdateTaskGuiStuff();
            }
        }

        private void AddCalibrateTaskButton_Click(object sender, RoutedEventArgs e)
        {
            var dialog = new CalibrateTaskWindow();
            if (dialog.ShowDialog() == true)
            {
                staticTasksObservableCollection.Add(new PreRunTask(dialog.TheTask));
                UpdateTaskGuiStuff();
            }
        }

        private void AddGPTMDTaskButton_Click(object sender, RoutedEventArgs e)
        {
            var dialog = new GptmdTaskWindow();
            if (dialog.ShowDialog() == true)
            {
                staticTasksObservableCollection.Add(new PreRunTask(dialog.TheTask));
                UpdateTaskGuiStuff();
            }
        }

        private void RemoveLastTask_Click(object sender, RoutedEventArgs e)
        {
            staticTasksObservableCollection.RemoveAt(staticTasksObservableCollection.Count - 1);
            UpdateTaskGuiStuff();
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

                ForTreeView theEntityOnWhichToUpdateLabel = dynamicTasksObservableCollection.First(b => b.Id.Equals(s.nestedIDs[0]));

                for (int i = 1; i < s.nestedIDs.Count - 1; i++)
                {
                    var hm = s.nestedIDs[i];
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

                theEntityOnWhichToUpdateLabel.Children.Add(new CollectionForTreeView(s.nestedIDs.Last(), s.s));
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

                ForTreeView theEntityOnWhichToUpdateLabel = dynamicTasksObservableCollection.First(b => b.Id.Equals(s.nestedIDs[0]));

                foreach (var hm in s.nestedIDs.Skip(1))
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

                theEntityOnWhichToUpdateLabel.Status = s.s;
                theEntityOnWhichToUpdateLabel.IsIndeterminate = true;
            }
        }

        private void NewoutProgressBar(object sender, ProgressEventArgs s)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => NewoutProgressBar(sender, s)));
            }
            else
            {
                // Find the task or the collection!!!

                ForTreeView theEntityOnWhichToUpdateLabel = dynamicTasksObservableCollection.First(b => b.Id.Equals(s.nestedIDs[0]));

                foreach (var hm in s.nestedIDs.Skip(1))
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

                theEntityOnWhichToUpdateLabel.Status = s.v;
                theEntityOnWhichToUpdateLabel.IsIndeterminate = false;
                theEntityOnWhichToUpdateLabel.Progress = s.new_progress;
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
                dataGridDatafiles.Items.Refresh();
                dataGridXMLs.Items.Refresh();
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
                //statusLabel.Content = "Starting all tasks...";
                //outProgressBar.IsIndeterminate = true;

                dataGridDatafiles.Items.Refresh();

                ClearTasksButton.IsEnabled = false;
                RemoveLastTaskButton.IsEnabled = false;
                RunTasksButton.IsEnabled = false;
                LoadTaskButton.IsEnabled = false;

                addCalibrateTaskButton.IsEnabled = false;
                addGPTMDTaskButton.IsEnabled = false;
                addSearchTaskButton.IsEnabled = false;

                AddXML.IsEnabled = false;
                ClearXML.IsEnabled = false;
                AddRaw.IsEnabled = false;
                ClearRaw.IsEnabled = false;

                proteinDatabasesGroupBox.IsEnabled = false;
                datafilesGroupBox.IsEnabled = false;
            }
        }

        private void NewSuccessfullyFinishedAllTasks(object sender, EventArgs e)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => NewSuccessfullyFinishedAllTasks(sender, e)));
            }
            else
            {
                ResetTasksButton.IsEnabled = true;

                dataGridDatafiles.Items.Refresh();
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
                ForTreeView theEntityOnWhichToUpdateLabel = dynamicTasksObservableCollection.First(b => b.Id.Equals(v.nestedIDs[0]));

                foreach (var hm in v.nestedIDs.Skip(1))
                {
                    try
                    {
                        theEntityOnWhichToUpdateLabel = theEntityOnWhichToUpdateLabel.Children.First(b => b.Id.Equals(hm));
                    }
                    catch { }
                }
                theEntityOnWhichToUpdateLabel.Children.Add(new OutputFileForTreeView(v.writtenFile));
            }
        }

        private void ClearXML_Click(object sender, RoutedEventArgs e)
        {
            proteinDbObservableCollection.Clear();
        }

        private void ClearOutput_Click(object sender, RoutedEventArgs e)
        {
            finishedFileObservableCollection.Clear();
        }

        private void ResetTasksButton_Click(object sender, RoutedEventArgs e)
        {
            tasksGroupBox.IsEnabled = true;
            ClearTasksButton.IsEnabled = true;
            RemoveLastTaskButton.IsEnabled = true;
            RunTasksButton.IsEnabled = true;
            addCalibrateTaskButton.IsEnabled = true;
            addGPTMDTaskButton.IsEnabled = true;
            addSearchTaskButton.IsEnabled = true;
            ResetTasksButton.IsEnabled = false;

            proteinDatabasesGroupBox.IsEnabled = true;
            datafilesGroupBox.IsEnabled = true;

            AddXML.IsEnabled = true;
            ClearXML.IsEnabled = true;
            AddRaw.IsEnabled = true;
            ClearRaw.IsEnabled = true;

            LoadTaskButton.IsEnabled = true;

            tasksTreeView.DataContext = staticTasksObservableCollection;
        }

        private void TasksTreeView_MouseDoubleClick(object sender, MouseButtonEventArgs e)
        {
            var a = sender as TreeView;
            if (a.SelectedItem is PreRunTask preRunTask)
                switch (preRunTask.metaMorpheusTask.TaskType)
                {
                    case MyTask.Search:
                        var searchDialog = new SearchTaskWindow(preRunTask.metaMorpheusTask as SearchTask);
                        searchDialog.ShowDialog();
                        return;

                    case MyTask.Gptmd:
                        var gptmddialog = new GptmdTaskWindow(preRunTask.metaMorpheusTask as GptmdTask);
                        gptmddialog.ShowDialog();
                        return;

                    case MyTask.Calibrate:
                        var calibratedialog = new CalibrateTaskWindow(preRunTask.metaMorpheusTask as CalibrationTask);
                        calibratedialog.ShowDialog();
                        return;
                }

            if (a.SelectedItem is InRunTask inRunTask)
            {
                // Display the params (always) and the results (if done)
            }

            if (a.SelectedItem is OutputFileForTreeView fileThing)
            {
                System.Diagnostics.Process.Start(fileThing.Id);
            }
        }

        private void LoadTaskButton_Click(object sender, RoutedEventArgs e)
        {
            Microsoft.Win32.OpenFileDialog openFileDialog1 = new Microsoft.Win32.OpenFileDialog()
            {
                Filter = "TOML files(*.toml)|*.toml",
                FilterIndex = 1,
                RestoreDirectory = true,
                Multiselect = true
            };
            if (openFileDialog1.ShowDialog() == true)
                foreach (var tomlFromSelected in openFileDialog1.FileNames)
                {
                    var uhum = Toml.ReadFile(tomlFromSelected, MetaMorpheusTask.tomlConfig);
                    switch (uhum.Get<string>("TaskType"))
                    {
                        case "Search":
                            var ye1 = Toml.ReadFile<SearchTask>(tomlFromSelected, MetaMorpheusTask.tomlConfig);
                            staticTasksObservableCollection.Add(new PreRunTask(ye1));
                            break;

                        case "Calibrate":
                            var ye2 = Toml.ReadFile<CalibrationTask>(tomlFromSelected, MetaMorpheusTask.tomlConfig);
                            staticTasksObservableCollection.Add(new PreRunTask(ye2));
                            break;

                        case "Gptmd":
                            var ye3 = Toml.ReadFile<GptmdTask>(tomlFromSelected, MetaMorpheusTask.tomlConfig);
                            staticTasksObservableCollection.Add(new PreRunTask(ye3));
                            break;
                    }
                }
            UpdateTaskGuiStuff();
        }

        #endregion Private Methods

    }
}