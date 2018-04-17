using EngineLayer;
using MzLibUtil;
using Nett;
using Newtonsoft.Json.Linq;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.IO;
using System.Linq;
using System.Net.Http;
using System.Threading.Tasks;
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
        private readonly ObservableCollection<PreRunTask> staticTasksObservableCollection = new ObservableCollection<PreRunTask>();
        private readonly ObservableCollection<RawDataForDataGrid> SelectedRawFiles = new ObservableCollection<RawDataForDataGrid>();
        private ObservableCollection<InRunTask> dynamicTasksObservableCollection;

        #endregion Private Fields

        #region Public Constructors

        public MainWindow()
        {
            InitializeComponent();

            Title = "MetaMorpheus: version " + GlobalVariables.MetaMorpheusVersion;

            dataGridXMLs.DataContext = proteinDbObservableCollection;
            dataGridDatafiles.DataContext = rawDataObservableCollection;
            tasksTreeView.DataContext = staticTasksObservableCollection;

            EverythingRunnerEngine.NewDbsHandler += AddNewDB;
            EverythingRunnerEngine.NewSpectrasHandler += AddNewSpectra;
            EverythingRunnerEngine.NewFileSpecificTomlHandler += AddNewFileSpecificToml;
            EverythingRunnerEngine.StartingAllTasksEngineHandler += NewSuccessfullyStartingAllTasks;
            EverythingRunnerEngine.FinishedAllTasksEngineHandler += NewSuccessfullyFinishedAllTasks;
            EverythingRunnerEngine.WarnHandler += GuiWarnHandler;
            EverythingRunnerEngine.FinishedWritingAllResultsFileHandler += EverythingRunnerEngine_FinishedWritingAllResultsFileHandler;

            MetaMorpheusTask.StartingSingleTaskHander += Po_startingSingleTaskHander;
            MetaMorpheusTask.FinishedSingleTaskHandler += Po_finishedSingleTaskHandler;
            MetaMorpheusTask.FinishedWritingFileHandler += NewSuccessfullyFinishedWritingFile;
            MetaMorpheusTask.StartingDataFileHandler += MyTaskEngine_StartingDataFileHandler;
            MetaMorpheusTask.FinishedDataFileHandler += MyTaskEngine_FinishedDataFileHandler;
            MetaMorpheusTask.OutLabelStatusHandler += NewoutLabelStatus;
            MetaMorpheusTask.NewCollectionHandler += NewCollectionHandler;
            MetaMorpheusTask.OutProgressHandler += NewoutProgressBar;
            MetaMorpheusTask.WarnHandler += GuiWarnHandler;

            MetaMorpheusEngine.OutProgressHandler += NewoutProgressBar;
            MetaMorpheusEngine.OutLabelStatusHandler += NewoutLabelStatus;
            MetaMorpheusEngine.WarnHandler += GuiWarnHandler;

            MyFileManager.WarnHandler += GuiWarnHandler;

            UpdateRawFileGuiStuff();
            UpdateTaskGuiStuff();
            UpdateOutputFolderTextbox();
            FileSpecificParameters.ValidateFileSpecificVariableNames();

            // LOAD GUI SETTINGS
            GuiGlobalParams = Toml.ReadFile<GuiGlobalParams>(Path.Combine(GlobalVariables.DataDir, @"GUIsettings.toml"));

            if (GlobalVariables.MetaMorpheusVersion.Contains("Not a release version"))
                GuiGlobalParams.AskAboutUpdating = false;

            try
            {
                GetVersionNumbersFromWeb();
            }
            catch (Exception e)
            {
                GuiWarnHandler(null, new StringEventArgs("Could not get newest MM version from web: " + e.Message, null));
            }
        }

        #endregion Public Constructors

        #region Public Properties

        public static string NewestKnownVersion { get; private set; }

        #endregion Public Properties

        #region Internal Properties

        internal GuiGlobalParams GuiGlobalParams { get; }

        #endregion Internal Properties

        #region Private Methods

        private static void GetVersionNumbersFromWeb()
        {
            // Attempt to get current MetaMorpheus version
            using (var client = new HttpClient())
            {
                client.DefaultRequestHeaders.Add("User-Agent", "Mozilla/5.0 (compatible; MSIE 10.0; Windows NT 6.2; WOW64; Trident/6.0)");

                using (var response = client.GetAsync("https://api.github.com/repos/smith-chem-wisc/MetaMorpheus/releases/latest").Result)
                {
                    var json = response.Content.ReadAsStringAsync().Result;
                    JObject deserialized = JObject.Parse(json);
                    var assets = deserialized["assets"].Select(b => b["name"].ToString()).ToList();
                    if (!assets.Contains("MetaMorpheusInstaller.msi"))
                        throw new MetaMorpheusException("Necessary files do not exist!");
                    NewestKnownVersion = deserialized["tag_name"].ToString();
                }
            }
        }

        private void MyWindow_Loaded(object sender, RoutedEventArgs e)
        {
            if (NewestKnownVersion != null && !GlobalVariables.MetaMorpheusVersion.Equals(NewestKnownVersion) && GuiGlobalParams.AskAboutUpdating)
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

            this.KeyDown += new KeyEventHandler(Window_KeyDown);

            // hide the "InProgress" column
            dataGridXMLs.Columns.Where(p => p.Header.Equals(nameof(ProteinDbForDataGrid.InProgress))).First().Visibility = Visibility.Hidden;
            dataGridDatafiles.Columns.Where(p => p.Header.Equals(nameof(RawDataForDataGrid.InProgress))).First().Visibility = Visibility.Hidden;
        }

        private void EverythingRunnerEngine_FinishedWritingAllResultsFileHandler(object sender, StringEventArgs e)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => EverythingRunnerEngine_FinishedWritingAllResultsFileHandler(sender, e)));
            }
            else
            {
                dynamicTasksObservableCollection.Add(new InRunTask("All Task Results", null));
                dynamicTasksObservableCollection.Last().Progress = 100;
                dynamicTasksObservableCollection.Last().Children.Add(new OutputFileForTreeView(e.S, Path.GetFileNameWithoutExtension(e.S)));
            }
        }

        private void GuiWarnHandler(object sender, StringEventArgs e)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => GuiWarnHandler(sender, e)));
            }
            else
            {
                outRichTextBox.AppendText(e.S);
                outRichTextBox.AppendText(Environment.NewLine);
            }
        }

        private void MyTaskEngine_FinishedDataFileHandler(object sender, StringEventArgs s)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => MyTaskEngine_FinishedDataFileHandler(sender, s)));
            }
            else
            {
                var huh = rawDataObservableCollection.First(b => b.FilePath.Equals(s.S));
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
                var huh = rawDataObservableCollection.First(b => b.FilePath.Equals(s.S));
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
                dataGridXMLs.Items.Refresh();
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
                {
                    uu.Use = false;
                }
                foreach (var newRawData in e.StringList)
                    rawDataObservableCollection.Add(new RawDataForDataGrid(newRawData));
                UpdateOutputFolderTextbox();
            }
        }

        private void AddNewFileSpecificToml(object sender, StringListEventArgs e)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => AddNewFileSpecificToml(sender, e)));
            }
            else
            {
                foreach (var path in e.StringList)
                {
                    UpdateFileSpecificParamsDisplayJustAdded(path);
                }
            }
        }

        private void UpdateOutputFolderTextbox()
        {
            if (rawDataObservableCollection.Any())
            {
                // if current output folder is blank and there is a spectra file, use the spectra file's path as the output path
                if (string.IsNullOrWhiteSpace(OutputFolderTextBox.Text))
                {
                    var pathOfFirstSpectraFile = Path.GetDirectoryName(rawDataObservableCollection.First().FilePath);
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

        private void Po_startingSingleTaskHander(object sender, SingleTaskEventArgs s)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => Po_startingSingleTaskHander(sender, s)));
            }
            else
            {
                var theTask = dynamicTasksObservableCollection.First(b => b.DisplayName.Equals(s.DisplayName));
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
                var theTask = dynamicTasksObservableCollection.First(b => b.DisplayName.Equals(s.DisplayName));
                theTask.IsIndeterminate = false;
                theTask.Progress = 100;
                theTask.Status = "Done!";

                dataGridDatafiles.Items.Refresh();
                dataGridXMLs.Items.Refresh();
            }
        }

        private void ClearRaw_Click(object sender, RoutedEventArgs e)
        {
            rawDataObservableCollection.Clear();
            UpdateOutputFolderTextbox();
        }

        private void AddXML_Click(object sender, RoutedEventArgs e)
        {
            Microsoft.Win32.OpenFileDialog openPicker = new Microsoft.Win32.OpenFileDialog()
            {
                Filter = "Database Files|*.xml;*.xml.gz;*.fasta;*.fa",
                FilterIndex = 1,
                RestoreDirectory = true,
                Multiselect = true
            };
            if (openPicker.ShowDialog() == true)
                foreach (var filepath in openPicker.FileNames.OrderBy(p => p))
                {
                    AddAFile(filepath);
                }
            dataGridXMLs.Items.Refresh();
        }

        private void AddRaw_Click(object sender, RoutedEventArgs e)
        {
            Microsoft.Win32.OpenFileDialog openFileDialog1 = new Microsoft.Win32.OpenFileDialog
            {
                Filter = "Spectra Files(*.raw;*.mzML)|*.raw;*.mzML",
                FilterIndex = 1,
                RestoreDirectory = true,
                Multiselect = true
            };
            if (openFileDialog1.ShowDialog() == true)
                foreach (var rawDataFromSelected in openFileDialog1.FileNames.OrderBy(p => p))
                {
                    AddAFile(rawDataFromSelected);
                }
            dataGridDatafiles.Items.Refresh();
        }

        private void Window_Drop(object sender, DragEventArgs e)
        {
            if (LoadTaskButton.IsEnabled)
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
                        dataGridDatafiles.CommitEdit(DataGridEditingUnit.Row, true);
                        dataGridXMLs.CommitEdit(DataGridEditingUnit.Row, true);
                        dataGridDatafiles.Items.Refresh();
                        dataGridXMLs.Items.Refresh();
                    }
                }
                UpdateTaskGuiStuff();
            }
        }

        private void AddAFile(string draggedFilePath)
        {
            // this line is NOT used because .xml.gz (extensions with two dots) mess up with Path.GetExtension
            //var theExtension = Path.GetExtension(draggedFilePath).ToLowerInvariant();

            // we need to get the filename before parsing out the extension because if we assume that everything after the dot
            // is the extension and there are dots in the file path (i.e. in a folder name), this will mess up
            var filename = Path.GetFileName(draggedFilePath);
            var theExtension = filename.Substring(filename.IndexOf(".")).ToLowerInvariant();

            switch (theExtension)
            {
                case ".raw":
                    // check for MSFileReader and display a warning if the expected DLLs are not found
                    var versionCheckerResult = MyFileManager.ValidateThermoMsFileReaderVersion();

                    if (versionCheckerResult.Equals(MyFileManager.ThermoMsFileReaderVersionCheck.IncorrectVersion))
                    {
                        GuiWarnHandler(null, new StringEventArgs("Warning! Thermo MSFileReader is not version 3.0 SP2; a crash may result from searching this .raw file", null));
                    }
                    else if (versionCheckerResult.Equals(MyFileManager.ThermoMsFileReaderVersionCheck.DllsNotFound))
                    {
                        GuiWarnHandler(null, new StringEventArgs("Warning! Cannot find Thermo MSFileReader (v3.0 SP2 is preferred); a crash may result from searching this .raw file", null));
                    }

                    // check for ManagedThermoHelperLayer.dll and display a warning if it's not found
                    // this is one hacky way of checking if the user has C++ redistributable installed
                    string assumedManagedThermoHelperLayerDllPath = Path.Combine(Environment.CurrentDirectory, "ManagedThermoHelperLayer.dll");
                    if(!File.Exists(assumedManagedThermoHelperLayerDllPath))
                    {
                        GuiWarnHandler(null, new StringEventArgs("Warning! You may not have Microsoft Visual C++ Redistributable installed; MetaMorpheus needs this to read .raw files and may crash", null));
                    }

                    goto case ".mzml";

                case ".mzml":
                    RawDataForDataGrid zz = new RawDataForDataGrid(draggedFilePath);
                    if (!ExistRaw(rawDataObservableCollection, zz)) { rawDataObservableCollection.Add(zz); }
                    UpdateFileSpecificParamsDisplayJustAdded(Path.ChangeExtension(draggedFilePath, ".toml"));
                    UpdateOutputFolderTextbox();
                    break;

                case ".mzml.gz":  // not implemented yet
                case ".fasta.gz": // not implemented yet
                    GuiWarnHandler(null, new StringEventArgs("Cannot read, try uncompressing: " + draggedFilePath, null));
                    break;

                case ".xml":
                case ".xml.gz":
                case ".fasta":
                case ".fa":
                    ProteinDbForDataGrid uu = new ProteinDbForDataGrid(draggedFilePath);

                    if (!ExistDa(proteinDbObservableCollection, uu))
                    {
                        proteinDbObservableCollection.Add(uu);
                        if (theExtension.Equals(".xml") || theExtension.Equals(".xml.gz"))
                        {
                            try
                            {
                                GlobalVariables.AddMods(UsefulProteomicsDatabases.ProteinDbLoader.GetPtmListFromProteinXml(draggedFilePath).OfType<ModificationWithLocation>());
                            }
                            catch (Exception ee)
                            {
                                MessageBox.Show(ee.ToString());
                                GuiWarnHandler(null, new StringEventArgs("Cannot read: " + draggedFilePath, null));
                                proteinDbObservableCollection.Remove(uu);
                            }
                        }
                    }
                    break;

                case ".toml":
                    TomlTable tomlFile = null;
                    try
                    {
                        tomlFile = Toml.ReadFile(draggedFilePath, MetaMorpheusTask.tomlConfig);
                    }
                    catch (Exception)
                    {
                        GuiWarnHandler(null, new StringEventArgs("Cannot read: " + draggedFilePath, null));
                        break;
                    }
                    if (tomlFile.Keys.Contains("PrecursorMassTolerance") && tomlFile.Keys.Contains("ProductMassTolerance") && tomlFile.Keys.Count == 2)
                    {
                        // do nothing; it's a ppm suggested tolerance toml from calibration, this gets read in elsewhere
                    }
                    // need to check for file-specific parameters toml here... right now it's displaying a warning
                    else
                    {
                        try
                        {
                            switch (tomlFile.Get<string>("TaskType"))
                            {
                                case "Search":
                                    var ye1 = Toml.ReadFile<SearchTask>(draggedFilePath, MetaMorpheusTask.tomlConfig);
                                    AddTaskToCollection(ye1);
                                    break;

                                case "Calibrate":
                                    var ye2 = Toml.ReadFile<CalibrationTask>(draggedFilePath, MetaMorpheusTask.tomlConfig);
                                    AddTaskToCollection(ye2);
                                    break;

                                case "Gptmd":
                                    var ye3 = Toml.ReadFile<GptmdTask>(draggedFilePath, MetaMorpheusTask.tomlConfig);
                                    AddTaskToCollection(ye3);
                                    break;

                                case "XLSearch":
                                    var ye4 = Toml.ReadFile<XLSearchTask>(draggedFilePath, MetaMorpheusTask.tomlConfig);
                                    AddTaskToCollection(ye4);
                                    break;

                                case "Neo":
                                    var ye5 = Toml.ReadFile<NeoSearchTask>(draggedFilePath, MetaMorpheusTask.tomlConfig);
                                    foreach (MetaMorpheusTask task in NeoLoadTomls.LoadTomls(ye5))
                                        AddTaskToCollection(task);
                                    break;
                            }
                        }
                        catch (Exception e)
                        {
                            GuiWarnHandler(null, new StringEventArgs("Could not parse .toml: " + e.Message, null));
                        }
                    }
                    break;
                default:
                    GuiWarnHandler(null, new StringEventArgs("Unrecognized file type: " + theExtension, null));
                    break;
            }
        }

        private void AddTaskToCollection(MetaMorpheusTask ye)
        {
            PreRunTask te = new PreRunTask(ye);
            staticTasksObservableCollection.Add(te);
            staticTasksObservableCollection.Last().DisplayName = "Task" + (staticTasksObservableCollection.IndexOf(te) + 1) + "-" + ye.CommonParameters.TaskDescriptor;
        }

        private void InsertTaskToCollectionAtPosition(MetaMorpheusTask ye, int index, out int actualIndexInserted)
        {
            PreRunTask te = new PreRunTask(ye);
            if (index < 0)
            {
                index = 0;
            }

            if (staticTasksObservableCollection.Count == 0)
            {
                AddTaskToCollection(ye);
                actualIndexInserted = 0;
                return;
            }

            if (index > staticTasksObservableCollection.Count - 1)
            {
                index = staticTasksObservableCollection.Count - 1;
            }

            actualIndexInserted = index;
            staticTasksObservableCollection.Insert(index, te);
            staticTasksObservableCollection[index].DisplayName = "Task" + (staticTasksObservableCollection.IndexOf(te) + 1) + "-" + ye.CommonParameters.TaskDescriptor;
        }

        private void Row_DoubleClick(object sender, MouseButtonEventArgs e)
        {
            var ye = sender as DataGridCell;
            
            // prevent opening protein DB or spectra files if a run is in progress
            if ((ye.DataContext is ProteinDbForDataGrid || ye.DataContext is RawDataForDataGrid) && !LoadTaskButton.IsEnabled)
            {
                return;
            }

            if (ye.Content is TextBlock hm && hm != null && !string.IsNullOrEmpty(hm.Text))
            {
                try
                {
                    System.Diagnostics.Process.Start(hm.Text);
                }
                catch (Exception)
                {

                }
            }
        }

        private void RunAllTasks_Click(object sender, RoutedEventArgs e)
        {
            dynamicTasksObservableCollection = new ObservableCollection<InRunTask>();

            for (int i = 0; i < staticTasksObservableCollection.Count; i++)
            {
                dynamicTasksObservableCollection.Add(new InRunTask("Task" + (i + 1) + "-" + staticTasksObservableCollection[i].metaMorpheusTask.CommonParameters.TaskDescriptor, staticTasksObservableCollection[i].metaMorpheusTask));
            }
            tasksTreeView.DataContext = dynamicTasksObservableCollection;

            outRichTextBox.Document.Blocks.Clear();

            EverythingRunnerEngine a = new EverythingRunnerEngine(dynamicTasksObservableCollection.Select(b => (b.DisplayName, b.task)).ToList(), rawDataObservableCollection.Where(b => b.Use).Select(b => b.FilePath).ToList(), proteinDbObservableCollection.Where(b => b.Use).Select(b => new DbForTask(b.FilePath, b.Contaminant)).ToList(), OutputFolderTextBox.Text);
            var t = new Task(a.Run);
            t.ContinueWith(EverythingRunnerExceptionHandler, TaskContinuationOptions.OnlyOnFaulted);
            t.Start();
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
                while (e.InnerException != null) e = e.InnerException;
                var message = "Run failed, Exception: " + e.Message;
                var messageBoxResult = System.Windows.MessageBox.Show(message + "\n\nWould you like to report this crash?", "Runtime Error", MessageBoxButton.YesNo);
                outRichTextBox.AppendText(message + Environment.NewLine);
                Exception exception = e;
                //Find Output Folder
                string outputFolder = e.Data["folder"].ToString();
                string tomlText = "";
                if (Directory.Exists(outputFolder))
                {
                    var tomls = Directory.GetFiles(outputFolder, "*.toml");
                    //will only be 1 toml per task
                    foreach (var tomlFile in tomls)
                        tomlText += "\n" + File.ReadAllText(tomlFile);
                    if (!tomls.Any())
                        tomlText = "TOML not found";
                }
                else
                    tomlText = "Directory not found";
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
                    System.Diagnostics.Process.Start(mailto);
                    Console.WriteLine(body);
                }
            }
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
                DeleteSelectedTaskButton.IsEnabled = false;
                ClearTasksButton.IsEnabled = false;
                ResetTasksButton.IsEnabled = false;
            }
            else
            {
                RunTasksButton.IsEnabled = true;
                DeleteSelectedTaskButton.IsEnabled = true;
                ClearTasksButton.IsEnabled = true;

                // this exists so that when a task is deleted, the remaining tasks are renamed to keep the task numbers correct
                for(int i = 0; i < staticTasksObservableCollection.Count; i++)
                {
                    string newName = "Task" + (i + 1) + "-" + staticTasksObservableCollection[i].metaMorpheusTask.CommonParameters.TaskDescriptor;
                    staticTasksObservableCollection[i].DisplayName = newName;
                }
                tasksTreeView.Items.Refresh();
            }
        }

        private void UpdateRawFileGuiStuff()
        {
            ChangeFileParameters.IsEnabled = SelectedRawFiles.Count > 0 && LoadTaskButton.IsEnabled;
        }

        private void AddSearchTaskButton_Click(object sender, RoutedEventArgs e)
        {
            var dialog = new SearchTaskWindow();
            if (dialog.ShowDialog() == true)
            {
                AddTaskToCollection(dialog.TheTask);
                UpdateTaskGuiStuff();
            }
        }

        private void AddCalibrateTaskButton_Click(object sender, RoutedEventArgs e)
        {
            var dialog = new CalibrateTaskWindow();
            if (dialog.ShowDialog() == true)
            {
                AddTaskToCollection(dialog.TheTask);
                UpdateTaskGuiStuff();
            }
        }

        private void AddGPTMDTaskButton_Click(object sender, RoutedEventArgs e)
        {
            var dialog = new GptmdTaskWindow();
            if (dialog.ShowDialog() == true)
            {
                AddTaskToCollection(dialog.TheTask);
                UpdateTaskGuiStuff();
            }
        }

        private void BtnAddCrosslinkSearch_Click(object sender, RoutedEventArgs e)
        {
            var dialog = new XLSearchTaskWindow();
            if (dialog.ShowDialog() == true)
            {
                AddTaskToCollection(dialog.TheTask);
                UpdateTaskGuiStuff();
            }
        }

        private void AddNeoTaskButton_Click(object sender, RoutedEventArgs e)
        {
            var dialog = new NeoSearchTaskWindow();
            if (dialog.ShowDialog() == true)
            {
                var ye5 = dialog.TheTask;
                foreach (MetaMorpheusTask task in NeoLoadTomls.LoadTomls(ye5))
                    AddTaskToCollection(task);
                UpdateTaskGuiStuff();
            }
        }

        private void DeleteSelectedTask(object sender, RoutedEventArgs e)
        {
            var selectedTask = (PreRunTask)tasksTreeView.SelectedItem;
            if (selectedTask != null)
            {
                staticTasksObservableCollection.Remove(selectedTask);
                UpdateTaskGuiStuff();
            }
        }

        private void MoveSelectedTask(object sender, RoutedEventArgs e, bool moveTaskUp)
        {
            var selectedTask = (PreRunTask)tasksTreeView.SelectedItem;
            if (selectedTask == null)
            {
                return;
            }
            
            int indexOfSelected = staticTasksObservableCollection.IndexOf(selectedTask);
            int indexToMoveTo = indexOfSelected - 1;
            if(moveTaskUp)
            {
                indexToMoveTo = indexOfSelected + 1;
            }

            if (indexToMoveTo >= 0 && indexToMoveTo < staticTasksObservableCollection.Count)
            {
                var temp = staticTasksObservableCollection[indexToMoveTo];
                staticTasksObservableCollection[indexToMoveTo] = selectedTask;
                staticTasksObservableCollection[indexOfSelected] = temp;

                UpdateTaskGuiStuff();

                var item = tasksTreeView.ItemContainerGenerator.ContainerFromItem(selectedTask);
                ((TreeViewItem)item).IsSelected = true;
            }
        }

        // handles keyboard input in the main window
        private void Window_KeyDown(object sender, KeyEventArgs e)
        {
            if (LoadTaskButton.IsEnabled)
            {
                // delete selected task
                if (e.Key == Key.Delete || e.Key == Key.Back)
                {
                    DeleteSelectedTask(sender, e);
                    e.Handled = true;
                }

                // move task up
                if(e.Key == Key.Add)
                {
                    MoveSelectedTask(sender, e, true);
                    e.Handled = true;
                }

                // move task down
                if (e.Key == Key.Subtract)
                {
                    MoveSelectedTask(sender, e, false);
                    e.Handled = true;
                }
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

                theEntityOnWhichToUpdateLabel.Children.Add(new CollectionForTreeView(s.S, s.nestedIDs.Last()));
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

                theEntityOnWhichToUpdateLabel.Status = s.S;
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
                dataGridDatafiles.Items.Refresh();

                ClearTasksButton.IsEnabled = false;
                DeleteSelectedTaskButton.IsEnabled = false;
                RunTasksButton.IsEnabled = false;
                LoadTaskButton.IsEnabled = false;

                addCalibrateTaskButton.IsEnabled = false;
                addGPTMDTaskButton.IsEnabled = false;
                addSearchTaskButton.IsEnabled = false;
                btnAddCrosslinkSearch.IsEnabled = false;
                //addNeoTaskButton.IsEnabled = false;

                AddXML.IsEnabled = false;
                ClearXML.IsEnabled = false;
                AddRaw.IsEnabled = false;
                ClearRaw.IsEnabled = false;

                OutputFolderTextBox.IsEnabled = false;

                dataGridDatafiles.IsReadOnly = true;
                dataGridXMLs.IsReadOnly = true;
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
                ForTreeView AddWrittenFileToThisOne = dynamicTasksObservableCollection.First(b => b.Id.Equals(v.nestedIDs[0]));

                foreach (var hm in v.nestedIDs.Skip(1))
                {
                    try
                    {
                        AddWrittenFileToThisOne = AddWrittenFileToThisOne.Children.First(b => b.Id.Equals(hm));
                    }
                    catch
                    {
                    }
                }
                AddWrittenFileToThisOne.Children.Add(new OutputFileForTreeView(v.writtenFile, Path.GetFileName(v.writtenFile)));
            }
        }

        private void ClearXML_Click(object sender, RoutedEventArgs e)
        {
            proteinDbObservableCollection.Clear();
        }

        private void ResetTasksButton_Click(object sender, RoutedEventArgs e)
        {
            tasksGroupBox.IsEnabled = true;
            ClearTasksButton.IsEnabled = true;
            DeleteSelectedTaskButton.IsEnabled = true;
            RunTasksButton.IsEnabled = true;
            addCalibrateTaskButton.IsEnabled = true;
            addGPTMDTaskButton.IsEnabled = true;
            addSearchTaskButton.IsEnabled = true;
            btnAddCrosslinkSearch.IsEnabled = true;
            //addNeoTaskButton.IsEnabled = true;
            ResetTasksButton.IsEnabled = false;
            OutputFolderTextBox.IsEnabled = true;

            dataGridDatafiles.IsReadOnly = false;
            dataGridXMLs.IsReadOnly = false;

            AddXML.IsEnabled = true;
            ClearXML.IsEnabled = true;
            AddRaw.IsEnabled = true;
            ClearRaw.IsEnabled = true;

            LoadTaskButton.IsEnabled = true;

            tasksTreeView.DataContext = staticTasksObservableCollection;
            UpdateRawFileGuiStuff();
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

                        preRunTask.DisplayName = "Task" + (staticTasksObservableCollection.IndexOf(preRunTask) + 1) + "-" + searchDialog.TheTask.CommonParameters.TaskDescriptor;
                        tasksTreeView.Items.Refresh();

                        return;

                    case MyTask.Gptmd:
                        var gptmddialog = new GptmdTaskWindow(preRunTask.metaMorpheusTask as GptmdTask);
                        gptmddialog.ShowDialog();
                        preRunTask.DisplayName = "Task" + (staticTasksObservableCollection.IndexOf(preRunTask) + 1) + "-" + gptmddialog.TheTask.CommonParameters.TaskDescriptor;
                        tasksTreeView.Items.Refresh();

                        return;

                    case MyTask.Calibrate:
                        var calibratedialog = new CalibrateTaskWindow(preRunTask.metaMorpheusTask as CalibrationTask);
                        calibratedialog.ShowDialog();
                        preRunTask.DisplayName = "Task" + (staticTasksObservableCollection.IndexOf(preRunTask) + 1) + "-" + calibratedialog.TheTask.CommonParameters.TaskDescriptor;
                        tasksTreeView.Items.Refresh();
                        return;

                    case MyTask.XLSearch:
                        var XLSearchdialog = new XLSearchTaskWindow(preRunTask.metaMorpheusTask as XLSearchTask);
                        XLSearchdialog.ShowDialog();
                        preRunTask.DisplayName = "Task" + (staticTasksObservableCollection.IndexOf(preRunTask) + 1) + "-" + XLSearchdialog.TheTask.CommonParameters.TaskDescriptor;
                        tasksTreeView.Items.Refresh();
                        return;

                    case MyTask.Neo:
                        var Neodialog = new NeoSearchTaskWindow(preRunTask.metaMorpheusTask as NeoSearchTask);
                        Neodialog.ShowDialog();
                        preRunTask.DisplayName = "Task" + (staticTasksObservableCollection.IndexOf(preRunTask) + 1) + "-" + Neodialog.TheTask.CommonParameters.TaskDescriptor;
                        tasksTreeView.Items.Refresh();
                        return;
                }

            if (a.SelectedItem is OutputFileForTreeView fileThing)
            {
                System.Diagnostics.Process.Start(fileThing.FullPath);
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
                foreach (var tomlFromSelected in openFileDialog1.FileNames.OrderBy(p => p))
                {
                    AddAFile(tomlFromSelected);
                }
            UpdateTaskGuiStuff();
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
                                GuiWarnHandler(null, new StringEventArgs("Problem parsing the file-specific toml " + Path.GetFileName(fullPathofTomls[j]) + "; " + e.Message + "; is the toml from an older version of MetaMorpheus?", null));
                            }
                        }
                        else
                        {
                            file.SetParametersText(null);
                        }
                    }
                }
            }
            UpdateRawFileGuiStuff();
            dataGridDatafiles.Items.Refresh();
        }

        //run if data file has just been added with and checks for Existing fileSpecficParams
        private void UpdateFileSpecificParamsDisplayJustAdded(string tomlLocation)
        {
            string assumedPathOfSpectraFileWithoutExtension = Path.Combine(Directory.GetParent(tomlLocation).ToString(), Path.GetFileNameWithoutExtension(tomlLocation));

            for (int i = 0; i < rawDataObservableCollection.Count; i++)
            {
                string thisFilesPathWihoutExtension = Path.Combine(Directory.GetParent(rawDataObservableCollection[i].FilePath).ToString(), Path.GetFileNameWithoutExtension(rawDataObservableCollection[i].FilePath));
                if (File.Exists(tomlLocation) && assumedPathOfSpectraFileWithoutExtension.Equals(thisFilesPathWihoutExtension))
                {
                    TomlTable fileSpecificSettings = Toml.ReadFile(tomlLocation, MetaMorpheusTask.tomlConfig);
                    try
                    {
                        // parse to make sure toml is readable
                        var temp = new FileSpecificParameters(fileSpecificSettings);

                        // toml is ok; display the file-specific settings in the gui
                        rawDataObservableCollection[i].SetParametersText(File.ReadAllText(tomlLocation));
                    }
                    catch (MetaMorpheusException e)
                    {
                        GuiWarnHandler(null, new StringEventArgs("Problem parsing the file-specific toml " + Path.GetFileName(tomlLocation) + "; " + e.Message + "; is the toml from an older version of MetaMorpheus?", null));
                    }
                }
            }
            UpdateRawFileGuiStuff();
            dataGridDatafiles.Items.Refresh();
        }

        private void AddSelectedRaw(object sender, RoutedEventArgs e)
        {
            DataGridRow obj = (DataGridRow)sender;

            RawDataForDataGrid ok = (RawDataForDataGrid)obj.DataContext;
            SelectedRawFiles.Add(ok);
            UpdateRawFileGuiStuff();
        }

        private void RemoveSelectedRaw(object sender, RoutedEventArgs e)
        {
            DataGridRow obj = (DataGridRow)sender;
            RawDataForDataGrid ok = (RawDataForDataGrid)obj.DataContext;
            SelectedRawFiles.Remove(ok);
            UpdateRawFileGuiStuff();
        }

        private void MenuItem_Click(object sender, RoutedEventArgs e)
        {
            System.Diagnostics.Process.Start(@"https://github.com/smith-chem-wisc/MetaMorpheus/wiki");
        }

        private void MenuItem_Click_1(object sender, RoutedEventArgs e)
        {
            var globalSettingsDialog = new GlobalSettingsWindow();
            globalSettingsDialog.ShowDialog();
        }

        private bool ExistDa(ObservableCollection<ProteinDbForDataGrid> pDOC, ProteinDbForDataGrid uuu)
        {
            foreach (ProteinDbForDataGrid pdoc in pDOC)
                if (pdoc.FilePath == uuu.FilePath) { return true; }
            return false;
        }

        private bool ExistRaw(ObservableCollection<RawDataForDataGrid> rDOC, RawDataForDataGrid zzz)
        {
            foreach (RawDataForDataGrid rdoc in rDOC)
                if (rdoc.FileName == zzz.FileName) { return true; }
            return false;
        }

        private void ChangeFileParameters_Click(object sender, RoutedEventArgs e)
        {
            var dialog = new FileSpecificParametersWindow(SelectedRawFiles);
            if (dialog.ShowDialog() == true)
            {
                var tomlPathsForSelectedFiles = SelectedRawFiles.Select(p => Path.Combine(Directory.GetParent(p.FilePath).ToString(), Path.GetFileNameWithoutExtension(p.FileName)) + ".toml");
                UpdateFileSpecificParamsDisplay(tomlPathsForSelectedFiles.ToArray());
            }
        }

        private void MenuItem_Click_2(object sender, RoutedEventArgs e)
        {
            try
            {
                GetVersionNumbersFromWeb();
            }
            catch (Exception ex)
            {
                GuiWarnHandler(null, new StringEventArgs("Could not get newest MM version from web: " + ex.Message, null));
                return;
            }

            if (GlobalVariables.MetaMorpheusVersion.Equals(NewestKnownVersion))
                MessageBox.Show("You have the most updated version!");
            else
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
        }

        private void MenuItem_Click_4(object sender, RoutedEventArgs e)
        {
            string mailto = string.Format("mailto:{0}?Subject=MetaMorpheus. Issue:", "mm_support@chem.wisc.edu");
            System.Diagnostics.Process.Start(mailto);
        }

        private void MenuItem_Click_5(object sender, RoutedEventArgs e)
        {
            System.Diagnostics.Process.Start(@"https://github.com/smith-chem-wisc/MetaMorpheus/issues/new");
        }

        private void MenuItem_Click_6(object sender, RoutedEventArgs e)
        {
            System.Diagnostics.Process.Start(GlobalVariables.DataDir);
        }

        private void MenuItem_Click_3(object sender, RoutedEventArgs e)
        {
            System.Diagnostics.Process.Start(Path.Combine(GlobalVariables.DataDir, @"settings.toml"));
            Application.Current.Shutdown();
        }

        private void MenuItem_Click_7(object sender, RoutedEventArgs e)
        {
            System.Diagnostics.Process.Start(Path.Combine(GlobalVariables.DataDir, @"GUIsettings.toml"));
            Application.Current.Shutdown();
        }

        #endregion Private Methods
    }
}