using EngineLayer;
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
        private bool WarnedAboutThermoAlready = false;

        public MainWindow()
        {
            InitializeComponent();

            Title = "MetaMorpheus: version " + GlobalVariables.MetaMorpheusVersion;

            dataGridProteinDatabases.DataContext = ProteinDbObservableCollection;
            dataGridSpectraFiles.DataContext = SpectraFilesObservableCollection;
            tasksTreeView.DataContext = StaticTasksObservableCollection;

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

            UpdateSpectraFileGuiStuff();
            UpdateTaskGuiStuff();
            UpdateOutputFolderTextbox();
            FileSpecificParameters.ValidateFileSpecificVariableNames();
            SearchModifications.SetUpModSearchBoxes();

            // LOAD GUI SETTINGS

            if (File.Exists(Path.Combine(GlobalVariables.DataDir, @"GUIsettings.toml")))
            {
                GuiGlobalParams = Toml.ReadFile<GuiGlobalParams>(Path.Combine(GlobalVariables.DataDir, @"GUIsettings.toml"));
            }
            else
            {
                Toml.WriteFile(GuiGlobalParams, Path.Combine(GlobalVariables.DataDir, @"GUIsettings.toml"), MetaMorpheusTask.tomlConfig);
                notificationsTextBox.Document = YoutubeWikiNotification();
            }

            if (GlobalVariables.MetaMorpheusVersion.Contains("Not a release version"))
                GuiGlobalParams.AskAboutUpdating = false;

            try
            {
                GetVersionNumbersFromWeb();
            }
            catch (Exception e)
            {
                GuiWarnHandler(null, new StringEventArgs("Could not get newest version from web: " + e.Message, null));
            }

            Application.Current.MainWindow.Closing += new CancelEventHandler(MainWindow_Closing);
        }

        private FlowDocument YoutubeWikiNotification()
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

            var links = new List<Hyperlink> {wikiLink, youtubeLink};

            p.Inlines.Add(run1);
            p.Inlines.Add(wikiLink);
            p.Inlines.Add(run3);
            p.Inlines.Add(youtubeLink);
            p.Inlines.Add(run5);

            foreach (Hyperlink link in links)
            {
                link.RequestNavigate += (sender, e) =>
                {
                    System.Diagnostics.Process.Start(e.Uri.ToString());
                };
            }

            doc.Blocks.Add(p);
            return doc;
        }

        public static string NewestKnownVersion { get; private set; }

        internal GuiGlobalParams GuiGlobalParams = new GuiGlobalParams();

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
                        throw new MetaMorpheusException("A new version of MetaMorpheus was detected, but the files haven't been" +
                            " uploaded yet. Try again in a few minutes.");
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
            dataGridProteinDatabases.Columns.Where(p => p.Header.Equals(nameof(ProteinDbForDataGrid.InProgress))).First().Visibility = Visibility.Hidden;
            dataGridSpectraFiles.Columns.Where(p => p.Header.Equals(nameof(RawDataForDataGrid.InProgress))).First().Visibility = Visibility.Hidden;

            PrintErrorsReadingMods();
        }

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

        private void GuiWarnHandler(object sender, StringEventArgs e)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => GuiWarnHandler(sender, e)));
            }
            else
            {
                notificationsTextBox.AppendText(e.S);
                notificationsTextBox.AppendText(Environment.NewLine);
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
                var huh = SpectraFilesObservableCollection.First(b => b.FilePath.Equals(s.S));
                huh.SetInProgress(false);

                dataGridSpectraFiles.Items.Refresh();
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
                var huh = SpectraFilesObservableCollection.First(b => b.FilePath.Equals(s.S));
                huh.SetInProgress(true);
                dataGridSpectraFiles.Items.Refresh();
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
                foreach (var uu in ProteinDbObservableCollection)
                    uu.Use = false;
                foreach (var uu in e.NewDatabases)
                    ProteinDbObservableCollection.Add(new ProteinDbForDataGrid(uu));
                dataGridProteinDatabases.Items.Refresh();
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

        private void ClearSpectraFiles_Click(object sender, RoutedEventArgs e)
        {
            SpectraFilesObservableCollection.Clear();
            UpdateOutputFolderTextbox();
        }

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
                    GuiWarnHandler(null, new StringEventArgs("Error opening directory: " + ex.Message, null));
                }
            }

            if (Directory.Exists(outputFolder))
            {
                // open the directory
                System.Diagnostics.Process.Start(new System.Diagnostics.ProcessStartInfo()
                {
                    FileName = outputFolder,
                    UseShellExecute = true,
                    Verb = "open"
                });
            }
            else
            {
                // this should only happen if the file path is empty or something unexpected happened
                GuiWarnHandler(null, new StringEventArgs("Output folder does not exist", null));
            }
        }

        private void AddProteinDatabase_Click(object sender, RoutedEventArgs e)
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
            dataGridProteinDatabases.Items.Refresh();
        }

        private void AddSpectraFile_Click(object sender, RoutedEventArgs e)
        {
            Microsoft.Win32.OpenFileDialog openFileDialog1 = new Microsoft.Win32.OpenFileDialog
            {
                Filter = "Spectra Files(*.raw;*.mzML;*.mgf)|*.raw;*.mzML;*.mgf",
                FilterIndex = 1,
                RestoreDirectory = true,
                Multiselect = true
            };
            if (openFileDialog1.ShowDialog() == true)
                foreach (var rawDataFromSelected in openFileDialog1.FileNames.OrderBy(p => p))
                {
                    AddAFile(rawDataFromSelected);
                }
            dataGridSpectraFiles.Items.Refresh();
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
                        dataGridSpectraFiles.CommitEdit(DataGridEditingUnit.Row, true);
                        dataGridProteinDatabases.CommitEdit(DataGridEditingUnit.Row, true);
                        dataGridSpectraFiles.Items.Refresh();
                        dataGridProteinDatabases.Items.Refresh();
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
            var theExtension = Path.GetExtension(filename).ToLowerInvariant();
            bool compressed = theExtension.EndsWith("gz"); // allows for .bgz and .tgz, too which are used on occasion
            theExtension = compressed ? Path.GetExtension(Path.GetFileNameWithoutExtension(filename)).ToLowerInvariant() : theExtension;

            switch (theExtension)
            {
                case ".raw":
                    if (!WarnedAboutThermoAlready)
                    {
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
                        else if (versionCheckerResult.Equals(MyFileManager.ThermoMsFileReaderVersionCheck.SomeDllsMissing))
                        {
                            GuiWarnHandler(null, new StringEventArgs("Warning! Found only some of the expected Thermo MSFileReader .dll files; a crash may result from searching this .raw file", null));
                        }

                        // check for ManagedThermoHelperLayer.dll and display a warning if it's not found
                        // this is one hacky way of checking if the user has C++ redistributable installed
                        string assumedManagedThermoHelperLayerDllPath = Path.Combine(Environment.CurrentDirectory, "ManagedThermoHelperLayer.dll");
                        if (!File.Exists(assumedManagedThermoHelperLayerDllPath))
                        {
                            GuiWarnHandler(null, new StringEventArgs("Warning! Cannot find Microsoft Visual C++ Redistributable; " +
                                "a crash may result from searching this .raw file. If you have just installed the C++ redistributable, " +
                                "please uninstall and reinstall MetaMorpheus", null));
                        }
                    }

                    WarnedAboutThermoAlready = true;
                    goto case ".mzml";

                case ".mgf":
                    GuiWarnHandler(null, new StringEventArgs(".mgf files lack MS1 spectra, which are needed for quantification and searching for coisolated peptides. All other features of MetaMorpheus will function.", null));
                    goto case ".mzml";

                case ".mzml":
                    if (compressed) // not implemented yet
                    {
                        GuiWarnHandler(null, new StringEventArgs("Cannot read, try uncompressing: " + draggedFilePath, null));
                        break;
                    }
                    RawDataForDataGrid zz = new RawDataForDataGrid(draggedFilePath);
                    if (!SpectraFileExists(SpectraFilesObservableCollection, zz))
                    {
                        SpectraFilesObservableCollection.Add(zz);
                    }
                    UpdateFileSpecificParamsDisplayJustAdded(Path.ChangeExtension(draggedFilePath, ".toml"));
                    UpdateOutputFolderTextbox();
                    break;

                case ".xml":
                case ".fasta":
                case ".fa":
                    ProteinDbForDataGrid uu = new ProteinDbForDataGrid(draggedFilePath);
                    if (!DatabaseExists(ProteinDbObservableCollection, uu))
                    {
                        ProteinDbObservableCollection.Add(uu);
                        if (theExtension.Equals(".xml"))
                        {
                            try
                            {
                                GlobalVariables.AddMods(UsefulProteomicsDatabases.ProteinDbLoader.GetPtmListFromProteinXml(draggedFilePath).OfType<Modification>(), true);

                                PrintErrorsReadingMods();
                            }
                            catch (Exception ee)
                            {
                                MessageBox.Show(ee.ToString());
                                GuiWarnHandler(null, new StringEventArgs("Cannot parse modification info from: " + draggedFilePath, null));
                                ProteinDbObservableCollection.Remove(uu);
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
                        GuiWarnHandler(null, new StringEventArgs("Cannot read toml: " + draggedFilePath, null));
                        break;
                    }

                    if (tomlFile.ContainsKey("TaskType"))
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
                            }
                        }
                        catch (Exception e)
                        {
                            GuiWarnHandler(null, new StringEventArgs("Cannot read task toml: " + e.Message, null));
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
            StaticTasksObservableCollection.Add(te);
            StaticTasksObservableCollection.Last().DisplayName = "Task" + (StaticTasksObservableCollection.IndexOf(te) + 1) + "-" + ye.CommonParameters.TaskDescriptor;
        }

        // handles double-clicking on a data grid row
        private void Row_DoubleClick(object sender, MouseButtonEventArgs e)
        {
            var ye = sender as DataGridCell;

            // prevent opening protein DB or spectra files if a run is in progress
            if ((ye.DataContext is ProteinDbForDataGrid || ye.DataContext is RawDataForDataGrid) && !LoadTaskButton.IsEnabled)
            {
                return;
            }

            // open the file with the default process for this file format
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
            GlobalVariables.StopLoops = false;
            CancelButton.IsEnabled = true;

            // check for valid tasks/spectra files/protein databases
            if (!StaticTasksObservableCollection.Any())
            {
                GuiWarnHandler(null, new StringEventArgs("You need to add at least one task!", null));
                return;
            }
            if (!SpectraFilesObservableCollection.Any())
            {
                GuiWarnHandler(null, new StringEventArgs("You need to add at least one spectra file!", null));
                return;
            }
            if (!ProteinDbObservableCollection.Any())
            {
                GuiWarnHandler(null, new StringEventArgs("You need to add at least one protein database!", null));
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
            BtnQuantSet.IsEnabled = false;

            // everything is OK to run
            EverythingRunnerEngine a = new EverythingRunnerEngine(DynamicTasksObservableCollection.Select(b => (b.DisplayName, b.Task)).ToList(),
                SpectraFilesObservableCollection.Where(b => b.Use).Select(b => b.FilePath).ToList(),
                ProteinDbObservableCollection.Where(b => b.Use).Select(b => new DbForTask(b.FilePath, b.Contaminant)).ToList(),
                outputFolder);

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
                ResetTasksButton.IsEnabled = true;
            }
        }

        private void ClearTasks_Click(object sender, RoutedEventArgs e)
        {
            StaticTasksObservableCollection.Clear();
            UpdateTaskGuiStuff();
        }

        private void UpdateTaskGuiStuff()
        {
            if (StaticTasksObservableCollection.Count == 0)
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

        // deletes the selected task
        private void DeleteSelectedTask(object sender, RoutedEventArgs e)
        {
            var selectedTask = (PreRunTask)tasksTreeView.SelectedItem;
            if (selectedTask != null)
            {
                StaticTasksObservableCollection.Remove(selectedTask);
                UpdateTaskGuiStuff();
            }
        }

        // move the task up or down in the GUI
        private void MoveSelectedTask(object sender, RoutedEventArgs e, bool moveTaskUp)
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
                if (e.Key == Key.Add || e.Key == Key.OemPlus)
                {
                    MoveSelectedTask(sender, e, true);
                    e.Handled = true;
                }

                // move task down
                if (e.Key == Key.Subtract || e.Key == Key.OemMinus)
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

                ClearTasksButton.IsEnabled = false;
                DeleteSelectedTaskButton.IsEnabled = false;
                RunTasksButton.IsEnabled = false;
                LoadTaskButton.IsEnabled = false;

                addCalibrateTaskButton.IsEnabled = false;
                addGPTMDTaskButton.IsEnabled = false;
                addSearchTaskButton.IsEnabled = false;
                btnAddCrosslinkSearch.IsEnabled = false;

                AddXML.IsEnabled = false;
                ClearXML.IsEnabled = false;
                AddRaw.IsEnabled = false;
                ClearRaw.IsEnabled = false;

                OutputFolderTextBox.IsEnabled = false;

                dataGridSpectraFiles.IsReadOnly = true;
                dataGridProteinDatabases.IsReadOnly = true;
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

        private void ClearXML_Click(object sender, RoutedEventArgs e)
        {
            ProteinDbObservableCollection.Clear();
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
            ResetTasksButton.IsEnabled = false;
            OutputFolderTextBox.IsEnabled = true;

            dataGridSpectraFiles.IsReadOnly = false;
            dataGridProteinDatabases.IsReadOnly = false;

            AddXML.IsEnabled = true;
            ClearXML.IsEnabled = true;
            AddRaw.IsEnabled = true;
            ClearRaw.IsEnabled = true;
            BtnQuantSet.IsEnabled = true;

            LoadTaskButton.IsEnabled = true;

            tasksTreeView.DataContext = StaticTasksObservableCollection;
            UpdateSpectraFileGuiStuff();

            var pathOfFirstSpectraFile = Path.GetDirectoryName(SpectraFilesObservableCollection.First().FilePath);
            OutputFolderTextBox.Text = Path.Combine(pathOfFirstSpectraFile, @"$DATETIME");
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
                }

            if (a.SelectedItem is OutputFileForTreeView fileThing)
            {
                if (File.Exists(fileThing.FullPath))
                {
                    System.Diagnostics.Process.Start(fileThing.FullPath);
                }
                else
                {
                    MessageBox.Show("File " + Path.GetFileName(fileThing.FullPath) + " does not exist");
                }
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
                        GuiWarnHandler(null, new StringEventArgs("Problem parsing the file-specific toml " + Path.GetFileName(tomlLocation) + "; " + e.Message + "; is the toml from an older version of MetaMorpheus?", null));
                    }
                }
            }
            UpdateSpectraFileGuiStuff();
            dataGridSpectraFiles.Items.Refresh();
        }

        private void AddSelectedSpectra(object sender, RoutedEventArgs e)
        {
            DataGridRow obj = (DataGridRow)sender;

            RawDataForDataGrid ok = (RawDataForDataGrid)obj.DataContext;
            SelectedRawFiles.Add(ok);
            UpdateSpectraFileGuiStuff();
        }

        private void RemoveSelectedSpectra(object sender, RoutedEventArgs e)
        {
            DataGridRow obj = (DataGridRow)sender;
            RawDataForDataGrid ok = (RawDataForDataGrid)obj.DataContext;
            SelectedRawFiles.Remove(ok);
            UpdateSpectraFileGuiStuff();
        }

        private void MenuItem_Click(object sender, RoutedEventArgs e)
        {
            System.Diagnostics.Process.Start(@"https://github.com/smith-chem-wisc/MetaMorpheus/wiki");
        }

        private void MenuItem_YouTube(object sender, RoutedEventArgs e)
        {
            System.Diagnostics.Process.Start(@"https://www.youtube.com/playlist?list=PLVk5tTSZ1aWlhNPh7jxPQ8pc0ElyzSUQb");
        }

        private void MenuItem_ProteomicsNewsBlog(object sender, RoutedEventArgs e)
        {
            System.Diagnostics.Process.Start(@"https://proteomicsnews.blogspot.com/");
        }

        private void MenuItem_Click_1(object sender, RoutedEventArgs e)
        {
            var globalSettingsDialog = new GlobalSettingsWindow();
            globalSettingsDialog.ShowDialog();
        }

        private bool DatabaseExists(ObservableCollection<ProteinDbForDataGrid> pDOC, ProteinDbForDataGrid uuu)
        {
            foreach (ProteinDbForDataGrid pdoc in pDOC)
                if (pdoc.FilePath == uuu.FilePath) { return true; }
            return false;
        }

        private bool SpectraFileExists(ObservableCollection<RawDataForDataGrid> rDOC, RawDataForDataGrid zzz)
        {
            foreach (RawDataForDataGrid rdoc in rDOC)
                if (rdoc.FileName == zzz.FileName) { return true; }
            return false;
        }

        private void CancelButton_Click(object sender, RoutedEventArgs e)
        {
            string grammar = StaticTasksObservableCollection.Count <= 1 ? "this task" : "these tasks";
            if (MessageBox.Show("Are you sure you want to cancel " + grammar + "?", "Cancel Tasks", MessageBoxButton.OKCancel) == MessageBoxResult.OK)
            {
                GlobalVariables.StopLoops = true;
                CancelButton.IsEnabled = false;
                notificationsTextBox.AppendText("Canceling...\n");
            }
        }

        private void ChangeFileParameters_Click(object sender, RoutedEventArgs e)
        {
            var dialog = new FileSpecificParametersWindow(SelectedRawFiles);
            if (dialog.ShowDialog() == true)
            {
                var tomlPathsForSelectedFiles = SelectedRawFiles.Select(p => Path.Combine(Directory.GetParent(p.FilePath).ToString(), Path.GetFileNameWithoutExtension(p.FileName)) + ".toml").ToList();
                UpdateFileSpecificParamsDisplay(tomlPathsForSelectedFiles.ToArray());
            }
        }

        private void BtnQuantSet_Click(object sender, RoutedEventArgs e)
        {
            var dialog = new ExperimentalDesignWindow(SpectraFilesObservableCollection);
            dialog.ShowDialog();
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

        private void MenuItem_Twitter(object sender, RoutedEventArgs e)
        {
            System.Diagnostics.Process.Start(@"https://twitter.com/Smith_Chem_Wisc");
        }

        private void MenuItem_Slack(object sender, RoutedEventArgs e)
        {
            System.Diagnostics.Process.Start(@"https://join.slack.com/t/smith-chem-public/shared_invite/enQtNDYzNTM5Mzg5NzY0LTRiYWQ5MzVmYmExZWIyMTcyZmNlODJjMWI0YjVhNGM2MmQ2NjE4ZDAzNmM4NWYxMDFhNTQyNDBiM2E0MWE0NGU");
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

        private void MetaDrawMenu_Click(object sender, RoutedEventArgs e)
        {
            MetaDraw metaDrawGui = new MetaDraw();
            metaDrawGui.Show();
        }

        private void PrintErrorsReadingMods()
        {
            // print any error messages reading the mods to the notifications area
            foreach (var error in GlobalVariables.ErrorsReadingMods)
            {
                GuiWarnHandler(null, new StringEventArgs(error, null));
            }
            GlobalVariables.ErrorsReadingMods.Clear();
        }

        private void AddContaminantXML_Click(object sender, RoutedEventArgs e)
        {
            string[] contaminantFiles = Directory.GetFiles(Path.Combine(GlobalVariables.DataDir, "Contaminants"));
            foreach (string contaminantFile in contaminantFiles)
            {
                AddAFile(contaminantFile);
            }
            dataGridProteinDatabases.Items.Refresh();
        }

        private void AddCustomMod_Click(object sender, RoutedEventArgs e)
        {
            var dialog = new CustomModButtonWindow();
            dialog.ShowDialog();
        }

        // handle window closing
        private void MainWindow_Closing(object sender, CancelEventArgs e)
        {
            if (!GuiGlobalParams.DisableCloseWindow && !GlobalVariables.MetaMorpheusVersion.Contains("DEBUG"))
            {
                e.Cancel = true;
                var exit = CustomMsgBox.Show("Exit MetaMorpheus", "Are you sure you want to exit MetaMorpheus?", "Yes", "No", "Yes and don't ask me again");

                if (exit == MessageBoxResult.Yes)
                {
                    e.Cancel = false;
                }
                else if (exit == MessageBoxResult.OK)
                {
                    GuiGlobalParams.DisableCloseWindow = true;
                    Toml.WriteFile(GuiGlobalParams, Path.Combine(GlobalVariables.DataDir, @"GUIsettings.toml"), MetaMorpheusTask.tomlConfig);
                    e.Cancel = false;
                }
            }
        }
    }
}