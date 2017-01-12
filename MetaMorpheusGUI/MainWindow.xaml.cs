using InternalLogicEngineLayer;
using InternalLogicTaskLayer;
using Spectra;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.IO;
using System.Linq;
using System.Threading;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        private readonly ObservableCollection<RawData> rawDataObservableCollection = new ObservableCollection<RawData>();
        private readonly ObservableCollection<XMLdb> xmlDBobservableCollection = new ObservableCollection<XMLdb>();
        private readonly ObservableCollection<ModList> modListObservableCollection = new ObservableCollection<ModList>();
        private readonly ObservableCollection<SearchMode> searchModeObservableCollection = new ObservableCollection<SearchMode>();
        private readonly ObservableCollection<FinishedFile> finishedFileObservableCollection = new ObservableCollection<FinishedFile>();
        private readonly ObservableCollection<MyTaskEngine> taskEngineObservableCollection = new ObservableCollection<MyTaskEngine>();

        public const string elementsLocation = @"elements.dat";
        public const string unimodLocation = @"unimod_tables.xml";
        public const string uniprotLocation = @"ptmlist.txt";

        public MainWindow()
        {
            InitializeComponent();

            UsefulProteomicsDatabases.Loaders.LoadElements(elementsLocation);
            MyEngine.unimodDeserialized = UsefulProteomicsDatabases.Loaders.LoadUnimod(unimodLocation);
            MyEngine.uniprotDeseralized = UsefulProteomicsDatabases.Loaders.LoadUniprot(uniprotLocation);

            dataGridXMLs.DataContext = xmlDBobservableCollection;
            dataGridDatafiles.DataContext = rawDataObservableCollection;
            tasksDataGrid.DataContext = taskEngineObservableCollection;
            outputFilesDataGrid.DataContext = finishedFileObservableCollection;

            modListObservableCollection.Add(new ModList("f.txt"));
            modListObservableCollection.Add(new ModList("v.txt"));
            modListObservableCollection.Add(new ModList("p.txt"));
            modListObservableCollection.Add(new ModList("m.txt"));
            modListObservableCollection.Add(new ModList("r.txt"));
            modListObservableCollection.Add(new ModList("s.txt"));

            LoadSearchModesFromFile();

            xmlDBobservableCollection.Add(new XMLdb(@"C:\Users\stepa\Data\CalibrationPaperData\OrigData\uniprot-mouse-reviewed-12-23-2016.xml"));

            rawDataObservableCollection.Add(new RawData(@"C:\Users\stepa\Data\CalibrationPaperData\Step2\Mouse\Calib-0.1.2\04-29-13_B6_Frac9_9p5uL-Calibrated.mzML"));

            EverythingRunnerEngine.newDbsHandler += AddNewDB;
            EverythingRunnerEngine.newSpectrasHandler += AddNewSpectra;

            EverythingRunnerEngine.startingAllTasksEngineHandler += NewSuccessfullyStartingAllTasks;
            EverythingRunnerEngine.finishedAllTasksEngineHandler += NewSuccessfullyFinishedAllTasks;

            MyTaskEngine.startingSingleTaskHander += Po_startingSingleTaskHander;
            MyTaskEngine.finishedSingleTaskHandler += Po_finishedSingleTaskHandler;
            MyTaskEngine.finishedWritingFileHandler += NewSuccessfullyFinishedWritingFile;

            MyEngine.outProgressHandler += NewoutProgressBar;
            MyEngine.outLabelStatusHandler += NewoutLabelStatus;
            MyEngine.finishedSingleEngineHandler += MyEngine_finishedSingleEngineHandler;

            UpdateTaskGuiStuff();
        }

        private void MyEngine_finishedSingleEngineHandler(object sender, SingleEngineFinishedEventArgs e)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => MyEngine_finishedSingleEngineHandler(sender, e)));
            }
            else
            {
                outRichTextBox.AppendText(e.ToString() + Environment.NewLine);
                outRichTextBox.ScrollToEnd();
            }
        }

        private void LoadSearchModesFromFile()
        {
            searchModeObservableCollection.Add(new DotSearchMode("5ppm", new double[] { 0 }, new Tolerance(ToleranceUnit.PPM, 5)));
            searchModeObservableCollection.Add(new DotSearchMode("10ppm", new double[] { 0 }, new Tolerance(ToleranceUnit.PPM, 10)));
            searchModeObservableCollection.Add(new IntervalSearchMode("twoPointOneDalton", new List<DoubleRange>() { new DoubleRange(-2.1, 2.1) }));
        }

        private void AddNewDB(object sender, List<string> e)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => AddNewDB(sender, e)));
            }
            else
            {
                foreach (var uu in xmlDBobservableCollection)
                    uu.Use = false;
                foreach (var uu in e)
                    xmlDBobservableCollection.Add(new XMLdb(uu));
            }
        }

        private void AddNewSpectra(object sender, List<string> e)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => AddNewSpectra(sender, e)));
            }
            else
            {
                foreach (var uu in rawDataObservableCollection)
                    uu.Use = false;
                foreach (var uu in e)
                    rawDataObservableCollection.Add(new RawData(uu));
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
                s.theTask.IsMySelected = true;
                statusLabel.Content = "Running " + s.theTask.taskType + " task";
                outProgressBar.IsIndeterminate = true;

                tasksDataGrid.Items.Refresh();
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
                s.theTask.IsMySelected = false;
                statusLabel.Content = "Finished " + s.theTask.taskType + " task";
                outProgressBar.Value = 100;

                tasksDataGrid.Items.Refresh();
                dataGridDatafiles.Items.Refresh();
                dataGridXMLs.Items.Refresh();
            }
        }
        private void addFinishedFile(string filepath)
        {
            finishedFileObservableCollection.Add(new FinishedFile(filepath));
            outputFilesDataGrid.Items.Refresh();
        }

        private RawData GetCorrespondingRawDataAndResultsEntry(string filepath)
        {
            var fileNameNoExtension = Path.GetFileNameWithoutExtension(filepath);
            foreach (var a in rawDataObservableCollection)
            {
                if (a.FileName != null)
                {
                    var aNoExtension = Path.GetFileNameWithoutExtension(a.FileName);
                    if (aNoExtension.Equals(fileNameNoExtension))
                        return a;
                }
            }
            return null;
        }

        private void ClearRaw_Click(object sender, RoutedEventArgs e)
        {
            rawDataObservableCollection.Clear();
        }


        private void AddXML_Click(object sender, RoutedEventArgs e)
        {
            // Create the OpenFIleDialog object
            Microsoft.Win32.OpenFileDialog openPicker = new Microsoft.Win32.OpenFileDialog();
            openPicker.Filter = "XML Files|*.xml";
            openPicker.FilterIndex = 1;
            openPicker.RestoreDirectory = true;
            if (openPicker.ShowDialog() == true)
            {
                xmlDBobservableCollection.Add(new XMLdb(openPicker.FileName));
            }
            dataGridXMLs.Items.Refresh();
        }

        private void AddRaw_Click(object sender, RoutedEventArgs e)
        {
            Microsoft.Win32.OpenFileDialog openFileDialog1 = new Microsoft.Win32.OpenFileDialog();

            openFileDialog1.Filter = "Spectra Files(*.raw;*.mzML;*.mzid;*.psmtsv)|*.raw;*.mzML;*.mzid;*.psmtsv";
            openFileDialog1.FilterIndex = 1;
            openFileDialog1.RestoreDirectory = true;
            openFileDialog1.Multiselect = true;

            if (openFileDialog1.ShowDialog() == true)
                foreach (var filepath in openFileDialog1.FileNames)
                    rawDataObservableCollection.Add(new RawData(filepath));
            dataGridDatafiles.Items.Refresh();
        }

        private void Window_Drop(object sender, DragEventArgs e)
        {
            string[] files = (string[])e.Data.GetData(DataFormats.FileDrop);
            foreach (var file in files)
            {
                var theExtension = Path.GetExtension(file);
                switch (theExtension)
                {
                    case ".raw":
                    case ".mzML":
                        rawDataObservableCollection.Add(new RawData(file));
                        break;

                    case ".xml":
                        xmlDBobservableCollection.Add(new XMLdb(file));
                        break;
                }
                dataGridDatafiles.Items.Refresh();
            }
        }

        private void Row_DoubleClick(object sender, MouseButtonEventArgs e)
        {
            var ye = sender as DataGridCell;
            var hm = ye.Content as TextBlock;
            if (hm != null && !hm.Text.Equals(""))
            {
                System.Diagnostics.Process.Start(hm.Text);
            }
        }

        private void RunAllTasks_Click(object sender, RoutedEventArgs e)
        {
            EverythingRunnerEngine a = new EverythingRunnerEngine(taskEngineObservableCollection.ToList(), rawDataObservableCollection.Where(b => b.Use).Select(b => b.FileName).ToList(), xmlDBobservableCollection.Where(b => b.Use).Select(b => b.FileName).ToList());
            var t = new Thread(() => a.Run());
            t.IsBackground = true;
            t.Start();
        }

        private void ClearTasks_Click(object sender, RoutedEventArgs e)
        {
            taskEngineObservableCollection.Clear();
            UpdateTaskGuiStuff();
        }

        private void UpdateTaskGuiStuff()
        {
            if (taskEngineObservableCollection.Count == 0)
            {
                RunTasksButton.IsEnabled = false;
                RemoveLastTaskButton.IsEnabled = false;
                ClearTasksButton.IsEnabled = false;
            }
            else
            {
                RunTasksButton.IsEnabled = true;
                RemoveLastTaskButton.IsEnabled = true;
                ClearTasksButton.IsEnabled = true;
            }
        }

        private void addSearchTaskButton_Click(object sender, RoutedEventArgs e)
        {
            var dialog = new SearchTaskWindow(modListObservableCollection, searchModeObservableCollection);
            if (dialog.ShowDialog() == true)
            {
                taskEngineObservableCollection.Add(dialog.TheTask);
                UpdateTaskGuiStuff();
            }
            else
            {
            }
        }

        private void addCalibrateTaskButton_Click(object sender, RoutedEventArgs e)
        {
            var dialog = new CalibrateTaskWindow(modListObservableCollection);
            if (dialog.ShowDialog() == true)
            {
                taskEngineObservableCollection.Add(dialog.TheTask);
                UpdateTaskGuiStuff();
            }
            else
            {
            }
        }

        private void addGPTMDTaskButton_Click(object sender, RoutedEventArgs e)
        {
            var dialog = new GPTMDTaskWindow(modListObservableCollection);
            if (dialog.ShowDialog() == true)
            {
                taskEngineObservableCollection.Add(dialog.TheTask);
                UpdateTaskGuiStuff();
            }
            else
            {
            }
        }

        private void RemoveLastTask_Click(object sender, RoutedEventArgs e)
        {
            taskEngineObservableCollection.RemoveAt(taskEngineObservableCollection.Count - 1);
            UpdateTaskGuiStuff();
        }

        private void tasksDataGrid_MouseDoubleClick(object sender, MouseButtonEventArgs e)
        {
            var a = sender as DataGrid;
            var ok = (MyTaskEngine)a.SelectedItem;
            if (ok != null)
                switch (ok.taskType)
                {
                    case MyTaskEnum.Search:
                        var searchDialog = new SearchTaskWindow(ok as SearchTask, modListObservableCollection, searchModeObservableCollection);
                        searchDialog.ShowDialog();
                        break;

                    case MyTaskEnum.GPTMD:
                        var gptmddialog = new GPTMDTaskWindow(ok as GPTMDTask, modListObservableCollection);
                        gptmddialog.ShowDialog();
                        break;

                    case MyTaskEnum.Calibrate:
                        var calibratedialog = new CalibrateTaskWindow(ok as CalibrationTask, modListObservableCollection);
                        calibratedialog.ShowDialog();
                        break;
                }
        }
        private void NewoutLabelStatus(object sender, string s)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => NewoutLabelStatus(sender, s)));
            }
            else
            {
                outProgressBar.IsIndeterminate = true;
                statusLabel.Content = s;
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
                outProgressBar.IsIndeterminate = false;
                outProgressBar.Value = s.new_progress;
                statusLabel.Content = s.v;
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
                //TODO: Check those
                XMLdbPanel.IsEnabled = false;
                DatafilesStackPanel.IsEnabled = false;
                addSearchTaskButton.IsEnabled = false;
                addCalibrateTaskButton.IsEnabled = false;
                addGPTMDTaskButton.IsEnabled = false;
                tasksPanel.IsEnabled = false;

                statusLabel.Content = "Starting all tasks...";
                outProgressBar.IsIndeterminate = true;

                dataGridDatafiles.Items.Refresh();
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
                //TODO: Check those
                XMLdbPanel.IsEnabled = true;
                DatafilesStackPanel.IsEnabled = true;
                addSearchTaskButton.IsEnabled = true;
                addCalibrateTaskButton.IsEnabled = true;
                addGPTMDTaskButton.IsEnabled = true;
                tasksPanel.IsEnabled = true;

                statusLabel.Content = "Finished all tasks!";
                outProgressBar.IsIndeterminate = false;
                outProgressBar.Value = 100;

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
                addFinishedFile(v.writtenFile);
            }
        }

    }
}
