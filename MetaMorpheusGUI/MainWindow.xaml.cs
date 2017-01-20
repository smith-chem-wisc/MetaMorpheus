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

        #region Private Fields

        private readonly ObservableCollection<RawData> rawDataObservableCollection = new ObservableCollection<RawData>();
        private readonly ObservableCollection<XMLdb> xmlDBobservableCollection = new ObservableCollection<XMLdb>();
        private readonly ObservableCollection<ModList> modListObservableCollection = new ObservableCollection<ModList>();
        private readonly ObservableCollection<SearchMode> searchModeObservableCollection = new ObservableCollection<SearchMode>();
        private readonly ObservableCollection<FinishedFile> finishedFileObservableCollection = new ObservableCollection<FinishedFile>();
        private readonly ObservableCollection<MyTaskEngine> taskEngineObservableCollection = new ObservableCollection<MyTaskEngine>();

        #endregion Private Fields

        #region Public Constructors

        public MainWindow()
        {
            InitializeComponent();

            if (MyEngine.MetaMorpheusVersion.Equals("1.0.0.0"))
                this.Title = "MetaMorpheus: Not a release version";
            else
                this.Title = "MetaMorpheus: version " + MyEngine.MetaMorpheusVersion;

            dataGridXMLs.DataContext = xmlDBobservableCollection;
            dataGridDatafiles.DataContext = rawDataObservableCollection;
            tasksDataGrid.DataContext = taskEngineObservableCollection;
            outputFilesDataGrid.DataContext = finishedFileObservableCollection;

            modListObservableCollection.Add(new ModList("f.txt"));
            modListObservableCollection.Add(new ModList("v.txt"));
            modListObservableCollection.Add(new ModList("ptmlist.txt"));
            modListObservableCollection.Add(new ModList("m.txt"));
            modListObservableCollection.Add(new ModList("r.txt"));
            modListObservableCollection.Add(new ModList("s.txt"));

            LoadSearchModesFromFile();

            //xmlDBobservableCollection.Add(new XMLdb(@"C:\Users\stepa\Data\CalibrationPaperData\OrigData\uniprot-human-reviewed-12-15-2016.xml"));
            //xmlDBobservableCollection.Add(new XMLdb(@"C:\Users\stepa\Data\CalibrationPaperData\OrigData\uniprot-mouse-reviewed-1-17-2017.xml.gz"));
            //xmlDBobservableCollection.Add(new XMLdb(@"C:\Users\stepa\Data\CalibrationPaperData\OrigData\cRAP-11-11-2016.xml"));

            //rawDataObservableCollection.Add(new RawData(@"C:\Users\stepa\Data\CalibrationPaperData\OrigData\Jurkat\120426_Jurkat_highLC_Frac17.raw"));
            //rawDataObservableCollection.Add(new RawData(@"C:\Users\stepa\Data\CalibrationPaperData\Step2\Mouse\Calib-0.1.2\04-29-13_B6_Frac9_9p5uL-Calibrated.mzML"));

            //rawDataObservableCollection.Add(new RawData(@"C:\Users\stepa\Data\CalibrationPaperData\OrigData\Mouse\04-29-13_B6_Frac1_9uL.raw"));
            //rawDataObservableCollection.Add(new RawData(@"C:\Users\stepa\Data\CalibrationPaperData\OrigData\Mouse\04-29-13_B6_Frac2_9p5uL.raw"));
            //rawDataObservableCollection.Add(new RawData(@"C:\Users\stepa\Data\CalibrationPaperData\OrigData\Mouse\04-29-13_B6_Frac3_9p5uL.raw"));
            //rawDataObservableCollection.Add(new RawData(@"C:\Users\stepa\Data\CalibrationPaperData\OrigData\Mouse\04-29-13_B6_Frac4_8uL.raw"));
            //rawDataObservableCollection.Add(new RawData(@"C:\Users\stepa\Data\CalibrationPaperData\OrigData\Mouse\04-29-13_B6_Frac5_4uL.raw"));
            //rawDataObservableCollection.Add(new RawData(@"C:\Users\stepa\Data\CalibrationPaperData\OrigData\Mouse\04-29-13_B6_Frac6_5uL.raw"));
            //rawDataObservableCollection.Add(new RawData(@"C:\Users\stepa\Data\CalibrationPaperData\OrigData\Mouse\04-29-13_B6_Frac7_5uL.raw"));
            //rawDataObservableCollection.Add(new RawData(@"C:\Users\stepa\Data\CalibrationPaperData\OrigData\Mouse\04-29-13_B6_Frac8_9p5uL.raw"));
            //rawDataObservableCollection.Add(new RawData(@"C:\Users\stepa\Data\CalibrationPaperData\OrigData\Mouse\04-29-13_B6_Frac9_9p5uL.raw"));
            //rawDataObservableCollection.Add(new RawData(@"C:\Users\stepa\Data\CalibrationPaperData\OrigData\Mouse\04-30-13_CAST_Frac1_9uL.raw"));
            //rawDataObservableCollection.Add(new RawData(@"C:\Users\stepa\Data\CalibrationPaperData\OrigData\Mouse\04-30-13_CAST_Frac2_9uL.raw"));
            //rawDataObservableCollection.Add(new RawData(@"C:\Users\stepa\Data\CalibrationPaperData\OrigData\Mouse\04-30-13_CAST_Frac3_6uL.raw"));
            //rawDataObservableCollection.Add(new RawData(@"C:\Users\stepa\Data\CalibrationPaperData\OrigData\Mouse\04-30-13_CAST_Frac4_6uL.raw"));
            //rawDataObservableCollection.Add(new RawData(@"C:\Users\stepa\Data\CalibrationPaperData\OrigData\Mouse\04-30-13_CAST_Frac5_4uL.raw"));
            //rawDataObservableCollection.Add(new RawData(@"C:\Users\stepa\Data\CalibrationPaperData\OrigData\Mouse\04-30-13_CAST_Frac6_5uL.raw"));
            //rawDataObservableCollection.Add(new RawData(@"C:\Users\stepa\Data\CalibrationPaperData\OrigData\Mouse\04-30-13_CAST_Frac7_6uL.raw"));
            //rawDataObservableCollection.Add(new RawData(@"C:\Users\stepa\Data\CalibrationPaperData\OrigData\Mouse\04-30-13_CAST_Frac8_9p5uL.raw"));
            //rawDataObservableCollection.Add(new RawData(@"C:\Users\stepa\Data\CalibrationPaperData\OrigData\Mouse\04-30-13_CAST_Frac9_9p5uL.raw"));

            //rawDataObservableCollection.Add(new RawData(@"C:\Users\stepa\Data\CalibrationPaperData\OrigData\Mouse\2017-01-17-13-30-41\Task1Calibrate\04-30-13_CAST_Frac5_4uL-Calibrated.mzML"));

            //rawDataObservableCollection.Add(new RawData(@"C:\Users\stepa\Data\CalibrationPaperData\OrigData\Mouse\04-29-13_B6_Frac9_9p5uL-Calibrated.mzML"));

            EverythingRunnerEngine.newDbsHandler += AddNewDB;
            EverythingRunnerEngine.newSpectrasHandler += AddNewSpectra;

            EverythingRunnerEngine.startingAllTasksEngineHandler += NewSuccessfullyStartingAllTasks;
            EverythingRunnerEngine.finishedAllTasksEngineHandler += NewSuccessfullyFinishedAllTasks;

            MyTaskEngine.StartingSingleTaskHander += Po_startingSingleTaskHander;
            MyTaskEngine.FinishedSingleTaskHandler += Po_finishedSingleTaskHandler;
            MyTaskEngine.FinishedWritingFileHandler += NewSuccessfullyFinishedWritingFile;

            MyEngine.OutProgressHandler += NewoutProgressBar;
            MyEngine.OutLabelStatusHandler += NewoutLabelStatus;
            MyEngine.StartingSingleEngineHander += MyEngine_startingSingleEngineHander;
            MyEngine.FinishedSingleEngineHandler += MyEngine_finishedSingleEngineHandler;

            UpdateTaskGuiStuff();
        }

        #endregion Public Constructors

        #region Private Methods

        private void MyEngine_startingSingleEngineHander(object sender, SingleEngineEventArgs e)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => MyEngine_startingSingleEngineHander(sender, e)));
            }
            else
            {
                statusLabel.Content = "Running " + e.myEngine.GetType().Name + " engine...";
                outProgressBar.IsIndeterminate = true;
            }
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
            searchModeObservableCollection.Add(new SinglePpmAroundZeroSearchMode("5ppmAroundZero", 5));
            searchModeObservableCollection.Add(new DotSearchMode("5ppm", new double[] { 0 }, new Tolerance(ToleranceUnit.PPM, 5)));
            searchModeObservableCollection.Add(new DotSearchMode("10ppm", new double[] { 0 }, new Tolerance(ToleranceUnit.PPM, 10)));
            searchModeObservableCollection.Add(new IntervalSearchMode("twoPointOneDalton", new List<DoubleRange>() { new DoubleRange(-2.1, 2.1) }));
            searchModeObservableCollection.Add(new OpenSearchMode("Open"));
            searchModeObservableCollection.Add(new SingleAbsoluteAroundZeroSearchMode("0.05daltonsaroundzero", 0.05));
        }

        private void AddNewDB(object sender, XmlForTaskListEventArgs e)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => AddNewDB(sender, e)));
            }
            else
            {
                foreach (var uu in xmlDBobservableCollection)
                    uu.Use = false;
                foreach (var uu in e.newDatabases)
                    xmlDBobservableCollection.Add(new XMLdb(uu.FileName));
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
                foreach (var uu in e.StringList)
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
                s.TheTask.IsMySelected = true;
                statusLabel.Content = "Running " + s.TheTask.TaskType + " task";
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
                s.TheTask.IsMySelected = false;
                statusLabel.Content = "Finished " + s.TheTask.TaskType + " task";
                outProgressBar.Value = 100;

                tasksDataGrid.Items.Refresh();
                dataGridDatafiles.Items.Refresh();
                dataGridXMLs.Items.Refresh();
            }
        }

        private void AddFinishedFile(string filepath)
        {
            finishedFileObservableCollection.Add(new FinishedFile(filepath));
            outputFilesDataGrid.Items.Refresh();
        }

        private void ClearRaw_Click(object sender, RoutedEventArgs e)
        {
            rawDataObservableCollection.Clear();
        }

        private void AddXML_Click(object sender, RoutedEventArgs e)
        {
            // Create the OpenFIleDialog object
            Microsoft.Win32.OpenFileDialog openPicker = new Microsoft.Win32.OpenFileDialog();
            openPicker.Filter = "XML Files|*.xml;*.xml.gz";
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
                    case ".gz":
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
            if (hm != null && !string.IsNullOrEmpty(hm.Text))
            {
                System.Diagnostics.Process.Start(hm.Text);
            }
        }

        private void RunAllTasks_Click(object sender, RoutedEventArgs e)
        {
            EverythingRunnerEngine a = new EverythingRunnerEngine(taskEngineObservableCollection.ToList(), rawDataObservableCollection.Where(b => b.Use).Select(b => b.FileName).ToList(), xmlDBobservableCollection.Where(b => b.Use).Select(b => new XmlForTask(b.FileName, b.Contaminant)).ToList());
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
        }

        private void addCalibrateTaskButton_Click(object sender, RoutedEventArgs e)
        {
            var dialog = new CalibrateTaskWindow(modListObservableCollection);
            if (dialog.ShowDialog() == true)
            {
                taskEngineObservableCollection.Add(dialog.TheTask);
                UpdateTaskGuiStuff();
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
                switch (ok.TaskType)
                {
                    case MyTask.Search:
                        var searchDialog = new SearchTaskWindow(ok as SearchTask, modListObservableCollection, searchModeObservableCollection);
                        searchDialog.ShowDialog();
                        break;

                    case MyTask.Gptmd:
                        var gptmddialog = new GPTMDTaskWindow(ok as GPTMDTask, modListObservableCollection);
                        gptmddialog.ShowDialog();
                        break;

                    case MyTask.Calibrate:
                        var calibratedialog = new CalibrateTaskWindow(ok as CalibrationTask, modListObservableCollection);
                        calibratedialog.ShowDialog();
                        break;
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
                outProgressBar.IsIndeterminate = true;
                statusLabel.Content = s.s;
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
                AddFinishedFile(v.writtenFile);
            }
        }

        private void ClearXML_Click(object sender, RoutedEventArgs e)
        {
            xmlDBobservableCollection.Clear();
        }

        #endregion Private Methods

    }
}