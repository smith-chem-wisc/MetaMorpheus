using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Threading;
using System.Globalization;
using System.Collections.ObjectModel;
using TaskLayer;
using System.IO;
using EngineLayer;
using Proteomics;


[assembly: log4net.Config.XmlConfigurator(Watch = true)]

namespace RealTimeGUI
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        private static readonly log4net.ILog log = log4net.LogManager.GetLogger(System.Reflection.MethodBase.GetCurrentMethod().DeclaringType);

        private readonly ObservableCollection<ProteinDbForDataGrid> proteinDbObservableCollection = new ObservableCollection<ProteinDbForDataGrid>();
        private LogWatcher logWatcher;

        public DataReceiver DataReceiver { get; set; }
        public RealTimeTask RealTimeTask { get; set; }


        public MainWindow()
        {
            InitializeComponent();
            RealTimeTask = new RealTimeTask();

            DataReceiver = new DataReceiver();
            dataGridProteinDatabases.DataContext = proteinDbObservableCollection;

            DataReceiver.DataReceiverNotificationEventHandler += UpdateTbNotification;
            MyFileManager.WarnHandler += GuiWarnHandler;
            RealTimeSearch.WarnHandler += GuiWarnHandler;

            // Create a LogFileWatcher to display the log and bind the log textbox to it
            logWatcher = new LogWatcher();
            logWatcher.Updated += logWatcher_Updated;

        }

        private void UpdateTbNotification(object sender, NotificationEventArgs e)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => UpdateTbNotification(sender, e)));
            }
            else
            {
                RtbNotifications.AppendText(e.Notification);
            }    
        }

        private void BtnConnection_Click(object sender, RoutedEventArgs e)
        {
            log.Debug("Start log");
            //DataReceiver.TestLog();
            DataReceiver.InstrumentAccess = Connection.GetFirstInstrument();
            DataReceiver.ScanContainer = DataReceiver.InstrumentAccess.GetMsScanContainer(0);
            RtbNotifications.AppendText(DataReceiver.InstrumentAccess.InstrumentName);
        }

        private void BtnDisConnection_Click(object sender, RoutedEventArgs e)
        {
            //Put the following code into a function 
            DataReceiver.ScanContainer = null;
            DataReceiver.InstrumentAccess = null;
        }

        private void BtnRealTimeData_Click(object sender, RoutedEventArgs e)
        {
            log.Debug("Start Receive Data");
            DataReceiver.ReceiveData();
            Thread.CurrentThread.Join(DataReceiver.RTParameters.TimeScale);
            DataReceiver.StopReceiveData();
            //DataReceiver.TestLog();
        }

        public void logWatcher_Updated(object sender, EventArgs e)
        {
            UpdateLogTextbox(logWatcher.LogContent);
        }

        public void UpdateLogTextbox(string value)
        {
            // Check whether invoke is required and then invoke as necessary
            if (!Dispatcher.CheckAccess())
            {

                Dispatcher.BeginInvoke(new Action(() => UpdateLogTextbox(value)));
                return;
            }

            // Set the textbox value
            RtbNotifications.AppendText(value);
        }

        private void UpdateParametersFromTask()
        {
            TbTimeScale.Text = DataReceiver.RTParameters.TimeScale.ToString(CultureInfo.InvariantCulture);
        }

        private void UpdateParametersFromGUI()
        {
            DataReceiver.RTParameters.TimeScale = int.Parse(TbTimeScale.Text);
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
                dataGridProteinDatabases.Items.Refresh();
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
                case ".xml":
                case ".fasta":
                case ".fa":
                    ProteinDbForDataGrid uu = new ProteinDbForDataGrid(draggedFilePath);
                    if (!DatabaseExists(proteinDbObservableCollection, uu))
                    {
                        proteinDbObservableCollection.Add(uu);
                        if (theExtension.Equals(".xml"))
                        {
                            try
                            {
                                GlobalVariables.AddMods(UsefulProteomicsDatabases.ProteinDbLoader.GetPtmListFromProteinXml(draggedFilePath).OfType<ModificationWithLocation>());
                            }
                            catch (Exception ee)
                            {
                                MessageBox.Show(ee.ToString());
                                GuiWarnHandler(null, new StringEventArgs("Cannot parse modification info from: " + draggedFilePath, null));
                                proteinDbObservableCollection.Remove(uu);
                            }
                        }
                    }
                    break;

                default:
                    GuiWarnHandler(null, new StringEventArgs("Unrecognized file type: " + theExtension, null));
                    break;
            }
        }

        private void Window_Drop(object sender, DragEventArgs e)
        {
            //if (LoadTaskButton.IsEnabled)
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
                        dataGridProteinDatabases.CommitEdit(DataGridEditingUnit.Row, true);
                        dataGridProteinDatabases.Items.Refresh();
                    }
                }
            }
        }

        private bool DatabaseExists(ObservableCollection<ProteinDbForDataGrid> pDOC, ProteinDbForDataGrid uuu)
        {
            foreach (ProteinDbForDataGrid pdoc in pDOC)
                if (pdoc.FilePath == uuu.FilePath) { return true; }
            return false;
        }

        private void GuiWarnHandler(object sender, StringEventArgs e)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => GuiWarnHandler(sender, e)));
            }
            else
            {
                RtbNotifications.AppendText(e.S);
                RtbNotifications.AppendText(Environment.NewLine);
            }
        }

        private void BtnAddXML_Click(object sender, RoutedEventArgs e)
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

        private void BtnClearXML_Click(object sender, RoutedEventArgs e)
        {
            proteinDbObservableCollection.Clear();
        }

        private void BtnDatabasePrepare_Click(object sender, RoutedEventArgs e)
        {
            if (!proteinDbObservableCollection.Any())
            {
                GuiWarnHandler(null, new StringEventArgs("You need to add at least one protein database!", null));
                return;
            }
            string outputFolder = TbOutputFolder.Text;

            RealTimeSearch realTimeSearch = new RealTimeSearch(proteinDbObservableCollection.Where(b=>b.Use).Select(b=>new DbForTask(b.FilePath, b.Contaminant)).ToList(), RealTimeTask, outputFolder);
            var t = new Task(realTimeSearch.Run);
            t.Start();
        }

        private void BtnOpenOutputFolder_Click(object sender, RoutedEventArgs e)
        {

            string outputFolder = TbOutputFolder.Text;

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

        private void BtnDataBaseParameter_Click(object sender, RoutedEventArgs e)
        {
            var dialog = new DatabaseParameterWindow();
            if (dialog.ShowDialog() == true)
            {
                RealTimeTask = dialog.TheTask;
            }
        }
    }
}
