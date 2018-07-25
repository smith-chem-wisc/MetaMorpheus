using System.Linq;
using System.Windows;
using System.Windows.Controls;
using System.Collections.ObjectModel;
using System.IO;
using ViewModels;
using EngineLayer;
using System.Collections.Generic;
using MassSpectrometry;
using TaskLayer;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for MetaDraw.xaml
    /// </summary>
    public partial class MetaDraw : Window
    {
        private MainViewModel mainViewModel;
        private MyFileManager spectraFileManager;
        private MsDataFile MsDataFile;
        private readonly ObservableCollection<MetaDrawPsm> peptideSpectralMatches;
        private DrawParams DrawParameters;
        private string spectraFilePath;
        private string tsvResultsFilePath;

        public MetaDraw()
        {
            InitializeComponent();

            mainViewModel = new MainViewModel();
            plotView.DataContext = mainViewModel;

            peptideSpectralMatches = new ObservableCollection<MetaDrawPsm>();
            dataGridScanNums.DataContext = peptideSpectralMatches;

            Title = "MetaDraw: version " + GlobalVariables.MetaMorpheusVersion;
            DrawParameters = new DrawParams();
            spectraFileManager = new MyFileManager(true);
        }
        
        private void Window_Drop(object sender, DragEventArgs e)
        {
            string[] files = ((string[])e.Data.GetData(DataFormats.FileDrop)).OrderBy(p => p).ToArray();

            if (files != null)
            {
                foreach (var draggedFilePath in files)
                {
                    if (File.Exists(draggedFilePath))
                    {
                        LoadFile(draggedFilePath);
                    }
                }
            }
        }
        
        private void LoadFile(string filePath)
        {
            var theExtension = Path.GetExtension(filePath).ToLower();

            switch (theExtension)
            {
                case ".raw":
                case ".mzml":
                    spectraFilePath = filePath;
                    spectraFileNameLabel.Text = filePath;
                    break;
                case ".psmtsv":
                    tsvResultsFilePath = filePath;
                    psmFileNameLabel.Text = filePath;
                    break;
                default:
                    MessageBox.Show("Cannot read file type: " + theExtension);
                    break;
            }
        }

        private void LoadPsms(string filePath)
        {
            tsvResultsFilePath = filePath;

            string fileNameWithExtension = Path.GetFileName(spectraFilePath);
            string fileNameWithoutExtension = Path.GetFileNameWithoutExtension(spectraFilePath);

            var psms = TsvResultReader.ReadTsv(filePath);

            foreach (var psm in psms)
            {
                if (psm.FileName == fileNameWithExtension || psm.FileName == fileNameWithoutExtension)
                {
                    peptideSpectralMatches.Add(psm);
                }
            }
        }
        
        private void DrawPsm(int oneBasedScanNumber, string fullSequence = null)
        {
            MsDataScan msDataScanToDraw = MsDataFile.GetOneBasedScan(oneBasedScanNumber);
            IEnumerable<MetaDrawPsm> scanPsms = peptideSpectralMatches.Where(p => p.ScanNum == oneBasedScanNumber);

            if (fullSequence != null)
            {
                scanPsms = scanPsms.Where(p => p.FullSequence == fullSequence);
            }

            mainViewModel.DrawPeptideSpectralMatch(msDataScanToDraw, scanPsms.FirstOrDefault(), DrawParameters);
        }
        
        /// <summary>
        /// Event triggers when a different cell is selected in the PSM data grid
        /// </summary>
        private void dataGridScanNums_SelectedCellsChanged(object sender, SelectedCellsChangedEventArgs e)
        {
            // draw the selected PSM
            MetaDrawPsm row = (MetaDrawPsm)dataGridScanNums.SelectedItem;
            DrawPsm(row.ScanNum, row.FullSequence);
        }

        private void selectSpectraFileButton_Click(object sender, RoutedEventArgs e)
        {
            Microsoft.Win32.OpenFileDialog openFileDialog1 = new Microsoft.Win32.OpenFileDialog
            {
                Filter = "Spectra Files(*.raw;*.mzML)|*.raw;*.mzML",
                FilterIndex = 1,
                RestoreDirectory = true,
                Multiselect = false
            };
            if (openFileDialog1.ShowDialog() == true)
            {
                foreach (var filePath in openFileDialog1.FileNames.OrderBy(p => p))
                {
                    LoadFile(filePath);
                }
            }
        }
        
        private void selectPsmFileButton_Click(object sender, RoutedEventArgs e)
        {
            Microsoft.Win32.OpenFileDialog openFileDialog1 = new Microsoft.Win32.OpenFileDialog
            {
                Filter = "Result Files(*.psmtsv)|*.psmtsv",
                FilterIndex = 1,
                RestoreDirectory = true,
                Multiselect = false
            };
            if (openFileDialog1.ShowDialog() == true)
            {
                foreach (var filePath in openFileDialog1.FileNames.OrderBy(p => p))
                {
                    LoadFile(filePath);
                }
            }
        }
        
        private void loadFilesButton_Click(object sender, RoutedEventArgs e)
        {
            if (spectraFilePath == null)
            {
                MessageBox.Show("Please add a spectra file.");
                return;
            }

            if (tsvResultsFilePath == null)
            {
                MessageBox.Show("Please add a search result file.");
                return;
            }

            // load the spectra file
            MsDataFile = spectraFileManager.LoadFile(spectraFilePath, null, null, false, false, new CommonParameters());

            // load the PSMs
            LoadPsms(tsvResultsFilePath);
            dataGridScanNums.Items.Refresh();

            // hide the filename and ion columns
            dataGridScanNums.Columns[2].Visibility = Visibility.Hidden;
            dataGridScanNums.Columns[3].Visibility = Visibility.Hidden;
        }
    }
}
