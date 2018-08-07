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
using System.ComponentModel;
using System.Threading.Tasks;
using System.Windows.Threading;
using System;
using System.Threading;
using System.Data;
using System.Windows.Data;
using System.Windows.Media;

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
        ICollectionView peptideSpectralMatchesView;
        private readonly DataTable propertyView;
        private DrawParams DrawParameters;
        private string spectraFilePath;
        private string tsvResultsFilePath;
        
        //private readonly BackgroundWorker psmBack = new BackgroundWorker();
        public MetaDraw()
        {
            InitializeComponent();

            mainViewModel = new MainViewModel();
            plotView.DataContext = mainViewModel;
            MessageBox.Show("1234".ToLower());
            peptideSpectralMatches = new ObservableCollection<MetaDrawPsm>(); //Already specified headers according to properties
            propertyView = new DataTable();
            propertyView.Columns.Add("Name", typeof(string));
            propertyView.Columns.Add("Value", typeof(string));
            peptideSpectralMatchesView = CollectionViewSource.GetDefaultView(peptideSpectralMatches);
            dataGridScanNums.DataContext = peptideSpectralMatchesView;
            dataGridProperties.DataContext = propertyView.DefaultView;
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

        private void LoadPsms(string filename)
        {
            string fileNameWithExtension = Path.GetFileName(spectraFilePath);
            string fileNameWithoutExtension = Path.GetFileNameWithoutExtension(spectraFilePath);
                        
            var psms = TsvResultReader.ReadTsv(filename);

            foreach (var psm in psms)
                {
                    
                    if (psm.FileName == fileNameWithExtension || psm.FileName == fileNameWithoutExtension)
                    {
                        Dispatcher.BeginInvoke(new Action(() =>
                        {
                            peptideSpectralMatches.Add(psm);
                        }));
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
            if (dataGridScanNums.SelectedItem == null)
                return;
            // draw the selected PSM
            propertyView.Clear();
            MetaDrawPsm row = (MetaDrawPsm)dataGridScanNums.SelectedItem;
            System.Reflection.PropertyInfo[] temp = row.GetType().GetProperties();
            
            for (int i = 4; i < temp.Length; i++)
            {
                propertyView.Rows.Add(temp[i].Name, temp[i].GetValue(row,null));
            }
            dataGridProperties.Items.Refresh();
            DrawPsm(row.ScanNum, row.FullSequence);

            //BaseDraw.clearCanvas(cav);
            drawBaseSeq();
            
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

        private async void loadFilesButton_Click(object sender, RoutedEventArgs e)
        {
            propertyView.Clear();
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
            (sender as Button).IsEnabled = false;
            selectSpectraFileButton.IsEnabled = false;
            selectPsmFileButton.IsEnabled = false;
            prgsFeed.IsOpen = true;
            prgsText.Content = "Loading Spectrum Files...";

            //initial Progress
            var slowProcess = Task<MsDataFile>.Factory.StartNew(() =>spectraFileManager.LoadFile(spectraFilePath, null, null, false, false, new CommonParameters()) );
            //start process
            await slowProcess;
            MsDataFile = slowProcess.Result;
            
            //update Progress indicator
            this.prgsText.Content = "Loading PSM...";
            
            // load the PSMs
            await Task.Run(() => LoadPsms(tsvResultsFilePath));

           

            //restore controls
            this.prgsFeed.IsOpen = false;
            (sender as Button).IsEnabled = true;
            selectSpectraFileButton.IsEnabled = true;
            selectPsmFileButton.IsEnabled = true;
        }
        
        private void TextBox_TextChanged(object sender, TextChangedEventArgs e)
        {
            string txt = (sender as TextBox).Text;
            if (txt == "")
                peptideSpectralMatchesView.Filter = null;
            else
            {
                peptideSpectralMatchesView.Filter = obj => { MetaDrawPsm mdp = obj as MetaDrawPsm; return (("" + mdp.ScanNum).StartsWith(txt) || mdp.FullSequence.ToUpper().Contains(txt.ToUpper())); };
            }            
            //dataGridScanNums.Items.Refresh();
        }

        private void drawBaseSeq()
        {
            string peptide = "PEPTIDE";
            double aaPos = 0;
            double fragPos = 0;
            double modPos = 0;

            for (int r = 0; r < peptide.Length; r++)
            {
                aaPos += 3;
                aaPos += BaseDraw.txtDrawing(cav, new Point(aaPos, 10), peptide[r].ToString(), Brushes.Black);

                if (r == 0)
                {
                    fragPos = aaPos + 2;
                    BaseDraw.topSplittingDrawing(cav, new Point(fragPos, 2), Colors.Blue, "y2");
                }
                if (r == 2)
                {
                    fragPos = aaPos + 2;
                    BaseDraw.botSplittingDrawing(cav, new Point(fragPos, 50), Colors.Purple, "y1");
                    BaseDraw.topSplittingDrawing(cav, new Point(fragPos, 2), Colors.Blue, "y2");
                }
                if (r == 4)
                {
                    modPos = aaPos - 2;
                    BaseDraw.circledTxtDraw(cav, new Point(modPos, 11), Brushes.Blue);
                }
            }
        }

        private void dataGridProperties_SelectedCellsChanged(object sender, SelectedCellsChangedEventArgs e)
        {
            (sender as DataGrid).UnselectAll();
        }
    }
}
