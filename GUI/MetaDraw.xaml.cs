using System;
using System.Linq;
using System.Windows;
using System.Windows.Controls;
using System.Collections.ObjectModel;
using System.IO;
using ViewModels;
using EngineLayer;
using System.Collections.Generic;
using MzLibUtil;
using MassSpectrometry;
using OxyPlot;
using System.Globalization;
using TaskLayer;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for MetaDraw.xaml
    /// </summary>
    public partial class MetaDraw : Window
    {
        private readonly ObservableCollection<RawDataForDataGrid> spectraFilesObservableCollection = new ObservableCollection<RawDataForDataGrid>();
        private readonly ObservableCollection<RawDataForDataGrid> resultFilesObservableCollection = new ObservableCollection<RawDataForDataGrid>();
        private MainViewModel mainViewModel;
        private MyFileManager fileManager;
        private MsDataFile MsDataFile;
        private List<MetaDrawPsm> peptideSpectralMatches;
        private readonly ObservableCollection<MetaDrawPsm> spectrumNumsObservableCollection = new ObservableCollection<MetaDrawPsm>();
        private DrawParams DrawParameters;

        public MetaDraw()
        {
            InitializeComponent();

            mainViewModel = new MainViewModel();

            plotView.DataContext = mainViewModel;

            dataGridMassSpectraFiles.DataContext = spectraFilesObservableCollection;

            dataGridResultFiles.DataContext = resultFilesObservableCollection;

            dataGridScanNums.DataContext = spectrumNumsObservableCollection;
            
            Title = "MetaDraw: version " + GlobalVariables.MetaMorpheusVersion;

            DrawParameters = new DrawParams();
            productMassToleranceComboBox.Items.Add("Da");
            productMassToleranceComboBox.Items.Add("ppm");
            UpdateFieldsFromPanel();
            fileManager = new MyFileManager(true);
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
                        AddAFile(draggedFilePath);
                    }

                    dataGridMassSpectraFiles.CommitEdit(DataGridEditingUnit.Row, true);
                    dataGridMassSpectraFiles.Items.Refresh();
                }
            }
        }

        private void btnClearFiles_Click(object sender, RoutedEventArgs e)
        {
            if ((sender as Button).Name == "btnClearResultFiles")
                resultFilesObservableCollection.Clear();
            else
                spectraFilesObservableCollection.Clear();
        }

        private void btnAddFiles_Click(object sender, RoutedEventArgs e)
        {
            btnReset.IsEnabled = true;
            Microsoft.Win32.OpenFileDialog openFileDialog1 = new Microsoft.Win32.OpenFileDialog
            {
                Filter = "Spectra Files(*.raw;*.mzML)|*.raw;*.mzML",
                FilterIndex = 1,
                RestoreDirectory = true,
                Multiselect = true
            };
            if (openFileDialog1.ShowDialog() == true)
            {
                foreach (var rawDataFromSelected in openFileDialog1.FileNames.OrderBy(p => p))
                {
                    AddAFile(rawDataFromSelected);
                }
            }
            dataGridMassSpectraFiles.Items.Refresh();
        }

        private void btnAddResults_Click(object sender, RoutedEventArgs e)
        {
            btnReset.IsEnabled = true;
            Microsoft.Win32.OpenFileDialog openFileDialog1 = new Microsoft.Win32.OpenFileDialog
            {
                Filter = "Result Files(*.csv;*..psmtsv)|*.csv;*.psmtsv",
                FilterIndex = 1,
                RestoreDirectory = true,
                Multiselect = true
            };
            if (openFileDialog1.ShowDialog() == true)
            {
                foreach (var rawDataFromSelected in openFileDialog1.FileNames.OrderBy(p => p))
                {
                    AddAFile(rawDataFromSelected);
                }
            }
            dataGridResultFiles.Items.Refresh();
        }

        private void AddAFile(string draggedFilePath)
        {
            var filename = Path.GetFileName(draggedFilePath);
            var theExtension = filename.Substring(filename.IndexOf(".")).ToLowerInvariant();

            switch (theExtension)
            {
                case ".raw":
                case ".mzml":
                    RawDataForDataGrid file = new RawDataForDataGrid(draggedFilePath);
                    if (!SpectraFileExists(spectraFilesObservableCollection, file))
                    {
                        spectraFilesObservableCollection.Add(file);
                    }
                    break;
                case ".pep.xml":
                    break;
                case ".psmtsv":
                case ".tsv":
                    RawDataForDataGrid resultFileDataGrid = new RawDataForDataGrid(draggedFilePath);
                    if (!SpectraFileExists(resultFilesObservableCollection, resultFileDataGrid))
                    {
                        resultFilesObservableCollection.Add(resultFileDataGrid);
                    }
                    break;
                default:
                    break;
            }
        }

        private bool SpectraFileExists(ObservableCollection<RawDataForDataGrid> data, RawDataForDataGrid file)
        {
            return data.Select(p => p.FileName).Contains(file.FileName);
        }

        /*private void UpdateOutputFolderTextbox()
        {
            if (spectraFilesObservableCollection.Any())
            {
                // if current output folder is blank and there is a spectra file, use the spectra file's path as the output path
                if (string.IsNullOrWhiteSpace(txtBoxOutputFolder.Text))
                {
                    var pathOfFirstSpectraFile = Path.GetDirectoryName(spectraFilesObservableCollection.Where(p => p.Use).First().FilePath);
                    txtBoxOutputFolder.Text = Path.Combine(pathOfFirstSpectraFile, @"$DATETIME");
                }
                // else do nothing (do not override if there is a path already there; might clear user-defined path)
            }
            else
            {
                // no spectra files; clear the output folder from the GUI
                txtBoxOutputFolder.Clear();
            }
        }*/

        private void btnDraw_Click(object sender, RoutedEventArgs e)
        {
            mainViewModel.Model.InvalidatePlot(true);

            if (int.TryParse(txtScanNum.Text, out int scanNumber))
            {
                DrawPsm(scanNumber);
            }
            else
            {
                MessageBox.Show("The input scan number must be an integer.");
                return;
            }
        }

        private void btnReadResultFile_Click(object sender, RoutedEventArgs e)
        {
            btnReset.IsEnabled = true;

            if (!spectraFilesObservableCollection.Any())
            {
                MessageBox.Show("Please add spectra files.");
                return;
            }

            if (!resultFilesObservableCollection.Any())
            {
                MessageBox.Show("Please add search result files.");
                return;
            }

            MsDataFile = fileManager.LoadFile(spectraFilesObservableCollection.Where(b => b.Use).First().FilePath, null, null, false, false);

            btnReadResultFile.IsEnabled = true;

            var resultFilePath = resultFilesObservableCollection.Where(b => b.Use).First().FilePath;
            peptideSpectralMatches = TsvResultReader.ReadTsv(resultFilePath);
            foreach (var psm in peptideSpectralMatches)
            {
                spectrumNumsObservableCollection.Add(psm);
            }
            dataGridScanNums.Items.Refresh();

            // hide the "peptide" column because the PeptideWithSetModifications string is too large
            dataGridScanNums.Columns[2].Visibility = Visibility.Hidden;

            btnReadResultFile.IsEnabled = false;
            btnDraw.IsEnabled = true;
        }

        private void DrawPsm(int oneBasedScanNumber, string fullSequence = null)
        {
            if (MsDataFile == null)
            {
                MessageBox.Show("Please check that the mass spectral data is loaded.");
                return;
            }

            MsDataScan msDataScanToDraw = MsDataFile.GetOneBasedScan(oneBasedScanNumber);
            IEnumerable<MetaDrawPsm> scanPsms = peptideSpectralMatches.Where(p => p.ScanNum == oneBasedScanNumber);

            if (fullSequence != null)
            {
                scanPsms = scanPsms.Where(p => p.Peptide.Sequence == fullSequence);
            }

            mainViewModel.DrawPeptideSpectralMatch(msDataScanToDraw, scanPsms.FirstOrDefault(), DrawParameters);
        }

        private void btnReset_Click(object sender, RoutedEventArgs e)
        {
            btnReset.IsEnabled = false;

            resultFilesObservableCollection.Clear();

            spectraFilesObservableCollection.Clear();

            spectrumNumsObservableCollection.Clear();

            btnReadResultFile.IsEnabled = true;

            mainViewModel.Model = new PlotModel { Title = "Spectrum Annotation", Subtitle = "using OxyPlot" };
        }

        private void clearText(object sender, RoutedEventArgs e)
        {
            TextBox tb = (TextBox)sender;
            if (tb.Text.Equals("Scan Number"))
                tb.Text = string.Empty;
        }

        private void restoreText(object sender, RoutedEventArgs e)
        {
            TextBox tb = (TextBox)sender;
            if (tb.Text.Equals(string.Empty))
                tb.Text = "Scan Number";
        }

        private void bCheckBox_Checked(object sender, RoutedEventArgs e)
        {
            DrawParameters.BIons = true;
        }

        private void yCheckBox_Checked(object sender, RoutedEventArgs e)
        {
            DrawParameters.YIons = true;
        }

        private void cCheckBox_Checked(object sender, RoutedEventArgs e)
        {
            DrawParameters.CIons = true;
        }

        private void zdotCheckBox_Checked(object sender, RoutedEventArgs e)
        {
            DrawParameters.ZdotIons = true;
        }

        private void productMassToleranceTextBox_TextChanged(object sender, TextChangedEventArgs e)
        {
            if (double.TryParse(productMassToleranceTextBox.Text, out double tol))
            {
                if (productMassToleranceComboBox.SelectedIndex == 0)
                {
                    DrawParameters.ProductMassTolerance = new AbsoluteTolerance(tol);
                }
                else
                {
                    DrawParameters.ProductMassTolerance = new PpmTolerance(tol);
                }
            }
        }

        private void UpdateFieldsFromPanel()
        {
            bCheckBox.IsChecked = DrawParameters.BIons;
            yCheckBox.IsChecked = DrawParameters.YIons;
            cCheckBox.IsChecked = DrawParameters.CIons;
            zdotCheckBox.IsChecked = DrawParameters.ZdotIons;
            productMassToleranceTextBox.Text = DrawParameters.ProductMassTolerance.Value.ToString(CultureInfo.InvariantCulture);
            productMassToleranceComboBox.SelectedIndex = DrawParameters.ProductMassTolerance is AbsoluteTolerance ? 0 : 1;
        }

        private void dataGridScanNums_SelectedCellsChanged(object sender, SelectedCellsChangedEventArgs e)
        {
            // draw the selected PSM
            MetaDrawPsm row = (MetaDrawPsm)dataGridScanNums.SelectedItem;
            DrawPsm(row.ScanNum, row.Peptide.Sequence);
        }
    }
}
