using System;
using System.Linq;
using System.Windows;
using System.Windows.Controls;
using System.Collections.ObjectModel;
using System.IO;
using ViewModels;
using EngineLayer;
using EngineLayer.CrosslinkSearch;
using System.Collections.Generic;
using MzLibUtil;
using System.Text.RegularExpressions;
using MassSpectrometry;
using OxyPlot;
using System.Globalization;

namespace MetaDrawGUI
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        private readonly ObservableCollection<RawDataForDataGrid> spectraFilesObservableCollection = new ObservableCollection<RawDataForDataGrid>();
        private readonly ObservableCollection<RawDataForDataGrid> resultFilesObservableCollection = new ObservableCollection<RawDataForDataGrid>();
        private MainViewModel mainViewModel = null;
        private MsDataFile MsDataFile = null;   
        private List<PsmDraw> PSMs = null;
        private readonly ObservableCollection<SpectrumForDataGrid> spectrumNumsObservableCollection = new ObservableCollection<SpectrumForDataGrid>();
        private CommonParameters CommonParameters;

        public MainWindow()
        {

            InitializeComponent();

            mainViewModel = new MainViewModel();

            plotView.DataContext = mainViewModel;

            dataGridMassSpectraFiles.DataContext = spectraFilesObservableCollection;

            dataGridResultFiles.DataContext = resultFilesObservableCollection;

            dataGridScanNums.DataContext = spectrumNumsObservableCollection;

            Title = "MetaDraw: version " + GlobalVariables.MetaMorpheusVersion;

            CommonParameters = new CommonParameters();
            productMassToleranceComboBox.Items.Add("Da");
            productMassToleranceComboBox.Items.Add("ppm");
            UpdateFieldsFromPanel();
        }

        private void Window_Drop(object sender, DragEventArgs e)
        {
            if (true)
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
                        dataGridMassSpectraFiles.CommitEdit(DataGridEditingUnit.Row, true);
                        dataGridMassSpectraFiles.Items.Refresh();
                    }
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
                foreach (var rawDataFromSelected in openFileDialog1.FileNames.OrderBy(p => p))
                {
                    AddAFile(rawDataFromSelected);
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
                foreach (var rawDataFromSelected in openFileDialog1.FileNames.OrderBy(p => p))
                {
                    AddAFile(rawDataFromSelected);
                }
            dataGridResultFiles.Items.Refresh();
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
                case ".mzml":
                    RawDataForDataGrid zz = new RawDataForDataGrid(draggedFilePath);
                    if (!SpectraFileExists(spectraFilesObservableCollection, zz)) { spectraFilesObservableCollection.Add(zz); }
                    break;
                case ".pep.XML":
                case ".pep.xml":
                    break;
                case ".psmtsv":
                case ".tsv":
                    RawDataForDataGrid resultFileDataGrid = new RawDataForDataGrid(draggedFilePath);
                    if (!SpectraFileExists(resultFilesObservableCollection, resultFileDataGrid)) { resultFilesObservableCollection.Add(resultFileDataGrid); }
                    break;
                default:
                    break;
            }
        }

        private bool SpectraFileExists(ObservableCollection<RawDataForDataGrid> rDOC, RawDataForDataGrid zzz)
        {
            foreach (RawDataForDataGrid rdoc in rDOC)
                if (rdoc.FileName == zzz.FileName) { return true; }
            return false;
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

            int x = Convert.ToInt32(txtScanNum.Text);

            UpdateModel(x);

        }

        private void btnReadResultFile_Click(object sender, RoutedEventArgs e)
        {
            btnReset.IsEnabled = true;

            if (!spectraFilesObservableCollection.Any())
            {
                return;
            }
            
            LoadScans loadScans = new LoadScans(spectraFilesObservableCollection.Where(b => b.Use).First().FilePath,null);

            MsDataFile = loadScans.Run();
            
            btnReadResultFile.IsEnabled = true;

            if (resultFilesObservableCollection.Count == 0)
            {
                MessageBox.Show("Please add result files.");
                return;
            }
            var resultFilePath = resultFilesObservableCollection.Where(b => b.Use).First().FilePath;
            PSMs = TsvResultReader.ReadTsv(resultFilePath);
            foreach (var item in PSMs)
            {
                spectrumNumsObservableCollection.Add(new SpectrumForDataGrid(item.ScanNumber));
            }
            dataGridScanNums.Items.Refresh();
            

            btnReadResultFile.IsEnabled = false;
            btnDraw.IsEnabled = true;
        }

        private void Row_DoubleClick(object sender, System.Windows.Input.MouseButtonEventArgs e)
        {

            if (sender != null)
            {
                try
                {
                    Regex regex = new Regex(@"\d+");
                    int x = Convert.ToInt32(regex.Match(sender.ToString()).Value);
                    UpdateModel(x);

                }
                catch (Exception)
                {
                    MessageBox.Show("Please check the data loaded.");
                }
            }
        }

        private void UpdateModel(int x)
        {
            if (MsDataFile == null)
            {
                MessageBox.Show("Please check the MS data loaded.");
                return;
            }

            var msScanForDraw = MsDataFile.GetAllScansList().Where(p => p.OneBasedScanNumber == x).First();

            PsmDraw psmDraw = PSMs.Where(p => p.ScanNumber == x).First();

            var lp = new List<ProductType>();
            if (CommonParameters.BIons)
            {
                lp.Add(ProductType.BnoB1ions);
            }
            if (CommonParameters.YIons)
            {
                lp.Add(ProductType.Y);
            }
            if (CommonParameters.CIons)
            {
                lp.Add(ProductType.C);
            }
            if (CommonParameters.ZdotIons)
            {
                lp.Add(ProductType.Zdot);
            }

            var pmm = PsmDraw.XlCalculateTotalProductMassesForSingle(psmDraw, lp, false);

            var matchedIonMassesListPositiveIsMatch = new MatchedIonInfo(pmm.ProductMz.Length);

            double pmmScore = PsmCross.XlMatchIons(msScanForDraw, CommonParameters.ProductMassTolerance, pmm.ProductMz, pmm.ProductName, matchedIonMassesListPositiveIsMatch);

            psmDraw.MatchedIonInfo = matchedIonMassesListPositiveIsMatch;

            mainViewModel.UpdateForSingle(msScanForDraw, psmDraw);

        }

        private void btnReset_Click(object sender, RoutedEventArgs e)
        {
            btnReset.IsEnabled = false;

            resultFilesObservableCollection.Clear();

            spectraFilesObservableCollection.Clear();

            spectrumNumsObservableCollection.Clear();

            btnReadResultFile.IsEnabled = true;

            mainViewModel = new MainViewModel();
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
            if(tb.Text.Equals(string.Empty))
                tb.Text = "Scan Number";
        }

        private void bCheckBox_Checked(object sender, RoutedEventArgs e)
        {
            CommonParameters.BIons = true;
        }

        private void yCheckBox_Checked(object sender, RoutedEventArgs e)
        {
            CommonParameters.YIons = true;
        }

        private void cCheckBox_Checked(object sender, RoutedEventArgs e)
        {
            CommonParameters.CIons = true;
        }

        private void zdotCheckBox_Checked(object sender, RoutedEventArgs e)
        {
            CommonParameters.ZdotIons = true;
        }

        private void productMassToleranceTextBox_TextChanged(object sender, TextChangedEventArgs e)
        {
            if (productMassToleranceComboBox.SelectedIndex == 0)
            {
                CommonParameters.ProductMassTolerance = new AbsoluteTolerance(double.Parse(productMassToleranceTextBox.Text, CultureInfo.InvariantCulture));
            }
            else
            {
                CommonParameters.ProductMassTolerance = new PpmTolerance(double.Parse(productMassToleranceTextBox.Text, CultureInfo.InvariantCulture));
            }
        }

        private void UpdateFieldsFromPanel()
        {
            bCheckBox.IsChecked = CommonParameters.BIons;
            yCheckBox.IsChecked = CommonParameters.YIons;
            cCheckBox.IsChecked = CommonParameters.CIons;
            zdotCheckBox.IsChecked = CommonParameters.ZdotIons;
            productMassToleranceTextBox.Text = CommonParameters.ProductMassTolerance.Value.ToString(CultureInfo.InvariantCulture);
            productMassToleranceComboBox.SelectedIndex = CommonParameters.ProductMassTolerance is AbsoluteTolerance ? 0 : 1;
        }
    }
}
