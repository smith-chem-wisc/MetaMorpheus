using EngineLayer;
using MassSpectrometry;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Data;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Media;
using TaskLayer;
using ViewModels;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for MetaDraw.xaml
    /// </summary>
    public partial class MetaDraw : Window
    {
        private ItemsControlSampleViewModel itemsControlSampleViewModel;
        private PsmAnnotationViewModel mainViewModel;
        private MyFileManager spectraFileManager;
        private MsDataFile MsDataFile;
        private readonly ObservableCollection<PsmFromTsv> peptideSpectralMatches;
        ICollectionView peptideSpectralMatchesView;
        private readonly DataTable propertyView;
        private string spectraFilePath;
        private string tsvResultsFilePath;
        private Dictionary<ProductType, double> productTypeToYOffset;
        private Dictionary<ProductType, Color> productTypeToColor;
        private SolidColorBrush modificationAnnotationColor;
        private Regex illegalInFileName = new Regex(@"[\\/:*?""<>|]");

        public MetaDraw()
        {
            InitializeComponent();

            itemsControlSampleViewModel = new ItemsControlSampleViewModel();
            DataContext = itemsControlSampleViewModel;
            mainViewModel = new PsmAnnotationViewModel();
            plotView.DataContext = mainViewModel;
            peptideSpectralMatches = new ObservableCollection<PsmFromTsv>();
            propertyView = new DataTable();
            propertyView.Columns.Add("Name", typeof(string));
            propertyView.Columns.Add("Value", typeof(string));
            peptideSpectralMatchesView = CollectionViewSource.GetDefaultView(peptideSpectralMatches);
            dataGridScanNums.DataContext = peptideSpectralMatchesView;
            dataGridProperties.DataContext = propertyView.DefaultView;
            Title = "MetaDraw: version " + GlobalVariables.MetaMorpheusVersion;
            spectraFileManager = new MyFileManager(true);
            SetUpDictionaries();
            modificationAnnotationColor = Brushes.Yellow;

            ParentChildScanView.Visibility = Visibility.Collapsed;
            ParentScanView.Visibility = Visibility.Collapsed;
        }

        private void SetUpDictionaries()
        {
            // colors of each fragment to annotate on base sequence
            productTypeToColor = ((ProductType[])Enum.GetValues(typeof(ProductType))).ToDictionary(p => p, p => Colors.Aqua);
            productTypeToColor[ProductType.b] = Colors.Blue;
            productTypeToColor[ProductType.y] = Colors.Purple;
            productTypeToColor[ProductType.zDot] = Colors.Orange;
            productTypeToColor[ProductType.c] = Colors.Gold;

            // offset for annotation on base sequence
            productTypeToYOffset = ((ProductType[])Enum.GetValues(typeof(ProductType))).ToDictionary(p => p, p => 0.0);
            productTypeToYOffset[ProductType.b] = 50;
            productTypeToYOffset[ProductType.y] = 0;
            productTypeToYOffset[ProductType.c] = 50;
            productTypeToYOffset[ProductType.zDot] = 0;
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
                case ".mgf":
                    spectraFilePath = filePath;
                    spectraFileNameLabel.Text = filePath;
                    break;
                case ".psmtsv":
                case ".tsv":
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

            try
            {
                // TODO: print warnings
                foreach (var psm in PsmTsvReader.ReadTsv(filename, out List<string> warnings))
                {
                    if (psm.Filename == fileNameWithExtension || psm.Filename == fileNameWithoutExtension || psm.Filename.Contains(fileNameWithoutExtension))
                    {
                        Dispatcher.BeginInvoke(new Action(() =>
                        {
                            peptideSpectralMatches.Add(psm);
                        }));
                    }
                }
            }
            catch (Exception e)
            {
                MessageBox.Show("Could not open PSM file:\n" + e.Message);
            }
        }

        private void DrawPsm(int oneBasedScanNumber, string fullSequence = null, string fileName = null)
        {
            MsDataScan msDataScanToDraw = MsDataFile.GetOneBasedScan(oneBasedScanNumber);
            IEnumerable<PsmFromTsv> scanPsms = peptideSpectralMatches.Where(p => p.Ms2ScanNumber == oneBasedScanNumber);

            if (fullSequence != null)
            {
                scanPsms = scanPsms.Where(p => p.FullSequence == fullSequence);
            }

            PsmFromTsv psmToDraw = scanPsms.FirstOrDefault();

            // if this spectrum has child scans, draw them in the "advanced" tab
            if ((psmToDraw.ChildScanMatchedIons != null && psmToDraw.ChildScanMatchedIons.Count > 0)
                || (psmToDraw.BetaPeptideChildScanMatchedIons != null && psmToDraw.BetaPeptideChildScanMatchedIons.Count > 0))
            {
                ParentChildScanView.Visibility = Visibility.Visible;
                ParentScanView.Visibility = Visibility.Visible;

                // draw parent scans
                var parentPsmModel = new PsmAnnotationViewModel();
                MsDataScan parentScan = MsDataFile.GetOneBasedScan(psmToDraw.Ms2ScanNumber);

                parentPsmModel.DrawPeptideSpectralMatch(parentScan, psmToDraw);

                string parentAnnotation = "Scan: " + parentScan.OneBasedScanNumber.ToString()
                        + " Dissociation Type: " + parentScan.DissociationType.ToString()
                        + " MsOrder: " + parentScan.MsnOrder.ToString()
                        + " Selected Mz: " + parentScan.SelectedIonMZ.Value.ToString("0.##")
                        + " Retention Time: " + parentScan.RetentionTime.ToString("0.##");

                itemsControlSampleViewModel.AddNewRow(parentPsmModel, parentAnnotation);

                // draw child scans
                HashSet<int> scansDrawn = new HashSet<int>();
                foreach (var childScanMatchedIons in psmToDraw.ChildScanMatchedIons.Concat(psmToDraw.BetaPeptideChildScanMatchedIons))
                {
                    int scanNumber = childScanMatchedIons.Key;

                    if (scansDrawn.Contains(scanNumber))
                    {
                        continue;
                    }
                    scansDrawn.Add(scanNumber);

                    List<MatchedFragmentIon> matchedIons = childScanMatchedIons.Value;

                    var childPsmModel = new PsmAnnotationViewModel();
                    MsDataScan childScan = MsDataFile.GetOneBasedScan(scanNumber);

                    childPsmModel.DrawPeptideSpectralMatch(childScan, psmToDraw);

                    string childAnnotation = "Scan: " + scanNumber.ToString()
                        + " Dissociation Type: " + childScan.DissociationType.ToString()
                        + " MsOrder: " + childScan.MsnOrder.ToString()
                        + " Selected Mz: " + childScan.SelectedIonMZ.Value.ToString("0.##")
                        + " RetentionTime: " + childScan.RetentionTime.ToString("0.##");

                    itemsControlSampleViewModel.AddNewRow(childPsmModel, childAnnotation);
                }
            }
            else
            {
                ParentChildScanView.Visibility = Visibility.Collapsed;
                ParentScanView.Visibility = Visibility.Collapsed;
            }

            // if this is a crosslink spectrum match, there are two base sequence annotations to draw
            // this makes the canvas taller to fit both of these peptide sequences
            if (psmToDraw.BetaPeptideBaseSequence != null)
            {
                int height = 150;
                
                canvas.Height = height;
                PsmAnnotationGrid.RowDefinitions[1].Height = new GridLength(height);
            }
            else
            {
                int height = 60;
                
                canvas.Height = height;
                PsmAnnotationGrid.RowDefinitions[1].Height = new GridLength(height);
            }

            // draw annotated spectrum
            mainViewModel.DrawPeptideSpectralMatch(msDataScanToDraw, psmToDraw);

            // draw annotated base sequence
            DrawAnnotatedBaseSequence(psmToDraw);
        }

        /// <summary>
        /// Event triggers when a different cell is selected in the PSM data grid
        /// </summary>
        private void dataGridScanNums_SelectedCellsChanged(object sender, SelectedCellsChangedEventArgs e)
        {
            if (dataGridScanNums.SelectedItem == null)
            {
                return;
            }

            // draw the selected PSM
            propertyView.Clear();
            PsmFromTsv row = (PsmFromTsv)dataGridScanNums.SelectedItem;
            System.Reflection.PropertyInfo[] temp = row.GetType().GetProperties();

            for (int i = 0; i < temp.Length; i++)
            {
                if (temp[i].Name == nameof(row.MatchedIons))
                {
                    propertyView.Rows.Add(temp[i].Name, string.Join(", ", row.MatchedIons.Select(p => p.Annotation)));
                }
                else
                {
                    propertyView.Rows.Add(temp[i].Name, temp[i].GetValue(row, null));
                }
            }
            dataGridProperties.Items.Refresh();
            itemsControlSampleViewModel.Data.Clear();
            DrawPsm(row.Ms2ScanNumber, row.FullSequence);
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
            peptideSpectralMatches.Clear();

            // check for validity
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
            prgsText.Content = "Loading spectra file...";

            var slowProcess = Task<MsDataFile>.Factory.StartNew(() => spectraFileManager.LoadFile(spectraFilePath, new CommonParameters()));
            await slowProcess;
            MsDataFile = slowProcess.Result;

            // load the PSMs
            this.prgsText.Content = "Loading PSMs...";
            await Task.Run(() => LoadPsms(tsvResultsFilePath));

            // done loading - restore controls
            this.prgsFeed.IsOpen = false;
            (sender as Button).IsEnabled = true;
            selectSpectraFileButton.IsEnabled = true;
            selectPsmFileButton.IsEnabled = true;
        }

        private void TextBox_TextChanged(object sender, TextChangedEventArgs e)
        {
            string txt = (sender as TextBox).Text;
            if (txt == "")
            {
                peptideSpectralMatchesView.Filter = null;
            }
            else
            {
                peptideSpectralMatchesView.Filter = obj =>
                {
                    PsmFromTsv psm = obj as PsmFromTsv;
                    return ((psm.Ms2ScanNumber.ToString()).StartsWith(txt) || psm.FullSequence.ToUpper().Contains(txt.ToUpper()));
                };
            }
        }

        private void DrawAnnotatedBaseSequence(PsmFromTsv psm)
        {
            double spacing = 22;
            BaseDraw.clearCanvas(canvas);

            // don't draw ambiguous sequences
            if (psm.FullSequence.Contains("|"))
            {
                return;
            }

            // draw base sequence
            for (int r = 0; r < psm.BaseSeq.Length; r++)
            {
                BaseDraw.txtDrawing(canvas, new Point(r * spacing + 10, 10), psm.BaseSeq[r].ToString(), Brushes.Black);
            }

            // draw the fragment ion annotations on the base sequence
            foreach (var ion in psm.MatchedIons)
            {
                int residue = ion.NeutralTheoreticalProduct.TerminusFragment.AminoAcidPosition;
                string annotation = ion.NeutralTheoreticalProduct.ProductType + "" + ion.NeutralTheoreticalProduct.TerminusFragment.FragmentNumber;

                if (ion.NeutralTheoreticalProduct.NeutralLoss != 0)
                {
                    annotation += "-" + ion.NeutralTheoreticalProduct.NeutralLoss;
                }

                if (ion.NeutralTheoreticalProduct.TerminusFragment.Terminus == FragmentationTerminus.C)
                {
                    BaseDraw.topSplittingDrawing(canvas, new Point(residue * spacing + 8,
                        productTypeToYOffset[ion.NeutralTheoreticalProduct.ProductType]), productTypeToColor[ion.NeutralTheoreticalProduct.ProductType], annotation);
                }
                else if (ion.NeutralTheoreticalProduct.TerminusFragment.Terminus == FragmentationTerminus.N)
                {
                    BaseDraw.botSplittingDrawing(canvas, new Point(residue * spacing + 8,
                        productTypeToYOffset[ion.NeutralTheoreticalProduct.ProductType]), productTypeToColor[ion.NeutralTheoreticalProduct.ProductType], annotation);
                }
                // don't draw diagnostic ions, precursor ions, etc
            }

            // draw modifications
            var peptide = new PeptideWithSetModifications(psm.FullSequence, GlobalVariables.AllModsKnownDictionary);
            foreach (var mod in peptide.AllModsOneIsNterminus)
            {
                BaseDraw.circledTxtDraw(canvas, new Point((mod.Key - 1) * spacing - 17, 12), modificationAnnotationColor);
            }

            if (psm.BetaPeptideBaseSequence != null)
            {
                for (int r = 0; r < psm.BetaPeptideBaseSequence.Length; r++)
                {
                    BaseDraw.txtDrawing(canvas, new Point(r * spacing + 10, 100), psm.BetaPeptideBaseSequence[r].ToString(), Brushes.Black);
                }

                foreach (var ion in psm.BetaPeptideMatchedIons)
                {
                    int residue = ion.NeutralTheoreticalProduct.TerminusFragment.AminoAcidPosition;
                    string annotation = ion.NeutralTheoreticalProduct.ProductType + "" + ion.NeutralTheoreticalProduct.TerminusFragment.FragmentNumber;

                    if (ion.NeutralTheoreticalProduct.NeutralLoss != 0)
                    {
                        annotation += "-" + ion.NeutralTheoreticalProduct.NeutralLoss;
                    }

                    if (ion.NeutralTheoreticalProduct.TerminusFragment.Terminus == FragmentationTerminus.C)
                    {
                        BaseDraw.topSplittingDrawing(canvas, new Point(residue * spacing + 8,
                            productTypeToYOffset[ion.NeutralTheoreticalProduct.ProductType] + 90), productTypeToColor[ion.NeutralTheoreticalProduct.ProductType], annotation);
                    }
                    else if (ion.NeutralTheoreticalProduct.TerminusFragment.Terminus == FragmentationTerminus.N)
                    {
                        BaseDraw.botSplittingDrawing(canvas, new Point(residue * spacing + 8,
                            productTypeToYOffset[ion.NeutralTheoreticalProduct.ProductType] + 90), productTypeToColor[ion.NeutralTheoreticalProduct.ProductType], annotation);
                    }
                    // don't draw diagnostic ions, precursor ions, etc
                }

                var betaPeptide = new PeptideWithSetModifications(psm.BetaPeptideFullSequence, GlobalVariables.AllModsKnownDictionary);
                foreach (var mod in betaPeptide.AllModsOneIsNterminus)
                {
                    BaseDraw.circledTxtDraw(canvas, new Point((mod.Key - 1) * spacing - 17, 12 + 90), modificationAnnotationColor);
                }

                int alphaSite = Int32.Parse(Regex.Match(psm.FullSequence, @"\d+").Value);
                int betaSite = Int32.Parse(Regex.Match(psm.BetaPeptideFullSequence, @"\d+").Value);
                //BaseDraw.linkDrawing(canvas, new Point(alphaSite * spacing, 50), new Point(betaSite * spacing, 90), Colors.Black);
            }
        }

        private void dataGridProperties_SelectedCellsChanged(object sender, SelectedCellsChangedEventArgs e)
        {
            (sender as DataGrid).UnselectAll();
        }

        private void PDFButton_Click(object sender, RoutedEventArgs e)
        {
            PsmFromTsv tempPsm = null;
            if (dataGridScanNums.SelectedCells.Count == 0)
            {
                MessageBox.Show("Please select at least one scan to export");
            }

            int num = dataGridScanNums.SelectedItems.Count;

            foreach (object selectedItem in dataGridScanNums.SelectedItems)
            {
                PsmFromTsv psm = (PsmFromTsv)selectedItem;

                if (tempPsm == null)
                {
                    tempPsm = psm;
                }

                MsDataScan msDataScanToDraw = MsDataFile.GetOneBasedScan(psm.Ms2ScanNumber);

                string myString = illegalInFileName.Replace(psm.FullSequence, "");

                if (myString.Length > 30)
                {
                    myString = myString.Substring(0, 30);
                }

                string filePath = Path.Combine(Path.GetDirectoryName(tsvResultsFilePath), "MetaDrawExport", psm.Ms2ScanNumber + "_" + myString + ".pdf");

                mainViewModel.DrawPeptideSpectralMatchPdf(msDataScanToDraw, psm, filePath, num > 1);
            }

            dataGridScanNums.SelectedItem = dataGridScanNums.SelectedItem;

            DrawPsm(tempPsm.Ms2ScanNumber, tempPsm.FullSequence);

            MessageBox.Show(string.Format("{0} PDFs exported", num));
        }

        private void BtnChangeGridColumns_Click(object sender, RoutedEventArgs e)
        {
            itemsControlSampleViewModel.MyColumnCount++;
            if (itemsControlSampleViewModel.MyColumnCount > itemsControlSampleViewModel.Data.Count / 3)
            {
                itemsControlSampleViewModel.MyColumnCount = 1;
            }
        }
    }
}