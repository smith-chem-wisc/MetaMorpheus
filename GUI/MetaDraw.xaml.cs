using EngineLayer;
using MassSpectrometry;
using OxyPlot;
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
using System.Windows.Media.Imaging;
using TaskLayer;
using ViewModels;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for MetaDraw.xaml
    /// </summary>
    public partial class MetaDraw : Window
    {
        private MetaDrawGraphicalSettings metaDrawGraphicalSettings;
        private MetaDrawFilterSettings metaDrawFilterSettings;
        private ItemsControlSampleViewModel itemsControlSampleViewModel;
        private PsmAnnotationViewModel mainViewModel;
        private MyFileManager spectraFileManager;
        private MsDataFile MsDataFile;
        private readonly ObservableCollection<PsmFromTsv> allPsms; // all loaded PSMs
        private readonly ObservableCollection<PsmFromTsv> filteredListOfPsms; // this is the filtered list of PSMs to display (after q-value filter, etc.)
        private readonly Dictionary<string, ObservableCollection<PsmFromTsv>> psmsBySourceFile;  // filtered PSMs grouped by the "file name" field
        private readonly ObservableCollection<string> sourceFilesList;
        ICollectionView peptideSpectralMatchesView;
        private readonly DataTable propertyView;
        private string spectraFilePath;
        private string tsvResultsFilePath;
        private Dictionary<ProductType, double> productTypeToYOffset;
        private Dictionary<ProductType, Color> productTypeToColor;
        private Color variantCrossColor;
        private SolidColorBrush modificationAnnotationColor;
        private Regex illegalInFileName = new Regex(@"[\\/:*?""<>|]");
        private ObservableCollection<string> plotTypes;

        public MetaDraw()
        {
            UsefulProteomicsDatabases.Loaders.LoadElements();

            InitializeComponent();

            itemsControlSampleViewModel = new ItemsControlSampleViewModel();
            DataContext = itemsControlSampleViewModel;
            mainViewModel = new PsmAnnotationViewModel();
            plotView.DataContext = mainViewModel;
            allPsms = new ObservableCollection<PsmFromTsv>();
            filteredListOfPsms = new ObservableCollection<PsmFromTsv>();
            psmsBySourceFile = new Dictionary<string, ObservableCollection<PsmFromTsv>>();
            sourceFilesList = new ObservableCollection<string>();
            selectSourceFileListBox.DataContext = sourceFilesList;
            propertyView = new DataTable();
            propertyView.Columns.Add("Name", typeof(string));
            propertyView.Columns.Add("Value", typeof(string));
            peptideSpectralMatchesView = CollectionViewSource.GetDefaultView(filteredListOfPsms);
            dataGridScanNums.DataContext = peptideSpectralMatchesView;
            dataGridProperties.DataContext = propertyView.DefaultView;
            Title = "MetaDraw: version " + GlobalVariables.MetaMorpheusVersion;
            spectraFileManager = new MyFileManager(true);
            SetUpDictionaries();
            variantCrossColor = Colors.Green;
            modificationAnnotationColor = Brushes.Orange;
            metaDrawGraphicalSettings = new MetaDrawGraphicalSettings();
            metaDrawFilterSettings = new MetaDrawFilterSettings();
            base.Closing += this.OnClosing;

            ParentChildScanView.Visibility = Visibility.Collapsed;
            ParentScanView.Visibility = Visibility.Collapsed;

            plotTypes = new ObservableCollection<string>();
            SetUpPlots();
            plotsListBox.ItemsSource = plotTypes;
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
            productTypeToYOffset[ProductType.c] = 53.6;
            productTypeToYOffset[ProductType.zDot] = -3.6;
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
                    if (spectraFilePath == null || !spectraFilePath.Equals(filePath))
                    {
                        resetFilesButton_Click(resetSpectraFileButton, null);
                        spectraFilePath = filePath;
                        spectraFileNameLabel.Text = filePath;
                        spectraFileNameLabel.ToolTip = filePath;
                        resetSpectraFileButton.IsEnabled = true;
                    }
                    break;               
                case ".tsv":
                case ".psmtsv":
                    if (tsvResultsFilePath == null || !tsvResultsFilePath.Equals(filePath))
                    {
                        resetFilesButton_Click(resetPsmFileButton, null);
                        tsvResultsFilePath = filePath;
                        psmFileNameLabel.Text = filePath;
                        psmFileNameLabel.ToolTip = filePath;
                        psmFileNameLabelStat.Text = filePath;
                        psmFileNameLabelStat.ToolTip = filePath;
                        resetPsmFileButton.IsEnabled = true;
                        resetPsmFileButtonStat.IsEnabled = true;
                    }
                    break;
                default:
                    MessageBox.Show("Cannot read file type: " + theExtension);
                    break;
            }
        }

        private void LoadPsms(string filename, bool loadPsmsOfAllSpectraFiles)
        {
            allPsms.Clear();

            string fileNameWithExtension = Path.GetFileName(spectraFilePath);
            string fileNameWithoutExtension = Path.GetFileNameWithoutExtension(spectraFilePath);

            try
            {
                // TODO: print warnings
                foreach (var psm in PsmTsvReader.ReadTsv(filename, out List<string> warnings))
                {
                    if (loadPsmsOfAllSpectraFiles || Path.GetFileName(psm.Filename) == fileNameWithExtension || psm.Filename == fileNameWithExtension || psm.Filename == fileNameWithoutExtension || psm.Filename.Equals(fileNameWithoutExtension + "-calib")) // in case results are from a calibrated file
                    {
                        allPsms.Add(psm);
                    }
                }
            }
            catch (Exception e)
            {
                MessageBox.Show("Could not open PSM file:\n" + e.Message);
            }
        }

        private void DisplayLoadedAndFilteredPsms()
        {
            filteredListOfPsms.Clear();

            var filteredList = allPsms.Where(p =>
                p.QValue <= metaDrawFilterSettings.QValueFilter
                && (p.QValueNotch < metaDrawFilterSettings.QValueFilter || p.QValueNotch == null)
                && (p.DecoyContamTarget == "T" || (p.DecoyContamTarget == "D" && metaDrawFilterSettings.ShowDecoys) || (p.DecoyContamTarget == "C" && metaDrawFilterSettings.ShowContaminants))
                && (p.GlycanLocalizationLevel == null || p.GlycanLocalizationLevel >= metaDrawFilterSettings.LocalizationLevelStart && p.GlycanLocalizationLevel <= metaDrawFilterSettings.LocalizationLevelEnd)
                );

            foreach (PsmFromTsv psm in filteredList)
            {
                filteredListOfPsms.Add(psm);
            }
        }

        private void GroupPsmsBySourceFile()
        {
            psmsBySourceFile.Clear();
            sourceFilesList.Clear();

            foreach (PsmFromTsv psm in filteredListOfPsms)
            {
                if (psmsBySourceFile.Keys.Contains(psm.Filename))
                {
                    psmsBySourceFile[psm.Filename].Add(psm);
                }
                else
                {
                    psmsBySourceFile.Add(psm.Filename, new ObservableCollection<PsmFromTsv> { psm });
                    sourceFilesList.Add(psm.Filename);
                }
            }
        }

        private void DrawPsm(int oneBasedScanNumber, string fullSequence = null, string fileName = null)
        {
            MsDataScan msDataScanToDraw = MsDataFile.GetOneBasedScan(oneBasedScanNumber);
            IEnumerable<PsmFromTsv> scanPsms = filteredListOfPsms.Where(p => p.Ms2ScanNumber == oneBasedScanNumber);

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

                parentPsmModel.DrawPeptideSpectralMatch(parentScan, psmToDraw, metaDrawGraphicalSettings.ShowMzValues,
                        metaDrawGraphicalSettings.ShowAnnotationCharges, metaDrawGraphicalSettings.AnnotatedFontSize, metaDrawGraphicalSettings.BoldText);

                string parentAnnotation = "Scan: " + parentScan.OneBasedScanNumber.ToString()
                        + " Dissociation Type: " + parentScan.DissociationType.ToString()
                        + " MsOrder: " + parentScan.MsnOrder.ToString()
                        + " Selected Mz: " + parentScan.SelectedIonMZ.Value.ToString("0.##")
                        + " Retention Time: " + parentScan.RetentionTime.ToString("0.##");

                itemsControlSampleViewModel.AddNewRow(parentPsmModel, parentAnnotation, null);

                // draw child scans
                HashSet<int> scansDrawn = new HashSet<int>();
                var allChildScanMatchedIons = psmToDraw.ChildScanMatchedIons;
                if (psmToDraw.BetaPeptideChildScanMatchedIons != null)
                {
                    allChildScanMatchedIons = allChildScanMatchedIons.Concat(psmToDraw.BetaPeptideChildScanMatchedIons).ToDictionary(p=>p.Key, q=>q.Value);
                }
                foreach (var childScanMatchedIons in allChildScanMatchedIons)
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

                    childPsmModel.DrawPeptideSpectralMatch(childScan, psmToDraw, metaDrawGraphicalSettings.ShowMzValues,
                        metaDrawGraphicalSettings.ShowAnnotationCharges, metaDrawGraphicalSettings.AnnotatedFontSize, metaDrawGraphicalSettings.BoldText);

                    string childAnnotation = "Scan: " + scanNumber.ToString()
                        + " Dissociation Type: " + childScan.DissociationType.ToString()
                        + " MsOrder: " + childScan.MsnOrder.ToString()
                        + " Selected Mz: " + childScan.SelectedIonMZ.Value.ToString("0.##")
                        + " RetentionTime: " + childScan.RetentionTime.ToString("0.##");

                    Canvas aCanvas = new Canvas { Height = 60  };
                    DrawAnnotatedBaseSequence(aCanvas, psmToDraw, true);
          
                    itemsControlSampleViewModel.AddNewRow(childPsmModel, childAnnotation, aCanvas);
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
            mainViewModel.DrawPeptideSpectralMatch(msDataScanToDraw, psmToDraw, metaDrawGraphicalSettings.ShowMzValues,
                metaDrawGraphicalSettings.ShowAnnotationCharges, metaDrawGraphicalSettings.AnnotatedFontSize, metaDrawGraphicalSettings.BoldText);

            // draw annotated base sequence
            DrawAnnotatedBaseSequence(canvas, psmToDraw);
        }

        /// <summary>
        /// Event triggers when a different cell is selected in the PSM data grid
        /// </summary>
        private void dataGridScanNums_SelectedCellsChanged(object sender, SelectedCellsChangedEventArgs e)
        {
            OnSelectionChanged();
        }

        private void OnSelectionChanged()
        {
            if (dataGridScanNums.SelectedItem == null)
            {
                return;
            }

            if (spectraFilePath == null)
            {
                MessageBox.Show("Please add a spectra file.");
                return;
            }

            if (MsDataFile == null)
            {
                MessageBox.Show("Please load spectra file.");
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
                else if(temp[i].Name == nameof(row.VariantCrossingIons))
                {
                    propertyView.Rows.Add(temp[i].Name, string.Join(", ", row.VariantCrossingIons.Select(p => p.Annotation)));
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

        private void resetFilesButton_Click(object sender, RoutedEventArgs e)
        {
            // clear spectra file path
            if ((sender as Button).Name.Equals(resetSpectraFileButton.Name))
            {
                spectraFilePath = null;
                spectraFileNameLabel.Text = "None Selected";
                spectraFileNameLabel.ToolTip = null;
                resetSpectraFileButton.IsEnabled = false;
            }

            // clear psm file path, psm data, and stat plot
            else
            {
                tsvResultsFilePath = null;
                psmFileNameLabel.Text = "None Selected";
                psmFileNameLabel.ToolTip = null;
                psmFileNameLabelStat.Text = "None Selected";
                psmFileNameLabelStat.ToolTip = null;
                resetPsmFileButton.IsEnabled = false;
                resetPsmFileButtonStat.IsEnabled = false;

                allPsms.Clear();
                filteredListOfPsms.Clear();
                psmsBySourceFile.Clear();
                sourceFilesList.Clear();
                plotViewStat.DataContext = null;
            }

            // clear main plot, base sequence, and properties table
            mainViewModel.Model = null;
            BaseDraw.clearCanvas(canvas);
            propertyView.Clear();
            MsDataFile = null;
        }

        private void OnClosing(object sender, CancelEventArgs e)
        {
            metaDrawGraphicalSettings.Close();
            metaDrawFilterSettings.Close();
        }

        private void graphicalSettings_Click(object sender, RoutedEventArgs e)
        {
            metaDrawGraphicalSettings.MZCheckBox.IsChecked = metaDrawGraphicalSettings.ShowMzValues;
            metaDrawGraphicalSettings.ChargesCheckBox.IsChecked = metaDrawGraphicalSettings.ShowAnnotationCharges;
            metaDrawGraphicalSettings.BoldTextCheckBox.IsChecked = metaDrawGraphicalSettings.BoldText;
            metaDrawGraphicalSettings.TextSizeBox.Text = metaDrawGraphicalSettings.AnnotatedFontSize.ToString();

            metaDrawGraphicalSettings.ShowDialog();

            OnSelectionChanged();
        }

        private void filterSettings_Click(object sender, RoutedEventArgs e)
        {
            metaDrawFilterSettings.DecoysCheckBox.IsChecked = metaDrawFilterSettings.ShowDecoys;
            metaDrawFilterSettings.ContaminantsCheckBox.IsChecked = metaDrawFilterSettings.ShowContaminants;
            metaDrawFilterSettings.qValueBox.Text = metaDrawFilterSettings.QValueFilter.ToString();
            metaDrawFilterSettings.CmbGlycanLocalizationLevelStart.SelectedIndex = (int)metaDrawFilterSettings.LocalizationLevelStart;
            metaDrawFilterSettings.CmbGlycanLocalizationLevelEnd.SelectedIndex = (int)metaDrawFilterSettings.LocalizationLevelEnd;
            var result = metaDrawFilterSettings.ShowDialog();

            DisplayLoadedAndFilteredPsms();
        }

        private async void loadFilesButton_Click(object sender, RoutedEventArgs e)
        {
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
            resetSpectraFileButton.IsEnabled = false;
            resetPsmFileButton.IsEnabled = false;
            prgsFeed.IsOpen = true;
            prgsText.Content = "Loading spectra file...";

            // Add EventHandlers for popup click-in/click-out behaviour
            Deactivated += new EventHandler(prgsFeed_Deactivator);
            Activated += new EventHandler(prgsFeed_Reactivator);

            var slowProcess = Task<MsDataFile>.Factory.StartNew(() => spectraFileManager.LoadFile(spectraFilePath, new CommonParameters(trimMsMsPeaks: false)));
            await slowProcess;
            MsDataFile = slowProcess.Result;

            // load the PSMs
            this.prgsText.Content = "Loading PSMs...";
            LoadPsms(tsvResultsFilePath, false);
            DisplayLoadedAndFilteredPsms();

            // done loading - restore controls
            this.prgsFeed.IsOpen = false;

            // Remove added EventHandlers
            Deactivated -= new EventHandler(prgsFeed_Deactivator);
            Activated -= new EventHandler(prgsFeed_Reactivator);

            (sender as Button).IsEnabled = true;
            selectSpectraFileButton.IsEnabled = true;
            selectPsmFileButton.IsEnabled = true;
            resetSpectraFileButton.IsEnabled = true;
            resetPsmFileButton.IsEnabled = true;
        }

        /// <summary>
        /// Deactivates the popup if one clicks out of the main window
        /// </summary>
        /// <param name="sender">unused</param>
        /// <param name="e">unused</param>
        private void prgsFeed_Deactivator(object sender, EventArgs e)
        {
            prgsFeed.IsOpen = false;
        }

        /// <summary>
        /// Reacctivates the popup if one clicks into the main window
        /// </summary>
        /// <param name="sender">unused</param>
        /// <param name="e">unused</param>
        private void prgsFeed_Reactivator(object sender, EventArgs e)
        {
            prgsFeed.IsOpen = true;
        }   
        
        private void loadFilesButtonStat_Click(object sender, RoutedEventArgs e)
        {
            // check for validity
            if (tsvResultsFilePath == null)
            {
                MessageBox.Show("Please add a search result file.");
                return;
            }

            (sender as Button).IsEnabled = false;
            selectPsmFileButtonStat.IsEnabled = false;
            resetPsmFileButtonStat.IsEnabled = false;
            prgsFeedStat.IsOpen = true;

            // load the PSMs
            this.prgsTextStat.Content = "Loading PSMs...";
            LoadPsmsStat(tsvResultsFilePath);

            // done loading - restore controls
            this.prgsFeedStat.IsOpen = false;
            (sender as Button).IsEnabled = true;
            selectPsmFileButtonStat.IsEnabled = true;
            resetPsmFileButtonStat.IsEnabled = true;
        }

        private void LoadPsmsStat(string filepath)
        {

            LoadPsms(filepath, true);
            DisplayLoadedAndFilteredPsms();
            GroupPsmsBySourceFile();
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

        private void DrawAnnotatedBaseSequence(Canvas canvas, PsmFromTsv psm, bool ToAnnotateChildScan = false)
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

            var matchIons = psm.MatchedIons;

            if (ToAnnotateChildScan)
            {
                matchIons = new List<MatchedFragmentIon>();
                foreach (var x in psm.ChildScanMatchedIons)
                {
                    matchIons.AddRange(x.Value);
                }
            }
            // draw the fragment ion annotations on the base sequence
            foreach (var ion in matchIons)
            {
                int residue = ion.NeutralTheoreticalProduct.AminoAcidPosition;
                string annotation = ion.NeutralTheoreticalProduct.ProductType + "" + ion.NeutralTheoreticalProduct.FragmentNumber;
                Color color = psm.VariantCrossingIons.Contains(ion) ? variantCrossColor : productTypeToColor[ion.NeutralTheoreticalProduct.ProductType];

                if (ion.NeutralTheoreticalProduct.NeutralLoss != 0)
                {
                    annotation += "-" + ion.NeutralTheoreticalProduct.NeutralLoss;
                }

                if (ion.NeutralTheoreticalProduct.Terminus == FragmentationTerminus.C)
                {
                    BaseDraw.topSplittingDrawing(canvas, new Point(residue * spacing + 8,
                        productTypeToYOffset[ion.NeutralTheoreticalProduct.ProductType]), color, annotation);
                }
                else if (ion.NeutralTheoreticalProduct.Terminus == FragmentationTerminus.N)
                {
                    BaseDraw.botSplittingDrawing(canvas, new Point(residue * spacing + 8,
                        productTypeToYOffset[ion.NeutralTheoreticalProduct.ProductType]), color, annotation);
                }
                // don't draw diagnostic ions, precursor ions, etc
            }

            // draw modifications
            var peptide = new PeptideWithSetModifications(psm.FullSequence, GlobalVariables.AllModsKnownDictionary);
            if (psm.GlycanLocalizationLevel!=null)
            {
                var localGlycans = PsmFromTsv.ReadLocalizedGlycan(psm.LocalizedGlycan);

                foreach (var mod in peptide.AllModsOneIsNterminus)
                {
                    if (mod.Value.ModificationType == "O-Glycosylation")
                    {
                        if (localGlycans.Where(p=>p.Item1 + 1 == mod.Key).Count() > 0)
                        {
                            BaseDraw.circledTxtDraw(canvas, new Point((mod.Key - 1) * spacing - 17, 12), modificationAnnotationColor);
                        }
                        else
                        {
                            BaseDraw.circledTxtDraw(canvas, new Point((mod.Key - 1) * spacing - 17, 12), Brushes.Gray);
                        }
                    }

                   
                }
            }
            else
            {
                foreach (var mod in peptide.AllModsOneIsNterminus)
                {
                    BaseDraw.circledTxtDraw(canvas, new Point((mod.Key - 1) * spacing - 17, 12), modificationAnnotationColor);
                }
            }


            if (psm.BetaPeptideBaseSequence != null)
            {
                for (int r = 0; r < psm.BetaPeptideBaseSequence.Length; r++)
                {
                    BaseDraw.txtDrawing(canvas, new Point(r * spacing + 10, 100), psm.BetaPeptideBaseSequence[r].ToString(), Brushes.Black);
                }

                foreach (var ion in psm.BetaPeptideMatchedIons)
                {
                    int residue = ion.NeutralTheoreticalProduct.AminoAcidPosition;
                    string annotation = ion.NeutralTheoreticalProduct.ProductType + "" + ion.NeutralTheoreticalProduct.FragmentNumber;

                    if (ion.NeutralTheoreticalProduct.NeutralLoss != 0)
                    {
                        annotation += "-" + ion.NeutralTheoreticalProduct.NeutralLoss;
                    }

                    if (ion.NeutralTheoreticalProduct.Terminus == FragmentationTerminus.C)
                    {
                        BaseDraw.topSplittingDrawing(canvas, new Point(residue * spacing + 8,
                            productTypeToYOffset[ion.NeutralTheoreticalProduct.ProductType] + 90), productTypeToColor[ion.NeutralTheoreticalProduct.ProductType], annotation);
                    }
                    else if (ion.NeutralTheoreticalProduct.Terminus == FragmentationTerminus.N)
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
                BaseDraw.DrawCrosslinker(canvas, new Point(alphaSite * spacing, 50), new Point(betaSite * spacing, 90), Colors.Black);
            }
        }

        private void SetUpPlots()
        {
            foreach (var plot in PlotModelStat.PlotNames)
            {
                plotTypes.Add(plot);
            }
        }

        private void dataGridProperties_SelectedCellsChanged(object sender, SelectedCellsChangedEventArgs e)
        {
            (sender as DataGrid).UnselectAll();
        }

        private void CreatePlotPdf_Click(object sender, RoutedEventArgs e)
        {
            var selectedItem = plotsListBox.SelectedItem;

            if (selectedItem == null)
            {
                MessageBox.Show("Select a plot type to export!");
                return;
            }

            if (!filteredListOfPsms.Any())
            {
                MessageBox.Show("No PSMs are loaded!");
                return;
            }

            if (selectSourceFileListBox.SelectedItems.Count == 0)
            {
                MessageBox.Show("Please select a source file.");
                return;
            }

            var plotName = selectedItem as string;
            var fileDirectory = Directory.GetParent(tsvResultsFilePath).ToString();
            var fileName = String.Concat(plotName, ".pdf");

            // update font sizes to exported PDF's size
            double tmpW = plotViewStat.Width;
            double tmpH = plotViewStat.Height;
            plotViewStat.Width = 1000;
            plotViewStat.Height = 700;
            plotViewStat.UpdateLayout();
            PlotViewStat_SizeChanged(plotViewStat, null);

            using (Stream writePDF = File.Create(Path.Combine(fileDirectory, fileName)))
            {
                PdfExporter.Export(plotViewStat.Model, writePDF, 1000, 700);
            }
            plotViewStat.Width = tmpW;
            plotViewStat.Height = tmpH;
            MessageBox.Show("PDF Created at " + Path.Combine(fileDirectory, fileName) + "!");
        }

        private void PDFButton_Click(object sender, RoutedEventArgs e)
        {
            PsmFromTsv tempPsm = null;
            if (dataGridScanNums.SelectedCells.Count == 0)
            {
                MessageBox.Show("Please select at least one scan to export");
            }
            else
            {
                int numberOfScansToExport = dataGridScanNums.SelectedItems.Count;

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
                    string dir = Path.GetDirectoryName(filePath);

                    if (!Directory.Exists(dir))
                    {
                        Directory.CreateDirectory(dir);
                    }

                    DrawPdfAnnotatedBaseSequence(psm, canvas, filePath); // captures the annotation for the pdf
                    mainViewModel.DrawPeptideSpectralMatchPdf(msDataScanToDraw, psm, filePath, numberOfScansToExport > 1);

                    // if this spectrum has child scans, draw them in the "advanced" tab
                    if ((psm.ChildScanMatchedIons != null && psm.ChildScanMatchedIons.Count > 0)
                        || (psm.BetaPeptideChildScanMatchedIons != null && psm.BetaPeptideChildScanMatchedIons.Count > 0))
                    {
                        // draw child scans
                        HashSet<int> scansDrawn = new HashSet<int>();
                        var allChildScanMatchedIons = psm.ChildScanMatchedIons;
                        if (psm.BetaPeptideChildScanMatchedIons != null)
                        {
                            allChildScanMatchedIons = allChildScanMatchedIons.Concat(psm.BetaPeptideChildScanMatchedIons).ToDictionary(p => p.Key, q => q.Value);
                        }
                        foreach (var childScanMatchedIons in allChildScanMatchedIons)
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

                            childPsmModel.DrawPeptideSpectralMatch(childScan, psm, metaDrawGraphicalSettings.ShowMzValues,
                                metaDrawGraphicalSettings.ShowAnnotationCharges, metaDrawGraphicalSettings.AnnotatedFontSize, metaDrawGraphicalSettings.BoldText);

                            string childfilePath = Path.Combine(Path.GetDirectoryName(tsvResultsFilePath), "MetaDrawExport", childScan.OneBasedScanNumber + "_" + myString + ".pdf");
                            
                            DrawPdfAnnotatedBaseSequence(psm, canvas, childfilePath); // captures the annotation for the pdf
                            childPsmModel.DrawPeptideSpectralMatchPdf(childScan, psm, childfilePath, false);
                        }
                    }
                }

                dataGridScanNums.SelectedItem = dataGridScanNums.SelectedItem;

                itemsControlSampleViewModel.Data.Clear();
                DrawPsm(tempPsm.Ms2ScanNumber, tempPsm.FullSequence);

                MessageBox.Show(string.Format("{0} PDFs exported", numberOfScansToExport));
            }
        }

        private void DrawPdfAnnotatedBaseSequence(PsmFromTsv psm, Canvas canvas, string path)
        {
            if (psm.CrossType == null)
            {
                DrawAnnotatedBaseSequence(canvas, psm);
            }

            canvas.Measure(new Size((int)canvas.Width, 600));
            canvas.Arrange(new Rect(new Size((int)canvas.Width, 600)));

            RenderTargetBitmap renderBitmap = new RenderTargetBitmap((int)(canvas.Width), 600, 96, 96, PixelFormats.Pbgra32);

            renderBitmap.Render(canvas);
            PngBitmapEncoder encoder = new PngBitmapEncoder();
            encoder.Frames.Add(BitmapFrame.Create(renderBitmap));

            string tempPath = Path.Combine(Path.GetDirectoryName(tsvResultsFilePath), "MetaDrawExport", "annotation.png");

            using (FileStream file = File.Create(tempPath))
            {
                encoder.Save(file);
            }
        }

        private async void PlotSelected(object sender, SelectionChangedEventArgs e)
        {
            var listview = sender as ListView;
            var plotName = listview.SelectedItem as string;

            if (filteredListOfPsms.Count == 0)
            {
                MessageBox.Show("There are no PSMs to analyze.\n\nLoad the current file or choose a new file.");
                return;
            }
            if (selectSourceFileListBox.SelectedItems.Count == 0)
            {
                MessageBox.Show("Please select a source file.");
                return;
            }

            // get psms from selected source files
            ObservableCollection<PsmFromTsv> psms = new ObservableCollection<PsmFromTsv>();
            Dictionary<string, ObservableCollection<PsmFromTsv>> psmsBSF = new Dictionary<string, ObservableCollection<PsmFromTsv>>();
            foreach (string fileName in selectSourceFileListBox.SelectedItems)
            {
                psmsBSF.Add(fileName, psmsBySourceFile[fileName]);
                foreach (PsmFromTsv psm in psmsBySourceFile[fileName])
                {
                    psms.Add(psm);
                }
            }
            PlotModelStat plot = await Task.Run(() => new PlotModelStat(plotName, psms, psmsBSF));
            plotViewStat.DataContext = plot;
            PlotViewStat_SizeChanged(plotViewStat, null);
        }

        private void BtnChangeGridColumns_Click(object sender, RoutedEventArgs e)
        {
            itemsControlSampleViewModel.MyColumnCount++;
            if (itemsControlSampleViewModel.MyColumnCount > itemsControlSampleViewModel.Data.Count / 3)
            {
                itemsControlSampleViewModel.MyColumnCount = 1;
            }
        }

        private void selectSourceFileListBox_SelectionChanged(object sender, EventArgs e)
        {
            // refreshes the plot using the new source file
            if (plotsListBox.SelectedIndex > -1 && selectSourceFileListBox.SelectedItems.Count != 0)
            {
                PlotSelected(plotsListBox, null);
            }
        }

        private void selectAllSourceFiles_Click(object sender, RoutedEventArgs e)
        {
            selectSourceFileListBox.SelectAll();
        }

        private void deselectAllSourceFiles_Click(object sender, RoutedEventArgs e)
        {
            selectSourceFileListBox.SelectedIndex = -1;
        }

        // scales the font size down for the x axis labels of the PTM histogram when the window gets too small
        private void PlotViewStat_SizeChanged(object sender, SizeChangedEventArgs e)
        {
            if(plotsListBox.SelectedItem == null || !plotsListBox.SelectedItem.ToString().Equals("Histogram of PTM Spectral Counts"))
            {
                return;
            }
            OxyPlot.Wpf.PlotView plot = sender as OxyPlot.Wpf.PlotView;
            if(plot != null && plot.Model != null)
            {
                plot.Model.DefaultXAxis.TitleFontSize = plot.Model.DefaultFontSize; // stops the title from being scaled
                int count = (int)plot.Model.DefaultXAxis.ActualMaximum;
                int widthCountRatio = 23;   // maintains this width:number of PTM types ratio
                if (plot.ActualWidth / count < widthCountRatio) 
                {
                    plot.Model.DefaultXAxis.FontSize = plot.Model.DefaultFontSize * (plot.ActualWidth / (count * widthCountRatio));
                }
                else
                {
                    plot.Model.DefaultXAxis.FontSize = plot.Model.DefaultFontSize;
                }
            }
        }
    }
}