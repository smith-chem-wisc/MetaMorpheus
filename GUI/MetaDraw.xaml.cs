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
        private MetaDrawGraphicalSettings metaDrawGraphicalSettings;
        private MetaDrawFilterSettings metaDrawFilterSettings;
        private ItemsControlSampleViewModel itemsControlSampleViewModel;
        private PsmAnnotationViewModel mainViewModel;
        private MyFileManager spectraFileManager;
        private MsDataFile MsDataFile;
        private readonly ObservableCollection<PsmFromTsv> allPsms; // all loaded PSMs
        private readonly ObservableCollection<PsmFromTsv> filteredListOfPsms; // this is the filtered list of PSMs to display (after q-value filter, etc.)
        private ObservableCollection<ProteinForTreeView> proteinTree;
        private ObservableCollection<ProteinForTreeView> filteredTree;
        //ICollectionView proteinView;
        private readonly DataTable propertyView;
        private string spectraFilePath;
        private string psmFilePath;
        private string proteinFilePath;
        private Dictionary<ProductType, double> productTypeToYOffset;
        private Dictionary<ProductType, Color> productTypeToColor;
        private SolidColorBrush modificationAnnotationColor;
        private Dictionary<string, Color> proteaseByColor;
        private Dictionary<string, string> proteinGroups;
        private Dictionary<string, ProteinForTreeView> proteinGroupsForTreeView;
        private Regex illegalInFileName = new Regex(@"[\\/:*?""<>|]");
        private ObservableCollection<string> plotTypes;
        private List<PsmFromTsv> psmsWithMatch;
        private string[] proteases;

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
            proteinTree = new ObservableCollection<ProteinForTreeView>();
            filteredTree = new ObservableCollection<ProteinForTreeView>();
            propertyView = new DataTable();
            propertyView.Columns.Add("Name", typeof(string));
            propertyView.Columns.Add("Value", typeof(string));
            psmsWithMatch = new List<PsmFromTsv>();
            proteinTreeView.DataContext = proteinTree;
            dataGridProperties.DataContext = propertyView.DefaultView;
            Title = "MetaDraw: version " + GlobalVariables.MetaMorpheusVersion;
            spectraFileManager = new MyFileManager(true);
            SetUpDictionaries();
            modificationAnnotationColor = Brushes.Yellow;
            metaDrawGraphicalSettings = new MetaDrawGraphicalSettings();
            metaDrawFilterSettings = new MetaDrawFilterSettings();
            SearchTimer.Timer.Tick += new EventHandler(searchBox_TextChangedHandler);
            base.Closing += this.OnClosing;

            ParentChildScanView.Visibility = Visibility.Collapsed;
            ParentScanView.Visibility = Visibility.Collapsed;
            mapViewer.Visibility = Visibility.Collapsed;
            legend.Visibility = Visibility.Collapsed;

            plotTypes = new ObservableCollection<string>();
            proteases = new string[1] { "trypsin" };

            SetUpPlots();
            //plotsListBox.ItemsSource = plotTypes;
            ChangeMapScrollViewSize();
        }
        
        private void ChangeMapScrollViewSize()
        {
            mapViewer.Height = 0.75 * PsmAnnotationGrid.ActualHeight;
            mapViewer.Width = 0.75 * PsmAnnotationGrid.ActualWidth;

            ChangeMapScrollViewVisibility();
        }

        private void ChangeMapScrollViewVisibility()
        {
            if (mapViewer.Width != 0 && mapGrid.Width > mapViewer.Width)
            {
                mapViewer.HorizontalScrollBarVisibility = ScrollBarVisibility.Visible;
            }
            else
            {
                mapViewer.HorizontalScrollBarVisibility = ScrollBarVisibility.Hidden;
            }


            if (mapViewer.Height != 0 && mapGrid.Height > mapViewer.Height)
            {
                mapViewer.VerticalScrollBarVisibility = ScrollBarVisibility.Visible;
            }
            else
            {
                mapViewer.VerticalScrollBarVisibility = ScrollBarVisibility.Hidden;
            }
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

            // colors for peptide lines by protease
            proteaseByColor = new Dictionary<string, Color>();
            proteaseByColor["trypsin"] = Colors.Blue;

            // protein groups according to base sequence
            proteinGroups = new Dictionary<string, string>();
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
                case ".tsv":
                    proteinFilePath = filePath;
                    proteinGroupFileNameLabel.Text = filePath;
                    break;
                case ".psmtsv":
                    psmFilePath = filePath;
                    psmFileNameLabel.Text = filePath;
                    //psmFileNameLabelStat.Text = filePath;
                    break;
                default:
                    MessageBox.Show("Cannot read file type: " + theExtension);
                    break;
            }
        }

        private void LoadProteinGroups(string filename)
        {
            try
            {
                proteinGroups.Clear();
                proteinGroupsForTreeView = new Dictionary<string, ProteinForTreeView>();

                int lineCount = 0;
                foreach (string line in File.ReadLines(filename))
                {
                    ++lineCount;
                    if (lineCount == 1)
                    {
                        continue;
                    }

                    var spl = line.Split('\t');
                    string accession = spl[0];
                    string name = spl[3];
                    string sequence = spl[11];
                    var uniquePeptides = spl[6].Trim().Length > 0 ? spl[6].Split('|') : new string[0];
                    var sharedPeptides = spl[7].Trim().Length > 0 ? spl[7].Split('|') : new string[0];

                    var ptv = CreateNewProtein(accession, sequence, name, new Dictionary<string, PeptideForTreeView>(), new Dictionary<string, PeptideForTreeView>());

                    // store unique peptides
                    foreach (string pep in uniquePeptides)
                    {
                        var pepTV = new PeptideForTreeView(pep, ptv);
                        ptv.Children.Add(pepTV);
                        ptv.UniquePeptides.Add(pep, pepTV);
                        ptv.AllPeptides.Add(pep, pepTV);
                    }

                    // store shared peptides
                    foreach (string pep in sharedPeptides)
                    {
                        var pepTV = new PeptideForTreeView(pep, ptv);
                        ptv.Children.Add(pepTV);
                        ptv.SharedPeptides.Add(pep, pepTV);
                        ptv.AllPeptides.Add(pep, pepTV);
                    }
                }
            }
            catch (Exception e)
            {
                MessageBox.Show("Could not open protein groups file:\n" + e.Message);
            }
        }

        private void LoadAndDisplayPsms(string filename)
        {
            allPsms.Clear();

            string fileNameWithExtension = Path.GetFileName(spectraFilePath);
            string fileNameWithoutExtension = Path.GetFileNameWithoutExtension(spectraFilePath);

            try
            {
                // TODO: print warnings
                foreach (var psm in PsmTsvReader.ReadTsv(filename, out List<string> warnings))
                {
                    if (spectraFilePath == null || psm.Filename == fileNameWithExtension || psm.Filename == fileNameWithoutExtension || psm.Filename.Contains(fileNameWithoutExtension))
                    {
                        allPsms.Add(psm);
                        FilterPsm(psm);
                    }
                }

                GroupAmbiguousPsms();
            }
            catch (Exception e)
            {
                MessageBox.Show("Could not open PSM file:\n" + e.Message);
            }
        }

        // called when filter settings are changed
        private void Refilter()
        {
            // for each peptide, delete all psms - what's the best way to do thissss
            foreach (ProteinForTreeView protein in proteinTree)
            {
                protein.AllPeptides.Values.ToList().ForEach(pep => pep.Children.Clear());
            }

            foreach (PsmFromTsv psm in allPsms)
            {
                FilterPsm(psm);
            }

            GroupAmbiguousPsms();
        }

        // group remaining psms that do not match any protein groups
        private void GroupAmbiguousPsms()
        {            
            foreach (PsmFromTsv psm in filteredListOfPsms.Except(psmsWithMatch.Distinct()))
            {
                if (proteinGroupsForTreeView.ContainsKey(psm.ProteinAccession))
                {
                    AddPsmToTreeView(proteinGroupsForTreeView[psm.ProteinAccession], psm);
                }
                else
                {
                    // create new protein group
                    var protein = CreateNewProtein(psm.ProteinAccession, null, psm.ProteinName);
                    AddPsmToTreeView(protein, psm);
                }
            }
        }

        // create a new instance of a protein group
        private ProteinForTreeView CreateNewProtein(string accession, string sequence, string name, Dictionary<string, PeptideForTreeView> unique = null, Dictionary<string, PeptideForTreeView> shared = null)
        {
            proteinGroups.Add(accession, sequence);
            var ptv = new ProteinForTreeView(name, accession, unique, shared, new Dictionary<string, PeptideForTreeView>());
            proteinGroupsForTreeView.Add(accession, ptv);
            proteinTree.Add(ptv);

            return ptv;
        }

        // add psm to corresponding peptide
        private void AddPsmToTreeView(ProteinForTreeView protein, PsmFromTsv psm)
        {
            PeptideForTreeView peptide = null;

            if (protein.AllPeptides.ContainsKey(psm.BaseSeq)) // retrieve corresponding peptide 
            {
                peptide = protein.AllPeptides[psm.BaseSeq];
            }
            else // psm doesnt match any peptide, create new peptide
            {
                peptide = new PeptideForTreeView(psm.BaseSeq, protein);
                protein.Children.Add(peptide);
                protein.AllPeptides.Add(psm.BaseSeq, peptide);
            }

            int i = 0;
            while (i < peptide.Children.Count) // O(N) 
            {
                // add to sorted collection
                if (psm.Ms2ScanNumber < peptide.Children[i].ScanNo)
                {
                    peptide.Children.Insert(i, new PsmForTreeView(psm.Ms2ScanNumber, psm, peptide));
                    return;
                }
                ++i;
            }

            peptide.Children.Add(new PsmForTreeView(psm.Ms2ScanNumber, psm, peptide));
        }

        // filter psms based on settings and group psms by protein
        private void FilterPsm(PsmFromTsv psm)
        {
            if (psm.QValue <= metaDrawFilterSettings.QValueFilter && (psm.QValueNotch < metaDrawFilterSettings.QValueFilter || psm.QValueNotch == null)
                    && (psm.DecoyContamTarget == "T" || (psm.DecoyContamTarget == "D" && metaDrawFilterSettings.ShowDecoys) || (psm.DecoyContamTarget == "C" && metaDrawFilterSettings.ShowContaminants)))
            {
                filteredListOfPsms.Add(psm);

                // group psms by protein
                foreach (var protein in proteinGroupsForTreeView.Where(x => psm.ProteinAccession.Contains(x.Key) || x.Key.Contains(psm.ProteinAccession)).Select(x => x.Value))
                {
                    AddPsmToTreeView(protein, psm);
                    psmsWithMatch.Add(psm);
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

                    childPsmModel.DrawPeptideSpectralMatch(childScan, psmToDraw, metaDrawGraphicalSettings.ShowMzValues,
                        metaDrawGraphicalSettings.ShowAnnotationCharges, metaDrawGraphicalSettings.AnnotatedFontSize, metaDrawGraphicalSettings.BoldText);

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
            mainViewModel.DrawPeptideSpectralMatch(msDataScanToDraw, psmToDraw, metaDrawGraphicalSettings.ShowMzValues,
                metaDrawGraphicalSettings.ShowAnnotationCharges, metaDrawGraphicalSettings.AnnotatedFontSize, metaDrawGraphicalSettings.BoldText);

            // draw annotated base sequence
            DrawAnnotatedBaseSequence(psmToDraw);
        }

        /// <summary>
        /// Event triggers when a different cell is selected in the PSM data grid
        /// </summary>
        private void proteinTreeView_SelectedCellsChanged(object sender, RoutedPropertyChangedEventArgs<Object> e)
        {
            OnSelectionChanged();
        }

        private void psmAnnotationSizeChanged(object sender, SizeChangedEventArgs e) {
            ChangeMapScrollViewSize();
        }

        private void OnSelectionChanged()
        {
            // can select protein or peptide only at a time (not both)
            if (proteinTreeView.SelectedItem == null)
            {
                return;
            }

            // draw the selected PSM
            propertyView.Clear();

            // differentiate item click protein VS psm
            var selection = proteinTreeView.SelectedItem.GetType().Name;
            PsmForTreeView psmToDisplay = null; 

            switch (selection)
            {
                case "ProteinForTreeView":
                    {
                        ProteinForTreeView selectedItem = (ProteinForTreeView)proteinTreeView.SelectedItem;
                        itemsControlSampleViewModel.Data.Clear();

                        if (proteinGroups[selectedItem.Accession] != null)
                        {
                            DrawSequenceCoverageMap(selectedItem);
                            ChangeMapScrollViewVisibility();
                        }
                    }
                    break;

                case "PeptideForTreeView":
                    {
                        // find psm with best score
                        PeptideForTreeView selectedItem = (PeptideForTreeView)proteinTreeView.SelectedItem;
                        psmToDisplay = selectedItem.Children.OrderByDescending(psm => psm.Psm.Score).First();
                        goto case "display";
                    }

                case "PsmForTreeView":
                    {
                        psmToDisplay = (PsmForTreeView)proteinTreeView.SelectedItem;
                        goto case "display";
                    }

                case "display":
                    {
                        mapViewer.Visibility = Visibility.Collapsed;
                        legend.Visibility = Visibility.Collapsed;
                        displayPsm(psmToDisplay);
                    }
                    break;
            }

            
        }

        // displays psm spectrum annotation
        private void displayPsm(PsmForTreeView selection)
        {
            PsmFromTsv row = selection.Psm;
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

        private void selectProteinGroupFileButton_Click(object sender, RoutedEventArgs e)
        {
            Microsoft.Win32.OpenFileDialog openFileDialog1 = new Microsoft.Win32.OpenFileDialog
            {
                Filter = "Result Files(*.tsv)|*.tsv",
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

            var result = metaDrawFilterSettings.ShowDialog();

            Refilter();
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

            if (psmFilePath == null)
            {
                MessageBox.Show("Please add a search result file.");
                return;
            }

            if (proteinFilePath == null)
            {
                MessageBox.Show("Please add a protein search result file."); //fixme
                return;
            }

            // load the spectra file
            (sender as Button).IsEnabled = false;
            selectSpectraFileButton.IsEnabled = false;
            selectPsmFileButton.IsEnabled = false;
            selectProteinGroupFileButton.IsEnabled = false;
            prgsFeed.IsOpen = true;
            prgsText.Content = "Loading spectra file...";

            var slowProcess = Task<MsDataFile>.Factory.StartNew(() => spectraFileManager.LoadFile(spectraFilePath, new CommonParameters(trimMsMsPeaks: false)));
            await slowProcess;
            MsDataFile = slowProcess.Result;

            // load protein groups
            this.prgsText.Content = "Loading protein groups...";
            LoadProteinGroups(proteinFilePath);

            // load the PSMs
            this.prgsText.Content = "Loading PSMs...";
            LoadAndDisplayPsms(psmFilePath);

            // done loading - restore controls
            prgsFeed.IsOpen = false;
            (sender as Button).IsEnabled = true;
            selectSpectraFileButton.IsEnabled = true;
            selectPsmFileButton.IsEnabled = true;
            selectProteinGroupFileButton.IsEnabled = true;
        }

        //private void loadFilesButtonStat_Click(object sender, RoutedEventArgs e)
        //{
        //    // check for validity
        //    if (tsvResultsFilePath == null)
        //    {
        //        MessageBox.Show("Please add a search result file.");
        //        return;
        //    }

        //    (sender as Button).IsEnabled = false;
        //    selectPsmFileButtonStat.IsEnabled = false;
        //    prgsFeedStat.IsOpen = true;

        //    // load the PSMs
        //    this.prgsTextStat.Content = "Loading PSMs...";
        //    LoadPsmsStat(tsvResultsFilePath);

        //    // done loading - restore controls
        //    this.prgsFeedStat.IsOpen = false;
        //    (sender as Button).IsEnabled = true;
        //    selectPsmFileButtonStat.IsEnabled = true;
        //}

        private void LoadPsmsStat(string filepath)
        {
            LoadAndDisplayPsms(filepath);
        }

        private void searchBox_TextChanged(object sender, TextChangedEventArgs e)
        {
            SearchTimer.Set();
        }

        // handler for searching through tree
        private void searchBox_TextChangedHandler(object sender, EventArgs e)
        {
            string userInput = searchBox.Text;

            if (string.IsNullOrEmpty(userInput))
            {
                proteinTreeView.DataContext = proteinTree;
                return;
            }

            searchProtein(userInput);
            proteinTreeView.DataContext = filteredTree;
            SearchTimer.Timer.Stop();
        }

        // search through protein list based on user input
        private void searchProtein(string txt)
        {
            filteredTree.Clear();
            foreach (var protein in proteinTree)
            {
                if (protein.DisplayName.ToUpper().Contains(txt.ToUpper()))
                {
                    filteredTree.Add(protein);
                }
                else
                {
                    searchPeptide(protein, txt);
                }
            }
        }
        
        // search through peptide list based on user input
        private void searchPeptide(ProteinForTreeView protein, string txt)
        {
            // create copy of protein
            var prot = new ProteinForTreeView(protein.DisplayName, protein.Accession, protein.UniquePeptides, protein.SharedPeptides, protein.AllPeptides);

            foreach (var peptide in protein.Children)
            {
                if (peptide.DisplayName.ToUpper().Contains(txt.ToUpper()))
                {
                    prot.Children.Add(peptide);
                }
                else
                {
                    prot = searchPsm(prot, peptide, txt);
                }
            }

            if (prot.Children.Count() > 0)
            {
                prot.Expanded = true;
                filteredTree.Add(prot);
            }
        }

        // search through psm list based on user input
        private ProteinForTreeView searchPsm(ProteinForTreeView prot, PeptideForTreeView peptide, string txt)
        {
            // create copy of peptide
            var pep = new PeptideForTreeView(peptide.DisplayName, peptide.Parent);

            foreach (var psm in peptide.Children)
            {
                if (psm.ScanNo.ToString().StartsWith(txt))
                {
                    pep.Expanded = true;
                    pep.Children.Add(psm);
                }
            }

            if (pep.Children.Count() > 0)
            {
                prot.Children.Add(pep);
            }

            return prot;
        }

        // split sequence into fixed amino acids per line
        private List<string> Split(string sequence, double spacing)
        {
            int size = Convert.ToInt32(mapGrid.Width / spacing);
            var splitSequence = Enumerable.Range(0, sequence.Length / size).Select(i => sequence.Substring(i * size, size)).ToList();
            splitSequence.Add(sequence.Substring(splitSequence.Count() * size));

            return splitSequence;
        }

        private void DrawSequenceCoverageMap(ProteinForTreeView protein) //string accession, Dictionary<string, PeptideForTreeView> uniquePeptides, Dictionary<string, PeptideForTreeView> sharedPeptides)
        {
            string protease = "trypsin"; // only works for single protease for now
            string seqCoverage = proteinGroups[protein.Accession];
            mapViewer.Visibility = Visibility.Visible;

            BaseDraw.clearCanvas(map);
            mainViewModel.Model = null; // clear plot
            BaseDraw.clearCanvas(canvas);

            double spacing = 22;
            int height = 10;
            int totalHeight = 0;
            int accumIndex = 0;

            foreach (string seq in seqCoverage.Split('|'))
            {
                var splitSeq = Split(seq, spacing);

                var allPeptides = new List<string>(protein.AllPeptides.Keys);
                foreach (var line in splitSeq)
                {
                    List<int> indices = new List<int>();

                    // draw sequence
                    for (int r = 0; r < line.Length; r++)
                    {
                        SequenceCoverageMap.txtDrawing(map, new Point(r * spacing + 10, height), line[r].ToString().ToUpper(), Brushes.Black);
                    }

                    // highlight partial peptide sequences (broken off into multiple lines)
                    if (partialPeptideMatches.Count > 0)
                    {
                        var temp = new Dictionary<string, int>(partialPeptideMatches);
                        partialPeptideMatches.Clear();

                        foreach (var peptide in temp)
                        {
                            if (MatchPeptideSequence(peptide.Key, line, 0, peptide.Value, accumIndex - peptide.Value == seq.IndexOf(peptide.Key)))
                            {
                                int start = 0;
                                int end = Math.Min(start + peptide.Key.Length - peptide.Value - 1, line.Length - 1);
                                SequenceCoverageMap.Highlight(start, end, map, indices, height, proteaseByColor[protease], protein.UniquePeptides.Keys.Any(u => u.Contains(peptide.Key))); // draw line for peptide sequence
                            }
                        }
                    }

                    // highlight full peptide sequences
                    for (int i = 0; i < line.Length; ++i)
                    {
                        var temp = new List<string>(allPeptides);
                        foreach (string peptide in temp)
                        {
                            if (MatchPeptideSequence(peptide, line, i, 0, accumIndex + i == seq.IndexOf(peptide)))
                            {
                                int start = i;
                                int end = Math.Min(start + peptide.Length - 1, line.Length - 1);
                                SequenceCoverageMap.Highlight(start, end, map, indices, height, proteaseByColor[protease], protein.UniquePeptides.Keys.Any(u => u.Contains(peptide)));
                                allPeptides.Remove(peptide);
                            }
                        }
                    }

                    height += 50;
                    accumIndex += line.Length;
                }

                totalHeight += splitSeq.Count() * 50;
                height += 50;
            }

            mapGrid.Height = totalHeight + 50 * seqCoverage.Split('|').Count();
            SequenceCoverageMap.drawLegend(legend, proteaseByColor, proteases, legendGrid); 
        }

        Dictionary<string, int> partialPeptideMatches = new Dictionary<string, int>();

        private bool MatchPeptideSequence(string peptide, string line, int proteinStartIndex, int peptideStartIndex, bool fits)
        {
            bool match = true;
            char current;
            int m;

            // compare protein sequence and peptide
            for (m = 0; peptideStartIndex < peptide.Length && match; ++m, ++peptideStartIndex)
            {
                if (proteinStartIndex + m >= line.Length)
                {
                    if (!fits)
                    {
                        match = false;
                    }
                    else
                    {
                        partialPeptideMatches.Add(peptide, peptideStartIndex);
                    }
                    break;
                }

                current = line[proteinStartIndex + m];
                if (current != peptide[peptideStartIndex])
                {
                    match = false;
                }
            }

            return match;
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

        private void PDFButton_Click(object sender, RoutedEventArgs e)
        {
            object temp = null;
            var selections = new List<object>();

            // select best scoring psm for any peptide selections
            foreach (object selectedItem in proteinTreeView.SelectedItems)
            {
                if (selectedItem.GetType().Name.Equals("PeptideForTreeView"))
                {
                    PeptideForTreeView peptide = (PeptideForTreeView)selectedItem;
                    PsmForTreeView bestScoringPsm = peptide.Children.OrderByDescending(psm => psm.Psm.Score).First();
                    selections.Add(bestScoringPsm);
                    continue;
                }

                selections.Add(selectedItem);
            }

            if (proteinTreeView.SelectedItems.Count == 0)
            {
                MessageBox.Show("Please select at least one scan or coverage map to export");
            }
            else
            {
                int numberOfScansToExport = selections.Distinct().Count();

                foreach (object selectedItem in selections.Distinct())
                {
                    if (temp == null)
                    {
                        temp = selectedItem;
                    }

                    var selection = selectedItem.GetType().Name;
                    switch (selection)
                    {
                        case "ProteinForTreeView":
                            {
                                ProteinForTreeView protein = (ProteinForTreeView)selectedItem;
                                string myString = illegalInFileName.Replace(protein.DisplayName, "");
                                myString = myString.Length > 30 ? myString.Substring(0, 30) : myString;

                                string filePath = Path.Combine(Path.GetDirectoryName(psmFilePath), "MetaDrawExport", myString.Trim(), myString + ".pdf");
                                string dir = Path.GetDirectoryName(filePath);

                                if (!Directory.Exists(dir))
                                {
                                    Directory.CreateDirectory(dir);
                                }
                                
                                DrawPdfCoverageMap(protein, mapGrid, filePath);
                            }
                            break;

                        case "PsmForTreeView":
                            {
                                PsmForTreeView psm = (PsmForTreeView)selectedItem;
                                var proteinFolder = illegalInFileName.Replace(psm.Parent.Parent.DisplayName, "");
                                proteinFolder = proteinFolder.Length > 30 ? proteinFolder.Substring(0, 30) : proteinFolder;

                                string peptideFolder = illegalInFileName.Replace(psm.Parent.DisplayName, "");
                                peptideFolder = peptideFolder.Length > 30 ? peptideFolder.Substring(0, 30) : peptideFolder;

                                string filePath = Path.Combine(Path.GetDirectoryName(psmFilePath), "MetaDrawExport", proteinFolder.Trim(), peptideFolder.Trim(), psm.ScanNo + ".pdf");
                                string dir = Path.GetDirectoryName(filePath);

                                MsDataScan msDataScanToDraw = MsDataFile.GetOneBasedScan(psm.ScanNo);

                                if (!Directory.Exists(dir))
                                {
                                    Directory.CreateDirectory(dir);
                                }

                                DrawPdfAnnotatedBaseSequence(psm.Psm, canvas, filePath); // captures the annotation for the pdf
                                mainViewModel.DrawPeptideSpectralMatchPdf(msDataScanToDraw, psm.Psm, filePath, numberOfScansToExport > 1);
                            }
                            break;
                    }
                }

                RestoreDefault(temp);
                MessageBox.Show(string.Format("{0} PDFs exported", numberOfScansToExport));
            }
        }

        // redraw display based on first selection
        private void RestoreDefault(object defaultSelection)
        {
            mapGrid.Width = 485;
            var selection = defaultSelection.GetType().Name;
            switch (selection)
            {
                case "ProteinForTreeView":
                    {
                        ProteinForTreeView protein = (ProteinForTreeView)defaultSelection;
                        DrawSequenceCoverageMap(protein);
                    }
                    break;

                case "PsmForTreeView":
                    {
                        PsmForTreeView psm = (PsmForTreeView)defaultSelection;
                        DrawPsm(psm.Psm.Ms2ScanNumber, psm.Psm.FullSequence);
                    }
                    break;
            }
        }
        
        private void DrawPdfCoverageMap(ProteinForTreeView protein, Grid mapGrid, string path)
        {
            // draw coverage map
            mapGrid.Width = 1000;
            DrawSequenceCoverageMap(protein);

            string pathToCoverageMap = Path.Combine(Path.GetDirectoryName(path), "map.png");
            CustomPdfWriter.RenderImage((int)mapGrid.Width, (int)mapGrid.Height, pathToCoverageMap, map);

            // draw legend
            SequenceCoverageMap.drawLegend(legend, proteaseByColor, proteases, legendGrid);
            string pathToLegend = Path.Combine(Path.GetDirectoryName(path), "legend.png");
            CustomPdfWriter.RenderImage((int)(legend.Width), 50, pathToLegend, legend);

            CustomPdfWriter.WriteToPdfMetaDraw(mapGrid.Width, mapGrid.Height, pathToCoverageMap, pathToLegend, path);
        }

        private void DrawPdfAnnotatedBaseSequence(PsmFromTsv psm, Canvas canvas, string path)
        {
            if (psm.CrossType == null)
            {
                DrawAnnotatedBaseSequence(psm);
            }

            string pathToBaseSeq = Path.Combine(Path.GetDirectoryName(path), "annotation.png");
            CustomPdfWriter.RenderImage((int)canvas.Width, 600, pathToBaseSeq, canvas);
        }

        //private async void PlotSelected(object sender, SelectionChangedEventArgs e)
        //{
        //    var listview = sender as ListView;
        //    var plotName = listview.SelectedItem as string;

        //    if (filteredListOfPsms.Count == 0)
        //    {
        //        MessageBox.Show("There are no PSMs to analyze.\n\nLoad the current file or choose a new file.");
        //        return;
        //    }
        //    PlotModelStat plot = await Task.Run(() => new PlotModelStat(plotName, filteredListOfPsms));
        //    plotViewStat.DataContext = plot;
        //}

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
