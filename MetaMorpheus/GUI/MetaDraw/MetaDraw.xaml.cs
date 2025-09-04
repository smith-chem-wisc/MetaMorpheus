using Easy.Common.Extensions;
using EngineLayer;
using GuiFunctions;
using GuiFunctions.MetaDraw;
using Omics.Fragmentation;
using OxyPlot;
using Readers;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Data;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Input;
using System.Windows.Media;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for MetaDraw.xaml
    /// </summary>
    public partial class MetaDraw : Window
    {
        private ParentChildScanPlotsView itemsControlSampleViewModel;
        private MetaDrawLogic MetaDrawLogic;
        private readonly DataTable propertyView;
        private ObservableCollection<string> plotTypes;
        private ObservableCollection<string> PsmStatPlotFiles;
        public PtmLegendViewModel PtmLegend;
        private ObservableCollection<ModTypeForTreeViewModel> Modifications = new ObservableCollection<ModTypeForTreeViewModel>();
        private static List<string> AcceptedSpectraFormats => SpectrumMatchFromTsvHeader.AcceptedSpectraFormats.Concat(new List<string> { ".msalign", ".tdf", ".tdf_bin" }).Select(format => format.ToLower()).ToList();
        private static List<string> AcceptedResultsFormats = new List<string> { ".psmtsv", ".tsv" };
        private static List<string> AcceptedSpectralLibraryFormats = new List<string> { ".msp" };
        private FragmentationReanalysisViewModel FragmentationReanalysisViewModel;
        public ChimeraAnalysisTabViewModel ChimeraAnalysisTabViewModel { get; set; }
        public DeconExplorationTabViewModel DeconExplorationViewModel { get; set; } = new();
        public BioPolymerTabViewModel BioPolymerTabViewModel { get; set; } 

        public MetaDraw(string[]? filesToLoad = null)
        {
            InitializeComponent();

            MetaDrawLogic = new MetaDrawLogic();
            SettingsButtonControl.SettingsChanged += RefreshPlotsAfterSettingsChange;
            BindingOperations.EnableCollectionSynchronization(MetaDrawLogic.SpectralMatchResultFilePaths, MetaDrawLogic.ThreadLocker);
            BindingOperations.EnableCollectionSynchronization(MetaDrawLogic.SpectraFilePaths, MetaDrawLogic.ThreadLocker);
            BindingOperations.EnableCollectionSynchronization(MetaDrawLogic.FilteredListOfPsms, MetaDrawLogic.ThreadLocker);
            BindingOperations.EnableCollectionSynchronization(MetaDrawLogic.SpectralMatchesGroupedByFile, MetaDrawLogic.ThreadLocker);

            itemsControlSampleViewModel = new ParentChildScanPlotsView();
            DeconExplorationTabView.DataContext = DeconExplorationViewModel;
            ParentChildScanViewPlots.DataContext = itemsControlSampleViewModel;
            AdditionalFragmentIonControl.DataContext = FragmentationReanalysisViewModel ??= new FragmentationReanalysisViewModel();
            AdditionalFragmentIonControl.LinkMetaDraw(this);
            BioPolymerTabViewModel = new BioPolymerTabViewModel(MetaDrawLogic);
            BioPolymerCoverageTabView.DataContext = BioPolymerTabViewModel;
            

            propertyView = new DataTable();
            propertyView.Columns.Add("Name", typeof(string));
            propertyView.Columns.Add("Value", typeof(string));
            dataGridProperties.DataContext = propertyView.DefaultView;

            dataGridScanNums.DataContext = MetaDrawLogic.PeptideSpectralMatchesView;


            Title = "MetaDraw: version " + GlobalVariables.MetaMorpheusVersion;
            base.Closing += this.OnClosing;

            ParentChildScanView.Visibility = Visibility.Collapsed;

            PsmStatPlotFiles = new ObservableCollection<string>();
            selectSourceFileListBox.DataContext = PsmStatPlotFiles;
            plotTypes = new ObservableCollection<string>();
            SetUpPlots();
            plotsListBox.ItemsSource = plotTypes;

            // if files piped in from a MetaMorpheus search, load them immediately
            if (filesToLoad is not null)
            {
                filesToLoad.ForEach(AddFile);
            }
        }

        private void Window_Drop(object sender, DragEventArgs e)
        {
            string[] files = ((string[])e.Data.GetData(DataFormats.FileDrop))?.OrderBy(p => p).ToArray();

            if (files != null)
            {
                foreach (var draggedFilePath in files)
                {
                    if (File.Exists(draggedFilePath) | (Directory.Exists(draggedFilePath) && Regex.IsMatch(draggedFilePath, @".d$")) )
                    {
                        AddFile(draggedFilePath);
                    }
                }
            }
        }

        private void AddFile(string filePath)
        {
            var theExtension = GlobalVariables.GetFileExtension(filePath).ToLowerInvariant();

            if (AcceptedSpectraFormats.Contains(theExtension))
            {
                // If a bruker timsTof file was selected, we actually want the parent folder
                if(theExtension == ".tdf" || theExtension == ".tdf_bin")
                {
                    filePath = Path.GetDirectoryName(filePath);
                }
                if (!MetaDrawLogic.SpectraFilePaths.Contains(filePath))
                {
                    MetaDrawLogic.SpectraFilePaths.Add(filePath);

                    if (MetaDrawLogic.SpectraFilePaths.Count == 1)
                    {
                        spectraFileNameLabel.Text = filePath;
                    }
                    else
                    {
                        spectraFileNameLabel.Text = "[Mouse over to view files]";
                    }

                    spectraFileNameLabel.ToolTip = string.Join("\n", MetaDrawLogic.SpectraFilePaths);
                    resetSpectraFileButton.IsEnabled = true;
                }
            }
            else if (AcceptedResultsFormats.Contains(theExtension))
            {
                if (!MetaDrawLogic.SpectralMatchResultFilePaths.Contains(filePath))
                {
                    MetaDrawLogic.SpectralMatchResultFilePaths.Add(filePath);

                    if (MetaDrawLogic.SpectralMatchResultFilePaths.Count == 1)
                    {
                        psmFileNameLabel.Text = filePath;
                        psmFileNameLabelStat.Text = filePath;
                    }
                    else
                    {
                        psmFileNameLabel.Text = "[Mouse over to view files]";
                        psmFileNameLabelStat.Text = "[Mouse over to view files]";
                    }

                    psmFileNameLabel.ToolTip = string.Join("\n", MetaDrawLogic.SpectralMatchResultFilePaths);
                    resetPsmFileButton.IsEnabled = true;

                    psmFileNameLabelStat.ToolTip = string.Join("\n", MetaDrawLogic.SpectralMatchResultFilePaths);
                    resetPsmFileButtonStat.IsEnabled = true;
                }
            }
            else if (AcceptedSpectralLibraryFormats.Contains(theExtension))
            {
                // TODO: display this somewhere in the GUI
                if (!MetaDrawLogic.SpectralLibraryPaths.Contains(filePath))
                {
                    MetaDrawLogic.SpectralLibraryPaths.Add(filePath);
                    specLibraryLabel.Text = filePath;
                    specLibraryLabel.ToolTip = string.Join("\n", MetaDrawLogic.SpectralLibraryPaths);
                    resetSpecLibraryButton.IsEnabled = true;
                }
            }
            else if (GlobalVariables.AcceptedDatabaseFormats.Contains(theExtension))
            {
                BioPolymerTabViewModel.DatabasePath = filePath;
            }
            else
            {
                MessageBox.Show("Cannot read file type: " + theExtension);
            }
        }

        /// <summary>
        /// Event triggers when a different cell is selected in the PSM data grid
        /// <remarks>
        ///  if sender is FragmentationReanalysisViewModel, then this method was run by clicking the search button on the FragmentationReanalysisViewModel
        /// </remarks>
        /// </summary>
        private void dataGridScanNums_SelectedCellsChanged(object sender, SelectedCellsChangedEventArgs e)
        {
            if (dataGridScanNums.SelectedItem == null || sender == null)
            {
                ClearPresentationArea();
                return;
            }

            var strings = sender.ToString();

            MetaDrawLogic.CleanUpCurrentlyDisplayedPlots();
            wholeSequenceCoverageHorizontalScroll.ScrollToLeftEnd();
            plotView.Visibility = Visibility.Visible;
            SpectrumMatchFromTsv psm = (SpectrumMatchFromTsv)dataGridScanNums.SelectedItem;

            List<MatchedFragmentIon> oldMatchedIons = null;
            if (FragmentationReanalysisViewModel.Persist && sender is DataGrid)
            {
                oldMatchedIons = psm.MatchedIons;
                ReplaceFragmentIonsOnPsmFromFragmentReanalysisViewModel(psm);
            }

            wholeSequenceCoverageHorizontalScroll.Visibility = Visibility.Visible;
            

            SetSequenceDrawingPositionSettings(true);
            // Psm selected from ambiguous dropdown => adjust the psm to be drawn
            // Clicking the research button on an ambiguous psm => research with new ions
            if (psm.FullSequence.Contains('|') && (sender.ToString() == "System.Object" || sender is FragmentationReanalysisViewModel))
            {
                // From chimeric scan view to child scan view with ambiguous selected
                if (AmbiguousSequenceOptionBox.SelectedItem == null) 
                {
                    // set to first and break the loop if casted successfully
                    foreach (var ambiguousResult in AmbiguousSequenceOptionBox.Items)
                    {
                        psm = ambiguousResult as SpectrumMatchFromTsv;
                        if (psm != null)
                            break;
                    }

                    AmbiguousWarningTextBlocks.Visibility = Visibility.Collapsed;
                    AmbiguousSequenceOptionBox.Visibility = Visibility.Visible;
                    AmbiguousSequenceOptionBox.SelectedItem = psm;
                }
                // selecting a different ambiguous result from the combobox in child scan view
                else
                {
                    psm = (SpectrumMatchFromTsv)AmbiguousSequenceOptionBox.SelectedItem;
                }

                if (FragmentationReanalysisViewModel.Persist || sender is FragmentationReanalysisViewModel)
                {
                    oldMatchedIons = psm.MatchedIons;
                    ReplaceFragmentIonsOnPsmFromFragmentReanalysisViewModel(psm);
                }
            }
            // Selection of ambiguous psm => clean up the canvases and show the option box
            else if(psm.FullSequence.Contains('|') && sender.ToString() != "System.Object")
            {
                // clear all drawings of the previous non-ambiguous psm
                ClearPresentationArea();

                AmbiguousWarningTextBlocks.Visibility = Visibility.Visible;
                AmbiguousSequenceOptionBox.Visibility = Visibility.Visible;

                // create a psm object for each ambiguous option and add it to the dropdown box
                var fullSeqs = psm.FullSequence.Split('|');
                foreach (var fullSeq in fullSeqs)
                {
                    SpectrumMatchFromTsv oneAmbiguousPsm = psm.ReplaceFullSequence(fullSeq);
                    AmbiguousSequenceOptionBox.Items.Add(oneAmbiguousPsm);
                }
                return;
            }
            
            // Selection of non-ambiguous psm => clear psms in the drop down
            else if (!psm.FullSequence.Contains('|'))
            {
                AmbiguousSequenceOptionBox.Items.Clear();
                AmbiguousWarningTextBlocks.Visibility = Visibility.Collapsed;
                AmbiguousSequenceOptionBox.Visibility = Visibility.Collapsed;
                wholeSequenceCoverageHorizontalScroll.Visibility = Visibility.Visible;
                if (psm.BaseSeq.Length > MetaDrawSettings.NumberOfAAOnScreen)
                {
                    GrayBox.Opacity = 0.7;
                }
                else
                {
                    GrayBox.Opacity = 0;
                }
            }

            // display the ion and elements correctly
            MetaDrawSettings.DrawMatchedIons = true;


            // define initial limits for sequence annotation
            double maxDisplayedPerRow = (int)Math.Round((UpperSequenceAnnotaiton.ActualWidth - 10) / MetaDrawSettings.AnnotatedSequenceTextSpacing, 0) + 7;
            MetaDrawSettings.SequenceAnnotationSegmentPerRow = (int)Math.Floor(maxDisplayedPerRow / (double)(MetaDrawSettings.SequenceAnnotaitonResiduesPerSegment + 1));

            // draw the annotated spectrum
            MetaDrawLogic.DisplaySequences(stationarySequenceCanvas, scrollableSequenceCanvas, sequenceAnnotationCanvas, psm);
            MetaDrawLogic.DisplaySpectrumMatch(plotView, psm, itemsControlSampleViewModel, out var errors);

            // add ptm legend if desired
            if (MetaDrawSettings.ShowLegend)
            {
                int descriptionLineCount = MetaDrawSettings.SpectrumDescription.Count(p => p.Value);
                if (psm.Name.IsNotNullOrEmptyOrWhiteSpace())
                {
                    descriptionLineCount += (int)Math.Floor((psm.Name.Length - 20) / (double)SpectrumMatchPlot.MaxCharactersPerDescriptionLine);
                }
                if (psm.Accession.Length > 10)
                    descriptionLineCount++;
                double verticalOffset = descriptionLineCount * 1.4 * MetaDrawSettings.SpectrumDescriptionFontSize;
                
                PtmLegend = new PtmLegendViewModel(psm, verticalOffset);
                ChildScanPtmLegendControl.DataContext = PtmLegend;
                SequenceCoveragePtmLegendControl.DataContext = PtmLegend;
            }

            //draw the sequence coverage if not crosslinked
            if (psm.ChildScanMatchedIons == null)
            {
                MetaDrawLogic.DrawSequenceCoverageMap(psm, sequenceText, map); //TODO: figure out how to show coverage on crosslinked peptides
                ParentChildScanView.Visibility = Visibility.Collapsed;
                SequenceCoverageAnnotationView.Visibility = Visibility.Visible;
            }
            else
            {
                ParentChildScanView.Visibility = Visibility.Visible;
                SequenceCoverageAnnotationView.Visibility = Visibility.Collapsed;
            }

            mapViewer.Width = map.Width;

            if (errors != null && errors.Any())
            {
                MessageBox.Show(errors.First());
                return;
            }

            // display PSM properties
            propertyView.Clear();
            System.Reflection.PropertyInfo[] temp = psm.GetType().GetProperties();

            for (int i = 0; i < temp.Length; i++)
            {
                if (temp[i].Name == nameof(psm.MatchedIons))
                {
                    propertyView.Rows.Add(temp[i].Name, string.Join(", ", psm.MatchedIons.Select(p => p.Annotation)));
                }
                else if (temp[i].Name == nameof(psm.VariantCrossingIons))
                {
                    propertyView.Rows.Add(temp[i].Name, string.Join(", ", psm.VariantCrossingIons.Select(p => p.Annotation)));
                }
                else
                {
                    // Hacky fix for some of the properties in IQuantifiableRecord that are only populated as needed by FlashLFQ
                    try
                    {
                        propertyView.Rows.Add(temp[i].Name, temp[i].GetValue(psm, null));
                    }
                    catch
                    {
                        // do nothing
                    }
                }
            }

            // put the original ions back in place if they were altered
            if (oldMatchedIons != null && !psm.MatchedIons.SequenceEqual(oldMatchedIons))
                psm.MatchedIons = oldMatchedIons;
        }

        private void selectSpectraFileButton_Click(object sender, RoutedEventArgs e)
        {
            string filterString = string.Join(";", AcceptedSpectraFormats.Select(p => "*" + p));

            Microsoft.Win32.OpenFileDialog openFileDialog1 = new Microsoft.Win32.OpenFileDialog
            {
                Filter = "Spectra Files(" + filterString + ")|" + filterString,
                FilterIndex = 1,
                RestoreDirectory = true,
                Multiselect = true
            };
            if (openFileDialog1.ShowDialog() == true)
            {
                foreach (var filePath in openFileDialog1.FileNames.OrderBy(p => p))
                {
                    AddFile(filePath);
                }
            }
        }

        private void selectPsmFileButton_Click(object sender, RoutedEventArgs e)
        {
            string filterString = string.Join(";", AcceptedResultsFormats.Concat(AcceptedSpectralLibraryFormats).Select(p => "*" + p));

            Microsoft.Win32.OpenFileDialog openFileDialog1 = new Microsoft.Win32.OpenFileDialog
            {
                Filter = "Result Files(" + filterString + ")|" + filterString,
                FilterIndex = 1,
                RestoreDirectory = true,
                Multiselect = false
            };
            if (openFileDialog1.ShowDialog() == true)
            {
                foreach (var filePath in openFileDialog1.FileNames.OrderBy(p => p))
                {
                    AddFile(filePath);
                }
            }
        }

        private void selectSpecLibraryButton_Click(object sender, RoutedEventArgs e)
        {
            string filterString = string.Join(";", AcceptedResultsFormats.Concat(AcceptedSpectralLibraryFormats).Select(p => "*" + p));

            Microsoft.Win32.OpenFileDialog openFileDialog1 = new Microsoft.Win32.OpenFileDialog
            {
                Filter = "Result Files(" + filterString + ")|" + filterString,
                FilterIndex = 1,
                RestoreDirectory = true,
                Multiselect = false
            };
            if (openFileDialog1.ShowDialog() == true)
            {
                foreach (var filePath in openFileDialog1.FileNames.OrderBy(p => p))
                {
                    AddFile(filePath);
                }
            }
        }

        private void resetFilesButton_Click(object sender, RoutedEventArgs e)
        {
            if (((Button)sender).Name.Equals("resetSpectraFileButton"))
            {
                spectraFileNameLabel.Text = "None Selected";
                MetaDrawLogic.CleanUpSpectraFiles();
            }

            else if (((Button)sender).Name.Equals("resetPsmFileButton"))
            {
                psmFileNameLabel.Text = "None Selected";
                MetaDrawLogic.CleanUpPSMFiles();
            }

            else if (((Button)sender).Name.Equals("resetSpecLibraryButton"))
            {
                specLibraryLabel.Text = "None Selected";
                MetaDrawLogic.CleanUpSpectralLibraryFiles();
            }
            else
            {
                MetaDrawLogic.CleanUpResources();
            }

            // if a psm is selected
            if (MetaDrawLogic.ScrollableSequence != null)
            {
                ClearPresentationArea();
                MetaDrawLogic.FilteredListOfPsms.Clear();
            }
        }

        private void OnClosing(object sender, CancelEventArgs e)
        {
            MetaDrawLogic.CleanUpResources();
        }

        /// <summary>
        /// Will be fired by settings button control if settings change. 
        /// </summary>
        private void RefreshPlotsAfterSettingsChange(object sender, MetaDrawSettingsChangedEventArgs e)
        {
            if (MetaDrawLogic.FilteredListOfPsms.Count == 0)
                return;

            // save current selected PSM
            var selectedItem = dataGridScanNums.SelectedItem as SpectrumMatchFromTsv;
            var selectedChimeraGroup = ChimeraAnalysisTabViewModel?.SelectedChimeraGroup;

            // filter based on new settings
            if (e.FilterChanged)
            {
                MetaDrawLogic.FilterPsms();

                ChimeraAnalysisTabViewModel.ChimeraGroupViewModels.Clear();
                foreach (var chimeraGroup in ChimeraAnalysisTabViewModel.ConstructChimericPsms(MetaDrawLogic.FilteredListOfPsms.ToList(), MetaDrawLogic.MsDataFiles)
                             .OrderByDescending(p => p.Count)
                             .ThenByDescending(p => p.UniqueFragments)
                             .ThenByDescending(p => p.TotalFragments))
                {
                    ChimeraAnalysisTabViewModel.ChimeraGroupViewModels.Add(chimeraGroup);
                }

                foreach (var group in BioPolymerTabViewModel.AllGroups)
                {
                    group.UpdatePropertiesAfterFilter();
                }
            }

            if (e.DataVisualizationChanged && (string)((TabItem)MainTabControl.SelectedItem).Header == "Data Visualization")
            {
                PlotSelected(plotsListBox, null);
            }

            // Reselect items and refresh plots
            if (selectedItem != null)
            {
                dataGridScanNums.SelectedItem = selectedItem;
            }
            if (selectedChimeraGroup != null)
            {
                ChimeraAnalysisTabViewModel.SelectedChimeraGroup = selectedChimeraGroup;
                ChimeraAnalysisTabViewModel.Ms1ChimeraPlot = new Ms1ChimeraPlot(ChimeraAnalysisTabView.ms1ChimeraOverlaPlot, selectedChimeraGroup);
                ChimeraAnalysisTabViewModel.ChimeraSpectrumMatchPlot = new ChimeraSpectrumMatchPlot(ChimeraAnalysisTabView.ms2ChimeraPlot, selectedChimeraGroup);
                ChimeraAnalysisTabViewModel.ChimeraDrawnSequence = new ChimeraDrawnSequence(ChimeraAnalysisTabView.chimeraSequenceCanvas, selectedChimeraGroup, ChimeraAnalysisTabViewModel);
            }
        }

        private async void loadFilesButton_Click(object sender, RoutedEventArgs e)
        {
            // check for validity
            propertyView.Clear();
            if (!MetaDrawLogic.SpectraFilePaths.Any())
            {
                MessageBox.Show("Please add a spectra file.");
                return;
            }

            if (!MetaDrawLogic.SpectralMatchResultFilePaths.Any())
            {
                MessageBox.Show("Please add a search result file.");
                return;
            }

            // load the spectra file
            ToggleButtonsEnabled(false);
 
            prgsFeed.IsOpen = true;
            prgsText.Content = "Loading data...";

            // Add EventHandlers for popup click-in/click-out behaviour
            Deactivated += new EventHandler(prgsFeed_Deactivator);
            Activated += new EventHandler(prgsFeed_Reactivator);

            var slowProcess = Task<List<string>>.Factory.StartNew(() => MetaDrawLogic.LoadFiles(loadSpectra: true, loadPsms: true));
            await slowProcess;

            string directoryPath = Path.Combine(Path.GetDirectoryName(MetaDrawLogic.SpectralMatchResultFilePaths.First()), "MetaDrawExport",
                DateTime.Now.ToString("yyyy-MM-dd", CultureInfo.InvariantCulture));
            //ChimeraAnalysisTabViewModel = new ChimeraAnalysisTabViewModel(MetaDrawLogic.FilteredListOfPsms.ToList(), MetaDrawLogic.MsDataFiles, directoryPath);
            //ChimeraTab.DataContext = ChimeraAnalysisTabViewModel;
            DeconExplorationViewModel.MsDataFiles.Clear();
            foreach (var dataFile in MetaDrawLogic.MsDataFiles)
            {
                DeconExplorationViewModel.MsDataFiles.Add(dataFile.Value);
            }


            BioPolymerTabViewModel.ExportDirectory = directoryPath;
            if (BioPolymerTabViewModel.IsDatabaseLoaded)
                BioPolymerTabViewModel.ProcessSpectralMatches(MetaDrawLogic.AllSpectralMatches);

            var errors = slowProcess.Result;

            if (errors.Any())
            {
                string errorList = string.Join("\n", errors);
                MessageBox.Show(errorList);
            }

            PsmStatPlotFiles.Clear();
            foreach (var item in MetaDrawLogic.SpectralMatchesGroupedByFile)
            {
                PsmStatPlotFiles.Add(item.Key);
            }

            // done loading - restore controls
            this.prgsFeed.IsOpen = false;

            // Remove added EventHandlers
            Deactivated -= new EventHandler(prgsFeed_Deactivator);
            Activated -= new EventHandler(prgsFeed_Reactivator);

            ToggleButtonsEnabled(true);
        }

        /// <summary>
        /// Deactivates the "loading data" popup if one clicks out of the main window
        /// </summary>
        private void prgsFeed_Deactivator(object sender, EventArgs e)
        {
            prgsFeed.IsOpen = false;
        }

        /// <summary>
        /// Reactivates the "loading data" popup if one clicks into the main window
        /// </summary>
        private void prgsFeed_Reactivator(object sender, EventArgs e)
        {
            prgsFeed.IsOpen = true;
        }

        private void TextBox_TextChanged(object sender, TextChangedEventArgs e)
        {
            string txt = (sender as TextBox).Text;
            MetaDrawLogic.FilterPsmsByString(txt);
        }

        /// <summary>
        /// Exports images of the parent and child scan
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        private void PDFButton_Click(object sender, RoutedEventArgs e)
        {
            if (dataGridScanNums.SelectedCells.Count == 0)
            {
                MessageBox.Show("Please select at least one scan to export");
                return;
            }

            SetSequenceDrawingPositionSettings();
            List<SpectrumMatchFromTsv> items = new();

            foreach (var cell in dataGridScanNums.SelectedItems)
            {
                var psm = (SpectrumMatchFromTsv)cell;
                items.Add(psm);
            }

            string directoryPath = Path.Combine(Path.GetDirectoryName(MetaDrawLogic.SpectralMatchResultFilePaths.First()), "MetaDrawExport",
                    DateTime.Now.ToString("yyyy-MM-dd", CultureInfo.InvariantCulture));


            Canvas legendCanvas = null;
            Vector ptmLegendLocationVector = new();
            List<string> errors = new();
            if (((Grid)MetaDrawTabControl.SelectedContent).Name == "PsmAnnotationGrid")
            {

                if (MetaDrawSettings.ShowLegend)
                {
                    PtmLegendControl ptmLegendCopy = new();
                    ptmLegendCopy.DataContext = ChildScanPtmLegendControl.DataContext;
                    legendCanvas = new();
                    legendCanvas.Children.Add(ptmLegendCopy);
                    Size ptmLegendSize = new Size((int)ChildScanPtmLegendControl.ActualWidth, (int)ChildScanPtmLegendControl.ActualHeight);
                    legendCanvas.Measure(ptmLegendSize);
                    legendCanvas.Arrange(new Rect(ptmLegendSize));
                    legendCanvas.UpdateLayout();
                    ptmLegendLocationVector = (Vector)ChildScanPtmLegendControl.GetType().GetProperty("VisualOffset", BindingFlags.NonPublic | BindingFlags.Instance).GetValue(ChildScanPtmLegendControl);
                    ptmLegendLocationVector.X = PsmAnnotationGrid.ActualWidth - ChildScanPtmLegendControl.ActualWidth;
                }

                // If psm and GUI display have different number of matched ions, send in the refragmenter for exporting. 
                FragmentationReanalysisViewModel toPlotForExport = null;
                if (MetaDrawLogic.SpectrumAnnotation.SpectrumMatch.MatchedIons.Count != MetaDrawLogic.SpectrumAnnotation.MatchedFragmentIons.Count)
                {
                    toPlotForExport = FragmentationReanalysisViewModel;
                }

                MetaDrawLogic.ExportPlot(plotView, stationarySequenceCanvas, items, itemsControlSampleViewModel,
                    directoryPath, out errors, legendCanvas, ptmLegendLocationVector, toPlotForExport);
            }

            if (errors != null && errors.Any())
            {
                MessageBox.Show(errors.First());
            }
            else
            {
                MessageBoxHelper.Show(MetaDrawSettings.ExportType + "(s) exported to: " + directoryPath);
            }
        }

        private void ExportSpectrumLibraryButton_Click(object sender, RoutedEventArgs e)
        {
            if (dataGridScanNums.SelectedCells.Count == 0)
            {
                MessageBox.Show("Please select at least one scan to export");
                return;
            }

            List<SpectrumMatchFromTsv> psms = new();

            foreach (var cell in dataGridScanNums.SelectedItems)
            {
                var psm = (SpectrumMatchFromTsv)cell;
                psms.Add(psm);
            }

            string directoryPath = Path.Combine(Path.GetDirectoryName(MetaDrawLogic.SpectralMatchResultFilePaths.First()),
                "MetaDrawExport",    
                DateTime.Now.ToString("yyyy-MM-dd", CultureInfo.InvariantCulture));

            if(!Directory.Exists(directoryPath)) 
            { 
                Directory.CreateDirectory(directoryPath);
            }

            var libraryPath = Path.Combine(directoryPath, "spectrumLibrary.msp");

            // research all identifications with the newly matched ion types
            using (var sw = new StreamWriter(File.Create(libraryPath)))
            {
                foreach (var psm in psms)
                {
                    var oldIons = psm.MatchedIons;
                    ReplaceFragmentIonsOnPsmFromFragmentReanalysisViewModel(psm);
                    sw.WriteLine(psm.ToLibrarySpectrum().ToString());
                    psm.MatchedIons = oldIons;
                }
            }

            MessageBoxHelper.Show("Spectral Library exported to: " + libraryPath);
        }

        private void SequenceCoverageExportButton_Click(object sender, RoutedEventArgs e)
        {
            if (dataGridScanNums.SelectedItems.Count == 0 || dataGridScanNums.SelectedItems.Count > 1)
            {
                MessageBox.Show("Please select one psm to export sequence coverage");
                return;
            }

            string directoryPath = Path.Combine(Path.GetDirectoryName(MetaDrawLogic.SpectralMatchResultFilePaths.First()), "MetaDrawExport",
                    DateTime.Now.ToString("yyyy-MM-dd", CultureInfo.InvariantCulture));
            SpectrumMatchFromTsv psm = (SpectrumMatchFromTsv)dataGridScanNums.SelectedItem;
            MetaDrawLogic.ExportSequenceCoverage(sequenceText, map, directoryPath, psm);
            
            if (Directory.Exists(directoryPath))
            {
                MessageBoxHelper.Show(MetaDrawSettings.ExportType + " exported to: " + directoryPath);
            }
        }    

        private void SequenceAnnotationExportButton_Click(object sender, RoutedEventArgs e)
        {
            if (dataGridScanNums.SelectedItems.Count == 0 || dataGridScanNums.SelectedItems.Count > 1)
            {
                MessageBox.Show("Please select one psm to export sequence annotation");
                return;
            }

            string directoryPath = Path.Combine(Path.GetDirectoryName(MetaDrawLogic.SpectralMatchResultFilePaths.First()), "MetaDrawExport",
                    DateTime.Now.ToString("yyyy-MM-dd", CultureInfo.InvariantCulture));
            SpectrumMatchFromTsv psm = (SpectrumMatchFromTsv)dataGridScanNums.SelectedItem;

            int width = (int)SequenceAnnotationGrid.ActualWidth;
            MetaDrawLogic.ExportAnnotatedSequence(sequenceAnnotationCanvas, SequenceCoveragePtmLegendControl, psm, directoryPath, width);
            
            if (Directory.Exists(directoryPath))
            {
                MessageBoxHelper.Show(MetaDrawSettings.ExportType + " exported to: " + directoryPath);
            }
        }

        #region Data Visualization Tab

        private void SetUpPlots()
        {
            foreach (var plot in PlotModelStat.PlotNames)
            {
                plotTypes.Add(plot);
            }
        }

        private void loadFilesButtonStat_Click(object sender, RoutedEventArgs e)
        {
            // check for validity
            if (!MetaDrawLogic.SpectralMatchResultFilePaths.Any())
            {
                MessageBox.Show("Please add a search result file.");
                return;
            }

            (sender as Button).IsEnabled = false;
            selectPsmFileButtonStat.IsEnabled = false;
            resetPsmFileButtonStat.IsEnabled = false;
            prgsFeedStat.IsOpen = true;

            // load the PSMs
            this.prgsTextStat.Content = "Loading data...";
            MetaDrawLogic.LoadFiles(loadSpectra: false, loadPsms: true);

            PsmStatPlotFiles.Clear();
            foreach (var item in MetaDrawLogic.SpectralMatchesGroupedByFile)
            {
                PsmStatPlotFiles.Add(item.Key);
            }

            // done loading - restore controls
            this.prgsFeedStat.IsOpen = false;
            (sender as Button).IsEnabled = true;
            selectPsmFileButtonStat.IsEnabled = true;
            resetPsmFileButtonStat.IsEnabled = true;
        }

        /// <summary>
        /// Export for Data Visualization tab
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        private void CreatePlotPdf_Click(object sender, RoutedEventArgs e)
        {
            var selectedItem = plotsListBox.SelectedItem;

            if (selectedItem == null)
            {
                MessageBox.Show("Select a plot type to export!");
                return;
            }

            if (!MetaDrawLogic.SpectralMatchResultFilePaths.Any())
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
            var fileDirectory = Path.Combine(Path.GetDirectoryName(MetaDrawLogic.SpectralMatchResultFilePaths.First()), "MetaDrawExport",
                    DateTime.Now.ToString("yyyy-MM-dd", CultureInfo.InvariantCulture));
            var fileName = String.Concat(plotName, ".pdf");

            // update font sizes to exported PDF's size
            double tmpW = plotViewStat.Width;
            double tmpH = plotViewStat.Height;
            plotViewStat.Width = 1000;
            plotViewStat.Height = 700;
            plotViewStat.UpdateLayout();
            PlotViewStat_SizeChanged(plotViewStat, null);

            if (!Directory.Exists(fileDirectory))
            {
                Directory.CreateDirectory(fileDirectory);
            }

            using (Stream writePDF = File.Create(Path.Combine(fileDirectory, fileName)))
            {
                PdfExporter.Export(plotViewStat.Model, writePDF, 1000, 700);
            }
            plotViewStat.Width = tmpW;
            plotViewStat.Height = tmpH;
            MessageBoxHelper.Show(MetaDrawSettings.ExportType + " Created at " + Path.Combine(fileDirectory, fileName) + "!");
        }

        private async void PlotSelected(object sender, SelectionChangedEventArgs e)
        {
            var listview = sender as ListView;
            var plotName = listview.SelectedItem as string;

            if (MetaDrawLogic.FilteredListOfPsms.Count == 0)
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
            ObservableCollection<SpectrumMatchFromTsv> psms = new();
            Dictionary<string, ObservableCollection<SpectrumMatchFromTsv>> psmsBSF = new();
            foreach (string fileName in selectSourceFileListBox.SelectedItems)
            {
                psmsBSF.Add(fileName, new ObservableCollection<SpectrumMatchFromTsv>());
                foreach (SpectrumMatchFromTsv psm in MetaDrawLogic.SpectralMatchesGroupedByFile[fileName])
                {
                    if (!MetaDrawSettings.DisplayFilteredOnly)
                    {
                        psms.Add(psm);
                        psmsBSF[fileName].Add(psm);
                    }
                    else if (MetaDrawSettings.FilterAcceptsPsm(psm))
                    {
                        psms.Add(psm);
                        psmsBSF[fileName].Add(psm);
                    }
                }
            }

            PlotModelStat plot = null;
            try
            {
                plot = await Task.Run(() => new PlotModelStat(plotName, psms, psmsBSF));
            }
            catch (Exception ex)
            {
                MessageBox.Show($"An error occurred while generating the plot '{plotName}':\n{ex.Message}", "Plot Generation Error", MessageBoxButton.OK, MessageBoxImage.Error);
                return;
            }
            plotViewStat.DataContext = plot;
            PlotViewStat_SizeChanged(plotViewStat, null);
        }

        private void selectSourceFileListBox_SelectionChanged(object sender, EventArgs e)
        {
            // refreshes the plot using the new source file
            if (plotsListBox?.SelectedIndex > -1 && selectSourceFileListBox?.SelectedItems.Count != 0)
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
            if (plotsListBox.SelectedItem == null || !plotsListBox.SelectedItem.ToString().Equals("Histogram of PTM Spectral Counts"))
            {
                return;
            }
            OxyPlot.Wpf.PlotView plot = sender as OxyPlot.Wpf.PlotView;
            if (plot != null && plot.Model != null)
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


        private Point _dragStartPoint;
        private object _draggedData;
        private void selectSourceFileListBox_PreviewMouseLeftButtonDown(object sender, MouseButtonEventArgs e)
        {
            _dragStartPoint = e.GetPosition(null);
            _draggedData = GetDataFromListBoxItemUnderMouse(e.GetPosition(selectSourceFileListBox));
        }

        private void selectSourceFileListBox_PreviewMouseMove(object sender, MouseEventArgs e)
        {
            if (e.LeftButton == MouseButtonState.Pressed && _draggedData != null)
            {
                Point currentPosition = e.GetPosition(null);
                if (Math.Abs(currentPosition.X - _dragStartPoint.X) > SystemParameters.MinimumHorizontalDragDistance ||
                    Math.Abs(currentPosition.Y - _dragStartPoint.Y) > SystemParameters.MinimumVerticalDragDistance)
                {
                    DragDrop.DoDragDrop(selectSourceFileListBox, _draggedData, DragDropEffects.Move);
                    _draggedData = null;
                }
            }
        }

        private void selectSourceFileListBox_Drop(object sender, DragEventArgs e)
        {
            var droppedData = e.Data.GetData(typeof(string)) as string;
            var targetData = GetDataFromListBoxItemUnderMouse(e.GetPosition(selectSourceFileListBox));
            if (droppedData == null || targetData == null || droppedData == targetData)
                return;

            int removedIdx = PsmStatPlotFiles.IndexOf(droppedData);
            int targetIdx = PsmStatPlotFiles.IndexOf(targetData as string);
            if (removedIdx < 0 || targetIdx < 0)
                return;

            if (removedIdx < targetIdx)
            {
                PsmStatPlotFiles.Insert(targetIdx + 1, droppedData);
                PsmStatPlotFiles.RemoveAt(removedIdx);
            }
            else
            {
                int remIdx = removedIdx + 1;
                if (PsmStatPlotFiles.Count + 1 <= remIdx) return;
                PsmStatPlotFiles.Insert(targetIdx, droppedData);
                PsmStatPlotFiles.RemoveAt(remIdx);
            }
        }

        // Helper to get the data object from the ListBoxItem under the mouse
        private object GetDataFromListBoxItemUnderMouse(Point point)
        {
            var element = selectSourceFileListBox.InputHitTest(point) as DependencyObject;
            while (element != null && !(element is ListBoxItem))
                element = VisualTreeHelper.GetParent(element);

            return (element as ListBoxItem)?.DataContext;
        }

        private void DataVisualizationFilters_OnChecked(object sender, RoutedEventArgs e)
        {
            if (plotsListBox?.SelectedIndex > -1 && selectSourceFileListBox?.SelectedItems.Count != 0)
            {
                PlotSelected(plotsListBox, null);
            }
        }

        #endregion

        /// <summary>
        /// Redraws the Stationary Sequence whenever the scrolling sequence is scrolled
        /// </summary>
        /// <param name="sender"> 
        /// <param name="e"></param>
        private void wholeSequenceCoverageHorizontalScroll_Scroll(object sender, ScrollChangedEventArgs e)
        {
            SpectrumMatchFromTsv psm = (SpectrumMatchFromTsv)dataGridScanNums.SelectedItem;
            if (AmbiguousSequenceOptionBox.Items.Count > 1 && AmbiguousSequenceOptionBox.SelectedItem != null)
            {
                psm = (SpectrumMatchFromTsv)AmbiguousSequenceOptionBox.SelectedItem;

                // Draw the matched ions for the first ambiguous sequence only
                if (AmbiguousSequenceOptionBox.SelectedIndex == 0)
                {
                    MetaDrawSettings.DrawMatchedIons = true;
                }
                else
                {
                    MetaDrawSettings.DrawMatchedIons = false;
                }
            }
            SetSequenceDrawingPositionSettings();
            if (MetaDrawLogic.StationarySequence != null && !psm.FullSequence.Contains('|'))
            {
                DrawnSequence.DrawStationarySequence(psm, MetaDrawLogic.StationarySequence, 10);
            }
                
        }

        /// <summary>
        /// Redraws the Stationary and Scrollable sequences upon the selection of an Ambiguous PSM
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        private void AmbiguousSequenceOptionBox_SelectionChanged(object sender, SelectionChangedEventArgs e)
        {
            if (AmbiguousWarningTextBlocks.Visibility != Visibility.Collapsed)
            {
                AmbiguousWarningTextBlocks.Visibility = Visibility.Collapsed;
            }

            wholeSequenceCoverageHorizontalScroll.ScrollToLeftEnd();
            SpectrumMatchFromTsv psm = (SpectrumMatchFromTsv)dataGridScanNums.SelectedItem;
            if (AmbiguousSequenceOptionBox.Items.Count > 1 && AmbiguousSequenceOptionBox.SelectedItem != null)
            {
                psm = (SpectrumMatchFromTsv)AmbiguousSequenceOptionBox.SelectedItem;
                wholeSequenceCoverageHorizontalScroll.Visibility = Visibility.Visible;

                // Draw the matched ions for the first ambiguous sequence only
                if (AmbiguousSequenceOptionBox.SelectedIndex == 0)
                {
                    MetaDrawSettings.DrawMatchedIons = true;
                }
                else
                {
                    MetaDrawSettings.DrawMatchedIons = false;
                }
            }
            SetSequenceDrawingPositionSettings(true);
            object obj = new object();
            if (AmbiguousSequenceOptionBox.Items.Count > 0)
            {
                dataGridScanNums_SelectedCellsChanged(obj, null);
                MetaDrawLogic.DisplaySequences(stationarySequenceCanvas, scrollableSequenceCanvas,
                    sequenceAnnotationCanvas, psm);
            }
            else
            {
                AmbiguousSequenceOptionBox.Visibility = Visibility.Hidden;
            }
        }

        /// <summary>
        /// Method to set the MetaDrawSettings fields FirstAAOnScreen and NumberofAAonScreen to the current scrolling sequence position
        /// </summary>
        private void SetSequenceDrawingPositionSettings(bool reset = false)
        {
            if (dataGridScanNums.SelectedItem == null)
                return;

            // Get the total width of the sequence annotation area
            double totalWidth = SequenceAnnotationArea.ActualWidth;

            // Get the width of the description area or take a guess
            double descriptionWidth;
            if (plotView.Model != null)
            {
                var description =
                    plotView.ActualModel.Annotations.First(p =>
                        p is PlotTextAnnotation anno && anno.Text.Contains("\r\n")) as PlotTextAnnotation;
                descriptionWidth = -description!.X - 60;
            }
            else
            {
                descriptionWidth = 160;
            }

            // Define the offset (gap) you want between the sequence and the description
            double rightOffset = 0.0; 

            // Calculate the available width for the sequence
            double availableWidth = totalWidth - descriptionWidth - rightOffset;
            if (availableWidth < 0) 
                availableWidth = 0;

            // Use the scrolling sequence offset to determine where to start
            double offset = wholeSequenceCoverageHorizontalScroll.HorizontalOffset;
            if (reset)
                offset = 0;

            SpectrumMatchFromTsv psm = (SpectrumMatchFromTsv)dataGridScanNums.SelectedItem;
            if (AmbiguousSequenceOptionBox.Items.Count > 1 && AmbiguousSequenceOptionBox.SelectedItem != null)
                psm = (SpectrumMatchFromTsv)AmbiguousSequenceOptionBox.SelectedItem;

            int lettersOnScreen = (int)Math.Round((availableWidth) / MetaDrawSettings.AnnotatedSequenceTextSpacing, 0);
            int firstLetterOnScreen = (int)Math.Round((offset) / MetaDrawSettings.AnnotatedSequenceTextSpacing, 0);
            if ((firstLetterOnScreen + lettersOnScreen) > psm.BaseSeq.Length)
            {
                lettersOnScreen = psm.BaseSeq.Length - firstLetterOnScreen;
            }
            MetaDrawSettings.FirstAAonScreenIndex = firstLetterOnScreen;
            MetaDrawSettings.NumberOfAAOnScreen = lettersOnScreen;
        }

        /// <summary>
        /// Fires the command in the PtmLegend to decrease the residues per segment by one
        /// </summary>
        /// <param name="sender">Button Object</param>
        /// <param name="e"></param>
        private void residuesPerSegmentcmdDown_Click(object sender, RoutedEventArgs e)
        {
            if (PtmLegend.ResiduesPerSegment == 1)
            {
                MessageBox.Show("Value must be greater than 0");
                return;
            }

            PtmLegend.DecreaseResiduesPerSegment();
            SpectrumMatchFromTsv psm = (SpectrumMatchFromTsv)dataGridScanNums.SelectedItem;
            MetaDrawLogic.DisplaySequences(null, null, sequenceAnnotationCanvas, psm);
        }

        /// <summary>
        /// Fires the command in the PtmLegend to increase the residues per segment of the annotated sequence by one
        /// </summary>
        /// <param name="sender">Button Object</param>
        /// <param name="e"></param>
        private void residuesPerSegmentcmdUp_Click(object sender, RoutedEventArgs e)
        {
            PtmLegend.IncreaseResiduesPerSegment();
            SpectrumMatchFromTsv psm = (SpectrumMatchFromTsv)dataGridScanNums.SelectedItem;
            MetaDrawLogic.DisplaySequences(null, null, sequenceAnnotationCanvas, psm);
        }

        /// <summary>
        /// Fires the command in the PtmLegend to decrease the segments per row of the annotated sequence by one
        /// </summary>
        /// <param name="sender">Button Object</param>
        /// <param name="e"></param>
        private void segmentsPerRowcmdDown_Click(object sender, RoutedEventArgs e)
        {
            if (PtmLegend.SegmentsPerRow == 1)
            {
                MessageBox.Show("Value must be greater than 0");
                return;
            }

            PtmLegend.DecreaseSegmentsPerRow();
            SpectrumMatchFromTsv psm = (SpectrumMatchFromTsv)dataGridScanNums.SelectedItem;
            MetaDrawLogic.DisplaySequences(null, null, sequenceAnnotationCanvas, psm);
        }

        /// <summary>
        /// Fires the command in the PtmLegend to increase the segments per row of the annotated sequence by one
        /// </summary>
        /// <param name="sender">Button Object</param>
        /// <param name="e"></param>
        private void segmentsPerRowcmdUp_Click(object sender, RoutedEventArgs e)
        {
            PtmLegend.IncreaseSegmentsPerRow();
            SpectrumMatchFromTsv psm = (SpectrumMatchFromTsv)dataGridScanNums.SelectedItem;
            MetaDrawLogic.DisplaySequences(null, null, sequenceAnnotationCanvas, psm);
        }


        private void MetaDrawTabControl_OnSelectionChanged(object sender, SelectionChangedEventArgs e)
        {
            SpectrumMatchFromTsv selectedPsm = (SpectrumMatchFromTsv)dataGridScanNums.SelectedItem;

            if (e.OriginalSource is not TabControl) // only clicking on different MetaDrawTabs will trigger this event
                return;

            // switch from chimera to other views
            if (e.RemovedItems.Count > 0 && ((TabItem)e.RemovedItems[0]).Name == "ChimeraScanPlot")
            {
                MetaDrawLogic.FilterPsms();
                ClearPresentationArea();

                // reselect what was selected
                if (selectedPsm != null && MetaDrawLogic.FilteredListOfPsms.Contains(selectedPsm))
                {
                    int psmIndex = MetaDrawLogic.FilteredListOfPsms.IndexOf(selectedPsm);
                    dataGridScanNums.SelectedIndex = psmIndex;
                    dataGridScanNums_SelectedCellsChanged(new object(), null);
                }
            }
        }

        /// <summary>
        /// Clears and resets the presentation area
        /// </summary>
        private void ClearPresentationArea()
        {
            DrawnSequence.ClearCanvas(scrollableSequenceCanvas);
            DrawnSequence.ClearCanvas(stationarySequenceCanvas);
            DrawnSequence.ClearCanvas(map);
            DrawnSequence.ClearCanvas(sequenceText);
            DrawnSequence.ClearCanvas(sequenceAnnotationCanvas);
            GrayBox.Opacity = 0;
            wholeSequenceCoverageHorizontalScroll.Visibility = Visibility.Collapsed;
            AmbiguousSequenceOptionBox.Items.Clear();
            plotView.Visibility = Visibility.Hidden;

            if (PtmLegend != null)
                PtmLegend.Visibility = false;
        }

        /// <summary>
        /// Enables and disables the buttons on the main MetaDraw view
        /// </summary>
        /// <param name="value">true = enabled, false = disable</param>
        private void ToggleButtonsEnabled(bool value)
        {
            loadFiles.IsEnabled = value;
            selectSpectraFileButton.IsEnabled = value;
            selectPsmFileButton.IsEnabled = value;
            selectSpecLibraryButton.IsEnabled = value;
            resetPsmFileButton.IsEnabled = value;
            resetSpectraFileButton.IsEnabled = value;
            resetSpectraFileButton.IsEnabled = value;
            exportPdfs.IsEnabled = value;
            exportSpectrumLibrary.IsEnabled = value;
        }

        /// <summary>
        /// Method to fire the plotting method with new fragment ions
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        /// <exception cref="NotImplementedException"></exception>
        internal void SearchWithNewIons_OnClick(object sender, RoutedEventArgs e)
        {
            // find currently selected psm
            var psm = dataGridScanNums.SelectedItem as SpectrumMatchFromTsv;
            if (psm is null)
                return;
            
            // replace the ions and replot
            var oldIons = psm.MatchedIons;
            ReplaceFragmentIonsOnPsmFromFragmentReanalysisViewModel(psm);
            dataGridScanNums.SelectedItem = psm;
            dataGridScanNums_SelectedCellsChanged(FragmentationReanalysisViewModel, null);

            // put the old ions back
            psm.MatchedIons = oldIons;
        }

        /// <summary>
        /// Replaces matched fragment ions on a psm with new ion types after a quick search
        /// </summary>
        /// <param name="psm"></param>
        private void ReplaceFragmentIonsOnPsmFromFragmentReanalysisViewModel(SpectrumMatchFromTsv psm)
        {
            var scan = MetaDrawLogic.GetMs2ScanFromPsm(psm);
            var newIons = FragmentationReanalysisViewModel.MatchIonsWithNewTypes(scan, psm);
            psm.MatchedIons = newIons;
        }

        private void MetaDraw_OnClosing(object sender, CancelEventArgs e)
        {
            MetaDrawLogic.CleanUpResources();
        }
        
        #region Fragment Plot Click Effects 

        // Copy the entire m/z spectrum (all peaks)
        private void CopyMzSpectrum_Click(object sender, RoutedEventArgs e)
        {
            var plot = MetaDrawLogic.SpectrumAnnotation;
            if (plot?.Scan == null)
                return;

            var sb = new StringBuilder();
            var mzs = plot.Scan.MassSpectrum.XArray;
            var intensities = plot.Scan.MassSpectrum.YArray;
            for (int i = 0; i < mzs.Length; i++)
                sb.AppendLine($"{mzs[i]:F6}\t{intensities[i]:F6}");

            Clipboard.SetText(sb.ToString());
        }

        // Copy only annotated peaks (matched ions)
        private void CopyAnnotatedMzSpectrum_Click(object sender, RoutedEventArgs e)
        {
            var plot = MetaDrawLogic.SpectrumAnnotation;
            if (plot?.Scan == null || plot.SpectrumMatch == null)
                return;

            var matched = plot.MatchedFragmentIons;

            var sb = new StringBuilder();
            foreach (var ion in matched)
                sb.AppendLine($"{ion.Mz:F6}\t{ion.Intensity:F0}");
            Clipboard.SetText(sb.ToString());
        }

        // Copy matched ions with details
        private void CopyMatchedIons_Click(object sender, RoutedEventArgs e)
        {
            var plot = MetaDrawLogic.SpectrumAnnotation;
            if (plot?.SpectrumMatch == null)
                return;

            var matched = plot.MatchedFragmentIons;

            var sb = new StringBuilder();
            sb.AppendLine("Annotation\tm/z\tIntensity\tType\tFragmentNumber");
            foreach (var ion in matched)
            {
                sb.AppendLine($"{ion.Annotation}\t{ion.Mz:F6}\t{ion.Intensity:F0}\t{ion.NeutralTheoreticalProduct.ProductType}\t{ion.NeutralTheoreticalProduct.FragmentNumber}");
            }

            Clipboard.SetText(sb.ToString());
        }

        #endregion
    }
}
