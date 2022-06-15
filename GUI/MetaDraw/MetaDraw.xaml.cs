using EngineLayer;
using GuiFunctions;
using GuiFunctions.MetaDraw;
using Nett;
using OxyPlot;
using Proteomics.Fragmentation;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Data;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;

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
        private static List<string> AcceptedSpectraFormats = new List<string> { ".mzml", ".raw", ".mgf" };
        private static List<string> AcceptedResultsFormats = new List<string> { ".psmtsv", ".tsv" };
        private static List<string> AcceptedSpectralLibraryFormats = new List<string> { ".msp" };

        public MetaDraw()
        {
            UsefulProteomicsDatabases.Loaders.LoadElements();

            InitializeComponent();

            MetaDrawLogic = new MetaDrawLogic();
            BindingOperations.EnableCollectionSynchronization(MetaDrawLogic.PsmResultFilePaths, MetaDrawLogic.ThreadLocker);
            BindingOperations.EnableCollectionSynchronization(MetaDrawLogic.SpectraFilePaths, MetaDrawLogic.ThreadLocker);
            BindingOperations.EnableCollectionSynchronization(MetaDrawLogic.FilteredListOfPsms, MetaDrawLogic.ThreadLocker);
            BindingOperations.EnableCollectionSynchronization(MetaDrawLogic.PsmsGroupedByFile, MetaDrawLogic.ThreadLocker);

            itemsControlSampleViewModel = new ParentChildScanPlotsView();
            ParentChildScanViewPlots.DataContext = itemsControlSampleViewModel;

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

            // checks to see if default settings have been saved, and loads them for the first opening of the window
            MetaDrawSettingsSnapshot settings = null;
            string settingsPath = Path.Combine(GlobalVariables.DataDir, "DefaultParameters", @"MetaDrawSettingsDefault.toml");
            if (File.Exists(settingsPath))
            {
                settings = Toml.ReadFile<MetaDrawSettingsSnapshot>(settingsPath);
                MetaDrawSettings.LoadSettings(settings);
            }
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
                if (!MetaDrawLogic.PsmResultFilePaths.Contains(filePath))
                {
                    MetaDrawLogic.PsmResultFilePaths.Add(filePath);

                    if (MetaDrawLogic.PsmResultFilePaths.Count == 1)
                    {
                        psmFileNameLabel.Text = filePath;
                        psmFileNameLabelStat.Text = filePath;
                    }
                    else
                    {
                        psmFileNameLabel.Text = "[Mouse over to view files]";
                        psmFileNameLabelStat.Text = "[Mouse over to view files]";
                    }

                    psmFileNameLabel.ToolTip = string.Join("\n", MetaDrawLogic.PsmResultFilePaths);
                    resetPsmFileButton.IsEnabled = true;

                    psmFileNameLabelStat.ToolTip = string.Join("\n", MetaDrawLogic.PsmResultFilePaths);
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
            else
            {
                MessageBox.Show("Cannot read file type: " + theExtension);
            }
        }

        /// <summary>
        /// Event triggers when a different cell is selected in the PSM data grid
        /// </summary>
        private void dataGridScanNums_SelectedCellsChanged(object sender, SelectedCellsChangedEventArgs e)
        {
            if (dataGridScanNums.SelectedItem == null || sender == null)
            {
                return;
            }

            wholeSequenceCoverageHorizontalScroll.ScrollToLeftEnd();
            AmbiguousSequenceOptionBox.Items.Clear();
            plotView.Visibility = Visibility.Visible;
            PsmFromTsv psm = (PsmFromTsv)dataGridScanNums.SelectedItem;
            SetSequenceDrawingPositionSettings(true);
            // If psm is ambiguous, split into several psm objects and add each as an option to see
            if (psm.FullSequence.Contains('|'))
            {
                MetaDrawSettings.DrawMatchedIons = false;
                DrawnSequence.ClearCanvas(scrollableSequenceCanvas);
                DrawnSequence.ClearCanvas(stationarySequenceCanvas);
                GrayBox.Opacity = 0;
                var fullSeqs = psm.FullSequence.Split('|');
                foreach (var fullSeq in fullSeqs)
                {
                    PsmFromTsv oneAmbiguousPsm = new(psm, fullSeq);
                    AmbiguousSequenceOptionBox.Items.Add(oneAmbiguousPsm);
                }
                AmbiguousWarningTextBlocks.Visibility = Visibility.Visible;
                AmbiguousSequenceOptionBox.Visibility = Visibility.Visible;
                wholeSequenceCoverageHorizontalScroll.Visibility = Visibility.Collapsed;
            }
            else
            {
                MetaDrawSettings.DrawMatchedIons = true;
                AmbiguousWarningTextBlocks.Visibility = Visibility.Collapsed;
                AmbiguousSequenceOptionBox.Visibility = Visibility.Collapsed;
                wholeSequenceCoverageHorizontalScroll.Visibility = Visibility.Visible;
                if (psm.BaseSeq.Length > MetaDrawSettings.NumberOfAAOnScreen)
                    GrayBox.Opacity = 0.7;
                else
                    GrayBox.Opacity = 0;
            }

            // draw the annotated spectrum
            MetaDrawLogic.DisplaySequences(stationarySequenceCanvas, scrollableSequenceCanvas, psm);
            MetaDrawLogic.DisplaySpectrumMatch(plotView, psm, itemsControlSampleViewModel, out var errors);

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
                    propertyView.Rows.Add(temp[i].Name, temp[i].GetValue(psm, null));
                }
            }
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
                DrawnSequence.ClearCanvas(MetaDrawLogic.ScrollableSequence.SequenceDrawingCanvas);
                DrawnSequence.ClearCanvas(MetaDrawLogic.StationarySequence.SequenceDrawingCanvas);
                plotView.Visibility = Visibility.Hidden;
                MetaDrawLogic.FilteredListOfPsms.Clear();
            }
        }

        private void OnClosing(object sender, CancelEventArgs e)
        {
            MetaDrawLogic.CleanUpResources();
        }

        private void settings_Click(object sender, RoutedEventArgs e)
        {
            
            // save current selected PSM
            var selectedItem = dataGridScanNums.SelectedItem;
            var settingsWindow = new MetaDrawSettingsWindow();
            var result = settingsWindow.ShowDialog();

            // re-select selected PSM
            

            if (result == true)
            {
                // refresh chart
                dataGridScanNums_SelectedCellsChanged(null, null);

                // filter based on new settings
                MetaDrawLogic.FilterPsms();
            }

            if (selectedItem != null)
            {
                dataGridScanNums.SelectedItem = selectedItem;
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

            if (!MetaDrawLogic.PsmResultFilePaths.Any())
            {
                MessageBox.Show("Please add a search result file.");
                return;
            }

            // load the spectra file
            (sender as Button).IsEnabled = false;
            selectSpectraFileButton.IsEnabled = false;
            selectPsmFileButton.IsEnabled = false;
            selectSpecLibraryButton.IsEnabled = false;
 
            prgsFeed.IsOpen = true;
            prgsText.Content = "Loading data...";

            // Add EventHandlers for popup click-in/click-out behaviour
            Deactivated += new EventHandler(prgsFeed_Deactivator);
            Activated += new EventHandler(prgsFeed_Reactivator);

            var slowProcess = Task<List<string>>.Factory.StartNew(() => MetaDrawLogic.LoadFiles(loadSpectra: true, loadPsms: true));
            await slowProcess;
            var errors = slowProcess.Result;

            if (errors.Any())
            {
                string errorList = string.Join("\n", errors);
                MessageBox.Show(errorList);
            }

            PsmStatPlotFiles.Clear();
            foreach (var item in MetaDrawLogic.PsmsGroupedByFile)
            {
                PsmStatPlotFiles.Add(item.Key);
            }

            // done loading - restore controls
            this.prgsFeed.IsOpen = false;

            // Remove added EventHandlers
            Deactivated -= new EventHandler(prgsFeed_Deactivator);
            Activated -= new EventHandler(prgsFeed_Reactivator);

            (sender as Button).IsEnabled = true;
            selectSpectraFileButton.IsEnabled = true;
            selectPsmFileButton.IsEnabled = true;
            selectSpecLibraryButton.IsEnabled = true;
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

        private void PDFButton_Click(object sender, RoutedEventArgs e)
        {
            if (dataGridScanNums.SelectedCells.Count == 0)
            {
                MessageBox.Show("Please select at least one scan to export");
                return;
            }

            List<PsmFromTsv> items = new List<PsmFromTsv>();

            foreach (var cell in dataGridScanNums.SelectedItems)
            {
                var psm = (PsmFromTsv)cell;
                items.Add(psm);
            }

            string directoryPath = Path.Combine(Path.GetDirectoryName(MetaDrawLogic.PsmResultFilePaths.First()), "MetaDrawExport",
                    DateTime.Now.ToString("yyyy-MM-dd-HH-mm-ss", CultureInfo.InvariantCulture));

            MetaDrawLogic.ExportToPdf(plotView, stationarySequenceCanvas, items, itemsControlSampleViewModel, directoryPath, out var errors) ;

            if (errors.Any())
            {
                MessageBox.Show(errors.First());
            }
            else
            {
                MessageBox.Show("PDFs exported to: " + directoryPath);
            }
        }

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
            if (!MetaDrawLogic.PsmResultFilePaths.Any())
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
            foreach (var item in MetaDrawLogic.PsmsGroupedByFile)
            {
                PsmStatPlotFiles.Add(item.Key);
            }

            // done loading - restore controls
            this.prgsFeedStat.IsOpen = false;
            (sender as Button).IsEnabled = true;
            selectPsmFileButtonStat.IsEnabled = true;
            resetPsmFileButtonStat.IsEnabled = true;
        }

        private void CreatePlotPdf_Click(object sender, RoutedEventArgs e)
        {
            var selectedItem = plotsListBox.SelectedItem;

            if (selectedItem == null)
            {
                MessageBox.Show("Select a plot type to export!");
                return;
            }

            if (!MetaDrawLogic.PsmResultFilePaths.Any())
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
            var fileDirectory = Directory.GetParent(MetaDrawLogic.PsmResultFilePaths.First()).ToString();
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
            ObservableCollection<PsmFromTsv> psms = new ObservableCollection<PsmFromTsv>();
            Dictionary<string, ObservableCollection<PsmFromTsv>> psmsBSF = new Dictionary<string, ObservableCollection<PsmFromTsv>>();
            foreach (string fileName in selectSourceFileListBox.SelectedItems)
            {
                psmsBSF.Add(fileName, MetaDrawLogic.PsmsGroupedByFile[fileName]);
                foreach (PsmFromTsv psm in MetaDrawLogic.PsmsGroupedByFile[fileName])
                {
                    psms.Add(psm);
                }
            }
            PlotModelStat plot = await Task.Run(() => new PlotModelStat(plotName, psms, psmsBSF));
            plotViewStat.DataContext = plot;
            PlotViewStat_SizeChanged(plotViewStat, null);
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

        private void AnnotationSizeChanged(object sender, SizeChangedEventArgs e)
        {
            mapViewer.Height = .8 * SequenceAnnotationGrid.ActualHeight;
            mapViewer.Width = .99 * SequenceAnnotationGrid.ActualWidth;
        }

        /// <summary>
        /// Redraws the Stationary Sequence whenever the scrolling sequence is scrolled
        /// </summary>
        /// <param name="sender"> 
        /// <param name="e"></param>
        private void wholeSequenceCoverageHorizontalScroll_Scroll(object sender, ScrollChangedEventArgs e)
        {
            PsmFromTsv psm = (PsmFromTsv)dataGridScanNums.SelectedItem;
            if (AmbiguousSequenceOptionBox.Items.Count > 1 && AmbiguousSequenceOptionBox.SelectedItem != null)
            {
                psm = (PsmFromTsv)AmbiguousSequenceOptionBox.SelectedItem;

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
            if (MetaDrawLogic.StationarySequence != null)
                DrawnSequence.DrawStationarySequence(psm, MetaDrawLogic.StationarySequence);
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
            PsmFromTsv psm = (PsmFromTsv)dataGridScanNums.SelectedItem;
            if (AmbiguousSequenceOptionBox.Items.Count > 1 && AmbiguousSequenceOptionBox.SelectedItem != null)
            {
                psm = (PsmFromTsv)AmbiguousSequenceOptionBox.SelectedItem;
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
            MetaDrawLogic.DisplaySequences(stationarySequenceCanvas, scrollableSequenceCanvas, psm);
        }

        /// <summary>
        /// Method to set the MetaDrawSettings fields FirstAAOnScreen and NumberofAAonScreen to the current scrolling sequence position
        /// </summary>
        private void SetSequenceDrawingPositionSettings(bool reset = false)
        {
            double width = SequenceAnnotationArea.ActualWidth;
            double offset = wholeSequenceCoverageHorizontalScroll.HorizontalOffset;
            if (reset)
            {
                offset = 0;
            }
            PsmFromTsv psm = (PsmFromTsv)dataGridScanNums.SelectedItem;
            if (AmbiguousSequenceOptionBox.Items.Count > 1 && AmbiguousSequenceOptionBox.SelectedItem != null)
            {
                psm = (PsmFromTsv)AmbiguousSequenceOptionBox.SelectedItem;
            }

            int lettersOnScreen = (int)Math.Round((width - 10) / MetaDrawSettings.AnnotatedSequenceTextSpacing, 0);
            int firstLetterOnScreen = (int)Math.Round((offset) / MetaDrawSettings.AnnotatedSequenceTextSpacing, 0);
            if ((firstLetterOnScreen + lettersOnScreen) > psm.BaseSeq.Length)
            {
                lettersOnScreen = psm.BaseSeq.Length - firstLetterOnScreen;
            }
            MetaDrawSettings.FirstAAonScreenIndex = firstLetterOnScreen;
            MetaDrawSettings.NumberOfAAOnScreen = lettersOnScreen;
        }
    }
}