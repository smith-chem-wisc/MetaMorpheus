using EngineLayer;
using OxyPlot;
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
            ParentScanView.Visibility = Visibility.Collapsed;

            PsmStatPlotFiles = new ObservableCollection<string>();
            selectSourceFileListBox.DataContext = PsmStatPlotFiles;
            plotTypes = new ObservableCollection<string>();
            SetUpPlots();
            plotsListBox.ItemsSource = plotTypes;
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
            if (dataGridScanNums.SelectedItem == null)
            {
                return;
            }

            PsmFromTsv psm = (PsmFromTsv)dataGridScanNums.SelectedItem;

            // draw the PSM
            MetaDrawLogic.DisplaySpectrumMatch(plotView, canvas, psm, itemsControlSampleViewModel, out var errors);

            if (psm.ChildScanMatchedIons != null)
            {
                ParentChildScanView.Visibility = Visibility.Visible;
                ParentScanView.Visibility = Visibility.Visible;
            }
            else
            {
                ParentChildScanView.Visibility = Visibility.Collapsed;
                ParentScanView.Visibility = Visibility.Collapsed;
            }

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

        private void resetFilesButton_Click(object sender, RoutedEventArgs e)
        {
            MetaDrawLogic.CleanUpResources();
            spectraFileNameLabel.Text = "None Selected";
            psmFileNameLabel.Text = "None Selected";
        }

        private void OnClosing(object sender, CancelEventArgs e)
        {
            MetaDrawLogic.CleanUpResources();
        }

        private void settings_Click(object sender, RoutedEventArgs e)
        {
            var settingsWindow = new MetaDrawSettingsWindow();
            var result = settingsWindow.ShowDialog();

            // save current selected PSM
            var selectedItem = dataGridScanNums.SelectedItem;

            if (result == true)
            {
                // refresh chart
                dataGridScanNums_SelectedCellsChanged(null, null);

                // filter based on new settings
                MetaDrawLogic.FilterPsms();
            }

            // re-select selected PSM
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
            resetSpectraFileButton.IsEnabled = false;
            resetPsmFileButton.IsEnabled = false;
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
            resetSpectraFileButton.IsEnabled = true;
            resetPsmFileButton.IsEnabled = true;
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

            MetaDrawLogic.ExportToPdf(plotView, canvas, items, itemsControlSampleViewModel, directoryPath, out var errors);

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
    }
}