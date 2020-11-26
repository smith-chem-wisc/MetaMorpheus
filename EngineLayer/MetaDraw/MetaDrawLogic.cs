using IO.Mgf;
using IO.MzML;
using IO.ThermoRawFileReader;
using MassSpectrometry;
using mzPlot;
using OxyPlot;
using OxyPlot.Wpf;
using Proteomics.Fragmentation;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.ComponentModel;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Text;
using System.Text.RegularExpressions;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Media;

namespace EngineLayer
{
    public class MetaDrawLogic
    {
        public ObservableCollection<string> PsmResultFilePaths { get; private set; }
        public ObservableCollection<string> SpectraFilePaths { get; private set; }
        public ObservableCollection<string> SpectralLibraryPaths { get; private set; }
        public ObservableCollection<PsmFromTsv> FilteredListOfPsms { get; private set; } // filtered list of PSMs after q-value filter, etc.
        public Dictionary<string, ObservableCollection<PsmFromTsv>> PsmsGroupedByFile { get; private set; }
        public object ThreadLocker;
        public ICollectionView PeptideSpectralMatchesView;

        private List<PsmFromTsv> AllPsms; // all loaded PSMs
        private Dictionary<string, DynamicDataConnection> MsDataFiles; // key is file name without extension
        private List<PeptideSpectrumMatchPlot> CurrentlyDisplayedPlots;
        private Regex illegalInFileName = new Regex(@"[\\/:*?""<>|]");
        private SpectralLibrary SpectralLibrary;

        public MetaDrawLogic()
        {
            PsmResultFilePaths = new ObservableCollection<string>();
            SpectraFilePaths = new ObservableCollection<string>();
            SpectralLibraryPaths = new ObservableCollection<string>();
            FilteredListOfPsms = new ObservableCollection<PsmFromTsv>();
            PsmsGroupedByFile = new Dictionary<string, ObservableCollection<PsmFromTsv>>();
            AllPsms = new List<PsmFromTsv>();
            MsDataFiles = new Dictionary<string, DynamicDataConnection>();
            PeptideSpectralMatchesView = CollectionViewSource.GetDefaultView(FilteredListOfPsms);
            ThreadLocker = new object();
            CurrentlyDisplayedPlots = new List<PeptideSpectrumMatchPlot>();
        }

        public List<string> LoadFiles(bool loadSpectra, bool loadPsms)
        {
            List<string> errors = new List<string>();

            lock (ThreadLocker)
            {
                FilteredListOfPsms.Clear();
                PsmsGroupedByFile.Clear();
                AllPsms.Clear();
                MsDataFiles.Clear();
            }

            // load MS data files
            if (loadSpectra)
            {
                LoadSpectraFiles(out var errors1);
                errors.AddRange(errors1);
            }

            // load PSMs
            if (loadPsms)
            {
                LoadPsms(out var errors2, loadSpectra);
                errors.AddRange(errors2);
            }

            // load spectral libraries
            LoadSpectralLibraries(out var errors3);
            errors.AddRange(errors3);

            return errors;
        }

        public void DisplaySpectrumMatch(PlotView plotView, Canvas canvas, PsmFromTsv psm, ParentChildScanPlotsView parentChildScanPlotsView, out List<string> errors)
        {
            errors = null;

            // clear old parent/child scans
            parentChildScanPlotsView.Plots.Clear();
            CurrentlyDisplayedPlots.Clear();

            // get the scan
            if (!MsDataFiles.TryGetValue(psm.FileNameWithoutExtension, out DynamicDataConnection spectraFile))
            {
                errors = new List<string>();
                errors.Add("The spectra file could not be found for this PSM: " + psm.FileNameWithoutExtension);
                return;
            }

            MsDataScan scan = spectraFile.GetOneBasedScanFromDynamicConnection(psm.Ms2ScanNumber);

            LibrarySpectrum librarySpectrum = null;

            // plot the annotated spectrum match
            PeptideSpectrumMatchPlot plot;
            if (psm.BetaPeptideBaseSequence == null)
            {
                // get the library spectrum if relevant
                if (SpectralLibrary != null)
                {
                    SpectralLibrary.TryGetSpectrum(psm.FullSequence, psm.PrecursorCharge, out var librarySpectrum1);
                    librarySpectrum = librarySpectrum1;
                }

                plot = new PeptideSpectrumMatchPlot(plotView, canvas, psm, scan, psm.MatchedIons, librarySpectrum: librarySpectrum);
            }
            else
            {
                plot = new CrosslinkSpectrumMatchPlot(plotView, canvas, psm, scan);
            }

            CurrentlyDisplayedPlots.Add(plot);

            // plot parent/child scans
            if (psm.ChildScanMatchedIons != null)
            {
                // draw parent scan
                string parentAnnotation = "Scan: " + scan.OneBasedScanNumber
                        + " Dissociation Type: " + scan.DissociationType
                        + " MsOrder: " + scan.MsnOrder
                        + " Selected Mz: " + scan.SelectedIonMZ.Value.ToString("0.##")
                        + " Retention Time: " + scan.RetentionTime.ToString("0.##");

                var parentPlotView = new PlotView(); // placeholder
                var parentCanvas = new Canvas();
                var item = new ParentChildScanPlotTemplate()
                {
                    Plot = new PeptideSpectrumMatchPlot(parentPlotView, parentCanvas, psm, scan, psm.MatchedIons),
                    SpectrumLabel = parentAnnotation,
                    TheCanvas = parentCanvas
                };

                parentChildScanPlotsView.Plots.Add(item);

                // remove model from placeholder (the model can only be referenced by 1 plotview at a time)
                parentPlotView.Model = null;

                // draw child scans
                HashSet<int> scansDrawn = new HashSet<int>();
                var allChildScanMatchedIons = psm.ChildScanMatchedIons;

                if (psm.BetaPeptideChildScanMatchedIons != null)
                {
                    allChildScanMatchedIons = allChildScanMatchedIons.Concat(psm.BetaPeptideChildScanMatchedIons)
                        .GroupBy(p => p.Key)
                        .ToDictionary(p => p.Key, q => q.SelectMany(p => p.Value).ToList());
                }

                foreach (var childScanMatchedIons in allChildScanMatchedIons)
                {
                    int childScanNumber = childScanMatchedIons.Key;

                    if (scansDrawn.Contains(childScanNumber))
                    {
                        continue;
                    }

                    scansDrawn.Add(childScanNumber);

                    List<MatchedFragmentIon> matchedIons = childScanMatchedIons.Value;

                    MsDataScan childScan = spectraFile.GetOneBasedScanFromDynamicConnection(childScanNumber);

                    string childAnnotation = "Scan: " + childScan.OneBasedScanNumber
                        + " Dissociation Type: " + childScan.DissociationType
                        + " MsOrder: " + childScan.MsnOrder
                        + " Selected Mz: " + childScan.SelectedIonMZ.Value.ToString("0.##")
                        + " RetentionTime: " + childScan.RetentionTime.ToString("0.##");

                    Canvas childCanvas = new Canvas();
                    PlotView childPlotView = new PlotView(); // placeholder

                    // make the plot
                    var childPlot = new PeptideSpectrumMatchPlot(childPlotView, childCanvas, psm, childScan, matchedIons, annotateProperties: false);
                    childPlot.Model.Title = null;
                    childPlot.Model.Subtitle = null;

                    item = new ParentChildScanPlotTemplate() { Plot = childPlot, SpectrumLabel = childAnnotation, TheCanvas = childCanvas };

                    // remove model from placeholder (the model can only be referenced by 1 plotview at a time)
                    childPlotView.Model = null;

                    parentChildScanPlotsView.Plots.Add(item);

                    CurrentlyDisplayedPlots.Add(childPlot);
                }
            }
        }

        public void ExportToPdf(PlotView plotView, Canvas canvas, List<PsmFromTsv> spectrumMatches, ParentChildScanPlotsView parentChildScanPlotsView, string directory, out List<string> errors)
        {
            errors = new List<string>();

            if (!Directory.Exists(directory))
            {
                Directory.CreateDirectory(directory);
            }

            foreach (var psm in spectrumMatches)
            {
                DisplaySpectrumMatch(plotView, canvas, psm, parentChildScanPlotsView, out var displayErrors);

                if (displayErrors != null)
                {
                    errors.AddRange(displayErrors);
                }

                string sequence = illegalInFileName.Replace(psm.FullSequence, string.Empty);

                if (sequence.Length > 30)
                {
                    sequence = sequence.Substring(0, 30);
                }

                foreach (var plot in CurrentlyDisplayedPlots)
                {
                    string filePath = Path.Combine(directory, plot.Scan.OneBasedScanNumber + "_" + sequence + ".pdf");

                    int i = 2;
                    while (File.Exists(filePath))
                    {
                        filePath = Path.Combine(directory, plot.Scan.OneBasedScanNumber + "_" + sequence + "_" + i + ".pdf");
                        i++;
                    }

                    plot.ExportToPdf(filePath, plotView.ActualWidth, plotView.ActualHeight);
                }
            }

            DisplaySpectrumMatch(plotView, canvas, spectrumMatches.First(), parentChildScanPlotsView, out var moreDisplayErrors);
        }

        public void FilterPsms()
        {
            lock (ThreadLocker)
            {
                FilteredListOfPsms.Clear();

                foreach (var psm in AllPsms.Where(p => MetaDrawSettings.FilterAcceptsPsm(p)))
                {
                    FilteredListOfPsms.Add(psm);
                }
            }
        }

        public void FilterPsmsByString(string searchString)
        {
            if (searchString == "")
            {
                PeptideSpectralMatchesView.Filter = null;
            }
            else
            {
                PeptideSpectralMatchesView.Filter = obj =>
                {
                    PsmFromTsv psm = obj as PsmFromTsv;
                    return ((psm.Ms2ScanNumber.ToString()).StartsWith(searchString) || psm.FullSequence.ToUpper().Contains(searchString.ToUpper()));
                };
            }
        }

        public void CleanUpResources()
        {
            lock (ThreadLocker)
            {
                AllPsms.Clear();
                FilteredListOfPsms.Clear();
                PsmResultFilePaths.Clear();
                SpectraFilePaths.Clear();
                SpectralLibraryPaths.Clear();

                foreach (var connection in MsDataFiles)
                {
                    connection.Value.CloseDynamicConnection();
                }

                MsDataFiles.Clear();

                if (SpectralLibrary != null)
                {
                    SpectralLibrary.CloseConnections();
                }
            }
        }

        private void LoadPsms(out List<string> errors, bool haveLoadedSpectra)
        {
            errors = new List<string>();

            HashSet<string> fileNamesWithoutExtension = new HashSet<string>(
                SpectraFilePaths.Select(p => Path.GetFileName(p.Replace(GlobalVariables.GetFileExtension(p), string.Empty))));
            List<PsmFromTsv> psmsThatDontHaveMatchingSpectraFile = new List<PsmFromTsv>();

            try
            {
                foreach (var resultsFile in PsmResultFilePaths)
                {
                    lock (ThreadLocker)
                    {
                        foreach (PsmFromTsv psm in PsmTsvReader.ReadTsv(resultsFile, out List<string> warnings))
                        {
                            if (fileNamesWithoutExtension.Contains(psm.FileNameWithoutExtension) || !haveLoadedSpectra)
                            {
                                AllPsms.Add(psm);
                            }
                            else
                            {
                                psmsThatDontHaveMatchingSpectraFile.Add(psm);
                            }

                            if (PsmsGroupedByFile.TryGetValue(psm.FileNameWithoutExtension, out var psmsForThisFile))
                            {
                                psmsForThisFile.Add(psm);
                            }
                            else
                            {
                                PsmsGroupedByFile.Add(psm.FileNameWithoutExtension, new ObservableCollection<PsmFromTsv> { psm });
                            }
                        }
                    }
                }
            }
            catch (Exception e)
            {
                errors.Add("Error reading PSM file:\n" + e.Message);
            }

            if (psmsThatDontHaveMatchingSpectraFile.Any())
            {
                foreach (var file in psmsThatDontHaveMatchingSpectraFile.GroupBy(p => p.FileNameWithoutExtension))
                {
                    errors.Add(file.Count() + " PSMs from " + file.Key + " were not loaded because this spectra file was not found");
                }
            }

            FilterPsms();
        }

        private void LoadSpectraFiles(out List<string> errors)
        {
            errors = new List<string>();

            foreach (var filepath in SpectraFilePaths)
            {
                lock (ThreadLocker)
                {
                    var fileNameWithoutExtension = filepath.Replace(GlobalVariables.GetFileExtension(filepath), string.Empty);
                    fileNameWithoutExtension = Path.GetFileName(fileNameWithoutExtension);

                    DynamicDataConnection spectraFile = null;
                    string extension = GlobalVariables.GetFileExtension(filepath);

                    if (extension.Equals(".mzML", StringComparison.OrdinalIgnoreCase))
                    {
                        spectraFile = new MzmlDynamicData(filepath);
                    }
                    else if (extension.Equals(".mgf", StringComparison.OrdinalIgnoreCase))
                    {
                        spectraFile = new MgfDynamicData(filepath);
                    }
                    else if (extension.Equals(".raw", StringComparison.OrdinalIgnoreCase))
                    {
                        spectraFile = new ThermoDynamicData(filepath);
                    }
                    else
                    {
                        errors.Add("Unrecognized spectra file type: " + extension);
                        continue;
                    }

                    if (!MsDataFiles.TryAdd(fileNameWithoutExtension, spectraFile))
                    {
                        spectraFile.CloseDynamicConnection();
                        // print warning? but probably unnecessary. this means the data file was loaded twice. 
                        // which is an error but not an important one because the data is loaded
                    }
                }
            }
        }

        private void LoadSpectralLibraries(out List<string> errors)
        {
            errors = new List<string>();

            try
            {
                SpectralLibrary = new SpectralLibrary(SpectralLibraryPaths.ToList());
            }
            catch (Exception e)
            {
                SpectralLibrary = null;
                errors.Add("Problem loading spectral library: " + e.Message);
            }
        }
    }
}
