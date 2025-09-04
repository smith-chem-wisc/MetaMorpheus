using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Diagnostics;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Windows;
using System.Windows.Data;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using Easy.Common.Extensions;
using MzLibUtil;
using Omics;
using Readers;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace GuiFunctions.MetaDraw;

public class BioPolymerTabViewModel : MetaDrawTabViewModel
{
    private Dictionary<string, IBioPolymer> _allBioPolymers;
    private MetaDrawLogic _metaDrawLogic;

    public BioPolymerTabViewModel(MetaDrawLogic metaDrawLogic, string exportDirectory = null)
    {
        IsDatabaseLoaded = false;
        _metaDrawLogic = metaDrawLogic;
        _allBioPolymers = new Dictionary<string, IBioPolymer>();
        AllGroups = new ObservableCollection<BioPolymerGroupViewModel>();
        ExportDirectory = exportDirectory ?? Environment.GetFolderPath(Environment.SpecialFolder.MyDocuments);

        LoadDatabaseCommand = new RelayCommand(LoadDatabase);
        ResetDatabaseCommand = new RelayCommand(ResetDatabase);
        ExportImageCommand = new RelayCommand(ExportImage);

        BindingOperations.EnableCollectionSynchronization(AllGroups, ThreadLocker);
        BindingOperations.EnableCollectionSynchronization(FilteredGroups, ThreadLocker);
    }

    #region Database Loading Handling

    private bool _isDatabaseLoaded;

    public bool IsDatabaseLoaded
    {
        get => _isDatabaseLoaded;
        set
        {
            _isDatabaseLoaded = value;
            OnPropertyChanged(nameof(IsDatabaseLoaded));
        }
    }

    private string _databasePath;

    public string DatabasePath
    {
        get => _databasePath;
        set
        {
            _databasePath = value;
            OnPropertyChanged(nameof(DatabasePath));
            DatabaseName =PeriodTolerantFilenameWithoutExtension.GetPeriodTolerantFilenameWithoutExtension(_databasePath);
        }
    }

    private string _databaseName;

    public string DatabaseName
    {
        get => _databaseName;
        set
        {
            _databaseName = value;
            OnPropertyChanged(nameof(DatabaseName));
        }
    }

    public ICommand LoadDatabaseCommand { get; set; }
    public ICommand ResetDatabaseCommand { get; set; }

    private void LoadDatabase()
    {
        if (DatabasePath is null or "")
            return;
        try
        {
            _allBioPolymers = new SearchTask().LoadBioPolymers("", new()
                { new DbForTask(DatabasePath, false) }, true, DecoyType.None, new(), new())
                .ToDictionary(p => p.Accession, p => p);

            if (_allBioPolymers.Count == 0) 
                return;

            IsDatabaseLoaded = true;
            if (_metaDrawLogic.AllSpectralMatches is not null or { Count: 0 })
                ProcessSpectralMatches(_metaDrawLogic.AllSpectralMatches);
        }
        catch (Exception e)
        {
            // do nothing
        }
    }

    private void ResetDatabase()
    {
        DatabasePath = string.Empty;
        DatabaseName = string.Empty;
        _allBioPolymers.Clear();
        AllGroups.Clear();
        FilteredGroups.Clear();
    }

    #endregion

    private string _searchText;
    public string SearchText
    {
        get => _searchText;
        set
        {
            if (_searchText != value)
            {
                _searchText = value;
                OnPropertyChanged(nameof(SearchText));
                UpdateFilteredGroups();
            }
        }
    }

    public ObservableCollection<BioPolymerGroupViewModel> FilteredGroups { get; } = new();

    public ObservableCollection<BioPolymerGroupViewModel> AllGroups { get; set; }
    public BioPolymerCoverageMapViewModel CoverageMapViewModel { get; } = new();

    private BioPolymerGroupViewModel _selectedGroup;
    public BioPolymerGroupViewModel SelectedGroup
    {
        get => _selectedGroup;
        set
        {
            _selectedGroup = value;
            OnPropertyChanged(nameof(SelectedGroup));
            CoverageMapViewModel.Group = value;
        }
    }

    public void ProcessSpectralMatches(IList<SpectrumMatchFromTsv> matches)
    {
        AllGroups.Clear();

        // Build a mapping: accession => list of best matches for that accession
        var bestMatchesBySequence = matches
            .GroupBy(p => p.BaseSequence)
            .Select(group => group.MinBy(sm => sm.QValue))
            .ToList();

        // Map: accession => List<SpectrumMatchFromTsv>
        var accessionToMatches = new Dictionary<string, List<SpectrumMatchFromTsv>>();

        foreach (var match in bestMatchesBySequence)
        {
            foreach (var acc in match.Accession.Split('|'))
            {
                if (!accessionToMatches.TryGetValue(acc, out var list))
                {
                    list = new List<SpectrumMatchFromTsv>();
                    accessionToMatches[acc] = list;
                }
                list.Add(match);
            }
        }

        var processedResults = new List<BioPolymerCoverageResultModel>();
        foreach (var accessionToGroup in accessionToMatches.Keys) 
        {
            if (!_allBioPolymers.TryGetValue(accessionToGroup, out var bioPolymer))
                continue;

            var relevantMatches = accessionToMatches[accessionToGroup];
            processedResults.Clear();

            foreach (var match in relevantMatches)
            {
                BioPolymerCoverageType covType;

                var accessionSplits = match.Accession.Split("|");
                var siteSplits = match.GetStartAndEndPosition().ToArray();
                bool singleAccession = accessionSplits.Length == 1;
                bool singleLocalization = siteSplits.Length == 1;

                // TODO: Better ambiguity handling for missed cleavages. 
                var missedCleavageSplits = match.MissedCleavage.Split("|").Select(int.Parse).ToArray();
                int missedCleavages = missedCleavageSplits[0];

                // belongs to single biopolymer and single localization 
                if (singleAccession && singleLocalization)
                {
                    if (missedCleavages == 0)
                        covType = BioPolymerCoverageType.Unique;
                    else
                        covType = BioPolymerCoverageType.UniqueMissedCleavage;

                    processedResults.Add(new(match, match.BaseSeq, siteSplits[0].Start, siteSplits[0].End, covType));
                }
                // belongs to multiple biopolymers
                else if (accessionSplits.Length == siteSplits.Length) 
                {
                    int relevantIndex = accessionSplits.IndexOf(accessionToGroup);

                    if (missedCleavages == 0)
                        covType = BioPolymerCoverageType.Shared;
                    else
                        covType = BioPolymerCoverageType.SharedMissedCleavage;

                    processedResults.Add(new(match, match.BaseSeq, siteSplits[relevantIndex].Start, siteSplits[relevantIndex].End, covType));
                }
                // belongs to multiple biopolymers but shares same localization
                else if (!singleAccession && singleLocalization)
                {
                    if (missedCleavages == 0)
                        covType = BioPolymerCoverageType.Shared;
                    else
                        covType = BioPolymerCoverageType.SharedMissedCleavage;

                    processedResults.Add(new(match, match.BaseSeq, siteSplits[0].Start, siteSplits[0].End, covType));
                }
                // belongs to single biopolymer with multiple localizations
                else if (singleAccession) 
                {
                    if (missedCleavages == 0)
                        covType = BioPolymerCoverageType.TandemRepeat;
                    else
                        covType = BioPolymerCoverageType.TandemRepeatMissedCleavage;

                    for (int i = 0; i < siteSplits.Length; i++)
                    {
                        processedResults.Add(new(match, match.BaseSeq, siteSplits[i].Start, siteSplits[i].End, covType));
                    }
                }
                // Mismatch b/n accessions and start/end indexes. This occurs when one of the fields does not have a distinct value (e.g. peptide is in 3 protein sequences but has the same start and end in two of the sequences). 
                else
                {
                    var startIndicesInBioPolymer = bioPolymer.BaseSequence.IndexOfAll(match.BaseSeq)
                        .Select(p => p + 1)
                        .ToArray();
                    foreach (var potentialSite in siteSplits)
                    {
                        if (startIndicesInBioPolymer.Contains(potentialSite.Start))
                        {
                            // Found a match
                            if (missedCleavages == 0)
                                covType = BioPolymerCoverageType.Shared;
                            else
                                covType = BioPolymerCoverageType.SharedMissedCleavage;
                            processedResults.Add(new(match, match.BaseSeq, potentialSite.Start, potentialSite.End, covType));
                        }
                    }
                }
            }

            if (processedResults.Count == 0)
                continue;

            var group = new BioPolymerGroupViewModel(accessionToGroup, bioPolymer.Name, bioPolymer.BaseSequence, processedResults);
            AllGroups.Add(group);
        }

        UpdateFilteredGroups();
    }

    private void UpdateFilteredGroups()
    {
        FilteredGroups.Clear();
        var query = string.IsNullOrWhiteSpace(SearchText)
            ? AllGroups
            : AllGroups.Where(g =>
                (g.Accession?.IndexOf(SearchText, StringComparison.OrdinalIgnoreCase) >= 0) ||
                (g.ProteinName?.IndexOf(SearchText, StringComparison.OrdinalIgnoreCase) >= 0));
        foreach (var group in query.OrderByDescending(p => p.UniqueSequenceCoverage).ThenByDescending(p => p.MaximumSequenceCoverage))
            FilteredGroups.Add(group);
    }

    #region Image Export

    public ICommand ExportImageCommand { get; set; }

    private void ExportImage()
    {
        // Get the DrawingImage from the CoverageMapViewModel
        var drawingImage = CoverageMapViewModel.CoverageDrawing;
        if (drawingImage == null)
            return;

        // Use the DrawingImage's dimensions for export
        double width = 800;
        double height = 480;
        if (drawingImage.Drawing != null)
        {
            Rect bounds = drawingImage.Drawing.Bounds;
            if (!bounds.IsEmpty && bounds.Width > 0 && bounds.Height > 0)
            {
                width = bounds.Width;
                height = bounds.Height;
            }
        }

        var drawingVisual = new DrawingVisual();
        using (var dc = drawingVisual.RenderOpen())
        {
            dc.DrawImage(drawingImage, new Rect(0, 0, width, height));
        }

        // Use MetaDrawSettings.CanvasPdfExportDpi for DPI
        double dpi = MetaDrawSettings.CanvasPdfExportDpi;
        var rtb = new RenderTargetBitmap(
            (int)Math.Ceiling(width * dpi / 96.0),
            (int)Math.Ceiling(height * dpi / 96.0),
            dpi, dpi, PixelFormats.Pbgra32);
        rtb.Render(drawingVisual);

        string path = Path.Combine(ExportDirectory, $"{SelectedGroup.Accession}_SequenceCoverage.{MetaDrawSettings.ExportType}");

        // Convert RenderTargetBitmap to System.Drawing.Bitmap
        using (var ms = new MemoryStream())
        {
            var encoder = new PngBitmapEncoder();
            encoder.Frames.Add(BitmapFrame.Create(rtb));
            encoder.Save(ms);
            ms.Seek(0, SeekOrigin.Begin);
            using (var bmp = new System.Drawing.Bitmap(ms))
            {
                MetaDrawLogic.ExportBitmap(bmp, path);
            }
        }

        MessageBoxHelper.Show($"Exported Coverage Plot to {path}");
    }

    #endregion
}


[ExcludeFromCodeCoverage]
public class BioPolymerTabModel() : BioPolymerTabViewModel(new())
{
    public static BioPolymerTabModel Instance => new();
}