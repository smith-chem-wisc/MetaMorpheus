using EngineLayer;
using OxyPlot;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using Readers;
using Omics;
using MzLibUtil;
using System.Collections.ObjectModel;
using System.Windows.Input;
using TaskLayer;
using UsefulProteomicsDatabases;
using System.Diagnostics;
using Easy.Common.Extensions;

namespace GuiFunctions;

public class BioPolymerTabViewModel : BaseViewModel
{
    private Dictionary<string, IBioPolymer> _allBioPolymers;
    private MetaDrawLogic _metaDrawLogic;

    public BioPolymerTabViewModel(MetaDrawLogic metaDrawLogic)
    {
        IsDatabaseLoaded = false;
        _metaDrawLogic = metaDrawLogic;
        _allBioPolymers = new Dictionary<string, IBioPolymer>();
        AllGroups = new ObservableCollection<BioPolymerGroupViewModel>();

        LoadDatabaseCommand = new RelayCommand(LoadDatabase);
        ResetDatabaseCommand = new RelayCommand(ResetDatabase);
        ExportImageCommand = new RelayCommand(ExportImage);
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
    }

    #endregion

    public ObservableCollection<BioPolymerGroupViewModel> AllGroups { get; set; }

    private BioPolymerGroupViewModel _selectedGroup;

    public BioPolymerGroupViewModel SelectedGroup
    {
        get => _selectedGroup;
        set
        {
            _selectedGroup = value;
            OnPropertyChanged(nameof(SelectedGroup));
        }
    }

    public void ProcessSpectralMatches(IList<SpectrumMatchFromTsv> matches)
    {
        var allAccessions = matches.SelectMany(p => p.Accession.Split('|'))
            .Distinct();
        var bestMatchesBySequence = matches.GroupBy(p => p.BaseSequence)
            .Select(group => group.MinBy(sm => sm.QValue))
            .ToList();

        foreach (var accessionToGroup in allAccessions) 
        {
            var relevantMatches = bestMatchesBySequence.Where(p => p.Accession.Contains(accessionToGroup));
            var processedResults = new List<BioPolymerCoverageResultModel>();

            foreach (var match in relevantMatches)
            {
                BioPolymerCoverageType covType;

                var accessionSplits = match.Accession.Split("|");
                var siteSplits = match.GetStartAndEndPosition().ToArray();

                // TODO: Better ambiguity handling for missed cleavages. 
                var missedCleavageSplits = match.MissedCleavage.Split("|").Select(int.Parse).ToArray();
                int missedCleavages = missedCleavageSplits[0];

                // belongs to single biopolymer and single localization 
                if (accessionSplits.Length == 1 && siteSplits.Length == 1)
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
                // belongs to single biopolymer with multiple localizations
                else if (accessionToGroup.Length == 1 && siteSplits.Length > 1) 
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
                else
                {
                    Debugger.Break();
                }
            }

            if (processedResults.Count == 0)
                continue;

            var bioPolymer = _allBioPolymers[accessionToGroup];
            var group = new BioPolymerGroupViewModel(accessionToGroup, bioPolymer.Name, bioPolymer.BaseSequence, processedResults);
            AllGroups.Add(group);
        }
    }

    #region Image Export

    public ICommand ExportImageCommand { get; set; }

    private void ExportImage()
    {
    }

    #endregion
}