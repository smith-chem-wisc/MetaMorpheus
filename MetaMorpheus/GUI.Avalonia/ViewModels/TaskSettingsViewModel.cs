using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using CommunityToolkit.Mvvm.ComponentModel;
using EngineLayer;
using EngineLayer.GlycoSearch;
using MassSpectrometry;
using MzLibUtil;
using Omics.Digestion;
using Proteomics.ProteolyticDigestion;
using SpectralAveraging;
using TaskLayer;
using Transcriptomics.Digestion;
using UsefulProteomicsDatabases;

namespace MetaMorpheusAvalonia.ViewModels;

/// <summary>
/// Editable settings for one task, shared by every task type.
///
/// Two different shapes have to be handled, and this is why the WPF task windows run to thousands of
/// lines each:
///
///   * CommonParameters exposes most of its properties with a private setter, so a change cannot be
///     bound onto an existing instance. It has to be rebuilt through its 42-argument constructor.
///     Apply() goes through CommonParameters.CloneWithNewValues, which copies every setting and
///     overrides only the ones shown here - the constructor would reset the rest to their defaults.
///     DigestionParams is the same shape, handled by BuildDigestionParams.
///   * Task-specific options - SearchParameters and friends - are plain settable properties, so those
///     are mutated in place.
///
/// Loading is the reverse: CommonParameters' getters are public, so the current values are read
/// straight out of the task.
/// </summary>
internal sealed partial class TaskSettingsViewModel : ObservableObject
{
    private readonly MetaMorpheusTask _task;

    public string TaskKind { get; }

    // --- shared across all task types -------------------------------------------------------
    [ObservableProperty] private string _precursorTolerance = "5";
    [ObservableProperty] private bool _precursorToleranceIsPpm = true;
    [ObservableProperty] private string _productTolerance = "20";
    [ObservableProperty] private bool _productToleranceIsPpm = true;
    [ObservableProperty] private DissociationType _dissociationType = DissociationType.HCD;
    [ObservableProperty] private string _protease = "trypsin";
    [ObservableProperty] private int _maxMissedCleavages = 2;
    [ObservableProperty] private int _minPeptideLength = 7;
    [ObservableProperty] private int _maxPeptideLength = int.MaxValue;
    [ObservableProperty] private int _maxModificationIsoforms = 1024;
    [ObservableProperty] private int _maxModsPerPeptide = 2;
    [ObservableProperty] private double _scoreCutoff = 5;
    [ObservableProperty] private double _qValueThreshold = 0.01;
    [ObservableProperty] private int _maxThreadsPerFile = -1;
    [ObservableProperty] private bool _doPrecursorDeconvolution = true;
    [ObservableProperty] private bool _useProvidedPrecursorInfo = true;
    [ObservableProperty] private bool _trimMsMsPeaks = true;
    [ObservableProperty] private bool _reportAllAmbiguity = true;

    // --- modifications, shared by all task types ---------------------------------------------
    public ModificationSelectionViewModel Modifications { get; private set; }

    // --- crosslink search only ---------------------------------------------------------------
    [ObservableProperty] private bool _isXlSearchTask;
    [ObservableProperty] private string _crosslinker = "DSSO";
    [ObservableProperty] private int _crosslinkSearchTopNum = 300;
    [ObservableProperty] private bool _crosslinkAtCleavageSite;
    [ObservableProperty] private bool _writePepXml = true;

    // --- glyco search only -------------------------------------------------------------------
    [ObservableProperty] private bool _isGlycoSearchTask;
    [ObservableProperty] private GlycoSearchType _glycoSearchType = GlycoSearchType.OGlycanSearch;
    [ObservableProperty] private int _glycoSearchTopNum = 50;
    [ObservableProperty] private int _maximumOGlycanAllowed = 4;
    [ObservableProperty] private bool _oxoniumIonFilt = true;

    // --- spectral averaging only -------------------------------------------------------------
    [ObservableProperty] private bool _isAveragingTask;
    [ObservableProperty] private int _numberOfScansToAverage = 5;
    [ObservableProperty] private int _scanOverlap = 4;

    // --- GPTMD only --------------------------------------------------------------------------
    [ObservableProperty] private bool _isGptmdTask;

    // --- search only -------------------------------------------------------------------------
    [ObservableProperty] private bool _isSearchTask;
    [ObservableProperty] private bool _doParsimony = true;
    [ObservableProperty] private bool _noOneHitWonders;
    [ObservableProperty] private bool _doQuantification = true;
    [ObservableProperty] private bool _matchBetweenRuns;
    [ObservableProperty] private bool _writePrunedDatabase;
    [ObservableProperty] private bool _doLabelFreeQuantification = true;
    [ObservableProperty] private DecoyType _decoyType = DecoyType.Reverse;
    [ObservableProperty] private MassDiffAcceptorType _massDiffAcceptorType = MassDiffAcceptorType.OneMM;

    /// <summary>
    /// Only meaningful when MassDiffAcceptorType is Custom, which is the whole point: offering Custom
    /// in the combo box without carrying this meant the choice could be made and then silently ignored.
    /// SearchTaskWindow.xaml.cs sets it; this did not.
    /// </summary>
    [ObservableProperty] private string _customMassDiffAcceptor = string.Empty;

    public IReadOnlyList<DissociationType> DissociationTypes { get; } =
        Enum.GetValues<DissociationType>().ToList();

    public IReadOnlyList<DecoyType> DecoyTypes { get; } = Enum.GetValues<DecoyType>().ToList();

    public IReadOnlyList<MassDiffAcceptorType> MassDiffAcceptorTypes { get; } =
        Enum.GetValues<MassDiffAcceptorType>().ToList();

    /// <summary>Protease names come from GlobalVariables, which loads them from the data directory.</summary>
    public IReadOnlyList<string> Proteases { get; } = ProteaseDictionary.Dictionary is null
        ? new List<string> { "trypsin" }
        : ProteaseDictionary.Dictionary.Keys.OrderBy(k => k).ToList();

    public TaskSettingsViewModel(MetaMorpheusTask task, string taskKind)
    {
        _task = task;
        TaskKind = taskKind;
        LoadFrom(task);
    }

    private void LoadFrom(MetaMorpheusTask task)
    {
        CommonParameters common = task.CommonParameters ?? new CommonParameters();

        (PrecursorTolerance, PrecursorToleranceIsPpm) = Describe(common.PrecursorMassTolerance);
        (ProductTolerance, ProductToleranceIsPpm) = Describe(common.ProductMassTolerance);
        DissociationType = common.DissociationType;
        ScoreCutoff = common.ScoreCutoff;
        QValueThreshold = common.QValueThreshold;
        MaxThreadsPerFile = common.MaxThreadsToUsePerFile;
        DoPrecursorDeconvolution = common.DoPrecursorDeconvolution;
        UseProvidedPrecursorInfo = common.UseProvidedPrecursorInfo;
        TrimMsMsPeaks = common.TrimMsMsPeaks;
        ReportAllAmbiguity = common.ReportAllAmbiguity;

        if (common.DigestionParams is { } digestion)
        {
            // SpecificProtease, not Protease: a non-specific search holds singleN/singleC in Protease
            // and the user's actual choice in SpecificProtease.
            Protease = (digestion as DigestionParams)?.SpecificProtease?.Name
                ?? digestion.DigestionAgent?.Name
                ?? Protease;
            MaxMissedCleavages = digestion.MaxMissedCleavages;
            MinPeptideLength = digestion.MinLength;
            MaxPeptideLength = digestion.MaxLength;
            MaxModificationIsoforms = digestion.MaxModificationIsoforms;
            MaxModsPerPeptide = digestion.MaxMods;
        }

        Modifications = new ModificationSelectionViewModel(
            common.ListOfModsFixed,
            common.ListOfModsVariable,
            (task as GptmdTask)?.GptmdParameters?.ListOfModsGptmd,
            isRna: common.DigestionParams is RnaDigestionParams);

        switch (task)
        {
            case GptmdTask:
                IsGptmdTask = true;
                break;

            case SearchTask search:
                IsSearchTask = true;
                DoParsimony = search.SearchParameters.DoParsimony;
                NoOneHitWonders = search.SearchParameters.NoOneHitWonders;
                DoQuantification = search.SearchParameters.DoLabelFreeQuantification;
                DoLabelFreeQuantification = search.SearchParameters.DoLabelFreeQuantification;
                MatchBetweenRuns = search.SearchParameters.MatchBetweenRuns;
                WritePrunedDatabase = search.SearchParameters.WritePrunedDatabase;
                DecoyType = search.SearchParameters.DecoyType;
                MassDiffAcceptorType = search.SearchParameters.MassDiffAcceptorType;
                CustomMassDiffAcceptor = search.SearchParameters.CustomMdac ?? string.Empty;
                break;

            case XLSearchTask xl:
                IsXlSearchTask = true;
                DecoyType = xl.XlSearchParameters.DecoyType;
                Crosslinker = xl.XlSearchParameters.Crosslinker?.CrosslinkerName ?? Crosslinker;
                CrosslinkSearchTopNum = xl.XlSearchParameters.CrosslinkSearchTopNum;
                CrosslinkAtCleavageSite = xl.XlSearchParameters.CrosslinkAtCleavageSite;
                WritePepXml = xl.XlSearchParameters.WritePepXml;
                break;

            case GlycoSearchTask glyco:
                IsGlycoSearchTask = true;
                DecoyType = glyco._glycoSearchParameters.DecoyType;
                GlycoSearchType = glyco._glycoSearchParameters.GlycoSearchType;
                GlycoSearchTopNum = glyco._glycoSearchParameters.GlycoSearchTopNum;
                MaximumOGlycanAllowed = glyco._glycoSearchParameters.MaximumOGlycanAllowed;
                OxoniumIonFilt = glyco._glycoSearchParameters.OxoniumIonFilt;
                DoParsimony = glyco._glycoSearchParameters.DoParsimony;
                NoOneHitWonders = glyco._glycoSearchParameters.NoOneHitWonders;
                break;

            case SpectralAveragingTask averaging:
                IsAveragingTask = true;
                // SpectralAveragingTask's parameterless constructor leaves Parameters null - its own
                // comment says it exists only for reading TOML - so a task built any other way has to
                // be given defaults before it can be edited or run.
                averaging.Parameters ??= new SpectralAveragingParameters();
                NumberOfScansToAverage = averaging.Parameters.NumberOfScansToAverage;
                ScanOverlap = averaging.Parameters.ScanOverlap;
                break;
        }
    }

    /// <summary>Crosslinker names come from GlobalVariables, which loads them from the data directory.</summary>
    public IReadOnlyList<string> Crosslinkers { get; } = GlobalVariables.Crosslinkers
        .Select(c => c.CrosslinkerName)
        .Where(n => !string.IsNullOrEmpty(n))
        .Distinct()
        .OrderBy(n => n, StringComparer.Ordinal)
        .ToList();

    public IReadOnlyList<GlycoSearchType> GlycoSearchTypes { get; } = Enum.GetValues<GlycoSearchType>().ToList();

    /// <summary>
    /// Which digestion parameters the default modifications should be derived from. CommonParameters
    /// picks a different default set for RNA, so this decides which one "reset" means.
    /// </summary>
    internal IDigestionParams DigestionParametersForDefaults => _task.CommonParameters?.DigestionParams;

    /// <summary>
    /// Every double here round-trips through InvariantCulture, matching the 157 uses in the WPF task
    /// windows. Current culture would be actively dangerous rather than merely inconsistent: on a
    /// de-DE machine "." is the group separator, so double.TryParse("0.05") succeeds with the value
    /// 5.0 and a user typing a tolerance silently gets one a hundred times looser.
    /// </summary>
    private static (string Value, bool IsPpm) Describe(Tolerance tolerance) => tolerance switch
    {
        PpmTolerance ppm => (Format(ppm.Value), true),
        AbsoluteTolerance absolute => (Format(absolute.Value), false),
        null => ("5", true),
        _ => (Format(tolerance.Value), true),
    };

    private static string Format(double value) => value.ToString(CultureInfo.InvariantCulture);

    private static bool TryParse(string value, out double parsed) =>
        double.TryParse(value, NumberStyles.Float, CultureInfo.InvariantCulture, out parsed);

    private static Tolerance BuildTolerance(string value, bool isPpm, double fallback)
    {
        if (!TryParse(value, out double parsed) || parsed <= 0)
        {
            parsed = fallback;
        }
        return isPpm ? new PpmTolerance(parsed) : new AbsoluteTolerance(parsed);
    }

    /// <summary>Validates the editable fields, returning the reasons rather than throwing.</summary>
    public IReadOnlyList<string> Validate()
    {
        var problems = new List<string>();
        if (!TryParse(PrecursorTolerance, out double precursor) || precursor <= 0)
        {
            problems.Add("Precursor mass tolerance must be a number greater than zero.");
        }
        if (!TryParse(ProductTolerance, out double product) || product <= 0)
        {
            problems.Add("Product mass tolerance must be a number greater than zero.");
        }
        if (MinPeptideLength < 1)
        {
            problems.Add("Minimum peptide length must be at least 1.");
        }
        if (MaxPeptideLength < MinPeptideLength)
        {
            problems.Add("Maximum peptide length cannot be below the minimum.");
        }
        if (MaxMissedCleavages < 0)
        {
            problems.Add("Missed cleavages cannot be negative.");
        }
        // >= 1, matching TaskValidator.CheckMinScoreAllowed. Same field, so the two GUIs must not
        // disagree about what is valid.
        if (ScoreCutoff < 1)
        {
            problems.Add("Score cutoff must be at least 1.");
        }

        // Custom is only a real choice if the expression parses, so reject it here rather than letting
        // the search fail later. SearchTaskWindow.xaml.cs:732 validates the same way.
        if (MassDiffAcceptorType == MassDiffAcceptorType.Custom)
        {
            try
            {
                SearchTask.GetMassDiffAcceptor(null, MassDiffAcceptorType.Custom, CustomMassDiffAcceptor);
            }
            catch (Exception exception)
            {
                problems.Add($"The custom mass difference acceptor could not be read: {exception.Message}");
            }
        }
        return problems;
    }

    /// <summary>
    /// Rebuilds the digestion parameters, carrying forward the six settings this dialog does not show.
    /// Most of DigestionParams' properties have private setters, so the object has to go back through
    /// its constructor; every argument not edited here is read off <paramref name="existing"/>.
    /// </summary>
    private IDigestionParams BuildDigestionParams(IDigestionParams existing) => existing switch
    {
        // The enzyme box is filled from ProteaseDictionary, so it has nothing to offer an RNA task and
        // the rnase is carried through untouched. Without this branch an RNA task came back proteolytic.
        RnaDigestionParams rna => new RnaDigestionParams(
            rnase: rna.Rnase.Name,
            maxMissedCleavages: MaxMissedCleavages,
            minLength: MinPeptideLength,
            maxLength: MaxPeptideLength,
            maxModificationIsoforms: MaxModificationIsoforms,
            maxMods: MaxModsPerPeptide,
            fragmentationTerminus: rna.FragmentationTerminus)
        {
            SearchModeType = rna.SearchModeType,
        },

        // searchModeType and fragmentationTerminus have to be passed for a semi or non-specific search
        // to survive; RecordSpecificProtease() then re-derives Protease from them, so the name given
        // here is the specific protease either way.
        DigestionParams proteolytic => new DigestionParams(
            protease: Protease,
            maxMissedCleavages: MaxMissedCleavages,
            minPeptideLength: MinPeptideLength,
            maxPeptideLength: MaxPeptideLength,
            maxModificationIsoforms: MaxModificationIsoforms,
            initiatorMethionineBehavior: proteolytic.InitiatorMethionineBehavior,
            maxModsForPeptides: MaxModsPerPeptide,
            searchModeType: proteolytic.SearchModeType,
            fragmentationTerminus: proteolytic.FragmentationTerminus,
            generateUnlabeledProteinsForSilac: proteolytic.GeneratehUnlabeledProteinsForSilac,
            keepNGlycopeptide: proteolytic.KeepNGlycopeptide,
            keepOGlycopeptide: proteolytic.KeepOGlycopeptide),

        _ => new DigestionParams(
            protease: Protease,
            maxMissedCleavages: MaxMissedCleavages,
            minPeptideLength: MinPeptideLength,
            maxPeptideLength: MaxPeptideLength,
            maxModificationIsoforms: MaxModificationIsoforms,
            maxModsForPeptides: MaxModsPerPeptide),
    };

    /// <summary>
    /// Writes the edited values back onto the task. CommonParameters and the digestion parameters are
    /// cloned-and-overridden because their setters are private; the search options are assigned in
    /// place because theirs are not.
    /// </summary>
    public void Apply()
    {
        CommonParameters existing = _task.CommonParameters ?? new CommonParameters();
        IDigestionParams digestion = BuildDigestionParams(existing.DigestionParams);

        // CloneWithNewValues, not the constructor: this view model shows 13 of the settings
        // CommonParameters holds, and the constructor would default the rest back.
        _task.CommonParameters = existing.CloneWithNewValues(
            dissociationType: DissociationType,
            doPrecursorDeconvolution: DoPrecursorDeconvolution,
            useProvidedPrecursorInfo: UseProvidedPrecursorInfo,
            reportAllAmbiguity: ReportAllAmbiguity,
            trimMsMsPeaks: TrimMsMsPeaks,
            qValueThreshold: QValueThreshold,
            scoreCutoff: ScoreCutoff,
            precursorMassTolerance: BuildTolerance(PrecursorTolerance, PrecursorToleranceIsPpm, 5),
            productMassTolerance: BuildTolerance(ProductTolerance, ProductToleranceIsPpm, 20),
            maxThreadsToUsePerFile: MaxThreadsPerFile,
            digestionParams: digestion,
            listOfModsVariable: Modifications.VariableSelection,
            listOfModsFixed: Modifications.FixedSelection);

        switch (_task)
        {
            case GptmdTask gptmd:
                gptmd.GptmdParameters.ListOfModsGptmd = Modifications.GptmdSelection.ToList();
                break;

            case SearchTask search:
                search.SearchParameters.DoParsimony = DoParsimony;
                search.SearchParameters.NoOneHitWonders = NoOneHitWonders;
                search.SearchParameters.DoLabelFreeQuantification = DoLabelFreeQuantification;
                search.SearchParameters.MatchBetweenRuns = MatchBetweenRuns;
                search.SearchParameters.WritePrunedDatabase = WritePrunedDatabase;
                search.SearchParameters.DecoyType = DecoyType;
                search.SearchParameters.MassDiffAcceptorType = MassDiffAcceptorType;
                search.SearchParameters.CustomMdac = CustomMassDiffAcceptor;
                break;

            case XLSearchTask xl:
                xl.XlSearchParameters.DecoyType = DecoyType;
                xl.XlSearchParameters.CrosslinkSearchTopNum = CrosslinkSearchTopNum;
                xl.XlSearchParameters.CrosslinkAtCleavageSite = CrosslinkAtCleavageSite;
                xl.XlSearchParameters.WritePepXml = WritePepXml;
                // Look the crosslinker up by name rather than keeping the object, so the choice
                // survives being written to and read back from TOML.
                Crosslinker chosen = GlobalVariables.Crosslinkers
                    .FirstOrDefault(c => c.CrosslinkerName == Crosslinker);
                if (chosen is not null)
                {
                    xl.XlSearchParameters.Crosslinker = chosen;
                }
                break;

            case GlycoSearchTask glyco:
                glyco._glycoSearchParameters.DecoyType = DecoyType;
                glyco._glycoSearchParameters.GlycoSearchType = GlycoSearchType;
                glyco._glycoSearchParameters.GlycoSearchTopNum = GlycoSearchTopNum;
                glyco._glycoSearchParameters.MaximumOGlycanAllowed = MaximumOGlycanAllowed;
                glyco._glycoSearchParameters.OxoniumIonFilt = OxoniumIonFilt;
                glyco._glycoSearchParameters.DoParsimony = DoParsimony;
                glyco._glycoSearchParameters.NoOneHitWonders = NoOneHitWonders;
                break;

            case SpectralAveragingTask averaging:
                averaging.Parameters ??= new SpectralAveragingParameters();
                averaging.Parameters.NumberOfScansToAverage = NumberOfScansToAverage;
                averaging.Parameters.ScanOverlap = ScanOverlap;
                break;
        }
    }
}
