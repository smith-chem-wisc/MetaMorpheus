using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Linq;
using System.Windows;
using Chemistry;
using Easy.Common.Extensions;
using EngineLayer;
using GuiFunctions.ViewModels.Legends;
using MassSpectrometry;
using MzLibUtil;
using Omics;
using Omics.Fragmentation;
using Omics.Modifications;
using OxyPlot;
using Readers;

namespace GuiFunctions;

/// <summary>
/// All data and user triggered manipulations for the Chimera Analysis tab
/// </summary>
public class ChimeraAnalysisTabViewModel : BaseViewModel
{
    #region Displayed in GUI
    public List<ChimeraGroupViewModel> ChimeraGroupViewModels { get; set; }

    private ChimeraLegendViewModel _chimeraLegendViewModel;
    public ChimeraLegendViewModel ChimeraLegendViewModel
    {
        get => _chimeraLegendViewModel;
        set
        {
            _chimeraLegendViewModel = value;
            OnPropertyChanged(nameof(ChimeraLegendViewModel));
        }
    }

    private ChimeraDrawnSequence _chimeraDrawnSequence;
    public ChimeraDrawnSequence ChimeraDrawnSequence
    {
        get => _chimeraDrawnSequence;
        set
        {
            _chimeraDrawnSequence = value;
            OnPropertyChanged(nameof(ChimeraDrawnSequence));
        }
    }

    #endregion

    #region Settings that change behavior

    private bool _groupProteinsInSequenceAnnotation;
    public bool GroupProteinsInSequenceAnnotation
    {
        get => _groupProteinsInSequenceAnnotation;
        set
        {
            _groupProteinsInSequenceAnnotation = value;
            OnPropertyChanged(nameof(GroupProteinsInSequenceAnnotation));
        }
    }

    #endregion

    public ChimeraAnalysisTabViewModel(List<SpectrumMatchFromTsv> allPsms, Dictionary<string, MsDataFile> dataFiles)
    {
        ChimeraGroupViewModels = ConstructChimericPsms(allPsms, dataFiles)
            .OrderByDescending(p => p.Count)
            .ToList();
    }

    private static IEnumerable<ChimeraGroupViewModel> ConstructChimericPsms(List<SpectrumMatchFromTsv> psms, Dictionary<string, MsDataFile> dataFiles)
    {
        return psms
            .Where(p => p.QValue <= 0.01 && p.DecoyContamTarget == "T")
            .GroupBy(p => (p.FileNameWithoutExtension, p.Ms2ScanNumber))
            .Select(group =>
            {
                if (!dataFiles.TryGetValue(group.First().FileNameWithoutExtension, out MsDataFile spectraFile))
                    return null;

                var ms1Scan = spectraFile.GetOneBasedScanFromDynamicConnection(group.First().PrecursorScanNum);
                var ms2Scan = spectraFile.GetOneBasedScanFromDynamicConnection(group.First().Ms2ScanNumber);

                if (ms1Scan == null || ms2Scan == null)
                    return null;

                var groupVm = new ChimeraGroupViewModel(group.OrderBy(p => p.PrecursorMz), ms1Scan, ms2Scan);
                return groupVm.ChimericPsms.Count > 0 ? groupVm : null;
            })
            .Where(groupVm => groupVm != null);
    }
}

/// <summary>
/// View model for a group of chimeric IDs from a single MS2 scan
/// </summary>
public class ChimeraGroupViewModel : BaseViewModel
{
    public string FileNameWithoutExtension { get; set; }
    public int OneBasedPrecursorScanNumber { get; set; }
    public int Ms2ScanNumber { get; set; }
    public int Count => ChimericPsms.Count;
    public int ProteinCount { get; }
    public MsDataScan Ms1Scan { get; set; }
    public MsDataScan Ms2Scan { get; set; }
    public ObservableCollection<ChimericSpectralMatchModel> ChimericPsms { get; set; }

    #region Plotting 

    private List<string> _letters = new List<string> { "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z" };
    public Queue<string> Letters { get; }
    public Dictionary<string, List<ChimeraLegendItemViewModel>> LegendItems { get; set; }

    private bool IsColorInitialized { get; set; } = false;
    private Dictionary<OxyColor, List<(MatchedFragmentIon, string)>> _matchedFragmentIonsByColor;
    public Dictionary<OxyColor, List<(MatchedFragmentIon, string)>> MatchedFragmentIonsByColor
    {
        get
        {
            if (IsColorInitialized) return _matchedFragmentIonsByColor;
            AssignIonColors();
            IsColorInitialized = true;
            return _matchedFragmentIonsByColor;
        }
        set
        {
            _matchedFragmentIonsByColor = value;
            OnPropertyChanged(nameof(MatchedFragmentIonsByColor));
        }
    }

    private Dictionary<OxyColor, List<(MatchedFragmentIon, string)>> _precursorIonsByColor;

    public Dictionary<OxyColor, List<(MatchedFragmentIon, string)>> PrecursorIonsByColor
    {
        get
        {
            if (IsColorInitialized) return _precursorIonsByColor;
            AssignIonColors();
            IsColorInitialized = true;
            return _precursorIonsByColor;
        }
        set
        {
            _precursorIonsByColor = value;
            OnPropertyChanged(nameof(PrecursorIonsByColor));
        }
    }

    /// <summary>
    /// Lazy loading of colors. Ensures we dont spend time parsing out what should be what color to never end up plotting it. 
    /// </summary>
    /// <param name="useLetterOnly"></param>
    internal void AssignIonColors(bool useLetterOnly = false)
    {
        // precursor peaks
        foreach (var group in ChimericPsms.SelectMany(psm => psm.PrecursorEnvelope.Peaks.Select(peak =>
        {
            var neutralTheoreticalProduct = new Product(ProductType.M, FragmentationTerminus.None,
                peak.mz.ToMass(psm.PrecursorEnvelope.Charge),
                0, 0, 0);
            return (psm, new MatchedFragmentIon(
                neutralTheoreticalProduct,
                peak.mz,
                peak.intensity,
                psm.PrecursorEnvelope.Charge));
        }))
                     .GroupBy(p => p.Item2))
        {

            // distinct ions
            if (group.Count() == 1)
            {
                var psm = group.First().psm;
                var maxIntensityPrecursorIon = psm.PrecursorEnvelope.Peaks.MaxBy(p => p.intensity);
                if (Math.Abs(group.Key.Intensity - maxIntensityPrecursorIon.intensity) < 0.00001)
                {
                    string annotation = "";

                    if (useLetterOnly)
                    {
                        annotation += psm.Letter;
                    }
                    else
                    {
                        annotation += $"Charge = {group.Key.Charge}";
                        annotation += $"\nm/z = {group.Key.Mz:0.00}";
                        annotation += $"\nMono Mass = {psm.PrecursorEnvelope.MonoisotopicMass:0.00}";
                        annotation += $"\nProtein = {psm.Psm.Name}";


                        //PeptideWithSetModifications pepWithSetMods = new(psm.Psm.FullSequence.Split("|")[0], GlobalVariables.AllModsKnownDictionary);
                        //foreach (var mod in pepWithSetMods.AllModsOneIsNterminus)
                        //{
                        //    annotation += $"\n{mod.Value.IdWithMotif}{mod.Key}";
                        //}
                    }


                    _precursorIonsByColor.AddOrReplace(psm.Color, group.Key, annotation);
                }
                else
                    _precursorIonsByColor.AddOrReplace(psm.Color, group.Key, "");
            }
            // shared ions
            else
            {
                if (group.Select(p => p.psm.Psm.Accession).Distinct().Count() == 1)
                {
                    _precursorIonsByColor.AddOrReplace(group.First().psm.ProteinColor, group.Key, "");
                }
                else
                {
                    _precursorIonsByColor.AddOrReplace(ChimeraSpectrumMatchPlot.MultipleProteinSharedColor, group.Key, "");
                }
            }
        }

        // matched fragment ions
        var accessionDict = ChimericPsms.Select(p => p.Psm.Accession)
            .Distinct()
            .ToDictionary(p => p, p => ChimericPsms.Count(psm => psm.Psm.Accession == p));
        foreach (var annotationGroup in ChimericPsms
                     .SelectMany(psm => psm.Psm.MatchedIons
                         .Select(ion => (psm.Psm.Accession, psm.Color, psm.ProteinColor, ion)))
                         .GroupBy(g => g.ion.Annotation))
        {
            //TODO: Group by mz
            // if only one ion has the same annotation, unique proteoform color
            if (annotationGroup.Count() == 1)
                _matchedFragmentIonsByColor.AddOrReplace(annotationGroup.First().Color, annotationGroup.First().ion, "");
            else
            {
                foreach (var mzGroup in annotationGroup.GroupBy(p => p.ion.Mz))
                {
                    if (mzGroup.Count() == 1)
                    {
                        _matchedFragmentIonsByColor.AddOrReplace(mzGroup.First().Color, mzGroup.First().ion, "");
                    }
                    // if only one protein present
                    else if (mzGroup.Select(p => p.Accession).Distinct().Count() == 1)
                    {
                        // if all proteoforms of the protein have the ion, protein shared color
                        if (mzGroup.Count() == accessionDict[mzGroup.First().Accession])
                            _matchedFragmentIonsByColor.AddOrReplace(mzGroup.First().ProteinColor, mzGroup.First().ion, "");
                        // if not all proteoforms have the same ion, their unique color
                        else
                            foreach (var item in mzGroup)
                                _matchedFragmentIonsByColor.AddOrReplace(item.Color, item.ion, "");

                    }
                    // if only one mz value and multiple proteins, shared color
                    else
                    {
                        _matchedFragmentIonsByColor.AddOrReplace(ChimeraSpectrumMatchPlot.MultipleProteinSharedColor, mzGroup.First().ion, "");
                    }
                }
            }
        }
    }

    #endregion

    public ChimeraGroupViewModel(IEnumerable<SpectrumMatchFromTsv> chimericSpectrumMatches, MsDataScan precursorSpectrum, MsDataScan fragmentationSpectrum)
    {
        Ms1Scan = precursorSpectrum;
        Ms2Scan = fragmentationSpectrum;
        ChimericPsms = [..ConstructChimericPsmModels(chimericSpectrumMatches)];
        LegendItems = new();
        Letters = new Queue<string>(_letters);

        var representative = ChimericPsms.FirstOrDefault()!.Psm;
        FileNameWithoutExtension = representative.FileNameWithoutExtension;
        OneBasedPrecursorScanNumber = representative.PrecursorScanNum;
        Ms2ScanNumber = representative.Ms2ScanNumber;
    }

    private IEnumerable<ChimericSpectralMatchModel> ConstructChimericPsmModels(IEnumerable<SpectrumMatchFromTsv> psms)
    {

        // Deconvolute the isolation window of the MS1 scan. 
        var deconParams = new DeconHostViewModel().PrecursorDeconvolutionParameters.Parameters; // Use default for current AnalyteType
        List<IsotopicEnvelope> envelopes = Ms2Scan.GetIsolatedMassesAndCharges(Ms1Scan, deconParams).ToList();

        var tolerance = new PpmTolerance(100);

        // TODO: polish this. 
        // match each envelope with a SpectrumMatch based upon the spectrumMatches peptidemonoMass + massdiffda
        // considering each scan needs to be matched with teh closes spectrum match
        List<(SpectrumMatchFromTsv, IsotopicEnvelope)> matchedPsms = new();
        foreach (var scan in envelopes.Where(p => p.Peaks.Count >= 3 && Ms2Scan.IsolationRange.MajorityWithin(p.Peaks.Select(m => m.mz))))
        {
            var psm = psms
                .Where(p => tolerance.Within(p.PrecursorMass, scan.MonoisotopicMass)
                //|| (MetaDrawSettings.CheckMzForChimeras && mzTolerance.Within(p.PrecursorMz, scan.PrecursorMonoisotopicPeakMz))
                //|| scan.PrecursorEnvelope.Peaks.Any(peak => tolerance.Within(peak.mz, p.PrecursorMz)))
                )
                .OrderBy(p => Math.Abs(p.PrecursorMass - scan.MonoisotopicMass))
                .FirstOrDefault();

            if (psm != null)
                matchedPsms.Add((psm, scan));
        }

        var distinct = matchedPsms
            .GroupBy(p => p.Item1.QValue + p.Item1.FullSequence)
                .Select(g => g.First())
                .OrderBy(p => p.Item1.PrecursorMz)
                .ToList();


        int proteinIndex = 0;
        foreach (var group in distinct.GroupBy(p => p.Item1.Name)
                        .OrderByDescending(p => p.Count()))
        {
            var proteinColor = ChimeraSpectrumMatchPlot.ColorByProteinDictionary[proteinIndex][0];
            LegendItems.Add(group.Key, new List<ChimeraLegendItemViewModel>());

            //if (annotationGroup.Count( ) > 1)
            //LegendItems[annotationGroup.Key].Add(new ChimeraLegendItemViewModel("Shared Ions", proteinColor));
            for (int i = 0; i < group.Count(); i++)
            {
                var color = ChimeraSpectrumMatchPlot.ColorByProteinDictionary[proteinIndex][i + 1];
                var chimericPsm = new ChimericSpectralMatchModel(group.ElementAt(i).Item1, group.ElementAt(i).Item2,
                        color, proteinColor)
                    { Letter = Letters.Dequeue() };
                LegendItems[group.Key].Add(new(chimericPsm.ModString, color));
                yield return chimericPsm;
            }
            proteinIndex++;
        }
    }
}



/// <summary>
/// View model for a single chimeric spectral match
/// </summary>
public class ChimericSpectralMatchModel : BaseViewModel
{
    public SpectrumMatchFromTsv Psm { get; set; }
    public IsotopicEnvelope PrecursorEnvelope { get; set; }
    public OxyColor Color { get; set; }
    public OxyColor ProteinColor { get; set; }
    public string Letter { get; set; }
    public string ModString { get; set; }
    public Dictionary<int, Modification> AllModsOneIsNterminus { get; set; }

    public ChimericSpectralMatchModel(SpectrumMatchFromTsv psm, IsotopicEnvelope precursorEnvelope, OxyColor color, OxyColor proteinColor)
    {
        Psm = psm;
        PrecursorEnvelope = precursorEnvelope;
        Color = color;
        ProteinColor = proteinColor;

        var modDict = GlobalVariables.AnalyteType switch {
            AnalyteType.Oligo => GlobalVariables.AllRnaModsKnownDictionary,
            _ => GlobalVariables.AllModsKnownDictionary
        };
        AllModsOneIsNterminus = IBioPolymerWithSetMods.GetModificationDictionaryFromFullSequence(psm.FullSequence, modDict);
        ModString = string.Join(", ", AllModsOneIsNterminus.Select(m => $"{m.Key} - {m.Value.IdWithMotif}"));
    }
}


internal static class Extensions
{
    /// <summary>
    /// Determines if a majority of values are within a range
    /// </summary>
    /// <param name="range"></param>
    /// <param name="values"></param>
    /// <returns></returns>
    internal static bool MajorityWithin(this MzRange range, IEnumerable<double> values)
    {
        int within = values.Count(p => p >= range.Minimum && p <= range.Maximum);
        return within > values.Count() / 2;
    }

    // method to add a value to a list in a dictionary if the key is present, and craete a new list if the key is not present
    public static void AddOrReplace<TKey, TValue, TValue2>(this Dictionary<TKey, List<(TValue, TValue2)>> dictionary, TKey key,
        TValue value, TValue2 value2)
    {
        if (dictionary.ContainsKey(key))
        {
            var previousVersion = dictionary[key].FirstOrDefault(p => p.Item1.Equals(value));
            if (!previousVersion.GetType().IsDefault())
            {
                dictionary[key].Remove(previousVersion);
            }
            dictionary[key].Add((value, value2));

        }
        else
            dictionary.Add(key, new List<(TValue, TValue2)> { (value, value2) });

    }
}