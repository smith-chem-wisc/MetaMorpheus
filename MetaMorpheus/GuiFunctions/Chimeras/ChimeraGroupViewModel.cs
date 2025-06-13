using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Linq;
using Chemistry;
using EngineLayer;
using GuiFunctions.ViewModels.Legends;
using MassSpectrometry;
using MzLibUtil;
using Omics.Fragmentation;
using OxyPlot;
using Readers;
using ClassExtensions = Chemistry.ClassExtensions;

namespace GuiFunctions;

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
    public int TotalFragments { get; }
    public int UniqueFragments { get; }
    public MsDataScan Ms1Scan { get; set; }
    public MsDataScan Ms2Scan { get; set; }
    public ObservableCollection<ChimericSpectralMatchModel> ChimericPsms { get; set; }

    #region Plotting 

    private List<string> _letters = new List<string> { "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z" };
    public Queue<string> Letters { get; }
    public Dictionary<string, List<ChimeraLegendItemViewModel>> LegendItems { get; set; }

    private bool IsColorInitialized { get; set; } = false;
    private Dictionary<OxyColor, List<(MatchedFragmentIon, string)>> _matchedFragmentIonsByColor = [];
    public Dictionary<OxyColor, List<(MatchedFragmentIon, string)>> MatchedFragmentIonsByColor
    {
        get
        {
            if (IsColorInitialized) return _matchedFragmentIonsByColor;
            AssignIonColors();
            return _matchedFragmentIonsByColor;
        }
        set
        {
            _matchedFragmentIonsByColor = value;
            OnPropertyChanged(nameof(MatchedFragmentIonsByColor));
        }
    }

    private Dictionary<OxyColor, List<(MatchedFragmentIon, string)>> _precursorIonsByColor = [];

    public Dictionary<OxyColor, List<(MatchedFragmentIon, string)>> PrecursorIonsByColor
    {
        get
        {
            if (IsColorInitialized) return _precursorIonsByColor;
            AssignIonColors();
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
                             ClassExtensions.ToMass(peak.mz, psm.PrecursorEnvelope.Charge),
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
                    }


                    _precursorIonsByColor.AddOrReplace(psm.Color, group.Key, annotation);
                }
                else
                    _precursorIonsByColor.AddOrReplace(psm.Color, group.Key, "Miso");
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

        IsColorInitialized = true;
    }

    #endregion

    public ChimeraGroupViewModel(IEnumerable<SpectrumMatchFromTsv> chimericSpectrumMatches, MsDataScan precursorSpectrum, MsDataScan fragmentationSpectrum)
    {
        Ms1Scan = precursorSpectrum;
        Ms2Scan = fragmentationSpectrum;
        LegendItems = new();
        Letters = new (_letters);

        ChimericPsms = [.. ConstructChimericPsmModels(chimericSpectrumMatches)];
        var representative = chimericSpectrumMatches.FirstOrDefault()!;
        FileNameWithoutExtension = representative.FileNameWithoutExtension;
        OneBasedPrecursorScanNumber = representative.PrecursorScanNum;
        Ms2ScanNumber = representative.Ms2ScanNumber;

        var terminalFrags = ChimericPsms.SelectMany(p => p.Psm.MatchedIons)
            .Where(p => !p.IsInternalFragment).ToList();
        TotalFragments = terminalFrags.Count;
        UniqueFragments = terminalFrags.Distinct().Count();
        ProteinCount = ChimericPsms.GroupBy(p => p.Psm.Accession).Count();
    }

    private IEnumerable<ChimericSpectralMatchModel> ConstructChimericPsmModels(IEnumerable<SpectrumMatchFromTsv> psms)
    {
        // Deconvolute the isolation window of the MS1 scan. 
        var deconParams = new DeconHostViewModel().PrecursorDeconvolutionParameters.Parameters; // Use default for current AnalyteType
        List<IsotopicEnvelope> envelopes = Ms2Scan.GetIsolatedMassesAndCharges(Ms1Scan, deconParams).ToList();

        // Prepare a list of available psms (not yet assigned)
        var availablePsms = psms.ToList();
        var matchedPsms = new List<(SpectrumMatchFromTsv, IsotopicEnvelope)>();

        // For each envelope, find the best matching PSM by minimizing a combined score of mz and neutral mass differences. 
        foreach (var envelope in envelopes.Where(p => p.Peaks.Count >= 2 && Ms2Scan.IsolationRange.MajorityWithin(p.Peaks.Select(m => m.mz)))
            .OrderByDescending(p => p.Peaks.Count)
            .ThenByDescending(p => p.Score))
        {
            SpectrumMatchFromTsv? bestPsm = null;
            double bestScore = double.MaxValue;

            foreach (var psm in availablePsms)
            {
                // Experimental mass is precursor mass + mass diff in daltons
                double experimentalMass = psm.MonoisotopicMass + double.Parse(psm.MassDiffDa.Split('|')[0]);
                double massDiff = Math.Abs(experimentalMass - envelope.MonoisotopicMass);

                // Calculate m/z for the envelope's monoisotopic peak
                double envelopeMz = envelope.MonoisotopicMass.ToMz(envelope.Charge);
                double mzDiff = Math.Abs(psm.PrecursorMz - envelopeMz);

                // Combine mass and m/z difference into a score (weights can be adjusted)
                double score = massDiff + mzDiff;

                if (score < bestScore)
                {
                    bestScore = score;
                    bestPsm = psm;
                }
            }

            if (bestPsm == null) continue;
            matchedPsms.Add((bestPsm, envelope));
            availablePsms.Remove(bestPsm); // Prevent double assignment
        }

        int proteinIndex = 0;
        foreach (var group in matchedPsms.GroupBy(p => p.Item1.Name)
                     .OrderByDescending(p => p.Count()))
        {
            var proteinColor = ChimeraSpectrumMatchPlot.ColorByProteinDictionary[proteinIndex][0];
            LegendItems.Add(group.Key, new List<ChimeraLegendItemViewModel>());

            if (group.Count() > 1)
                LegendItems[group.Key].Add(new ChimeraLegendItemViewModel("Shared Ions", proteinColor));

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