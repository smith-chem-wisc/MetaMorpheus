using System;
using System.Collections;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Linq;
using Chemistry;
using EngineLayer;
using MassSpectrometry;
using Omics.Fragmentation;
using OxyPlot;
using Readers;

namespace GuiFunctions;

/// <summary>
/// View model for a group of chimeric IDs from a single MS2 scan
/// </summary>
public class ChimeraGroupViewModel : BaseViewModel, IEnumerable<ChimericSpectralMatchModel>
{
    public string FileNameWithoutExtension { get; }
    public int OneBasedPrecursorScanNumber { get;}
    public int Ms2ScanNumber { get; }
    public int Count => ChimericPsms.Count;
    public int ProteinCount { get; }
    public int PrimarySequenceCount { get; }
    public int TotalFragments { get; }
    public int UniqueFragments { get; }
    public MsDataScan Ms1Scan { get; }
    public MsDataScan Ms2Scan { get; }
    public ObservableCollection<ChimericSpectralMatchModel> ChimericPsms { get; }

    #region Plotting 

    private List<string> _letters = new List<string> { "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z" };
    public Queue<string> Letters { get; }

    private bool IsColorInitialized { get; set; } = false;
    private Dictionary<OxyColor, List<(MatchedFragmentIon, string)>> _matchedFragmentIonsByColor = [];

    /// <summary>
    /// Key: Color, Value is a list of tuples where the first item is the matched fragment ion and the second item is an annotation string.
    /// </summary>
    public Dictionary<OxyColor, List<(MatchedFragmentIon, string)>> MatchedFragmentIonsByColor
    {
        get
        {
            if (IsColorInitialized) return _matchedFragmentIonsByColor;
            AssignIonColors();
            return _matchedFragmentIonsByColor;
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
    }

    /// <summary>
    /// Lazy loading of colors. Ensures we dont spend time parsing out what should be what color to never end up plotting it. 
    /// </summary>
    /// <param name="useLetterOnly"></param>
    internal void AssignIonColors(bool useLetterOnly = false)
    {
        // precursor peaks
        foreach (var group in ChimericPsms
                     .SelectMany(psm => psm.PrecursorEnvelope.Peaks
                         .Select(peak =>
                         {
                             var neutralTheoreticalProduct = new Product(ProductType.M, FragmentationTerminus.None,
                                 peak.mz.ToMass(psm.PrecursorEnvelope.Charge), 0, 0, 0);
                             return (psm, new MatchedFragmentIon(
                                 neutralTheoreticalProduct,
                                 peak.mz,
                                 peak.intensity,
                                 psm.PrecursorEnvelope.Charge));
                         }))
                     .GroupBy(p => p.Item2))
        {
            // distinct ions
            if (group.Count() != 1) 
                continue;

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

        _matchedFragmentIonsByColor = AssignFragmentIonColors();
        IsColorInitialized = true;
    }

    public Dictionary<OxyColor, List<(MatchedFragmentIon, string)>> AssignFragmentIonColors()
    {
        var toReturn = new Dictionary<OxyColor, List<(MatchedFragmentIon, string)>>();
        // matched fragment ions
        var accessionDict = ChimericPsms.Select(p => p.Psm.Accession)
            .Distinct()
            .ToDictionary(p => p, p => ChimericPsms.Count(psm => psm.Psm.Accession == p));
        foreach (var mzGroup in ChimericPsms
                     .SelectMany(psm => psm.Psm.MatchedIons
                         .Select(ion => (psm.Psm.Accession, psm.Color, psm.ProteinColor, ion)))
                     .GroupBy(g => g.ion.Mz))
        {
            var first = mzGroup.First();
            if (mzGroup.Count() == 1)
            {
                toReturn.AddOrReplace(first.Color, first.ion, "");
            }
            // if only one protein present
            else if (mzGroup.Select(p => p.Accession).Distinct().Count() == 1)
            {
                // if all proteoforms of the protein have the ion, protein shared color
                if (mzGroup.Count() == accessionDict[first.Accession])
                    toReturn.AddOrReplace(first.ProteinColor, first.ion, "");
                // if not all proteoforms have the same ion, their unique color
                else
                    foreach (var item in mzGroup)
                        toReturn.AddOrReplace(item.Color, item.ion, "");

            }
            // if only one mz value and multiple proteins, shared color
            else
            {
                toReturn.AddOrReplace(ChimeraSpectrumMatchPlot.MultipleProteinSharedColor, first.ion, "");
            }
        }
        return toReturn;
    }

    #endregion

    public ChimeraGroupViewModel(IEnumerable<SpectrumMatchFromTsv> chimericSpectrumMatches, MsDataScan precursorSpectrum, MsDataScan fragmentationSpectrum)
    {
        Ms1Scan = precursorSpectrum;
        Ms2Scan = fragmentationSpectrum;
        Letters = new (_letters);

        ChimericPsms = [.. ConstructChimericPsmModels(chimericSpectrumMatches)];
        var representative = ChimericPsms.FirstOrDefault()!.Psm;
        FileNameWithoutExtension = representative.FileNameWithoutExtension;
        OneBasedPrecursorScanNumber = representative.PrecursorScanNum;
        Ms2ScanNumber = representative.Ms2ScanNumber;

        var terminalFrags = ChimericPsms.SelectMany(p => p.Psm.MatchedIons)
            .Where(p => !p.IsInternalFragment).ToList();
        TotalFragments = terminalFrags.Count;
        UniqueFragments = terminalFrags.Select(p => (p.NeutralTheoreticalProduct.Annotation, p.Mz.RoundedDouble(1)))
            .Distinct()
            .Count();
        ProteinCount = ChimericPsms.GroupBy(p => p.Psm.Accession).Count();
        PrimarySequenceCount = ChimericPsms.Select(p => p.Psm.BaseSeq.Split('|')[0]).Distinct().Count();
    }

    private IEnumerable<ChimericSpectralMatchModel> ConstructChimericPsmModels(IEnumerable<SpectrumMatchFromTsv> psms)
    {
        // Deconvolute the isolation window of the MS1 scan. 
        var deconParams = MetaDrawSettingsViewModel.Instance.DeconHostViewModel.PrecursorDeconvolutionParameters.Parameters;
        List<IsotopicEnvelope> envelopes = Ms2Scan.GetIsolatedMassesAndCharges(Ms1Scan, deconParams)
            .Where(p => p.Peaks.Count >= 2) // Magic number for quality filtering. 
            .ToList();

        var psmList = psms.ToList();
        int n = psmList.Count;
        int m = envelopes.Count;
        double[,] costMatrix = new double[n, m];

        // Build cost matrix
        for (int i = 0; i < n; i++)
        {
            var psm = psmList[i];
            int psmCharge = psm.PrecursorCharge;
            double psmExperimentalMass = psm.PrecursorMass;

            for (int j = 0; j < m; j++)
            {
                double cost = 0;
                var env = envelopes[j];
                if (psmCharge != env.Charge)
                    cost += 1000; // Charge mismatch penalty, can be adjusted

                cost += Math.Abs(psmExperimentalMass - env.MonoisotopicMass); // Neutral Mass Diff Cost
                cost += Math.Abs(psmExperimentalMass.ToMz(psm.PrecursorCharge) - env.MonoisotopicMass.ToMz(env.Charge)); // Mz Diff Cost

                costMatrix[i, j] = cost;
            }
        }

        // Run Hungarian algorithm
        int[] assignment = HungarianAlgorithm.FindAssignments(costMatrix);

        // Build matched list
        var matchedPsms = new List<(SpectrumMatchFromTsv, IsotopicEnvelope)>();
        for (int i = 0; i < assignment.Length; i++)
        {
            int j = assignment[i];
            if (j >= 0 && j < m) 
            {
                matchedPsms.Add((psmList[i], envelopes[j]));
            }
        }

        int proteinIndex = 0;
        foreach (var group in matchedPsms.GroupBy(p => p.Item1.Name)
                     .OrderByDescending(p => p.Count()))
        {
            var proteinColor = ChimeraSpectrumMatchPlot.ColorByProteinDictionary[proteinIndex][0];

            if (group.Count() > 1)
            {
                for (int i = 0; i < group.Count(); i++)
                {
                    var color = ChimeraSpectrumMatchPlot.ColorByProteinDictionary[proteinIndex][i + 1];
                    var chimericPsm = new ChimericSpectralMatchModel(group.ElementAt(i).Item1, group.ElementAt(i).Item2,
                            color, proteinColor)
                        { Letter = Letters.Dequeue() };

                    yield return chimericPsm;
                }
            }
            else
            {
                var color = ChimeraSpectrumMatchPlot.ColorByProteinDictionary[proteinIndex][1];
                var chimericPsm = new ChimericSpectralMatchModel(group.First().Item1, group.First().Item2,
                        color, proteinColor)
                    { Letter = Letters.Dequeue() };
                yield return chimericPsm;
            }
            proteinIndex++;
        }
    }

    public IEnumerator<ChimericSpectralMatchModel> GetEnumerator() => ChimericPsms.GetEnumerator();
    IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();
}