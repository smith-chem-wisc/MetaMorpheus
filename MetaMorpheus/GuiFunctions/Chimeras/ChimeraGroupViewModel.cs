using System;
using System.Collections;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Linq;
using Chemistry;
using EngineLayer;
using GuiFunctions.ViewModels.Legends;
using MassSpectrometry;
using MzLibUtil;
using Nett;
using Omics.Fragmentation;
using OxyPlot;
using OxyPlot.Reporting;
using Readers;
using ClassExtensions = Chemistry.ClassExtensions;

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
    public int TotalFragments { get; }
    public int UniqueFragments { get; }
    public MsDataScan Ms1Scan { get; }
    public MsDataScan Ms2Scan { get; }
    public ObservableCollection<ChimericSpectralMatchModel> ChimericPsms { get; }

    #region Plotting 

    private List<string> _letters = new List<string> { "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z" };
    public Queue<string> Letters { get; }
    public Dictionary<string, List<ChimeraLegendItemViewModel>> LegendItems { get; }

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
        UniqueFragments = terminalFrags.Select(p => (p.NeutralTheoreticalProduct.Annotation, p.Mz.RoundedDouble(1)))
            .Distinct()
            .Count();
        ProteinCount = ChimericPsms.GroupBy(p => p.Psm.Accession).Count();
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

    public IEnumerator<ChimericSpectralMatchModel> GetEnumerator() => ChimericPsms.GetEnumerator();
    IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();
}

/// <summary>
/// The Hungarian algorithm, also known as the Kuhn-Munkres algorithm or Munkres assignment algorithm, is a powerful combinatorial optimization algorithm used to solve the assignment problem. The goal is to find the optimal assignment between two sets of equal size(e.g., psms and precursor envelopes) to minimize the total cost or maximize the total value, given a cost matrix representing the cost of assigning each psms to each envelopes where cost is mass error.
/// </summary>
public static class HungarianAlgorithm
{
    /// <summary>
    /// Solves the assignment problem (also known as the minimum cost bipartite matching problem) using the Hungarian algorithm.
    /// Given a cost matrix where each entry represents the cost of assigning a row (e.g., a PSM) to a column (e.g., an isotopic envelope),
    /// this method finds the optimal one-to-one assignment that minimizes the total cost.
    /// If the input matrix is not square, it is padded with zeros to make it square.
    /// Returns an array where each element at index i indicates the column assigned to row i, or -1 if no assignment is possible.
    /// </summary>
    /// <param name="costMatrix">
    /// A 2D array of assignment costs. Rows typically represent items to assign (e.g., PSMs), columns represent targets (e.g., envelopes).
    /// Use a large value (e.g., double.MaxValue) to indicate forbidden assignments.
    /// </param>
    /// <returns>
    /// An array of assignments: result[i] = j means row i is assigned to column j; -1 means row i is unassigned.
    /// </returns>
    public static int[] FindAssignments(double[,] costMatrix)
    {
        int rows = costMatrix.GetLength(0);
        int cols = costMatrix.GetLength(1);
        int dim = Math.Max(rows, cols);

        // Pad to square matrix if needed
        double[,] cost = new double[dim, dim];
        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++)
                cost[i, j] = (i < rows && j < cols) ? costMatrix[i, j] : 0;

        // Step 1: Subtract row minima
        for (int i = 0; i < dim; i++)
        {
            double min = double.MaxValue;
            for (int j = 0; j < dim; j++)
                if (cost[i, j] < min) min = cost[i, j];
            for (int j = 0; j < dim; j++)
                cost[i, j] -= min;
        }

        // Step 2: Subtract column minima
        for (int j = 0; j < dim; j++)
        {
            double min = double.MaxValue;
            for (int i = 0; i < dim; i++)
                if (cost[i, j] < min) min = cost[i, j];
            for (int i = 0; i < dim; i++)
                cost[i, j] -= min;
        }

        // Marked zeros
        int[] rowCover = new int[dim];
        int[] colCover = new int[dim];
        int[,] mask = new int[dim, dim];

        // Step 3: Cover zeros with minimum number of lines
        int step = 1;
        int[] pathRow = new int[dim * 2];
        int[] pathCol = new int[dim * 2];
        int pathCount = 0;

        while (true)
        {
            switch (step)
            {
                case 1:
                    // Star zeros
                    for (int i = 0; i < dim; i++)
                        for (int j = 0; j < dim; j++)
                            if (cost[i, j] == 0 && rowCover[i] == 0 && colCover[j] == 0)
                            {
                                mask[i, j] = 1;
                                rowCover[i] = 1;
                                colCover[j] = 1;
                            }
                    Array.Clear(rowCover, 0, dim);
                    Array.Clear(colCover, 0, dim);
                    step = 2;
                    break;

                case 2:
                    // Cover columns with starred zeros
                    for (int i = 0; i < dim; i++)
                        for (int j = 0; j < dim; j++)
                            if (mask[i, j] == 1)
                                colCover[j] = 1;
                    if (colCover.Sum() >= dim)
                        step = 7;
                    else
                        step = 3;
                    break;

                case 3:
                    // Prime uncovered zeros
                    bool done = false;
                    int zRow = -1, zCol = -1;
                    while (!done)
                    {
                        (zRow, zCol) = FindZero(cost, rowCover, colCover, dim);
                        if (zRow == -1)
                        {
                            step = 5;
                            done = true;
                        }
                        else
                        {
                            mask[zRow, zCol] = 2;
                            // Find the column in this row where mask[zRow, col] == 1
                            int starCol = -1;
                            for (int j = 0; j < dim; j++)
                            {
                                if (mask[zRow, j] == 1)
                                {
                                    starCol = j;
                                    break;
                                }
                            }
                            if (starCol != -1)
                            {
                                rowCover[zRow] = 1;
                                colCover[starCol] = 0;
                            }
                            else
                            {
                                step = 4;
                                pathRow[0] = zRow;
                                pathCol[0] = zCol;
                                pathCount = 1;
                                done = true;
                            }
                        }
                    }
                    break;

                case 4:
                    // Augmenting path
                    bool found;
                    do
                    {
                        int row = FindStarInCol(mask, pathCol[pathCount - 1], dim);
                        if (row != -1)
                        {
                            pathRow[pathCount] = row;
                            pathCol[pathCount] = pathCol[pathCount - 1];
                            pathCount++;
                        }
                        else
                        {
                            found = false;
                            break;
                        }
                        int col = FindPrimeInRow(mask, pathRow[pathCount - 1], dim);
                        pathRow[pathCount] = pathRow[pathCount - 1];
                        pathCol[pathCount] = col;
                        pathCount++;
                        found = true;
                    } while (found);

                    // Convert path
                    for (int i = 0; i < pathCount; i++)
                        mask[pathRow[i], pathCol[i]] = mask[pathRow[i], pathCol[i]] == 1 ? 0 : 1;

                    // Clear covers and primes
                    Array.Clear(rowCover, 0, dim);
                    Array.Clear(colCover, 0, dim);
                    for (int i = 0; i < dim; i++)
                        for (int j = 0; j < dim; j++)
                            if (mask[i, j] == 2)
                                mask[i, j] = 0;
                    step = 2;
                    break;

                case 5:
                    // Add minimum uncovered value to every element of covered rows, subtract from uncovered columns
                    double minVal = double.MaxValue;
                    for (int i = 0; i < dim; i++)
                        if (rowCover[i] == 0)
                            for (int j = 0; j < dim; j++)
                                if (colCover[j] == 0 && cost[i, j] < minVal)
                                    minVal = cost[i, j];
                    for (int i = 0; i < dim; i++)
                        for (int j = 0; j < dim; j++)
                        {
                            if (rowCover[i] == 1)
                                cost[i, j] += minVal;
                            if (colCover[j] == 0)
                                cost[i, j] -= minVal;
                        }
                    step = 3;
                    break;

                case 7:
                    // Done
                    int[] result = new int[rows];
                    for (int i = 0; i < rows; i++)
                    {
                        result[i] = -1;
                        for (int j = 0; j < cols; j++)
                            if (mask[i, j] == 1)
                                result[i] = j;
                    }
                    return result;
            }
        }
    }

    private static (int, int) FindZero(double[,] cost, int[] rowCover, int[] colCover, int dim)
    {
        for (int i = 0; i < dim; i++)
            if (rowCover[i] == 0)
                for (int j = 0; j < dim; j++)
                    if (cost[i, j] == 0 && colCover[j] == 0)
                        return (i, j);
        return (-1, -1);
    }

    private static int FindStarInCol(int[,] mask, int col, int dim)
    {
        for (int i = 0; i < dim; i++)
            if (mask[i, col] == 1)
                return i;
        return -1;
    }

    private static int FindPrimeInRow(int[,] mask, int row, int dim)
    {
        for (int j = 0; j < dim; j++)
            if (mask[row, j] == 2)
                return j;
        return -1;
    }
}
