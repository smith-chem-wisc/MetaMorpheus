using EngineLayer.SpectrumMatch;
using EngineLayer.FdrAnalysis;
using Omics;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;

namespace EngineLayer.ParallelSearch;

/// <summary>
/// Runs parsimony on a neighbourhood around the scans the transient search changed: those scans, the baseline scans
/// of the same peptides, and every baseline scan of each baseline protein those PSMs map to, so a baseline protein
/// keeps all its evidence and is not dropped for a transient protein that only shares one of its peptides.
/// </summary>
public sealed class TransientProteinParsimonyEngine : ProteinParsimonyEngine
{
    private readonly HashSet<IBioPolymer> _transientProteins;

    public List<SpectralMatch> NeighborhoodPsms { get; }

    public static string GetParsimonyPeptideKey(IBioPolymerWithSetMods peptide, bool modPeptidesAreDifferent)
    {
        string sequence = modPeptidesAreDifferent ? peptide.FullSequence : peptide.BaseSequence;
        return $"{peptide.DigestionParams.DigestionAgent}|{sequence}";
    }

    public TransientProteinParsimonyEngine(
        SpectralMatch[] allPsms,
        SpectralMatch[] baseSearchPsms,
        HashSet<int> updatedPsmIndexes,
        List<IBioPolymer> transientProteins,
        Dictionary<string, List<int>> baselinePeptideKeyToScanIndexes,
        Dictionary<IBioPolymer, List<int>> baselineProteinToScanIndexes,
        bool modPeptidesAreDifferent,
        CommonParameters commonParameters,
        List<(string fileName, CommonParameters fileSpecificParameters)> fileSpecificParameters,
        List<string> nestedIds)
        : this(
            BuildState(
                allPsms,
                baseSearchPsms,
                updatedPsmIndexes,
                transientProteins,
                baselinePeptideKeyToScanIndexes,
                baselineProteinToScanIndexes,
                modPeptidesAreDifferent),
            modPeptidesAreDifferent,
            commonParameters,
            fileSpecificParameters,
            nestedIds)
    {
    }

    private TransientProteinParsimonyEngine(
        NeighborhoodBuildState state,
        bool modPeptidesAreDifferent,
        CommonParameters commonParameters,
        List<(string fileName, CommonParameters fileSpecificParameters)> fileSpecificParameters,
        List<string> nestedIds)
        : base(state.NeighborhoodPsms, state.NeighborhoodPsms, modPeptidesAreDifferent, commonParameters, fileSpecificParameters, nestedIds)
    {
        _transientProteins = state.TransientProteins;
        NeighborhoodPsms = state.NeighborhoodPsms;
    }

    protected override MetaMorpheusEngineResults RunSpecific()
    {
        Status($"Transient parsimony neighborhood: {NeighborhoodPsms.Count} PSMs ({_transientProteins.Count} transient proteins)");

        ProteinParsimonyResults results = (ProteinParsimonyResults)base.RunSpecific();
        results.ProteinGroups = FilterProteinGroupsToTransientProteins(results.ProteinGroups, _transientProteins).ToList();

        return results;
    }

    private static NeighborhoodBuildState BuildState(
        SpectralMatch[] allPsms,
        SpectralMatch[] baseSearchPsms,
        HashSet<int> updatedPsmIndexes,
        List<IBioPolymer> transientProteins,
        Dictionary<string, List<int>> baselinePeptideKeyToScanIndexes,
        Dictionary<IBioPolymer, List<int>> baselineProteinToScanIndexes,
        bool modPeptidesAreDifferent)
    {
        ArgumentNullException.ThrowIfNull(allPsms);
        ArgumentNullException.ThrowIfNull(baseSearchPsms);
        ArgumentNullException.ThrowIfNull(updatedPsmIndexes);
        ArgumentNullException.ThrowIfNull(transientProteins);
        ArgumentNullException.ThrowIfNull(baselinePeptideKeyToScanIndexes);
        ArgumentNullException.ThrowIfNull(baselineProteinToScanIndexes);

        HashSet<IBioPolymer> transientProteinSet = new(transientProteins, ReferenceComparer<IBioPolymer>.Instance);
        List<SpectralMatch> neighborhoodPsms = BuildParsimonyNeighborhoodPsms(
            allPsms,
            baseSearchPsms,
            updatedPsmIndexes,
            transientProteinSet,
            baselinePeptideKeyToScanIndexes,
            baselineProteinToScanIndexes,
            modPeptidesAreDifferent);

        return new NeighborhoodBuildState(transientProteinSet, neighborhoodPsms);
    }

    private static List<SpectralMatch> BuildParsimonyNeighborhoodPsms(
        SpectralMatch[] allPsms,
        SpectralMatch[] baseSearchPsms,
        HashSet<int> updatedPsmIndexes,
        HashSet<IBioPolymer> transientProteins,
        Dictionary<string, List<int>> baselinePeptideKeyToScanIndexes,
        Dictionary<IBioPolymer, List<int>> baselineProteinToScanIndexes,
        bool modPeptidesAreDifferent)
    {
        List<SpectralMatch> neighborhoodPsms = new(updatedPsmIndexes.Count * 2);
        HashSet<int> includedScanIndexes = new();
        HashSet<string> neighborhoodPeptideKeys = new(StringComparer.Ordinal);
        // baseline proteins in discovery order, so the neighbourhood is built the same way every run
        List<IBioPolymer> neighborhoodBaselineProteins = new();
        HashSet<IBioPolymer> seenBaselineProteins = new(ReferenceComparer<IBioPolymer>.Instance);

        void AddBaselineProteins(SpectralMatch psm)
        {
            foreach (var hypothesis in psm.BestMatchingBioPolymersWithSetMods)
            {
                IBioPolymer protein = hypothesis.SpecificBioPolymer.Parent;
                if (!transientProteins.Contains(protein) && seenBaselineProteins.Add(protein))
                {
                    neighborhoodBaselineProteins.Add(protein);
                }
            }
        }

        void AddBaselineScan(int scanIndex)
        {
            if (scanIndex < 0 || scanIndex >= baseSearchPsms.Length || !includedScanIndexes.Add(scanIndex))
            {
                return;
            }

            var baselinePsm = baseSearchPsms[scanIndex];
            if (baselinePsm is null)
            {
                return;
            }

            var clonedBaselinePsm = ClonePsmForParsimony(baselinePsm);
            if (clonedBaselinePsm is not null)
            {
                neighborhoodPsms.Add(clonedBaselinePsm);
            }
        }

        foreach (int scanIndex in updatedPsmIndexes.OrderBy(p => p))
        {
            if (scanIndex < 0 || scanIndex >= allPsms.Length)
            {
                continue;
            }

            SpectralMatch psm = allPsms[scanIndex];
            if (psm is null)
            {
                continue;
            }

            SpectralMatch originalBaselinePsm = scanIndex >= 0 && scanIndex < baseSearchPsms.Length
                ? baseSearchPsms[scanIndex]
                : null;

            foreach (var hypothesis in psm.BestMatchingBioPolymersWithSetMods)
            {
                if (!transientProteins.Contains(hypothesis.SpecificBioPolymer.Parent) && originalBaselinePsm is not null)
                {
                    foreach (var originalHypothesis in originalBaselinePsm.BestMatchingBioPolymersWithSetMods)
                    {
                        neighborhoodPeptideKeys.Add(GetParsimonyPeptideKey(originalHypothesis.SpecificBioPolymer, modPeptidesAreDifferent));
                    }
                }
                else
                {
                    neighborhoodPeptideKeys.Add(GetParsimonyPeptideKey(hypothesis.SpecificBioPolymer, modPeptidesAreDifferent));
                }
            }

            neighborhoodPsms.Add(psm);
            includedScanIndexes.Add(scanIndex);
            AddBaselineProteins(psm);
            if (originalBaselinePsm is not null)
            {
                AddBaselineProteins(originalBaselinePsm);
            }
        }

        int firstBaselineClone = neighborhoodPsms.Count;
        foreach (string peptideKey in neighborhoodPeptideKeys)
        {
            if (!baselinePeptideKeyToScanIndexes.TryGetValue(peptideKey, out var baselineScanIndexes))
            {
                continue;
            }

            foreach (int scanIndex in baselineScanIndexes)
            {
                AddBaselineScan(scanIndex);
            }
        }

        // every baseline protein matched so far brings all its baseline scans, so it competes in parsimony with its full
        // support rather than with only the peptide it shares. Proteins reached only through those added scans are not
        // expanded again; they bear on how baseline peptides are shared among baseline proteins, and groups without a
        // transient protein are filtered out after parsimony.
        for (int i = firstBaselineClone; i < neighborhoodPsms.Count; i++)
        {
            AddBaselineProteins(neighborhoodPsms[i]);
        }

        foreach (IBioPolymer baselineProtein in neighborhoodBaselineProteins)
        {
            if (!baselineProteinToScanIndexes.TryGetValue(baselineProtein, out var baselineScanIndexes))
            {
                continue;
            }

            foreach (int scanIndex in baselineScanIndexes)
            {
                AddBaselineScan(scanIndex);
            }
        }

        return neighborhoodPsms;
    }

    private static SpectralMatch ClonePsmForParsimony(SpectralMatch source)
    {
        var bestMatches = source.BestMatchingBioPolymersWithSetMods.ToList();
        if (bestMatches.Count == 0)
        {
            return null;
        }

        SpectralMatch clone = source switch
        {
            PeptideSpectralMatch peptidePsm => peptidePsm.Clone(bestMatches),
            _ => null
        };

        if (clone is null)
        {
            return null;
        }

        clone.PsmFdrInfo = source.PsmFdrInfo?.Clone() ?? new FdrInfo();
        clone.PeptideFdrInfo = source.PeptideFdrInfo?.Clone() ?? new FdrInfo();
        clone.ResolveAllAmbiguities();

        return clone;
    }

    private static IEnumerable<ProteinGroup> FilterProteinGroupsToTransientProteins(List<ProteinGroup> proteinGroups,
        HashSet<IBioPolymer> transientProteins)
    {
        return proteinGroups.Where(pg => pg.Proteins.Any(transientProteins.Contains));
    }

    private sealed class NeighborhoodBuildState(HashSet<IBioPolymer> transientProteins, List<SpectralMatch> neighborhoodPsms)
    {
        public HashSet<IBioPolymer> TransientProteins { get; } = transientProteins;
        public List<SpectralMatch> NeighborhoodPsms { get; } = neighborhoodPsms;
    }

    private sealed class ReferenceComparer<T> : IEqualityComparer<T> where T : class
    {
        public static ReferenceComparer<T> Instance { get; } = new();

        public bool Equals(T x, T y)
        {
            return ReferenceEquals(x, y);
        }

        public int GetHashCode(T obj)
        {
            return RuntimeHelpers.GetHashCode(obj);
        }
    }
}
