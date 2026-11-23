using Chromatography.RetentionTimeCalibration;
using System;
using System.Collections.Generic;
using System.Linq;
using EngineLayer;
using EngineLayer.FdrAnalysis;
using MassSpectrometry;
using Omics;
using Omics.Fragmentation;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using TaskLayer.ParallelSearch.Analysis;
using EngineLayer.SpectrumMatch;

namespace Test.ParallelSearchTask.Utility;

internal static class ParallelSearchTestContextFactory
{
    /// <summary>The search's identification filter at its threshold; the q-value notch at 0.01 unless asked otherwise.</summary>
    public static CommonParameters CreateCommonParameters(FilterType identificationFilter = FilterType.QValueNotch, double filterThreshold = 0.01)
    {
        return new CommonParameters(identificationFilter: identificationFilter, filterThreshold: filterThreshold);
    }

    public static SpectralMatch CreateSpectralMatch(
        CommonParameters commonParameters,
        bool isDecoy,
        double score,
        double psmQValue,
        double peptideQValue,
        int scanNumber = 1,
        string baseSequence = "PEPTIDE",
        double retentionTime = double.NaN)
    {
        var protein = new Protein(baseSequence, isDecoy ? "DECOY" : "TARGET", isDecoy: isDecoy);
        var peptide = protein.Digest(new DigestionParams(minPeptideLength: 1), [], []).First();
        var scan = CreateScan(scanNumber, commonParameters, retentionTime);

        var psm = new PeptideSpectralMatch(peptide, 0, score, scanNumber, scan, commonParameters, new List<MatchedFragmentIon>());
        psm.ResolveAllAmbiguities();
        psm.SetFdrValues(0, 0, psmQValue, 0, 0, psmQValue, 0, psmQValue);
        psm.PeptideFdrInfo = new FdrInfo
        {
            QValue = peptideQValue,
            QValueNotch = peptideQValue,
            PEP_QValue = peptideQValue,
        };
        return psm;
    }

    public static TransientDatabaseContext CreateContext(
        CommonParameters commonParameters,
        List<SpectralMatch> allPsms,
        List<SpectralMatch> transientPsms,
        List<SpectralMatch> allPeptides,
        List<SpectralMatch> transientPeptides,
        List<ProteinGroup>? proteinGroups = null,
        List<ProteinGroup>? transientProteinGroups = null,
        int totalProteins = 0,
        int transientPeptideCount = 0,
        IReadOnlyDictionary<string, IrtCalibrationModel>? retentionTimeCalibrations = null,
        HashSet<string>? transientProteinAccessions = null,
        BaselineTies? baselineTies = null)
    {
        // by default every protein behind the transient matches and groups is the transient database's own, so none is
        // shared with the baseline; pass the accessions to make some of them baseline proteins
        // (some tests pass null lists on purpose)
        transientProteinAccessions ??= (transientPsms ?? []).Concat(transientPeptides ?? [])
            .SelectMany(m => m.BestMatchingBioPolymersWithSetMods.Select(h => h.SpecificBioPolymer.Parent.Accession))
            .Concat((transientProteinGroups ?? []).SelectMany(g => g.Proteins.Select(p => p.Accession)))
            .ToHashSet();

        return new TransientDatabaseContext
        {
            DatabaseName = "TransientDb",
            TransientDatabase = null!,
            TransientProteins = [],
            TransientProteinAccessions = transientProteinAccessions,
            BaselineTies = baselineTies ?? BaselineTies.Empty,
            AllPsms = allPsms,
            TransientPsms = transientPsms,
            AllPeptides = allPeptides,
            TransientPeptides = transientPeptides,
            ProteinGroups = proteinGroups,
            TransientProteinGroups = transientProteinGroups,
            RetentionTimeCalibrations = retentionTimeCalibrations,
            CommonParameters = commonParameters,
            TotalProteins = totalProteins,
            TransientPeptideCount = transientPeptideCount,
            OutputFolder = string.Empty,
            NestedIds = [],
        };
    }

    public static ProteinGroup CreateProteinGroup(bool isDecoy, double qValue, int peptideCount)
    {
        var proteins = new HashSet<Omics.IBioPolymer>();
        var peptides = new HashSet<Omics.IBioPolymerWithSetMods>();

        for (int i = 0; i < peptideCount; i++)
        {
            var protein = new Protein($"PEPTIDE{i}", isDecoy ? $"DECOY_PG_{i}" : $"TARGET_PG_{i}", isDecoy: isDecoy);
            proteins.Add(protein);
            var digestedPeptides = protein.Digest(new DigestionParams(minPeptideLength: 1), [], []).ToList();
            peptides.Add(digestedPeptides[0]);
        }

        var proteinGroup = new ProteinGroup(proteins, peptides, new HashSet<Omics.IBioPolymerWithSetMods>(peptides))
        {
            QValue = qValue,
        };

        // one PSM on the group's own first peptide, so the group has evidence of its own
        var commonParameters = CreateCommonParameters();
        var psm = new PeptideSpectralMatch(peptides.First(), 0, 10, 1, CreateScan(1, commonParameters), commonParameters,
            new List<MatchedFragmentIon>());
        psm.ResolveAllAmbiguities();
        proteinGroup.AllPsmsBelowOnePercentFDR.Add(psm);

        return proteinGroup;
    }

    internal static Ms2ScanWithSpecificMass CreateScan(int scanNumber, CommonParameters commonParameters, double retentionTime = double.NaN)
    {
        var msDataScan = new MassSpectrometry.MsDataScan(
            new MzSpectrum([1.0], [1.0], false),
            0,
            scanNumber,
            true,
            Polarity.Positive,
            retentionTime,
            null,
            null,
            MZAnalyzerType.Orbitrap,
            double.NaN,
            null,
            null,
            $"scan={scanNumber}",
            double.NaN,
            null,
            null,
            double.NaN,
            null,
            MassSpectrometry.DissociationType.AnyActivationType,
            0,
            null);

        return new Ms2ScanWithSpecificMass(msDataScan, 2, 0, "TestFile", commonParameters);
    }
}
