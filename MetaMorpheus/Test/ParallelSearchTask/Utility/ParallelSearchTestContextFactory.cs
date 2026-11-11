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

namespace Test.ParallelSearchTask.Utility;

internal static class ParallelSearchTestContextFactory
{
    public static CommonParameters CreateCommonParameters(double qValueThreshold = 0.01, double pepQValueThreshold = 0.01)
    {
        return new CommonParameters(qValueThreshold: qValueThreshold, pepQValueThreshold: pepQValueThreshold);
    }

    public static SpectralMatch CreateSpectralMatch(
        CommonParameters commonParameters,
        bool isDecoy,
        double score,
        double psmQValue,
        double peptideQValue,
        int scanNumber = 1)
    {
        var protein = new Protein("PEPTIDE", isDecoy ? "DECOY" : "TARGET", isDecoy: isDecoy);
        var peptide = protein.Digest(new DigestionParams(minPeptideLength: 1), [], []).First();
        var scan = CreateScan(scanNumber, commonParameters);

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
        int transientPeptideCount = 0)
    {
        return new TransientDatabaseContext
        {
            DatabaseName = "TransientDb",
            TransientDatabase = null!,
            TransientProteins = [],
            TransientProteinAccessions = new HashSet<string> { "A", "B" },
            AllPsms = allPsms,
            TransientPsms = transientPsms,
            AllPeptides = allPeptides,
            TransientPeptides = transientPeptides,
            ProteinGroups = proteinGroups,
            TransientProteinGroups = transientProteinGroups,
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

        return proteinGroup;
    }

    private static Ms2ScanWithSpecificMass CreateScan(int scanNumber, CommonParameters commonParameters)
    {
        var msDataScan = new MassSpectrometry.MsDataScan(
            new MzSpectrum([1.0], [1.0], false),
            0,
            scanNumber,
            true,
            Polarity.Positive,
            double.NaN,
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
