using System;
using System.Collections.Generic;
using System.Linq;
using EngineLayer;
using EngineLayer.ParallelSearch;
using EngineLayer.SpectrumMatch;
using MassSpectrometry;
using NUnit.Framework;
using Omics;
using Omics.Fragmentation;
using Proteomics;
using Proteomics.ProteolyticDigestion;

namespace Test.ParallelSearchTask.Engine;

[TestFixture]
public class TransientProteinParsimonyEngineTests
{
    [Test]
    public void NeighborhoodPsms_IncludesUpdatedAndClonedBaselineOverlap()
    {
        CommonParameters commonParameters = new(digestionParams: new DigestionParams(minPeptideLength: 1));

        Protein transientProtein = new("MPEPTIDEK", "TRANSIENT_1");
        Protein baselineProtein0 = new("MPEPTIDEK", "BASE_0");
        Protein baselineProtein1 = new("MPEPTIDEK", "BASE_1");

        SpectralMatch updatedTransientPsm = CreatePsm(transientProtein, commonParameters, scanNumber: 1001);
        SpectralMatch baselinePsm0 = CreatePsm(baselineProtein0, commonParameters, scanNumber: 2001);
        SpectralMatch baselinePsm1 = CreatePsm(baselineProtein1, commonParameters, scanNumber: 2002);

        SpectralMatch[] allPsms = [updatedTransientPsm, null];
        SpectralMatch[] baseSearchPsms = [baselinePsm0, baselinePsm1];
        HashSet<int> updatedIndexes = [0];

        string overlapKey = TransientProteinParsimonyEngine.GetParsimonyPeptideKey(
            baselinePsm1.BestMatchingBioPolymersWithSetMods.First().SpecificBioPolymer,
            modPeptidesAreDifferent: true);

        Dictionary<string, List<int>> baselineLookup = new(StringComparer.Ordinal)
        {
            { overlapKey, [1] }
        };

        var engine = new TransientProteinParsimonyEngine(
            allPsms,
            baseSearchPsms,
            updatedIndexes,
            [transientProtein],
            baselineLookup,
            new Dictionary<IBioPolymer, List<int>>(ReferenceEqualityComparer.Instance),
            modPeptidesAreDifferent: true,
            commonParameters,
            fileSpecificParameters: null,
            nestedIds: []);

        List<SpectralMatch> neighborhood = engine.NeighborhoodPsms;
        SpectralMatch clonedBaseline = neighborhood.Single(p => !ReferenceEquals(p, updatedTransientPsm));

        Assert.Multiple(() =>
        {
            Assert.That(neighborhood.Count, Is.EqualTo(2));
            Assert.That(neighborhood.Any(p => ReferenceEquals(p, updatedTransientPsm)), Is.True);
            Assert.That(ReferenceEquals(clonedBaseline, baselinePsm1), Is.False);
            Assert.That(clonedBaseline.ScanNumber, Is.EqualTo(baselinePsm1.ScanNumber));
        });
    }

    [Test]
    public void RunSpecific_FiltersProteinGroupsToTransientReferenceSet()
    {
        CommonParameters commonParameters = new(digestionParams: new DigestionParams(minPeptideLength: 1));

        Protein baselineProtein = new("MPEPTIDEK", "SHARED_ACCESSION");
        Protein transientProtein = new("MPEPTIDEK", "SHARED_ACCESSION");

        SpectralMatch transientPsm = CreatePsm(transientProtein, commonParameters, scanNumber: 1101);
        var baselinePeptide = baselineProtein.Digest(commonParameters.DigestionParams, [], []).First();

        transientPsm.AddOrReplace(
            baselinePeptide,
            transientPsm.Score,
            notch: 0,
            reportAllAmbiguity: true,
            new List<MatchedFragmentIon>());
        transientPsm.ResolveAllAmbiguities();

        SpectralMatch baselinePsm = CreatePsm(baselineProtein, commonParameters, scanNumber: 2101);

        string sharedKey = TransientProteinParsimonyEngine.GetParsimonyPeptideKey(
            baselinePeptide,
            modPeptidesAreDifferent: true);

        var engine = new TransientProteinParsimonyEngine(
            allPsms: [transientPsm],
            baseSearchPsms: [baselinePsm],
            updatedPsmIndexes: [0],
            transientProteins: [transientProtein],
            baselinePeptideKeyToScanIndexes: new Dictionary<string, List<int>>(StringComparer.Ordinal)
            {
                { sharedKey, [0] }
            },
            baselineProteinToScanIndexes: new Dictionary<IBioPolymer, List<int>>(ReferenceEqualityComparer.Instance)
            {
                { baselineProtein, [0] }
            },
            modPeptidesAreDifferent: true,
            commonParameters,
            fileSpecificParameters: null,
            nestedIds: []);

        ProteinParsimonyResults results = (ProteinParsimonyResults)engine.Run();

        Assert.Multiple(() =>
        {
            Assert.That(results.ProteinGroups.Count, Is.GreaterThan(0));
            Assert.That(results.ProteinGroups.Any(pg => pg.Proteins.Any(p => ReferenceEquals(p, baselineProtein))), Is.False);
            Assert.That(results.ProteinGroups.All(pg => pg.Proteins.All(p => ReferenceEquals(p, transientProtein))), Is.True);
        });
    }

    [Test]
    public void SharedPeptide_StaysAmbiguousWhenTheBaselineProteinHasOtherEvidence()
    {
        // human H = AAAAAAK + GGGGGGK + LLLLLLK, transient T = AAAAAAK + VVVVVVK; the AAAAAAK scan ties between H and T.
        // Full parsimony keeps both proteins and leaves that scan on H;T. The neighbourhood must bring in H's other scans
        // for the same result, or T's unique VVVVVVK makes H redundant and the shared scan reads as T only.
        CommonParameters commonParameters = new(digestionParams: new DigestionParams(minPeptideLength: 1));
        Protein human = new("AAAAAAKGGGGGGKLLLLLLK", "HUMAN_H");
        Protein transient = new("AAAAAAKVVVVVVK", "TRANSIENT_T");

        SpectralMatch baselineShared = CreatePsm(human, commonParameters, scanNumber: 1, "AAAAAAK");
        SpectralMatch[] baseSearchPsms =
        [
            baselineShared,
            CreatePsm(human, commonParameters, scanNumber: 2, "GGGGGGK"),
            CreatePsm(human, commonParameters, scanNumber: 3, "LLLLLLK"),
            null
        ];

        // the transient search ties the shared scan with T and finds T's unique peptide on scan 4
        SpectralMatch transientShared = CreatePsm(human, commonParameters, scanNumber: 1, "AAAAAAK");
        transientShared.AddOrReplace(Digest(transient, commonParameters, "AAAAAAK"), transientShared.Score, notch: 0,
            reportAllAmbiguity: true, new List<MatchedFragmentIon>());
        transientShared.ResolveAllAmbiguities();
        SpectralMatch[] allPsms =
        [
            transientShared,
            baseSearchPsms[1],
            baseSearchPsms[2],
            CreatePsm(transient, commonParameters, scanNumber: 4, "VVVVVVK")
        ];

        Dictionary<string, List<int>> peptideLookup = new(StringComparer.Ordinal);
        for (int i = 0; i < 3; i++)
        {
            peptideLookup[TransientProteinParsimonyEngine.GetParsimonyPeptideKey(
                baseSearchPsms[i].BestMatchingBioPolymersWithSetMods.First().SpecificBioPolymer, false)] = [i];
        }

        var engine = new TransientProteinParsimonyEngine(allPsms, baseSearchPsms, updatedPsmIndexes: [0, 3], [transient],
            peptideLookup,
            new Dictionary<IBioPolymer, List<int>>(ReferenceEqualityComparer.Instance) { { human, [0, 1, 2] } },
            modPeptidesAreDifferent: false, commonParameters, fileSpecificParameters: null, nestedIds: []);
        ProteinParsimonyResults results = (ProteinParsimonyResults)engine.Run();

        Assert.Multiple(() =>
        {
            // both updated scans plus H's GGGGGGK and LLLLLLK scans
            Assert.That(engine.NeighborhoodPsms, Has.Count.EqualTo(4));
            Assert.That(transientShared.BestMatchingBioPolymersWithSetMods.Select(m => m.SpecificBioPolymer.Parent.Accession),
                Is.EquivalentTo(new[] { "HUMAN_H", "TRANSIENT_T" }));
            Assert.That(results.ProteinGroups.SelectMany(pg => pg.Proteins).Select(p => p.Accession),
                Is.EquivalentTo(new[] { "TRANSIENT_T" }));
        });
    }

    private static PeptideWithSetModifications Digest(Protein protein, CommonParameters commonParameters, string baseSequence) =>
        protein.Digest(commonParameters.DigestionParams, [], []).Cast<PeptideWithSetModifications>()
            .First(p => p.BaseSequence == baseSequence);

    private static SpectralMatch CreatePsm(Protein protein, CommonParameters commonParameters, int scanNumber, string baseSequence)
    {
        var peptide = Digest(protein, commonParameters, baseSequence);
        Ms2ScanWithSpecificMass scan = CreateScan(scanNumber, commonParameters);
        SpectralMatch psm = new PeptideSpectralMatch(peptide, notch: 0, score: 10, scanNumber, scan, commonParameters,
            new List<MatchedFragmentIon>());
        psm.ResolveAllAmbiguities();
        psm.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0);
        return psm;
    }

    private static SpectralMatch CreatePsm(Protein protein, CommonParameters commonParameters, int scanNumber)
    {
        var peptide = protein.Digest(commonParameters.DigestionParams, [], []).First();
        Ms2ScanWithSpecificMass scan = CreateScan(scanNumber, commonParameters);

        SpectralMatch psm = new PeptideSpectralMatch(
            peptide,
            notch: 0,
            score: 10,
            scanNumber,
            scan,
            commonParameters,
            new List<MatchedFragmentIon>());

        psm.ResolveAllAmbiguities();
        psm.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0);
        return psm;
    }

    private static Ms2ScanWithSpecificMass CreateScan(int scanNumber, CommonParameters commonParameters)
    {
        MsDataScan msDataScan = new(
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
