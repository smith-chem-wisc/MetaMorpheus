using System.Collections.Generic;
using System.Linq;
using EngineLayer;
using NUnit.Framework;
using Omics.Fragmentation;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using TaskLayer.ParallelSearch.Analysis.Collectors;
using Test.ParallelSearchTask.Utility;

namespace Test.ParallelSearchTask.Analysis;

/// <summary>
/// A transient match that ties with a protein outside the transient database is a baseline identification, not evidence
/// for the database. On PXD038928 one phage peptide, an I/L isomer of a human/cRAP peptide, gave 748 such PSMs; the
/// contaminant has no "Homo sapiens" organism, so only the transient-accession check catches it.
/// </summary>
[TestFixture]
public class BaselineSharedEvidenceTests
{
    private static readonly HashSet<string> TransientAccessions = ["TARGET", "DECOY"];

    /// <summary>A confident transient PSM that ties with a contaminant protein of no stated organism.</summary>
    private static SpectralMatch SharedPsm(CommonParameters commonParameters, int scanNumber)
    {
        var psm = ParallelSearchTestContextFactory.CreateSpectralMatch(commonParameters, isDecoy: false, score: 120,
            psmQValue: 0.001, peptideQValue: 0.001, scanNumber: scanNumber, baseSequence: "IVTDITK");
        var contaminant = new Protein("LVTDLTK", "CONTAMINANT");
        var baselinePeptide = contaminant.Digest(new DigestionParams(minPeptideLength: 1), [], []).First();
        psm.AddOrReplace(baselinePeptide, psm.Score, 0, true, new List<MatchedFragmentIon>());
        return psm;
    }

    [Test]
    public void PsmCollector_CountsAndDistributions_LeaveOutMatchesSharedWithTheBaseline()
    {
        var commonParameters = ParallelSearchTestContextFactory.CreateCommonParameters();
        var own = ParallelSearchTestContextFactory.CreateSpectralMatch(commonParameters, isDecoy: false, score: 110,
            psmQValue: 0.005, peptideQValue: 0.005, scanNumber: 1);
        var shared = SharedPsm(commonParameters, scanNumber: 2);
        Assert.That(shared.BestMatchingBioPolymersWithSetMods.Count(), Is.EqualTo(2), "the contaminant ties");

        var context = ParallelSearchTestContextFactory.CreateContext(commonParameters,
            [own, shared], [own, shared], [own, shared], [own, shared],
            transientProteinAccessions: TransientAccessions);

        var results = new PsmPeptideSearchCollector().CollectData(context);

        Assert.Multiple(() =>
        {
            Assert.That(context.IsSharedWithBaseline(shared), Is.True);
            Assert.That(context.IsSharedWithBaseline(own), Is.False);
            Assert.That(context.ConfidentUnambiguousPsms, Is.EqualTo(new[] { own }));
            Assert.That(context.ConfidentUnambiguousPeptides, Is.EqualTo(new[] { own }));

            Assert.That(results[PsmPeptideSearchCollector.PsmBacterialTargets], Is.EqualTo(2));
            Assert.That(results[PsmPeptideSearchCollector.PsmBacterialUnambiguousTargets], Is.EqualTo(1));
            Assert.That(results[PsmPeptideSearchCollector.PsmBacterialAmbiguous], Is.EqualTo(1));
            Assert.That(results[PsmPeptideSearchCollector.PeptideBacterialUnambiguousTargets], Is.EqualTo(1));
            Assert.That((double[])results[PsmPeptideSearchCollector.PsmBacterialUnambiguousTargetScores], Is.EqualTo(new[] { 110d }));
            Assert.That((double[])results[PsmPeptideSearchCollector.PsmBacteriaTargetlDeltaScores], Has.Length.EqualTo(1));
            Assert.That((double[])results[PsmPeptideSearchCollector.PsmPrecursorMassErrors], Has.Length.EqualTo(1));
        });
    }

    [Test]
    public void ProteinGroupCollector_UsesOnlyConfidentGroupsOfTheDatabaseItself()
    {
        var commonParameters = ParallelSearchTestContextFactory.CreateCommonParameters();
        var confident = ParallelSearchTestContextFactory.CreateProteinGroup(isDecoy: false, qValue: 0.005, peptideCount: 2);
        var unconfident = ParallelSearchTestContextFactory.CreateProteinGroup(isDecoy: false, qValue: 0.5, peptideCount: 3);
        var withBaselineProtein = ParallelSearchTestContextFactory.CreateProteinGroup(isDecoy: false, qValue: 0.001, peptideCount: 4);
        withBaselineProtein.Proteins.Add(new Protein("PEPTIDEK", "BASELINE"));
        // the database's own proteins, but its only PSM ties with a contaminant
        var onlySharedEvidence = ParallelSearchTestContextFactory.CreateProteinGroup(isDecoy: false, qValue: 0.002, peptideCount: 3);
        onlySharedEvidence.AllPsmsBelowOnePercentFDR.Clear();
        onlySharedEvidence.AllPsmsBelowOnePercentFDR.Add(SharedPsm(commonParameters, scanNumber: 3));

        var groups = new List<ProteinGroup> { confident, unconfident, withBaselineProtein, onlySharedEvidence };
        var context = ParallelSearchTestContextFactory.CreateContext(commonParameters, [], [], [], [],
            proteinGroups: groups, transientProteinGroups: groups,
            transientProteinAccessions: ["TARGET", "TARGET_PG_0", "TARGET_PG_1", "TARGET_PG_2", "TARGET_PG_3"]);

        var results = new ProteinGroupCollector().CollectData(context);

        Assert.Multiple(() =>
        {
            Assert.That(context.IsSharedWithBaseline(onlySharedEvidence), Is.True);
            Assert.That(context.IsSharedWithBaseline(confident), Is.False);
            Assert.That(results[ProteinGroupCollector.TargetProteinGroupsFromTransientDbAtQValueThreshold], Is.EqualTo(3));
            Assert.That(results[ProteinGroupCollector.ProteinGroupBacterialUnambiguousTargets], Is.EqualTo(1));
            // only the confident group of the database's own proteins, with its 2 peptides
            Assert.That((double[])results[ProteinGroupCollector.AllPeptidesPerProteinGroup], Is.EqualTo(new[] { 2d }));
        });
    }
}
