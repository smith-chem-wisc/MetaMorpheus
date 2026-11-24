using System.Collections.Generic;
using System.Linq;
using EngineLayer;
using NUnit.Framework;
using Omics.Fragmentation;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using TaskLayer.ParallelSearch.Analysis;
using TaskLayer.ParallelSearch.Analysis.Collectors;
using Test.ParallelSearchTask.Utility;

namespace Test.ParallelSearchTask.Analysis;

/// <summary>
/// Transient parsimony trims hypotheses in place. When a transient protein explains a peptide that the only baseline
/// protein holding it also explains, parsimony can keep the transient protein and drop the baseline one, and the match
/// then looks like the database's own evidence (PXD061094 sample 103: HLA-A's FIAVGYVDDTQFVR, 22 PSMs). The ties are read
/// before parsimony, so the match still counts as shared.
/// </summary>
[TestFixture]
public class BaselineTiesTests
{
    private static readonly HashSet<string> TransientAccessions = ["TARGET", "DECOY"];

    /// <summary>A transient PSM (IVTDITK) that ties with a baseline protein's I/L isomer (LVTDLTK).</summary>
    private static SpectralMatch TiedPsm(CommonParameters commonParameters, int scanNumber)
    {
        var psm = ParallelSearchTestContextFactory.CreateSpectralMatch(commonParameters, isDecoy: false, score: 120,
            psmQValue: 0.001, peptideQValue: 0.001, scanNumber: scanNumber, baseSequence: "IVTDITK");
        var host = new Protein("LVTDLTK", "HOST");
        var hostPeptide = host.Digest(new DigestionParams(minPeptideLength: 1), [], []).First();
        psm.AddOrReplace(hostPeptide, psm.Score, 0, true, new List<MatchedFragmentIon>());
        return psm;
    }

    /// <summary>What transient parsimony does when it keeps the transient protein: drops the baseline hypothesis.</summary>
    private static void DropBaselineHypotheses(SpectralMatch psm)
    {
        foreach (var hypothesis in psm.BestMatchingBioPolymersWithSetMods
                     .Where(h => !TransientAccessions.Contains(h.SpecificBioPolymer.Parent.Accession)).ToList())
            psm.RemoveThisAmbiguousPeptide(hypothesis);
    }

    [Test]
    public void TieTrimmedByParsimony_StillCountsAsShared()
    {
        var commonParameters = ParallelSearchTestContextFactory.CreateCommonParameters();
        var own = ParallelSearchTestContextFactory.CreateSpectralMatch(commonParameters, isDecoy: false, score: 110,
            psmQValue: 0.005, peptideQValue: 0.005, scanNumber: 1);
        var tied = TiedPsm(commonParameters, scanNumber: 2);

        var ties = BaselineTies.Read([own, tied], TransientAccessions);
        DropBaselineHypotheses(tied);
        Assert.That(tied.BestMatchingBioPolymersWithSetMods.Select(h => h.SpecificBioPolymer.Parent.Accession),
            Is.EqualTo(new[] { "TARGET" }), "parsimony left only the transient protein");

        var withoutTies = ParallelSearchTestContextFactory.CreateContext(commonParameters,
            [own, tied], [own, tied], [own, tied], [own, tied], transientProteinAccessions: TransientAccessions);
        var withTies = ParallelSearchTestContextFactory.CreateContext(commonParameters,
            [own, tied], [own, tied], [own, tied], [own, tied], transientProteinAccessions: TransientAccessions,
            baselineTies: ties);

        var results = new PsmPeptideSearchCollector().CollectData(withTies);

        Assert.Multiple(() =>
        {
            Assert.That(withoutTies.IsSharedWithBaseline(tied), Is.False, "the trimmed match alone looks like the database's own");
            Assert.That(withTies.IsSharedWithBaseline(tied), Is.True);
            Assert.That(withTies.IsSharedWithBaseline(own), Is.False);
            Assert.That(withTies.ConfidentUnambiguousPsms, Is.EqualTo(new[] { own }));
            Assert.That(withTies.ConfidentUnambiguousPeptides, Is.EqualTo(new[] { own }));
            Assert.That(results[PsmPeptideSearchCollector.PsmBacterialUnambiguousTargets], Is.EqualTo(1));
            Assert.That(results[PsmPeptideSearchCollector.PeptideBacterialUnambiguousTargets], Is.EqualTo(1));
        });
    }

    [Test]
    public void BaselineIsobaricPeptide_IsSharedInItsOtherMatchesToo()
    {
        var commonParameters = ParallelSearchTestContextFactory.CreateCommonParameters();
        var tied = TiedPsm(commonParameters, scanNumber: 1);
        // the same transient peptide on another scan, where no baseline hypothesis was kept
        var sameSequence = ParallelSearchTestContextFactory.CreateSpectralMatch(commonParameters, isDecoy: false, score: 100,
            psmQValue: 0.002, peptideQValue: 0.002, scanNumber: 2, baseSequence: "IVTDITK");
        var other = ParallelSearchTestContextFactory.CreateSpectralMatch(commonParameters, isDecoy: false, score: 90,
            psmQValue: 0.003, peptideQValue: 0.003, scanNumber: 3, baseSequence: "SAMPLEK");

        var ties = BaselineTies.Read([tied, sameSequence, other], TransientAccessions);

        Assert.Multiple(() =>
        {
            Assert.That(ties.Contains(sameSequence, TransientAccessions), Is.True);
            Assert.That(ties.Contains(other, TransientAccessions), Is.False);
        });
    }

    [Test]
    public void TieWithADifferentBaselinePeptide_DoesNotMarkTheSequence()
    {
        var commonParameters = ParallelSearchTestContextFactory.CreateCommonParameters();
        // a transient peptide tying by chance with an unrelated baseline peptide: that match is shared, but the sequence
        // is not baseline-isobaric, so its other matches stay the database's own
        var psm = ParallelSearchTestContextFactory.CreateSpectralMatch(commonParameters, isDecoy: false, score: 120,
            psmQValue: 0.001, peptideQValue: 0.001, scanNumber: 1, baseSequence: "IVTDITK");
        var host = new Protein("TKDTIVK", "HOST");
        psm.AddOrReplace(host.Digest(new DigestionParams(minPeptideLength: 1), [], []).First(), psm.Score, 0, true,
            new List<MatchedFragmentIon>());
        var elsewhere = ParallelSearchTestContextFactory.CreateSpectralMatch(commonParameters, isDecoy: false, score: 100,
            psmQValue: 0.002, peptideQValue: 0.002, scanNumber: 2, baseSequence: "IVTDITK");

        var ties = BaselineTies.Read([psm, elsewhere], TransientAccessions);

        Assert.Multiple(() =>
        {
            Assert.That(ties.Contains(psm, TransientAccessions), Is.True);
            Assert.That(ties.Contains(elsewhere, TransientAccessions), Is.False);
        });
    }
}
