using NUnit.Framework;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using Chemistry;
using EngineLayer;
using MassSpectrometry;
using Omics.Modifications;
using System.Linq;
using EngineLayer.FdrAnalysis;

namespace Test;

[TestFixture]
public class GptmdFilterTests
{
    // Helper to create a dummy PeptideWithSetModifications
    public static PeptideWithSetModifications DummyPeptide(bool isDecoy = false) => new("PEPTIDE", [], p: new Proteomics.Protein("PEPTIDE", "a", isDecoy: isDecoy));
    private static Modification DummyMod() => new("a", "a", "a", "a", DummyMotif(), "Anywhere.", null, 20);
    private static ModificationMotif DummyMotif() => ModificationMotif.TryGetMotif("X", out var motif) ? motif : null;
    // Helper to create a dummy SpectralMatch
    public static SpectralMatch DummySpectralMatch(double score = 5, int notch = 0, bool isDecoy=false) => new PeptideSpectralMatch(DummyPeptide(isDecoy), notch, score, 0,
        new Ms2ScanWithSpecificMass(
            new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true,
                MassSpectrometry.Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN,
                null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType,
                0, null),
            (new Proteomics.AminoAcidPolymer.Peptide(DummyPeptide().BaseSequence).MonoisotopicMass + 21.981943)
            .ToMz(1), 1, "filepath", new CommonParameters())
        ,
        new(), new List<MatchedFragmentIon>())
    {
        PsmFdrInfo = new FdrInfo(), PeptideFdrInfo = new FdrInfo()
    };

    // Helper to create a dummy MatchedFragmentIon
    private MatchedFragmentIon CreateIon(FragmentationTerminus terminus, int fragmentNumber, int residuePosition)
    {
        var product = new Product(ProductType.b, terminus, 0, fragmentNumber, residuePosition, 0, 0);
        return new MatchedFragmentIon(product, 100, 100, 1);
    }

    [Test]
    public void ImprovedScoreFilter_Passes_ReturnsTrueIfNewScoreGreater()
    {
        var filter = new ImprovedScoreFilter();
        bool result = filter.Passes(
            DummyPeptide(),
            DummySpectralMatch(),
            newScore: 2.0,
            originalScore: 1.0,
            matchedIons: null,
            peptideOneBasedModSite: 1,
            peptideLength: 7,
            DummyMod());

        Assert.That(result, Is.True);
        Assert.That(IGptmdFilter.GetFilterTypeName(filter), Is.EqualTo("ImprovedScoreFilter"));
    }

    [Test]
    public void ImprovedScoreFilter_Passes_ReturnsFalseIfNewScoreNotGreater()
    {
        var filter = new ImprovedScoreFilter();
        bool result = filter.Passes(
            DummyPeptide(),
            DummySpectralMatch(),
            newScore: 1.0,
            originalScore: 2.0,
            matchedIons: null,
            peptideOneBasedModSite: 1,
            peptideLength: 7,
            DummyMod());

        Assert.That(result, Is.False);
    }

    [Test]
    public void DualDirectionalIonCoverageFilter_Passes_ReturnsFalseIfNoMatchedIons()
    {
        var filter = new DualDirectionalIonCoverageFilter();
        bool result = filter.Passes(
            DummyPeptide(),
            DummySpectralMatch(),
            0, 0,
            matchedIons: null,
            peptideOneBasedModSite: 3,
            peptideLength: 7,
            DummyMod());

        Assert.That(result, Is.False);
        Assert.That(IGptmdFilter.GetFilterTypeName(filter), Is.EqualTo("DualDirectionalIonCoverageFilter"));
    }

    [Test]
    public void DualDirectionalIonCoverageFilter_Passes_ReturnsTrueIfBothDirectionsCovered()
    {
        var filter = new DualDirectionalIonCoverageFilter();
        var ions = new List<MatchedFragmentIon>
        {
            CreateIon(FragmentationTerminus.N, fragmentNumber: 5, residuePosition: 5), // covers N-term
            CreateIon(FragmentationTerminus.C, fragmentNumber: 6, residuePosition: 1)  // covers C-term
        };

        bool result = filter.Passes(
            DummyPeptide(),
            DummySpectralMatch(),
            0, 0,
            ions,
            peptideOneBasedModSite: 3,
            peptideLength: 7,
            DummyMod());

        Assert.That(result, Is.True);
    }

    [Test]
    public void DualDirectionalIonCoverageFilter_Passes_ReturnsFalseIfOnlyOneDirectionCovered()
    {
        var filter = new DualDirectionalIonCoverageFilter();
        var ions = new List<MatchedFragmentIon>
        {
            CreateIon(FragmentationTerminus.N, fragmentNumber: 3, residuePosition: 2)
        };

        bool result = filter.Passes(
            DummyPeptide(),
            DummySpectralMatch(),
            0, 0,
            ions,
            peptideOneBasedModSite: 3,
            peptideLength: 7,
            DummyMod());

        Assert.That(result, Is.False);
    }

    [Test]
    public void FlankingIonCoverageFilter_Passes_ReturnsFalseIfNoMatchedIons()
    {
        var filter = new FlankingIonCoverageFilter();
        bool result = filter.Passes(
            DummyPeptide(),
            DummySpectralMatch(),
            0, 0,
            matchedIons: null,
            peptideOneBasedModSite: 3,
            peptideLength: 7,
            DummyMod());

        Assert.That(result, Is.False);
        Assert.That(IGptmdFilter.GetFilterTypeName(filter), Is.EqualTo("FlankingIonCoverageFilter"));
    }

    [Test]
    public void FlankingIonCoverageFilter_Passes_ReturnsTrueIfBothFlanksCovered()
    {
        var filter = new FlankingIonCoverageFilter();
        var ions = new List<MatchedFragmentIon>
        {
            CreateIon(FragmentationTerminus.N, fragmentNumber: 2, residuePosition: 2), // left flank
            CreateIon(FragmentationTerminus.C, fragmentNumber: 3, residuePosition: 3)  // right flank
        };

        bool result = filter.Passes(
            DummyPeptide(),
            DummySpectralMatch(),
            0, 0,
            ions,
            peptideOneBasedModSite: 3,
            peptideLength: 7,
            DummyMod());

        Assert.That(result, Is.True);
    }

    [Test]
    public void FlankingIonCoverageFilter_Passes_ReturnsFalseIfOnlyLeftFlankCovered()
    {
        var filter = new FlankingIonCoverageFilter();
        var ions = new List<MatchedFragmentIon>
        {
            CreateIon(FragmentationTerminus.N, fragmentNumber: 2, residuePosition: 2)
        };

        bool result = filter.Passes(
            DummyPeptide(),
            DummySpectralMatch(),
            0, 0,
            ions,
            peptideOneBasedModSite: 3,
            peptideLength: 7,
            DummyMod());

        Assert.That(result, Is.False);
    }

    [Test]
    public void FlankingIonCoverageFilter_Passes_ReturnsFalseIfOnlyRightFlankCovered()
    {
        var filter = new FlankingIonCoverageFilter();
        var ions = new List<MatchedFragmentIon>
        {
            CreateIon(FragmentationTerminus.C, fragmentNumber: 3, residuePosition: 3)
        };

        bool result = filter.Passes(
            DummyPeptide(),
            DummySpectralMatch(),
            0, 0,
            ions,
            peptideOneBasedModSite: 3,
            peptideLength: 7,
            DummyMod());

        Assert.That(result, Is.False);
    }

    [Test]
    public void ImprovedScoreFilter_Passes_ReturnsFalseIfScoresAreEqual()
    {
        var filter = new ImprovedScoreFilter();
        bool result = filter.Passes(
            DummyPeptide(),
            DummySpectralMatch(),
            newScore: 1.0,
            originalScore: 1.0,
            matchedIons: null,
            peptideOneBasedModSite: 1,
            peptideLength: 7,
            DummyMod());

        Assert.That(result, Is.False);
    }

    [Test]
    public void DualDirectionalIonCoverageFilter_Passes_ReturnsFalseIfEmptyMatchedIons()
    {
        var filter = new DualDirectionalIonCoverageFilter();
        bool result = filter.Passes(
            DummyPeptide(),
            DummySpectralMatch(),
            0, 0,
            new List<MatchedFragmentIon>(),
            peptideOneBasedModSite: 3,
            peptideLength: 7,
            DummyMod());

        Assert.That(result, Is.False);
        Assert.That(IGptmdFilter.GetFilterTypeName(filter), Is.EqualTo("DualDirectionalIonCoverageFilter"));
    }

    [Test]
    public void DualDirectionalIonCoverageFilter_Passes_ReturnsTrueIfFivePrimeAndThreePrimeCovered()
    {
        var filter = new DualDirectionalIonCoverageFilter();
        var ions = new List<MatchedFragmentIon>
        {
            CreateIon(FragmentationTerminus.FivePrime, fragmentNumber: 4, residuePosition: 4), // N-term
            CreateIon(FragmentationTerminus.ThreePrime, fragmentNumber: 5, residuePosition: 1) // C-term
        };

        bool result = filter.Passes(
            DummyPeptide(),
            DummySpectralMatch(),
            0, 0,
            ions,
            peptideOneBasedModSite: 3,
            peptideLength: 7,
            DummyMod());

        Assert.That(result, Is.True);
    }

    [Test]
    public void FlankingIonCoverageFilter_Passes_ReturnsFalseIfEmptyMatchedIons()
    {
        var filter = new FlankingIonCoverageFilter();
        bool result = filter.Passes(
            DummyPeptide(),
            DummySpectralMatch(),
            0, 0,
            new List<MatchedFragmentIon>(),
            peptideOneBasedModSite: 3,
            peptideLength: 7,
            DummyMod());

        Assert.That(result, Is.False);
    }

    [Test]
    public void FlankingIonCoverageFilter_Passes_ReturnsFalseIfBothFlanksAreSameIon()
    {
        var filter = new FlankingIonCoverageFilter();
        var ions = new List<MatchedFragmentIon>
        {
            CreateIon(FragmentationTerminus.N, fragmentNumber: 2, residuePosition: 2), // left flank only
            CreateIon(FragmentationTerminus.N, fragmentNumber: 2, residuePosition: 2)  // duplicate left flank
        };

        bool result = filter.Passes(
            DummyPeptide(),
            DummySpectralMatch(),
            0, 0,
            ions,
            peptideOneBasedModSite: 3,
            peptideLength: 7,
            DummyMod());

        Assert.That(result, Is.False);
        Assert.That(IGptmdFilter.GetFilterTypeName(filter), Is.EqualTo("FlankingIonCoverageFilter"));
    }

    [Test]
    public void UniDirectionalIonCoverageFilter_Passes_ReturnsFalseIfNoMatchedIons()
    {
        var filter = new UniDirectionalIonCoverageFilter();
        bool result = filter.Passes(
            DummyPeptide(),
            DummySpectralMatch(),
            0, 0,
            matchedIons: null,
            peptideOneBasedModSite: 3,
            peptideLength: 7,
            DummyMod());

        Assert.That(result, Is.False);
    }

    [Test]
    public void UniDirectionalIonCoverageFilter_Passes_ReturnsFalseIfEmptyMatchedIons()
    {
        var filter = new UniDirectionalIonCoverageFilter();
        bool result = filter.Passes(
            DummyPeptide(),
            DummySpectralMatch(),
            0, 0,
            new List<MatchedFragmentIon>(),
            peptideOneBasedModSite: 3,
            peptideLength: 7,
            DummyMod());

        Assert.That(result, Is.False);
    }

    [Test]
    public void UniDirectionalIonCoverageFilter_Passes_ReturnsTrueIfCoveredFromNTerm()
    {
        var filter = new UniDirectionalIonCoverageFilter();
        var ions = new List<MatchedFragmentIon>
        {
            CreateIon(FragmentationTerminus.N, fragmentNumber: 4, residuePosition: 3) // covers N-term, residuePosition >= site
        };

        bool result = filter.Passes(
            DummyPeptide(),
            DummySpectralMatch(),
            0, 0,
            ions,
            peptideOneBasedModSite: 3,
            peptideLength: 7,
            DummyMod());

        Assert.That(result, Is.True);
    }

    [Test]
    public void UniDirectionalIonCoverageFilter_Passes_ReturnsTrueIfCoveredFromCTerm()
    {
        var filter = new UniDirectionalIonCoverageFilter();
        var ions = new List<MatchedFragmentIon>
        {
            CreateIon(FragmentationTerminus.C, fragmentNumber: 2, residuePosition: 2) // covers C-term, residuePosition < site
        };

        bool result = filter.Passes(
            DummyPeptide(),
            DummySpectralMatch(),
            0, 0,
            ions,
            peptideOneBasedModSite: 3,
            peptideLength: 7,
            DummyMod());

        Assert.That(result, Is.True);
    }

    [Test]
    public void UniDirectionalIonCoverageFilter_Passes_ReturnsTrueIfCoveredFromFivePrime()
    {
        var filter = new UniDirectionalIonCoverageFilter();
        var ions = new List<MatchedFragmentIon>
        {
            CreateIon(FragmentationTerminus.FivePrime, fragmentNumber: 5, residuePosition: 4) // covers N-term, residuePosition >= site
        };

        bool result = filter.Passes(
            DummyPeptide(),
            DummySpectralMatch(),
            0, 0,
            ions,
            peptideOneBasedModSite: 3,
            peptideLength: 7,
            DummyMod());

        Assert.That(result, Is.True);
    }

    [Test]
    public void UniDirectionalIonCoverageFilter_Passes_ReturnsTrueIfCoveredFromThreePrime()
    {
        var filter = new UniDirectionalIonCoverageFilter();
        var ions = new List<MatchedFragmentIon>
        {
            CreateIon(FragmentationTerminus.ThreePrime, fragmentNumber: 2, residuePosition: 2) // covers C-term, residuePosition < site
        };

        bool result = filter.Passes(
            DummyPeptide(),
            DummySpectralMatch(),
            0, 0,
            ions,
            peptideOneBasedModSite: 3,
            peptideLength: 7,
            DummyMod());

        Assert.That(result, Is.True);
    }

    [Test]
    public void UniDirectionalIonCoverageFilter_Passes_ReturnsFalseIfNeitherDirectionCovered()
    {
        var filter = new UniDirectionalIonCoverageFilter();
        var ions = new List<MatchedFragmentIon>
        {
            CreateIon(FragmentationTerminus.N, fragmentNumber: 2, residuePosition: 1), // residuePosition < site, not covered
            CreateIon(FragmentationTerminus.C, fragmentNumber: 4, residuePosition: 5)  // residuePosition >= site, not covered
        };

        bool result = filter.Passes(
            DummyPeptide(),
            DummySpectralMatch(),
            0, 0,
            ions,
            peptideOneBasedModSite: 3,
            peptideLength: 7,
            DummyMod());

        Assert.That(result, Is.False);
    }

    [Test]
    public void DualDirectionalIonCoverageFilter_Passes_NTerminalAcetylation_PassesWithOneDirection()
    {
        var filter = new DualDirectionalIonCoverageFilter();
        // Only N-term ion present
        var ions = new List<MatchedFragmentIon>
        {
            CreateIon(FragmentationTerminus.N, fragmentNumber: 1, residuePosition: 1)
        };
        var nTermAcetyl = GlobalVariables.AllModsKnown.First(p => p.IdWithMotif.Contains("Acetylation") && p.LocationRestriction.Contains("terminal"));

        bool result = filter.Passes(
            DummyPeptide(),
            DummySpectralMatch(),
            0, 0,
            ions,
            peptideOneBasedModSite: 1,
            peptideLength: 7,
            nTermAcetyl);

        Assert.That(result, Is.True);
    }

    [Test]
    public void FlankingIonCoverageFilter_Passes_NTerminalAcetylation_PassesWithOneFlank()
    {
        var filter = new FlankingIonCoverageFilter();

        // Only right flank present (site-1)
        var ions = new List<MatchedFragmentIon>
        {
            CreateIon(FragmentationTerminus.N, fragmentNumber: 1, residuePosition: 1)
        };

        // N-terminal acetylation mod with "terminal" in LocationRestriction
        var nTermAcetyl = GlobalVariables.AllModsKnown.First(p => p.IdWithMotif.Contains("Acetylation") && p.LocationRestriction.Contains("terminal"));

        bool result = filter.Passes(
            DummyPeptide(),
            DummySpectralMatch(),
            0, 0,
            ions,
            peptideOneBasedModSite: 1,
            peptideLength: 7,
            nTermAcetyl);

        Assert.That(result, Is.True);
    }
}


[TestFixture]
public class CleavageSiteFilterTests
{
    private static Modification Substitution(char original, char substituted)
    {
        ModificationMotif.TryGetMotif(original.ToString(), out var motif);
        return new Modification(
            _originalId: $"{original}->{substituted}",
            _modificationType: "1 nucleotide substitution",
            _target: motif,
            _locationRestriction: "Anywhere.",
            _monoisotopicMass: 1);
    }

    /// <summary>
    /// Digests a protein and returns the peptide with the given base sequence, so the peptide
    /// carries the digestion parameters and the parent the filter reads.
    /// </summary>
    private static PeptideWithSetModifications Peptide(string proteinSequence, string peptideSequence,
        string protease = "trypsin")
    {
        var digestionParams = new DigestionParams(protease: protease, maxMissedCleavages: 0, minPeptideLength: 1);
        return new Proteomics.Protein(proteinSequence, "accession")
            .Digest(digestionParams, new List<Modification>(), new List<Modification>())
            .Single(p => p.BaseSequence == peptideSequence);
    }

    /// <summary>
    /// GptmdEngine passes the site in "one is N-terminus" space, so residue one arrives as two.
    /// </summary>
    private static bool Passes(PeptideWithSetModifications peptide, int oneBasedResidueInPeptide, Modification mod)
        => new CleavageSiteFilter().Passes(
            peptide,
            GptmdFilterTests.DummySpectralMatch(),
            newScore: 10,
            originalScore: 1,
            matchedIons: new List<MatchedFragmentIon>(),
            peptideOneBasedModSite: oneBasedResidueInPeptide + 1,
            peptideLength: peptide.BaseSequence.Length,
            mod);

    [Test]
    public void Rejects_SubstitutionThatCreatesACleavageSiteInsideThePeptide()
    {
        // PEPTIDEK | AAAR. T at residue 4 becomes R, which trypsin would cut after.
        var peptide = Peptide("PEPTIDEKAAAR", "PEPTIDEK");
        Assert.That(Passes(peptide, 4, Substitution('T', 'R')), Is.False);
    }

    [Test]
    public void Rejects_SubstitutionThatDestroysTheCleavageSiteEndingThePeptide()
    {
        // The K at residue 8 is why this peptide ends here; without it the peptide runs on.
        var peptide = Peptide("PEPTIDEKAAAR", "PEPTIDEK");
        Assert.That(Passes(peptide, 8, Substitution('K', 'A')), Is.False);
    }

    [Test]
    public void Rejects_SubstitutionThatBlocksCleavageBeforeProline()
    {
        // trypsin|P does not cut before proline, so A at the first residue becoming P removes the
        // cut after the preceding K that produced this peptide.
        var peptide = Peptide("PEPTIDEKAAAR", "AAAR", protease: "trypsin|P");
        Assert.That(Passes(peptide, 1, Substitution('A', 'P')), Is.False);
    }

    [Test]
    public void Accepts_SubstitutionThatLeavesCleavageUnchanged()
    {
        // T at residue 4 becomes S: no cut site is created, destroyed or blocked.
        var peptide = Peptide("PEPTIDEKAAAR", "PEPTIDEK");
        Assert.That(Passes(peptide, 4, Substitution('T', 'S')), Is.True);
    }

    [Test]
    public void Accepts_SubstitutionBetweenTwoCleavableResidues()
    {
        // K to R at the C terminus keeps the site, so the peptide's boundaries do not move.
        var peptide = Peptide("PEPTIDEKAAAR", "PEPTIDEK");
        Assert.That(Passes(peptide, 8, Substitution('K', 'R')), Is.True);
    }

    [Test]
    public void Accepts_SubstitutionWhenProlineRuleDoesNotApply()
    {
        // Plain trypsin cuts after K regardless of a following proline, so the same substitution
        // that fails under trypsin|P is inert here.
        var peptide = Peptide("PEPTIDEKAAAR", "AAAR");
        Assert.That(Passes(peptide, 1, Substitution('A', 'P')), Is.True);
    }

    [Test]
    public void Accepts_ModificationsThatAreNotSubstitutions()
    {
        var peptide = Peptide("PEPTIDEKAAAR", "PEPTIDEK");
        ModificationMotif.TryGetMotif("T", out var motif);
        var phospho = new Modification(_originalId: "Phospho", _modificationType: "Common Biological",
            _target: motif, _locationRestriction: "Anywhere.", _monoisotopicMass: 79.96633);

        Assert.That(Passes(peptide, 4, phospho), Is.True);
    }

    [Test]
    public void Accepts_WhenTheSiteIsNotAResidue()
    {
        // Site one is the peptide N terminus, not a residue, so there is nothing to substitute.
        // GptmdEngine never sends it for a residue modification; accept rather than guess.
        var peptide = Peptide("PEPTIDEKAAAR", "PEPTIDEK");

        Assert.That(
            new CleavageSiteFilter().Passes(
                peptide,
                GptmdFilterTests.DummySpectralMatch(),
                10, 1, new List<MatchedFragmentIon>(),
                peptideOneBasedModSite: 1,
                peptideLength: peptide.BaseSequence.Length,
                Substitution('K', 'A')),
            Is.True);
    }

    [Test]
    public void Accepts_WhenTheCandidateHasNoDigestionAgent()
    {
        // The dummy peptide is constructed without digestion parameters; the filter cannot reason
        // about cleavage and must not reject on that basis.
        Assert.That(
            new CleavageSiteFilter().Passes(
                GptmdFilterTests.DummyPeptide(),
                GptmdFilterTests.DummySpectralMatch(),
                10, 1, new List<MatchedFragmentIon>(),
                peptideOneBasedModSite: 5,
                peptideLength: 7,
                Substitution('T', 'R')),
            Is.True);
    }
}
