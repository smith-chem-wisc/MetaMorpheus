using NUnit.Framework;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using Chemistry;
using EngineLayer;
using MassSpectrometry;

namespace Test
{
    [TestFixture]
    public class IGptmdFilterTests
    {
        // Helper to create a dummy PeptideWithSetModifications
        private PeptideWithSetModifications DummyPeptide() => new("PEPTIDE", []);

        // Helper to create a dummy SpectralMatch
        private SpectralMatch DummySpectralMatch() => new PeptideSpectralMatch(DummyPeptide(), 0, 0, 0,
            new Ms2ScanWithSpecificMass(
                new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true,
                    MassSpectrometry.Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN,
                    null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType,
                    0, null),
                (new Proteomics.AminoAcidPolymer.Peptide(DummyPeptide().BaseSequence).MonoisotopicMass + 21.981943)
                .ToMz(1), 1, "filepath", new CommonParameters())
            ,
            new(), new List<MatchedFragmentIon>());

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
                peptideLength: 7);

            Assert.That(result, Is.True);
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
                peptideLength: 7);

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
                peptideLength: 7);

            Assert.That(result, Is.False);
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
                peptideLength: 7);

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
                peptideLength: 7);

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
                peptideLength: 7);

            Assert.That(result, Is.False);
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
                peptideLength: 7);

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
                peptideLength: 7);

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
                peptideLength: 7);

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
                peptideLength: 7);

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
                peptideLength: 7);

            Assert.That(result, Is.False);
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
                peptideLength: 7);

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
                peptideLength: 7);

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
                peptideLength: 7);

            Assert.That(result, Is.False);
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
                peptideLength: 7);

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
                peptideLength: 7);

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
                peptideLength: 7);

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
                peptideLength: 7);

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
                peptideLength: 7);

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
                peptideLength: 7);

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
                peptideLength: 7);

            Assert.That(result, Is.False);
        }
    }
}

