using System.Collections.Generic;
using System.Linq;
using EngineLayer.Truncation;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;

namespace Test
{
    /// <summary>
    /// Tests for the sequence-tag candidate-filtering foundation (doc §11.2.4): de-novo tag extraction from
    /// fragment-mass gaps (<see cref="SequenceTagExtractor"/>) and the protein k-mer index
    /// (<see cref="ProteinTagIndex"/>) that turns a scan's tags into a small candidate-protein set.
    /// </summary>
    [TestFixture]
    public class SequenceTagTests
    {
        // Monoisotopic residue masses used to synthesize fragment ladders (match SequenceTagExtractor's table).
        private const double P = 97.05276, E = 129.04259, T = 101.04768, I = 113.08406, D = 115.02694;

        /// <summary>Builds an ascending fragment-mass ladder whose consecutive gaps spell the given residues.</summary>
        private static List<double> Ladder(double baseMass, params double[] residues)
        {
            var masses = new List<double> { baseMass };
            double m = baseMass;
            foreach (double r in residues) { m += r; masses.Add(m); }
            return masses;
        }

        [Test]
        public void ExtractTags_RecoversContiguousTags_FromCleanLadder()
        {
            // Gaps spell P,E,P,T,I,D,E -> contiguous 4-mers (I normalized to L).
            List<double> ladder = Ladder(100.0, P, E, P, T, I, D, E);
            HashSet<string> tags = SequenceTagExtractor.ExtractTags(ladder, new AbsoluteTolerance(0.02), tagLength: 4);

            Assert.That(tags, Does.Contain("PEPT"));
            Assert.That(tags, Does.Contain("TLDE"));     // T,I,D,E with I->L
            Assert.That(tags, Does.Contain("EPTL"));     // I->L normalization
            Assert.That(tags, Does.Not.Contain("EPTI")); // never emits 'I'
            Assert.That(tags.Count, Is.EqualTo(4));      // exactly the four contiguous windows
        }

        [Test]
        public void ExtractTags_TooFewPeaks_ReturnsEmpty()
        {
            List<double> ladder = Ladder(100.0, P, E); // only 3 peaks, need >= tagLength+1 = 5
            HashSet<string> tags = SequenceTagExtractor.ExtractTags(ladder, new AbsoluteTolerance(0.02), tagLength: 4);
            Assert.That(tags, Is.Empty);
        }

        [Test]
        public void ExtractTags_ShorterTagLength_Works()
        {
            List<double> ladder = Ladder(200.0, P, E, T);
            HashSet<string> tags = SequenceTagExtractor.ExtractTags(ladder, new AbsoluteTolerance(0.02), tagLength: 3);
            Assert.That(tags, Does.Contain("PET"));
            Assert.That(tags.Count, Is.EqualTo(1));
        }

        [Test]
        public void TagIndex_SelectsProteinContainingTag()
        {
            var proteins = new List<Protein>
            {
                new Protein("PEPTIDESEQK", "A"),  // normalized PEPTLDESEQK -> has PEPT and TLDE
                new Protein("GGPEPTGGGGG", "B"),  // has PEPT, not TLDE
                new Protein("WYACFGHMNRK", "C"),  // neither
            };
            var index = new ProteinTagIndex(proteins, tagLength: 4, maxThreads: 2);

            List<int> hitPept = index.GetCandidateProteinIds(new[] { "PEPT" }, minTagHits: 1);
            Assert.That(hitPept.OrderBy(x => x), Is.EqualTo(new[] { 0, 1 }));   // A and B contain PEPT

            List<int> hitC = index.GetCandidateProteinIds(new[] { "WYAC" }, minTagHits: 1);
            Assert.That(hitC, Is.EqualTo(new[] { 2 }));
        }

        [Test]
        public void TagIndex_MinTagHits_RequiresMultipleSupportingTags()
        {
            var proteins = new List<Protein>
            {
                new Protein("PEPTIDESEQK", "A"),  // has both PEPT and TLDE
                new Protein("GGPEPTGGGGG", "B"),  // has PEPT only
            };
            var index = new ProteinTagIndex(proteins, tagLength: 4, maxThreads: 2);

            List<int> both = index.GetCandidateProteinIds(new[] { "PEPT", "TLDE" }, minTagHits: 2);
            Assert.That(both, Is.EqualTo(new[] { 0 }));   // only A is supported by >= 2 tags
        }

        [Test]
        public void TagIndex_IsoleucineLeucineNormalized()
        {
            var proteins = new List<Protein> { new Protein("AAAILKAAA", "X") }; // contains "AILK" -> normalized "ALLK"
            var index = new ProteinTagIndex(proteins, tagLength: 4, maxThreads: 1);

            Assert.That(index.GetCandidateProteinIds(new[] { "ALLK" }, 1), Is.EqualTo(new[] { 0 }));
            Assert.That(index.GetCandidateProteinIds(new[] { "AILK" }, 1), Is.Empty); // raw 'I' never indexed
        }
    }
}
