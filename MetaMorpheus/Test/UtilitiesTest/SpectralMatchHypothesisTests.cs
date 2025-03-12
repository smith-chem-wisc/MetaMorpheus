using EngineLayer.SpectrumMatch;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics;
using System;
using System.Collections.Generic;
using EngineLayer;
using Proteomics.ProteolyticDigestion;
using System.Diagnostics.CodeAnalysis;
using Proteomics;

namespace Test.UtilitiesTest
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    internal class SpectralMatchHypothesisTests
    {
        IBioPolymerWithSetMods testPeptide1 = new PeptideWithSetModifications("PEPTIDE", GlobalVariables.AllModsKnownDictionary, p: new Protein("PEPTIDE", "protein"));
        IBioPolymerWithSetMods testPeptide2 = new PeptideWithSetModifications("PE[UniProt:4-carboxyglutamate on E]PTIDE", GlobalVariables.AllModsKnownDictionary, p: new Protein("PEPTIDE", "protein"));

        [Test]
        public void TestEquals_SameObject()
        {
            var matchedIons = new List<MatchedFragmentIon>();
            var tsm = new SpectralMatchHypothesis(1, testPeptide1, matchedIons, 0);

            Assert.That(tsm.Equals(tsm), Is.True);
            Assert.That(tsm.Equals((object)tsm), Is.True);
        }

        [Test]
        public void TestEquals_NullObject()
        {
            var matchedIons = new List<MatchedFragmentIon>();
            var tsm = new SpectralMatchHypothesis(1, testPeptide1, matchedIons, 0);

            Assert.That(tsm.Equals((SpectralMatchHypothesis)null), Is.False);
            Assert.That(tsm.Equals((object)null), Is.False);
        }

        [Test]
        public void TestEquals_DifferentType()
        {
            var matchedIons = new List<MatchedFragmentIon>();
            var tsm = new SpectralMatchHypothesis(1, testPeptide1, matchedIons, 0);

            Assert.That(tsm.Equals(new object()), Is.False);
        }

        [Test]
        public void TestEquals_DifferentValues()
        {
            var matchedIons1 = new List<MatchedFragmentIon>();
            var matchedIons2 = new List<MatchedFragmentIon> { new MatchedFragmentIon(default, 1, 1, 1) };
            var tsm1 = new SpectralMatchHypothesis(1, testPeptide1, matchedIons1, 0);
            var tsm2 = new SpectralMatchHypothesis(2, testPeptide2, matchedIons2, 0);

            Assert.That(tsm1.Equals(tsm2), Is.False);
            Assert.That(tsm1.Equals((object)tsm2), Is.False);
        }

        [Test]
        public void TestEquals_SameValues()
        {
            var matchedIons = new List<MatchedFragmentIon>();
            var tsm1 = new SpectralMatchHypothesis(1, testPeptide1, matchedIons, 0);
            var tsm2 = new SpectralMatchHypothesis(1, testPeptide1, matchedIons, 0);

            Assert.That(tsm1.Equals(tsm2), Is.True);
            Assert.That(tsm1.Equals((object)tsm2), Is.True);
        }

        [Test]
        public void TestGetHashCode()
        {
            var matchedIons = new List<MatchedFragmentIon>();
            var tsm = new SpectralMatchHypothesis(1, testPeptide1, matchedIons, 0);

            Assert.That(tsm.GetHashCode(), Is.EqualTo(HashCode.Combine(testPeptide1, matchedIons, 0, 1)));
        }

        [Test]
        public void TestEquals_ISearchAttempt_SameValues()
        {
            var matchedIons = new List<MatchedFragmentIon>();
            var tsm1 = new SpectralMatchHypothesis(1, testPeptide1, matchedIons, 0);
            var tsm2 = new SpectralMatchHypothesis(1, testPeptide1, matchedIons, 0);

            Assert.That(tsm1.Equals((ISearchAttempt)tsm2), Is.True);
        }

        [Test]
        public void TestEquals_ISearchAttempt_DifferentValues()
        {
            var matchedIons1 = new List<MatchedFragmentIon>();
            var matchedIons2 = new List<MatchedFragmentIon> { new MatchedFragmentIon(default, 1, 1, 1) };
            var tsm1 = new SpectralMatchHypothesis(1, testPeptide1, matchedIons1, 0);

            // different notch
            var tsm2 = new SpectralMatchHypothesis(2, testPeptide2, matchedIons2, 0);
            Assert.That(tsm1.Equals((ISearchAttempt)tsm2), Is.False);

            // different score
            tsm2 = new SpectralMatchHypothesis(1, testPeptide2, matchedIons2, 20);
            Assert.That(tsm1.Equals((ISearchAttempt)tsm2), Is.False);

            // decoy vs target
            var decoyProtein = new Protein("PEPTIDE", "decoy", isDecoy: true);
            var decoyPeptide = new PeptideWithSetModifications("PEPTIDE", null, p: decoyProtein);
            tsm2 = new SpectralMatchHypothesis(1, decoyPeptide, matchedIons2, 0);
            Assert.That(tsm1.Equals((ISearchAttempt)tsm2), Is.False);
        }

        [Test]
        public void TestEquals_ISearchAttempt_Null()
        {
            var matchedIons = new List<MatchedFragmentIon>();
            var tsm = new SpectralMatchHypothesis(1, testPeptide1, matchedIons, 0);

            Assert.That(tsm.Equals((ISearchAttempt)null), Is.False);
        }
    }
}
