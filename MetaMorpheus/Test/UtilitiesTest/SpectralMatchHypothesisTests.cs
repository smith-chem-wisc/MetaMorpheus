using EngineLayer.SpectrumMatch;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics;
using System;
using System.Collections.Generic;
using EngineLayer;
using Proteomics.ProteolyticDigestion;
using System.Diagnostics.CodeAnalysis;

namespace Test.UtilitiesTest
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    internal class SpectralMatchHypothesisTests
    {
        IBioPolymerWithSetMods testPeptide1 = new PeptideWithSetModifications("PEPTIDE", GlobalVariables.AllModsKnownDictionary);
        IBioPolymerWithSetMods testPeptide2 = new PeptideWithSetModifications("PE[UniProt:4-carboxyglutamate on E]PTIDE", GlobalVariables.AllModsKnownDictionary);

        [Test]
        public void TestEquals_SameObject()
        {
            var matchedIons = new List<MatchedFragmentIon>();
            var tsm = new SpectralMatchHypothesis(1, testPeptide1, matchedIons);

            Assert.That(tsm.Equals(tsm), Is.True);
        }

        [Test]
        public void TestEquals_NullObject()
        {
            var matchedIons = new List<MatchedFragmentIon>();
            var tsm = new SpectralMatchHypothesis(1, testPeptide1, matchedIons);

            Assert.That(tsm.Equals(null), Is.False);
        }

        [Test]
        public void TestEquals_DifferentType()
        {
            var matchedIons = new List<MatchedFragmentIon>();
            var tsm = new SpectralMatchHypothesis(1, testPeptide1, matchedIons);

            Assert.That(tsm.Equals(new object()), Is.False);
        }

        [Test]
        public void TestEquals_DifferentValues()
        {
            var matchedIons1 = new List<MatchedFragmentIon>();
            var matchedIons2 = new List<MatchedFragmentIon> { new MatchedFragmentIon(default, 1, 1, 1) };
            var tsm1 = new SpectralMatchHypothesis(1, testPeptide1, matchedIons1);
            var tsm2 = new SpectralMatchHypothesis(2, testPeptide2, matchedIons2);

            Assert.That(tsm1.Equals(tsm2), Is.False);
        }

        [Test]
        public void TestEquals_SameValues()
        {
            var matchedIons = new List<MatchedFragmentIon>();
            var tsm1 = new SpectralMatchHypothesis(1, testPeptide1, matchedIons);
            var tsm2 = new SpectralMatchHypothesis(1, testPeptide1, matchedIons);

            Assert.That(tsm1.Equals(tsm2), Is.True);
        }

        [Test]
        public void TestGetHashCode()
        {
            var matchedIons = new List<MatchedFragmentIon>();
            var tsm = new SpectralMatchHypothesis(1, testPeptide1, matchedIons);

            Assert.That(tsm.GetHashCode(), Is.EqualTo(HashCode.Combine(1, testPeptide1, matchedIons)));
        }
    }
}
