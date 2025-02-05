using EngineLayer;
using NUnit.Framework;
using Omics;
using Omics.Fragmentation;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;

namespace Test
{
    [TestFixture]
    public static class BioPolymerNotchFragmentIonComparerTest
    {
        private static BioPolymerNotchFragmentIonComparer<(int, IBioPolymerWithSetMods, List<MatchedFragmentIon>)> comparer;

        [SetUp]
        public static void Setup()
        {
            comparer = new BioPolymerNotchFragmentIonComparer<(int, IBioPolymerWithSetMods, List<MatchedFragmentIon>)>();
        }

        private static Protein _exampleProtein = new Protein("PEPTIDEK", "accession");
        private static PeptideWithSetModifications _examplePwsm = new PeptideWithSetModifications("PEPTIDEK", null, p: _exampleProtein);
        private static MatchedFragmentIon _exampleIon = new MatchedFragmentIon(new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 0), 100, 100, 2);

        [Test]
        public static void Compare_DifferentNotches()
        {
            var x = (0, _examplePwsm, new List<MatchedFragmentIon>());
            var y = (2, _examplePwsm, new List<MatchedFragmentIon>());
            Assert.That(comparer.Compare(x, y), Is.GreaterThan(0));
        }

        [Test]
        public static void Compare_DifferentFragmentIonCounts()
        {
            var x = (1, _examplePwsm, new List<MatchedFragmentIon> { _exampleIon });
            var y = (1, _examplePwsm, new List<MatchedFragmentIon>());
            Assert.That(comparer.Compare(x, y), Is.GreaterThan(0));
        }

        [Test]
        public static void Compare_DifferentSequenceLengths()
        {
            var modifiedPwsm = new PeptideWithSetModifications("PEP[Oxidation]T[Reduction]IDEK", null, p: _exampleProtein);
            var x = (0, _examplePwsm, new List<MatchedFragmentIon>());
            var y = (0, modifiedPwsm, new List<MatchedFragmentIon>());
            Assert.That(comparer.Compare(x, y), Is.GreaterThan(0));
        }

        [Test]
        public static void Compare_DifferentFullSequences()
        {
            var modifiedPwsmFirst = new PeptideWithSetModifications("PEP[Oxidation]TIDEK", null, p: _exampleProtein);
            var modifiedPwsmSecond = new PeptideWithSetModifications("P[Oxidation]EPTIDEK", null, p: _exampleProtein);
            var x = (0, modifiedPwsmFirst, new List<MatchedFragmentIon>());
            var y = (0, modifiedPwsmSecond, new List<MatchedFragmentIon>());
            Assert.That(comparer.Compare(x, y), Is.GreaterThan(0));
        }

        //[Test]
        //public static void Compare_DifferentAccessions()
        //{
        //    var x = (1, new MockBioPolymerWithSetMods("SEQ", "ACC1", 1), new List<MatchedFragmentIon>());
        //    var y = (1, new MockBioPolymerWithSetMods("SEQ", "ACC2", 1), new List<MatchedFragmentIon>());
        //    Assert.That(comparer.Compare(x, y), Is.LessThan(0));
        //}

        //[Test]
        //public static void Compare_DifferentStartResidues()
        //{
        //    var x = (1, new MockBioPolymerWithSetMods("SEQ", "ACC", 1), new List<MatchedFragmentIon>());
        //    var y = (1, new MockBioPolymerWithSetMods("SEQ", "ACC", 2), new List<MatchedFragmentIon>());
        //    Assert.That(comparer.Compare(x, y), Is.LessThan(0));
        //}
    }
}

