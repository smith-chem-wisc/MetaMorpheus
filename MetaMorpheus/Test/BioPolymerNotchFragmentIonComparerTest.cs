using NUnit.Framework;
using Omics.Fragmentation;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.Reflection;
using EngineLayer;
using Omics.Modifications;

namespace Test
{
    [TestFixture]
    public static class BioPolymerNotchFragmentIonComparerTest
    {
        private static BioPolymerNotchFragmentIonComparer comparer;
        private static Protein exampleProtein;
        private static PeptideWithSetModifications examplePwsm;
        private static MatchedFragmentIon exampleIon;
        private static PropertyInfo fullSequenceProperty;
        private static FieldInfo modDictField;

        [SetUp]
        public static void Setup()
        {
            comparer = new BioPolymerNotchFragmentIonComparer();
            exampleProtein = new Protein("PEPTIDEK", "accession");
            examplePwsm = new PeptideWithSetModifications("PEPTIDEK", null, p: exampleProtein);
            exampleIon = new MatchedFragmentIon(new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 0), 100, 100, 1);
            fullSequenceProperty = examplePwsm.GetType().GetProperty("FullSequence", BindingFlags.Instance | BindingFlags.Public | BindingFlags.NonPublic);
            modDictField = examplePwsm.GetType().GetField("_allModsOneIsNterminus", BindingFlags.Instance | BindingFlags.Public | BindingFlags.NonPublic);
        }

        [Test]
        public static void Compare_DifferentNotches()
        {
            var x = (0, _examplePwsm: examplePwsm, new List<MatchedFragmentIon>());
            var y = (2, _examplePwsm: examplePwsm, new List<MatchedFragmentIon>());
            Assert.That(comparer.Compare(x, y), Is.LessThan(0));
        }

        [Test]
        public static void Compare_DifferentFragmentIonCounts()
        {
            var x = (1, _examplePwsm: examplePwsm, new List<MatchedFragmentIon> { exampleIon });
            var y = (1, _examplePwsm: examplePwsm, new List<MatchedFragmentIon>());
            Assert.That(comparer.Compare(x, y), Is.LessThan(0));
        }

        [Test]
        public static void Compare_DifferentNumberOfMods()
        {
            var modifiedPwsm = new PeptideWithSetModifications("PEPTIDEK", null, p: exampleProtein);
            fullSequenceProperty.SetValue(modifiedPwsm, "P[Oxidation]EPT[Reduction]IDEK", null); 
            modDictField.SetValue(modifiedPwsm,
                new Dictionary<int, Modification>
                {
                    {1, new Modification() },
                    {4, new Modification() }
                });
            var x = (0, _examplePwsm: examplePwsm, new List<MatchedFragmentIon>());
            var y = (0, modifiedPwsm, new List<MatchedFragmentIon>());
            Assert.That(comparer.Compare(x, y), Is.LessThan(0));

            // double check that mods are considered before sequence
            fullSequenceProperty.SetValue(modifiedPwsm, "AAAAAAAA", null);
            Assert.That(comparer.Compare(x, y), Is.LessThan(0));
        }

        [Test]
        public static void Compare_DifferentFullSequences()
        {
            var modifiedPwsmFirst = new PeptideWithSetModifications("PEPTIDEK", null, p: exampleProtein);
            var modifiedPwsmSecond = new PeptideWithSetModifications("PEPTIDEK", null, p: exampleProtein);
            fullSequenceProperty.SetValue(modifiedPwsmFirst, "P[Oxidation]EPTIDEK");
            fullSequenceProperty.SetValue(modifiedPwsmSecond, "PEP[Oxidation]TIDEK");

            // Full sequences are compared alphabetically, and '[' comes before 'E'
            var x = (0, modifiedPwsmFirst, new List<MatchedFragmentIon>());
            var y = (0, modifiedPwsmSecond, new List<MatchedFragmentIon>());
            Assert.That(comparer.Compare(x, y), Is.LessThan(0));
        }

        [Test]
        public static void Compare_DifferentAccessions()
        {
            var protein1 = new Protein("PEPTIDEK", "accession1");
            var protein2 = new Protein("PEPTIDEK", "accession2");
            var x = (1, new PeptideWithSetModifications("PEPTIDEK", null, p: protein1), new List<MatchedFragmentIon>());
            var y = (1, new PeptideWithSetModifications("PEPTIDEK", null, p: protein2), new List<MatchedFragmentIon>());
            Assert.That(comparer.Compare(x, y), Is.LessThan(0));
        }

        [Test]
        public static void Compare_DifferentStartResidues()
        {
            var x = (1, new PeptideWithSetModifications("PEPTIDEK", null, p: exampleProtein, oneBasedStartResidueInProtein: 1), new List<MatchedFragmentIon>());
            var y = (1, new PeptideWithSetModifications("PEPTIDEK", null, p: exampleProtein, oneBasedStartResidueInProtein: 5), new List<MatchedFragmentIon>());
            Assert.That(comparer.Compare(x, y), Is.LessThan(0));
        }

        [Test]
        public static void Compare_NullPwsm()
        {
            var x = (0, (PeptideWithSetModifications)null, new List<MatchedFragmentIon>());
            var y = (0, examplePwsm, new List<MatchedFragmentIon>());
            Assert.That(comparer.Compare(x, y), Is.GreaterThan(0));

            x = (0, examplePwsm, new List<MatchedFragmentIon>());
            y = (0, null, new List<MatchedFragmentIon>());
            Assert.That(comparer.Compare(x, y), Is.LessThan(0));

            x = (0, null, new List<MatchedFragmentIon>());
            y = (0, null, new List<MatchedFragmentIon>());
            Assert.That(comparer.Compare(x, y), Is.EqualTo(0));
        }
    }
}

