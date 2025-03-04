using NUnit.Framework;
using Omics.Fragmentation;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.Reflection;
using EngineLayer.SpectrumMatch;
using Omics.Modifications;

namespace Test.UtilitiesTest
{
    [TestFixture]
    public static class BioPolymerNotchFragmentIonComparerFromSpectralMatchHypothesisTest
    {
        private static BioPolymerNotchFragmentIonComparer comparer;
        private static Protein exampleProtein;
        private static PeptideWithSetModifications examplePwsm;
        private static MatchedFragmentIon exampleIon;
        private static PropertyInfo fullSequenceProperty;
        private static FieldInfo modDictField;
        private static List<MatchedFragmentIon> emptyList;

        [SetUp]
        public static void Setup()
        {
            comparer = new BioPolymerNotchFragmentIonComparer();
            exampleProtein = new Protein("PEPTIDEK", "accession");
            examplePwsm = new PeptideWithSetModifications("PEPTIDEK", null, p: exampleProtein);
            exampleIon = new MatchedFragmentIon(new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 0), 100, 100, 1);
            fullSequenceProperty = examplePwsm.GetType().GetProperty("FullSequence", BindingFlags.Instance | BindingFlags.Public | BindingFlags.NonPublic);
            modDictField = examplePwsm.GetType().GetField("_allModsOneIsNterminus", BindingFlags.Instance | BindingFlags.Public | BindingFlags.NonPublic);
            emptyList = new();
        }

        [Test]
        public static void Compare_DifferentNotches()
        {
            var x = new SpectralMatchHypothesis(0, examplePwsm, emptyList);
            var y = new SpectralMatchHypothesis(2, examplePwsm, emptyList);
            Assert.That(comparer.Compare(x, y), Is.GreaterThan(0));
        }

        [Test]
        public static void Compare_DifferentFragmentIonCounts()
        {
            var x = new SpectralMatchHypothesis(1, examplePwsm, new List<MatchedFragmentIon> { exampleIon });
            var y = new SpectralMatchHypothesis(1, examplePwsm, emptyList);
            Assert.That(comparer.Compare(x, y), Is.GreaterThan(0));
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
            var x = new SpectralMatchHypothesis(0, examplePwsm, emptyList);
            var y = new SpectralMatchHypothesis(0, modifiedPwsm, emptyList);
            Assert.That(comparer.Compare(x, y), Is.GreaterThan(0));

            // double check that mods are considered before sequence
            fullSequenceProperty.SetValue(modifiedPwsm, "AAAAAAAA", null);
            Assert.That(comparer.Compare(x, y), Is.GreaterThan(0));
        }

        [Test]
        public static void Compare_DifferentFullSequences()
        {
            var modifiedPwsmFirst = new PeptideWithSetModifications("PEPTIDEK", null, p: exampleProtein);
            var modifiedPwsmSecond = new PeptideWithSetModifications("PEPTIDEK", null, p: exampleProtein);
            fullSequenceProperty.SetValue(modifiedPwsmFirst, "P[Oxidation]EPTIDEK");
            fullSequenceProperty.SetValue(modifiedPwsmSecond, "PEP[Oxidation]TIDEK");

            // Full sequences are compared alphabetically, and '[' comes before 'E'
            var x = new SpectralMatchHypothesis(0, modifiedPwsmFirst, emptyList);
            var y = new SpectralMatchHypothesis(0, modifiedPwsmSecond, emptyList);
            Assert.That(comparer.Compare(x, y), Is.GreaterThan(0));
        }

        [Test]
        public static void Compare_DifferentAccessions()
        {
            var protein1 = new Protein("PEPTIDEK", "accession1");
            var protein2 = new Protein("PEPTIDEK", "accession2");
            var x = new SpectralMatchHypothesis(1, new PeptideWithSetModifications("PEPTIDEK", null, p: protein1), emptyList);
            var y = new SpectralMatchHypothesis(1, new PeptideWithSetModifications("PEPTIDEK", null, p: protein2), emptyList);
            Assert.That(comparer.Compare(x, y), Is.GreaterThan(0));
        }

        [Test]
        public static void Compare_DifferentStartResidues()
        {
            var x = new SpectralMatchHypothesis(1, new PeptideWithSetModifications("PEPTIDEK", null, p: exampleProtein, oneBasedStartResidueInProtein: 1), emptyList);
            var y = new SpectralMatchHypothesis(1, new PeptideWithSetModifications("PEPTIDEK", null, p: exampleProtein, oneBasedStartResidueInProtein: 5), emptyList);
            Assert.That(comparer.Compare(x, y), Is.GreaterThan(0));
        }

        [Test]
        public static void Compare_NullPwsm()
        {
            var x = new SpectralMatchHypothesis(0, null, emptyList);
            var y = new SpectralMatchHypothesis(0, examplePwsm, emptyList);
            Assert.That(comparer.Compare(x, y), Is.LessThan(0));

            x = new SpectralMatchHypothesis(0, examplePwsm, emptyList);
            y = new SpectralMatchHypothesis(0, null, emptyList);
            Assert.That(comparer.Compare(x, y), Is.GreaterThan(0));

            x = new SpectralMatchHypothesis(0, null, emptyList);
            y = new SpectralMatchHypothesis(0, null, emptyList);
            Assert.That(comparer.Compare(x, y), Is.EqualTo(0));
        }

        [Test]
        public static void Compare_NullFragmentIons_DoesNotCrash()
        {
            var protein1 = new Protein("PEPTIDEK", "accession1");
            var peptide = new PeptideWithSetModifications("PEPTIDEK", null, p: protein1);
            List<MatchedFragmentIon> nullList = null;

            var x = new SpectralMatchHypothesis(1, peptide, nullList);
            var y = new SpectralMatchHypothesis(1, peptide, emptyList);
            Assert.That(comparer.Compare(x, y), Is.LessThan(0));

            x = new SpectralMatchHypothesis(1, peptide, emptyList);
            y = new SpectralMatchHypothesis(1, peptide, nullList);
            Assert.That(comparer.Compare(x, y), Is.GreaterThan(0));

            x = new SpectralMatchHypothesis(1, peptide, nullList);
            y = new SpectralMatchHypothesis(1, peptide, nullList);
            Assert.That(comparer.Compare(x, y), Is.EqualTo(0));
        }

        [Test]
        public static void Compare_NullParentAccessions()
        {
            var protein1 = new Protein("PEPTIDEK", null);
            var protein2 = new Protein("PEPTIDEK", "accession");

            var nullAccessionPeptide = new PeptideWithSetModifications("PEPTIDEK", null, p: protein1);
            var normalPeptide = new PeptideWithSetModifications("PEPTIDEK", null, p: protein2);

            var x = new SpectralMatchHypothesis(1, nullAccessionPeptide, emptyList);
            var y = new SpectralMatchHypothesis(1, normalPeptide, emptyList);
            Assert.That(comparer.Compare(x, y), Is.LessThan(0));

            x = new SpectralMatchHypothesis(1, normalPeptide, emptyList);
            y = new SpectralMatchHypothesis(1, nullAccessionPeptide, emptyList);
            Assert.That(comparer.Compare(x, y), Is.GreaterThan(0));

            x = new SpectralMatchHypothesis(1, nullAccessionPeptide, emptyList);
            y = new SpectralMatchHypothesis(1, nullAccessionPeptide, emptyList);
            Assert.That(comparer.Compare(x, y), Is.EqualTo(0));
        }

        [Test]
        public static void Compare_NullParent()
        {
            var nullProteinPeptide = new PeptideWithSetModifications("PEPTIDEK", null, p: null);
            var normalPeptide = new PeptideWithSetModifications("PEPTIDEK", null, p: exampleProtein);

            var x = new SpectralMatchHypothesis(1, nullProteinPeptide, emptyList);
            var y = new SpectralMatchHypothesis(1, normalPeptide, emptyList);
            Assert.That(comparer.Compare(x, y), Is.LessThan(0));

            x = new SpectralMatchHypothesis(1, normalPeptide, emptyList);
            y = new SpectralMatchHypothesis(1, nullProteinPeptide, emptyList);
            Assert.That(comparer.Compare(x, y), Is.GreaterThan(0));

            x = new SpectralMatchHypothesis(1, nullProteinPeptide, emptyList);
            y = new SpectralMatchHypothesis(1, nullProteinPeptide, emptyList);
            Assert.That(comparer.Compare(x, y), Is.EqualTo(0));
        }
    }
}
