using EngineLayer;
using MzLibUtil;
using NUnit.Framework;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace Test.CircularSearch
{
    /// <summary>
    /// Tests for <see cref="CircularSearchParameters"/> default values and the
    /// <see cref="CircularSearchTask.GetMassDiffAcceptor"/> factory method.
    ///
    /// Mirrors the MassDiffAceptorTest / ParseSearchModeTest pattern from
    /// SearchTaskTest, adapted for circular-specific acceptor types.
    /// </summary>
    [TestFixture]
    public static class CircularSearchParametersTests
    {
        // ── Default-value tests ───────────────────────────────────────────────

        [Test]
        public static void DefaultParameters_MassDiffAcceptorType_IsOneMM()
        {
            var p = new CircularSearchParameters();
            Assert.That(p.MassDiffAcceptorType, Is.EqualTo(MassDiffAcceptorType.OneMM));
        }

        [Test]
        public static void DefaultParameters_MinInternalFragmentLength_IsThree()
        {
            var p = new CircularSearchParameters();
            Assert.That(p.MinInternalFragmentLength, Is.EqualTo(3));
        }

        [Test]
        public static void DefaultParameters_SearchTarget_IsTrue()
        {
            var p = new CircularSearchParameters();
            Assert.That(p.SearchTarget, Is.True);
        }

        [Test]
        public static void DefaultParameters_OutputFlags_AreAllTrue()
        {
            var p = new CircularSearchParameters();
            Assert.That(p.WriteMzId, Is.True);
            Assert.That(p.WriteHighQValuePsms, Is.True);
            Assert.That(p.WriteDecoys, Is.True);
            Assert.That(p.WriteContaminants, Is.True);
        }

        [Test]
        public static void DefaultParameters_Parsimony_IsEnabledWithNoOneHitWondersFalse()
        {
            var p = new CircularSearchParameters();
            Assert.That(p.DoParsimony, Is.True);
            Assert.That(p.NoOneHitWonders, Is.False);
        }

        [Test]
        public static void DefaultParameters_MaxFragmentSize_IsPositive()
        {
            var p = new CircularSearchParameters();
            Assert.That(p.MaxFragmentSize, Is.GreaterThan(0));
        }

        [Test]
        public static void DefaultParameters_CustomMdac_IsEmptyString()
        {
            var p = new CircularSearchParameters();
            Assert.That(p.CustomMdac, Is.EqualTo(string.Empty));
        }

        [Test]
        public static void DefaultParameters_ModsToWriteSelection_ContainsExpectedCategories()
        {
            var p = new CircularSearchParameters();
            Assert.That(p.ModsToWriteSelection, Contains.Key("Common Biological"));
            Assert.That(p.ModsToWriteSelection, Contains.Key("Less Common"));
            Assert.That(p.ModsToWriteSelection, Contains.Key("Metal"));
            Assert.That(p.ModsToWriteSelection, Contains.Key("UniProt"));
        }

        // ── GetMassDiffAcceptor factory tests ─────────────────────────────────

        [Test]
        public static void GetMassDiffAcceptor_DefaultOneMM_ReturnsOneMmAcceptor()
        {
            var task = new CircularSearchTask();
            var result = CircularSearchTask.GetMassDiffAcceptor(
                task.CommonParameters.PrecursorMassTolerance,
                task.CircularSearchParameters.MassDiffAcceptorType,
                task.CircularSearchParameters.CustomMdac);

            Assert.That(result.FileNameAddition, Is.EqualTo("1mm"));
        }

        [Test]
        public static void GetMassDiffAcceptor_TwoMM_ReturnsTwoMmAcceptor()
        {
            var task = new CircularSearchTask();
            var result = CircularSearchTask.GetMassDiffAcceptor(
                task.CommonParameters.PrecursorMassTolerance,
                MassDiffAcceptorType.TwoMM,
                task.CircularSearchParameters.CustomMdac);

            Assert.That(result.FileNameAddition, Is.EqualTo("2mm"));
        }

        [Test]
        public static void GetMassDiffAcceptor_ThreeMM_ReturnsThreeMmAcceptor()
        {
            var task = new CircularSearchTask();
            var result = CircularSearchTask.GetMassDiffAcceptor(
                task.CommonParameters.PrecursorMassTolerance,
                MassDiffAcceptorType.ThreeMM,
                task.CircularSearchParameters.CustomMdac);

            Assert.That(result.FileNameAddition, Is.EqualTo("3mm"));
        }

        [Test]
        public static void GetMassDiffAcceptor_Open_ReturnsOpenSearchAcceptor()
        {
            var task = new CircularSearchTask();
            var result = CircularSearchTask.GetMassDiffAcceptor(
                task.CommonParameters.PrecursorMassTolerance,
                MassDiffAcceptorType.Open,
                task.CircularSearchParameters.CustomMdac);

            Assert.That(result.FileNameAddition, Is.EqualTo("OpenSearch"));
        }

        [Test]
        public static void GetMassDiffAcceptor_ModOpen_ReturnsModOpenAcceptor()
        {
            var task = new CircularSearchTask();
            var result = CircularSearchTask.GetMassDiffAcceptor(
                task.CommonParameters.PrecursorMassTolerance,
                MassDiffAcceptorType.ModOpen,
                task.CircularSearchParameters.CustomMdac);

            Assert.That(result.FileNameAddition, Is.EqualTo("-187andUp"));
        }

        [Test]
        public static void GetMassDiffAcceptor_Exact_WithPpmTolerance_ReturnsPpmAcceptor()
        {
            var task = new CircularSearchTask();
            var result = CircularSearchTask.GetMassDiffAcceptor(
                new PpmTolerance(5),
                MassDiffAcceptorType.Exact,
                task.CircularSearchParameters.CustomMdac);

            Assert.That(result.FileNameAddition, Is.EqualTo("5ppmAroundZero"));
        }

        [Test]
        public static void GetMassDiffAcceptor_Exact_WithAbsoluteTolerance_ReturnsDaltonAcceptor()
        {
            var task = new CircularSearchTask();
            var result = CircularSearchTask.GetMassDiffAcceptor(
                new AbsoluteTolerance(0.01),
                MassDiffAcceptorType.Exact,
                task.CircularSearchParameters.CustomMdac);

            Assert.That(result.FileNameAddition, Is.EqualTo("0.01daltonsAroundZero"));
        }

        [Test]
        public static void GetMassDiffAcceptor_CustomPpmAroundZero_ParsesCorrectly()
        {
            var task = new CircularSearchTask();
            var result = CircularSearchTask.GetMassDiffAcceptor(
                task.CommonParameters.PrecursorMassTolerance,
                MassDiffAcceptorType.Custom,
                "custom ppmAroundZero 4");

            Assert.That(result.FileNameAddition, Is.EqualTo("4ppmAroundZero"));
        }

        [Test]
        public static void GetMassDiffAcceptor_CustomDot_ParsesCorrectNotchCount()
        {
            var task = new CircularSearchTask();
            var result = CircularSearchTask.GetMassDiffAcceptor(
                task.CommonParameters.PrecursorMassTolerance,
                MassDiffAcceptorType.Custom,
                "TestCustom dot 5 ppm 0,1.0029,2.0052");

            Assert.That(result.NumNotches, Is.EqualTo(3));
        }

        [Test]
        public static void GetMassDiffAcceptor_CustomInterval_ParsesCorrectly()
        {
            var task = new CircularSearchTask();
            var result = CircularSearchTask.GetMassDiffAcceptor(
                task.CommonParameters.PrecursorMassTolerance,
                MassDiffAcceptorType.Custom,
                "TestCustom interval [0,5];[0,5]");

            Assert.That(result.NumNotches, Is.EqualTo(1));
        }

        [Test]
        public static void GetMassDiffAcceptor_UnknownType_ThrowsMetaMorpheusException()
        {
            var task = new CircularSearchTask();
            Assert.That(() => CircularSearchTask.GetMassDiffAcceptor(
                    task.CommonParameters.PrecursorMassTolerance,
                    (MassDiffAcceptorType)999,
                    task.CircularSearchParameters.CustomMdac),
                Throws.TypeOf<MetaMorpheusException>());
        }

        [Test]
        public static void GetMassDiffAcceptor_CustomInvalidString_ThrowsMetaMorpheusException()
        {
            var task = new CircularSearchTask();
            Assert.That(() => CircularSearchTask.GetMassDiffAcceptor(
                    task.CommonParameters.PrecursorMassTolerance,
                    MassDiffAcceptorType.Custom,
                    "TestCustom bogusMode 5"),
                Throws.TypeOf<MetaMorpheusException>());
        }
    }
}