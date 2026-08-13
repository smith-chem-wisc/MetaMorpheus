using EngineLayer;
using NUnit.Framework;
using Omics.Digestion;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using Transcriptomics.Digestion;

namespace Test
{
    /// <summary>
    /// The digestion agent name is what FlashLFQ compares to decide whether match-between-runs may transfer
    /// an identification between two files (smith-chem-wisc/mzLib#1148). It has to name the enzyme that
    /// digested the sample, since that is what determines which analytes a file can contain.
    /// </summary>
    [TestFixture]
    public static class DigestionAgentNameTests
    {
        private static DigestionParams Params(string protease, CleavageSpecificity specificity) =>
            new DigestionParams(protease, searchModeType: specificity,
                fragmentationTerminus: FragmentationTerminus.N);

        [Test]
        [TestCase("trypsin", CleavageSpecificity.Full)]
        [TestCase("trypsin", CleavageSpecificity.Semi)]
        [TestCase("Glu-C", CleavageSpecificity.Full)]
        [TestCase("Glu-C", CleavageSpecificity.Semi)]
        public static void ReportsTheProteaseForSpecificSearches(string protease, CleavageSpecificity specificity)
        {
            Assert.That(Params(protease, specificity).DigestionAgentName(), Is.EqualTo(protease));
        }

        /// <summary>
        /// The regression. A fully non-specific search replaces the protease with singleN or singleC, so
        /// DigestionAgent names the search strategy rather than the enzyme. Reading it directly made two
        /// files digested with different proteases report the same agent, which silently disabled the
        /// match-between-runs restriction for exactly the datasets it exists to protect.
        /// </summary>
        [Test]
        public static void NonSpecificSearchReportsTheProteaseNotTheSearchStrategy()
        {
            var trypsin = Params("trypsin", CleavageSpecificity.None);
            var gluC = Params("Glu-C", CleavageSpecificity.None);

            Assert.That(trypsin.DigestionAgent.Name, Is.EqualTo("singleN"),
                "premise of this test: the protease really is replaced, so the naive read collapses");
            Assert.That(trypsin.DigestionAgent.Name, Is.EqualTo(gluC.DigestionAgent.Name));

            Assert.That(trypsin.DigestionAgentName(), Is.EqualTo("trypsin"));
            Assert.That(gluC.DigestionAgentName(), Is.EqualTo("Glu-C"));
        }

        /// <summary>
        /// The non-specific search clones its digestion params per terminus, which flips singleN to singleC.
        /// The enzyme has to survive that, or the same file reports two different agents across the passes.
        /// </summary>
        [Test]
        public static void SurvivesThePerTerminusCloneOfANonSpecificSearch()
        {
            var cloned = Params("trypsin", CleavageSpecificity.None).Clone(FragmentationTerminus.C);

            Assert.That(cloned.DigestionAgent.Name, Is.EqualTo("singleC"), "premise: the clone flips the strategy");
            Assert.That(cloned.DigestionAgentName(), Is.EqualTo("trypsin"));
        }

        [Test]
        public static void ReportsTheRnaseForOligos()
        {
            Assert.That(new RnaDigestionParams("top-down").DigestionAgentName(), Is.EqualTo("top-down"));
        }

        [Test]
        public static void UnknownDigestionParamsReportNoAgent()
        {
            Assert.That(((IDigestionParams)null).DigestionAgentName(), Is.Null);
        }
    }
}
