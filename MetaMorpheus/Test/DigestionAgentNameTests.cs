using EngineLayer;
using NUnit.Framework;
using System.Collections.Generic;
using System.Linq;
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

        /// <summary>
        /// A digestion-parameter type that is neither proteolytic nor RNA and names no agent has to come
        /// back null rather than throw. The call site is PostSearchAnalysisTask.QuantificationAnalysis,
        /// which reaches it once per identification while assembling the FlashLFQ input, so a
        /// NullReferenceException here would abort quantification for the whole run rather than degrade
        /// it. Null is the answer FlashLfqEngine already handles: it reads an unknown agent as "do not
        /// restrict match-between-runs", which is the safe direction.
        /// </summary>
        [Test]
        public static void ParametersThatNameNoAgentReportNullRatherThanThrowing()
        {
            IDigestionParams namesNoAgent = new AgentlessDigestionParams();

            Assert.That(namesNoAgent.DigestionAgent, Is.Null, "premise: this is the case the null check exists for");
            Assert.That(() => namesNoAgent.DigestionAgentName(), Throws.Nothing);
            Assert.That(namesNoAgent.DigestionAgentName(), Is.Null);
        }

        /// <summary>
        /// DigestionAgentName reads SpecificProtease unguarded, which is only safe because DigestionParams
        /// populates it unconditionally -- RecordSpecificProtease assigns it from a ProteaseDictionary
        /// lookup that throws rather than returning null, and the property is private-set. mzLib relies on
        /// the same invariant in ToString and Equals. This pins it across every way the parameters are
        /// built, so that if it ever stops holding the failure lands here instead of as a
        /// NullReferenceException inside quantification.
        /// </summary>
        [Test]
        public static void ProteolyticParametersAlwaysCarryASpecificProtease()
        {
            var built = new List<DigestionParams>
            {
                new DigestionParams(),  // the parameterless constructor the toml reader needs
                Params("trypsin", CleavageSpecificity.Full),
                Params("trypsin", CleavageSpecificity.Semi),
                Params("Glu-C", CleavageSpecificity.None),
            };
            built.Add((DigestionParams)built.Last().Clone(FragmentationTerminus.C));

            foreach (var digestionParams in built)
            {
                Assert.That(digestionParams.SpecificProtease, Is.Not.Null);
                Assert.That(digestionParams.DigestionAgentName(), Is.EqualTo(digestionParams.SpecificProtease.Name));
            }
        }

        /// <summary>
        /// Minimal stand-in for a digestion-parameter implementation that names no agent. ParameterTest
        /// carries a similar double for its own unknown-type case, but that one is private to its fixture.
        /// </summary>
        private class AgentlessDigestionParams : IDigestionParams
        {
            public int MaxMissedCleavages { get; set; }
            public int MinLength { get; set; }
            public int MaxLength { get; set; }
            public int MaxModificationIsoforms { get; set; }
            public int MaxMods { get; set; }
            public DigestionAgent DigestionAgent => null;
            public FragmentationTerminus FragmentationTerminus => FragmentationTerminus.Both;
            public CleavageSpecificity SearchModeType => CleavageSpecificity.Full;
            public IDigestionParams Clone(FragmentationTerminus? newTerminus = null) => this;
            public bool Equals(IDigestionParams other) => ReferenceEquals(this, other);
        }
    }
}
