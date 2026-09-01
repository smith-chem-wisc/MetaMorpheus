using System.Collections.Generic;
using System.Linq;
using EngineLayer;
using EngineLayer.Truncation;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Omics.Fragmentation;

namespace Test
{
    /// <summary>
    /// Tests for the fragmentation-propensity scoring foundation (the borrowable C-Score core): the relative
    /// cleavage-weight table and the fragment-ion -> backbone-bond mapping that the weighted score relies on.
    /// </summary>
    [TestFixture]
    public class FragmentationPropensityTests
    {
        [Test]
        public void CleavageWeight_EnhancedSites()
        {
            Assert.That(FragmentationPropensity.CleavageWeight('A', 'A'), Is.EqualTo(1.0)); // average bond
            Assert.That(FragmentationPropensity.CleavageWeight('D', 'A'), Is.EqualTo(3.0)); // C-terminal to Asp
            Assert.That(FragmentationPropensity.CleavageWeight('E', 'A'), Is.EqualTo(2.0)); // C-terminal to Glu
            Assert.That(FragmentationPropensity.CleavageWeight('A', 'P'), Is.EqualTo(3.0)); // N-terminal to Pro
            Assert.That(FragmentationPropensity.CleavageWeight('D', 'P'), Is.EqualTo(9.0)); // both effects
            Assert.That(FragmentationPropensity.CleavageWeight('E', 'P'), Is.EqualTo(6.0));
        }

        [Test]
        public void TryGetCleavageBond_NSeries_MapsToFlankingResidues()
        {
            // "PEPTIDE" (length 7). b2 reports the bond between residue 2 (idx1) and residue 3 (idx2).
            Assert.That(FragmentationPropensity.TryGetCleavageBond(FragmentationTerminus.N, 2, 7, out int before, out int after), Is.True);
            Assert.That(before, Is.EqualTo(1));
            Assert.That(after, Is.EqualTo(2));

            // b7 spans the whole peptide -> no internal bond.
            Assert.That(FragmentationPropensity.TryGetCleavageBond(FragmentationTerminus.N, 7, 7, out _, out _), Is.False);
        }

        [Test]
        public void TryGetCleavageBond_CSeries_MapsToFlankingResidues()
        {
            // "PEPTIDE" (length 7). y2 reports the bond between residue 5 (idx4) and residue 6 (idx5).
            Assert.That(FragmentationPropensity.TryGetCleavageBond(FragmentationTerminus.C, 2, 7, out int before, out int after), Is.True);
            Assert.That(before, Is.EqualTo(4));
            Assert.That(after, Is.EqualTo(5));

            // y7 spans the whole peptide -> no internal bond.
            Assert.That(FragmentationPropensity.TryGetCleavageBond(FragmentationTerminus.C, 7, 7, out _, out _), Is.False);
        }

        /// <summary>
        /// The claim that licenses comparing an internal-fragment score with a terminal one: with uniform
        /// weights this reproduces MetaMorpheusEngine.CalculatePeptideScore EXACTLY, not nearly. Uniform is
        /// what any non-collisional dissociation gets, so an ETD scan is scored on the standard score.
        ///
        /// Both details that make it exact are exercised here: the whole-length b ion, which reports no
        /// internal bond, still contributes (it takes the average weight rather than being dropped), and the
        /// diagnostic ion contributes nothing, as it does there.
        /// </summary>
        [Test]
        public void Score_OutsideCollisionalDissociation_IsExactlyTheStandardScore()
        {
            const string sequence = "PEPTIDE";
            List<MatchedFragmentIon> ions = MatchedIons(sequence.Length);
            MsDataScan scan = ScanWithTotalIonCurrent(1000.0);

            double standard = MetaMorpheusEngine.CalculatePeptideScore(scan, ions);
            double weighted = FragmentationPropensity.Score(ions, sequence, scan.TotalIonCurrent, DissociationType.ETD);

            Assert.That(FragmentationPropensity.ModelApplies(DissociationType.ETD), Is.False,
                "premise: the collisional model must not claim to describe ETD");
            Assert.That(weighted, Is.EqualTo(standard).Within(1e-12),
                "uniform weights must reproduce the standard matched-ion score exactly");
        }

        /// <summary>
        /// And under HCD the weights actually apply, so the same ions score higher -- otherwise the test
        /// above would pass for the trivial reason that nothing is ever weighted.
        /// </summary>
        [Test]
        public void Score_UnderCollisionalDissociation_WeightsTheEnhancedBonds()
        {
            const string sequence = "PEPTIDE";
            List<MatchedFragmentIon> ions = MatchedIons(sequence.Length);
            MsDataScan scan = ScanWithTotalIonCurrent(1000.0);

            double standard = MetaMorpheusEngine.CalculatePeptideScore(scan, ions);
            double weighted = FragmentationPropensity.Score(ions, sequence, scan.TotalIonCurrent, DissociationType.HCD);

            Assert.That(FragmentationPropensity.ModelApplies(DissociationType.HCD), Is.True);
            Assert.That(weighted, Is.GreaterThan(standard),
                "PEPTIDE carries E-P and D-E bonds, so a weighted score must exceed the uniform one");
        }

        [Test]
        public void ModelApplies_OnlyToCollisionalDissociation()
        {
            foreach (DissociationType collisional in new[]
                     { DissociationType.HCD, DissociationType.CID, DissociationType.LowCID, DissociationType.ISCID })
            {
                Assert.That(FragmentationPropensity.ModelApplies(collisional), Is.True, collisional.ToString());
            }

            // EThcD is excluded on purpose: it yields b/y and c/z together, and one scalar per bond cannot be
            // right for both halves of that ladder.
            foreach (DissociationType other in new[]
                     { DissociationType.ETD, DissociationType.ECD, DissociationType.EThcD, DissociationType.UVPD,
                       DissociationType.NETD, DissociationType.Unknown })
            {
                Assert.That(FragmentationPropensity.ModelApplies(other), Is.False, other.ToString());
            }
        }

        /// <summary>
        /// b1..bL and y1..yL over a peptide of the given length, plus one diagnostic ion. bL and yL report no
        /// internal bond; the diagnostic ion must count nothing on either side of the comparison.
        /// </summary>
        private static List<MatchedFragmentIon> MatchedIons(int length)
        {
            var ions = new List<MatchedFragmentIon>();
            for (int i = 1; i <= length; i++)
            {
                ions.Add(new MatchedFragmentIon(
                    new Product(ProductType.b, FragmentationTerminus.N, 100.0 * i, i, i, 0), 100.0 * i, 50.0, 1));
                ions.Add(new MatchedFragmentIon(
                    new Product(ProductType.y, FragmentationTerminus.C, 100.0 * i, i, i, 0), 100.0 * i, 50.0, 1));
            }
            ions.Add(new MatchedFragmentIon(
                new Product(ProductType.D, FragmentationTerminus.None, 204.0, 1, 1, 0), 204.0, 50.0, 1));
            return ions;
        }

        private static MsDataScan ScanWithTotalIonCurrent(double totalIonCurrent) => new MsDataScan(
            massSpectrum: new MzSpectrum(new[] { 100.0 }, new[] { totalIonCurrent }, false),
            oneBasedScanNumber: 1, msnOrder: 2, isCentroid: true, polarity: Polarity.Positive,
            retentionTime: 1, scanWindowRange: new MzRange(0, 2000), scanFilter: "f",
            mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: totalIonCurrent, injectionTime: 1.0,
            noiseData: null, nativeId: "scan=1");
    }
}
