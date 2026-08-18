using EngineLayer.Truncation;
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
    }
}
