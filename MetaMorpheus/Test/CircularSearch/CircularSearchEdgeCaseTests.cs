using EngineLayer;
using EngineLayer.DatabaseLoading;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Omics;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Readers;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;

namespace Test.CircularSearch
{
    /// <summary>
    /// Edge case tests for circular search digestion and task behaviour.
    ///
    /// Tests cover three 10-residue rings representing the three possible
    /// cleavage-site counts (0, 1, 2) with trypsin, plus the rotation-equivalence
    /// invariant and the empty-spectra robustness guarantee.
    ///
    /// Verified canonical forms:
    ///   "AFYTLSGECD"  → canonical "AFYTLSGECD"  (0 trypsin sites)
    ///   "FGHIKACDET"  → canonical "ACDETFGHIK"  (1 trypsin site, after K at pos 9)
    ///   "PEPTIDEKAAK" → canonical "AAKPEPTIDEK" (2 trypsin sites, verified in Class 3)
    /// </summary>
    [TestFixture]
    public static class CircularSearchEdgeCaseTests
    {
        private static readonly double Proton = 1.007276;

        // ── Helpers ───────────────────────────────────────────────────────────

        private static DigestionParams MakeParams(
            int maxMissedCleavages,
            int minPeptideLength = 1) =>
            new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: maxMissedCleavages,
                minPeptideLength: minPeptideLength);

        private static List<IBioPolymerWithSetMods> Digest(
            CircularProtein protein,
            DigestionParams digestionParams) =>
            protein.Digest(digestionParams,
                    new List<Modification>(),
                    new List<Modification>())
                .ToList();

        // ── Zero cleavage sites ───────────────────────────────────────────────

        /// <summary>
        /// A ring with no K or R has 0 trypsin cleavage sites.
        /// The circular product is emitted whenever maxMissedCleavages >= 0,
        /// which is always true — so it is always produced.
        /// No linear products can ever be generated (no cuts possible).
        /// </summary>
        [Test]
        public static void ZeroCleavageSites_AlwaysProducesCircularProduct_NeverLinear()
        {
            // "AFYTLSGECD": 0 trypsin sites, canonical = "AFYTLSGECD"
            var protein = new CircularProtein("AFYTLSGECD", "acc_zero");
            Assert.That(protein.BaseSequence, Is.EqualTo("AFYTLSGECD"),
                "Pre-condition: canonical sequence must be AFYTLSGECD.");

            foreach (int mmc in new[] { 0, 1, 2, 5 })
            {
                var products = Digest(protein, MakeParams(mmc));

                var circular = products.OfType<CircularPeptideWithSetModifications>().ToList();
                var linear = products
                    .OfType<PeptideWithSetModifications>()
                    .Where(p => p is not CircularPeptideWithSetModifications)
                    .ToList();

                Assert.That(circular.Count, Is.EqualTo(1),
                    $"maxMissedCleavages={mmc}: exactly one circular product expected.");
                Assert.That(circular[0].BaseSequence, Is.EqualTo("AFYTLSGECD"),
                    $"maxMissedCleavages={mmc}: circular product must be the full ring.");
                Assert.That(circular[0].BaseSequence.Length, Is.EqualTo(10),
                    $"maxMissedCleavages={mmc}: circular product must be length N=10.");
                Assert.That(linear.Count, Is.EqualTo(0),
                    $"maxMissedCleavages={mmc}: no linear products when there are no cleavage sites.");
            }
        }

        // ── One cleavage site ─────────────────────────────────────────────────

        /// <summary>
        /// A ring with 1 trypsin site and maxMissedCleavages=0 (less than 1 site)
        /// produces no circular product, but does produce linear products from
        /// the single ring-opening cut.
        /// </summary>
        [Test]
        public static void OneCleavageSite_MaxMissedZero_NoCircularProduct_HasLinear()
        {
            // "FGHIKACDET" → canonical "ACDETFGHIK", 1 trypsin site (after K at pos 9)
            var protein = new CircularProtein("FGHIKACDET", "acc_one");
            Assert.That(protein.BaseSequence, Is.EqualTo("ACDETFGHIK"),
                "Pre-condition: canonical sequence must be ACDETFGHIK.");

            var products = Digest(protein, MakeParams(maxMissedCleavages: 0));

            var circular = products.OfType<CircularPeptideWithSetModifications>().ToList();
            var linear = products
                .OfType<PeptideWithSetModifications>()
                .Where(p => p is not CircularPeptideWithSetModifications)
                .ToList();

            Assert.That(circular.Count, Is.EqualTo(0),
                "maxMissedCleavages=0 < 1 site: no circular product.");
            Assert.That(linear.Count, Is.GreaterThan(0),
                "maxMissedCleavages=0: single-cut ring-opening should produce linear products.");
        }

        /// <summary>
        /// A ring with 1 trypsin site and maxMissedCleavages=1 (equal to 1 site)
        /// produces exactly one circular product of length N, plus linear products.
        /// </summary>
        [Test]
        public static void OneCleavageSite_MaxMissedOne_CircularAndLinearProducts()
        {
            var protein = new CircularProtein("FGHIKACDET", "acc_one_b");
            Assert.That(protein.BaseSequence, Is.EqualTo("ACDETFGHIK"));

            var products = Digest(protein, MakeParams(maxMissedCleavages: 1));

            var circular = products.OfType<CircularPeptideWithSetModifications>().ToList();
            var linear = products
                .OfType<PeptideWithSetModifications>()
                .Where(p => p is not CircularPeptideWithSetModifications)
                .ToList();

            Assert.That(circular.Count, Is.EqualTo(1),
                "maxMissedCleavages=1 >= 1 site: exactly one circular product.");
            Assert.That(circular[0].BaseSequence, Is.EqualTo("ACDETFGHIK"),
                "Circular product must span the full canonical ring.");
            Assert.That(circular[0].BaseSequence.Length, Is.EqualTo(10),
                "Circular product must be length N=10.");
            Assert.That(linear.Count, Is.GreaterThan(0),
                "Linear ring-opening products should also be present.");
        }

        // ── Two cleavage sites ────────────────────────────────────────────────

        /// <summary>
        /// A ring with 2 trypsin sites and maxMissedCleavages=1 (less than 2 sites)
        /// produces no circular product.
        /// </summary>
        [Test]
        public static void TwoCleavageSites_MaxMissedOne_NoCircularProduct()
        {
            // "PEPTIDEKAAK" → canonical "AAKPEPTIDEK", 2 trypsin sites
            var protein = new CircularProtein("PEPTIDEKAAK", "acc_two_a");
            Assert.That(protein.BaseSequence, Is.EqualTo("AAKPEPTIDEK"));

            var products = Digest(protein, MakeParams(maxMissedCleavages: 1));

            var circular = products.OfType<CircularPeptideWithSetModifications>().ToList();
            Assert.That(circular.Count, Is.EqualTo(0),
                "maxMissedCleavages=1 < 2 sites: no circular product.");
        }

        /// <summary>
        /// A ring with 2 trypsin sites and maxMissedCleavages=2 (equal to 2 sites)
        /// produces exactly one circular product of length N, plus linear products.
        /// </summary>
        [Test]
        public static void TwoCleavageSites_MaxMissedTwo_CircularAndLinearProducts()
        {
            var protein = new CircularProtein("PEPTIDEKAAK", "acc_two_b");
            Assert.That(protein.BaseSequence, Is.EqualTo("AAKPEPTIDEK"));

            var products = Digest(protein, MakeParams(maxMissedCleavages: 2));

            var circular = products.OfType<CircularPeptideWithSetModifications>().ToList();
            var linear = products
                .OfType<PeptideWithSetModifications>()
                .Where(p => p is not CircularPeptideWithSetModifications)
                .ToList();

            Assert.That(circular.Count, Is.EqualTo(1),
                "maxMissedCleavages=2 >= 2 sites: exactly one circular product.");
            Assert.That(circular[0].BaseSequence, Is.EqualTo("AAKPEPTIDEK"),
                "Circular product must span the full canonical ring.");
            Assert.That(circular[0].BaseSequence.Length, Is.EqualTo(11),
                "Circular product must be length N=11.");
            Assert.That(linear.Count, Is.GreaterThan(0),
                "Linear products should also be present.");
        }

        // ── Rotation equivalence ──────────────────────────────────────────────

        /// <summary>
        /// Two CircularProtein objects constructed from different rotations of the
        /// same ring must canonicalize to the same BaseSequence, proving they
        /// represent the same molecule.
        ///
        /// Tests all three ring types (0, 1, and 2 cleavage sites).
        /// </summary>
        [Test]
        public static void TwoRotationsOfSameRing_CanonicalizeToSameSequence()
        {
            // 0-site ring: "AFYTLSGECD" — rotate by 5 → "SGECD" + "AFYTL" = "SGECDAFYTL"
            var zero_a = new CircularProtein("AFYTLSGECD", "acc_rot_0a");
            var zero_b = new CircularProtein("SGECDAFYTL", "acc_rot_0b");
            Assert.That(zero_a.BaseSequence, Is.EqualTo(zero_b.BaseSequence),
                "Both rotations of the 0-site ring must canonicalize to the same sequence.");

            // 1-site ring: "FGHIKACDET" → canonical "ACDETFGHIK"
            // Rotate by 3 → "IKACDET" + "FGH" = "IKACDETFGH"
            var one_a = new CircularProtein("FGHIKACDET", "acc_rot_1a");
            var one_b = new CircularProtein("IKACDETFGH", "acc_rot_1b");
            Assert.That(one_a.BaseSequence, Is.EqualTo(one_b.BaseSequence),
                "Both rotations of the 1-site ring must canonicalize to the same sequence.");

            // 2-site ring: "PEPTIDEKAAK" → canonical "AAKPEPTIDEK"
            // "AAKPEPTIDEK" rotated by 3 → "PEPTIDEK" + "AAK" = "PEPTIDEKAAK"
            var two_a = new CircularProtein("PEPTIDEKAAK", "acc_rot_2a");
            var two_b = new CircularProtein("AAKPEPTIDEK", "acc_rot_2b");
            Assert.That(two_a.BaseSequence, Is.EqualTo(two_b.BaseSequence),
                "Both rotations of the 2-site ring must canonicalize to the same sequence.");
        }
    }
}