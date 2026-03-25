using EngineLayer;
using NUnit.Framework;
using Omics;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.Linq;

namespace Test.CircularSearch
{
    /// <summary>
    /// Tests for the digestion logic of <see cref="CircularProtein"/>.
    ///
    /// Covers canonical rotation, product counts and types, parent references,
    /// and variable modification key positions.
    /// </summary>
    [TestFixture]
    public static class CircularDigestionLogicTests
    {
        // ── Helpers ───────────────────────────────────────────────────────────

        private static DigestionParams MakeParams(
            int maxMissedCleavages,
            int minPeptideLength = 1,
            int maxModsForPeptides = 2) =>
            new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: maxMissedCleavages,
                minPeptideLength: minPeptideLength,
                maxModsForPeptides: maxModsForPeptides);

        private static List<IBioPolymerWithSetMods> Digest(
            CircularProtein protein,
            DigestionParams digestionParams,
            List<Modification> variableMods = null) =>
            protein.Digest(
                    digestionParams,
                    new List<Modification>(),
                    variableMods ?? new List<Modification>())
                .ToList();

        // ── Canonical rotation tests ──────────────────────────────────────────

        /// <summary>
        /// Verifies that "PEPTIDEKAAK" is canonicalized to "AAKPEPTIDEK".
        /// A comes before P lexicographically, so the rotation starting at the
        /// first A is the smallest.
        /// </summary>
        [Test]
        public static void CanonicalRotation_PEPTIDEKAAK_BecomesAAKPEPTIDEK()
        {
            var protein = new CircularProtein("PEPTIDEKAAK", "acc1");
            Assert.That(protein.BaseSequence, Is.EqualTo("AAKPEPTIDEK"));
        }

        /// <summary>
        /// Verifies that a sequence that is already canonical is unchanged.
        /// "AAKPEPTIDEK" starts with A which is smallest — no rotation needed.
        /// </summary>
        [Test]
        public static void CanonicalRotation_AlreadyCanonical_Unchanged()
        {
            var protein = new CircularProtein("AAKPEPTIDEK", "acc2");
            Assert.That(protein.BaseSequence, Is.EqualTo("AAKPEPTIDEK"));
        }

        /// <summary>
        /// Verifies that two inputs that are rotations of each other canonicalize
        /// to the same sequence — they represent the same ring.
        /// </summary>
        [Test]
        public static void CanonicalRotation_TwoRotations_SameCanonical()
        {
            var p1 = new CircularProtein("PEPTIDEKAAK", "acc_a");
            var p2 = new CircularProtein("AAKPEPTIDEK", "acc_b");
            Assert.That(p1.BaseSequence, Is.EqualTo(p2.BaseSequence));
        }

        /// <summary>
        /// Verifies that "FGHIK" is already canonical (F is the smallest first
        /// character among all rotations).
        /// </summary>
        [Test]
        public static void CanonicalRotation_FGHIK_Unchanged()
        {
            var protein = new CircularProtein("FGHIK", "acc3");
            Assert.That(protein.BaseSequence, Is.EqualTo("FGHIK"));
        }

        // ── maxMissedCleavages < numCleavageSites → no circular product ───────

        /// <summary>
        /// Verifies that when maxMissedCleavages (0) is less than numCleavageSites (2),
        /// no CircularPeptideWithSetModifications is produced.
        ///
        /// Ring "AAKPEPTIDEK" has 2 trypsin cleavage sites (after K at position 3
        /// and after K at position 11). maxMissedCleavages=0 < 2 → circular product
        /// is never emitted.
        /// </summary>
        [Test]
        public static void Digestion_MaxMissedCleavagesLessThanSites_NoCircularProduct()
        {
            var protein = new CircularProtein("PEPTIDEKAAK", "acc4");
            Assert.That(protein.BaseSequence, Is.EqualTo("AAKPEPTIDEK"));

            var products = Digest(protein, MakeParams(maxMissedCleavages: 0));

            var circular = products.OfType<CircularPeptideWithSetModifications>().ToList();
            Assert.That(circular.Count, Is.EqualTo(0),
                "No CircularPeptideWithSetModifications when maxMissedCleavages < numCleavageSites.");
        }

        /// <summary>
        /// Same ring with maxMissedCleavages=1 < 2 cleavage sites → still no circular product.
        /// </summary>
        [Test]
        public static void Digestion_MaxMissedCleavagesOneLessThanSites_NoCircularProduct()
        {
            var protein = new CircularProtein("PEPTIDEKAAK", "acc5");
            var products = Digest(protein, MakeParams(maxMissedCleavages: 1));

            var circular = products.OfType<CircularPeptideWithSetModifications>().ToList();
            Assert.That(circular.Count, Is.EqualTo(0),
                "No CircularPeptideWithSetModifications when maxMissedCleavages=1 < numCleavageSites=2.");
        }

        // ── maxMissedCleavages >= numCleavageSites → exactly one circular product ──

        /// <summary>
        /// Verifies that when maxMissedCleavages (2) equals numCleavageSites (2),
        /// exactly one CircularPeptideWithSetModifications of length N is produced.
        ///
        /// Ring "AAKPEPTIDEK" has N=11. The circular product must:
        ///   - Be of type CircularPeptideWithSetModifications
        ///   - Have BaseSequence == "AAKPEPTIDEK" (length 11)
        ///   - Have OneBasedStartResidueInProtein == 1
        ///   - Have OneBasedEndResidueInProtein == 11
        /// </summary>
        [Test]
        public static void Digestion_MaxMissedCleavagesEqualsNumSites_ExactlyOneCircularProduct()
        {
            var protein = new CircularProtein("PEPTIDEKAAK", "acc6");
            Assert.That(protein.BaseSequence, Is.EqualTo("AAKPEPTIDEK"));

            var products = Digest(protein, MakeParams(maxMissedCleavages: 2));

            var circular = products.OfType<CircularPeptideWithSetModifications>().ToList();
            Assert.That(circular.Count, Is.EqualTo(1),
                "Exactly one CircularPeptideWithSetModifications when maxMissedCleavages >= numCleavageSites.");

            var cp = circular[0];
            Assert.That(cp.BaseSequence, Is.EqualTo("AAKPEPTIDEK"));
            Assert.That(cp.BaseSequence.Length, Is.EqualTo(11), "Circular product must be length N.");
            Assert.That(cp.OneBasedStartResidueInProtein, Is.EqualTo(1));
            Assert.That(cp.OneBasedEndResidueInProtein, Is.EqualTo(11));
        }

        /// <summary>
        /// Verifies that maxMissedCleavages > numCleavageSites also produces exactly one
        /// circular product — the count never exceeds 1 regardless of the budget.
        /// </summary>
        [Test]
        public static void Digestion_MaxMissedCleavagesGreaterThanNumSites_StillOneCircularProduct()
        {
            var protein = new CircularProtein("PEPTIDEKAAK", "acc7");
            var products = Digest(protein, MakeParams(maxMissedCleavages: 5));

            var circular = products.OfType<CircularPeptideWithSetModifications>().ToList();
            Assert.That(circular.Count, Is.EqualTo(1),
                "Still exactly one CircularPeptideWithSetModifications even when budget far exceeds sites.");
        }

        // ── Full 5-product verification for PEPTIDEKAAK ───────────────────────

        /// <summary>
        /// Verifies the complete digestion product table for "PEPTIDEKAAK" (canonical
        /// "AAKPEPTIDEK", N=11) with trypsin, maxMissedCleavages=2, minPeptideLength=1.
        ///
        /// Expected 5 products:
        ///   AAK              Linear  [1, 3]
        ///   PEPTIDEK         Linear  [4, 11]
        ///   AAKPEPTIDEK      Linear  [1, 11]   (single-cut ring-opening, +H2O)
        ///   PEPTIDEKAAK      Linear  [4, 14]   (wrapping single-cut product)
        ///   AAKPEPTIDEK      Circular [1, 11]  (intact ring, no termini)
        ///
        /// The two AAKPEPTIDEK products have the same base sequence but different
        /// types and different monoisotopic masses.
        /// </summary>
        [Test]
        public static void Digestion_PEPTIDEKAAK_ProducesFiveExpectedProducts()
        {
            var protein = new CircularProtein("PEPTIDEKAAK", "acc8");
            Assert.That(protein.BaseSequence, Is.EqualTo("AAKPEPTIDEK"));

            var products = Digest(protein, MakeParams(maxMissedCleavages: 2));

            Assert.That(products.Count, Is.EqualTo(5),
                "Expected exactly 5 digestion products for AAKPEPTIDEK with maxMissedCleavages=2.");

            // Exactly one circular product
            var circular = products.OfType<CircularPeptideWithSetModifications>().ToList();
            Assert.That(circular.Count, Is.EqualTo(1));

            // Exactly four linear products (including the ring-opening full-length one)
            var linear = products
                .OfType<PeptideWithSetModifications>()
                .Where(p => p is not CircularPeptideWithSetModifications)
                .ToList();
            Assert.That(linear.Count, Is.EqualTo(4));

            // Verify specific linear sequences are present
            var linearSequences = linear.Select(p => p.BaseSequence).OrderBy(s => s).ToList();
            Assert.That(linearSequences, Contains.Item("AAK"));
            Assert.That(linearSequences, Contains.Item("PEPTIDEK"));
            // Two full-length linear products (one non-wrapping, one wrapping)
            Assert.That(linearSequences.Count(s => s == "AAKPEPTIDEK" || s == "PEPTIDEKAAK"),
                Is.EqualTo(2),
                "Expected two full-length linear products: AAKPEPTIDEK and PEPTIDEKAAK.");

            // Circular product is full ring
            Assert.That(circular[0].BaseSequence, Is.EqualTo("AAKPEPTIDEK"));

            // Circular and linear full-length products share base sequence but differ in mass
            var linearFullLength = linear.Single(p => p.BaseSequence == "AAKPEPTIDEK");
            Assert.That(circular[0].MonoisotopicMass,
                Is.Not.EqualTo(linearFullLength.MonoisotopicMass).Within(0.001),
                "Circular and linear full-length products must have different masses (+H2O difference).");
        }

        // ── All products reference the originating CircularProtein ────────────

        /// <summary>
        /// Verifies that every digestion product — both circular and linear — has
        /// its Parent set to the originating CircularProtein instance.
        /// </summary>
        [Test]
        public static void Digestion_AllProducts_HaveCircularProteinAsParent()
        {
            var protein = new CircularProtein("PEPTIDEKAAK", "acc9");
            var products = Digest(protein, MakeParams(maxMissedCleavages: 2));

            foreach (var product in products)
            {
                Assert.That(product.Parent, Is.InstanceOf<CircularProtein>(),
                    $"Product '{product.BaseSequence}' should have a CircularProtein parent.");
                Assert.That(product.Parent, Is.SameAs(protein),
                    $"Product '{product.BaseSequence}' should reference the exact originating protein.");
            }
        }

        // ── Variable modifications ────────────────────────────────────────────

        /// <summary>
        /// Verifies that when MaxMods=0, no modified forms are emitted — the product
        /// count equals the unmodified count regardless of available variable mods.
        /// </summary>
        [Test]
        public static void Digestion_MaxModsZero_NoModifiedFormsEmitted()
        {
            var protein = new CircularProtein("PEPTIDEKAAK", "acc10");

            // Get a real "Anywhere." modification from GlobalVariables to use as variable mod.
            // Oxidation on M is the canonical example but M is not in this ring.
            // Use any Anywhere modification available.
            var anywheremod = GlobalVariables.AllModsKnown
                .FirstOrDefault(m => m.LocationRestriction == "Anywhere."
                                     && m.MonoisotopicMass.HasValue
                                     && protein.BaseSequence.Contains(m.Target?.ToString() ?? ""));

            if (anywheremod == null)
            {
                Assert.Ignore("No suitable Anywhere modification found in GlobalVariables; skipping.");
                return;
            }

            var paramsWithMods = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 2,
                minPeptideLength: 1,
                maxModsForPeptides: 0);   // MaxMods = 0

            var products = Digest(protein, paramsWithMods, new List<Modification> { anywheremod });
            var paramsNoMods = MakeParams(maxMissedCleavages: 2);
            var unmodifiedProducts = Digest(protein, paramsNoMods);

            Assert.That(products.Count, Is.EqualTo(unmodifiedProducts.Count),
                "With MaxMods=0, no additional modified forms should be emitted.");
        }

        /// <summary>
        /// Verifies that the circular peptide from "AAKPEPTIDEK" (N=11) with a
        /// modification on T (position 7 in canonical sequence) has the correct
        /// AllModsOneIsNterminus key.
        ///
        /// For the full-ring CircularPeptideWithSetModifications (start=1):
        ///   T is at canonical position 7 (1-based).
        ///   localIndex = 7 - 1 = 6  (0-based within peptide)
        ///   key = localIndex + 2 = 8
        /// </summary>
        [Test]
        public static void Digestion_CircularProduct_ModOnT_HasCorrectKeyPosition()
        {
            var protein = new CircularProtein("PEPTIDEKAAK", "acc11");
            Assert.That(protein.BaseSequence, Is.EqualTo("AAKPEPTIDEK"));
            // Confirm T is at position 7 in "AAKPEPTIDEK" (1-based)
            Assert.That(protein.BaseSequence[6], Is.EqualTo('T'),
                "Pre-condition: T must be at 1-based position 7 in canonical sequence.");

            // Find a phosphorylation or any Anywhere mod targeting T
            var tMod = GlobalVariables.AllModsKnown
                .FirstOrDefault(m => m.LocationRestriction == "Anywhere."
                                     && m.Target?.ToString() == "T"
                                     && m.MonoisotopicMass.HasValue);

            if (tMod == null)
            {
                Assert.Ignore("No Anywhere modification targeting T found; skipping.");
                return;
            }

            var digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 2,
                minPeptideLength: 1,
                maxModsForPeptides: 1);

            var products = Digest(protein, digestionParams, new List<Modification> { tMod });

            // Find modified circular products
            var modifiedCircular = products
                .OfType<CircularPeptideWithSetModifications>()
                .Where(p => p.AllModsOneIsNterminus.Count > 0)
                .ToList();

            Assert.That(modifiedCircular.Count, Is.GreaterThan(0),
                "At least one modified CircularPeptideWithSetModifications should be produced.");

            foreach (var cp in modifiedCircular)
            {
                // The circular peptide starts at position 1 in the ring.
                // T is at canonical position 7 → localIndex = 6 → key = 8.
                Assert.That(cp.AllModsOneIsNterminus.ContainsKey(8),
                    $"Modified circular peptide should have mod at key 8 (T at canonical pos 7). " +
                    $"Actual keys: [{string.Join(", ", cp.AllModsOneIsNterminus.Keys)}]");
            }
        }
        // ── Top-down protease ─────────────────────────────────────────────────

        /// <summary>
        /// Verifies that the "top-down" protease produces exactly one
        /// CircularPeptideWithSetModifications and no linear products,
        /// regardless of maxMissedCleavages.
        ///
        /// The top-down protease has no cleavage motif and specificity "none",
        /// meaning it identifies zero cleavage sites in any ring. Since
        /// numCleavageSites = 0, the condition maxMissedCleavages >= 0 is always
        /// satisfied and the circular product is always emitted. No linear
        /// PeptideWithSetModifications is ever produced because there are no
        /// cuts to make.
        ///
        /// This has a direct consequence for the CircularSearchEngine: because
        /// there are no linear ring-opening products, the search scores only
        /// internal fragment ions (via FragmentInternally). No b/y terminal
        /// ions are generated or matched.
        ///
        /// The ring used is "AAKPEPTIDEK" (canonical, 2 trypsin sites) to confirm
        /// that the choice of protease — not the ring sequence — determines
        /// whether cleavage sites are found.
        /// </summary>
        [Test]
        public static void Digestion_TopDownProtease_OnlyCircularProduct_NoLinearProducts()
        {
            // Use a ring that has 2 trypsin sites, confirming that it is the
            // protease choice — not the absence of K/R — that suppresses digestion.
            var protein = new CircularProtein("PEPTIDEKAAK", "acc_topdown");
            Assert.That(protein.BaseSequence, Is.EqualTo("AAKPEPTIDEK"),
                "Pre-condition: canonical sequence must be AAKPEPTIDEK.");

            var topDownParams = new DigestionParams(
                protease: "top-down",
                maxMissedCleavages: 0,
                minPeptideLength: 1);

            var products = Digest(protein, topDownParams);

            var circular = products.OfType<CircularPeptideWithSetModifications>().ToList();
            var linear = products
                .OfType<PeptideWithSetModifications>()
                .Where(p => p is not CircularPeptideWithSetModifications)
                .ToList();

            // Exactly one circular product — the intact ring
            Assert.That(circular.Count, Is.EqualTo(1),
                "top-down: exactly one CircularPeptideWithSetModifications expected.");
            Assert.That(circular[0].BaseSequence, Is.EqualTo("AAKPEPTIDEK"),
                "top-down: circular product must span the full canonical ring.");
            Assert.That(circular[0].BaseSequence.Length, Is.EqualTo(11),
                "top-down: circular product must be length N=11.");

            // No linear products — no cuts were made
            Assert.That(linear.Count, Is.EqualTo(0),
                "top-down: no linear PeptideWithSetModifications expected.");

            // Verify that the circular product's parent is the originating protein
            Assert.That(circular[0].Parent, Is.SameAs(protein),
                "top-down: circular product must reference the originating CircularProtein.");
        }

        /// <summary>
        /// Verifies that FragmentInternally on the top-down circular product
        /// produces internal ions, and that calling Fragment() (terminal ions)
        /// on it would be inappropriate — confirming the engine's branching
        /// logic is correct: isCircular == true → FragmentInternally only.
        ///
        /// For a ring of length 11 with minLength=2, we expect at least one
        /// internal fragment. No b/y terminal ion scoring applies.
        /// </summary>
        [Test]
        public static void Digestion_TopDownProtease_CircularProduct_ProducesInternalFragmentsOnly()
        {
            var protein = new CircularProtein("PEPTIDEKAAK", "acc_topdown_frag");
            Assert.That(protein.BaseSequence, Is.EqualTo("AAKPEPTIDEK"));

            var topDownParams = new DigestionParams(
                protease: "top-down",
                maxMissedCleavages: 0,
                minPeptideLength: 1);

            var products = Digest(protein, topDownParams);

            var circularPeptide = products
                .OfType<CircularPeptideWithSetModifications>()
                .Single();

            // Generate internal fragments (the only ion type used for circular products)
            var internalFragments = new List<Omics.Fragmentation.Product>();
            circularPeptide.FragmentInternally(
                MassSpectrometry.DissociationType.HCD,
                minLengthOfFragments: 2,
                internalFragments);

            Assert.That(internalFragments.Count, Is.GreaterThan(0),
                "top-down circular product must produce at least one internal fragment.");

            // All fragments must be internal (FragmentationTerminus.None)
            // and have no primary ion product type (b/y/a/c/z etc.)
            foreach (var frag in internalFragments)
            {
                Assert.That(frag.Terminus, Is.EqualTo(Omics.Fragmentation.FragmentationTerminus.None),
                    $"Internal fragment at [{frag.FragmentNumber}-{frag.SecondaryFragmentNumber}] " +
                    $"must have Terminus=None, not a terminal ion type.");
            }
        }
    }
}