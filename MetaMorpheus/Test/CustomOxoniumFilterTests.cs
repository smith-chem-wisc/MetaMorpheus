using Chemistry;
using EngineLayer;
using EngineLayer.GlycoSearch;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Test
{
    /// <summary>
    /// Strict custom-oxonium filter (option B): the diagnostic ions a user attaches to a custom
    /// monosaccharide (MonosaccharidesCustom.tsv column 4) feed GlycoPeptides.DiagonsticFilter with
    /// strict bidirectional semantics. All tests that register custom monosaccharides restore global
    /// state in a finally block via Glycan.ResetCustomMonosaccharides().
    /// </summary>
    [TestFixture]
    [NonParallelizable] // mutates process-wide Glycan monosaccharide statics
    public class CustomOxoniumFilterTests
    {
        private const int Ion512 = 51219700; // 512.197 m/z * 1e5
        private const int Ion733 = 73325000; // 733.25 m/z * 1e5, deliberately NOT a built-in oxonium ion

        // Build a GlycanBox with a chosen Kind[] without touching the glycan database, via the
        // internal test-only constructor. A real object, so these tests do not depend on which
        // members DiagonsticFilter happens to read today.
        private static GlycanBox BoxWithKind(byte[] kind)
        {
            return new GlycanBox(kind);
        }

        // Reset to *startup* state, not to empty: GlobalVariables.LoadGlycans() registers whatever the
        // shipped MonosaccharidesCustom.tsv contains, and resetting to empty would wipe that for every
        // test that runs afterwards the moment a real row is added to the shipped file.
        private static void RestoreStartupMonosaccharides()
        {
            Glycan.ResetCustomMonosaccharides();
            // GlobalVariables.CustomMonosaccharidePath, not a hand-built path: the file moved out of
            // Glycan_Mods to the data-directory root, and a stale literal here would silently restore
            // nothing (File.Exists false) instead of the startup registrations.
            string shipped = GlobalVariables.CustomMonosaccharidePath;
            if (File.Exists(shipped))
            {
                GlycanDatabase.LoadCustomMonosaccharides(shipped);
            }
        }

        [Test]
        public void CustomOxoniumIons_DefaultNoCustoms_EmptyAndArrayIsByteIdentical()
        {
            try
            {
                Glycan.ResetCustomMonosaccharides();
                Assert.Multiple(() =>
                {
                    Assert.That(Glycan.HasCustomOxoniumIons, Is.False);
                    Assert.That(Glycan.CustomOxoniumIons, Is.Empty);
                    // No customs => the combined array is the very same built-in array reference.
                    Assert.That(Glycan.AllOxoniumIonsIncludingCustoms, Is.SameAs(Glycan.AllOxoniumIons));
                });
            }
            finally
            {
                RestoreStartupMonosaccharides();
            }
        }

        [Test]
        public void RegisterCustomMonosaccharide_WithDiagnosticIons_AppearsInCustomOxoniumIonsAndCombinedArray()
        {
            try
            {
                Glycan.ResetCustomMonosaccharides();
                Glycan.RegisterCustomMonosaccharide("SugarU", 'U', 17603209, new[] { Ion512 });

                int customIndex = Glycan.NameCharDic["SugarU"].Item2; // Kind[] index of the new mono
                Assert.Multiple(() =>
                {
                    Assert.That(Glycan.HasCustomOxoniumIons, Is.True);
                    Assert.That(Glycan.CustomOxoniumIons.Count, Is.EqualTo(1));
                    Assert.That(Glycan.CustomOxoniumIons[0].MzScaled, Is.EqualTo(Ion512));
                    Assert.That(Glycan.CustomOxoniumIons[0].KindIndex, Is.EqualTo(customIndex));
                    Assert.That(Glycan.AllOxoniumIonsIncludingCustoms.Length,
                        Is.EqualTo(Glycan.AllOxoniumIons.Length + 1));
                    Assert.That(Glycan.AllOxoniumIonsIncludingCustoms[Glycan.AllOxoniumIons.Length],
                        Is.EqualTo(Ion512));
                });
            }
            finally
            {
                RestoreStartupMonosaccharides();
            }
        }

        [Test]
        public void ApplyStrictMonoFilter_AcceptsOnlyWhenSignalMatchesMonoPresence()
        {
            Assert.Multiple(() =>
            {
                Assert.That(GlycoPeptides.ApplyStrictMonoFilter(true, true), Is.True);
                Assert.That(GlycoPeptides.ApplyStrictMonoFilter(false, false), Is.True);
                Assert.That(GlycoPeptides.ApplyStrictMonoFilter(true, false), Is.False);
                Assert.That(GlycoPeptides.ApplyStrictMonoFilter(false, true), Is.False);
            });
        }

        [Test]
        public void CheckOxoniumPresence_UsesRelativeIntensityNotBarePresence()
        {
            // reference = 100, so the 0.02 threshold sits at intensity 2.
            double[] intensities = { 0.0, 12.5, 1.5, 2.0 };
            Assert.Multiple(() =>
            {
                Assert.That(GlycoPeptides.CheckOxoniumPresence(intensities, 0, 100), Is.False, "zero intensity");
                Assert.That(GlycoPeptides.CheckOxoniumPresence(intensities, 1, 100), Is.True, "well above threshold");
                Assert.That(GlycoPeptides.CheckOxoniumPresence(intensities, 2, 100), Is.False,
                    "a noise-level peak is below the threshold and must NOT count as observed");
                Assert.That(GlycoPeptides.CheckOxoniumPresence(intensities, 3, 100), Is.False, "exactly at the threshold is not above it");
                // No 138.055 reference: the ratio is undefined and the ion is treated as absent,
                // matching how the built-in NaN comparisons behave.
                Assert.That(GlycoPeptides.CheckOxoniumPresence(intensities, 1, 0), Is.False, "no HexNAc reference");
            });
        }

        [Test]
        public void DiagonsticFilter_NoCustoms_BuiltInRejectStillApplies()
        {
            // Built-in rule: both NeuAc ratios > 0.02 but glycan has no NeuAc (Kind[2]==0) -> reject.
            // Proves the legacy filter path is untouched when no customs are registered.
            try
            {
                Glycan.ResetCustomMonosaccharides();
                double[] intensities = new double[Glycan.AllOxoniumIons.Length];
                intensities[OxoniumIndex_R138] = 100;
                intensities[OxoniumIndex_NeuAc274] = 100;
                intensities[OxoniumIndex_NeuAc292] = 100;
                var box = BoxWithKind(new byte[Glycan.KindCapacity]); // all zero -> no NeuAc

                Assert.That(GlycoPeptides.DiagonsticFilter(intensities, box), Is.False);
            }
            finally
            {
                RestoreStartupMonosaccharides();
            }
        }

        [Test]
        public void DiagonsticFilter_CustomIonObservedAndMonoPresent_Accepts()
        {
            RunSingleCustom(ionObserved: true, monoPresent: true, expectedAccept: true);
        }

        [Test]
        public void DiagonsticFilter_CustomIonAbsentButMonoPresent_Rejects()
        {
            RunSingleCustom(ionObserved: false, monoPresent: true, expectedAccept: false);
        }

        [Test]
        public void DiagonsticFilter_CustomIonObservedButMonoAbsent_Rejects()
        {
            RunSingleCustom(ionObserved: true, monoPresent: false, expectedAccept: false);
        }

        [Test]
        public void DiagonsticFilter_CustomIonAbsentAndMonoAbsent_Accepts()
        {
            RunSingleCustom(ionObserved: false, monoPresent: false, expectedAccept: true);
        }

        private static void RunSingleCustom(bool ionObserved, bool monoPresent, bool expectedAccept)
        {
            try
            {
                Glycan.ResetCustomMonosaccharides();
                Glycan.RegisterCustomMonosaccharide("SugarU", 'U', 17603209, new[] { Ion512 });
                int customIndex = Glycan.NameCharDic["SugarU"].Item2;

                // 138.055 is the reference the relative-intensity test divides by; the NeuAc slots are
                // left at 0, which the built-in rules accept for a box with no NeuAc (Kind[2] == 0).
                double[] intensities = new double[Glycan.AllOxoniumIonsIncludingCustoms.Length];
                intensities[OxoniumIndex_R138] = 100;
                if (ionObserved)
                {
                    intensities[Glycan.AllOxoniumIons.Length] = 500; // first (only) custom slot, ratio 5.0
                }

                byte[] kind = new byte[Glycan.KindCapacity];
                if (monoPresent)
                {
                    kind[customIndex] = 1;
                }

                Assert.That(GlycoPeptides.DiagonsticFilter(intensities, BoxWithKind(kind)),
                    Is.EqualTo(expectedAccept));
            }
            finally
            {
                RestoreStartupMonosaccharides();
            }
        }

        [Test]
        public void DiagonsticFilter_TwoCustoms_OneMismatch_Rejects()
        {
            try
            {
                Glycan.ResetCustomMonosaccharides();
                Glycan.RegisterCustomMonosaccharide("SugarU", 'U', 17603209, new[] { Ion512 });
                Glycan.RegisterCustomMonosaccharide("SugarV", 'V', 18000000, new[] { Ion733 });
                int idxU = Glycan.NameCharDic["SugarU"].Item2;
                int idxV = Glycan.NameCharDic["SugarV"].Item2;

                // CustomOxoniumIons order is by Kind index: [Ion512@U, Ion733@V] at the two appended slots.
                double[] intensities = new double[Glycan.AllOxoniumIonsIncludingCustoms.Length];
                intensities[OxoniumIndex_R138] = 100;                // relative-intensity reference
                intensities[Glycan.AllOxoniumIons.Length] = 500;     // Ion512 observed (matches U)
                // Ion733 slot left 0 (absent) -> mismatches V which IS present below.

                byte[] kind = new byte[Glycan.KindCapacity];
                kind[idxU] = 1;
                kind[idxV] = 1;

                Assert.That(GlycoPeptides.DiagonsticFilter(intensities, BoxWithKind(kind)), Is.False);
            }
            finally
            {
                RestoreStartupMonosaccharides();
            }
        }

        [Test]
        public void DiagonsticFilter_UndersizedIntensityArray_Throws()
        {
            // A length mismatch is a programming error, not a data condition: the only producer of this
            // array sizes it from the same AllOxoniumIonsIncludingCustoms the filter indexes. Reading the
            // missing slot as "absent" would silently reject every candidate carrying the custom
            // monosaccharide, so it must fail loudly. This is the asymmetric case -- the array is short
            // but the monosaccharide IS present -- which a both-sides-undersized test cannot see.
            try
            {
                Glycan.ResetCustomMonosaccharides();
                Glycan.RegisterCustomMonosaccharide("SugarU", 'U', 17603209, new[] { Ion512 });
                int customIndex = Glycan.NameCharDic["SugarU"].Item2;

                double[] tooShortIntensities = new double[Glycan.AllOxoniumIons.Length]; // no custom slot
                tooShortIntensities[OxoniumIndex_R138] = 100;
                byte[] kind = new byte[Glycan.KindCapacity];
                kind[customIndex] = 1;                                                    // mono IS present

                var ex = Assert.Throws<ArgumentException>(
                    () => GlycoPeptides.DiagonsticFilter(tooShortIntensities, BoxWithKind(kind)));
                Assert.That(ex.Message, Does.Contain("slots but"));
            }
            finally
            {
                RestoreStartupMonosaccharides();
            }
        }

        [Test]
        public void DiagonsticFilter_UndersizedKindArray_Throws()
        {
            try
            {
                Glycan.ResetCustomMonosaccharides();
                Glycan.RegisterCustomMonosaccharide("SugarU", 'U', 17603209, new[] { Ion512 });

                double[] intensities = new double[Glycan.AllOxoniumIonsIncludingCustoms.Length];
                intensities[OxoniumIndex_R138] = 100;
                intensities[Glycan.AllOxoniumIons.Length] = 500;      // ion observed
                // Long enough for the built-in rules (which read Kind[0..2]) but one slot short of the
                // custom monosaccharide's index, so only the custom branch sees the mismatch.
                var boxShortKind = BoxWithKind(new byte[Glycan.KindCapacity - 1]);

                var ex = Assert.Throws<ArgumentException>(
                    () => GlycoPeptides.DiagonsticFilter(intensities, boxShortKind));
                Assert.That(ex.Message, Does.Contain("Kind"));
            }
            finally
            {
                RestoreStartupMonosaccharides();
            }
        }

        [Test]
        public void DiagonsticFilter_NoiseLevelCustomIon_IsNotObservedAndDoesNotReject()
        {
            // The footgun the relative-intensity threshold exists to close: a single low peak at the
            // custom ion's m/z must not reject every candidate that lacks the monosaccharide.
            try
            {
                Glycan.ResetCustomMonosaccharides();
                Glycan.RegisterCustomMonosaccharide("SugarU", 'U', 17603209, new[] { Ion512 });

                double[] intensities = new double[Glycan.AllOxoniumIonsIncludingCustoms.Length];
                intensities[OxoniumIndex_R138] = 1000;
                intensities[Glycan.AllOxoniumIons.Length] = 5;        // 0.5% of the reference: noise

                // Candidate has no SugarU. Under bare presence this would have been rejected.
                Assert.That(GlycoPeptides.DiagonsticFilter(intensities, BoxWithKind(new byte[Glycan.KindCapacity])),
                    Is.True);
            }
            finally
            {
                RestoreStartupMonosaccharides();
            }
        }

        [Test]
        public void DiagonsticFilter_NoHexNAcReference_MonoPresent_StrictRuleSkippedAndCandidateKept()
        {
            // The mirror of the noise footgun, arriving from the other side. With no 138.055 in the
            // spectrum the relative-intensity test has no denominator, so there is no evidence that the
            // custom ion is present or absent. Reporting "not observed" would not be neutral under
            // hasSignal == hasMono -- it would reject every candidate that CARRIES the custom
            // monosaccharide. A low-energy HCD spectrum with strong 204.087 and no 138.055 is ordinary,
            // so the strict rule has to stand down rather than pick a side.
            try
            {
                Glycan.ResetCustomMonosaccharides();
                Glycan.RegisterCustomMonosaccharide("SugarU", 'U', 17603209, new[] { Ion512 });
                int customIndex = Glycan.NameCharDic["SugarU"].Item2;

                double[] intensities = new double[Glycan.AllOxoniumIonsIncludingCustoms.Length];
                // A strong 204.087 with no 138.055 is the spectrum this test is about. Index 204 has to
                // be set explicitly: writing 0 into index 4 of a fresh array asserts nothing, and left
                // 204 at zero too, so GlycoSearchEngine's own "no 204 => not a glycopeptide" gate would
                // have dropped the scan before DiagonsticFilter ever saw it.
                intensities[OxoniumIndex_HexNAc204] = 1000;
                Assert.That(intensities[OxoniumIndex_R138], Is.Zero, "no 138.055 reference at all");
                byte[] kind = new byte[Glycan.KindCapacity];
                kind[customIndex] = 1;                               // candidate DOES carry the sugar

                Assert.That(GlycoPeptides.DiagonsticFilter(intensities, BoxWithKind(kind)), Is.True,
                    "with no 138.055 reference the strict custom rule must not reject a candidate carrying the monosaccharide");
            }
            finally
            {
                RestoreStartupMonosaccharides();
            }
        }

        [Test]
        public void DiagonsticFilter_NoHexNAcReference_IonPeakPresentMonoAbsent_StrictRuleSkippedAndCandidateKept()
        {
            // Same spectrum condition from the other direction: a peak sits at the custom ion's m/z but
            // there is no reference to size it against. Neither verdict is available, so neither
            // population may be rejected.
            try
            {
                Glycan.ResetCustomMonosaccharides();
                Glycan.RegisterCustomMonosaccharide("SugarU", 'U', 17603209, new[] { Ion512 });

                double[] intensities = new double[Glycan.AllOxoniumIonsIncludingCustoms.Length];
                intensities[OxoniumIndex_HexNAc204] = 1000;          // a glycopeptide spectrum, as above
                Assert.That(intensities[OxoniumIndex_R138], Is.Zero, "no 138.055 reference at all");
                intensities[Glycan.AllOxoniumIons.Length] = 500;     // custom ion slot has signal

                Assert.That(GlycoPeptides.DiagonsticFilter(intensities, BoxWithKind(new byte[Glycan.KindCapacity])),
                    Is.True);
            }
            finally
            {
                RestoreStartupMonosaccharides();
            }
        }

        [Test]
        public void DiagonsticFilter_NoHexNAcReference_UndersizedIntensityArray_StillThrows()
        {
            // The bounds guard is a programming-error check and must keep firing even on the spectra
            // where the strict comparison itself stands down.
            try
            {
                Glycan.ResetCustomMonosaccharides();
                Glycan.RegisterCustomMonosaccharide("SugarU", 'U', 17603209, new[] { Ion512 });
                int customIndex = Glycan.NameCharDic["SugarU"].Item2;

                double[] tooShort = new double[Glycan.AllOxoniumIons.Length]; // no slot for the custom ion
                byte[] kind = new byte[Glycan.KindCapacity];
                kind[customIndex] = 1;

                Assert.Throws<ArgumentException>(() => GlycoPeptides.DiagonsticFilter(tooShort, BoxWithKind(kind)));
            }
            finally
            {
                RestoreStartupMonosaccharides();
            }
        }

        [Test]
        public void RegisterCustomMonosaccharide_DiagnosticIonDuplicatingBuiltIn_Throws()
        {
            // 204.08720 is the HexNAc oxonium ion, present on essentially every glycopeptide spectrum.
            // Accepting it as a custom ion would make the strict filter reject every candidate lacking
            // the custom sugar -- a silent, near-total wipeout of results.
            try
            {
                Glycan.ResetCustomMonosaccharides();
                var ex = Assert.Throws<ArgumentException>(
                    () => Glycan.RegisterCustomMonosaccharide("SugarU", 'U', 17603209, new[] { 20408720 }));
                Assert.That(ex.Message, Does.Contain("built-in oxonium ion"));

                // Rejected before any state changed: the sugar is not registered.
                Assert.That(Glycan.NameCharDic.ContainsKey("SugarU"), Is.False);
                Assert.That(Glycan.HasCustomOxoniumIons, Is.False);
            }
            finally
            {
                RestoreStartupMonosaccharides();
            }
        }

        [Test]
        public void RegisterCustomMonosaccharide_DiagnosticIonWithinRoundingOfBuiltIn_Throws()
        {
            // "204.087" is not bit-identical to the built-in 204.08720 but is the same ion.
            try
            {
                Glycan.ResetCustomMonosaccharides();
                Assert.Throws<ArgumentException>(
                    () => Glycan.RegisterCustomMonosaccharide("SugarU", 'U', 17603209, new[] { 20408700 }));
            }
            finally
            {
                RestoreStartupMonosaccharides();
            }
        }

        [Test]
        public void RegisterCustomMonosaccharide_DiagnosticIonInsideProductToleranceOfHighMassBuiltIn_Throws()
        {
            // A fixed 0.01 Da window stops being conservative at the top of AllOxoniumIons. 20 ppm at
            // the built-in 657.23544 is 0.0131 Da, so an ion 0.012 Da away clears a fixed window and
            // still matches the same peak -- the exact collision the check exists to prevent. The
            // window has to widen to the product tolerance where the tolerance is the wider of the two.
            const int highMassBuiltIn = 65723544;  // 657.23544, second-largest built-in oxonium ion
            const int offsetScaled = 1200;         // 0.012 Da: outside the 0.01 Da floor, inside 20 ppm
            try
            {
                Glycan.ResetCustomMonosaccharides();
                var ex = Assert.Throws<ArgumentException>(
                    () => Glycan.RegisterCustomMonosaccharide("SugarU", 'U', 17603209, new[] { highMassBuiltIn + offsetScaled }));
                Assert.That(ex.Message, Does.Contain("657.23544"));
            }
            finally
            {
                RestoreStartupMonosaccharides();
            }
        }

        [Test]
        public void RegisterCustomMonosaccharide_SameOffsetFromALowMassBuiltIn_IsAccepted()
        {
            // The companion to the test above, so widening the window is not mistaken for widening it
            // everywhere. 20 ppm at the built-in 204.08720 is only 0.004 Da, so there the 0.01 Da floor
            // still governs and an ion 0.012 Da away is a genuinely distinct ion that must be allowed.
            const int lowMassBuiltIn = 20408720;   // 204.08720
            const int offsetScaled = 1200;         // the same 0.012 Da as above
            try
            {
                Glycan.ResetCustomMonosaccharides();
                Glycan.RegisterCustomMonosaccharide("SugarU", 'U', 17603209, new[] { lowMassBuiltIn + offsetScaled });
                Assert.That(Glycan.CustomOxoniumIons.Select(i => i.MzScaled),
                    Is.EquivalentTo(new[] { lowMassBuiltIn + offsetScaled }));
            }
            finally
            {
                RestoreStartupMonosaccharides();
            }
        }

        [Test]
        public void DiagnosticIonCollisionWindow_IsNeverNarrowerThanProductToleranceAtAnyBuiltIn()
        {
            // The property the two tests above sample at one mass each, asserted across the whole array:
            // for every built-in oxonium ion the window must be at least 20 ppm of that ion, otherwise
            // an ion that passes the duplicate check can still match the built-in's peak.
            foreach (int builtIn in Glycan.AllOxoniumIons)
            {
                int window = Glycan.DiagnosticIonCollisionWindowScaled(builtIn, builtIn);
                Assert.That(window, Is.GreaterThanOrEqualTo(builtIn * 20 / 1000000),
                    $"window at built-in {(double)builtIn / 1E5:F5} is narrower than a 20 ppm product tolerance");
                Assert.That(window, Is.GreaterThanOrEqualTo(1000),
                    $"window at built-in {(double)builtIn / 1E5:F5} dropped below the 0.01 Da floor");
            }
        }

        [Test]
        public void RegisterCustomMonosaccharide_DiagnosticIonDuplicatingAnotherCustom_Throws()
        {
            // Two monosaccharides claiming one ion cannot both satisfy the strict rule, so no candidate
            // carrying only one of them could ever pass.
            try
            {
                Glycan.ResetCustomMonosaccharides();
                Glycan.RegisterCustomMonosaccharide("SugarU", 'U', 17603209, new[] { Ion512 });
                var ex = Assert.Throws<ArgumentException>(
                    () => Glycan.RegisterCustomMonosaccharide("SugarV", 'V', 18000000, new[] { Ion512 }));
                Assert.That(ex.Message, Does.Contain("SugarU"));
            }
            finally
            {
                RestoreStartupMonosaccharides();
            }
        }

        [Test]
        public void GlycanDiagnosticIons_CustomIon_IsEmittedAsNeutralMassLikeTheBuiltIns()
        {
            // Column 4 holds observed singly-charged m/z, but GlycanDiagnosticIons becomes
            // Modification.DiagnosticIons, which mzLib consumes as NEUTRAL mass -- it assigns the value
            // straight to Product.NeutralMass. That is why every built-in literal in the property is
            // emitted minus a proton. A custom ion added raw was searched a proton high (at m/z 175,
            // ~287x a 20 ppm window), so it never matched and contributed nothing to the diagnostic-ion
            // score, even though the filter path read the same column correctly.
            int protonScaled = Convert.ToInt32(PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass * 1E5);
            try
            {
                Glycan.ResetCustomMonosaccharides();
                Glycan.RegisterCustomMonosaccharide("SugarU", 'U', 17603209, new[] { Ion512 });
                int customIndex = Glycan.NameCharDic["SugarU"].Item2;

                byte[] kind = new byte[Glycan.KindCapacity];
                kind[1] = 1;                 // one HexNAc, so a built-in ion is emitted alongside
                kind[customIndex] = 1;       // and one of the custom sugar
                var ions = new Glycan(kind, "Nxs", GlycanType.N_glycan).GlycanDiagnosticIons;

                Assert.Multiple(() =>
                {
                    Assert.That(ions, Does.Contain(20408720 - protonScaled),
                        "the built-in 204.087 anchors the convention: m/z minus a proton");
                    Assert.That(ions, Does.Contain(Ion512 - protonScaled),
                        "a custom ion must carry the same m/z-to-neutral-mass conversion as the built-ins");
                    Assert.That(ions, Does.Not.Contain(Ion512),
                        "emitting the raw m/z searches the ion a proton high, where it never matches");
                });
            }
            finally
            {
                RestoreStartupMonosaccharides();
            }
        }

        [Test]
        public void CustomOxoniumIons_AreCachedAndRebuiltOnlyOnRegistrationChange()
        {
            // Both are read on the hot path -- CustomOxoniumIons once per candidate, the combined array
            // once per scan -- so repeated reads must not re-allocate or re-sort.
            try
            {
                Glycan.ResetCustomMonosaccharides();
                Glycan.RegisterCustomMonosaccharide("SugarU", 'U', 17603209, new[] { Ion512 });

                var ionsFirst = Glycan.CustomOxoniumIons;
                var arrayFirst = Glycan.AllOxoniumIonsIncludingCustoms;
                Assert.Multiple(() =>
                {
                    Assert.That(Glycan.CustomOxoniumIons, Is.SameAs(ionsFirst));
                    Assert.That(Glycan.AllOxoniumIonsIncludingCustoms, Is.SameAs(arrayFirst));
                });

                // ...and a new registration does invalidate them.
                Glycan.RegisterCustomMonosaccharide("SugarV", 'V', 18000000, new[] { Ion733 });
                Assert.Multiple(() =>
                {
                    Assert.That(Glycan.CustomOxoniumIons, Is.Not.SameAs(ionsFirst));
                    Assert.That(Glycan.CustomOxoniumIons.Count, Is.EqualTo(2));
                    Assert.That(Glycan.AllOxoniumIonsIncludingCustoms.Length,
                        Is.EqualTo(Glycan.AllOxoniumIons.Length + 2));
                });

                // ...and reset returns the built-in array reference itself.
                Glycan.ResetCustomMonosaccharides();
                Assert.That(Glycan.AllOxoniumIonsIncludingCustoms, Is.SameAs(Glycan.AllOxoniumIons));
            }
            finally
            {
                RestoreStartupMonosaccharides();
            }
        }

        [Test]
        public void LoadCustomMonosaccharides_DiagnosticIonsColumn_RoundTripsIntoCustomOxoniumIons()
        {
            string tsv = string.Join("\n", new[]
            {
                "# comment",
                "Name\tSingleCharCode\tMonoisotopicMass\tDiagnosticIonMasses\tDescription",
                "HexA\tU\t176.03209\t512.197\tHexuronic acid"
            });
            string path = Path.GetTempFileName();
            try
            {
                Glycan.ResetCustomMonosaccharides();
                File.WriteAllText(path, tsv);
                GlycanDatabase.LoadCustomMonosaccharides(path);

                int idx = Glycan.NameCharDic["HexA"].Item2;
                Assert.Multiple(() =>
                {
                    Assert.That(Glycan.HasCustomOxoniumIons, Is.True);
                    Assert.That(Glycan.CustomOxoniumIons.Count, Is.EqualTo(1));
                    Assert.That(Glycan.CustomOxoniumIons[0].MzScaled, Is.EqualTo(Ion512));
                    Assert.That(Glycan.CustomOxoniumIons[0].KindIndex, Is.EqualTo(idx));
                });
            }
            finally
            {
                RestoreStartupMonosaccharides();
                File.Delete(path);
            }
        }

        [Test]
        public void LoadCustomMonosaccharides_DiagnosticIonDuplicatingBuiltIn_SkipsThatIonAndKeepsTheRest()
        {
            // A collision at LOAD time must not throw. DiagnosticIonMasses shipped in 1.1.8 without this
            // check, and EnsureCustomMonosaccharideFileExists only writes the template when the file is
            // MISSING, so an upgrading user's file can already hold a colliding ion. Throwing would be
            // fatal rather than corrective: LoadGlycans runs from SetUpGlobalVariables before the GUI's
            // InitializeComponent, with no handler, so the window never opens and the user cannot reach
            // "Open mods/data folder" to fix the file. Drop the ion, keep the sugar, and say so.
            string tsv = string.Join("\n", new[]
            {
                "Name\tSingleCharCode\tMonoisotopicMass\tDiagnosticIonMasses\tDescription",
                "HexA\tU\t176.03209\t204.08720,512.197\tHexuronic acid"
            });
            string path = Path.GetTempFileName();
            try
            {
                Glycan.ResetCustomMonosaccharides();
                File.WriteAllText(path, tsv);

                var warnings = CaptureLoadWarnings(() => GlycanDatabase.LoadCustomMonosaccharides(path));

                Assert.Multiple(() =>
                {
                    Assert.That(Glycan.NameCharDic.ContainsKey("HexA"), Is.True,
                        "the monosaccharide itself must still load -- glycan databases naming it would otherwise fail to parse");
                    Assert.That(Glycan.CustomOxoniumIons.Select(i => i.MzScaled), Is.EquivalentTo(new[] { Ion512 }),
                        "only the colliding ion is dropped; the good one survives");
                    Assert.That(warnings.Any(w => w.Contains("204.08720")), Is.True, "the warning must name the offending ion");
                    Assert.That(warnings.Any(w => w.Contains(Path.GetFileName(path))), Is.True, "and the file");
                    Assert.That(warnings.Any(w => w.Contains("line 2")), Is.True, "and the line");
                });
            }
            finally
            {
                RestoreStartupMonosaccharides();
                File.Delete(path);
            }
        }

        [Test]
        public void LoadCustomMonosaccharides_EveryIonClaimed_LoadsTheMonosaccharideWithNoIons()
        {
            // The degenerate case of the above: nothing survives the screen. The sugar must still
            // register, with no diagnostic ions, rather than the load failing.
            string tsv = string.Join("\n", new[]
            {
                "Name\tSingleCharCode\tMonoisotopicMass\tDiagnosticIonMasses\tDescription",
                "HexA\tU\t176.03209\t204.08720\tHexuronic acid"
            });
            string path = Path.GetTempFileName();
            try
            {
                Glycan.ResetCustomMonosaccharides();
                File.WriteAllText(path, tsv);

                CaptureLoadWarnings(() => GlycanDatabase.LoadCustomMonosaccharides(path));

                Assert.Multiple(() =>
                {
                    Assert.That(Glycan.NameCharDic.ContainsKey("HexA"), Is.True);
                    Assert.That(Glycan.HasCustomOxoniumIons, Is.False, "no ion survived, so there is no custom filter to run");
                });
            }
            finally
            {
                RestoreStartupMonosaccharides();
                File.Delete(path);
            }
        }

        [Test]
        public void LoadCustomMonosaccharides_DiagnosticIonsPresent_WarnsThatTheyActAsStrictGates()
        {
            // The other half of the upgrade problem, and the quieter one. An existing file keeps its
            // column-4 values, OxoniumIonFilt defaults to true, and those ions therefore become
            // accept/reject gates -- while the banner in the shipped template that explains this never
            // reaches the user, precisely because their file already exists and is never rewritten.
            string tsv = string.Join("\n", new[]
            {
                "Name\tSingleCharCode\tMonoisotopicMass\tDiagnosticIonMasses\tDescription",
                "HexA\tU\t176.03209\t512.197\tHexuronic acid"
            });
            string path = Path.GetTempFileName();
            try
            {
                Glycan.ResetCustomMonosaccharides();
                File.WriteAllText(path, tsv);

                var warnings = CaptureLoadWarnings(() => GlycanDatabase.LoadCustomMonosaccharides(path));

                Assert.Multiple(() =>
                {
                    Assert.That(warnings.Any(w => w.Contains("HexA") && w.Contains("OxoniumIonFilt")), Is.True,
                        "a user whose file predates the banner has to learn about the gate somewhere");
                    Assert.That(warnings.Any(w => w.Contains("512.19700")), Is.True);
                });
            }
            finally
            {
                RestoreStartupMonosaccharides();
                File.Delete(path);
            }
        }

        [Test]
        public void LoadCustomMonosaccharides_NoDiagnosticIons_WarnsAboutNothing()
        {
            // The warning above is emitted once per startup, so it must not fire for the common case of
            // a custom monosaccharide with an empty column 4 -- which includes the shipped template.
            string tsv = string.Join("\n", new[]
            {
                "Name\tSingleCharCode\tMonoisotopicMass\tDiagnosticIonMasses\tDescription",
                "HexA\tU\t176.03209\t\tHexuronic acid"
            });
            string path = Path.GetTempFileName();
            try
            {
                Glycan.ResetCustomMonosaccharides();
                File.WriteAllText(path, tsv);

                var warnings = CaptureLoadWarnings(() => GlycanDatabase.LoadCustomMonosaccharides(path));

                Assert.Multiple(() =>
                {
                    Assert.That(Glycan.NameCharDic.ContainsKey("HexA"), Is.True);
                    Assert.That(warnings, Is.Empty);
                });
            }
            finally
            {
                RestoreStartupMonosaccharides();
                File.Delete(path);
            }
        }

        // Non-fatal load problems go to GlobalVariables.ErrorsReadingMods, which is process-wide and
        // drained by the GUI rather than by tests. Snapshot the length and read only the tail, so these
        // tests neither clear messages another fixture put there nor see them as their own.
        private static List<string> CaptureLoadWarnings(Action load)
        {
            GlobalVariables.ErrorsReadingMods ??= new List<string>();
            int before = GlobalVariables.ErrorsReadingMods.Count;
            load();
            return GlobalVariables.ErrorsReadingMods.GetRange(before, GlobalVariables.ErrorsReadingMods.Count - before);
        }

        // Reserved built-in indices used by the legacy filter. Spelled as literals (rather than the
        // OxoniumIonReservedIndices constants) so this built-in-behavior test does not depend on the
        // very constants it is meant to be independent of.
        private const int OxoniumIndex_R138 = 4;
        private const int OxoniumIndex_HexNAc204 = 9;
        private const int OxoniumIndex_NeuAc274 = 10;
        private const int OxoniumIndex_NeuAc292 = 12;
    }
}
