using EngineLayer;
using EngineLayer.GlycoSearch;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.IO;

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
            string shipped = Path.Combine(GlobalVariables.DataDir, "Glycan_Mods", "MonosaccharidesCustom.tsv");
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
        public void LoadCustomMonosaccharides_DiagnosticIonDuplicatingBuiltIn_ReportsFileAndLine()
        {
            // The collision must reach the user as a MetaMorpheusException naming the file and line,
            // not as a bare ArgumentException from deep in the glycan code.
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

                var ex = Assert.Throws<MetaMorpheusException>(() => GlycanDatabase.LoadCustomMonosaccharides(path));
                Assert.Multiple(() =>
                {
                    Assert.That(ex.Message, Does.Contain(Path.GetFileName(path)));
                    Assert.That(ex.Message, Does.Contain("line 2"));
                    Assert.That(Glycan.HasCustomOxoniumIons, Is.False);
                });
            }
            finally
            {
                RestoreStartupMonosaccharides();
                File.Delete(path);
            }
        }

        // Reserved built-in indices used by the legacy filter. Spelled as literals (rather than the
        // OxoniumIonReservedIndices constants) so this built-in-behavior test does not depend on the
        // very constants it is meant to be independent of.
        private const int OxoniumIndex_R138 = 4;
        private const int OxoniumIndex_NeuAc274 = 10;
        private const int OxoniumIndex_NeuAc292 = 12;
    }
}
