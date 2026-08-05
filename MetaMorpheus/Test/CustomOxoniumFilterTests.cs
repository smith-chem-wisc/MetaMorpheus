using EngineLayer;
using EngineLayer.GlycoSearch;
using NUnit.Framework;
using System.Collections.Generic;
using System.IO;
using System.Reflection;
using System.Runtime.CompilerServices;

namespace Test
{
    /// <summary>
    /// Strict custom-oxonium filter (option B): the diagnostic ions a user attaches to a custom
    /// monosaccharide (MonosaccharidesCustom.tsv column 4) feed GlycoPeptides.DiagonsticFilter with
    /// strict bidirectional semantics. All tests that register custom monosaccharides restore global
    /// state in a finally block via Glycan.ResetCustomMonosaccharides().
    /// </summary>
    [TestFixture]
    public class CustomOxoniumFilterTests
    {
        private const int Ion512 = 51219700; // 512.197 m/z * 1e5
        private const int Ion657 = 65723544; // 657.23544 m/z * 1e5

        // Build a GlycanBox with a chosen Kind[] without touching the glycan database. DiagonsticFilter
        // only reads glycanBox.Kind, so an uninitialized instance with the private setter invoked is enough.
        private static GlycanBox BoxWithKind(byte[] kind)
        {
            var box = (GlycanBox)RuntimeHelpers.GetUninitializedObject(typeof(GlycanBox));
            MethodInfo setter = typeof(GlycanBox).GetProperty("Kind").GetSetMethod(nonPublic: true);
            Assert.That(setter, Is.Not.Null, "GlycanBox.Kind private setter not found via reflection.");
            setter.Invoke(box, new object[] { kind });
            return box;
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
                Glycan.ResetCustomMonosaccharides();
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
                Glycan.ResetCustomMonosaccharides();
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
        public void CheckOxoniumPresence_PositiveIntensityIsPresentZeroIsAbsent()
        {
            double[] intensities = { 0.0, 12.5 };
            Assert.Multiple(() =>
            {
                Assert.That(GlycoPeptides.CheckOxoniumPresence(intensities, 0), Is.False);
                Assert.That(GlycoPeptides.CheckOxoniumPresence(intensities, 1), Is.True);
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
                Glycan.ResetCustomMonosaccharides();
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

                // Built-in slots left at 0 so legacy rules do not reject; custom ion appended at the end.
                double[] intensities = new double[Glycan.AllOxoniumIonsIncludingCustoms.Length];
                if (ionObserved)
                {
                    intensities[Glycan.AllOxoniumIons.Length] = 500; // first (only) custom slot
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
                Glycan.ResetCustomMonosaccharides();
            }
        }

        [Test]
        public void DiagonsticFilter_TwoCustoms_OneMismatch_Rejects()
        {
            try
            {
                Glycan.ResetCustomMonosaccharides();
                Glycan.RegisterCustomMonosaccharide("SugarU", 'U', 17603209, new[] { Ion512 });
                Glycan.RegisterCustomMonosaccharide("SugarV", 'V', 18000000, new[] { Ion657 });
                int idxU = Glycan.NameCharDic["SugarU"].Item2;
                int idxV = Glycan.NameCharDic["SugarV"].Item2;

                // CustomOxoniumIons order is by Kind index: [Ion512@U, Ion657@V] at the two appended slots.
                double[] intensities = new double[Glycan.AllOxoniumIonsIncludingCustoms.Length];
                intensities[Glycan.AllOxoniumIons.Length] = 500;     // Ion512 observed (matches U)
                // Ion657 slot left 0 (absent) -> mismatches V which IS present below.

                byte[] kind = new byte[Glycan.KindCapacity];
                kind[idxU] = 1;
                kind[idxV] = 1;

                Assert.That(GlycoPeptides.DiagonsticFilter(intensities, BoxWithKind(kind)), Is.False);
            }
            finally
            {
                Glycan.ResetCustomMonosaccharides();
            }
        }

        [Test]
        public void DiagonsticFilter_UndersizedInputs_GuardsTreatAsAbsentAndDoNotThrow()
        {
            // Defensive guards: a custom ion index past the end of the intensity array, or a custom
            // monosaccharide index past the end of glycanBox.Kind, must be treated as "absent" rather
            // than throwing IndexOutOfRange. Here both are undersized; with neither observed nor present
            // the strict rule accepts.
            try
            {
                Glycan.ResetCustomMonosaccharides();
                Glycan.RegisterCustomMonosaccharide("SugarU", 'U', 17603209, new[] { Ion512 });

                double[] tooShortIntensities = new double[Glycan.AllOxoniumIons.Length]; // no custom slot
                var boxShortKind = BoxWithKind(new byte[2]);                              // shorter than custom index

                Assert.That(GlycoPeptides.DiagonsticFilter(tooShortIntensities, boxShortKind), Is.True);
            }
            finally
            {
                Glycan.ResetCustomMonosaccharides();
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
                Glycan.ResetCustomMonosaccharides();
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
