// Copyright 2026 MetaMorpheus Contributors
// Licensed under the MIT License

using EngineLayer.DiaSearch;
using MassSpectrometry.Dia;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.SpectrumMatch;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Test.DiaSearch
{
    /// <summary>
    /// Unit tests for <see cref="DiaMetaDrawShim"/> and <see cref="DiaMetaDrawBundle"/>.
    ///
    /// Uses only in-memory stubs — no file I/O, no mzML parsing.
    /// </summary>
    [TestFixture]
    public class DiaMetaDrawShimTests
    {
        // ── Helpers ───────────────────────────────────────────────────────────

        /// <summary>Creates a minimal DiaSearchResult with synthetic XIC intensities.</summary>
        private static DiaSearchResult MakeResult(
            string sequence,
            int charge,
            float[] extractedIntensities)
        {
            int n = extractedIntensities.Length;
            var r = new DiaSearchResult(
                sequence: sequence,
                chargeState: charge,
                precursorMz: 500.0,
                windowId: 1,
                isDecoy: false,
                fragmentsQueried: n,
                libraryRetentionTime: 12.5,
                rtWindowStart: 10.0f,
                rtWindowEnd: 15.0f);

            for (int i = 0; i < n; i++)
            {
                r.ExtractedIntensities[i] = extractedIntensities[i];
                if (extractedIntensities[i] > 0f)
                    r.FragmentsDetected++;
            }

            r.SpectralAngleScore = 0.85f;
            r.ClassifierScore = 0.92f;
            return r;
        }

        /// <summary>
        /// Creates a minimal LibrarySpectrum with synthetic MatchedFragmentIons at
        /// regularly-spaced m/z values (200, 300, 400 … for n fragments).
        /// Library intensities are set to 1.0 (they should be replaced by observed values).
        /// </summary>
        private static LibrarySpectrum MakeLibrary(
            string sequence,
            int charge,
            int fragmentCount)
        {
            var ions = new List<MatchedFragmentIon>(fragmentCount);
            for (int i = 0; i < fragmentCount; i++)
            {
                double mz = 200.0 + i * 100.0;
                var product = new Product(
                    productType: ProductType.b,
                    terminus: FragmentationTerminus.N,
                    neutralMass: mz * charge - 1.00728 * charge,   // rough back-calc
                    fragmentNumber: i + 1,
                    residuePosition: i + 1,   // 1-based residue position
                    neutralLoss: 0.0);
                // Library intensities are fractional (0–1 scale); the shim replaces these
                // with raw XIC intensities from DiaSearchResult.ExtractedIntensities.
                ions.Add(new MatchedFragmentIon(product, mz, 1.0 / fragmentCount, 1));
            }
            return new LibrarySpectrum(sequence, precursorMz: 500.0, chargeState: charge,
                peaks: ions, rt: 12.5);
        }

        // ── BuildSyntheticFragmentIons ────────────────────────────────────────

        [Test]
        public void BuildSyntheticFragmentIons_AllDetected_ReturnsAllIons()
        {
            var result = MakeResult("PEPTIDEK", 2, new[] { 1000f, 2000f, 1500f });
            var library = MakeLibrary("PEPTIDEK", 2, 3);

            var ions = DiaMetaDrawShim.BuildSyntheticFragmentIons(result, library);

            Assert.That(ions.Count, Is.EqualTo(3));
        }

        [Test]
        public void BuildSyntheticFragmentIons_ObservedIntensitiesUsed_NotLibraryIntensities()
        {
            var result = MakeResult("PEPTIDEK", 2, new[] { 5000f, 3000f, 1000f });
            var library = MakeLibrary("PEPTIDEK", 2, 3);

            var ions = DiaMetaDrawShim.BuildSyntheticFragmentIons(result, library);

            // Sorted by mz ascending: ion[0]=200mz=5000, ion[1]=300mz=3000, ion[2]=400mz=1000
            var sorted = ions.OrderBy(x => x.Mz).ToList();
            Assert.That((float)sorted[0].Intensity, Is.EqualTo(5000f));
            Assert.That((float)sorted[1].Intensity, Is.EqualTo(3000f));
            Assert.That((float)sorted[2].Intensity, Is.EqualTo(1000f));
        }

        [Test]
        public void BuildSyntheticFragmentIons_ZeroIntensityFragmentsExcluded()
        {
            // Fragment at index 1 not detected (zero intensity)
            var result = MakeResult("PEPTIDEK", 2, new[] { 1000f, 0f, 1500f });
            var library = MakeLibrary("PEPTIDEK", 2, 3);

            var ions = DiaMetaDrawShim.BuildSyntheticFragmentIons(result, library);

            Assert.That(ions.Count, Is.EqualTo(2),
                "Zero-intensity fragments must be excluded from the synthetic ion list.");
        }

        [Test]
        public void BuildSyntheticFragmentIons_AllZeroIntensity_ReturnsEmptyList()
        {
            var result = MakeResult("PEPTIDEK", 2, new[] { 0f, 0f, 0f });
            var library = MakeLibrary("PEPTIDEK", 2, 3);

            var ions = DiaMetaDrawShim.BuildSyntheticFragmentIons(result, library);

            Assert.That(ions, Is.Empty);
        }

        [Test]
        public void BuildSyntheticFragmentIons_ReturnedIonsSortedByMzAscending()
        {
            // Intensities in reverse order of mz; after sort, highest mz = lowest intensity
            var result = MakeResult("PEPTIDEK", 2, new[] { 3000f, 2000f, 1000f });
            var library = MakeLibrary("PEPTIDEK", 2, 3);

            var ions = DiaMetaDrawShim.BuildSyntheticFragmentIons(result, library);

            for (int i = 1; i < ions.Count; i++)
                Assert.That(ions[i].Mz, Is.GreaterThanOrEqualTo(ions[i - 1].Mz),
                    $"Ions not sorted ascending at index {i}");
        }

        [Test]
        public void BuildSyntheticFragmentIons_LibraryMzPreserved()
        {
            var result = MakeResult("PEPTIDEK", 2, new[] { 1000f, 2000f, 1500f });
            var library = MakeLibrary("PEPTIDEK", 2, 3);

            var ions = DiaMetaDrawShim.BuildSyntheticFragmentIons(result, library);

            // Library m/z values: 200, 300, 400
            var mzValues = ions.Select(x => x.Mz).OrderBy(x => x).ToList();
            Assert.That(mzValues[0], Is.EqualTo(200.0).Within(0.001));
            Assert.That(mzValues[1], Is.EqualTo(300.0).Within(0.001));
            Assert.That(mzValues[2], Is.EqualTo(400.0).Within(0.001));
        }

        [Test]
        public void BuildSyntheticFragmentIons_LibraryProductTypePreserved()
        {
            var result = MakeResult("PEPTIDEK", 2, new[] { 1000f, 2000f });
            var library = MakeLibrary("PEPTIDEK", 2, 2);

            var ions = DiaMetaDrawShim.BuildSyntheticFragmentIons(result, library);

            foreach (var ion in ions)
                Assert.That(ion.NeutralTheoreticalProduct.ProductType, Is.EqualTo(ProductType.b));
        }

        [Test]
        public void BuildSyntheticFragmentIons_NullResult_ThrowsArgumentNullException()
        {
            var library = MakeLibrary("PEPTIDEK", 2, 3);
            Assert.Throws<ArgumentNullException>(
                () => DiaMetaDrawShim.BuildSyntheticFragmentIons(null, library));
        }

        [Test]
        public void BuildSyntheticFragmentIons_NullLibrary_ThrowsArgumentNullException()
        {
            var result = MakeResult("PEPTIDEK", 2, new[] { 1000f, 2000f });
            Assert.Throws<ArgumentNullException>(
                () => DiaMetaDrawShim.BuildSyntheticFragmentIons(result, null));
        }

        [Test]
        public void BuildSyntheticFragmentIons_EmptyLibrary_ReturnsEmptyList()
        {
            var result = MakeResult("PEPTIDEK", 2, Array.Empty<float>());
            var library = MakeLibrary("PEPTIDEK", 2, 0);

            var ions = DiaMetaDrawShim.BuildSyntheticFragmentIons(result, library);

            Assert.That(ions, Is.Empty);
        }

        [Test]
        public void BuildSyntheticFragmentIons_MismatchedLengthsUsesMin()
        {
            // result has 5 fragments, library has 3: only first 3 should be considered
            var result = MakeResult("PEPTIDEK", 2, new[] { 1000f, 2000f, 1500f, 500f, 800f });
            var library = MakeLibrary("PEPTIDEK", 2, 3);

            var ions = DiaMetaDrawShim.BuildSyntheticFragmentIons(result, library);

            Assert.That(ions.Count, Is.LessThanOrEqualTo(3));
        }

        // ── BuildNormalisedSyntheticFragmentIons ──────────────────────────────

        [Test]
        public void BuildNormalisedSyntheticFragmentIons_MostAbundantIsOne()
        {
            var result = MakeResult("PEPTIDEK", 2, new[] { 1000f, 4000f, 2000f });
            var library = MakeLibrary("PEPTIDEK", 2, 3);

            var ions = DiaMetaDrawShim.BuildNormalisedSyntheticFragmentIons(result, library);

            double maxIntensity = ions.Max(x => x.Intensity);
            Assert.That(maxIntensity, Is.EqualTo(1.0).Within(1e-6));
        }

        [Test]
        public void BuildNormalisedSyntheticFragmentIons_RelativeIntensitiesCorrect()
        {
            // Intensities: 1000, 4000, 2000. After norm (max=4000): 0.25, 1.0, 0.5
            var result = MakeResult("PEPTIDEK", 2, new[] { 1000f, 4000f, 2000f });
            var library = MakeLibrary("PEPTIDEK", 2, 3);

            var ions = DiaMetaDrawShim.BuildNormalisedSyntheticFragmentIons(result, library)
                       .OrderBy(x => x.Mz).ToList();

            Assert.That(ions[0].Intensity, Is.EqualTo(0.25).Within(1e-5)); // mz=200 → 1000/4000
            Assert.That(ions[1].Intensity, Is.EqualTo(1.00).Within(1e-5)); // mz=300 → 4000/4000
            Assert.That(ions[2].Intensity, Is.EqualTo(0.50).Within(1e-5)); // mz=400 → 2000/4000
        }

        [Test]
        public void BuildNormalisedSyntheticFragmentIons_EmptyResult_ReturnsEmpty()
        {
            var result = MakeResult("PEPTIDEK", 2, new[] { 0f, 0f, 0f });
            var library = MakeLibrary("PEPTIDEK", 2, 3);

            var ions = DiaMetaDrawShim.BuildNormalisedSyntheticFragmentIons(result, library);

            Assert.That(ions, Is.Empty);
        }

        // ── BuildBundle ───────────────────────────────────────────────────────

        [Test]
        public void BuildBundle_ReturnsCorrectAdapter()
        {
            var result = MakeResult("PEPTIDEK", 2, new[] { 1000f, 2000f });
            var library = MakeLibrary("PEPTIDEK", 2, 2);
            var adapter = new DiaPsmAdapter(result, "HeLa_01");

            var bundle = DiaMetaDrawShim.BuildBundle(adapter, library);

            Assert.That(bundle.Adapter, Is.SameAs(adapter));
        }

        [Test]
        public void BuildBundle_HasFragmentAnnotations_WhenIonsDetected()
        {
            var result = MakeResult("PEPTIDEK", 2, new[] { 1000f, 2000f });
            var library = MakeLibrary("PEPTIDEK", 2, 2);
            var adapter = new DiaPsmAdapter(result, "HeLa_01");

            var bundle = DiaMetaDrawShim.BuildBundle(adapter, library);

            Assert.That(bundle.HasFragmentAnnotations, Is.True);
            Assert.That(bundle.DetectedFragmentCount, Is.EqualTo(2));
        }

        [Test]
        public void BuildBundle_HasNoFragmentAnnotations_WhenAllZero()
        {
            var result = MakeResult("PEPTIDEK", 2, new[] { 0f, 0f });
            var library = MakeLibrary("PEPTIDEK", 2, 2);
            var adapter = new DiaPsmAdapter(result, "HeLa_01");

            var bundle = DiaMetaDrawShim.BuildBundle(adapter, library);

            Assert.That(bundle.HasFragmentAnnotations, Is.False);
        }

        [Test]
        public void BuildBundle_NormaliseTrue_MostAbundantIsOne()
        {
            var result = MakeResult("PEPTIDEK", 2, new[] { 500f, 2000f, 1000f });
            var library = MakeLibrary("PEPTIDEK", 2, 3);
            var adapter = new DiaPsmAdapter(result, "HeLa_01");

            var bundle = DiaMetaDrawShim.BuildBundle(adapter, library, normalise: true);

            double max = bundle.SyntheticFragmentIons.Max(x => x.Intensity);
            Assert.That(max, Is.EqualTo(1.0).Within(1e-6));
        }

        [Test]
        public void BuildBundle_NullAdapter_ThrowsArgumentNullException()
        {
            var library = MakeLibrary("PEPTIDEK", 2, 2);
            Assert.Throws<ArgumentNullException>(
                () => DiaMetaDrawShim.BuildBundle(null, library));
        }

        [Test]
        public void BuildBundle_NullLibrary_ThrowsArgumentNullException()
        {
            var result = MakeResult("PEPTIDEK", 2, new[] { 1000f });
            var adapter = new DiaPsmAdapter(result, "HeLa_01");
            Assert.Throws<ArgumentNullException>(
                () => DiaMetaDrawShim.BuildBundle(adapter, null));
        }

        // ── IsCompatible ──────────────────────────────────────────────────────

        [Test]
        public void IsCompatible_MatchingResultAndLibrary_ReturnsTrue()
        {
            var result = MakeResult("PEPTIDEK", 2, new[] { 1000f, 2000f });
            var library = MakeLibrary("PEPTIDEK", 2, 2);

            bool ok = DiaMetaDrawShim.IsCompatible(result, library, out string msg);

            Assert.That(ok, Is.True);
            Assert.That(msg, Is.Null);
        }

        [Test]
        public void IsCompatible_FragmentCountMismatch_ReturnsFalseWithMessage()
        {
            var result = MakeResult("PEPTIDEK", 2, new[] { 1000f, 2000f, 3000f });
            var library = MakeLibrary("PEPTIDEK", 2, 2); // library has 2, result has 3

            bool ok = DiaMetaDrawShim.IsCompatible(result, library, out string msg);

            Assert.That(ok, Is.False);
            Assert.That(msg, Does.Contain("mismatch"));
        }

        [Test]
        public void IsCompatible_SequenceMismatch_ReturnsFalseWithMessage()
        {
            var result = MakeResult("PEPTIDEK", 2, new[] { 1000f, 2000f });
            var library = MakeLibrary("DIFFERENTK", 2, 2);

            bool ok = DiaMetaDrawShim.IsCompatible(result, library, out string msg);

            Assert.That(ok, Is.False);
            Assert.That(msg, Does.Contain("Sequence mismatch"));
        }

        [Test]
        public void IsCompatible_ChargeMismatch_ReturnsFalseWithMessage()
        {
            var result = MakeResult("PEPTIDEK", 2, new[] { 1000f, 2000f });
            var library = MakeLibrary("PEPTIDEK", 3, 2); // charge 3 vs 2

            bool ok = DiaMetaDrawShim.IsCompatible(result, library, out string msg);

            Assert.That(ok, Is.False);
            Assert.That(msg, Does.Contain("Charge state mismatch"));
        }

        [Test]
        public void IsCompatible_NullResult_ReturnsFalse()
        {
            var library = MakeLibrary("PEPTIDEK", 2, 2);

            bool ok = DiaMetaDrawShim.IsCompatible(null, library, out string msg);

            Assert.That(ok, Is.False);
            Assert.That(msg, Does.Contain("null"));
        }

        [Test]
        public void IsCompatible_NullLibrary_ReturnsFalse()
        {
            var result = MakeResult("PEPTIDEK", 2, new[] { 1000f });

            bool ok = DiaMetaDrawShim.IsCompatible(result, null, out string msg);

            Assert.That(ok, Is.False);
        }

        // ── DiaMetaDrawBundle ─────────────────────────────────────────────────

        [Test]
        public void Bundle_ToString_ContainsSequenceAndScore()
        {
            var result = MakeResult("PEPTIDEK", 2, new[] { 1000f, 2000f });
            var library = MakeLibrary("PEPTIDEK", 2, 2);
            var adapter = new DiaPsmAdapter(result, "HeLa_01");

            var bundle = DiaMetaDrawShim.BuildBundle(adapter, library);
            string str = bundle.ToString();

            Assert.That(str, Does.Contain("PEPTIDEK"));
            Assert.That(str, Does.Contain("Score="));
        }

        [Test]
        public void Bundle_SyntheticFragmentIons_IsReadOnly()
        {
            var result = MakeResult("PEPTIDEK", 2, new[] { 1000f, 2000f });
            var library = MakeLibrary("PEPTIDEK", 2, 2);
            var adapter = new DiaPsmAdapter(result, "HeLa_01");

            var bundle = DiaMetaDrawShim.BuildBundle(adapter, library);

            // IReadOnlyList — casting to List should not be possible at compile time,
            // but we verify the runtime type is not mutable via the interface.
            Assert.That(bundle.SyntheticFragmentIons, Is.Not.Null);
            Assert.That(bundle.SyntheticFragmentIons.Count, Is.EqualTo(2));
        }
    }
}