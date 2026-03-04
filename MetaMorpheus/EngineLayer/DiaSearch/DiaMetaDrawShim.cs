// Copyright 2026 MetaMorpheus Contributors
// Licensed under the MIT License

using MassSpectrometry.Dia;
using Omics.Fragmentation;
using Omics.SpectrumMatch;
using System;
using System.Collections.Generic;

namespace EngineLayer.DiaSearch
{
    /// <summary>
    /// Builds synthetic <see cref="MatchedFragmentIon"/> lists so that MetaDraw's
    /// fragment annotation panel can render DIA results.
    ///
    /// Background
    /// ----------
    /// DDA PSMs carry a <c>List&lt;MatchedFragmentIon&gt;</c> derived from a single MS2 scan.
    /// DIA results have no single representative scan; instead, evidence is accumulated across
    /// many MS2 scans as extracted ion chromatograms (XICs).  MetaDraw's fragment panel
    /// requires a <c>MatchedFragmentIon</c> list, so we synthesise one:
    ///
    ///   For each fragment ion i (0 … FragmentsQueried-1):
    ///     - Ion identity (ProductType, FragmentNumber, Charge, m/z) → LibrarySpectrum.MatchedFragmentIons[i]
    ///     - Observed intensity                                       → DiaSearchResult.ExtractedIntensities[i]
    ///
    /// Fragments that yielded zero XIC signal (ExtractedIntensities[i] == 0) are excluded
    /// from the returned list so MetaDraw does not annotate ghost peaks.
    ///
    /// The synthetic spectrum can be passed directly to
    /// <c>LibrarySpectrum.CalculateSpectralAngleOnTheFly</c> because that method only
    /// needs (Mz, Intensity) pairs.
    ///
    /// Thread-safety
    /// -------------
    /// All public methods are pure functions; no shared state is mutated.
    /// </summary>
    public static class DiaMetaDrawShim
    {
        // ── Public API ────────────────────────────────────────────────────────

        /// <summary>
        /// Builds a synthetic <see cref="MatchedFragmentIon"/> list for use in MetaDraw.
        ///
        /// Each returned ion carries:
        /// <list type="bullet">
        ///   <item>The library fragment's <see cref="Product"/> (ion type, number, neutral loss)</item>
        ///   <item>The library fragment's m/z</item>
        ///   <item>The <em>observed</em> XIC apex intensity from <paramref name="result"/></item>
        ///   <item>The library fragment's charge</item>
        /// </list>
        ///
        /// Fragments with zero observed intensity are omitted (they were not detected).
        /// </summary>
        /// <param name="result">
        /// The DIA search result.  <c>ExtractedIntensities</c> must be parallel to
        /// <paramref name="library"/>'s <c>MatchedFragmentIons</c>.
        /// </param>
        /// <param name="library">
        /// The <see cref="LibrarySpectrum"/> that was searched against.
        /// Must have the same fragment count as <paramref name="result"/>.<c>FragmentsQueried</c>.
        /// </param>
        /// <returns>
        /// Synthetic fragment ion list, sorted by m/z ascending (MetaDraw convention).
        /// Returns an empty list when the fragment arrays are mismatched or empty.
        /// </returns>
        /// <exception cref="ArgumentNullException">
        /// Thrown when <paramref name="result"/> or <paramref name="library"/> is null.
        /// </exception>
        public static List<MatchedFragmentIon> BuildSyntheticFragmentIons(
            DiaSearchResult result,
            LibrarySpectrum library)
        {
            if (result == null) throw new ArgumentNullException(nameof(result));
            if (library == null) throw new ArgumentNullException(nameof(library));

            var libraryIons = library.MatchedFragmentIons;

            // Guard: parallel arrays must agree in length.
            // A mismatch indicates a programming error (wrong library paired with result).
            if (libraryIons == null || libraryIons.Count == 0)
                return new List<MatchedFragmentIon>();

            int n = Math.Min(result.FragmentsQueried, libraryIons.Count);
            if (n == 0)
                return new List<MatchedFragmentIon>();

            var synthetic = new List<MatchedFragmentIon>(capacity: n);

            for (int i = 0; i < n; i++)
            {
                float observedIntensity = result.ExtractedIntensities[i];

                // Omit fragments with no observed signal — MetaDraw should only annotate
                // peaks that were actually detected in the XICs.
                if (observedIntensity <= 0f)
                    continue;

                var libIon = libraryIons[i];

                synthetic.Add(new MatchedFragmentIon(
                    libIon.NeutralTheoreticalProduct,
                    libIon.Mz,
                    observedIntensity,
                    libIon.Charge));
            }

            // MetaDraw expects ions sorted by m/z ascending.
            synthetic.Sort(static (a, b) => a.Mz.CompareTo(b.Mz));
            return synthetic;
        }

        /// <summary>
        /// Builds a normalised version of the synthetic fragment ion list where the most
        /// abundant observed fragment is set to intensity 1.0 and all others are scaled
        /// proportionally.
        ///
        /// Useful when the caller wants to render relative intensities (mirror plots, etc.)
        /// rather than raw XIC summed intensities.
        /// </summary>
        /// <param name="result">The DIA search result.</param>
        /// <param name="library">The corresponding library spectrum.</param>
        /// <returns>
        /// Normalised synthetic ion list, sorted by m/z ascending.
        /// Returns an empty list when no fragments were detected.
        /// </returns>
        public static List<MatchedFragmentIon> BuildNormalisedSyntheticFragmentIons(
            DiaSearchResult result,
            LibrarySpectrum library)
        {
            var ions = BuildSyntheticFragmentIons(result, library);

            if (ions.Count == 0)
                return ions;

            // Find the maximum observed intensity for normalisation.
            double maxIntensity = 0.0;
            foreach (var ion in ions)
            {
                if (ion.Intensity > maxIntensity)
                    maxIntensity = ion.Intensity;
            }

            if (maxIntensity <= 0.0)
                return ions;

            // Rebuild with normalised intensities.
            // MatchedFragmentIon is immutable after construction, so we create new instances.
            var normalised = new List<MatchedFragmentIon>(ions.Count);
            foreach (var ion in ions)
            {
                normalised.Add(new MatchedFragmentIon(
                    ion.NeutralTheoreticalProduct,
                    ion.Mz,
                    ion.Intensity / maxIntensity,
                    ion.Charge));
            }

            return normalised;
        }

        /// <summary>
        /// Convenience overload: builds the synthetic fragment list and attaches it to
        /// <paramref name="adapter"/> by returning both together as a
        /// <see cref="DiaMetaDrawBundle"/>.
        ///
        /// This is the primary entry point used by <c>PostDiaSearchAnalysisTask</c> when
        /// building the MetaDraw-compatible result set.
        /// </summary>
        /// <param name="adapter">The <see cref="DiaPsmAdapter"/> wrapping the result.</param>
        /// <param name="library">The corresponding library spectrum.</param>
        /// <param name="normalise">
        /// When <c>true</c>, intensities are normalised to the most abundant fragment.
        /// Default: <c>false</c> (raw XIC summed intensities).
        /// </param>
        /// <returns>A bundle pairing the adapter with its synthetic fragment ion list.</returns>
        public static DiaMetaDrawBundle BuildBundle(
            DiaPsmAdapter adapter,
            LibrarySpectrum library,
            bool normalise = false)
        {
            if (adapter == null) throw new ArgumentNullException(nameof(adapter));
            if (library == null) throw new ArgumentNullException(nameof(library));

            var ions = normalise
                ? BuildNormalisedSyntheticFragmentIons(adapter.UnderlyingResult, library)
                : BuildSyntheticFragmentIons(adapter.UnderlyingResult, library);

            return new DiaMetaDrawBundle(adapter, ions);
        }

        /// <summary>
        /// Verifies that a <see cref="DiaSearchResult"/> and a <see cref="LibrarySpectrum"/>
        /// are compatible for shim construction (same fragment count, matching sequence and
        /// charge state).
        ///
        /// Returns <c>true</c> if compatible. Populates <paramref name="diagnosticMessage"/>
        /// with a human-readable explanation when incompatible.
        /// </summary>
        public static bool IsCompatible(
            DiaSearchResult result,
            LibrarySpectrum library,
            out string diagnosticMessage)
        {
            if (result == null)
            {
                diagnosticMessage = "result is null";
                return false;
            }
            if (library == null)
            {
                diagnosticMessage = "library is null";
                return false;
            }
            if (library.MatchedFragmentIons == null)
            {
                diagnosticMessage = "library has null MatchedFragmentIons";
                return false;
            }

            int libCount = library.MatchedFragmentIons.Count;
            if (result.FragmentsQueried != libCount)
            {
                diagnosticMessage =
                    $"Fragment count mismatch: result.FragmentsQueried={result.FragmentsQueried} " +
                    $"but library has {libCount} ions. " +
                    $"Result: {result.Sequence}/{result.ChargeState}, " +
                    $"Library: {library.Name}";
                return false;
            }

            if (!string.Equals(result.Sequence, library.Sequence,
                    StringComparison.OrdinalIgnoreCase))
            {
                diagnosticMessage =
                    $"Sequence mismatch: result='{result.Sequence}' library='{library.Sequence}'";
                return false;
            }

            if (result.ChargeState != library.ChargeState)
            {
                diagnosticMessage =
                    $"Charge state mismatch: result={result.ChargeState} library={library.ChargeState}";
                return false;
            }

            diagnosticMessage = null;
            return true;
        }
    }

    // ── Supporting type ───────────────────────────────────────────────────────

    /// <summary>
    /// Pairs a <see cref="DiaPsmAdapter"/> with its synthetic fragment ion list,
    /// ready for MetaDraw rendering.
    ///
    /// Immutable value-like container; all fields are set at construction.
    /// </summary>
    public sealed class DiaMetaDrawBundle
    {
        /// <summary>The adapted DIA result (identification table data).</summary>
        public DiaPsmAdapter Adapter { get; }

        /// <summary>
        /// Synthetic fragment ions built from XIC intensities + library ion metadata.
        /// Sorted by m/z ascending. Empty when no fragments were detected.
        /// </summary>
        public IReadOnlyList<MatchedFragmentIon> SyntheticFragmentIons { get; }

        /// <summary>Number of synthetic ions (fragments detected in XICs).</summary>
        public int DetectedFragmentCount => SyntheticFragmentIons.Count;

        /// <summary>
        /// Whether there is at least one fragment ion to display in MetaDraw.
        /// </summary>
        public bool HasFragmentAnnotations => SyntheticFragmentIons.Count > 0;

        internal DiaMetaDrawBundle(
            DiaPsmAdapter adapter,
            List<MatchedFragmentIon> syntheticIons)
        {
            Adapter = adapter ?? throw new ArgumentNullException(nameof(adapter));
            SyntheticFragmentIons = syntheticIons
                ?? throw new ArgumentNullException(nameof(syntheticIons));
        }

        public override string ToString() =>
            $"{Adapter.FullSequence}/{Adapter.PrecursorCharge} " +
            $"RT={Adapter.ScanRetentionTime:F3} " +
            $"Fragments={DetectedFragmentCount} " +
            $"Score={Adapter.Score:F4}";
    }
}