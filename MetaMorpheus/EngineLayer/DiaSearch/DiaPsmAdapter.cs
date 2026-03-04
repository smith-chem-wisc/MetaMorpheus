// Copyright 2026 MetaMorpheus Contributors
// Licensed under the MIT License

using MassSpectrometry.Dia;
using System;

namespace EngineLayer.DiaSearch
{
    /// <summary>
    /// Adapts a <see cref="DiaSearchResult"/> to the flat-property surface expected by
    /// <see cref="DiaPsmTsvWriter"/> and MetaDraw's identification table.
    ///
    /// Design: Adapter pattern (not inheritance from SpectralMatch) because DIA results
    /// have no single representative scan — the scan number is the isolation window ID.
    ///
    /// Thread-safety: instances are read-only after construction (BaseSequence is
    /// lazily cached but behind a simple null-check; safe for concurrent reads after
    /// first access on a single thread).
    /// </summary>
    public class DiaPsmAdapter
    {
        private readonly DiaSearchResult _result;
        private string _baseSequence; // lazy cache

        /// <summary>
        /// Constructs an adapter wrapping the given result.
        /// </summary>
        /// <param name="result">The underlying DIA search result.</param>
        /// <param name="fileName">
        /// The raw file base name (no path, no extension) this result came from.
        /// </param>
        public DiaPsmAdapter(DiaSearchResult result, string fileName)
        {
            _result = result ?? throw new ArgumentNullException(nameof(result));
            FileName = fileName ?? string.Empty;
        }

        // ── Identity ──────────────────────────────────────────────────────────

        /// <summary>Base name of the raw file (no path, no extension).</summary>
        public string FileName { get; }

        /// <summary>
        /// Isolation window ID used as a proxy for scan number.
        /// MetaDraw will display this in the scan number column but cannot load
        /// the corresponding spectrum (window IDs are not real scan numbers).
        /// </summary>
        public int ScanNumber => _result.WindowId;

        /// <summary>
        /// Observed apex retention time in minutes.
        /// Falls back to the midpoint of the RT extraction window when apex RT
        /// was not resolved (e.g. fewer than 3 time points in the XIC).
        /// </summary>
        public double ScanRetentionTime
        {
            get
            {
                float apex = _result.ObservedApexRt;
                if (!float.IsNaN(apex) && apex > 0f)
                    return apex;
                return (_result.RtWindowStart + _result.RtWindowEnd) / 2.0;
            }
        }

        /// <summary>Precursor charge state.</summary>
        public int PrecursorCharge => _result.ChargeState;

        /// <summary>Precursor m/z (as stored; already charge-state adjusted).</summary>
        public double PrecursorMz => _result.PrecursorMz;

        // ── Sequence ──────────────────────────────────────────────────────────

        /// <summary>Full annotated sequence (may include modification brackets).</summary>
        public string FullSequence => _result.Sequence;

        /// <summary>
        /// Bare amino-acid sequence with modification annotations stripped.
        /// Lazily computed from <see cref="DiaSequenceParser"/> and cached.
        /// </summary>
        public string BaseSequence
        {
            get
            {
                if (_baseSequence == null)
                    _baseSequence = DiaSequenceParser.GetBaseSequence(_result.Sequence);
                return _baseSequence;
            }
        }

        // ── Scores and FDR ────────────────────────────────────────────────────

        /// <summary>
        /// Primary ranking score. Uses <see cref="DiaSearchResult.ClassifierScore"/>
        /// when available (post-FDR); falls back to <see cref="DiaSearchResult.SpectralAngle"/>.
        /// </summary>
        public double Score
        {
            get
            {
                float cs = _result.ClassifierScore;
                if (!float.IsNaN(cs)) return cs;
                return _result.SpectralAngle;
            }
        }

        /// <summary>
        /// Precursor-level q-value. Returns 1.0 when FDR has not yet been computed.
        /// </summary>
        public double QValue => _result.FdrInfo?.QValue ?? 1.0;

        /// <summary>
        /// Posterior error probability. Returns NaN when not computed.
        /// </summary>
        public double PEP => _result.FdrInfo?.PEP ?? double.NaN;

        /// <summary>
        /// PEP-based q-value. Returns NaN when not computed.
        /// </summary>
        public double PEP_QValue => _result.FdrInfo?.PEP_QValue ?? double.NaN;

        /// <summary>
        /// Peptide-level q-value (rollup across charge states).
        /// Returns null when not computed.
        /// </summary>
        public double? PeptideQValue => _result.FdrInfo?.PeptideQValue;

        /// <summary>"T" for target, "D" for decoy.</summary>
        public string TargetDecoyLabel => _result.IsDecoy ? "D" : "T";

        // ── Fragment evidence ────────────────────────────────────────────────

        /// <summary>Number of library fragments detected in the XIC.</summary>
        public int FragmentsDetected => _result.FragmentsDetected;

        /// <summary>Number of library fragments queried.</summary>
        public int FragmentsQueried => _result.FragmentsQueried;

        // ── DIA-specific scores ───────────────────────────────────────────────

        /// <summary>
        /// Spectral angle score. Coalesces the two storage properties:
        /// <see cref="DiaSearchResult.SpectralAngle"/> (temporal scoring path) and
        /// <see cref="DiaSearchResult.SpectralAngleScore"/> (simple assembler path).
        /// Returns NaN only when neither path populated the value.
        /// </summary>
        public float SpectralAngleScore
        {
            get
            {
                float sa = _result.SpectralAngle;
                if (!float.IsNaN(sa)) return sa;
                return _result.SpectralAngleScore;
            }
        }

        /// <summary>
        /// Dot product / temporal score. Coalesces
        /// <see cref="DiaSearchResult.TemporalScore"/> (temporal path) and
        /// <see cref="DiaSearchResult.DotProductScore"/> (simple assembler path).
        /// </summary>
        public float DotProductScore
        {
            get
            {
                float ts = _result.TemporalScore;
                if (!float.IsNaN(ts)) return ts;
                return _result.DotProductScore;
            }
        }

        /// <summary>Apex temporal correlation score.</summary>
        public float ApexScore => _result.ApexScore;

        /// <summary>Interference / co-isolation score (1.0 = no co-isolation).</summary>
        public float ChimericScore => _result.ChimericScore;

        // ── RT fields ────────────────────────────────────────────────────────

        /// <summary>Library-predicted retention time, if available.</summary>
        public double? LibraryRetentionTime => _result.LibraryRetentionTime;

        /// <summary>RT extraction window start (minutes).</summary>
        public float RtWindowStart => _result.RtWindowStart;

        /// <summary>RT extraction window end (minutes).</summary>
        public float RtWindowEnd => _result.RtWindowEnd;

        /// <summary>Absolute RT deviation from library prediction (minutes).</summary>
        public float RtDeviationMinutes => _result.RtDeviationMinutes;

        // ── Modification summary ──────────────────────────────────────────────

        /// <summary>
        /// Human-readable modification summary string, e.g. "Oxidation@5; Phospho@9".
        /// Empty string when no modifications are present.
        /// </summary>
        public string ModificationSummary =>
            DiaSequenceParser.FormatModificationSummary(_result.Sequence);

        // ── Pass-through ─────────────────────────────────────────────────────

        /// <summary>The underlying <see cref="DiaSearchResult"/> for callers that need raw access.</summary>
        public DiaSearchResult UnderlyingResult => _result;
    }
}