// Copyright 2026 MetaMorpheus Contributors
// Licensed under the MIT License

using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using MassSpectrometry;
using MassSpectrometry.Dia;
using Omics.Fragmentation;
using Omics.SpectrumMatch;

namespace EngineLayer.DiaSearch
{
    /// <summary>
    /// MetaMorpheus engine that orchestrates the mzLib DIA extraction pipeline.
    /// 
    /// This is the facade that bridges MetaMorpheus types (LibrarySpectrum, CommonParameters)
    /// to the mzLib DIA engine types (LibraryPrecursorInput, DiaScanIndex, etc.).
    /// 
    /// Pipeline:
    ///   1. Convert LibrarySpectrum[] → LibraryPrecursorInput[] (resolves Omics↔MassSpec boundary)
    ///   2. Build DiaScanIndex from raw scans via DiaScanIndexBuilder
    ///   3. Generate fragment queries via DiaLibraryQueryGenerator.Generate()
    ///   4. Create extractor via FragmentExtractorFactory (CPU or GPU)
    ///   5. Run parallel extraction via DiaExtractionOrchestrator.ExtractAll()
    ///   6. Assemble + score results via DiaLibraryQueryGenerator.AssembleResults()
    ///   7. Return DiaEngineResults
    /// </summary>
    public class DiaEngine : MetaMorpheusEngine
    {
        private readonly MsDataScan[] _allScans;
        private readonly List<LibrarySpectrum> _librarySpectra;
        private readonly DiaSearchParameters _diaParameters;

        public DiaEngine(
            MsDataScan[] allScans,
            List<LibrarySpectrum> librarySpectra,
            DiaSearchParameters diaParameters,
            CommonParameters commonParameters,
            List<(string FileName, CommonParameters Parameters)> fileSpecificParameters,
            List<string> nestedIds)
            : base(commonParameters, fileSpecificParameters, nestedIds)
        {
            _allScans = allScans ?? throw new ArgumentNullException(nameof(allScans));
            _librarySpectra = librarySpectra ?? throw new ArgumentNullException(nameof(librarySpectra));
            _diaParameters = diaParameters ?? throw new ArgumentNullException(nameof(diaParameters));
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            // ── Step 1: Convert LibrarySpectrum → LibraryPrecursorInput ─────
            Status("Converting library spectra...");
            var precursorInputs = ConvertLibrarySpectra(_librarySpectra);
            int targetCount = precursorInputs.Count(p => !p.IsDecoy);
            int decoyCount = precursorInputs.Length - targetCount;

            // ── Step 2: Build SoA scan index ────────────────────────────────
            Status("Building DIA scan index...");
            var sw = Stopwatch.StartNew();
            using var scanIndex = DiaScanIndexBuilder.Build(_allScans);
            sw.Stop();

            Status($"Indexed {scanIndex.ScanCount} MS2 scans across {scanIndex.WindowCount} windows " +
                   $"({scanIndex.TotalPeakCount:N0} peaks) in {sw.Elapsed.TotalSeconds:F2}s");

            // ── Step 3: Generate fragment queries ───────────────────────────
            Status("Generating fragment queries...");
            var generation = DiaLibraryQueryGenerator.Generate(precursorInputs, scanIndex, _diaParameters);

            Status($"Generated {generation.Queries.Length:N0} queries for " +
                   $"{generation.PrecursorGroups.Length:N0} precursors " +
                   $"(skipped: {generation.SkippedNoWindow} no-window, {generation.SkippedNoFragments} no-fragments)");

            if (generation.Queries.Length == 0)
            {
                Warn("No fragment queries generated. Check that library precursor m/z values " +
                     "fall within the DIA isolation windows.");
                return new DiaEngineResults(this, new List<DiaSearchResult>(),
                    targetCount, decoyCount, scanIndex.ScanCount, scanIndex.WindowCount);
            }

            // ── Step 4: Create extractor (CPU or GPU) ───────────────────────
            string backend = FragmentExtractorFactory.DescribeBackend(!_diaParameters.PreferGpu);
            Status($"Extraction backend: {backend}");

            var extractorFactory = FragmentExtractorFactory.CreateFactory(!_diaParameters.PreferGpu);

            // ── Step 5: Parallel extraction via orchestrator ────────────────
            Status("Extracting fragment ion chromatograms...");
            sw.Restart();

            using var orchestrator = new DiaExtractionOrchestrator(scanIndex, extractorFactory);
            var extraction = orchestrator.ExtractAll(generation.Queries, _diaParameters.EffectiveMaxThreads);

            sw.Stop();
            Status($"Extraction complete: {generation.Queries.Length:N0} queries, " +
                   $"{extraction.TotalDataPoints:N0} XIC data points in {sw.Elapsed.TotalSeconds:F2}s");

            // ── Step 6: Score and assemble results ──────────────────────────
            Status("Scoring and assembling results...");
            var dotScorer = new NormalizedDotProductScorer();
            var angleScorer = new SpectralAngleScorer();

            var results = DiaLibraryQueryGenerator.AssembleResults(
                precursorInputs,
                generation,
                extraction.Results,
                _diaParameters,
                dotProductScorer: dotScorer,
                spectralAngleScorer: angleScorer);

            int targetResults = results.Count(r => !r.IsDecoy);
            int decoyResults = results.Count - targetResults;
            Status($"DIA search complete: {targetResults} target + {decoyResults} decoy results " +
                   $"from {_librarySpectra.Count} library entries");

            return new DiaEngineResults(this, results,
                targetCount, decoyCount, scanIndex.ScanCount, scanIndex.WindowCount);
        }

        /// <summary>
        /// Converts MetaMorpheus LibrarySpectrum objects to mzLib LibraryPrecursorInput structs.
        /// This is where the Omics → MassSpectrometry.Dia boundary crossing happens.
        /// 
        /// MatchedFragmentIon.Mz and .Intensity are doubles in Omics;
        /// LibraryPrecursorInput uses float[] for cache/GPU efficiency.
        /// </summary>
        internal static LibraryPrecursorInput[] ConvertLibrarySpectra(List<LibrarySpectrum> spectra)
        {
            var inputs = new LibraryPrecursorInput[spectra.Count];

            for (int i = 0; i < spectra.Count; i++)
            {
                var lib = spectra[i];
                var ions = lib.MatchedFragmentIons;

                var fragmentMzs = new float[ions.Count];
                var fragmentIntensities = new float[ions.Count];

                for (int j = 0; j < ions.Count; j++)
                {
                    fragmentMzs[j] = (float)ions[j].Mz;
                    fragmentIntensities[j] = (float)ions[j].Intensity;
                }

                inputs[i] = new LibraryPrecursorInput(
                    sequence: lib.Sequence,
                    precursorMz: lib.PrecursorMz,           // double → double (no conversion needed)
                    chargeState: lib.ChargeState,
                    retentionTime: lib.RetentionTime,       // double? → double? (no conversion needed)
                    isDecoy: lib.IsDecoy,
                    fragmentMzs: fragmentMzs,
                    fragmentIntensities: fragmentIntensities
                );
            }

            return inputs;
        }
    }
}
