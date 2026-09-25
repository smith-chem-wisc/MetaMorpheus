using Chemistry;
﻿using NUnit.Framework;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using EngineLayer;
using MassSpectrometry;
using MzLibUtil;
using Omics.Fragmentation;
using Omics.SpectrumMatch;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Readers.SpectralLibrary;
using PredictionClients.Koina.AbstractClasses;
using TaskLayer;

namespace Test
{
    /// <summary>
    /// Covers the spectral-angle path that <see cref="SupplementalSpectralSimilarityTests"/> does
    /// not reach. Those tests build PSMs without calling ResolveAllAmbiguities, so FullSequence is
    /// null, so GetSpectralMatchesWithoutComputedSpectralAngle filters every one of them out and
    /// ComputeSpectrumSimilarity returns before it does anything. Every test here resolves its PSMs
    /// first, which is what gets past that gate.
    ///
    /// Nothing here touches the network. The one part that genuinely needs Koina is the Predict call
    /// itself; the logic around it is separated into BuildPredictionInputs and ResolveCollisionEnergy
    /// so it can be asserted on directly.
    /// </summary>
    [TestFixture]
    public static class SupplementalSpectralSimilarityCoverageTests
    {
        private static readonly CommonParameters CommonParams = new(
            digestionParams: new DigestionParams(protease: "trypsin"),
            scoreCutoff: 1,
            productMassTolerance: new PpmTolerance(20),
            precursorMassTolerance: new PpmTolerance(5));

        /// <summary>
        /// A PSM whose ambiguities are resolved, so FullSequence is populated and the PSM actually
        /// survives the filter. hcdEnergy null means the scan carries no recorded energy.
        /// </summary>
        private static SpectralMatch ResolvedPsm(string sequence, int charge, string hcdEnergy = null)
        {
            // digested from a real protein rather than constructed directly: ResolveAllAmbiguities
            // reaches through to the parent to decide decoy status, and a parentless peptide NREs
            var peptide = new Protein(sequence, "accession")
                .Digest(CommonParams.DigestionParams, new List<Omics.Modifications.Modification>(),
                        new List<Omics.Modifications.Modification>())
                .First(p => p.BaseSequence == sequence);

            var products = new List<Product>();
            peptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
            var matched = products.Take(5)
                .Select(p => new MatchedFragmentIon(p, p.ToMz(1), 100, 1))
                .ToList();

            var scan = new MsDataScan(
                massSpectrum: new MzSpectrum(matched.Select(m => m.Mz).ToArray(), matched.Select(_ => 100.0).ToArray(), false),
                oneBasedScanNumber: 1, msnOrder: 2, isCentroid: true, polarity: Polarity.Positive,
                retentionTime: 1.0, scanWindowRange: new MzRange(50, 2000), scanFilter: "FTMS + p NSI",
                mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: 1000, injectionTime: null,
                noiseData: null, nativeId: "scan=1", hcdEnergy: hcdEnergy);

            var psm = new PeptideSpectralMatch(peptide, 0, 50, 0,
                new Ms2ScanWithSpecificMass(scan, peptide.MonoisotopicMass.ToMz(charge), charge, "test.raw", CommonParams),
                CommonParams, matched);

            psm.ResolveAllAmbiguities();

            // Ms2Scan is opt-in on SpectralMatch and only some paths populate it, so set it
            // explicitly when the test is about the recorded collision energy.
            if (hcdEnergy != null)
            {
                psm.SetMs2Scan(scan);
            }

            psm.SpectralAngle = -1;
            return psm;
        }

        private static PostSearchAnalysisTask TaskWith(params SpectralMatch[] psms) =>
            new()
            {
                Parameters = new PostSearchAnalysisParameters
                {
                    SearchParameters = new SearchParameters(),
                    AllSpectralMatches = psms.ToList(),
                    OutputFolder = Path.GetTempPath()
                }
            };

        private static LibrarySpectrum LibrarySpectrumFor(SpectralMatch psm) =>
            new(psm.FullSequence, psm.ScanPrecursorMonoisotopicPeakMz, psm.ScanPrecursorCharge,
                psm.MatchedFragmentIons, 1.0);

        /// <summary>
        /// A real, file-backed library. SpectralLibrary.GetAllLibrarySpectra reads spectra from disk
        /// by byte offset, so setting Results on a default-constructed instance does not produce a
        /// readable library - it throws on enumeration.
        /// </summary>
        private static SpectralLibrary LibraryOf(params SpectralMatch[] psms)
        {
            string path = Path.Combine(Path.GetTempPath(), "spectralAngleCoverage_" + System.Guid.NewGuid().ToString("N") + ".msp");
            File.WriteAllLines(path, psms.Select(psm => LibrarySpectrumFor(psm).ToString()));
            return new SpectralLibrary(new List<string> { path });
        }

        // ---------- the scoring loop ----------

        /// <summary>
        /// A PSM the library covers must come out with a real angle. This is the feature working:
        /// before this test nothing asserted that a spectral angle is ever actually computed.
        /// </summary>
        [Test]
        public static void PsmCoveredByTheLibraryGetsAComputedAngle()
        {
            var psm = ResolvedPsm("PEPTIDEK", 2);
            var task = TaskWith(psm);

            task.ComputeSpectrumSimilarity(LibraryOf(psm));

            Assert.That(psm.SpectralAngle, Is.GreaterThan(-1),
                "a PSM present in the library must be scored, not left at the sentinel");
        }

        /// <summary>
        /// The other arm of the same loop: the lookup is non-empty, so the method does not return
        /// early, but this PSM is not in it. It must be left at the sentinel rather than scored
        /// against somebody else's spectrum.
        /// </summary>
        [Test]
        public static void PsmMissingFromANonEmptyLookupKeepsTheSentinel()
        {
            var covered = ResolvedPsm("PEPTIDEK", 2);
            var uncovered = ResolvedPsm("ELVISLIVESK", 3);
            var task = TaskWith(covered, uncovered);
            task.Parameters.SearchParameters.UsePredictedSpectraForSpectralAngle = false;

            task.ComputeSpectrumSimilarity(LibraryOf(covered));

            Assert.That(covered.SpectralAngle, Is.GreaterThan(-1));
            Assert.That(uncovered.SpectralAngle, Is.EqualTo(-1),
                "a PSM the lookup does not contain must not borrow another peptide's spectrum");
        }

        /// <summary>
        /// SpectralAngle 0 is a legitimate (terrible) score, and only a negative value means
        /// "not computed". A PSM that already has an angle must not be recomputed.
        /// </summary>
        [Test]
        public static void AlreadyScoredPsmsAreLeftAlone()
        {
            var psm = ResolvedPsm("PEPTIDEK", 2);
            psm.SpectralAngle = 0;
            var task = TaskWith(psm);

            task.ComputeSpectrumSimilarity(LibraryOf(psm));

            Assert.That(psm.SpectralAngle, Is.EqualTo(0), "zero is a score, not a sentinel");
        }

        // ---------- what would be sent for prediction ----------

        /// <summary>
        /// Only the PSMs the library did not cover are requested. Predicting one we already have a
        /// real spectrum for would be a wasted call and a worse answer.
        /// </summary>
        [Test]
        public static void OnlyPsmsTheLibraryMissedAreQueuedForPrediction()
        {
            var covered = ResolvedPsm("PEPTIDEK", 2);
            var uncovered = ResolvedPsm("ELVISLIVESK", 3);
            var lookup = new Dictionary<(string, int), LibrarySpectrum>
            {
                [(covered.FullSequence, covered.ScanPrecursorCharge)] = LibrarySpectrumFor(covered)
            };

            var inputs = PostSearchAnalysisTask.BuildPredictionInputs(
                new List<SpectralMatch> { covered, uncovered }, lookup);

            Assert.That(inputs.Keys, Is.EquivalentTo(new[] { (uncovered.FullSequence, 3) }));
        }

        /// <summary>
        /// The same peptide and charge seen twice is one request, not two. Duplicates are the normal
        /// case: the same peptide is routinely matched in many scans.
        /// </summary>
        [Test]
        public static void DuplicateSequenceAndChargeIsRequestedOnce()
        {
            var first = ResolvedPsm("PEPTIDEK", 2);
            var second = ResolvedPsm("PEPTIDEK", 2);

            var inputs = PostSearchAnalysisTask.BuildPredictionInputs(
                new List<SpectralMatch> { first, second },
                new Dictionary<(string, int), LibrarySpectrum>());

            Assert.That(inputs, Has.Count.EqualTo(1));
        }

        /// <summary>
        /// The same peptide at a different charge is a different spectrum and must be requested
        /// separately - the dedup key has to include the charge.
        /// </summary>
        [Test]
        public static void SameSequenceAtADifferentChargeIsRequestedSeparately()
        {
            var inputs = PostSearchAnalysisTask.BuildPredictionInputs(
                new List<SpectralMatch> { ResolvedPsm("PEPTIDEK", 2), ResolvedPsm("PEPTIDEK", 3) },
                new Dictionary<(string, int), LibrarySpectrum>());

            Assert.That(inputs, Has.Count.EqualTo(2));
        }

        [Test]
        public static void CollisionEnergyComesFromTheScanWhenItRecordedOne()
        {
            Assert.That(PostSearchAnalysisTask.ResolveCollisionEnergy(ResolvedPsm("PEPTIDEK", 2, hcdEnergy: "27")),
                Is.EqualTo(27));
        }

        /// <summary>
        /// 30 is a fallback, not a measurement. Pinned so a change to it is a deliberate decision
        /// about what energy unlabelled spectra are predicted at, rather than an accident.
        /// </summary>
        [Test]
        public static void CollisionEnergyFallsBackToThirtyWithoutARecordedEnergy()
        {
            Assert.That(PostSearchAnalysisTask.ResolveCollisionEnergy(ResolvedPsm("PEPTIDEK", 2, hcdEnergy: null)),
                Is.EqualTo(30));
        }

        [Test]
        public static void CollisionEnergyFallsBackWhenTheScanEnergyIsNotANumber()
        {
            Assert.That(PostSearchAnalysisTask.ResolveCollisionEnergy(ResolvedPsm("PEPTIDEK", 2, hcdEnergy: "not a number")),
                Is.EqualTo(30));
        }

        // ---------- the library lifecycle ----------

        /// <summary>
        /// SearchTask closes the spectral library before post-search analysis runs, unless the run is
        /// updating it. Reading it here therefore hits a closed file on an ordinary search with a
        /// supplied library. That must cost the angles, not the whole completed search.
        /// </summary>
        [Test]
        public static void ClosedLibraryDegradesInsteadOfLosingTheSearch()
        {
            var psm = ResolvedPsm("PEPTIDEK", 2);
            var task = TaskWith(psm);

            var library = LibraryOf(psm);
            library.CloseConnections();

            Assert.DoesNotThrow(() => task.ComputeSpectrumSimilarity(library));
            Assert.That(psm.SpectralAngle, Is.EqualTo(-1),
                "a closed library yields no spectra, so the PSM keeps the sentinel");
        }

        // ---------- filing predictions back under the key PSMs ask with ----------

        /// <summary>
        /// Predictions come back keyed by the UNIMOD form the service accepted, while PSMs look
        /// themselves up by their own FullSequence. Without the mapping the spectrum is in the
        /// lookup under a key nothing queries - present, and invisible.
        /// </summary>
        [Test]
        public static void PredictedSpectrumIsFiledUnderTheOriginalSequence()
        {
            var psm = ResolvedPsm("PEPTIDEK", 2);
            var prediction = new PeptideFragmentIntensityPrediction(
                FullSequence: "PEPTIDEK",
                ValidatedFullSequence: "PEPTIDEK[UNIMOD:1]",
                PrecursorCharge: 2,
                FragmentAnnotations: new List<string>(),
                FragmentMZs: new List<double>(),
                FragmentIntensities: new List<double>());

            var predicted = new LibrarySpectrum("PEPTIDEK[UNIMOD:1]", 500, 2, psm.MatchedFragmentIons, 10);
            var lookup = new Dictionary<(string, int), LibrarySpectrum>();

            PostSearchAnalysisTask.MergePredictedSpectra(new[] { prediction }, new[] { predicted }, lookup);

            Assert.That(lookup.ContainsKey(("PEPTIDEK", 2)), Is.True,
                "filed under the validated sequence, no PSM would ever find it");
        }

        /// <summary>A measured spectrum beats a predicted one for the same peptide and charge.</summary>
        [Test]
        public static void RealLibrarySpectraAreNotOverwrittenByPredictions()
        {
            var psm = ResolvedPsm("PEPTIDEK", 2);
            var real = LibrarySpectrumFor(psm);
            var lookup = new Dictionary<(string, int), LibrarySpectrum> { [("PEPTIDEK", 2)] = real };

            var prediction = new PeptideFragmentIntensityPrediction(
                FullSequence: "PEPTIDEK", ValidatedFullSequence: "PEPTIDEK", PrecursorCharge: 2,
                FragmentAnnotations: new List<string>(), FragmentMZs: new List<double>(),
                FragmentIntensities: new List<double>());
            var predicted = new LibrarySpectrum("PEPTIDEK", 500, 2, psm.MatchedFragmentIons, 10);

            PostSearchAnalysisTask.MergePredictedSpectra(new[] { prediction }, new[] { predicted }, lookup);

            Assert.That(lookup[("PEPTIDEK", 2)], Is.SameAs(real));
        }

        /// <summary>
        /// An unconvertible sequence comes back with a null ValidatedFullSequence, which must not
        /// blow up the mapping. Relevant right now: the model is configured to return null for these.
        /// </summary>
        [Test]
        public static void NullValidatedSequencesAreSkippedRatherThanThrowing()
        {
            var psm = ResolvedPsm("PEPTIDEK", 2);
            var prediction = new PeptideFragmentIntensityPrediction(
                FullSequence: "PEPTIDEK", ValidatedFullSequence: null, PrecursorCharge: 2,
                FragmentAnnotations: new List<string>(), FragmentMZs: new List<double>(),
                FragmentIntensities: new List<double>());
            var predicted = new LibrarySpectrum("PEPTIDEK", 500, 2, psm.MatchedFragmentIons, 10);
            var lookup = new Dictionary<(string, int), LibrarySpectrum>();

            Assert.DoesNotThrow(() => PostSearchAnalysisTask.MergePredictedSpectra(
                new[] { prediction }, new[] { predicted }, lookup));
            Assert.That(lookup.ContainsKey(("PEPTIDEK", 2)), Is.True);
        }

        // ---------- the opt-in gate ----------

        /// <summary>
        /// With predictions off, a PSM the library does not cover is simply not scored. The library
        /// half still works, which is the point of gating only the prediction.
        /// </summary>
        [Test]
        public static void PredictionsOffStillScoresFromTheLibrary()
        {
            var covered = ResolvedPsm("PEPTIDEK", 2);
            var uncovered = ResolvedPsm("ELVISLIVESK", 3);
            var task = TaskWith(covered, uncovered);
            task.Parameters.SearchParameters.UsePredictedSpectraForSpectralAngle = false;

            task.ComputeSpectrumSimilarity(LibraryOf(covered));

            Assert.That(covered.SpectralAngle, Is.GreaterThan(-1), "the library half must not be gated");
            Assert.That(uncovered.SpectralAngle, Is.EqualTo(-1));
        }

        /// <summary>
        /// No library and predictions off means there is nothing to score against at all, and the
        /// method must say so by leaving the sentinel rather than by throwing.
        /// </summary>
        [Test]
        public static void NoLibraryAndPredictionsOffLeavesEverySentinelInPlace()
        {
            var psm = ResolvedPsm("PEPTIDEK", 2);
            var task = TaskWith(psm);
            task.Parameters.SearchParameters.UsePredictedSpectraForSpectralAngle = false;

            Assert.DoesNotThrow(() => task.ComputeSpectrumSimilarity(null));
            Assert.That(psm.SpectralAngle, Is.EqualTo(-1));
        }

    }
}
