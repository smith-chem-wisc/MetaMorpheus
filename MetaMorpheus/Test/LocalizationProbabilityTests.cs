using Chemistry;
using EngineLayer;
using EngineLayer.FdrAnalysis;
using EngineLayer.Localization;
using MassSpectrometry;
using NUnit.Framework;
using Omics;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    [TestFixture]
    public static class LocalizationProbabilityTests
    {
        private static readonly Modification PhosphoS = Make("S");
        private static readonly Modification PhosphoT = Make("T");
        private static readonly Product AnyProduct = new Product(ProductType.b, FragmentationTerminus.N, 100, 1, 1, 0);
        private static readonly CommonParameters Parameters = new CommonParameters(digestionParams: new DigestionParams(minPeptideLength: 1));

        private static readonly LocalizationPlacement OnS3 = new(new[] { new LocalizationSite(2, 'S', PhosphoS) });
        private static readonly LocalizationPlacement OnT4 = new(new[] { new LocalizationSite(3, 'T', PhosphoT) });
        private static readonly LocalizationPlacement Unmodified = new(Array.Empty<LocalizationSite>());

        private static Modification Make(string residue)
        {
            ModificationMotif.TryGetMotif(residue, out ModificationMotif motif);
            return new Modification(_originalId: "Phosphorylation", _modificationType: "Common Biological", _target: motif,
                _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H1O3P1"));
        }

        /// <summary>
        /// Evidence whose peaks (relative intensity, and which candidates matched each) are given directly.
        /// </summary>
        private static LocalizationEvidence Evidence(LocalizationPlacement[] placements, params (double Intensity, int[] MatchedBy)[] peaks)
        {
            var builder = new LocalizationEvidenceBuilder(basePeakIntensity: 1);
            for (int c = 0; c < placements.Length; c++)
            {
                var ions = peaks
                    .Select((peak, i) => (peak, i))
                    .Where(x => x.peak.MatchedBy.Contains(c))
                    .Select(x => new MatchedFragmentIon(AnyProduct, 100 + x.i, x.peak.Intensity, 1))
                    .ToList();
                builder.AddCandidate(placements[c], 0, ions);
            }
            return builder.Build(massDifference: 0, peptideLength: 6, isProteinNTerminal: false, isProteinCTerminal: false);
        }

        private static Ms2ScanWithSpecificMass _scan;
        private static PeptideWithSetModifications _unmodified;
        private static PeptideWithSetModifications _phosphorylated;

        private static SpectralMatch Psm(LocalizationEvidence evidence, bool modified, double qValue = 0.001, bool decoy = false)
        {
            if (_scan == null)
            {
                _unmodified = (PeptideWithSetModifications)new Protein("PESTSK", "P1")
                    .Digest(Parameters.DigestionParams, new List<Modification>(), new List<Modification>()).First();
                _phosphorylated = (PeptideWithSetModifications)PositionalIsoforms.WithModifications(_unmodified, new() { { 4, PhosphoS } });
                MsDataFile file = new TestDataFile(new List<IBioPolymerWithSetMods> { _unmodified });
                _scan = new Ms2ScanWithSpecificMass(file.GetAllScansList().Last(), _unmodified.MonoisotopicMass.ToMz(1), 1, null, new CommonParameters());
            }

            PeptideWithSetModifications peptide = modified ? _phosphorylated : _unmodified;
            if (decoy)
            {
                var decoyProtein = new Protein("PESTSK", "DECOY_P1", isDecoy: true);
                peptide = (PeptideWithSetModifications)decoyProtein.Digest(Parameters.DigestionParams, new List<Modification>(), new List<Modification>()).First();
            }

            SpectralMatch psm = new PeptideSpectralMatch(peptide, 0, 10, 0, _scan, Parameters, new List<MatchedFragmentIon>());
            psm.ResolveAllAmbiguities();
            psm.PsmFdrInfo = new FdrInfo { QValue = qValue };
            psm.LocalizationEvidence = evidence;
            return psm;
        }

        /// <summary>
        /// Confident unmodified matches whose matched peaks spread evenly from 0.02 to 1, for F to be learned from.
        /// </summary>
        private static List<SpectralMatch> UnmodifiedMatches(int count = 20)
        {
            var peaks = Enumerable.Range(1, 50).Select(i => (i / 50.0, new[] { 0 })).ToArray();
            return Enumerable.Range(0, count).Select(_ => Psm(Evidence(new[] { Unmodified }, peaks), modified: false)).ToList();
        }

        private static void Run(List<SpectralMatch> psms, int minMatchesForPrior = LocalizationProbabilityEngine.DefaultMinMatchesForPrior)
        {
            new LocalizationProbabilityEngine(psms, Parameters, new List<(string, CommonParameters)>(), new List<string>())
            {
                MinMatchesForPrior = minMatchesForPrior
            }.Run();
        }

        [Test]
        public static void TheIntensityDistributionIsAnEmpiricalCdfStrictlyInsideZeroAndOne()
        {
            var model = new MatchedPeakIntensityModel(new[] { 0.005, 0.1, 0.2, 0.3, 0.4 });

            Assert.That(model.TrainingPeakCount, Is.EqualTo(4)); // 0.005 is fainter than the 0.01 floor
            Assert.That(model.CumulativeProbability(0.05), Is.EqualTo(0.5 / 5));
            Assert.That(model.CumulativeProbability(0.2), Is.EqualTo(2.5 / 5));
            Assert.That(model.CumulativeProbability(1), Is.EqualTo(4.5 / 5));
        }

        [Test]
        public static void OnlyPeaksThatTellPlacementsApartCount()
        {
            var model = new MatchedPeakIntensityModel(Enumerable.Range(1, 99).Select(i => i / 100.0));
            var evidence = Evidence(new[] { OnS3, OnT4 },
                (0.9, new[] { 0 }),       // intense, only S3 explains it
                (0.5, new[] { 0, 1 }),    // both explain it: cancels
                (0.005, new[] { 1 }));    // below the floor: ignored

            double[] logLikelihoods = LocalizationProbabilityEngine.LogLikelihoods(evidence, model);
            double f = model.CumulativeProbability((float)0.9); // evidence keeps intensities as floats

            Assert.That(logLikelihoods[0], Is.EqualTo(Math.Log(f)).Within(1e-12));
            Assert.That(logLikelihoods[1], Is.EqualTo(Math.Log(1 - f)).Within(1e-12));
        }

        [Test]
        public static void ThePriorLearnedFromClearMatchesDecidesTheUnclearOnes()
        {
            // 40 matches the spectrum puts on S3, 10 it cannot place (S3 and T4 match the same peaks)
            var clear = Enumerable.Range(0, 40).Select(_ => Psm(Evidence(new[] { OnS3, OnT4 }, (0.8, new[] { 0 }), (0.3, new[] { 0, 1 })), modified: true));
            var unclear = Enumerable.Range(0, 10).Select(_ => Psm(Evidence(new[] { OnS3, OnT4 }, (0.3, new[] { 0, 1 })), modified: true)).ToList();
            var psms = UnmodifiedMatches().Concat(clear).Concat(unclear).ToList();

            // with a flat prior an unclear match is a coin toss
            Run(psms, minMatchesForPrior: int.MaxValue);
            Assert.That(unclear[0].LocalizationProbabilities.PlacementProbabilities[0], Is.EqualTo(0.5).Within(1e-9));

            // with the learned prior, phosphorylation is far likelier on S than T, so S3 is favoured
            Run(psms);
            LocalizationProbabilities learned = unclear[0].LocalizationProbabilities;
            Assert.That(learned.BestPlacement, Is.SameAs(OnS3));
            Assert.That(learned.BestProbability, Is.GreaterThan(0.75));
            Assert.That(learned.SiteProbability(2), Is.EqualTo(learned.BestProbability));
            Assert.That(learned.PlacementProbabilities.Sum(), Is.EqualTo(1).Within(1e-9));
        }

        [Test]
        public static void ClearMatchesHaveLowLocalizationQValuesAndUnclearOnesHigher()
        {
            var clear = Enumerable.Range(0, 40).Select(_ => Psm(Evidence(new[] { OnS3, OnT4 }, (0.8, new[] { 0 }), (0.6, new[] { 0 })), modified: true)).ToList();
            var unclear = Enumerable.Range(0, 10).Select(_ => Psm(Evidence(new[] { OnS3, OnT4 }, (0.3, new[] { 0, 1 })), modified: true)).ToList();
            var unmodified = UnmodifiedMatches();
            Run(unmodified.Concat(clear).Concat(unclear).ToList(), minMatchesForPrior: int.MaxValue);

            double clearQ = clear[0].LocalizationProbabilities.QValue.Value;
            double unclearQ = unclear[0].LocalizationProbabilities.QValue.Value;
            Assert.That(clearQ, Is.LessThan(unclearQ));
            Assert.That(clear.All(p => p.LocalizationProbabilities.QValue == clearQ)); // tied probabilities share a q-value

            // q is the false localization rate of everything at least as probable: (40 x clear error + 10 x 0.5) / 50
            double clearError = 1 - clear[0].LocalizationProbabilities.BestProbability;
            Assert.That(clearQ, Is.EqualTo(clearError).Within(1e-12));
            Assert.That(unclearQ, Is.EqualTo((40 * clearError + 10 * 0.5) / 50).Within(1e-12));

            // a match with one placement is not in question and has no q-value
            Assert.That(unmodified[0].LocalizationProbabilities.IsLocalizationInQuestion, Is.False);
            Assert.That(unmodified[0].LocalizationProbabilities.QValue, Is.Null);
        }

        [Test]
        public static void DecoysAndUnconfidentMatchesGetProbabilitiesButNoQValue()
        {
            var placementEvidence = Evidence(new[] { OnS3, OnT4 }, (0.8, new[] { 0 }));
            var decoy = Psm(placementEvidence, modified: true, decoy: true);
            var unconfident = Psm(placementEvidence, modified: true, qValue: 0.2);
            var target = Psm(placementEvidence, modified: true);

            Run(UnmodifiedMatches().Concat(new[] { decoy, unconfident, target }).ToList());

            Assert.That(decoy.LocalizationProbabilities.BestPlacement, Is.SameAs(OnS3));
            Assert.That(decoy.LocalizationProbabilities.QValue, Is.Null);
            Assert.That(unconfident.LocalizationProbabilities.QValue, Is.Null);
            Assert.That(target.LocalizationProbabilities.QValue, Is.Not.Null);
        }

        [Test]
        public static void TooFewUnmodifiedPeaksLeavesProbabilitiesUncomputed()
        {
            var target = Psm(Evidence(new[] { OnS3, OnT4 }, (0.8, new[] { 0 })), modified: true);
            var psms = UnmodifiedMatches(count: 1).Concat(new[] { target }).ToList(); // 50 peaks of which 0.02 is the faintest

            new LocalizationProbabilityEngine(psms, Parameters, new List<(string, CommonParameters)>(), new List<string>()) { MinTrainingPeaks = 51 }.Run();

            Assert.That(target.LocalizationProbabilities, Is.Null);
        }
    }
}
