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

        private static LocalizationProbabilityResults Run(List<SpectralMatch> psms, int minMatchesForPrior = LocalizationProbabilityEngine.DefaultMinMatchesForPrior,
            LocalizationLikelihood likelihood = LocalizationLikelihood.PtminerCdf)
        {
            return (LocalizationProbabilityResults)new LocalizationProbabilityEngine(psms, Parameters, new List<(string, CommonParameters)>(), new List<string>())
            {
                MinMatchesForPrior = minMatchesForPrior,
                Likelihood = likelihood
            }.Run();
        }

        /// <summary>
        /// Evidence that also records how many fragments each candidate predicted and the chance-peak coverage of the
        /// spectrum, for the intensity-ratio likelihood.
        /// </summary>
        private static LocalizationEvidence RatioEvidence(LocalizationPlacement[] placements, int[] predicted, float[] noiseCoverage,
            params (double Intensity, int[] MatchedBy)[] peaks)
        {
            var builder = new LocalizationEvidenceBuilder(basePeakIntensity: 1);
            for (int c = 0; c < placements.Length; c++)
            {
                var ions = peaks
                    .Select((peak, i) => (peak, i))
                    .Where(x => x.peak.MatchedBy.Contains(c))
                    .Select(x => new MatchedFragmentIon(AnyProduct, 100 + x.i, x.peak.Intensity, 1))
                    .ToList();
                builder.AddCandidate(placements[c], 0, ions, predicted[c]);
            }
            if (noiseCoverage != null)
            {
                builder.SetNoiseCoverage(noiseCoverage);
            }
            return builder.Build(massDifference: 0, peptideLength: 6, isProteinNTerminal: false, isProteinCTerminal: false);
        }

        /// <summary>
        /// Chance-peak coverage typical of a fragment spectrum: unexplained peaks are mostly faint, and cover 5% of the
        /// fragment mass range in all.
        /// </summary>
        private static float[] FaintNoise()
        {
            var coverage = new float[PeakIntensityRatioModel.BinCount];
            int first = PeakIntensityRatioModel.BinOf(0.01);
            coverage[first] = 0.02f;      // 1-1.8%
            coverage[first + 1] = 0.015f; // 1.8-3.2%
            coverage[first + 2] = 0.008f; // 3.2-5.6%
            coverage[first + 3] = 0.004f; // 5.6-10%
            coverage[first + 4] = 0.002f; // 10-18%
            coverage[first + 5] = 0.001f; // 18-32%
            return coverage;
        }

        /// <summary>
        /// Confident unmodified matches whose 50 matched peaks spread evenly from 0.02 to 1, out of 100 predicted
        /// fragments (d = 0.5), and whose unexplained peaks are <see cref="FaintNoise"/>.
        /// </summary>
        private static List<SpectralMatch> UnmodifiedRatioMatches(int count = 20)
        {
            var peaks = Enumerable.Range(1, 50).Select(i => (i / 50.0, new[] { 0 })).ToArray();
            return Enumerable.Range(0, count).Select(_ => Psm(RatioEvidence(new[] { Unmodified }, new[] { 100 }, FaintNoise(), peaks), modified: false)).ToList();
        }

        /// <summary>The intensity-ratio model the matches of <see cref="UnmodifiedRatioMatches"/> train.</summary>
        private static PeakIntensityRatioModel TrainedRatioModel()
        {
            var matched = Enumerable.Range(1, 50).Select(i => (double)(float)(i / 50.0));
            return new PeakIntensityRatioModel(Enumerable.Repeat(matched, 20).SelectMany(m => m),
                Enumerable.Range(0, 20).Select(_ => (IReadOnlyList<float>)FaintNoise()), predictedFragments: 2000, observedFragments: 1000);
        }

        /// <summary>
        /// The pathology of PTMiner's likelihood on NMDPLNDNIATLLHQSSDK +15.995: one placement explains 27 faint
        /// differential peaks, another explains 9 just as faint; both predict 40 fragments.
        /// </summary>
        private static LocalizationEvidence ManyFaintAgainstFewFaint()
        {
            var peaks = new List<(double, int[])>();
            for (int i = 0; i < 27; i++)
            {
                peaks.Add((0.011 + 0.089 * i / 26, new[] { 0 })); // 1% to 10% of the base peak
            }
            for (int i = 0; i < 9; i++)
            {
                peaks.Add((0.011 + 0.089 * i / 8, new[] { 1 }));
            }
            peaks.Add((0.8, new[] { 0, 1 })); // shared peaks cancel
            peaks.Add((0.5, new[] { 0, 1 }));
            return RatioEvidence(new[] { OnS3, OnT4 }, new[] { 40, 40 }, FaintNoise(), peaks.ToArray());
        }

        [Test]
        public static void TheDensityRatioComparesMatchedToUnexplainedPeakIntensities()
        {
            // the same intensities matched as unexplained: no intensity favours a match
            var flat = new float[PeakIntensityRatioModel.BinCount];
            for (int b = PeakIntensityRatioModel.BinOf(0.01); b < flat.Length; b++)
            {
                flat[b] = 0.01f;
            }
            var evenlySpread = Enumerable.Range(PeakIntensityRatioModel.BinOf(0.01), 8)
                .Select(b => Math.Pow(10, PeakIntensityRatioModel.LogMin + (b + 0.5) * PeakIntensityRatioModel.BinWidth));
            var same = new PeakIntensityRatioModel(evenlySpread, new[] { flat }, predictedFragments: 10, observedFragments: 4);
            Assert.That(same.IsTrained);
            Assert.That(same.ObservationProbability, Is.EqualTo(0.4));
            foreach (double intensity in new[] { 0.01, 0.03, 0.2, 1 })
            {
                Assert.That(same.LogDensityRatio(intensity), Is.EqualTo(0).Within(1e-12));
                // a match counts (1 - d) + d / q: here q = 8 bins x 0.01
                Assert.That(same.LogExplainedPeak(intensity, 0.08), Is.EqualTo(Math.Log(0.6 + 0.4 / 0.08)).Within(1e-12));
            }
            Assert.That(same.LogMissingFragment, Is.EqualTo(Math.Log(0.6)).Within(1e-12));

            // matched peaks run intense, unexplained ones faint: the ratio climbs with intensity, through 0
            PeakIntensityRatioModel model = TrainedRatioModel();
            Assert.That(model.LogDensityRatio(0.015), Is.LessThan(0));
            Assert.That(model.LogDensityRatio(0.9), Is.GreaterThan(0));
            double[] ratios = new[] { 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1 }.Select(model.LogDensityRatio).ToArray();
            Assert.That(ratios, Is.Ordered);

            // halfway (in log intensity) between two bin centres the log ratio is halfway between theirs
            double centre10 = Math.Pow(10, PeakIntensityRatioModel.LogMin + 10.5 * PeakIntensityRatioModel.BinWidth);
            double centre11 = Math.Pow(10, PeakIntensityRatioModel.LogMin + 11.5 * PeakIntensityRatioModel.BinWidth);
            Assert.That(model.LogDensityRatio(Math.Sqrt(centre10 * centre11)),
                Is.EqualTo((model.LogDensityRatio(centre10) + model.LogDensityRatio(centre11)) / 2).Within(1e-9));

            // densities are floored: no ratio is infinite, however one-sided the training, and a match never counts
            // for less than a predicted fragment left unobserved
            Assert.That(ratios.All(double.IsFinite));
            Assert.That(new[] { 0.01, 0.02, 0.05, 1.0 }.Select(i => model.LogExplainedPeak(i, PeakIntensityRatioModel.MaxCoverage)),
                Has.All.GreaterThanOrEqualTo(model.LogMissingFragment));

            // the spectrum's own chance-peak coverage is q; one recorded without falls back to the training mean
            Assert.That(model.CoverageOf(ManyFaintAgainstFewFaint()), Is.EqualTo(0.05).Within(1e-6));
            Assert.That(model.CoverageOf(Evidence(new[] { OnS3 }, (0.5, new[] { 0 }))), Is.EqualTo(model.MeanCoverage));
            Assert.That(model.MeanCoverage, Is.EqualTo(0.05).Within(1e-6));
        }

        [Test]
        public static void ManyFaintPeaksBeatFewUnderTheIntensityRatioButLoseUnderPtminer()
        {
            LocalizationEvidence evidence = ManyFaintAgainstFewFaint();

            // PTMiner's likelihood counts every faint peak explained against a placement, so the one explaining fewer wins
            var cdf = new MatchedPeakIntensityModel(Enumerable.Range(1, 50).Select(i => i / 50.0));
            double[] ptminer = LocalizationProbabilityEngine.LogLikelihoods(evidence, cdf);
            Assert.That(ptminer[1], Is.GreaterThan(ptminer[0] + 10));

            double[] ratio = LocalizationProbabilityEngine.LogLikelihoodRatios(evidence, TrainedRatioModel());
            Assert.That(ratio[0], Is.GreaterThan(ratio[1] + 10));

            // and so through the engine, learning both models from the unmodified matches
            var modified = Psm(evidence, modified: true);
            var psms = UnmodifiedRatioMatches().Append(modified).ToList();

            Run(psms, likelihood: LocalizationLikelihood.PtminerCdf);
            Assert.That(modified.LocalizationProbabilities.BestPlacement, Is.SameAs(OnT4));

            LocalizationProbabilityResults results = Run(psms, likelihood: LocalizationLikelihood.IntensityRatio);
            Assert.That(results.Likelihood, Is.EqualTo(LocalizationLikelihood.IntensityRatio));
            Assert.That(results.FragmentObservationProbability, Is.EqualTo(0.5).Within(1e-12));
            Assert.That(modified.LocalizationProbabilities.BestPlacement, Is.SameAs(OnS3));
            Assert.That(modified.LocalizationProbabilities.BestProbability, Is.GreaterThan(0.99));
        }

        [Test]
        public static void AnIntenseUnexplainedPeakStillDecides()
        {
            PeakIntensityRatioModel model = TrainedRatioModel();

            // each placement explains one differential peak: S3 an intense one, T4 a faint one
            var oneEach = RatioEvidence(new[] { OnS3, OnT4 }, new[] { 40, 40 }, FaintNoise(), (0.9, new[] { 0 }), (0.02, new[] { 1 }));
            double[] oneEachRatios = LocalizationProbabilityEngine.LogLikelihoodRatios(oneEach, model);
            Assert.That(oneEachRatios[0], Is.GreaterThan(oneEachRatios[1] + 5));

            // the intense peak outweighs three faint ones
            var threeFaint = RatioEvidence(new[] { OnS3, OnT4 }, new[] { 40, 40 }, FaintNoise(),
                (0.9, new[] { 0 }), (0.02, new[] { 1 }), (0.03, new[] { 1 }), (0.015, new[] { 1 }));
            double[] threeFaintRatios = LocalizationProbabilityEngine.LogLikelihoodRatios(threeFaint, model);
            Assert.That(threeFaintRatios[0], Is.GreaterThan(threeFaintRatios[1]));

            // and the same peaks explained, the less intense one is worth less: explaining a peak grows with its intensity
            Assert.That(model.LogExplainedPeak(0.9, 0.05), Is.GreaterThan(model.LogExplainedPeak(0.2, 0.05)));
            Assert.That(model.LogExplainedPeak(0.2, 0.05), Is.GreaterThan(model.LogExplainedPeak(0.02, 0.05)));

            // a busier spectrum (more chance peaks) makes the same match weaker evidence
            Assert.That(model.LogExplainedPeak(0.2, 0.2), Is.LessThan(model.LogExplainedPeak(0.2, 0.05)));
        }

        [Test]
        public static void EachPredictedFragmentLeftUnobservedCostsLogOneMinusD()
        {
            PeakIntensityRatioModel model = TrainedRatioModel(); // d = 0.5
            Assert.That(model.ObservationProbability, Is.EqualTo(0.5));

            // the same peaks explained, but S3 predicts 10 more fragments within range than T4, none of them seen
            var evidence = RatioEvidence(new[] { OnS3, OnT4 }, new[] { 30, 20 }, FaintNoise(),
                (0.6, new[] { 0, 1 }), (0.3, new[] { 0, 1 }), (0.005, new[] { 0, 1 })); // the last is below the floor
            double[] ratios = LocalizationProbabilityEngine.LogLikelihoodRatios(evidence, model);
            Assert.That(ratios[0], Is.EqualTo((30 - 2) * Math.Log(0.5)).Within(1e-9)); // shared peaks count only as observed
            Assert.That(ratios[1], Is.EqualTo((20 - 2) * Math.Log(0.5)).Within(1e-9));
            Assert.That(ratios[0] - ratios[1], Is.EqualTo(10 * Math.Log(0.5)).Within(1e-9));

            // a differential peak adds its explained-peak term to its candidate and removes one missing fragment
            var oneMore = RatioEvidence(new[] { OnS3, OnT4 }, new[] { 20, 20 }, FaintNoise(), (0.6, new[] { 0, 1 }), (0.3, new[] { 0 }));
            double[] oneMoreRatios = LocalizationProbabilityEngine.LogLikelihoodRatios(oneMore, model);
            Assert.That(oneMoreRatios[0] - oneMoreRatios[1],
                Is.EqualTo(model.LogExplainedPeak((float)0.3, 0.05) - Math.Log(0.5)).Within(1e-6));

            // evidence that recorded no predicted counts has no missing-fragment term
            var uncounted = Evidence(new[] { OnS3, OnT4 }, (0.6, new[] { 0, 1 }));
            Assert.That(LocalizationProbabilityEngine.LogLikelihoodRatios(uncounted, model), Is.EqualTo(new[] { 0.0, 0.0 }));
        }

        [Test]
        public static void TheIntensityRatioFallsBackToPtminerWhenNothingRecordedUnexplainedPeaks()
        {
            var target = Psm(Evidence(new[] { OnS3, OnT4 }, (0.8, new[] { 0 })), modified: true);
            LocalizationProbabilityResults results = Run(UnmodifiedMatches().Append(target).ToList(), likelihood: LocalizationLikelihood.IntensityRatio);

            Assert.That(results.Likelihood, Is.EqualTo(LocalizationLikelihood.PtminerCdf));
            Assert.That(target.LocalizationProbabilities.BestPlacement, Is.SameAs(OnS3));
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
        public static void TheLogNormalDistributionIsFittedToLogIntensities()
        {
            // log intensities -3, -2, -1: mean -2, sample standard deviation 1
            var model = new MatchedPeakIntensityModel(new[] { Math.Exp(-3), Math.Exp(-2), Math.Exp(-1) },
                distribution: MatchedPeakIntensityDistribution.LogNormal);

            Assert.That(model.CumulativeProbability(Math.Exp(-2)), Is.EqualTo(0.5).Within(1e-12));
            Assert.That(model.CumulativeProbability(Math.Exp(-1)), Is.EqualTo(MathNet.Numerics.Distributions.Normal.CDF(0, 1, 1)).Within(1e-12));
            Assert.That(model.CumulativeProbability(1e-300), Is.GreaterThan(0)); // never exactly 0 or 1
            Assert.That(model.CumulativeProbability(1e300), Is.LessThan(1));
        }

        [Test]
        public static void WithoutAPseudocountTheLearnedPriorIsSharperAndStillFinite()
        {
            // phospho is seen on S3 only; with no pseudocount T4's prior shrinks toward zero, so an unclear match goes to S3
            var clear = Enumerable.Range(0, 40).Select(_ => Psm(Evidence(new[] { OnS3, OnT4 }, (0.8, new[] { 0 }), (0.6, new[] { 0 })), modified: true));
            var unclear = Psm(Evidence(new[] { OnS3, OnT4 }, (0.3, new[] { 0, 1 })), modified: true);
            var psms = UnmodifiedMatches().Concat(clear).Append(unclear).ToList();

            Run(psms);
            double withPseudocount = unclear.LocalizationProbabilities.BestProbability;

            new LocalizationProbabilityEngine(psms, Parameters, new List<(string, CommonParameters)>(), new List<string>()) { Pseudocount = 0 }.Run();

            Assert.That(unclear.LocalizationProbabilities.BestPlacement, Is.SameAs(OnS3));
            Assert.That(unclear.LocalizationProbabilities.BestProbability, Is.GreaterThan(withPseudocount));
            Assert.That(unclear.LocalizationProbabilities.PlacementProbabilities.All(p => !double.IsNaN(p)));
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
        [NonParallelizable] // listens to every engine that finishes
        [TestCase(true)]
        [TestCase(false)]
        public static void AClosedSearchComputesLocalizationProbabilitiesOnlyWhenAsked(bool doLocalizationProbabilities)
        {
            var results = new List<LocalizationProbabilityResults>();
            void Collect(object sender, SingleEngineFinishedEventArgs e)
            {
                if (e.MyResults is LocalizationProbabilityResults r)
                {
                    lock (results) { results.Add(r); }
                }
            }

            var task = new TaskLayer.SearchTask
            {
                SearchParameters = new TaskLayer.SearchParameters { DoLocalizationProbabilities = doLocalizationProbabilities, WriteMzId = false }
            };
            string folder = System.IO.Path.Combine(TestContext.CurrentContext.TestDirectory, $"LocalizationProbabilitySearch{doLocalizationProbabilities}");
            System.IO.Directory.CreateDirectory(folder);
            string spectra = System.IO.Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SmallCalibratible_Yeast.mzML");
            string database = System.IO.Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\smalldb.fasta");

            MetaMorpheusEngine.FinishedSingleEngineHandler += Collect;
            try
            {
                task.RunTask(folder, new List<EngineLayer.DatabaseLoading.DbForTask> { new(database, false) }, new List<string> { spectra }, "");
            }
            finally
            {
                MetaMorpheusEngine.FinishedSingleEngineHandler -= Collect;
                System.IO.Directory.Delete(folder, true);
            }

            if (!doLocalizationProbabilities)
            {
                Assert.That(results, Is.Empty);
                return;
            }
            Assert.That(results, Has.Count.EqualTo(1));
            TestContext.WriteLine(results[0]);
            // the unmodified matches of a closed search are enough to learn peak intensities from
            Assert.That(results[0].TrainingPeaks, Is.GreaterThanOrEqualTo(LocalizationProbabilityEngine.DefaultMinTrainingPeaks));
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
