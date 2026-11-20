using Chromatography;
using Chromatography.RetentionTimePrediction;
using EngineLayer;
using EngineLayer.ClassicSearch;
using EngineLayer.FdrAnalysis;
using NUnit.Framework;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace Test
{
    /// <summary>
    /// Top-down PEP uses CZE electrophoretic mobility as its separation feature, and nothing else: no
    /// retention-time predictor (SSRCalc, Chronologer, Prosit) is ever run on a proteoform.
    /// </summary>
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    internal class TopDownPepSeparationTests
    {
        /// <summary>Counts every call on both entry points, so a test can prove the predictor was never consulted.</summary>
        private sealed class CountingPredictor : IRetentionTimePredictor
        {
            private int _calls;
            internal int Calls => _calls;

            public string PredictorName => "Counting";
            public SeparationType SeparationType => SeparationType.HPLC;

            public double? PredictRetentionTimeEquivalent(IRetentionPredictable peptide, out RetentionTimeFailureReason? failureReason)
            {
                Interlocked.Increment(ref _calls);
                failureReason = null;
                return peptide.FullSequence.Length;
            }

            public IReadOnlyList<(double? PredictedValue, IRetentionPredictable Peptide, RetentionTimeFailureReason? FailureReason)>
                PredictRetentionTimeEquivalents(IEnumerable<IRetentionPredictable> peptides, int maxThreads = 1)
            {
                Interlocked.Increment(ref _calls);
                return peptides.Select(p => ((double?)p.FullSequence.Length, p, (RetentionTimeFailureReason?)null)).ToList();
            }

            public string GetFormattedSequence(IRetentionPredictable peptide, out RetentionTimeFailureReason? failureReason)
            {
                failureReason = null;
                return peptide.FullSequence;
            }

            public void Dispose() { }
        }

        [TestCase("HPLC")]
        [TestCase("CZE")]
        public static void TopDownRunsNoRetentionTimePredictor(string separationType)
        {
            var (psms, fsp) = SearchTopDown(separationType);
            var predictor = new CountingPredictor();

            var engine = new PepAnalysisEngine(psms, "top-down", fsp, PepOutputFolder(), predictor);
            engine.ComputePEPValuesForAllPSMs();

            Assert.That(predictor.Calls, Is.Zero);
            Assert.That(engine.RetentionTimePredictor, Is.SameAs(predictor), "the predictor must not be warmed");
            Assert.That(engine.FileSpecificTimeDependantHydrophobicityAverageAndDeviation_unmodified, Is.Null);
            Assert.That(engine.FileSpecificTimeDependantHydrophobicityAverageAndDeviation_modified, Is.Null);
        }

        /// <summary>
        /// Every match from an LC file carries the missing value, so the feature is constant there and the model
        /// cannot split on it.
        /// </summary>
        [Test]
        public static void TopDownLcFileCarriesTheMissingValue()
        {
            var (psms, fsp) = SearchTopDown("HPLC");
            var engine = new PepAnalysisEngine(psms, "top-down", fsp, PepOutputFolder());

            List<PsmData> features = FeaturesOfEveryHypothesis(engine, psms);

            Assert.That(features, Is.Not.Empty);
            Assert.That(features.Select(f => f.HydrophobicityZScore), Is.All.EqualTo(100));
            Assert.That(features.Select(f => f.HasHydrophobicity), Is.All.EqualTo(0));
        }

        [Test]
        public static void TopDownCzeFileTrainsOnMobility()
        {
            var (psms, fsp) = SearchTopDown("CZE");
            var engine = new PepAnalysisEngine(psms, "top-down", fsp, PepOutputFolder());

            Assert.That(engine.TrainingVariables, Does.Contain("HydrophobicityZScore").And.Contain("HasHydrophobicity"));
            Assert.That(engine.FileSpecificTimeDependantHydrophobicityAverageAndDeviation_CZE[DataFileName], Is.Not.Empty,
                "confident targets must build the mobility reference");

            List<PsmData> features = FeaturesOfEveryHypothesis(engine, psms);

            List<PsmData> withMobility = features.Where(f => f.HasHydrophobicity == 1).ToList();
            Assert.That(withMobility, Is.Not.Empty);
            Assert.That(withMobility.Select(f => f.HydrophobicityZScore), Is.All.InRange(0, 100));
        }

        private const string DataFileName = "TaGe_SA_HeLa_04_subset_longestSeq.mzML";

        private static string PepOutputFolder() => Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\");

        private static List<PsmData> FeaturesOfEveryHypothesis(PepAnalysisEngine engine, List<SpectralMatch> psms) => psms
            .SelectMany(psm => psm.BestMatchingBioPolymersWithSetMods
                .Select(h => engine.CreateOnePsmDataEntry("top-down", psm, h, !h.IsDecoy)))
            .ToList();

        /// <summary>The top-down search <see cref="FdrTest.TestComputePEPValueTopDown"/> runs, with the file's separation type set.</summary>
        private static (List<SpectralMatch> Psms, List<(string fileName, CommonParameters fileSpecificParameters)> Fsp) SearchTopDown(string separationType)
        {
            string origDataFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\" + DataFileName);
            var commonParameters = new CommonParameters(
                digestionParams: new DigestionParams(protease: "top-down"),
                scoreCutoff: 1,
                assumeOrphanPeaksAreZ1Fragments: false,
                separationType: separationType);
            var fsp = new List<(string fileName, CommonParameters fileSpecificParameters)>
            {
                (DataFileName, commonParameters)
            };

            var myMsDataFile = new MyFileManager(true).LoadFile(origDataFile, commonParameters);
            List<Protein> proteinList = ProteinDbLoader.LoadProteinFasta(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\HeLaFakeTopDown.fasta"), true, DecoyType.Reverse, false, out _,
                ProteinDbLoader.UniprotAccessionRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotGeneNameRegex,
                ProteinDbLoader.UniprotOrganismRegex, -1);
            var scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, origDataFile, new CommonParameters()).OrderBy(b => b.PrecursorMass).ToArray();
            SpectralMatch[] allPsms = new PeptideSpectralMatch[scans.Length];
            new ClassicSearchEngine(allPsms, scans, new List<Modification>(), new List<Modification>(), null, null, null,
                proteinList, new SinglePpmAroundZeroSearchMode(5), commonParameters, fsp, null, new List<string>(), false).Run();
            List<SpectralMatch> psms = allPsms.Where(p => p != null).ToList();
            new FdrAnalysisEngine(psms, 1, commonParameters, fsp, new List<string>(), doPEP: false).Run();
            return (psms, fsp);
        }
    }
}
