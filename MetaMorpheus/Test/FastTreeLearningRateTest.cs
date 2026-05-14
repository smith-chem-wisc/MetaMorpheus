using EngineLayer;
using EngineLayer.ClassicSearch;
using EngineLayer.FdrAnalysis;
using Microsoft.ML.Trainers.FastTree;
using NUnit.Framework;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace Test
{
    /// <summary>
    /// Integration test for the tuned FastTree PEP learning-rate default. Searches four ~2-minute
    /// snips cut from the retention-time heart of the BottomUp runs (a medium ~1,000-MS2 dataset),
    /// then runs the PEP pipeline twice - at the historical learning rate 0.2 and the tuned default
    /// 0.05 - and shows the tuned rate is markedly better calibrated (lower LogLoss) at no cost in
    /// target peptides at PEP q-value &lt; 0.01.
    ///
    /// The full hyperparameter sweep and the cross-dataset-size evidence behind the 0.2 -> 0.05
    /// change are written up in docs/PEP_FastTree_Options_Methodology.md.
    /// </summary>
    [TestFixture]
    public static class FastTreeLearningRateTest
    {
        private static readonly string[] SnipFiles =
        {
            "20100609_Velos1_TaGe_SA_293_2_snip_16173-16526.mzML",
            "20100609_Velos1_TaGe_SA_293_3_snip_14279-14594.mzML",
            "20100609_Velos1_TaGe_SA_293_4_snip_11989-12231.mzML",
            "20100609_Velos1_TaGe_SA_293_5_snip_11883-12145.mzML",
        };
        private const double PepQValueThreshold = 0.01;

        private static CommonParameters _commonParameters;
        private static List<(string, CommonParameters)> _fsp;
        private static List<MassSpectrometry.MsDataFile> _msDataFiles;
        private static List<string> _filePaths;
        private static List<Protein> _proteinList;

        [OneTimeSetUp]
        public static void SetUp()
        {
            string dir = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "FastTreeSnips");
            _commonParameters = new CommonParameters(digestionParams: new DigestionParams());
            _filePaths = SnipFiles.Select(f => Path.Combine(dir, f)).ToList();
            _fsp = _filePaths.Select(p => (p, _commonParameters)).ToList();

            var fileManager = new MyFileManager(true);
            _msDataFiles = _filePaths.Select(p => fileManager.LoadFile(p, _commonParameters)).ToList();

            _proteinList = ProteinDbLoader.LoadProteinXML(
                Path.Combine(dir, "BottomUpSnips_prunedDb.xml"),
                true, DecoyType.Reverse, new List<Modification>(), false, new List<string>(), out _);
        }

        [Test]
        public static void LowerLearningRate_ImprovesCalibration_AtNoCostInPeptides()
        {
            var high = EvaluatePep(learningRate: 0.2);   // the historical default
            var low = EvaluatePep(learningRate: 0.05);   // the tuned default (PepAnalysisEngine.BGDTreeOptions)

            TestContext.WriteLine(
                $"lr=0.20: {high.peptideCount} target peptides (PEP q<{PepQValueThreshold}), LogLoss {high.logLoss:F4}");
            TestContext.WriteLine(
                $"lr=0.05: {low.peptideCount} target peptides (PEP q<{PepQValueThreshold}), LogLoss {low.logLoss:F4}");

            // The tuned (gentler) learning rate produces a markedly better-calibrated PEP model.
            Assert.That(low.logLoss, Is.LessThan(high.logLoss),
                $"lr=0.05 LogLoss ({low.logLoss:F4}) should beat lr=0.2 LogLoss ({high.logLoss:F4})");

            // ...and it does so without costing target-peptide identifications on a realistically sized dataset.
            Assert.That(low.peptideCount, Is.GreaterThanOrEqualTo(high.peptideCount),
                $"lr=0.05 peptide count ({low.peptideCount}) should be >= lr=0.2 count ({high.peptideCount})");
        }

        /// <summary>
        /// Runs a fresh classic search of the four snips, then the PEP pipeline with the given FastTree
        /// learning rate, and returns the target-peptide count at PEP q-value &lt; threshold plus the
        /// trained model's LogLoss. A fresh search per call keeps the two learning-rate runs independent.
        /// </summary>
        private static (int peptideCount, double logLoss) EvaluatePep(double learningRate)
        {
            var psms = RunFreshSearch();
            new FdrAnalysisEngine(psms, 1, _commonParameters, _fsp, new List<string>(), doPEP: false).Run();

            var treeOptions = new FastTreeBinaryTrainer.Options
            {
                NumberOfTrees = 400,
                LearningRate = learningRate,
                NumberOfLeaves = 20,
                MinimumExampleCountPerLeaf = 10,
            };
            var pepEngine = new PepAnalysisEngine(psms, "standard", _fsp,
                TestContext.CurrentContext.TestDirectory, customTreeOptions: treeOptions);
            string metrics = pepEngine.ComputePEPValuesForAllPSMs();
            Assert.That(metrics, Does.Not.Contain("failed"), "PEP training failed for learning rate " + learningRate);

            var peptides = psms
                .OrderBy(p => p.FdrInfo.PEP)
                .ThenByDescending(p => p)
                .GroupBy(p => p.FullSequence)
                .Select(g => g.First())
                .ToList();
            FdrAnalysisEngine.CalculateQValue(peptides, peptideLevelCalculation: true, pepCalculation: true);

            int count = peptides.Count(p => !p.IsDecoy && !p.IsContaminant
                                            && p.PeptideFdrInfo.PEP_QValue < PepQValueThreshold);
            return (count, ParseLogLoss(metrics));
        }

        private static List<SpectralMatch> RunFreshSearch()
        {
            var allPsms = new List<SpectralMatch>();
            for (int i = 0; i < SnipFiles.Length; i++)
            {
                var ms2Scans = MetaMorpheusTask.GetMs2Scans(_msDataFiles[i], _filePaths[i], _commonParameters)
                    .OrderBy(b => b.PrecursorMass).ToArray();
                SpectralMatch[] fileResults = new PeptideSpectralMatch[ms2Scans.Length];
                new ClassicSearchEngine(fileResults, ms2Scans, new List<Modification>(), new List<Modification>(),
                    null, null, null, _proteinList, new SinglePpmAroundZeroSearchMode(5), _commonParameters, _fsp,
                    null, new List<string>(), false).Run();
                allPsms.AddRange(fileResults.Where(p => p != null));
            }
            return allPsms;
        }

        private static double ParseLogLoss(string metricsBlock)
        {
            foreach (string line in metricsBlock.Split('\n'))
            {
                int idx = line.IndexOf("LogLoss:", StringComparison.Ordinal);
                if (idx >= 0 && double.TryParse(line.Substring(idx + "LogLoss:".Length).Trim(), out double v))
                {
                    return v;
                }
            }
            return double.NaN;
        }
    }
}
