using EngineLayer;
using EngineLayer.FdrAnalysis;
using Microsoft.ML.Trainers.FastTree;
using NUnit.Framework;
using System.Collections.Generic;

namespace Test
{
    /// <summary>
    /// Tests for the FastTree options used by the PEP model: the tuned default learning rate, and
    /// the constructor seam that lets a caller inject custom options. The learning-rate default was
    /// lowered from 0.2 to 0.05 after a hyperparameter sweep - 0.2 left the model overconfident
    /// (stochastic infinite LogLoss) and iterative retraining degraded its calibration, while 0.05
    /// roughly halves LogLoss. See the fastTreeOptions methodology doc for the evidence.
    /// </summary>
    [TestFixture]
    public static class FastTreeOptionsTest
    {
        private static List<(string, CommonParameters)> Fsp() =>
            new List<(string, CommonParameters)> { ("test.mzML", new CommonParameters()) };

        [Test]
        public static void DefaultFastTreeOptions_UseTunedLearningRate()
        {
            var engine = new PepAnalysisEngine(new List<SpectralMatch>(), "standard", Fsp(), outputFolder: null);
            var options = engine.BGDTreeOptions;

            // the tuned default - lowered from 0.2 for calibration (see class summary)
            Assert.That(options.LearningRate, Is.EqualTo(0.05));
            Assert.That(options.NumberOfTrees, Is.EqualTo(400));
            Assert.That(options.NumberOfLeaves, Is.EqualTo(20));
            Assert.That(options.MinimumExampleCountPerLeaf, Is.EqualTo(10));

            // fixed plumbing fields are always enforced
            Assert.That(options.LabelColumnName, Is.EqualTo("Label"));
            Assert.That(options.FeatureColumnName, Is.EqualTo("Features"));
            Assert.That(options.NumberOfThreads, Is.EqualTo(1));
        }

        [Test]
        public static void CustomFastTreeOptions_AreUsed_WithFixedFieldsStillEnforced()
        {
            var custom = new FastTreeBinaryTrainer.Options
            {
                NumberOfTrees = 123,
                LearningRate = 0.07,
                NumberOfLeaves = 17,
                MinimumExampleCountPerLeaf = 9,
            };
            var engine = new PepAnalysisEngine(new List<SpectralMatch>(), "standard", Fsp(),
                outputFolder: null, customTreeOptions: custom);
            var options = engine.BGDTreeOptions;

            // the caller's tunable knobs come through
            Assert.That(options.NumberOfTrees, Is.EqualTo(123));
            Assert.That(options.LearningRate, Is.EqualTo(0.07));
            Assert.That(options.NumberOfLeaves, Is.EqualTo(17));
            Assert.That(options.MinimumExampleCountPerLeaf, Is.EqualTo(9));

            // the engine still enforces the fixed plumbing fields on injected options, so a caller
            // only has to supply the tunable knobs
            Assert.That(options.LabelColumnName, Is.EqualTo("Label"));
            Assert.That(options.FeatureColumnName, Is.EqualTo("Features"));
            Assert.That(options.NumberOfThreads, Is.EqualTo(1));
            Assert.That(options.Seed, Is.EqualTo(42));
            Assert.That(options.FeatureSelectionSeed, Is.EqualTo(42));
        }
    }
}
