using EngineLayer.Calibration;
using MassSpectrometry;
using NUnit.Framework;
using SharpLearning.AdaBoost.Learners;
using SharpLearning.Common.Interfaces;
using SharpLearning.Containers.Matrices;
using SharpLearning.DecisionTrees.Learners;
using SharpLearning.GradientBoost.Learners;
using SharpLearning.Metrics.Regression;
using SharpLearning.RandomForest.Learners;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    [TestFixture]
    public static class RfCalibrationFunctionTest
    {
        //[Test]
        //public static void TestRfCalibrationFunction()
        //{
        //    RandomForestCalibrationFunction randomForestCalibrationFunction = new RandomForestCalibrationFunction(1, 0, new[] { true });

        //    List<TestInputsOutputs> trainingList = GenerateTrainingData();
        //    randomForestCalibrationFunction.Learn(trainingList);
        //    Assert.AreEqual(1, randomForestCalibrationFunction.Predict(new[] { 5.0 }));
        //    Assert.AreEqual(2, randomForestCalibrationFunction.Predict(new[] { 15.0 }));
        //}

        //[Test]
        //public static void Test2dRfCalibrationFunction()
        //{
        //    RandomForestCalibrationFunction randomForestCalibrationFunction = new RandomForestCalibrationFunction(1, 0, new[] { true, true });

        //    List<TestInputsOutputs2> trainingList = GenerateTrainingData2();
        //    randomForestCalibrationFunction.Train(trainingList);
        //    Assert.AreEqual(1, randomForestCalibrationFunction.Predict(new[] { 1.5, 1.5 }));
        //    Assert.AreEqual(2, randomForestCalibrationFunction.Predict(new[] { 1.5, 4.5 }));
        //    Assert.AreEqual(2, randomForestCalibrationFunction.Predict(new[] { 4.5, 1.5 }));
        //    Assert.AreEqual(3, randomForestCalibrationFunction.Predict(new[] { 4.5, 4.5 }));
        //}

        //[Test]
        //public static void Test2dRfCalibrationFunction2()
        //{
        //    RandomForestCalibrationFunction randomForestCalibrationFunction = new RandomForestCalibrationFunction(1, 0, new[] { true, true });

        //    List<TestInputsOutputs2> trainingList = GenerateTrainingData3();
        //    randomForestCalibrationFunction.Train(trainingList);
        //    Assert.AreEqual(1, randomForestCalibrationFunction.Predict(new[] { 1.5, 1.5 }));
        //    Assert.AreEqual(2, randomForestCalibrationFunction.Predict(new[] { 1.5, 4.5 }));
        //    Assert.AreEqual(2, randomForestCalibrationFunction.Predict(new[] { 4.5, 1.5 }));
        //    Assert.AreEqual(1, randomForestCalibrationFunction.Predict(new[] { 4.5, 4.5 }));
        //}

        //[Test]
        //public static void TestQuadraticFunctionCalibration()
        //{
        //    Random rand = new Random(1);
        //    RandomForestCalibrationFunction randomForestCalibrationFunction = new RandomForestCalibrationFunction(100, 2, new[] { true });

        //    List<TestInputsOutputs> trainingList = GenerateQuadraticTrainingData().ToList();

        //    randomForestCalibrationFunction.Train(trainingList);

        //    foreach (var hah in trainingList)
        //        Console.WriteLine(hah.Inputs[0] + " , " + hah.Label + " , " + randomForestCalibrationFunction.Predict(hah.Inputs));
        //}

        #region Public Methods

        [Test]
        public static void TestQuadraticFunctionCalibrationNICE()
        {
            List<TestInputsOutputs> trainingList = GenerateQuadraticTrainingData().ToList();
            F64Matrix observations = new F64Matrix(trainingList.Select(b => b.Inputs[0]).ToArray(), trainingList.Count, 1);
            var targets = trainingList.Select(b => b.Label).ToArray();

            var learners = new List<ILearner<double>>
            {
                new RegressionAdaBoostLearner(),
                new RegressionDecisionTreeLearner(),
                new RegressionAbsoluteLossGradientBoostLearner(),
                new RegressionGradientBoostLearner(),
                new RegressionHuberLossGradientBoostLearner(),
                new RegressionQuantileLossGradientBoostLearner(),
                new RegressionSquareLossGradientBoostLearner(),
                new RegressionExtremelyRandomizedTreesLearner(),
                new RegressionRandomForestLearner(),
                //new RegressionEnsembleLearner(),
                //new RegressionForwardSearchModelSelectingEnsembleLearner(),
                //new RegressionModelSelectingEnsembleLearner(),
                //new RegressionRandomModelSelectingEnsembleLearner(),
                //new RegressionStackingEnsembleLearner(),
                //new RegressionNeuralNetLearner(),
            };
            foreach (var learner in learners)
            {
                var model = learner.Learn(observations, targets);
                var predictions = new double[targets.Length];
                for (int i = 0; i < targets.Length; i++)
                    predictions[i] = model.Predict(observations.Row(i));
                var metric = new MeanSquaredErrorRegressionMetric();
                var trainError = metric.Error(targets, predictions);
                Console.WriteLine("mse for " + learner.GetType() + ": " + trainError);
            }
        }

        #endregion Public Methods

        #region Private Methods

        private static IEnumerable<TestInputsOutputs> GenerateQuadraticTrainingData()
        {
            for (int i = -100; i <= 100; i++)
                yield return new TestInputsOutputs(i, i * i);
        }

        private static List<TestInputsOutputs2> GenerateTrainingData3()
        {
            List<TestInputsOutputs2> l = new List<TestInputsOutputs2>
            {
                new TestInputsOutputs2(1, 1, 1),
                new TestInputsOutputs2(1.1, 1.1, 1),
                new TestInputsOutputs2(1, 2, 1),
                new TestInputsOutputs2(2, 1, 1),
                new TestInputsOutputs2(2, 2, 1),

                new TestInputsOutputs2(4, 1, 2),
                new TestInputsOutputs2(5, 2, 2),
                new TestInputsOutputs2(5, 1, 2),
                new TestInputsOutputs2(4, 2, 2),

                new TestInputsOutputs2(1, 4, 2),
                new TestInputsOutputs2(2, 5, 2),
                new TestInputsOutputs2(1, 5, 2),
                new TestInputsOutputs2(2, 4, 2),

                new TestInputsOutputs2(4, 4, 1),
                new TestInputsOutputs2(4, 5, 1),
                new TestInputsOutputs2(5, 4, 1),
                new TestInputsOutputs2(4.9, 4.9, 1),
                new TestInputsOutputs2(5, 5, 1),
            };

            return l;
        }

        private static List<TestInputsOutputs2> GenerateTrainingData2()
        {
            List<TestInputsOutputs2> l = new List<TestInputsOutputs2>
            {
                new TestInputsOutputs2(1, 1, 1),
                new TestInputsOutputs2(1, 2, 1),
                new TestInputsOutputs2(2, 1, 1),
                new TestInputsOutputs2(2, 2, 1),

                new TestInputsOutputs2(4, 1, 2),
                new TestInputsOutputs2(5, 2, 2),
                new TestInputsOutputs2(5, 1, 2),
                new TestInputsOutputs2(4, 2, 2),

                new TestInputsOutputs2(1, 4, 2),
                new TestInputsOutputs2(2, 5, 2),
                new TestInputsOutputs2(1, 5, 2),
                new TestInputsOutputs2(2, 4, 2),

                new TestInputsOutputs2(4, 4, 3),
                new TestInputsOutputs2(4, 5, 3),
                new TestInputsOutputs2(5, 4, 3),
                new TestInputsOutputs2(5, 5, 3),
            };

            return l;
        }

        private static List<TestInputsOutputs> GenerateTrainingData()
        {
            List<TestInputsOutputs> l = new List<TestInputsOutputs>
            {
                new TestInputsOutputs(1, 1),
                new TestInputsOutputs(2, 1),
                new TestInputsOutputs(3, 1),
                new TestInputsOutputs(4, 1),
                new TestInputsOutputs(6, 1),
                new TestInputsOutputs(7, 1),
                new TestInputsOutputs(8, 1),
                new TestInputsOutputs(9, 1),
                new TestInputsOutputs(11, 2),
                new TestInputsOutputs(12, 2),
                new TestInputsOutputs(13, 2),
                new TestInputsOutputs(14, 2),
                new TestInputsOutputs(16, 2),
                new TestInputsOutputs(17, 2),
                new TestInputsOutputs(18, 2),
                new TestInputsOutputs(19, 2),
            };

            return l;
        }

        #endregion Private Methods

        #region Private Classes

        private class TestInputsOutputs : IHasInputsAndOutputs
        {
            #region Public Constructors

            public TestInputsOutputs(double v1, double v2)
            {
                Inputs = new[] { v1 };
                Label = v2;
            }

            #endregion Public Constructors

            #region Public Properties

            public double[] Inputs { get; }

            public double Label { get; }

            #endregion Public Properties

            #region Public Methods

            public override string ToString()
            {
                return "" + Inputs[0] + ", " + Label;
            }

            #endregion Public Methods
        }

        private class TestInputsOutputs2 : IHasInputsAndOutputs
        {
            #region Public Constructors

            public TestInputsOutputs2(double v1, double v2, double label)
            {
                Inputs = new[] { v1, v2 };
                Label = label;
            }

            #endregion Public Constructors

            #region Public Properties

            public double[] Inputs { get; }

            public double Label { get; }

            #endregion Public Properties

            #region Public Methods

            public override string ToString()
            {
                return "(" + Inputs[0] + ", " + Inputs[1] + ") " + Label;
            }

            #endregion Public Methods
        }

        private class CalibTestDataFile : IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>>
        {
            #region Private Fields

            private List<TestInputsOutputs2> list;

            #endregion Private Fields

            #region Public Constructors

            public CalibTestDataFile(List<TestInputsOutputs2> list)
            {
                this.list = list;
            }

            #endregion Public Constructors

            #region Public Properties

            public int NumSpectra => throw new NotImplementedException();

            #endregion Public Properties

            #region Public Methods

            public int GetClosestOneBasedSpectrumNumber(double retentionTime)
            {
                throw new NotImplementedException();
            }

            public IEnumerator<IMsDataScan<IMzSpectrum<IMzPeak>>> GetEnumerator()
            {
                throw new NotImplementedException();
            }

            public IEnumerable<IMsDataScan<IMzSpectrum<IMzPeak>>> GetMsScansInTimeRange(double firstRT, double lastRT)
            {
                throw new NotImplementedException();
            }

            public IMsDataScan<IMzSpectrum<IMzPeak>> GetOneBasedScan(int oneBasedScanNumber)
            {
                throw new NotImplementedException();
            }

            IEnumerator IEnumerable.GetEnumerator()
            {
                throw new NotImplementedException();
            }

            #endregion Public Methods
        }

        #endregion Private Classes
    }
}