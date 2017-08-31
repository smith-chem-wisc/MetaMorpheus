using EngineLayer.Calibration;
using NUnit.Framework;
using System;
using System.Collections.Generic;

namespace Test
{
    [TestFixture]
    public static class RfCalibrationFunctionTest
    {
        #region Public Methods

        [Test]
        public static void TestRfCalibrationFunction()
        {
            Random rand = new Random(1);
            RandomForestCalibrationFunction randomForestCalibrationFunction = new RandomForestCalibrationFunction(1, 0, new[] { true }, rand);

            List<TestInputsOutputs> trainingList = GenerateTrainingData();
            randomForestCalibrationFunction.Train(trainingList);
            Assert.AreEqual(1, randomForestCalibrationFunction.Predict(new[] { 5.0 }));
            Assert.AreEqual(2, randomForestCalibrationFunction.Predict(new[] { 15.0 }));
        }

        [Test]
        public static void Test2dRfCalibrationFunction()
        {
            Random rand = new Random(1);
            RandomForestCalibrationFunction randomForestCalibrationFunction = new RandomForestCalibrationFunction(1, 0, new[] { true, true }, rand);

            List<TestInputsOutputs2> trainingList = GenerateTrainingData2();
            randomForestCalibrationFunction.Train(trainingList);
            Assert.AreEqual(1, randomForestCalibrationFunction.Predict(new[] { 1.5, 1.5 }));
            Assert.AreEqual(2, randomForestCalibrationFunction.Predict(new[] { 1.5, 4.5 }));
            Assert.AreEqual(2, randomForestCalibrationFunction.Predict(new[] { 4.5, 1.5 }));
            Assert.AreEqual(3, randomForestCalibrationFunction.Predict(new[] { 4.5, 4.5 }));
        }

        [Test]
        public static void Test2dRfCalibrationFunction2()
        {
            Random rand = new Random(1);
            RandomForestCalibrationFunction randomForestCalibrationFunction = new RandomForestCalibrationFunction(1, 0, new[] { true, true }, rand);

            List<TestInputsOutputs2> trainingList = GenerateTrainingData3();
            randomForestCalibrationFunction.Train(trainingList);
            Assert.AreEqual(1, randomForestCalibrationFunction.Predict(new[] { 1.5, 1.5 }));
            Assert.AreEqual(2, randomForestCalibrationFunction.Predict(new[] { 1.5, 4.5 }));
            Assert.AreEqual(2, randomForestCalibrationFunction.Predict(new[] { 4.5, 1.5 }));
            Assert.AreEqual(1, randomForestCalibrationFunction.Predict(new[] { 4.5, 4.5 }));
        }

        #endregion Public Methods

        #region Private Methods

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

        #endregion Private Classes
    }
}